/* TeNeS - Massively parallel tensor network solver /
/ Copyright (C) 2019- The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify /
/ it under the terms of the GNU General Public License as published by /
/ the Free Software Foundation, either version 3 of the License, or /
/ (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, /
/ but WITHOUT ANY WARRANTY; without even the implied warranty of /
/ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the /
/ GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License /
/ along with this program. If not, see http://www.gnu.org/licenses/. */

/*! @file
 *  @brief Graded counterparts of the mptensor operations.
 *
 *  Each function mirrors its mptensor namesake — tensordot(), transpose(),
 *  trace(), conj(), qr(), svd(), svd_trunc(), slice(), extend(),
 *  reshape() — but acts on ::tenes::fermion::ftensor and generates the
 *  fermionic signs from the parity ledgers. The decompositions additionally
 *  preserve the grading: qr() and svd() sort the fused matrix legs
 *  even-first, factorize the even and odd diagonal blocks separately, and
 *  hand each factor a parity ledger for the new internal leg.
 *
 *  The header also hosts wrap_twosite_gate(), the shared loading adapter for
 *  two-site evolution gates and bundled-k measurement blobs. There are two
 *  loading conventions for measurement/evolution operators:
 *
 *  | path                       | loader              | swap applied |
 *  |----------------------------|---------------------|--------------|
 *  | two-site evolution / blob  | wrap_twosite_gate() | input legs   |
 *  | one-site operators         | plain ftensor wrap  | none         |
 */

#ifndef TENES_SRC_FERMION_FOPS_HPP_
#define TENES_SRC_FERMION_FOPS_HPP_

#include <algorithm>
#include <complex>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../mpi.hpp"
#include "finite_check.hpp"
#include "ftensor.hpp"
#include "sign_sweep.hpp"

namespace tenes {
namespace fermion {
namespace detail {

//! True iff ax appears in axes.
inline bool contains_axis(const mptensor::Axes& axes, std::size_t ax) {
  for (std::size_t i = 0; i < axes.size(); ++i) {
    if (axes[i] == ax) {
      return true;
    }
  }
  return false;
}

//! Diagnostic knob: TENES_FERMION_SVD_LOG_LIMIT caps how many svd_trunc()
//! calls dump their per-sector singular values to stderr (default 0: none).
inline int fermion_svd_log_limit() {
  const char* raw = std::getenv("TENES_FERMION_SVD_LOG_LIMIT");
  if (raw == nullptr) {
    return 0;
  }
  return std::max(0, std::atoi(raw));
}

//! Leg order mptensor::tensordot effectively gives its LEFT operand:
//! free legs first (in order), then the contracted legs in axes order.
inline mptensor::Axes tensordot_left_perm(std::size_t rank,
                                          const mptensor::Axes& axes) {
  mptensor::Axes perm;
  for (std::size_t ax = 0; ax < rank; ++ax) {
    if (!contains_axis(axes, ax)) {
      perm.push(ax);
    }
  }
  for (std::size_t i = 0; i < axes.size(); ++i) {
    perm.push(axes[i]);
  }
  return perm;
}

//! Leg order mptensor::tensordot effectively gives its RIGHT operand:
//! contracted legs first in REVERSED axes order, then the free legs.
//! The reversal is what makes the two operands' contracted legs meet
//! pairwise, and is the origin of the doubly-odd sign the gate wrappers
//! compensate.
inline mptensor::Axes tensordot_right_perm(std::size_t rank,
                                           const mptensor::Axes& axes) {
  mptensor::Axes perm;
  for (std::size_t i = axes.size(); i > 0; --i) {
    perm.push(axes[i - 1]);
  }
  for (std::size_t ax = 0; ax < rank; ++ax) {
    if (!contains_axis(axes, ax)) {
      perm.push(ax);
    }
  }
  return perm;
}

//! Copy of a with the Koszul sign of the permutation applied to the
//! elements, but WITHOUT permuting the legs — the sign half of a graded
//! transpose, used before handing the tensors to mptensor::tensordot.
template <class tensor>
ftensor<tensor> apply_transpose_sign_mask(const ftensor<tensor>& a,
                                          const mptensor::Axes& axes) {
  ftensor<tensor> ret = a;
  detail::validate_axes(axes, ret.t.shape().size(),
                        "fermion sign sweep: transpose axes out of range");
  if (detail::is_identity_axes(axes)) {
    return ret;
  }
  apply_swap_form(ret, detail::transpose_sign_form(axes));
  return ret;
}

//! Throw unless each contracted leg pair carries identical parity ledgers.
inline void validate_contracted_parity(const leg_parities& parity_a,
                                       const leg_parities& parity_b,
                                       const mptensor::Axes& axes_a,
                                       const mptensor::Axes& axes_b) {
  if (axes_a.size() != axes_b.size()) {
    throw std::runtime_error(
        "fermion tensordot parity check: axis count mismatch");
  }
  for (std::size_t i = 0; i < axes_a.size(); ++i) {
    if (parity_a[axes_a[i]] != parity_b[axes_b[i]]) {
      throw std::runtime_error(
          "fermion tensordot parity check: contracted parity mismatch");
    }
  }
}

//! Ledgers of the legs NOT in contracted, in their original order.
inline leg_parities free_leg_parities(const leg_parities& parity,
                                      const mptensor::Axes& contracted) {
  leg_parities ret;
  for (std::size_t ax = 0; ax < parity.size(); ++ax) {
    if (!contains_axis(contracted, ax)) {
      ret.push_back(parity[ax]);
    }
  }
  return ret;
}

//! Fused ledger of the listed legs, in axes order (trivial even ledger of
//! dimension 1 for an empty list).
inline parity_vector fuse_axes(const leg_parities& parity,
                               const mptensor::Axes& axes) {
  if (axes.size() == 0) {
    return parity_vector{false};
  }
  parity_vector ret = parity[axes[0]];
  for (std::size_t i = 1; i < axes.size(); ++i) {
    ret = fuse(ret, parity[axes[i]]);
  }
  return ret;
}

//! Dimensions of the listed legs, in axes order.
inline mptensor::Shape shape_from_axes(const mptensor::Shape& shape,
                                       const mptensor::Axes& axes) {
  mptensor::Shape ret;
  for (std::size_t i = 0; i < axes.size(); ++i) {
    ret.push(shape[axes[i]]);
  }
  return ret;
}

//! Product of all dimensions.
inline std::size_t product_shape(const mptensor::Shape& shape) {
  std::size_t ret = 1;
  for (std::size_t i = 0; i < shape.size(); ++i) {
    ret *= shape[i];
  }
  return ret;
}

//! Number of even index values in a ledger.
inline std::size_t count_even(const parity_vector& parity) {
  std::size_t ret = 0;
  for (bool p : parity) {
    if (!p) {
      ++ret;
    }
  }
  return ret;
}

/*!
 * @brief Debug check that a parity-sorted matrix is block diagonal.
 *
 * A parity-even tensor matricizes block-diagonally: an element is allowed
 * only when the fused row and column parities agree, so after the
 * even-first sort the (even row, odd col) and (odd row, even col) blocks
 * must vanish. qr()/svd() decompose the two diagonal blocks and never look
 * at the others, which would silently discard them if the input were not
 * even. Checked in debug builds only. Collective: the off-diagonal maximum
 * and the scale are reduced over the ranks before the comparison, so that
 * every rank throws or none does.
 *
 * @param[in] sorted Even-first sorted matricization.
 * @param[in] row_even Number of even rows.
 * @param[in] col_even Number of even columns.
 * @param[in] context Prefix of the error message.
 * @throw std::runtime_error If an off-diagonal block carries weight above
 *        1e-10 relative to the largest element (NDEBUG: never).
 */
template <class tensor>
void validate_block_diagonal(const tensor& sorted, std::size_t row_even,
                             std::size_t col_even, const char* context) {
#ifndef NDEBUG
  double off = 0.0;
  double scale = 0.0;
  mptensor::Index idx;
  idx.resize(sorted.shape().size());
  for (std::size_t n = 0; n < sorted.local_size(); ++n) {
    sorted.global_index_fast(n, idx);
    const double a = std::abs(sorted[n]);
    scale = std::max(scale, a);
    if ((idx[0] < row_even) != (idx[1] < col_even)) {
      off = std::max(off, a);
    }
  }
  std::vector<double> reduced{off, scale};
  tenes::allreduce_max(reduced, sorted.get_comm());
  off = reduced[0];
  scale = reduced[1];
  const double threshold = 1.0e-10 * std::max(1.0, scale);
  if (off > threshold) {
    std::stringstream ss;
    ss << context
       << ": input is not parity even; its off-diagonal parity blocks carry "
          "max_abs="
       << off << " (threshold " << threshold
       << ") and would be discarded silently";
    throw std::runtime_error(ss.str());
  }
#else
  static_cast<void>(sorted);
  static_cast<void>(row_even);
  static_cast<void>(col_even);
  static_cast<void>(context);
#endif
}

//! Complex conjugation that is a no-op for real scalars.
inline double scalar_conj(double v) { return v; }

//! @copybrief scalar_conj(double)
inline std::complex<double> scalar_conj(std::complex<double> v) {
  return std::conj(v);
}

}  // namespace detail

/*!
 * @brief Apply one swap-gate sign in place: @f$(-1)^{p_{ax1} p_{ax2}}@f$.
 *
 * For ax1 == ax2 the doubly-odd condition degenerates to "the leg is odd",
 * i.e. the parity operator @f$(-1)^p@f$ on that single leg (a case
 * SwapForm cannot express — its diagonal toggles are no-ops).
 *
 * @param[in,out] a Tensor to modify; parities are unchanged.
 * @param[in] ax1 First leg.
 * @param[in] ax2 Second leg (may equal @p ax1, see above).
 */
template <class tensor>
void apply_swap(ftensor<tensor>& a, int ax1, int ax2) {
  if (ax1 != ax2) {
    SwapForm form;
    form.toggle(ax1, ax2);
    apply_swap_form(a, form);
    return;
  }
  mptensor::Index idx;
  idx.resize(a.t.shape().size());
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    a.t.global_index_fast(n, idx);
    if (a.parity[ax1][idx[ax1]] && a.parity[ax2][idx[ax2]]) {
      a.t[n] = -a.t[n];
    }
  }
}

/*!
 * @brief Load a two-site evolution or bundled-k measurement operator.
 *
 * The gate is given as plain matrix elements
 * @f$\langle out_1\, out_2 | O | in_1\, in_2 \rangle@f$ in the ordered
 * two-site Fock basis; the result has legs
 * @f$(in_1, in_2, out_1, out_2)@f$ with parities (p1, p2, p1, p2) and the
 * input-leg swap mask pre-applied.
 *
 * The swap is a convention adapter, not extra physics: graded tensordot
 * contracts by moving the second operand's contracted legs to the front in
 * REVERSED order, which multiplies elements whose two input legs are both
 * odd by -1. Pre-applying the same mask cancels that factor, so the gate
 * that reaches theta is exactly the matrix the caller wrote. Without it
 * the doubly-odd input channel of every parity-conserving gate (e.g.
 * |11><11| in exp(-tau h)) is silently negated.
 *
 * @warning The bundled-k CTM blob and the simple-update kernel share this
 *          input-leg-only loading convention; see the file description.
 *
 * @param[in] op Gate matrix elements, legs (in1, in2, out1, out2).
 * @param[in] p1 Physical-leg ledger of the first site.
 * @param[in] p2 Physical-leg ledger of the second site.
 */
template <class tensor>
ftensor<tensor> wrap_twosite_gate(const tensor& op, const parity_vector& p1,
                                  const parity_vector& p2) {
  ftensor<tensor> fop{op, {p1, p2, p1, p2}};
  apply_swap(fop, 0, 1);
  return fop;
}

/*!
 * @brief Multiply the parity operator @f$(-1)^p@f$ onto one leg in place.
 */
template <class tensor>
void apply_parity(ftensor<tensor>& a, int ax) {
  std::vector<double> sign(a.parity[ax].size());
  for (std::size_t i = 0; i < sign.size(); ++i) {
    sign[i] = a.parity[ax][i] ? -1.0 : 1.0;
  }
  a.t.multiply_vector(sign, ax);
}

/*!
 * @brief Largest magnitude found in the parity-odd sector.
 *
 * A physical graded tensor must be parity even: every element whose index
 * selects an odd number of odd legs must vanish. Returns the worst
 * violation on the process-local slice (0 for a clean tensor); a
 * diagnostic, not a collective reduction.
 */
template <class tensor>
double parity_violation(const ftensor<tensor>& a) {
  double v = 0.0;
  mptensor::Index idx;
  idx.resize(a.t.shape().size());
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    a.t.global_index_fast(n, idx);
    if (count_odd(a.parity, idx) % 2 == 1) {
      v = std::max(v, std::abs(a.t[n]));
    }
  }
  return v;
}

/*!
 * @brief Graded transpose, returning a new tensor (see
 *        ftensor::transpose()).
 */
template <class tensor>
ftensor<tensor> transpose(const ftensor<tensor>& a,
                          const mptensor::Axes& axes) {
  ftensor<tensor> ret = a;
  ret.transpose(axes);
  return ret;
}

/*!
 * @brief Graded tensor contraction.
 *
 * Mirrors mptensor::tensordot: contracts axes_a of @p a against axes_b of
 * @p b pairwise; the result's legs are the free legs of @p a followed by
 * the free legs of @p b. The Koszul signs are those of the implied
 * reordering — free-then-contracted for @p a, reversed-contracted-then-free
 * for @p b (detail::tensordot_left_perm() / tensordot_right_perm()) — and
 * are applied as element masks before the plain contraction.
 *
 * @param[in] a Left operand.
 * @param[in] b Right operand.
 * @param[in] axes_a Legs of @p a to contract.
 * @param[in] axes_b Legs of @p b to contract, paired with @p axes_a in
 *            order.
 * @throw std::runtime_error If a contracted leg pair's ledgers differ.
 */
template <class tensor>
ftensor<tensor> tensordot(const ftensor<tensor>& a, const ftensor<tensor>& b,
                          const mptensor::Axes& axes_a,
                          const mptensor::Axes& axes_b) {
  detail::validate_contracted_parity(a.parity, b.parity, axes_a, axes_b);
  ftensor<tensor> a_masked = detail::apply_transpose_sign_mask(
      a, detail::tensordot_left_perm(a.parity.size(), axes_a));
  ftensor<tensor> b_masked = detail::apply_transpose_sign_mask(
      b, detail::tensordot_right_perm(b.parity.size(), axes_b));
  ftensor<tensor> ret;
  ret.t = mptensor::tensordot(a_masked.t, b_masked.t, axes_a, axes_b);
  ret.parity = detail::free_leg_parities(a.parity, axes_a);
  leg_parities b_free = detail::free_leg_parities(b.parity, axes_b);
  ret.parity.insert(ret.parity.end(), b_free.begin(), b_free.end());
  return ret;
}

/*!
 * @brief Graded full contraction of two tensors to a scalar.
 *
 * The all-legs-contracted case of tensordot(), evaluated without forming
 * the rank-0 intermediate (mirrors mptensor::trace).
 */
template <class tensor>
typename tensor::value_type trace(const ftensor<tensor>& a,
                                  const ftensor<tensor>& b,
                                  const mptensor::Axes& axes_a,
                                  const mptensor::Axes& axes_b) {
  detail::validate_contracted_parity(a.parity, b.parity, axes_a, axes_b);
  ftensor<tensor> a_masked = detail::apply_transpose_sign_mask(
      a, detail::tensordot_left_perm(a.parity.size(), axes_a));
  ftensor<tensor> b_masked = detail::apply_transpose_sign_mask(
      b, detail::tensordot_right_perm(b.parity.size(), axes_b));
  return mptensor::trace(a_masked.t, b_masked.t, axes_a, axes_b);
}

/*!
 * @brief Graded conjugate (the bra layer of a ket tensor).
 *
 * Elementwise complex conjugation times the sign
 * @f$(-1)^{m(m-1)/2}@f$, where @f$m@f$ is the number of odd legs the
 * element's index selects: reversing the @f$m@f$ fermionic factors of a
 * basis state costs @f$\binom{m}{2}@f$ exchanges. Parities are unchanged.
 */
template <class tensor>
ftensor<tensor> conj(const ftensor<tensor>& a) {
  ftensor<tensor> ret = a;
  mptensor::Index idx;
  idx.resize(ret.t.shape().size());
  for (std::size_t n = 0; n < ret.t.local_size(); ++n) {
    ret.t.global_index_fast(n, idx);
    const int m = count_odd(ret.parity, idx);
    const double sign = ((m * (m - 1) / 2) % 2 == 0) ? 1.0 : -1.0;
    ret.t[n] = sign * detail::scalar_conj(ret.t[n]);
  }
  return ret;
}

/*!
 * @brief Slice one leg to the index range [b, e), keeping its ledger slice.
 *
 * No signs arise: slicing does not reorder legs. Mirrors mptensor::slice.
 */
template <class tensor>
ftensor<tensor> slice(const ftensor<tensor>& a, std::size_t ax, std::size_t b,
                      std::size_t e) {
  ftensor<tensor> ret;
  ret.t = mptensor::slice(a.t, ax, b, e);
  ret.parity = a.parity;
  ret.parity[ax] =
      parity_vector(a.parity[ax].begin() + b, a.parity[ax].begin() + e);
  return ret;
}

/*!
 * @brief Zero-pad legs up to a larger shape; padded index values are even.
 *
 * Mirrors mptensor::extend.
 *
 * @throw std::runtime_error If any dimension of @p sh is smaller than the
 *        current one.
 */
template <class tensor>
ftensor<tensor> extend(const ftensor<tensor>& a, const mptensor::Shape& sh) {
  ftensor<tensor> ret;
  ret.t = mptensor::extend(a.t, sh);
  ret.parity = a.parity;
  for (std::size_t ax = 0; ax < ret.parity.size(); ++ax) {
    if (ret.parity[ax].size() > sh[ax]) {
      throw std::runtime_error("fermion extend: new shape is smaller");
    }
    ret.parity[ax].resize(sh[ax], false);
  }
  return ret;
}

/*!
 * @brief Graded reshape restricted to fusing groups of ADJACENT legs.
 *
 * Each output leg must be the product of one or more consecutive input
 * legs; their ledgers are combined with fuse() (column-major, matching
 * mptensor::reshape's flattening). Adjacent fusion generates no signs;
 * splitting a leg or fusing non-adjacent legs is not supported — reorder
 * with transpose() first, which is where the signs are accounted.
 *
 * @throw std::runtime_error If @p sh is not such an adjacent fusion of the
 *        input shape.
 */
template <class tensor>
ftensor<tensor> reshape(const ftensor<tensor>& a, const mptensor::Shape& sh) {
  ftensor<tensor> ret;
  ret.t = mptensor::reshape(a.t, sh);
  ret.parity.clear();
  std::size_t old_ax = 0;
  for (std::size_t new_ax = 0; new_ax < sh.size(); ++new_ax) {
    if (old_ax >= a.parity.size()) {
      throw std::runtime_error("fermion reshape: too many output axes");
    }
    std::size_t dim = a.parity[old_ax].size();
    parity_vector fused = a.parity[old_ax];
    ++old_ax;
    while (dim < sh[new_ax] && old_ax < a.parity.size()) {
      fused = fuse(fused, a.parity[old_ax]);
      dim *= a.parity[old_ax].size();
      ++old_ax;
    }
    if (dim != sh[new_ax]) {
      throw std::runtime_error(
          "fermion reshape: only adjacent leg fusion is supported");
    }
    ret.parity.push_back(fused);
  }
  if (old_ax != a.parity.size()) {
    throw std::runtime_error("fermion reshape: unused input axes");
  }
  return ret;
}

//! Largest element magnitude (collective; mirrors mptensor::max_abs).
template <class tensor>
double max_abs(const ftensor<tensor>& a) {
  return mptensor::max_abs(a.t);
}

/*!
 * @brief Collective check that a graded tensor is parity even.
 *
 * Every rank reduces the largest parity-odd magnitude (parity_violation(),
 * which is process-local) and the largest magnitude of @p a over all ranks,
 * so every rank throws or none does. A check that threw on the owning rank
 * only would leave the other ranks waiting in their next collective call.
 *
 * @param[in] a Tensor to check. Collective over its communicator.
 * @param[in] context Prefix of the error message.
 * @throw std::runtime_error On every rank, if the parity-odd sector carries
 *        weight above 1e-10 relative to max(1, largest element).
 */
template <class tensor>
void require_even_parity(const ftensor<tensor>& a, const char* context) {
  std::vector<double> v{parity_violation(a)};
  tenes::allreduce_max(v, a.t.get_comm());
  const double threshold = 1.0e-10 * std::max(1.0, max_abs(a));
  if (v[0] > threshold) {
    std::stringstream ss;
    ss << context
       << " is not parity even; odd-parity elements up to max_abs=" << v[0]
       << " (threshold " << threshold << ")";
    throw std::runtime_error(ss.str());
  }
}

/*!
 * @brief Guard of the simple update: refuse weight in the parity-odd
 *        sector, and clip the round-off that is left there.
 *
 * The graded update must keep the state parity even; an odd-sector element
 * above 1e-10 relative to the largest element indicates a sign-bookkeeping
 * bug upstream. Below that threshold the odd-sector elements are numerical
 * residue and are set to zero. Collective (see require_even_parity()).
 *
 * @param[in,out] a Updated site tensor.
 * @throw std::runtime_error On every rank, if the odd sector carries weight
 *        above the threshold.
 */
template <class tensor>
void enforce_even_parity(ftensor<tensor>& a) {
  require_even_parity(a, "fermion Simple_update_bond: the updated tensor");
  mptensor::Index index;
  index.resize(a.t.shape().size());
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    a.t.global_index_fast(n, index);
    if (count_odd(a.parity, index) % 2 == 1) {
      a.t[n] = typename tensor::value_type{};
    }
  }
}

//! Singular values only, of the plain matricization (no parity blocking;
//! mirrors the matrix overload of mptensor::svd).
template <class tensor>
int svd(const ftensor<tensor>& a, std::vector<double>& s) {
  return mptensor::svd(a.t, s);
}

//! Permutation matrix P with P[i][j] = 1 iff perm[i] == j, used to apply
//! the even-first sort to a matricization.
template <class tensor>
tensor make_perm_matrix(const std::vector<std::size_t>& perm) {
  tensor ret(mptensor::Shape(perm.size(), perm.size()));
  mptensor::Index idx;
  idx.resize(2);
  for (std::size_t n = 0; n < ret.local_size(); ++n) {
    ret.global_index_fast(n, idx);
    ret[n] = (perm[idx[0]] == idx[1]) ? 1.0 : 0.0;
  }
  return ret;
}

namespace detail {

/*! @brief Element budget below which a block is factorized on every rank.
 *
 * The graded decomposition splits every matrix by parity, so where the
 * bosonic path factorizes one D x D matrix this one factorizes two of
 * about half that side: a fermionic full update at D = 3, d = 4 hands
 * LAPACK a 2x2 and a 1x1. Spreading a matrix that small over a process
 * grid buys nothing - mptensor's block-cyclic block is 16x16, so it lives
 * on one process either way - and it costs pdgesvd's cross-rank singular
 * value check, which is what info_is_grid_heterogeneity() reports.
 *
 * So blocks up to this many elements are gathered onto every rank and
 * factorized there instead. The default, 4096, is where the largest blocks
 * a D = 8, d = 4 full update produces land: the QR of a site tensor is
 * D^3 x D*d, giving blocks of about 256 x 16, and the truncated SVD of
 * Theta is D*d^2 square, giving blocks of about 64 x 64. The replicated
 * copy is then at most 32 KB real, 64 KB complex.
 *
 * TENES_FERMION_LOCAL_DECOMP_MAX overrides it; 0 disables the local path
 * and sends every block to the grid, which is the way to tell whether this
 * routing is responsible for something without rebuilding.
 */
inline std::size_t local_decomposition_max_elements() {
  const char* raw = std::getenv("TENES_FERMION_LOCAL_DECOMP_MAX");
  if (raw == nullptr) {
    return 4096;
  }
  const long value = std::atol(raw);
  return value > 0 ? static_cast<std::size_t>(value) : 0;
}

//! True when a rows x cols block is within local_decomposition_max_elements().
//! Written as a division so a product cannot overflow on the way to the
//! comparison.
inline bool fits_local_decomposition(std::size_t rows, std::size_t cols) {
  const std::size_t budget = local_decomposition_max_elements();
  if (budget == 0) {
    return false;
  }
  if (rows == 0 || cols == 0) {
    return true;
  }
  return rows <= budget / cols;
}

//! The non-distributed tensor type over the same scalars as @p tensor.
template <class tensor>
using local_tensor_type =
    mptensor::Tensor<mptensor::lapack::Matrix<typename tensor::value_type>>;

/*!
 * @brief Copy of a rank-2 tensor held whole by every rank.
 *
 * mptensor::Tensor::gather() is meant for this and is what its own eig()
 * uses, but it hands the distributed communicator to the non-distributed
 * tensor's constructor, whose comm_type is int - so it does not compile
 * for a scalapack tensor. This does the same two steps (flatten, which is
 * an MPI_Allreduce and therefore gives every rank the same vector, then
 * reshape) without that conversion.
 *
 * @param[in,out] a The tensor; flatten() reorders its storage.
 */
template <class tensor>
local_tensor_type<tensor> gather_matrix(tensor& a) {
  const mptensor::Shape shape = a.shape();
  const std::vector<typename tensor::value_type> flat = a.flatten();
  local_tensor_type<tensor> local(shape);
  mptensor::Index idx;
  idx.resize(2);
  for (std::size_t n = 0; n < local.local_size(); ++n) {
    local.global_index_fast(n, idx);
    // The layout flatten() promises: row index fastest.
    local[n] = flat[idx[0] + idx[1] * shape[0]];
  }
  return local;
}

/*!
 * @brief SVD of one parity block, on one rank when it is small enough.
 *
 * gather_matrix() reduces over the ranks, so every rank starts the local
 * branch from the same matrix, and one routine given one matrix cannot
 * return different answers on different ranks - which is the failure this
 * branch exists to remove. The way back is mptensor's "all processes have
 * the same data" constructor, the one its own eig() uses.
 *
 * The local branch hands LAPACK a copy of @p block holding the same values,
 * so a diagnostic taken from @p block still describes what LAPACK saw.
 *
 * In a serial build the two branches call the same LAPACK routine on the
 * same values, so which one runs is not observable in the result.
 *
 * @param[in,out] block The block; gather_matrix() reorders its storage.
 * @param[out] u,s,vt The factors, in @p block's own tensor type.
 * @return The info of whichever routine ran.
 */
template <class tensor>
int block_svd(tensor& block, tensor& u, std::vector<double>& s, tensor& vt) {
  const mptensor::Shape shape = block.shape();
  if (!fits_local_decomposition(shape[0], shape[1])) {
    return mptensor::svd(block, mptensor::Axes(0), mptensor::Axes(1), u, s, vt);
  }
  local_tensor_type<tensor> local = gather_matrix(block);
  local_tensor_type<tensor> local_u, local_vt;
  const int info = mptensor::svd(local, mptensor::Axes(0), mptensor::Axes(1),
                                 local_u, s, local_vt);
  u = tensor(block.get_comm(), local_u);
  vt = tensor(block.get_comm(), local_vt);
  return info;
}

//! QR of one parity block, on one rank when it is small enough. See
//! block_svd() for why.
template <class tensor>
int block_qr(tensor& block, tensor& q, tensor& r) {
  const mptensor::Shape shape = block.shape();
  if (!fits_local_decomposition(shape[0], shape[1])) {
    return mptensor::qr(block, mptensor::Axes(0), mptensor::Axes(1), q, r);
  }
  local_tensor_type<tensor> local = gather_matrix(block);
  local_tensor_type<tensor> local_q, local_r;
  const int info = mptensor::qr(local, mptensor::Axes(0), mptensor::Axes(1),
                                local_q, local_r);
  q = tensor(block.get_comm(), local_q);
  r = tensor(block.get_comm(), local_r);
  return info;
}

/*! @brief True for ScaLAPACK's "the ranks disagreed" code.
 *
 * An MPI build decomposes through pdgesvd, not dgesvd, and pdgesvd.f
 * documents one positive info specially:
 *
 *   "> 0: if DBDSQR did not converge. If INFO = MIN(M,N) + 1, then PDGESVD
 *    has detected heterogeneity by finding that eigenvalues were not
 *    identical across the process grid. In this case, the accuracy of the
 *    results from PDGESVD cannot be guaranteed."
 *
 * The graded decomposition splits every matrix by parity, so its blocks
 * are small - 2x2 and 1x1 are routine - and their singular values often
 * nearly degenerate, which is where ranks stop agreeing. A run met exactly
 * this: info=3 on a 2x2 under four MPI processes, gone at one, with the
 * two singular values agreeing to 2e-5. A serial build never reaches
 * pdgesvd, so there the same value carries no such meaning.
 */
inline bool info_is_grid_heterogeneity(std::size_t rows, std::size_t cols,
                                       int info, bool from_svd) {
#ifdef _NO_MPI
  static_cast<void>(rows);
  static_cast<void>(cols);
  static_cast<void>(info);
  static_cast<void>(from_svd);
  return false;
#else
  // pdgeqrf documents INFO as 0 or -i only, so the code belongs to the SVD
  // alone; reading it into a QR would be the same mislabelling in reverse.
  return from_svd && info > 0 &&
         static_cast<std::size_t>(info) == std::min(rows, cols) + 1;
#endif
}

/*! @brief True when info cannot have come from a conforming library.
 *
 * A positive info from ?gesvd counts the superdiagonals of an intermediate
 * bidiagonal form that failed to converge; that form is min(rows, cols)
 * square, so it has min(rows, cols) - 1 of them. The bound used here is
 * the looser min(rows, cols), so only a value no implementation can
 * justify is called out - and ScaLAPACK's MIN(M,N)+1 is excluded first,
 * since it is a documented code rather than a count.
 */
inline bool info_out_of_range(std::size_t rows, std::size_t cols, int info,
                              bool from_svd) {
  if (info_is_grid_heterogeneity(rows, cols, info, from_svd)) {
    return false;
  }
  return info > 0 && static_cast<std::size_t>(info) > std::min(rows, cols);
}

//! "even sector 9x9 info=3 (LAPACK: did not converge)", or "odd sector
//! empty" for a sector the decomposition never had a block for.
inline std::string sector_phrase(const char* name, std::size_t rows,
                                 std::size_t cols, int info, bool from_svd) {
  std::stringstream ss;
  ss << name << " sector ";
  if (rows == 0 || cols == 0) {
    ss << "empty";
    return ss.str();
  }
  ss << rows << "x" << cols << " info=" << info;
  if (info_is_grid_heterogeneity(rows, cols, info, from_svd)) {
    ss << " (ScaLAPACK: MIN(M,N)+1 - pdgesvd detected heterogeneity, i.e. "
          "the process grid did not agree on this block's singular values; "
          "not a bad matrix)";
  } else if (info_out_of_range(rows, cols, info, from_svd)) {
    ss << " (LAPACK: out of range, at most " << std::min(rows, cols)
       << " for a " << rows << "x" << cols
       << " block - suspect the LAPACK/BLAS build, not the state)";
  } else if (info > 0) {
    ss << " (LAPACK: did not converge)";
  } else if (info < 0) {
    ss << " (LAPACK: illegal argument " << -info << ")";
  }
  return ss.str();
}

/*! @brief Write out a failed block, if it is small enough to be read.
 *
 * Small means "few enough numbers to retype into a standalone LAPACK
 * call": the point of the dump is to turn a failure that only reproduces
 * after two five-thousand-step simple updates into four numbers anyone can
 * hand to ?gesvd. @p m is the block as it was handed to LAPACK, not the
 * matrix it came from. Printed at 17 significant digits, or it would not
 * reproduce. Rank-local, like the rest of the scan: an element another
 * process holds is shown as "?".
 */
template <class tensor>
std::string format_block(const tensor& m, const char* name) {
  const std::size_t max_elements = 16;
  if (m.shape().size() != 2) {
    return std::string();
  }
  const std::size_t rows = m.shape()[0];
  const std::size_t cols = m.shape()[1];
  if (rows == 0 || cols == 0 || rows * cols > max_elements) {
    return std::string();
  }
  std::stringstream ss;
  ss << std::setprecision(17);
  ss << "; " << name << " block (" << rows << "x" << cols << ", row-major) [";
  for (std::size_t i = 0; i < rows; ++i) {
    if (i > 0) {
      ss << "; ";
    }
    for (std::size_t j = 0; j < cols; ++j) {
      if (j > 0) {
        ss << ", ";
      }
      typename tensor::value_type value;
      if (m.get_value(mptensor::Index(i, j), value)) {
        ss << value;
      } else {
        ss << "?";
      }
    }
  }
  ss << "]";
  return ss.str();
}

}  // namespace detail

/*!
 * @brief What a graded decomposition saw, for the error message.
 *
 * qr() and svd() factorize the even and the odd diagonal block separately
 * and return the FIRST nonzero LAPACK info of the two, which leaves the
 * caller unable to say which sector failed - and the odd sector is the
 * fermion-specific half, the one that goes empty or degenerate in the
 * known failure modes. Passing one of these to the decomposition records
 * both infos together with the block shapes; on failure it also records
 * whether the input was finite at all.
 */
struct decomposition_diagnostics {
  int info_even = 0;           //!< LAPACK info of the even block.
  int info_odd = 0;            //!< LAPACK info of the odd block.
  std::size_t row_even = 0;    //!< Rows of even parity.
  std::size_t col_even = 0;    //!< Columns of even parity.
  std::size_t row_odd = 0;     //!< Rows of odd parity.
  std::size_t col_odd = 0;     //!< Columns of odd parity.
  bool scanned = false;        //!< True once note_failed_block() has run.
  bool has_nonfinite = false;  //!< NaN or Inf found by that scan.
  double max_abs = 0.0;        //!< Largest finite |element| found by it.
  std::string block_dump;      //!< Elements of a small failed block.
  //! True when the infos came from svd(), whose ScaLAPACK counterpart has
  //! a documented code that qr()'s does not.
  bool from_svd = false;

  //! True iff either block factorization reported a failure.
  bool failed() const { return info_even != 0 || info_odd != 0; }

  //! True when a reported info cannot have come from a conforming LAPACK,
  //! i.e. when the build rather than the state is the thing to look at.
  bool suspect_library() const {
    return detail::info_out_of_range(row_even, col_even, info_even, from_svd) ||
           detail::info_out_of_range(row_odd, col_odd, info_odd, from_svd);
  }

  //! True when ScaLAPACK reported that the MPI ranks disagreed; then it is
  //! the process grid, not the state and not the library.
  bool grid_heterogeneity() const {
    return detail::info_is_grid_heterogeneity(row_even, col_even, info_even,
                                              from_svd) ||
           detail::info_is_grid_heterogeneity(row_odd, col_odd, info_odd,
                                              from_svd);
  }

  /*! @brief Record the block a failed decomposition was given.
   *
   * Called with the very tensor that went into LAPACK - qr() and svd()
   * decompose a slice() copy of each diagonal block, not the parity-sorted
   * matrix, and a block that was damaged on its way in would otherwise be
   * reported through a healthy source. The scan is an extra pass over the
   * block, so it is skipped while the factorization succeeds; a non-finite
   * element says what the sector shapes cannot, namely that the state had
   * already diverged. Process-local, like parity_violation().
   *
   * @param[in] m The block as handed to LAPACK.
   * @param[in] name "even" or "odd", for the message.
   */
  template <class tensor>
  void note_failed_block(const tensor& m, const char* name) {
    if (!failed()) {
      return;
    }
    scanned = true;
    for (std::size_t n = 0; n < m.local_size(); ++n) {
      const typename tensor::value_type v = m[n];
      if (!detail::is_finite_value(v)) {
        has_nonfinite = true;
        continue;
      }
      max_abs = std::max(max_abs, std::abs(v));
    }
    block_dump += detail::format_block(m, name);
  }

  //! One line naming each sector, its shape, and its LAPACK info, plus
  //! what the input scan found once it has run.
  std::string describe() const {
    std::stringstream ss;
    ss << detail::sector_phrase("even", row_even, col_even, info_even, from_svd)
       << ", "
       << detail::sector_phrase("odd", row_odd, col_odd, info_odd, from_svd);
    if (scanned) {
      ss << "; input max|.|=" << max_abs << ", "
         << (has_nonfinite ? "non-finite elements present"
                           : "all elements finite")
         << " (rank-local scan)";
    }
    ss << block_dump;
    return ss.str();
  }
};

/*!
 * @brief Grading-preserving QR decomposition.
 *
 * Matricizes @p a with @p rows fused as rows and @p cols as columns, sorts
 * both fused legs even-first, and QR-factorizes the even and odd diagonal
 * blocks separately (the input must be parity even, so the off-diagonal
 * blocks vanish — checked in debug builds by
 * detail::validate_block_diagonal()). The new internal leg concatenates
 * the two blocks' factors and gets an even-first ledger whose even count
 * is min(number of even rows, number of even columns).
 *
 * @param[in] a Parity-even tensor to decompose.
 * @param[in] rows Legs forming the rows.
 * @param[in] cols Legs forming the columns.
 * @param[out] q Orthogonal factor, legs (rows..., internal).
 * @param[out] r Triangular-block factor, legs (internal, cols...).
 * @param[out] diag Optional per-sector record; see
 *        decomposition_diagnostics.
 * @return The LAPACK-style info of the block factorizations (0 on
 *         success), the even block's first.
 */
template <class tensor>
int qr(const ftensor<tensor>& a, const mptensor::Axes& rows,
       const mptensor::Axes& cols, ftensor<tensor>& q, ftensor<tensor>& r,
       decomposition_diagnostics* diag = nullptr) {
  ftensor<tensor> a_ordered = transpose(a, rows + cols);
  mptensor::Shape row_shape = detail::shape_from_axes(a.shape(), rows);
  mptensor::Shape col_shape = detail::shape_from_axes(a.shape(), cols);
  const std::size_t drow = detail::product_shape(row_shape);
  const std::size_t dcol = detail::product_shape(col_shape);
  typename tensor::comm_type comm = a.t.get_comm();
  tensor mat = mptensor::reshape(a_ordered.t, mptensor::Shape(drow, dcol));

  parity_vector row_parity = detail::fuse_axes(a.parity, rows);
  parity_vector col_parity = detail::fuse_axes(a.parity, cols);
  std::vector<std::size_t> row_perm = parity_sort_perm(row_parity);
  std::vector<std::size_t> col_perm = parity_sort_perm(col_parity);
  tensor prow = make_perm_matrix<tensor>(row_perm);
  tensor pcol = make_perm_matrix<tensor>(col_perm);
  tensor sorted =
      mptensor::tensordot(prow, mat, mptensor::Axes(1), mptensor::Axes(0));
  sorted =
      mptensor::tensordot(sorted, pcol, mptensor::Axes(1), mptensor::Axes(1));

  const std::size_t row_even = detail::count_even(row_parity);
  const std::size_t col_even = detail::count_even(col_parity);
  detail::validate_block_diagonal(sorted, row_even, col_even, "fermion qr");
  const std::size_t row_odd = drow - row_even;
  const std::size_t col_odd = dcol - col_even;
  const std::size_t size_even = std::min(row_even, col_even);
  const std::size_t size_odd = std::min(row_odd, col_odd);
  const std::size_t size = size_even + size_odd;

  tensor q_sorted(comm, mptensor::Shape(drow, size));
  tensor r_sorted(comm, mptensor::Shape(size, dcol));
  if (diag != nullptr) {
    // Fully overwritten, so one instance can be reused across calls.
    *diag = decomposition_diagnostics{};
    diag->row_even = row_even;
    diag->col_even = col_even;
    diag->row_odd = row_odd;
    diag->col_odd = col_odd;
  }
  int info_even_block = 0;
  int info_odd_block = 0;
  if (size_even > 0) {
    tensor even_block = mptensor::slice(sorted, mptensor::Index(0, 0),
                                        mptensor::Index(row_even, col_even));
    tensor qe, re;
    info_even_block = detail::block_qr(even_block, qe, re);
    // Recorded here, where the block LAPACK actually saw is still in scope.
    if (diag != nullptr) {
      diag->info_even = info_even_block;
      if (info_even_block != 0) {
        diag->note_failed_block(even_block, "even");
      }
    }
    q_sorted.set_slice(qe, mptensor::Index(0, 0),
                       mptensor::Index(row_even, size_even));
    r_sorted.set_slice(re, mptensor::Index(0, 0),
                       mptensor::Index(size_even, col_even));
  }
  if (size_odd > 0) {
    tensor odd_block =
        mptensor::slice(sorted, mptensor::Index(row_even, col_even),
                        mptensor::Index(drow, dcol));
    tensor qo, ro;
    info_odd_block = detail::block_qr(odd_block, qo, ro);
    if (diag != nullptr) {
      diag->info_odd = info_odd_block;
      if (info_odd_block != 0) {
        diag->note_failed_block(odd_block, "odd");
      }
    }
    q_sorted.set_slice(qo, mptensor::Index(row_even, size_even),
                       mptensor::Index(drow, size));
    r_sorted.set_slice(ro, mptensor::Index(size_even, col_even),
                       mptensor::Index(size, dcol));
  }
  // The even block's failure is reported first, as it always has been.
  int info = info_even_block != 0 ? info_even_block : info_odd_block;

  tensor q_mat =
      mptensor::tensordot(prow, q_sorted, mptensor::Axes(0), mptensor::Axes(0));
  tensor r_mat =
      mptensor::tensordot(r_sorted, pcol, mptensor::Axes(1), mptensor::Axes(0));

  mptensor::Shape q_shape = row_shape;
  q_shape.push(size);
  q.t = mptensor::reshape(q_mat, q_shape);
  mptensor::Shape r_shape;
  r_shape.push(size);
  for (std::size_t i = 0; i < col_shape.size(); ++i) {
    r_shape.push(col_shape[i]);
  }
  r.t = mptensor::reshape(r_mat, r_shape);

  parity_vector internal(size, false);
  for (std::size_t i = size_even; i < size; ++i) {
    internal[i] = true;
  }
  q.parity.clear();
  for (std::size_t i = 0; i < rows.size(); ++i) {
    q.parity.push_back(a.parity[rows[i]]);
  }
  q.parity.push_back(internal);
  r.parity.clear();
  r.parity.push_back(internal);
  for (std::size_t i = 0; i < cols.size(); ++i) {
    r.parity.push_back(a.parity[cols[i]]);
  }
  return info;
}

/*!
 * @brief Grading-preserving singular value decomposition.
 *
 * Same block structure as qr(): the even-first sorted matricization is
 * SVD-ed per parity block, and the internal leg carries an even-first
 * ledger. The singular values are returned sector-concatenated — all even
 * ones first, then all odd ones — NOT globally sorted; use svd_trunc() to
 * truncate by magnitude.
 *
 * @param[in] a Parity-even tensor to decompose.
 * @param[in] rows Legs forming the rows.
 * @param[in] cols Legs forming the columns.
 * @param[out] u Left singular vectors, legs (rows..., internal).
 * @param[out] s Singular values, even sector then odd sector.
 * @param[out] vt Right singular vectors, legs (internal, cols...).
 * @param[out] diag Optional per-sector record; see
 *        decomposition_diagnostics.
 * @return The LAPACK-style info of the block factorizations (0 on
 *         success).
 */
template <class tensor>
int svd(const ftensor<tensor>& a, const mptensor::Axes& rows,
        const mptensor::Axes& cols, ftensor<tensor>& u, std::vector<double>& s,
        ftensor<tensor>& vt, decomposition_diagnostics* diag = nullptr) {
  ftensor<tensor> a_ordered = transpose(a, rows + cols);
  mptensor::Shape row_shape = detail::shape_from_axes(a.shape(), rows);
  mptensor::Shape col_shape = detail::shape_from_axes(a.shape(), cols);
  const std::size_t drow = detail::product_shape(row_shape);
  const std::size_t dcol = detail::product_shape(col_shape);
  typename tensor::comm_type comm = a.t.get_comm();
  tensor mat = mptensor::reshape(a_ordered.t, mptensor::Shape(drow, dcol));

  parity_vector row_parity = detail::fuse_axes(a.parity, rows);
  parity_vector col_parity = detail::fuse_axes(a.parity, cols);
  std::vector<std::size_t> row_perm = parity_sort_perm(row_parity);
  std::vector<std::size_t> col_perm = parity_sort_perm(col_parity);
  tensor prow = make_perm_matrix<tensor>(row_perm);
  tensor pcol = make_perm_matrix<tensor>(col_perm);
  tensor sorted =
      mptensor::tensordot(prow, mat, mptensor::Axes(1), mptensor::Axes(0));
  sorted =
      mptensor::tensordot(sorted, pcol, mptensor::Axes(1), mptensor::Axes(1));

  const std::size_t row_even = detail::count_even(row_parity);
  const std::size_t col_even = detail::count_even(col_parity);
  detail::validate_block_diagonal(sorted, row_even, col_even, "fermion svd");
  const std::size_t row_odd = drow - row_even;
  const std::size_t col_odd = dcol - col_even;
  const std::size_t size_even = std::min(row_even, col_even);
  const std::size_t size_odd = std::min(row_odd, col_odd);
  const std::size_t size = size_even + size_odd;

  tensor u_sorted(comm, mptensor::Shape(drow, size));
  tensor vt_sorted(comm, mptensor::Shape(size, dcol));
  std::vector<double> s_even, s_odd;
  if (diag != nullptr) {
    // Fully overwritten, so one instance can be reused across calls.
    *diag = decomposition_diagnostics{};
    diag->row_even = row_even;
    diag->col_even = col_even;
    diag->row_odd = row_odd;
    diag->col_odd = col_odd;
    diag->from_svd = true;
  }
  int info_even_block = 0;
  int info_odd_block = 0;
  if (size_even > 0) {
    tensor even_block = mptensor::slice(sorted, mptensor::Index(0, 0),
                                        mptensor::Index(row_even, col_even));
    tensor ue, vte;
    info_even_block = detail::block_svd(even_block, ue, s_even, vte);
    // Recorded here, where the block LAPACK actually saw is still in scope.
    if (diag != nullptr) {
      diag->info_even = info_even_block;
      if (info_even_block != 0) {
        diag->note_failed_block(even_block, "even");
      }
    }
    u_sorted.set_slice(ue, mptensor::Index(0, 0),
                       mptensor::Index(row_even, size_even));
    vt_sorted.set_slice(vte, mptensor::Index(0, 0),
                        mptensor::Index(size_even, col_even));
  }
  if (size_odd > 0) {
    tensor odd_block =
        mptensor::slice(sorted, mptensor::Index(row_even, col_even),
                        mptensor::Index(drow, dcol));
    tensor uo, vto;
    info_odd_block = detail::block_svd(odd_block, uo, s_odd, vto);
    if (diag != nullptr) {
      diag->info_odd = info_odd_block;
      if (info_odd_block != 0) {
        diag->note_failed_block(odd_block, "odd");
      }
    }
    u_sorted.set_slice(uo, mptensor::Index(row_even, size_even),
                       mptensor::Index(drow, size));
    vt_sorted.set_slice(vto, mptensor::Index(size_even, col_even),
                        mptensor::Index(size, dcol));
  }
  // The even block's failure is reported first, as it always has been.
  int info = info_even_block != 0 ? info_even_block : info_odd_block;

  tensor u_mat =
      mptensor::tensordot(prow, u_sorted, mptensor::Axes(0), mptensor::Axes(0));
  tensor vt_mat = mptensor::tensordot(vt_sorted, pcol, mptensor::Axes(1),
                                      mptensor::Axes(0));

  mptensor::Shape u_shape = row_shape;
  u_shape.push(size);
  u.t = mptensor::reshape(u_mat, u_shape);
  mptensor::Shape vt_shape;
  vt_shape.push(size);
  for (std::size_t i = 0; i < col_shape.size(); ++i) {
    vt_shape.push(col_shape[i]);
  }
  vt.t = mptensor::reshape(vt_mat, vt_shape);

  s = s_even;
  s.insert(s.end(), s_odd.begin(), s_odd.end());
  parity_vector internal(size, false);
  for (std::size_t i = size_even; i < size; ++i) {
    internal[i] = true;
  }
  u.parity.clear();
  for (std::size_t i = 0; i < rows.size(); ++i) {
    u.parity.push_back(a.parity[rows[i]]);
  }
  u.parity.push_back(internal);
  vt.parity.clear();
  vt.parity.push_back(internal);
  for (std::size_t i = 0; i < cols.size(); ++i) {
    vt.parity.push_back(a.parity[cols[i]]);
  }
  return info;
}

/*!
 * @brief Truncated grading-preserving SVD.
 *
 * Runs the full svd(), keeps the @p dc largest singular values across BOTH
 * parity sectors (ties resolved even-first, then by original position, so
 * the truncation is deterministic), and re-sorts the kept values
 * even-first for the internal leg's ledger. The per-sector spectra can be
 * dumped for diagnosis via TENES_FERMION_SVD_LOG_LIMIT (see
 * detail::fermion_svd_log_limit()).
 *
 * @param[in] a Parity-even tensor to decompose.
 * @param[in] rows Legs forming the rows.
 * @param[in] cols Legs forming the columns.
 * @param[out] u Left singular vectors, legs (rows..., internal).
 * @param[out] s Kept singular values, even sector then odd sector.
 * @param[out] vt Right singular vectors, legs (internal, cols...).
 * @param[in] dc Maximum number of singular values to keep.
 * @param[out] diag Optional per-sector record; see
 *        decomposition_diagnostics.
 * @return The LAPACK-style info of the underlying svd() (0 on success).
 */
template <class tensor>
int svd_trunc(const ftensor<tensor>& a, const mptensor::Axes& rows,
              const mptensor::Axes& cols, ftensor<tensor>& u,
              std::vector<double>& s, ftensor<tensor>& vt, int dc,
              decomposition_diagnostics* diag = nullptr) {
  ftensor<tensor> full_u, full_vt;
  std::vector<double> full_s;
  int info = svd(a, rows, cols, full_u, full_s, full_vt, diag);
  const int nkeep = std::min<int>(dc, full_s.size());
  std::vector<std::size_t> order(full_s.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    order[i] = i;
  }
  const parity_vector& full_internal = full_u.parity.back();
  static int svd_log_count = 0;
  const int svd_log_limit = detail::fermion_svd_log_limit();
  if (svd_log_count < svd_log_limit) {
    std::cerr << "TENES_FERMION_SVD call=" << svd_log_count << " dc=" << dc
              << " full_size=" << full_s.size() << " even=[";
    bool first = true;
    for (std::size_t i = 0; i < full_s.size(); ++i) {
      if (!full_internal[i]) {
        std::cerr << (first ? "" : ",") << full_s[i];
        first = false;
      }
    }
    std::cerr << "] odd=[";
    first = true;
    for (std::size_t i = 0; i < full_s.size(); ++i) {
      if (full_internal[i]) {
        std::cerr << (first ? "" : ",") << full_s[i];
        first = false;
      }
    }
    std::cerr << "]\n";
  }
  ++svd_log_count;
  std::stable_sort(order.begin(), order.end(),
                   [&](std::size_t lhs, std::size_t rhs) {
                     if (full_s[lhs] != full_s[rhs]) {
                       return full_s[lhs] > full_s[rhs];
                     }
                     if (full_internal[lhs] != full_internal[rhs]) {
                       return !full_internal[lhs];
                     }
                     return lhs < rhs;
                   });
  order.resize(nkeep);
  std::stable_sort(order.begin(), order.end(),
                   [&](std::size_t lhs, std::size_t rhs) {
                     if (full_internal[lhs] != full_internal[rhs]) {
                       return !full_internal[lhs];
                     }
                     return full_s[lhs] > full_s[rhs];
                   });

  mptensor::Shape row_shape = detail::shape_from_axes(a.shape(), rows);
  mptensor::Shape col_shape = detail::shape_from_axes(a.shape(), cols);
  const std::size_t drow = detail::product_shape(row_shape);
  const std::size_t dcol = detail::product_shape(col_shape);
  tensor full_u_mat =
      mptensor::reshape(full_u.t, mptensor::Shape(drow, full_s.size()));
  tensor full_vt_mat =
      mptensor::reshape(full_vt.t, mptensor::Shape(full_s.size(), dcol));
  tensor selector(a.t.get_comm(), mptensor::Shape(full_s.size(), nkeep));
  mptensor::Index sel_idx;
  sel_idx.resize(2);
  for (std::size_t n = 0; n < selector.local_size(); ++n) {
    selector.global_index_fast(n, sel_idx);
    selector[n] = (order[sel_idx[1]] == sel_idx[0]) ? 1.0 : 0.0;
  }
  tensor u_mat = mptensor::tensordot(full_u_mat, selector, mptensor::Axes(1),
                                     mptensor::Axes(0));
  tensor vt_mat = mptensor::tensordot(selector, full_vt_mat, mptensor::Axes(0),
                                      mptensor::Axes(0));

  mptensor::Shape u_shape = row_shape;
  u_shape.push(nkeep);
  u.t = mptensor::reshape(u_mat, u_shape);
  mptensor::Shape vt_shape;
  vt_shape.push(nkeep);
  for (std::size_t i = 0; i < col_shape.size(); ++i) {
    vt_shape.push(col_shape[i]);
  }
  vt.t = mptensor::reshape(vt_mat, vt_shape);

  s.resize(nkeep);
  parity_vector internal(nkeep, false);
  for (int i = 0; i < nkeep; ++i) {
    s[i] = full_s[order[i]];
    internal[i] = full_internal[order[i]];
  }
  u.parity.clear();
  for (std::size_t i = 0; i < rows.size(); ++i) {
    u.parity.push_back(a.parity[rows[i]]);
  }
  u.parity.push_back(internal);
  vt.parity.clear();
  vt.parity.push_back(internal);
  for (std::size_t i = 0; i < cols.size(); ++i) {
    vt.parity.push_back(a.parity[cols[i]]);
  }
  return info;
}

}  // namespace fermion
}  // namespace tenes

#endif  // TENES_SRC_FERMION_FOPS_HPP_
