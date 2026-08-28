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
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

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
 * even. Checked in debug builds only; like parity_violation() this
 * inspects the process-local slice.
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
 * @return The LAPACK-style info of the block factorizations (0 on
 *         success).
 */
template <class tensor>
int qr(const ftensor<tensor>& a, const mptensor::Axes& rows,
       const mptensor::Axes& cols, ftensor<tensor>& q, ftensor<tensor>& r) {
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
  int info = 0;
  if (size_even > 0) {
    tensor even_block = mptensor::slice(sorted, mptensor::Index(0, 0),
                                        mptensor::Index(row_even, col_even));
    tensor qe, re;
    info =
        mptensor::qr(even_block, mptensor::Axes(0), mptensor::Axes(1), qe, re);
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
    int info_odd =
        mptensor::qr(odd_block, mptensor::Axes(0), mptensor::Axes(1), qo, ro);
    if (info == 0) {
      info = info_odd;
    }
    q_sorted.set_slice(qo, mptensor::Index(row_even, size_even),
                       mptensor::Index(drow, size));
    r_sorted.set_slice(ro, mptensor::Index(size_even, col_even),
                       mptensor::Index(size, dcol));
  }

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
 * @return The LAPACK-style info of the block factorizations (0 on
 *         success).
 */
template <class tensor>
int svd(const ftensor<tensor>& a, const mptensor::Axes& rows,
        const mptensor::Axes& cols, ftensor<tensor>& u, std::vector<double>& s,
        ftensor<tensor>& vt) {
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
  int info = 0;
  if (size_even > 0) {
    tensor even_block = mptensor::slice(sorted, mptensor::Index(0, 0),
                                        mptensor::Index(row_even, col_even));
    tensor ue, vte;
    info = mptensor::svd(even_block, mptensor::Axes(0), mptensor::Axes(1), ue,
                         s_even, vte);
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
    int info_odd = mptensor::svd(odd_block, mptensor::Axes(0),
                                 mptensor::Axes(1), uo, s_odd, vto);
    if (info == 0) {
      info = info_odd;
    }
    u_sorted.set_slice(uo, mptensor::Index(row_even, size_even),
                       mptensor::Index(drow, size));
    vt_sorted.set_slice(vto, mptensor::Index(size_even, col_even),
                        mptensor::Index(size, dcol));
  }

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
 * @return The LAPACK-style info of the underlying svd() (0 on success).
 */
template <class tensor>
int svd_trunc(const ftensor<tensor>& a, const mptensor::Axes& rows,
              const mptensor::Axes& cols, ftensor<tensor>& u,
              std::vector<double>& s, ftensor<tensor>& vt, int dc) {
  ftensor<tensor> full_u, full_vt;
  std::vector<double> full_s;
  int info = svd(a, rows, cols, full_u, full_s, full_vt);
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
