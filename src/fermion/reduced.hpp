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
 *  @brief Double-layer (bra x ket) reduced tensors for fermionic
 *         measurement.
 *
 *  The fermionic CTM measurement strategy is to fold all fermionic signs
 *  into reduced tensors built here from open two-layer networks; after the
 *  fusion the network has no crossings left, so the existing single-layer
 *  bosonic CTM machinery contracts it unmodified. In particular the naive
 *  mechanical bookkeeping on CLOSED loops would pick up spurious supertrace
 *  @f$(-1)@f$ factors (the arrow/dual-flag design was tried and discarded
 *  for exactly that), which this open-network construction sidesteps.
 *
 *  build_reduced() / build_reduced_op() produce the one-site reduced
 *  tensor; build_reduced_pair_naive() / build_reduced_pair_direct()
 *  produce the two-site "blob" with the operator already inserted. Leg
 *  orders (each blob leg is a fused (ket, bra) pair of dimension
 *  @f$D^2@f$):
 *
 *  @verbatim
    horizontal pair (A left, B right):     vertical pair (A top, B bottom):

        T_A   T_B                                T_A
         |     |                                  |
    L_A- A --- B -R_B                        L_A- A -R_A
         |     |                                  |
        B_A   B_B                            L_B- B -R_B
                                                  |
    blob legs: (L_A, T_A, B_A,                   B_B
                T_B, R_B, B_B)
                                       blob legs: (L_A, T_A, R_A,
                                                   L_B, R_B, B_B)
    @endverbatim
 */

#ifndef TENES_SRC_FERMION_REDUCED_HPP_
#define TENES_SRC_FERMION_REDUCED_HPP_

#include <stdexcept>
#include <vector>

#include "fops.hpp"
#include "ftensor.hpp"

namespace tenes::fermion {

//! Orientation of a two-site window: horizontal (A left of B) or vertical
//! (A above B).
enum class reduced_pair_direction { horizontal, vertical };

namespace detail {

//! Bit position encoding the ORDERED direction-label pair (x, y),
//! x != y, in kDoubledJointMask (directions 0=left, 1=top, 2=right,
//! 3=bottom).
constexpr int joint_bit(int x, int y) { return x * 3 + (y < x ? y : y - 1); }

/*! @brief Which direction-label pairs generate joint-swap terms.
 *
 * Only the pairs (0,3), (1,0), (2,3), (3,0) are set (0=left, 1=top,
 * 2=right, 3=bottom). This joint-swap pattern is gauge-equivalent to
 * YASTN's fuse_layers() convention, empirically pinned down using the
 * oracle as arbiter — a measured convention, not an analytic derivation
 * (which is why the pair list looks arbitrary).
 */
constexpr unsigned kDoubledJointMask =
    (1u << joint_bit(0, 3)) | (1u << joint_bit(1, 0)) |
    (1u << joint_bit(2, 3)) | (1u << joint_bit(3, 0));

/*! @brief The two swap forms dressing a double-layer fuse.
 *
 * Each generated term is one half of the two-term expansion of "one bra
 * leg overtakes a fused (ket, bra) leg pair": @c cross collects the
 * ket x bra terms and @c bra the bra x bra terms. There is no @c ket
 * member — the ket layer is the resting reference frame, and a ket leg
 * never overtakes another ket leg in this convention. (Ket–ket sign terms
 * do exist, but they come from the Koszul side of the interleave
 * transpose, a separate bookkeeping.)
 */
struct JointSwapForms {
  //! ket x bra terms.
  SwapForm cross;
  //! bra x bra terms.
  SwapForm bra;
};

//! Throw std::runtime_error(context) if any axis appears twice.
inline void validate_unique_joint_axes(const std::vector<int>& axes,
                                       const char* context) {
  for (std::size_t i = 0; i < axes.size(); ++i) {
    for (std::size_t j = i + 1; j < axes.size(); ++j) {
      if (axes[i] == axes[j]) {
        throw std::runtime_error(context);
      }
    }
  }
}

/*!
 * @brief Build the joint-swap forms for fusing a bra and a ket layer.
 *
 * For every ordered pair (x, y) of direction labels enabled in
 * kDoubledJointMask, and every pair of open legs carrying those labels,
 * emits the two-term expansion described at JointSwapForms: a cross term
 * on (ket leg of x, bra leg of y) and a bra term on (bra leg of x, bra leg
 * of y). Axes are positions in whatever outer-product tensor the caller is
 * about to sweep; the leg ids attach the direction labels.
 *
 * @param[in] bra_axes Axis of each open leg in the bra layer.
 * @param[in] ket_axes Axis of the same open legs in the ket layer, in the
 *            same order.
 * @param[in] leg_ids Direction label (0=left, 1=top, 2=right, 3=bottom) of
 *            each open leg; the same label may appear more than once (a
 *            two-site window has e.g. two "top" legs).
 * @throw std::runtime_error On size mismatch, negative or duplicate axes.
 */
inline JointSwapForms joint_swap_forms(const std::vector<int>& bra_axes,
                                       const std::vector<int>& ket_axes,
                                       const std::vector<int>& leg_ids) {
  if (bra_axes.size() != leg_ids.size() || ket_axes.size() != leg_ids.size()) {
    throw std::runtime_error("joint_swap_forms: axis and leg id size mismatch");
  }
  for (std::size_t i = 0; i < leg_ids.size(); ++i) {
    if (bra_axes[i] < 0 || ket_axes[i] < 0) {
      throw std::runtime_error("joint_swap_forms: negative axis");
    }
  }
  validate_unique_joint_axes(bra_axes,
                             "fermion joint swap: duplicate bra axis");
  validate_unique_joint_axes(ket_axes,
                             "fermion joint swap: duplicate ket axis");

  JointSwapForms forms;
  for (int x = 0; x < 4; ++x) {
    for (int y = 0; y < 4; ++y) {
      if (x == y) {
        continue;
      }
      if ((kDoubledJointMask & (1u << joint_bit(x, y))) == 0) {
        continue;
      }
      for (std::size_t ix = 0; ix < leg_ids.size(); ++ix) {
        if (leg_ids[ix] != x) {
          continue;
        }
        for (std::size_t iy = 0; iy < leg_ids.size(); ++iy) {
          if (leg_ids[iy] != y) {
            continue;
          }
          forms.cross.toggle(ket_axes[ix], bra_axes[iy]);
          forms.bra.toggle(bra_axes[ix], bra_axes[iy]);
        }
      }
    }
  }
  return forms;
}

/*!
 * @brief Sign gauge on one fused (ket, bra) leg of a blob.
 *
 * Multiplies -1 onto the fused index values whose ket factor is odd and
 * bra factor even (@p ket_odd_bra_even true) or vice versa (false). This
 * aligns the pair blob's gauge with the single-site doubling convention of
 * doubled_pipeline() (see build_reduced_pair_impl()), so both can meet the
 * same CTM environment.
 *
 * @param[in,out] a Blob tensor (plain: the fused legs carry no ledger).
 * @param[in] leg_parity Single-layer ledger of the leg before doubling.
 * @param[in] ax Fused leg to act on.
 * @param[in] ket_odd_bra_even Which mixed-parity sector gets the -1.
 */
template <class tensor>
void apply_fused_leg_gauge(tensor& a, const parity_vector& leg_parity,
                           std::size_t ax, bool ket_odd_bra_even) {
  std::vector<double> sign(leg_parity.size() * leg_parity.size(), 1.0);
  for (std::size_t bra = 0; bra < leg_parity.size(); ++bra) {
    for (std::size_t ket = 0; ket < leg_parity.size(); ++ket) {
      const bool flip = ket_odd_bra_even
                            ? (leg_parity[ket] && !leg_parity[bra])
                            : (!leg_parity[ket] && leg_parity[bra]);
      sign[ket + leg_parity.size() * bra] = flip ? -1.0 : 1.0;
    }
  }
  a.multiply_vector(sign, ax);
}

/*!
 * @brief Shared single-site double-layer pipeline.
 *
 * Outer product of bra and ket layers, frozen joint-swap dressing,
 * (ket, bra) interleave, column-major fusion. Output legs:
 * ([l lb], [t tb], [r rb], [b bb], s_ket, s_bra) — the four virtual legs
 * fused (ket first and fastest-varying), the two physical legs left open.
 * The two arguments differ only where an operator has been inserted into
 * the ket layer beforehand; for the plain reduced tensor they are equal.
 *
 * @param[in] bra_Tn Rank-5 wrapped Tn to conjugate into the bra layer.
 * @param[in] ket_Tn Rank-5 wrapped Tn forming the ket layer.
 * @return Rank-6 plain tensor as above.
 */
template <class tensor>
tensor doubled_pipeline(const ftensor<tensor>& bra_Tn,
                        const ftensor<tensor>& ket_Tn) {
  const auto forms = joint_swap_forms({0, 1, 2, 3}, {5, 6, 7, 8}, {0, 1, 2, 3});
  ftensor<tensor> bra = conj(bra_Tn);
  // This pre-outer-product swap is the point of this task: it cuts rank-16
  // sweeps. Moving it after the outer product is numerically identical and
  // tests cannot detect that regression, but at D=4 it changes 16384 elements
  // into 2.68e8 elements, a 16384x larger pass. This relies on the empty
  // contraction axes here; do not generalize it to non-empty tensordot.
  apply_swap_form(bra, forms.bra);
  ftensor<tensor> doubled =
      tensordot(bra, ket_Tn, mptensor::Axes(), mptensor::Axes());
  mptensor::Axes interleaved;
  for (int ax = 0; ax < 4; ++ax) {
    interleaved.push(5 + ax);
    interleaved.push(ax);
  }
  interleaved.push(9);
  interleaved.push(4);
  transpose_with_swap_form(doubled, forms.cross, interleaved);
  mptensor::Shape sh;
  for (std::size_t ax = 0; ax < 4; ++ax) {
    sh.push(doubled.shape()[2 * ax] * doubled.shape()[2 * ax + 1]);
  }
  sh.push(doubled.shape()[8]);
  sh.push(doubled.shape()[9]);
  return mptensor::reshape(doubled.t, sh);
}

/*!
 * @brief Reference fuse of a two-site cluster: materializes the rank-16
 *        outer product.
 *
 * The two-site analogue of doubled_pipeline(): outer product of the bra
 * and ket pair states, joint-swap dressing, (ket, bra) interleave, fusion
 * of the six open legs, then the plain trace over the two physical leg
 * pairs. The rank-16 intermediate makes it the slow path — it exists as
 * the oracle fuse_doubled_cluster_direct() is pinned against.
 *
 * @param[in] bra_pair Rank-8 bra pair state, ALREADY conjugated by the
 *            caller (unlike doubled_pipeline(), which conjugates
 *            internally).
 * @param[in] ket_pair Rank-8 ket pair state with the operator applied.
 * @param[in] leg_ids Direction labels of the six open legs, in the order
 *            of build_reduced_pair_impl().
 * @return Rank-6 blob; leg order as in the file description.
 */
template <class tensor>
tensor fuse_doubled_cluster_naive(const ftensor<tensor>& bra_pair,
                                  const ftensor<tensor>& ket_pair,
                                  const std::vector<int>& leg_ids) {
  constexpr std::size_t kExternalLegs = 6;
  const std::vector<int> cluster_axes = {0, 1, 2, 4, 5, 6};
  std::vector<int> bra_axes;
  std::vector<int> ket_axes;
  for (const int ax : cluster_axes) {
    bra_axes.push_back(ax);
    ket_axes.push_back(ax + 8);
  }
  const auto forms = joint_swap_forms(bra_axes, ket_axes, leg_ids);
  ftensor<tensor> bra = bra_pair;
  // This pre-outer-product swap is the point of this task: it cuts rank-16
  // sweeps. Moving it after the outer product is numerically identical and
  // tests cannot detect that regression, but at D=4 it changes 16384 elements
  // into 2.68e8 elements, a 16384x larger pass. This relies on the empty
  // contraction axes here; do not generalize it to non-empty tensordot.
  apply_swap_form(bra, forms.bra);
  ftensor<tensor> doubled =
      tensordot(bra, ket_pair, mptensor::Axes(), mptensor::Axes());

  mptensor::Axes interleaved;
  for (std::size_t i = 0; i < kExternalLegs; ++i) {
    interleaved.push(ket_axes[i]);
    interleaved.push(bra_axes[i]);
  }
  interleaved.push(11);
  interleaved.push(15);
  interleaved.push(3);
  interleaved.push(7);
  transpose_with_swap_form(doubled, forms.cross, interleaved);

  mptensor::Shape sh;
  for (std::size_t ax = 0; ax < kExternalLegs; ++ax) {
    sh.push(doubled.shape()[2 * ax] * doubled.shape()[2 * ax + 1]);
  }
  sh.push(doubled.shape()[12]);
  sh.push(doubled.shape()[13]);
  sh.push(doubled.shape()[14]);
  sh.push(doubled.shape()[15]);

  tensor fused = mptensor::reshape(doubled.t, sh);
  return mptensor::contract(fused, mptensor::Axes(6, 7), mptensor::Axes(8, 9));
}

/*!
 * @brief Outer-product-free evaluation of fuse_doubled_cluster_naive().
 *
 * Mathematically the rank-16 outer product + sign-dressed interleave
 * transpose + fuse + physical trace equals a direct plain contraction of
 * the physical legs (a rank-12, blob-sized intermediate) once every
 * pairwise sign term of the cluster path — forms.cross plus the Koszul
 * terms of the interleave transpose — is redistributed by where its two
 * legs live:
 *   - both on the bra layer: mask on the bra pair tensor,
 *   - both on the ket layer: mask on the ket pair tensor,
 *   - one leg is traced: the delta trace equates a traced leg with its
 *     twin on the other layer, so the term is rewritten onto the twin
 *     (twin x twin collapses to a linear parity mask,
 *     @f$p(s)^2 = p(s)@f$),
 *   - open bra leg x open ket leg: mask on the rank-12 result.
 *
 * Pinned elementwise against build_reduced_pair_naive() by the
 * impurity_blob tests; derivation record in misc/fermion_twosite_blob/
 * (design.md rev. 2, and probe_direct.cpp, which spells the identity out
 * more compactly). Parameters and return as in
 * fuse_doubled_cluster_naive().
 */
template <class tensor>
tensor fuse_doubled_cluster_direct(const ftensor<tensor>& bra_pair,
                                   const ftensor<tensor>& ket_pair,
                                   const std::vector<int>& leg_ids) {
  constexpr std::size_t kExternalLegs = 6;
  const std::vector<int> cluster_axes = {0, 1, 2, 4, 5, 6};
  std::vector<int> bra_axes;
  std::vector<int> ket_axes;
  for (const int ax : cluster_axes) {
    bra_axes.push_back(ax);
    ket_axes.push_back(ax + 8);
  }
  const auto forms = joint_swap_forms(bra_axes, ket_axes, leg_ids);

  mptensor::Axes interleaved;
  for (std::size_t i = 0; i < kExternalLegs; ++i) {
    interleaved.push(ket_axes[i]);
    interleaved.push(bra_axes[i]);
  }
  interleaved.push(11);
  interleaved.push(15);
  interleaved.push(3);
  interleaved.push(7);

  SwapForm total = forms.cross;
  const SwapForm koszul = transpose_sign_form(interleaved);
  for (const auto& pr : koszul.terms()) {
    total.toggle(pr.first, pr.second);
  }

  auto is_bra = [](int ax) { return ax < 8; };
  auto trace_twin = [](int ax) {
    if (ax == 3) {
      return 11;
    }
    if (ax == 7) {
      return 15;
    }
    if (ax == 11) {
      return 3;
    }
    if (ax == 15) {
      return 7;
    }
    return -1;
  };
  auto to_ket_local = [](int ax) { return ax - 8; };
  auto post_pos = [&](int ax) {
    static const int bra_pos[8] = {0, 1, 2, -1, 3, 4, 5, -1};
    if (is_bra(ax)) {
      return bra_pos[ax];
    }
    return 6 + bra_pos[ax - 8];
  };

  // The twin rewrites below rely on the traced legs carrying the same
  // parity vector on both layers (bra s and its ket twin s'), which holds
  // by construction: conj() keeps parities and the wrapped operator maps
  // each physical leg to one of the same dimension and parity.
  if (bra_pair.parity[3] != ket_pair.parity[11 - 8] ||
      bra_pair.parity[7] != ket_pair.parity[15 - 8]) {
    throw std::runtime_error(
        "fuse_doubled_cluster_direct: traced-leg parity mismatch");
  }

  SwapForm bra_terms = forms.bra;
  SwapForm ket_terms;
  SwapForm post_terms;
  std::vector<LegGauge> ket_gauges;
  for (const auto& pr : total.terms()) {
    int a = pr.first;
    int b = pr.second;
    const bool a_traced = (a == 3 || a == 7 || a == 11 || a == 15);
    const bool b_traced = (b == 3 || b == 7 || b == 11 || b == 15);
    if (a_traced && b_traced && trace_twin(a) == b) {
      const int ax = to_ket_local(a > b ? a : b);
      std::vector<double> factor(ket_pair.parity[ax].size(), 1.0);
      for (std::size_t i = 0; i < factor.size(); ++i) {
        if (ket_pair.parity[ax][i]) {
          factor[i] = -1.0;
        }
      }
      ket_gauges.push_back({ax, factor});
      continue;
    }
    if (a_traced && is_bra(a)) {
      a = trace_twin(a);
    }
    if (b_traced && is_bra(b)) {
      b = trace_twin(b);
    }
    if (is_bra(a) && is_bra(b)) {
      bra_terms.toggle(a, b);
    } else if (!is_bra(a) && !is_bra(b)) {
      ket_terms.toggle(to_ket_local(a), to_ket_local(b));
    } else {
      const int ket_ax = is_bra(a) ? b : a;
      const int bra_ax = is_bra(a) ? a : b;
      if (ket_ax == 11 || ket_ax == 15) {
        bra_terms.toggle(trace_twin(ket_ax), bra_ax);
      } else {
        post_terms.toggle(post_pos(bra_ax), post_pos(ket_ax));
      }
    }
  }

  ftensor<tensor> bra = bra_pair;
  apply_swap_form(bra, bra_terms);
  ftensor<tensor> ket = ket_pair;
  apply_sign_sweep(ket, ket_terms, ket_gauges);

  tensor direct = mptensor::tensordot(bra.t, ket.t, mptensor::Axes(3, 7),
                                      mptensor::Axes(3, 7));
  leg_parities post_parity;
  for (const int ax : cluster_axes) {
    post_parity.push_back(bra.parity[ax]);
  }
  for (const int ax : cluster_axes) {
    post_parity.push_back(ket.parity[ax]);
  }
  ftensor<tensor> post;
  post.t = std::move(direct);
  post.parity = post_parity;
  apply_swap_form(post, post_terms);

  mptensor::Axes perm;
  for (int i = 0; i < static_cast<int>(kExternalLegs); ++i) {
    perm.push(6 + i);
    perm.push(i);
  }
  direct = mptensor::transpose(post.t, perm);
  mptensor::Shape sh;
  const mptensor::Shape ls = direct.shape();
  for (std::size_t ax = 0; ax < kExternalLegs; ++ax) {
    sh.push(ls[2 * ax] * ls[2 * ax + 1]);
  }
  return mptensor::reshape(direct, sh);
}

}  // namespace detail

/*!
 * @brief One-site reduced tensor with the physical legs left open.
 *
 * doubled_pipeline() of a wrapped Tn against itself: rank 6, legs
 * ([l lb], [t tb], [r rb], [b bb], s_ket, s_bra). Contracting a one-site
 * operator into (s_ket, s_bra) and handing the result to the bosonic CTM
 * is the fermionic one-site measurement.
 *
 * @param[in] Tn Rank-5 wrapped center tensor.
 * @throw std::runtime_error If @p Tn is not rank 5.
 */
template <class tensor>
tensor build_reduced_op(const ftensor<tensor>& Tn) {
  if (Tn.rank() != 5) {
    throw std::runtime_error("build_reduced_op expects a five-leg Tn ftensor");
  }
  return detail::doubled_pipeline(Tn, Tn);
}

/*!
 * @brief One-site reduced tensor with the physical legs traced: the rank-4
 *        double-layer tensor the CTM environment is built from.
 */
template <class tensor>
tensor build_reduced(const ftensor<tensor>& Tn) {
  return mptensor::contract(build_reduced_op(Tn), mptensor::Axes(4),
                            mptensor::Axes(5));
}

/*!
 * @brief Contract two neighboring wrapped Tn into a rank-8 pair state.
 *
 * Graded contraction of the shared bond. Output leg order:
 *   - horizontal: (l_A, t_A, b_A, s_A, t_B, r_B, b_B, s_B)
 *   - vertical:   (l_A, t_A, r_A, s_A, l_B, r_B, b_B, s_B)
 *
 * i.e. A's remaining legs in (l, t, r, b, s) order, then B's.
 *
 * @param[in] TnA First site (left, resp. top).
 * @param[in] TnB Second site (right, resp. bottom).
 * @param[in] direction Window orientation.
 * @throw std::runtime_error If either tensor is not rank 5.
 */
template <class tensor>
ftensor<tensor> build_pair_state(const ftensor<tensor>& TnA,
                                 const ftensor<tensor>& TnB,
                                 reduced_pair_direction direction) {
  if (TnA.rank() != 5 || TnB.rank() != 5) {
    throw std::runtime_error("build_pair_state expects five-leg Tn ftensors");
  }
  switch (direction) {
    case reduced_pair_direction::horizontal:
      // (l_A, t_A, b_A, s_A, t_B, r_B, b_B, s_B)
      return tensordot(TnA, TnB, mptensor::Axes(2), mptensor::Axes(0));
    case reduced_pair_direction::vertical:
      // (l_A, t_A, r_A, s_A, l_B, r_B, b_B, s_B)
      return tensordot(TnA, TnB, mptensor::Axes(3), mptensor::Axes(1));
  }
  throw std::runtime_error("build_pair_state: invalid direction");
}

/*!
 * @brief Apply a two-site operator to the physical legs of a pair state.
 *
 * @param[in] pair Rank-8 pair state from build_pair_state(); physical legs
 *            at positions 3 and 7.
 * @param[in] op12 Rank-4 operator with legs (in_A, in_B, out_A, out_B),
 *            loaded via one of the wrappers in fops.hpp.
 * @return Pair state with the operator applied and the original leg order
 *         restored.
 * @throw std::runtime_error On rank mismatch.
 */
template <class tensor>
ftensor<tensor> apply_pair_op(const ftensor<tensor>& pair,
                              const ftensor<tensor>& op12) {
  if (pair.rank() != 8) {
    throw std::runtime_error("apply_pair_op expects an eight-leg pair state");
  }
  if (op12.rank() != 4) {
    throw std::runtime_error("apply_pair_op expects a four-leg operator");
  }
  ftensor<tensor> applied =
      tensordot(pair, op12, mptensor::Axes(3, 7), mptensor::Axes(0, 1));
  return transpose(applied, mptensor::Axes(0, 1, 2, 6, 3, 4, 5, 7));
}

namespace detail {

/*!
 * @brief Shared body of the two blob builders.
 *
 * They differ in exactly one thing — which fuse implementation produces
 * the blob — so the surrounding sequence (pair state, operator insertion,
 * fuse, gauge alignment) lives here once and the equivalence test compares
 * the fuse step alone. The final apply_fused_leg_gauge() calls align the
 * blob with the single-site doubling convention; measured directly via
 * comparison against build_reduced_op() / build_reduced() based direct
 * composition, not derived analytically.
 *
 * @param[in] TnA First site (left, resp. top).
 * @param[in] TnB Second site (right, resp. bottom).
 * @param[in] op12 Two-site operator loaded with wrap_reduced_pair_op().
 * @param[in] direction Window orientation.
 * @param[in] fuse fuse_doubled_cluster_naive() or
 *            fuse_doubled_cluster_direct().
 * @return Rank-6 blob; leg order as in the file description.
 */
template <class tensor, class Fuse>
tensor build_reduced_pair_impl(const ftensor<tensor>& TnA,
                               const ftensor<tensor>& TnB,
                               const ftensor<tensor>& op12,
                               reduced_pair_direction direction, Fuse fuse) {
  const ftensor<tensor> ket_ab = build_pair_state(TnA, TnB, direction);
  const std::vector<int> leg_ids = [direction] {
    switch (direction) {
      case reduced_pair_direction::horizontal:
        return std::vector<int>{0, 1, 3, 1, 2, 3};
      case reduced_pair_direction::vertical:
        return std::vector<int>{0, 1, 2, 0, 2, 3};
    }
    throw std::runtime_error("build_reduced_pair: invalid direction");
  }();

  const ftensor<tensor> ket_op = apply_pair_op(ket_ab, op12);
  tensor ret = fuse(conj(ket_ab), ket_op, leg_ids);
  // Gauge alignment with the single-site doubling convention; measured
  // directly via comparison against build_reduced_op/build_reduced-based
  // direct composition, not derived analytically.
  if (direction == reduced_pair_direction::horizontal) {
    apply_fused_leg_gauge(ret, TnA.parity[3], 2, true);
    apply_fused_leg_gauge(ret, TnB.parity[3], 5, false);
  } else {
    apply_fused_leg_gauge(ret, TnA.parity[0], 0, true);
    apply_fused_leg_gauge(ret, TnB.parity[0], 3, false);
  }
  return ret;
}

}  // namespace detail

/*!
 * @brief Reference blob builder: materializes the rank-16 outer product of
 *        the two layers.
 *
 * Kept as the oracle the direct path is pinned against (and as the frozen
 * convention the Fock-verified tests check); the solver uses
 * build_reduced_pair_direct(). Parameters as in
 * detail::build_reduced_pair_impl().
 */
template <class tensor>
tensor build_reduced_pair_naive(const ftensor<tensor>& TnA,
                                const ftensor<tensor>& TnB,
                                const ftensor<tensor>& op12,
                                reduced_pair_direction direction) {
  return detail::build_reduced_pair_impl(
      TnA, TnB, op12, direction, &detail::fuse_doubled_cluster_naive<tensor>);
}

/*!
 * @brief Direct-contraction blob builder, used by the solver.
 *
 * Contracts the physical legs of the bra and ket pair layers directly (a
 * rank-12, blob-sized intermediate) after redistributing the frozen
 * joint-swap and transpose Koszul sign terms onto the two layers and the
 * contraction result (see detail::fuse_doubled_cluster_direct()). Returns
 * the same rank-6 blob as build_reduced_pair_naive() elementwise; the
 * rank-16 outer product is never materialized.
 * Design: misc/fermion_twosite_blob/design.md (rev. 2, "direct fuse").
 * Parameters as in detail::build_reduced_pair_impl().
 */
template <class tensor>
tensor build_reduced_pair_direct(const ftensor<tensor>& TnA,
                                 const ftensor<tensor>& TnB,
                                 const ftensor<tensor>& op12,
                                 reduced_pair_direction direction) {
  return detail::build_reduced_pair_impl(
      TnA, TnB, op12, direction, &detail::fuse_doubled_cluster_direct<tensor>);
}

}  // namespace tenes::fermion

#endif  // TENES_SRC_FERMION_REDUCED_HPP_
