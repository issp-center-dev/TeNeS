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
 *  One-site tensors are folded with a consistent planar routing: all six
 *  ordered direction pairs ending at the left-bottom corner dress the open
 *  bra x ket network before each virtual pair is fused. Two-site operator
 *  blobs use the bundled-k construction: a graded SVD of the operator,
 *  attachment of its factors to the sites, fusion of the channel k into the
 *  shared virtual bond, the required crossing mask, and the same asymmetric
 *  one-site fold on both sites.
 *
 *  Thus every fermionic crossing is resolved before the resulting reduced
 *  tensors are handed to the existing bosonic CTM contraction machinery.
 *  Keeping the construction open until its explicit physical contractions
 *  also avoids the spurious supertrace signs produced by mechanical graded
 *  bookkeeping on closed loops.
 *
 *  build_reduced() / build_reduced_op() produce the one-site reduced
 *  tensor. build_reduced_pair_halves() stops at the two folded sites and
 *  returns them as a reduced_pair_halves; that is what the production
 *  measurement path consumes, absorbing each half into its side of the CTM
 *  environment without ever forming their product (see
 *  reduced_measure.hpp). build_reduced_pair_direct() joins the halves over
 *  the shared bond with a plain contraction to produce the two-site "blob",
 *  and build_reduced_pair_naive() builds that same blob through an
 *  independent transcription of the construction; both exist so the tests
 *  have a materialized object to pin. Blob leg orders (each external leg is
 *  a fused (ket, bra) pair of dimension @f$D^2@f$):
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
#include <string>
#include <vector>

#include "../timer.hpp"
#include "fops.hpp"
#include "ftensor.hpp"

namespace tenes::fermion {

//! Orientation of a two-site window: horizontal (A left of B) or vertical
//! (A above B).
enum class reduced_pair_direction { horizontal, vertical };

/*! @brief The two folded site halves immediately before their shared-bond
 *         contraction.
 *
 * The bundled channel is fused into the shared virtual leg. Thus PA/PB have
 * leg orders (L_A,T_A,[R_A k],B_A) / ([L_B k],T_B,R_B,B_B) horizontally and
 * (L_A,T_A,R_A,[B_A k]) / (L_B,[T_B k],R_B,B_B) vertically. The shared axes
 * are derived from direction so inconsistent axis metadata cannot be stored.
 */
template <class tensor>
struct reduced_pair_halves {
  //! Folded first (left or top) site.
  tensor PA;
  //! Folded second (right or bottom) site.
  tensor PB;
  //! Orientation of the two-site window.
  reduced_pair_direction direction;

  //! Shared bundled-bond axis of PA.
  std::size_t axis_a() const {
    return direction == reduced_pair_direction::horizontal ? 2 : 3;
  }

  //! Shared bundled-bond axis of PB.
  std::size_t axis_b() const {
    return direction == reduced_pair_direction::horizontal ? 0 : 1;
  }
};

namespace detail {

//! Bit position encoding the ORDERED direction-label pair (x, y),
//! x != y, in kDoubledJointMask (directions 0=left, 1=top, 2=right,
//! 3=bottom).
constexpr int joint_bit(int x, int y) { return x * 3 + (y < x ? y : y - 1); }

/*! @brief Which direction-label pairs generate joint-swap terms.
 *
 * All six ordered pairs (x, y) with y in {0 = left, 3 = bottom} are set:
 * (1,0), (2,0), (3,0), (0,3), (1,3), (2,3). Geometrically this is one of
 * the four consistent planar routings of the bra legs around the site
 * (here: past the left-bottom corner); the four routings differ only by
 * gauge on the fused legs and any one of them, applied to EVERY pair,
 * folds the double layer exactly for arbitrary patches (verified against
 * the exact Fock oracle on 1xN, 2x2, 2x3, 3x2, 3x3 geometries in
 * work/fermion/ctm-fold-check).
 *
 * The previous mask {(0,3), (1,0), (2,3), (3,0)} was missing (1,3) and
 * (2,0): an inconsistent mixture that happens to be gauge-equivalent to a
 * consistent routing on 1xN chains and on the 2x2 plaquette — every
 * geometry the oracle pinning had used — but is wrong as soon as a site
 * has three or more nontrivial legs, i.e. on every true 2D network.
 */
constexpr unsigned kDoubledJointMask =
    (1u << joint_bit(1, 0)) | (1u << joint_bit(2, 0)) |
    (1u << joint_bit(3, 0)) | (1u << joint_bit(0, 3)) |
    (1u << joint_bit(1, 3)) | (1u << joint_bit(2, 3));

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
 * @brief Shared single-site double-layer pipeline.
 *
 * Outer product of bra and ket layers, frozen joint-swap dressing,
 * (ket, bra) interleave, column-major fusion. Output legs:
 * ([l lb], [t tb], [r rb], [b bb], s_ket, s_bra) — the four virtual legs
 * fused (ket first and fastest-varying), the two physical legs left open.
 *
 * The bra and ket arguments may be DIFFERENT tensors. In particular the
 * virtual legs may differ between the layers in both dimension and parity
 * ledger — every joint-swap and interleave step here is per-leg and
 * per-layer, so nothing requires the two layers to match. This asymmetry
 * is exactly what the two-site operator blob relies on: its ket layer
 * carries an operator factor with the SVD channel leg bundled into a bond
 * leg (fused ledger), while its bra layer is the bare site tensor. Only
 * the two physical legs must agree in dimension, because the caller
 * closes them against each other (build_reduced()) or hands them to a
 * physical operator as a matching (s_ket, s_bra) pair. Both layers must
 * be parity-even under their own ledgers.
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
 * @brief Physical-traced single-site double-layer pipeline.
 *
 * Contracts the matching physical legs before the frozen joint-swap
 * dressing, (ket, bra) interleave, and column-major fusion. Output legs:
 * ([l lb], [t tb], [r rb], [b bb]) — the four virtual legs fused with ket
 * first and fastest-varying. The bra and ket layers may have different
 * virtual dimensions and parity ledgers, but their physical legs must match.
 * Both layers must be parity-even under their own ledgers; the equivalence
 * to doubled_pipeline() followed by a plain physical trace relies on this
 * premise.
 *
 * @param[in] bra_Tn Rank-5 wrapped Tn to conjugate into the bra layer.
 * @param[in] ket_Tn Rank-5 wrapped Tn forming the ket layer.
 * @return Rank-4 plain tensor with the four fused virtual legs.
 * @throw std::runtime_error In debug builds, if either layer violates even
 *        parity above the numerical tolerance.
 */
template <class tensor>
tensor doubled_pipeline_traced(const ftensor<tensor>& bra_Tn,
                               const ftensor<tensor>& ket_Tn) {
#ifndef NDEBUG
  const double bra_violation = parity_violation(bra_Tn);
  const double bra_threshold = 1.0e-10 * std::max(1.0, max_abs(bra_Tn));
  if (bra_violation > bra_threshold) {
    throw std::runtime_error(
        "doubled_pipeline_traced: bra layer is not parity even");
  }
  const double ket_violation = parity_violation(ket_Tn);
  const double ket_threshold = 1.0e-10 * std::max(1.0, max_abs(ket_Tn));
  if (ket_violation > ket_threshold) {
    throw std::runtime_error(
        "doubled_pipeline_traced: ket layer is not parity even");
  }
#endif
  const auto forms = joint_swap_forms({0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 2, 3});
  ftensor<tensor> bra = conj(bra_Tn);
  apply_swap_form(bra, forms.bra);
  std::vector<double> phys_sign(bra.parity[4].size());
  for (std::size_t i = 0; i < phys_sign.size(); ++i) {
    phys_sign[i] = bra.parity[4][i] ? -1.0 : 1.0;
  }
  bra.multiply_vector(phys_sign, 4);
  ftensor<tensor> doubled =
      tensordot(bra, ket_Tn, mptensor::Axes(4), mptensor::Axes(4));
  transpose_with_swap_form(doubled, forms.cross,
                           mptensor::Axes(4, 0, 5, 1, 6, 2, 7, 3));
  mptensor::Shape sh;
  for (std::size_t ax = 0; ax < 4; ++ax) {
    sh.push(doubled.shape()[2 * ax] * doubled.shape()[2 * ax + 1]);
  }
  return mptensor::reshape(doubled.t, sh);
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

/*!
 * @brief Reference bundled-k two-site operator blob builder.
 *
 * Implements the SVD, attach, bundle, crossing-mask, asymmetric-fold, and
 * final-contraction construction literally, and folds each site with the
 * physical legs left open (detail::doubled_pipeline() followed by a plain
 * physical trace).
 *
 * The solver path has since diverged: build_reduced_pair_halves() folds with
 * the physical legs contracted first (detail::doubled_pipeline_traced()) and
 * stops before the join. That divergence is the point — the reference and the
 * solver must not share construction code, so T13 compares two independent
 * implementations rather than a function against itself. Do not deduplicate
 * them, and keep future optimizations out of this one. The operator must be
 * loaded with wrap_twosite_gate() (input swap only).
 */
template <class tensor>
tensor build_reduced_pair_naive(const ftensor<tensor>& TnA,
                                const ftensor<tensor>& TnB,
                                const ftensor<tensor>& op12,
                                reduced_pair_direction direction) {
  if (TnA.rank() != 5 || TnB.rank() != 5 || op12.rank() != 4) {
    throw std::runtime_error(
        "build_reduced_pair_naive expects rank-5 sites and a rank-4 gate");
  }
  if (direction != reduced_pair_direction::horizontal &&
      direction != reduced_pair_direction::vertical) {
    throw std::runtime_error("build_reduced_pair_naive: invalid direction");
  }

  ftensor<tensor> u, vt;
  std::vector<double> s;
  const int info =
      svd(op12, mptensor::Axes(0, 2), mptensor::Axes(1, 3), u, s, vt);
  if (info != 0) {
    throw std::runtime_error("build_reduced_pair_naive: gate SVD failed");
  }
  u.multiply_vector(s, 2);
  const ftensor<tensor> TA6 =
      tensordot(TnA, u, mptensor::Axes(4), mptensor::Axes(0));
  const ftensor<tensor> TB6 =
      tensordot(TnB, vt, mptensor::Axes(4), mptensor::Axes(1));

  const std::size_t nk = s.size();
  const std::size_t bond_axis =
      direction == reduced_pair_direction::horizontal ? 2 : 3;
  const parity_vector& bond_parity = TnA.parity[bond_axis];
  const parity_vector& k_parity = u.parity[2];
  std::vector<double> crossing_mask(bond_parity.size() * nk, 1.0);
  for (std::size_t k = 0; k < nk; ++k) {
    for (std::size_t b = 0; b < bond_parity.size(); ++b) {
      if (bond_parity[b] && k_parity[k]) {
        crossing_mask[b + bond_parity.size() * k] = -1.0;
      }
    }
  }

  ftensor<tensor> TA5;
  ftensor<tensor> TB5;
  std::size_t contract_a;
  std::size_t contract_b;
  if (direction == reduced_pair_direction::horizontal) {
    TA5 = reshape(
        transpose(TA6, mptensor::Axes(0, 1, 2, 5, 3, 4)),
        mptensor::Shape(TA6.shape()[0], TA6.shape()[1], TA6.shape()[2] * nk,
                        TA6.shape()[3], TA6.shape()[4]));
    TA5.multiply_vector(crossing_mask, 2);
    TB5 = reshape(
        transpose(TB6, mptensor::Axes(0, 4, 1, 2, 3, 5)),
        mptensor::Shape(TB6.shape()[0] * nk, TB6.shape()[1], TB6.shape()[2],
                        TB6.shape()[3], TB6.shape()[5]));
    contract_a = 2;
    contract_b = 0;
  } else {
    TA5 =
        reshape(transpose(TA6, mptensor::Axes(0, 1, 2, 3, 5, 4)),
                mptensor::Shape(TA6.shape()[0], TA6.shape()[1], TA6.shape()[2],
                                TA6.shape()[3] * nk, TA6.shape()[4]));
    TA5.multiply_vector(crossing_mask, 3);
    TB5 = reshape(
        transpose(TB6, mptensor::Axes(0, 1, 4, 2, 3, 5)),
        mptensor::Shape(TB6.shape()[0], TB6.shape()[1] * nk, TB6.shape()[2],
                        TB6.shape()[3], TB6.shape()[5]));
    contract_a = 3;
    contract_b = 1;
  }

  const tensor PA = mptensor::contract(detail::doubled_pipeline(TnA, TA5),
                                       mptensor::Axes(4), mptensor::Axes(5));
  const tensor PB = mptensor::contract(detail::doubled_pipeline(TnB, TB5),
                                       mptensor::Axes(4), mptensor::Axes(5));
  return mptensor::tensordot(PA, PB, mptensor::Axes(contract_a),
                             mptensor::Axes(contract_b));
}

/*!
 * @brief Build the folded halves of a bundled-k two-site operator network.
 *
 * Performs the direct builder's SVD, attachment, bundled-k crossing mask, and
 * asymmetric folds, but leaves the final shared bond open. No numerical
 * singular-value truncation is performed. The operator must be loaded with
 * wrap_twosite_gate().
 *
 * @return The two rank-4 folded halves and their orientation.
 * @throw std::runtime_error On rank mismatch, invalid direction, or failed
 *        gate SVD.
 */
template <class tensor>
reduced_pair_halves<tensor> build_reduced_pair_halves(
    const ftensor<tensor>& TnA, const ftensor<tensor>& TnB,
    const ftensor<tensor>& op12, reduced_pair_direction direction) {
  ::tenes::ScopedTimer scoped_timer("measure/twosite/halves");
  if (TnA.rank() != 5 || TnB.rank() != 5 || op12.rank() != 4) {
    throw std::runtime_error(
        "build_reduced_pair_halves expects rank-5 sites and a rank-4 gate");
  }
  if (direction != reduced_pair_direction::horizontal &&
      direction != reduced_pair_direction::vertical) {
    throw std::runtime_error("build_reduced_pair_halves: invalid direction");
  }

  ftensor<tensor> u, vt;
  std::vector<double> s;
  const int info =
      svd(op12, mptensor::Axes(0, 2), mptensor::Axes(1, 3), u, s, vt);
  if (info != 0) {
    throw std::runtime_error("build_reduced_pair_halves: gate SVD failed");
  }
  u.multiply_vector(s, 2);
  const ftensor<tensor> TA6 =
      tensordot(TnA, u, mptensor::Axes(4), mptensor::Axes(0));
  const ftensor<tensor> TB6 =
      tensordot(TnB, vt, mptensor::Axes(4), mptensor::Axes(1));

  const std::size_t nk = s.size();
  const std::size_t bond_axis =
      direction == reduced_pair_direction::horizontal ? 2 : 3;
  const parity_vector& bond_parity = TnA.parity[bond_axis];
  const parity_vector& k_parity = u.parity[2];
  std::vector<double> crossing_mask(bond_parity.size() * nk, 1.0);
  for (std::size_t k = 0; k < nk; ++k) {
    for (std::size_t b = 0; b < bond_parity.size(); ++b) {
      if (bond_parity[b] && k_parity[k]) {
        crossing_mask[b + bond_parity.size() * k] = -1.0;
      }
    }
  }

  ftensor<tensor> TA5;
  ftensor<tensor> TB5;
  if (direction == reduced_pair_direction::horizontal) {
    TA5 = reshape(
        transpose(TA6, mptensor::Axes(0, 1, 2, 5, 3, 4)),
        mptensor::Shape(TA6.shape()[0], TA6.shape()[1], TA6.shape()[2] * nk,
                        TA6.shape()[3], TA6.shape()[4]));
    TA5.multiply_vector(crossing_mask, 2);
    TB5 = reshape(
        transpose(TB6, mptensor::Axes(0, 4, 1, 2, 3, 5)),
        mptensor::Shape(TB6.shape()[0] * nk, TB6.shape()[1], TB6.shape()[2],
                        TB6.shape()[3], TB6.shape()[5]));
  } else {
    TA5 =
        reshape(transpose(TA6, mptensor::Axes(0, 1, 2, 3, 5, 4)),
                mptensor::Shape(TA6.shape()[0], TA6.shape()[1], TA6.shape()[2],
                                TA6.shape()[3] * nk, TA6.shape()[4]));
    TA5.multiply_vector(crossing_mask, 3);
    TB5 = reshape(
        transpose(TB6, mptensor::Axes(0, 1, 4, 2, 3, 5)),
        mptensor::Shape(TB6.shape()[0], TB6.shape()[1] * nk, TB6.shape()[2],
                        TB6.shape()[3], TB6.shape()[5]));
  }

  return {detail::doubled_pipeline_traced(TnA, TA5),
          detail::doubled_pipeline_traced(TnB, TB5), direction};
}

/*!
 * @brief Solver bundled-k two-site operator blob builder.
 *
 * Thin compatibility wrapper over build_reduced_pair_halves(). The final
 * bond contraction is deliberately plain mptensor::tensordot: a graded
 * contraction would introduce an extra closed-loop supertrace sign. The
 * returned leg order is unchanged from build_reduced_pair_naive().
 */
template <class tensor>
tensor build_reduced_pair_direct(const ftensor<tensor>& TnA,
                                 const ftensor<tensor>& TnB,
                                 const ftensor<tensor>& op12,
                                 reduced_pair_direction direction) {
  reduced_pair_halves<tensor> halves;
  try {
    halves = build_reduced_pair_halves(TnA, TnB, op12, direction);
  } catch (const std::runtime_error& error) {
    const std::string message = error.what();
    if (message ==
        "build_reduced_pair_halves expects rank-5 sites and a rank-4 gate") {
      throw std::runtime_error(
          "build_reduced_pair_direct expects rank-5 sites and a rank-4 gate");
    }
    if (message == "build_reduced_pair_halves: invalid direction") {
      throw std::runtime_error("build_reduced_pair_direct: invalid direction");
    }
    if (message == "build_reduced_pair_halves: gate SVD failed") {
      throw std::runtime_error("build_reduced_pair_direct: gate SVD failed");
    }
    throw;
  }
  return mptensor::tensordot(halves.PA, halves.PB,
                             mptensor::Axes(halves.axis_a()),
                             mptensor::Axes(halves.axis_b()));
}

}  // namespace tenes::fermion

#endif  // TENES_SRC_FERMION_REDUCED_HPP_
