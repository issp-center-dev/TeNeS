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
 *  @brief Measurement contractions over the fermionic reduced tensors.
 *
 *  Two environment flavors are served: the mean-field path
 *  (contract_pair_MF(), on lambda-dressed pair states with no environment
 *  tensors) and the CTM path
 *  (contract_reduced_pair_halves_density_CTM(), absorbing the two folded
 *  halves from reduced.hpp into the corner and edge tensors of the standard
 *  6-edge two-site CTM window). The blob-closing counterparts
 *  contract_reduced_pair_{horizontal,vertical}_density_CTM() take a
 *  materialized rank-6 blob instead; the solver no longer forms one, and
 *  they remain as the independent reference the tests close against.
 *  Helpers build the lambda-dressed site tensors and the per-site reduced
 *  density tensors the drivers hand to these contractions.
 */

#ifndef TENES_SRC_FERMION_REDUCED_MEASURE_HPP_
#define TENES_SRC_FERMION_REDUCED_MEASURE_HPP_

#include <stdexcept>
#include <vector>

#include "../timer.hpp"
#include "fermion_info.hpp"
#include "reduced.hpp"

namespace tenes::fermion {

/*!
 * @brief Multiply the mean-field weights of all four bonds onto one Tn.
 *
 * @param[in] Tn Plain rank-5 center tensor.
 * @param[in] lambda One weight vector per virtual leg (l, t, r, b).
 * @return Dressed copy of @p Tn.
 */
template <class tensor>
tensor lambda_dressed_tensor(const tensor& Tn,
                             const std::vector<std::vector<double>>& lambda) {
  tensor dressed = Tn;
  for (int leg = 0; leg < 4; ++leg) {
    dressed.multiply_vector(lambda[leg], leg);
  }
  return dressed;
}

//! lambda_dressed_tensor() applied to every site of the unit cell.
template <class tensor>
std::vector<tensor> lambda_dressed_tensors(
    const std::vector<tensor>& Tn,
    const std::vector<std::vector<std::vector<double>>>& lambda) {
  std::vector<tensor> dressed;
  dressed.reserve(Tn.size());
  for (int site = 0; site < static_cast<int>(Tn.size()); ++site) {
    dressed.push_back(lambda_dressed_tensor(Tn[site], lambda[site]));
  }
  return dressed;
}

/*!
 * @brief One-site reduced tensors (physical legs open) for every site.
 *
 * build_reduced_op() of the wrapped Tn, per site — the "density" tensors
 * one-site fermionic measurement contracts an operator into.
 */
template <class tensor>
std::vector<tensor> build_reduced_density_tensors(const std::vector<tensor>& Tn,
                                                  const FermionInfo& finfo) {
  std::vector<tensor> reduced;
  reduced.reserve(Tn.size());
  for (int site = 0; site < static_cast<int>(Tn.size()); ++site) {
    reduced.push_back(build_reduced_op(wrap_Tn(Tn[site], finfo, site)));
  }
  return reduced;
}

/*!
 * @brief Folded halves of a two-site identity network.
 *
 * Builds one physical-traced reduced tensor per site and leaves their shared
 * bond open. The rank-4 leg orders and derived shared axes match
 * build_reduced_pair_halves(), and so does the fold: both go through
 * detail::doubled_pipeline_traced(), which contracts the physical legs
 * before the fold rather than after.
 *
 * @param[in] TnA,TnB Rank-5 wrapped center tensors of the two sites.
 * @param[in] direction Orientation of the window.
 * @return The two folded identity halves and their orientation.
 * @throw std::runtime_error On rank mismatch or invalid direction; in debug
 *        builds, also if a site violates even parity (see
 *        detail::doubled_pipeline_traced()).
 */
template <class tensor>
reduced_pair_halves<tensor> build_reduced_identity_halves(
    const ftensor<tensor>& TnA, const ftensor<tensor>& TnB,
    reduced_pair_direction direction) {
  ::tenes::ScopedTimer scoped_timer("measure/twosite/halves");
  if (TnA.rank() != 5 || TnB.rank() != 5) {
    // Stated here because the fold below no longer runs through
    // build_reduced_op(), which used to supply this guard on the way past.
    // Without it a rank-4 site that happens to be parity even reaches the
    // fold and aborts the process instead of raising.
    throw std::runtime_error(
        "build_reduced_identity_halves expects rank-5 sites");
  }
  if (direction != reduced_pair_direction::horizontal &&
      direction != reduced_pair_direction::vertical) {
    // Same contract as build_reduced_pair_halves(): a direction outside the
    // enumeration must not be stored, because axis_a() / axis_b() and the
    // closure dispatcher all treat "not horizontal" as vertical and would
    // silently contract the wrong network.
    throw std::runtime_error(
        "build_reduced_identity_halves: invalid direction");
  }
  return {detail::doubled_pipeline_traced(TnA, TnA),
          detail::doubled_pipeline_traced(TnB, TnB), direction};
}

/*!
 * @brief Norm blob of a two-site window: the identity-operator counterpart
 *        of build_reduced_pair_direct().
 *
 * Thin wrapper over build_reduced_identity_halves() that joins the halves so
 * the tests have an object to pin; the solver stops at the halves. The final
 * shared-bond contraction remains a plain mptensor::tensordot and preserves
 * the established rank-6 blob leg order.
 *
 * @throw std::runtime_error Whatever build_reduced_identity_halves() throws,
 *        under that function's name: this wrapper adds no guards of its own.
 */
template <class tensor>
tensor build_reduced_identity_pair(const ftensor<tensor>& TnA,
                                   const ftensor<tensor>& TnB,
                                   reduced_pair_direction direction) {
  const auto halves = build_reduced_identity_halves(TnA, TnB, direction);
  return mptensor::tensordot(halves.PA, halves.PB,
                             mptensor::Axes(halves.axis_a()),
                             mptensor::Axes(halves.axis_b()));
}

namespace detail {

//! The axes list (0, 1, ..., rank-1).
inline mptensor::Axes all_axes(int rank) {
  mptensor::Axes axes;
  for (int ax = 0; ax < rank; ++ax) {
    axes.push(ax);
  }
  return axes;
}

}  // namespace detail

/*!
 * @brief Mean-field norm of a two-site pair state:
 *        @f$\langle pair | pair \rangle@f$ as a graded full contraction.
 *
 * No environment tensors: the lambda weights on the open legs are expected
 * to be multiplied into TnA / TnB beforehand (the same dressing the
 * bosonic mean-field path applies in measure_twosite).
 */
template <class tensor>
typename tensor::value_type contract_pair_MF(const ftensor<tensor>& pair) {
  const mptensor::Axes axes = detail::all_axes(pair.rank());
  return trace(conj(pair), pair, axes, axes);
}

/*!
 * @brief Mean-field expectation value (unnormalized) of op12 on a pair
 *        state.
 *
 * @warning op12 must be loaded with wrap_twosite_gate() (input-leg swap
 *          only), the same convention used by the bundled-k CTM blob.
 */
template <class tensor>
typename tensor::value_type contract_pair_MF(const ftensor<tensor>& pair,
                                             const ftensor<tensor>& op12) {
  const mptensor::Axes axes = detail::all_axes(pair.rank());
  return trace(conj(pair), apply_pair_op(pair, op12), axes, axes);
}

namespace detail {

//! Close the four remaining boundary legs pairwise: the sum over locally
//! stored elements with idx[0] == idx[1] and idx[2] == idx[3] (a double
//! delta trace).
template <class tensor>
typename tensor::value_type trace_boundary_pairs(const tensor& a) {
  if (a.rank() != 4) {
    throw std::runtime_error("expected four boundary legs");
  }
  typename tensor::value_type value = 0.0;
  mptensor::Index idx;
  idx.resize(4);
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    a.global_index_fast(n, idx);
    if (idx[0] == idx[1] && idx[2] == idx[3]) {
      value += a[n];
    }
  }
  return value;
}

}  // namespace detail

/*!
 * @brief Close a horizontal two-site blob with its CTM environment.
 *
 * The window (corners clockwise C1..C4 from the top left, edges eT1..eT6
 * clockwise from the top edge of the left site):
 *
 * @verbatim
   C1  eT1 eT2  C2
   eT6 [A   B]  eT3
   C4  eT5 eT4  C3
   @endverbatim
 *
 * @param[in] C1,C2,C3,C4 Corner transfer matrices as sketched above.
 * @param[in] eT1,eT2,eT3,eT4,eT5,eT6 Edge tensors as sketched above.
 * @param[in] blob Rank-6 blob (L_A, T_A, B_A, T_B, R_B, B_B) from
 *            build_reduced_pair_direct() or
 *            build_reduced_identity_pair().
 * @return The scalar value of the closed network (divide an operator
 *         blob's value by the identity blob's to normalize).
 */
template <class tensor>
typename tensor::value_type contract_reduced_pair_horizontal_density_CTM(
    const tensor& C1, const tensor& C2, const tensor& C3, const tensor& C4,
    const tensor& eT1, const tensor& eT2, const tensor& eT3, const tensor& eT4,
    const tensor& eT5, const tensor& eT6, const tensor& blob) {
  ::tenes::ScopedTimer scoped_timer("measure/twosite/absorb");
  using mptensor::Axes;
  const tensor left_lower = tensordot(
      eT5,
      tensordot(C1, tensordot(C4, eT6, Axes(1), Axes(0)), Axes(0), Axes(1)),
      Axes(1), Axes(1));
  tensor work = tensordot(blob, left_lower, Axes(0, 2), Axes(3, 1));
  work = tensordot(eT1, work, Axes(0, 2), Axes(5, 0));

  const tensor right_lower = tensordot(
      eT4,
      tensordot(C2, tensordot(C3, eT3, Axes(0), Axes(1)), Axes(1), Axes(1)),
      Axes(0), Axes(1));
  work = tensordot(work, right_lower, Axes(2, 3), Axes(3, 1));
  work = tensordot(eT2, work, Axes(1, 2), Axes(4, 1));
  return detail::trace_boundary_pairs(work);
}

/*!
 * @brief Close a vertical two-site blob with its CTM environment.
 *
 * The window (corners clockwise C1..C4 from the top left, edges eT1..eT6
 * clockwise from the top edge of the top site):
 *
 * @verbatim
   C1  eT1  C2
   eT6 [A]  eT2
   eT5 [B]  eT3
   C4  eT4  C3
   @endverbatim
 *
 * @param[in] C1,C2,C3,C4 Corner transfer matrices as sketched above.
 * @param[in] eT1,eT2,eT3,eT4,eT5,eT6 Edge tensors as sketched above.
 * @param[in] blob Rank-6 blob (L_A, T_A, R_A, L_B, R_B, B_B) from
 *            build_reduced_pair_direct() or
 *            build_reduced_identity_pair().
 * @return The scalar value of the closed network.
 */
template <class tensor>
typename tensor::value_type contract_reduced_pair_vertical_density_CTM(
    const tensor& C1, const tensor& C2, const tensor& C3, const tensor& C4,
    const tensor& eT1, const tensor& eT2, const tensor& eT3, const tensor& eT4,
    const tensor& eT5, const tensor& eT6, const tensor& blob) {
  ::tenes::ScopedTimer scoped_timer("measure/twosite/absorb");
  using mptensor::Axes;
  const tensor top_left = tensordot(
      eT6,
      tensordot(C1, tensordot(C2, eT1, Axes(0), Axes(1)), Axes(1), Axes(1)),
      Axes(1), Axes(0));
  tensor work = tensordot(blob, top_left, Axes(0, 1), Axes(1, 3));
  work = tensordot(eT2, work, Axes(0, 2), Axes(5, 0));

  const tensor bottom_right = tensordot(
      eT5,
      tensordot(C3, tensordot(C4, eT4, Axes(0), Axes(1)), Axes(1), Axes(1)),
      Axes(0), Axes(1));
  work = tensordot(work, bottom_right, Axes(1, 3), Axes(1, 3));
  work = tensordot(eT3, work, Axes(1, 2), Axes(4, 1));
  return detail::trace_boundary_pairs(work);
}

/*!
 * @brief Absorb horizontal pair halves into their CTM environment without
 *        materializing the rank-6 blob.
 *
 * Environment arguments and leg conventions are identical to
 * contract_reduced_pair_horizontal_density_CTM(). PA/PB carry the leg orders
 * documented by reduced_pair_halves.
 *
 * @return The scalar value of the closed network.
 */
template <class tensor>
typename tensor::value_type contract_reduced_pair_halves_horizontal_density_CTM(
    const tensor& C1, const tensor& C2, const tensor& C3, const tensor& C4,
    const tensor& eT1, const tensor& eT2, const tensor& eT3, const tensor& eT4,
    const tensor& eT5, const tensor& eT6,
    const reduced_pair_halves<tensor>& halves) {
  ::tenes::ScopedTimer scoped_timer("measure/twosite/absorb");
  using mptensor::Axes;
  const tensor left_lower = mptensor::tensordot(
      eT5,
      mptensor::tensordot(C1, mptensor::tensordot(C4, eT6, Axes(1), Axes(0)),
                          Axes(0), Axes(1)),
      Axes(1), Axes(1));
  tensor left =
      mptensor::tensordot(halves.PA, left_lower, Axes(0, 3), Axes(3, 1));
  left = mptensor::tensordot(eT1, left, Axes(0, 2), Axes(3, 0));

  const tensor right_lower = mptensor::tensordot(
      eT4,
      mptensor::tensordot(C2, mptensor::tensordot(C3, eT3, Axes(0), Axes(1)),
                          Axes(1), Axes(1)),
      Axes(0), Axes(1));
  tensor right =
      mptensor::tensordot(halves.PB, right_lower, Axes(2, 3), Axes(3, 1));
  right = mptensor::tensordot(eT2, right, Axes(1, 2), Axes(3, 1));

  tensor joined = mptensor::tensordot(left, right, Axes(1), Axes(1));
  joined = mptensor::transpose(joined, Axes(0, 2, 1, 3));
  return detail::trace_boundary_pairs(joined);
}

/*!
 * @brief Absorb vertical pair halves into their CTM environment without
 *        materializing the rank-6 blob.
 *
 * Environment arguments and leg conventions are identical to
 * contract_reduced_pair_vertical_density_CTM(). PA is the top-site half and
 * PB is the bottom-site half.
 *
 * @return The scalar value of the closed network.
 */
template <class tensor>
typename tensor::value_type contract_reduced_pair_halves_vertical_density_CTM(
    const tensor& C1, const tensor& C2, const tensor& C3, const tensor& C4,
    const tensor& eT1, const tensor& eT2, const tensor& eT3, const tensor& eT4,
    const tensor& eT5, const tensor& eT6,
    const reduced_pair_halves<tensor>& halves) {
  ::tenes::ScopedTimer scoped_timer("measure/twosite/absorb");
  using mptensor::Axes;
  const tensor top_left = mptensor::tensordot(
      eT6,
      mptensor::tensordot(C1, mptensor::tensordot(C2, eT1, Axes(0), Axes(1)),
                          Axes(1), Axes(1)),
      Axes(1), Axes(0));
  tensor top = mptensor::tensordot(halves.PA, top_left, Axes(0, 1), Axes(1, 3));
  top = mptensor::tensordot(eT2, top, Axes(0, 2), Axes(3, 0));

  const tensor bottom_right = mptensor::tensordot(
      eT5,
      mptensor::tensordot(C3, mptensor::tensordot(C4, eT4, Axes(0), Axes(1)),
                          Axes(1), Axes(1)),
      Axes(0), Axes(1));
  tensor bot =
      mptensor::tensordot(halves.PB, bottom_right, Axes(0, 3), Axes(1, 3));
  bot = mptensor::tensordot(eT3, bot, Axes(1, 2), Axes(3, 1));

  tensor joined = mptensor::tensordot(top, bot, Axes(1), Axes(1));
  joined = mptensor::transpose(joined, Axes(0, 2, 1, 3));
  return detail::trace_boundary_pairs(joined);
}

/*!
 * @brief Direction dispatcher for absorbing pair halves into a CTM
 *        environment.
 *
 * The orientation is taken solely from halves.direction; callers cannot pass
 * independent direction metadata that disagrees with the folded tensors.
 */
template <class tensor>
typename tensor::value_type contract_reduced_pair_halves_density_CTM(
    const tensor& C1, const tensor& C2, const tensor& C3, const tensor& C4,
    const tensor& eT1, const tensor& eT2, const tensor& eT3, const tensor& eT4,
    const tensor& eT5, const tensor& eT6,
    const reduced_pair_halves<tensor>& halves) {
  if (halves.direction == reduced_pair_direction::horizontal) {
    return contract_reduced_pair_halves_horizontal_density_CTM(
        C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6, halves);
  }
  return contract_reduced_pair_halves_vertical_density_CTM(
      C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6, halves);
}

}  // namespace tenes::fermion

#endif  // TENES_SRC_FERMION_REDUCED_MEASURE_HPP_
