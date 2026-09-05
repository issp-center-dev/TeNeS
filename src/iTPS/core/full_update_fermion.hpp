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
 *  @brief Fermionic two-site full-update interface.
 */

#ifndef TENES_SRC_ITPS_CORE_FULL_UPDATE_FERMION_HPP_
#define TENES_SRC_ITPS_CORE_FULL_UPDATE_FERMION_HPP_

#include "../../fermion/ftensor.hpp"
#include "../../fermion/reduced.hpp"

namespace tenes::itps {
class PEPS_Parameters;
namespace core {

/*! @brief Apply a fermionic full update to one nearest-neighbor bond.
 *
 *  @tparam tensor Tensor class in mptensor.
 *  @param[in] C1 First corner transfer matrix (see layout below).
 *  @param[in] C2 Second corner transfer matrix.
 *  @param[in] C3 Third corner transfer matrix.
 *  @param[in] C4 Fourth corner transfer matrix.
 *  @param[in] eT1 First direction-ordered edge tensor.
 *  @param[in] eT2 Second direction-ordered edge tensor.
 *  @param[in] eT3 Third direction-ordered edge tensor.
 *  @param[in] eT4 Fourth direction-ordered edge tensor.
 *  @param[in] eT5 Fifth direction-ordered edge tensor.
 *  @param[in] eT6 Sixth direction-ordered edge tensor.
 *  @param[in] Tn1 Center tensor at the left or top end of the bond.
 *  @param[in] Tn2 Center tensor at the right or bottom end of the bond.
 *  @param[in] wrapped_gate Graded, direction-normalized evolution gate.
 *  @param[in] direction Orientation of the updated bond.
 *  @param[in] peps_parameters Full-update hyperparameters.
 *  @param[out] Tn1_new Updated first center tensor.
 *  @param[out] Tn2_new Updated second center tensor.
 *
 *  The driver supplies C1..eT6 in direction-dependent order so the core sees
 *  the corresponding canonical window in both orientations.
 *
 *  @par Horizontal tensor layout
 *  @verbatim
      C1 eT1 eT2 C2
     eT6 Tn1 Tn2 eT3
      C4 eT5 eT4 C3
    @endverbatim
 *  For a vertical bond, Tn1 is above Tn2 and the direction-specific argument
 *  order maps the rotated window onto this canonical ordering.
 */
template <class tensor>
void Full_update_bond_fermion(
    const tensor& C1, const tensor& C2, const tensor& C3, const tensor& C4,
    const tensor& eT1, const tensor& eT2, const tensor& eT3, const tensor& eT4,
    const tensor& eT5, const tensor& eT6,
    const tenes::fermion::ftensor<tensor>& Tn1,
    const tenes::fermion::ftensor<tensor>& Tn2,
    const tenes::fermion::ftensor<tensor>& wrapped_gate,
    tenes::fermion::reduced_pair_direction direction,
    const PEPS_Parameters& peps_parameters,
    tenes::fermion::ftensor<tensor>& Tn1_new,
    tenes::fermion::ftensor<tensor>& Tn2_new);

}  // namespace core
}  // namespace tenes::itps

#endif  // TENES_SRC_ITPS_CORE_FULL_UPDATE_FERMION_HPP_
