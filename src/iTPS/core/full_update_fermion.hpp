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

#ifndef TENES_SRC_ITPS_CORE_FULL_UPDATE_FERMION_HPP_
#define TENES_SRC_ITPS_CORE_FULL_UPDATE_FERMION_HPP_

#include "../../fermion/ftensor.hpp"
#include "../../fermion/reduced.hpp"

namespace tenes::itps {
class PEPS_Parameters;
namespace core {

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
