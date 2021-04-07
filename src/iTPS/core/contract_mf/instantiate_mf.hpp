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

/*
 *
 Basic routines independent on unit cell structures.
 Using mptensor libraries
 (Test version)
 2015 Dec.  Tsuyoshi Okubo
*/

#include <vector>

namespace tenes {

/*! @brief contract tensors with mean field environment
 *
 *  @param[in] Tn center tensors
 *  @param[in] op onesite operators
 *
 *  @pre
 *  Tn should absorb MF environment
 */
template TENSOR_TYPE::value_type Contract_MF(
    const std::vector<std::vector<const TENSOR_TYPE *>> &Tn,
    const std::vector<std::vector<const TENSOR_TYPE *>> &op);

template TENSOR_TYPE::value_type Contract_one_site_MF(
    const TENSOR_TYPE &Tn_0_0, const TENSOR_TYPE &op_0_0);

template TENSOR_TYPE::value_type Contract_two_sites_vertical_MF(
    const TENSOR_TYPE &Tn_0_0, const TENSOR_TYPE &Tn_1_0,
    const TENSOR_TYPE &op_0_0, const TENSOR_TYPE &op_1_0);

template TENSOR_TYPE::value_type Contract_two_sites_horizontal_MF(
    const TENSOR_TYPE &Tn_0_0, const TENSOR_TYPE &Tn_0_1,
    const TENSOR_TYPE &op_0_0, const TENSOR_TYPE &op_0_1);

template TENSOR_TYPE::value_type Contract_two_sites_horizontal_op12_MF(
    const TENSOR_TYPE &Tn_0_0, const TENSOR_TYPE &Tn_0_1,
    const TENSOR_TYPE &op12);

template TENSOR_TYPE::value_type Contract_two_sites_vertical_op12_MF(
    const TENSOR_TYPE &Tn_0_0, const TENSOR_TYPE &Tn_1_0,
    const TENSOR_TYPE &op12);

template TENSOR_TYPE::value_type Contract_four_sites_MF(
    const TENSOR_TYPE &Tn_0_0, const TENSOR_TYPE &Tn_0_1,
    const TENSOR_TYPE &Tn_1_0, const TENSOR_TYPE &Tn_1_1,
    const TENSOR_TYPE &op_0_0, const TENSOR_TYPE &op_0_1,
    const TENSOR_TYPE &op_1_0, const TENSOR_TYPE &op_1_1);

}  // end of namespace tenes
