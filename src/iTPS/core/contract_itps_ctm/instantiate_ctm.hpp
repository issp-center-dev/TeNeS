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

//! @cond

namespace tenes {
namespace itps {
namespace core {

template TENSOR_TYPE::value_type Contract_iTPS_CTM(
    const std::vector<const TENSOR_TYPE *> &C,
    const std::vector<const TENSOR_TYPE *> &eTt,
    const std::vector<const TENSOR_TYPE *> &eTr,
    const std::vector<const TENSOR_TYPE *> &eTb,
    const std::vector<const TENSOR_TYPE *> &eTl,
    const std::vector<std::vector<const TENSOR_TYPE *>> &Tn,
    const std::vector<std::vector<const TENSOR_TYPE *>> &op);

template TENSOR_TYPE::value_type Contract_one_site_iTPS_CTM(
    const TENSOR_TYPE &C1, const TENSOR_TYPE &C2, const TENSOR_TYPE &C3,
    const TENSOR_TYPE &C4, const TENSOR_TYPE &eT1, const TENSOR_TYPE &eT2,
    const TENSOR_TYPE &eT3, const TENSOR_TYPE &eT4, const TENSOR_TYPE &Tn1,
    const TENSOR_TYPE &op1);

template TENSOR_TYPE::value_type Contract_two_sites_horizontal_iTPS_CTM(
    const TENSOR_TYPE &C1, const TENSOR_TYPE &C2, const TENSOR_TYPE &C3,
    const TENSOR_TYPE &C4, const TENSOR_TYPE &eT1, const TENSOR_TYPE &eT2,
    const TENSOR_TYPE &eT3, const TENSOR_TYPE &eT4, const TENSOR_TYPE &eT5,
    const TENSOR_TYPE &eT6, const TENSOR_TYPE &Tn1, const TENSOR_TYPE &Tn2,
    const TENSOR_TYPE &op1, const TENSOR_TYPE &op2);

template TENSOR_TYPE::value_type Contract_two_sites_vertical_iTPS_CTM(
    const TENSOR_TYPE &C1, const TENSOR_TYPE &C2, const TENSOR_TYPE &C3,
    const TENSOR_TYPE &C4, const TENSOR_TYPE &eT1, const TENSOR_TYPE &eT2,
    const TENSOR_TYPE &eT3, const TENSOR_TYPE &eT4, const TENSOR_TYPE &eT5,
    const TENSOR_TYPE &eT6, const TENSOR_TYPE &Tn1, const TENSOR_TYPE &Tn2,
    const TENSOR_TYPE &op1, const TENSOR_TYPE &op2);

template TENSOR_TYPE::value_type Contract_two_sites_horizontal_op12_iTPS_CTM(
    const TENSOR_TYPE &C1, const TENSOR_TYPE &C2, const TENSOR_TYPE &C3,
    const TENSOR_TYPE &C4, const TENSOR_TYPE &eT1, const TENSOR_TYPE &eT2,
    const TENSOR_TYPE &eT3, const TENSOR_TYPE &eT4, const TENSOR_TYPE &eT5,
    const TENSOR_TYPE &eT6, const TENSOR_TYPE &Tn1, const TENSOR_TYPE &Tn2,
    const TENSOR_TYPE &op12);

template TENSOR_TYPE::value_type Contract_two_sites_vertical_op12_iTPS_CTM(
    const TENSOR_TYPE &C1, const TENSOR_TYPE &C2, const TENSOR_TYPE &C3,
    const TENSOR_TYPE &C4, const TENSOR_TYPE &eT1, const TENSOR_TYPE &eT2,
    const TENSOR_TYPE &eT3, const TENSOR_TYPE &eT4, const TENSOR_TYPE &eT5,
    const TENSOR_TYPE &eT6, const TENSOR_TYPE &Tn1, const TENSOR_TYPE &Tn2,
    const TENSOR_TYPE &op12);

template TENSOR_TYPE::value_type Contract_four_sites_iTPS_CTM(
    const TENSOR_TYPE &C1, const TENSOR_TYPE &C2, const TENSOR_TYPE &C3,
    const TENSOR_TYPE &C4, const TENSOR_TYPE &eT1, const TENSOR_TYPE &eT2,
    const TENSOR_TYPE &eT3, const TENSOR_TYPE &eT4, const TENSOR_TYPE &eT5,
    const TENSOR_TYPE &eT6, const TENSOR_TYPE &eT7, const TENSOR_TYPE &eT8,
    const TENSOR_TYPE &Tn1, const TENSOR_TYPE &Tn2, const TENSOR_TYPE &Tn3,
    const TENSOR_TYPE &Tn4, const TENSOR_TYPE &op1, const TENSOR_TYPE &op2,
    const TENSOR_TYPE &op3, const TENSOR_TYPE &op4);
}  // namespace core
}  // namespace itps
}  // namespace tenes

//! @endcond
