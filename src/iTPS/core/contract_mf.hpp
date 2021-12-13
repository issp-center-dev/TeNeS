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

#ifndef TENES_SRC_ITPS_CORE_CONTRACT_MF_HPP_
#define TENES_SRC_ITPS_CORE_CONTRACT_MF_HPP_

#include <vector>
#include <cstddef>

namespace tenes {
namespace itps {
namespace core {

/*! @brief contract tensors with mean field environment
 *
 *  @param[in] Tn center tensors
 *  @param[in] op onesite operators
 *
 *  @pre
 *  Tn should absorb MF environment
 */
template <class tensor>
typename tensor::value_type Contract_MF(
    const std::vector<std::vector<const tensor *>> &Tn,
    const std::vector<std::vector<const tensor *>> &op);

template <class tensor>
typename tensor::value_type Contract_one_site_MF(const tensor &Tn_0_0,
                                                 const tensor &op_0_0);

template <class tensor>
typename tensor::value_type Contract_two_sites_vertical_MF(
    const tensor &Tn_0_0, const tensor &Tn_1_0, const tensor &op_0_0,
    const tensor &op_1_0);

template <class tensor>
typename tensor::value_type Contract_two_sites_horizontal_MF(
    const tensor &Tn_0_0, const tensor &Tn_0_1, const tensor &op_0_0,
    const tensor &op_0_1);

template <class tensor>
typename tensor::value_type Contract_two_sites_horizontal_op12_MF(
    const tensor &Tn_0_0, const tensor &Tn_0_1, const tensor &op12);

template <class tensor>
typename tensor::value_type Contract_two_sites_vertical_op12_MF(
    const tensor &Tn_0_0, const tensor &Tn_1_0, const tensor &op12);

template <class tensor>
typename tensor::value_type Contract_four_sites_MF(
    const tensor &Tn_0_0, const tensor &Tn_0_1, const tensor &Tn_1_0,
    const tensor &Tn_1_1, const tensor &op_0_0, const tensor &op_0_1,
    const tensor &op_1_0, const tensor &op_1_1);

/*
Declare template functions like

template <class tensor>
typename tensor::value_type
Contract_MF_1x1 (
  const std::vector<std::vector<const tensor*>> &Tn,
  const std::vector<std::vector<const tensor*>> &op
)
*/

#define DECLARE_CONTRACT(NROW, NCOL)                       \
  template <class tensor>                                  \
  typename tensor::value_type Contract_MF_##NROW##x##NCOL( \
      const std::vector<std::vector<const tensor *>> &Tn,  \
      const std::vector<std::vector<const tensor *>> &op)

//! @cond
DECLARE_CONTRACT(1, 1);
DECLARE_CONTRACT(2, 1);
DECLARE_CONTRACT(3, 1);
DECLARE_CONTRACT(4, 1);
DECLARE_CONTRACT(1, 2);
DECLARE_CONTRACT(2, 2);
DECLARE_CONTRACT(3, 2);
DECLARE_CONTRACT(4, 2);
DECLARE_CONTRACT(1, 3);
DECLARE_CONTRACT(2, 3);
DECLARE_CONTRACT(3, 3);
DECLARE_CONTRACT(4, 3);
DECLARE_CONTRACT(1, 4);
DECLARE_CONTRACT(2, 4);
DECLARE_CONTRACT(3, 4);
DECLARE_CONTRACT(4, 4);
//! @endcond

#undef DECLARE_CONTRACT

/* Lambda tensors should be absorbed before entering */
template <class tensor>
void StartCorrelation_MF(tensor &A, const tensor &Tn, const tensor &op,
                         size_t direction);

/* Lambda tensors should be absorbed before entering */
template <class tensor>
void Transfer_MF(tensor &A, const tensor &Tn, size_t direction);

/* Lambda tensors should be absorbed before entering */
template <class tensor>
typename tensor::value_type FinishCorrelation_MF(const tensor &A,
                                                 const tensor &Tn,
                                                 const tensor &op,
                                                 size_t direction);

}  // end of namespace core
}  // namespace itps
}  // namespace tenes

#endif  // TENES_SRC_ITPS_CORE_CONTRACT_MF_HPP_
