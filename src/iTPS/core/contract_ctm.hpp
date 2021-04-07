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

#ifndef TENES_SRC_ITPS_CORE_CONTRACT_CTM_HPP_
#define TENES_SRC_ITPS_CORE_CONTRACT_CTM_HPP_

#include <vector>
#include <cstddef>

namespace tenes {
namespace itps {
namespace core {

/*! @brief contract tensors with CTM
 *
 *  @param[in] C corner transfer matrix
 *  @param[in] eTt top edge tensors
 *  @param[in] eTr right edge tensors
 *  @param[in] eTb bottom edge tensors
 *  @param[in] eTl left edge tensors
 *  @param[in] Tn center tensors
 *  @param[in] op onesite operators
 */
template <class tensor>
typename tensor::value_type Contract(
    const std::vector<const tensor *> &C,
    const std::vector<const tensor *> &eTt,
    const std::vector<const tensor *> &eTr,
    const std::vector<const tensor *> &eTb,
    const std::vector<const tensor *> &eTl,
    const std::vector<std::vector<const tensor *>> &Tn,
    const std::vector<std::vector<const tensor *>> &op);

/*
Declare template functions like

template <class tensor>
typename tensor::value_type
Contract_1x1 (
  const std::vector<const tensor*> &C,
  const std::vector<const tensor*> &eTt,
  const std::vector<const tensor*> &eTr,
  const std::vector<const tensor*> &eTb,
  const std::vector<const tensor*> &eTl,
  const std::vector<std::vector<const tensor*>> &Tn,
  const std::vector<std::vector<const tensor*>> &op
)
*/
#define DECLARE_CONTRACT(NROW, NCOL)                      \
  template <class tensor>                                 \
  typename tensor::value_type Contract_##NROW##x##NCOL(   \
      const std::vector<const tensor *> &C,               \
      const std::vector<const tensor *> &eTt,             \
      const std::vector<const tensor *> &eTr,             \
      const std::vector<const tensor *> &eTb,             \
      const std::vector<const tensor *> &eTl,             \
      const std::vector<std::vector<const tensor *>> &Tn, \
      const std::vector<std::vector<const tensor *>> &op)

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

#undef DECLARE_CONTRACT

template <class tensor>
typename tensor::value_type Contract_one_site(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &Tn1, const tensor &op1);

template <class tensor>
typename tensor::value_type Contract_two_sites_horizontal(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op1, const tensor &op2);

template <class tensor>
typename tensor::value_type Contract_two_sites_vertical(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op1, const tensor &op2);

template <class tensor>
typename tensor::value_type Contract_two_sites_horizontal_op12(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op12);

template <class tensor>
typename tensor::value_type Contract_two_sites_vertical_op12(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op12);

template <class tensor>
typename tensor::value_type Contract_four_sites(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &eT7, const tensor &eT8,
    const tensor &Tn1, const tensor &Tn2, const tensor &Tn3, const tensor &Tn4,
    const tensor &op1, const tensor &op2, const tensor &op3, const tensor &op4);

template <class tensor>
void StartCorrelation(tensor &A, const tensor &C1, const tensor &C4,
                      const tensor &eT1, const tensor &eT3, const tensor &eT4,
                      const tensor &Tn1, const tensor &op);

template <class tensor>
void Transfer(tensor &A, const tensor &eT1, const tensor &eT3,
              const tensor &Tn1);

template <class tensor>
typename tensor::value_type FinishCorrelation(
    const tensor &A, const tensor &C2, const tensor &C3, const tensor &eT1,
    const tensor &eT2, const tensor &eT3, const tensor &Tn1, const tensor &op);

template <class tensor>
void TransferMatrix_MatVec(tensor &inoutvec, const tensor &eT1);

}  // end of namespace core
}  // namespace itps
}  // namespace tenes

#endif  // TENES_SRC_ITPS_CORE_CONTRACT_CTM_HPP_
