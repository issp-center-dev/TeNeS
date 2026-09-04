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

#ifndef TENES_SRC_ITPS_CORE_FULL_UPDATE_HPP_
#define TENES_SRC_ITPS_CORE_FULL_UPDATE_HPP_

namespace tenes::itps {
class PEPS_Parameters;
namespace core {

template <class tensor>
void prepare_environment(const tensor &Environment_in, const tensor &Theta_in,
                         const PEPS_Parameters &peps_parameters,
                         tensor &Environment_out, tensor &Theta_out,
                         tensor &LR1_inv, tensor &LR2_inv);

template <class tensor>
void als_iterate(const tensor &Environment, const tensor &Theta,
                 const PEPS_Parameters &peps_parameters, tensor &R1,
                 tensor &R2);

template <class tensor>
tensor Create_Environment_two_sites(const tensor &C1, const tensor &C2,
                                    const tensor &C3, const tensor &C4,
                                    const tensor &eT1, const tensor &eT2,
                                    const tensor &eT3, const tensor &eT4,
                                    const tensor &eT5, const tensor &eT6,
                                    const tensor &Q1, const tensor &Q2);

template <class tensor>
void Full_update_bond_horizontal(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op12, const PEPS_Parameters peps_parameters, tensor &Tn1_new,
    tensor &Tn2_new);

/*! @brief full update on a bond
 *
 *  @tparam    tensor tensor class in mptensor
 *  @param[in] C1 corner transfer matrix (see the layout below)
 *  @param[in] C2 corner transfer matrix
 *  @param[in] C3 corner transfer matrix
 *  @param[in] C4 corner transfer matrix
 *  @param[in] eT1 edge tensor
 *  @param[in] eT2 edge tensor
 *  @param[in] eT3 edge tensor
 *  @param[in] eT4 edge tensor
 *  @param[in] eT5 edge tensor
 *  @param[in] eT6 edge tensor
 *  @param[in] Tn1 center tensor on one end of the bond
 *  @param[in] Tn2 center tensor on the other end
 *  @param[in] op12 Imaginary time evolution operator
 *  @param[in] connect1  leg index from Tn1
 *  @param[in] peps_parameters  hyperparameters
 *  @param[out] Tn1_new   new tensor 1
 *  @param[out] Tn2_new   new tensor 2
 *
 * @par Tensors layout
 * @verbatim
  C1 eT1 eT2  C2
 eT6 Tn1 Tn2 eT3
  C4 eT5 eT6  C3
 @endverbatim
 *
 *
 *  @par Reference
 *  H. N. Phien, J. A. Bengua, H. D. Tuan, P. Corboz, and R. Orus, Phys. Rev. B
 @b 92, 035142 (2015)
 *
 */
template <class tensor>
void Full_update_bond(const tensor &C1, const tensor &C2, const tensor &C3,
                      const tensor &C4, const tensor &eT1, const tensor &eT2,
                      const tensor &eT3, const tensor &eT4, const tensor &eT5,
                      const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
                      const tensor &op12, const int connect1,
                      const PEPS_Parameters peps_parameters, tensor &Tn1_new,
                      tensor &Tn2_new);

}  // end of namespace core
}  // namespace tenes::itps

#endif  // TENES_SRC_ITPS_CORE_FULL_UPDATE_HPP_
