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
#include <sstream>

#include "../contract_mf.hpp"
#include "../../../tensor.hpp"

namespace tenes {
namespace itps {
namespace core {

using mptensor::Axes;

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
    const std::vector<std::vector<const tensor *>> &op) {
  const size_t nrow = Tn.size();
  const size_t ncol = Tn[0].size();

#define CALL_CONTRACT(NROW, NCOL)                 \
  do {                                            \
    if (nrow == NROW && ncol == NCOL) {           \
      return Contract_MF_##NROW##x##NCOL(Tn, op); \
    }                                             \
  } while (false)

  CALL_CONTRACT(1, 1);
  CALL_CONTRACT(2, 1);
  CALL_CONTRACT(1, 2);
  CALL_CONTRACT(2, 2);
  CALL_CONTRACT(3, 1);
  CALL_CONTRACT(1, 3);
  CALL_CONTRACT(3, 2);
  CALL_CONTRACT(2, 3);
  CALL_CONTRACT(3, 3);
  CALL_CONTRACT(1, 4);
  CALL_CONTRACT(4, 1);
  CALL_CONTRACT(2, 4);
  CALL_CONTRACT(4, 2);
  CALL_CONTRACT(4, 3);
  CALL_CONTRACT(3, 4);
  CALL_CONTRACT(4, 4);

#undef CALL_CONTRACT

  std::stringstream ss;
  ss << "Contract_MF_" << nrow << "_" << ncol << " is not implemented";
  throw std::runtime_error(ss.str());
}

template <class tensor>
typename tensor::value_type Contract_one_site_MF(const tensor &Tn_0_0,
                                                 const tensor &op_0_0) {
  ////////////////////////////////////////////////////////////
  // (op_0_0*(Tn_0_0*conj(Tn_0_0)))
  // cpu_cost= 262208  memory= 65664
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op_0_0,
      tensordot(Tn_0_0, conj(Tn_0_0), Axes(0, 1, 2, 3), Axes(0, 1, 2, 3)),
      Axes(0, 1), Axes(0, 1));
}

template <class tensor>
typename tensor::value_type Contract_two_sites_vertical_MF(
    const tensor &Tn_0_0, const tensor &Tn_1_0, const tensor &op_0_0,
    const tensor &op_1_0) {
  ////////////////////////////////////////////////////////////
  // (op_0_0*(Tn_0_0*(conj(Tn_0_0)*(Tn_1_0*(conj(Tn_1_0)*op_1_0)))))
  // cpu_cost= 1.04864e+06  memory= 163968
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op_0_0,
      tensordot(
          Tn_0_0,
          tensordot(conj(Tn_0_0),
                    tensordot(Tn_1_0,
                              tensordot(conj(Tn_1_0), op_1_0, Axes(4), Axes(1)),
                              Axes(0, 2, 3, 4), Axes(0, 2, 3, 4)),
                    Axes(3), Axes(1)),
          Axes(0, 1, 2, 3), Axes(0, 1, 2, 4)),
      Axes(0, 1), Axes(0, 1));
}

template <class tensor>
typename tensor::value_type Contract_two_sites_horizontal_MF(
    const tensor &Tn_0_0, const tensor &Tn_0_1, const tensor &op_0_0,
    const tensor &op_0_1) {
  ////////////////////////////////////////////////////////////
  // (op_0_0*(Tn_0_0*(conj(Tn_0_0)*(Tn_0_1*(conj(Tn_0_1)*op_0_1)))))
  // cpu_cost= 1.04864e+06  memory= 163968
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op_0_0,
      tensordot(
          Tn_0_0,
          tensordot(conj(Tn_0_0),
                    tensordot(Tn_0_1,
                              tensordot(conj(Tn_0_1), op_0_1, Axes(4), Axes(1)),
                              Axes(1, 2, 3, 4), Axes(1, 2, 3, 4)),
                    Axes(2), Axes(1)),
          Axes(0, 1, 2, 3), Axes(0, 1, 4, 2)),
      Axes(0, 1), Axes(0, 1));
}

template <class tensor>
typename tensor::value_type Contract_two_sites_horizontal_op12_MF(
    const tensor &Tn_0_0, const tensor &Tn_0_1, const tensor &op12) {
  ////////////////////////////////////////////////////////////
  // (op12*((Tn_0_0*conj(Tn_0_0))*(Tn_0_1*conj(Tn_0_1))))
  // cpu_cost= 4.46054e+06  memory= 139264
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op12,
      tensordot(tensordot(Tn_0_0, conj(Tn_0_0), Axes(0, 1, 3), Axes(0, 1, 3)),
                tensordot(Tn_0_1, conj(Tn_0_1), Axes(1, 2, 3), Axes(1, 2, 3)),
                Axes(0, 2), Axes(0, 2)),
      Axes(0, 1, 2, 3), Axes(0, 2, 1, 3));
}

template <class tensor>
typename tensor::value_type Contract_two_sites_vertical_op12_MF(
    const tensor &Tn_0_0, const tensor &Tn_1_0, const tensor &op12) {
  ////////////////////////////////////////////////////////////
  // (op12*((Tn_0_0*conj(Tn_0_0))*(Tn_1_0*conj(Tn_1_0))))
  // cpu_cost= 4.46054e+06  memory= 139264
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op12,
      tensordot(tensordot(Tn_0_0, conj(Tn_0_0), Axes(0, 1, 2), Axes(0, 1, 2)),
                tensordot(Tn_1_0, conj(Tn_1_0), Axes(0, 2, 3), Axes(0, 2, 3)),
                Axes(0, 2), Axes(0, 2)),
      Axes(0, 1, 2, 3), Axes(0, 2, 1, 3));
}

template <class tensor>
typename tensor::value_type Contract_four_sites_MF(
    const tensor &Tn_0_0, const tensor &Tn_0_1, const tensor &Tn_1_0,
    const tensor &Tn_1_1, const tensor &op_0_0, const tensor &op_0_1,
    const tensor &op_1_0, const tensor &op_1_1) {
  ////////////////////////////////////////////////////////////
  // (op_0_0*(Tn_0_0*(conj(Tn_0_0)*((Tn_0_1*(conj(Tn_0_1)*op_0_1))*((Tn_1_0*(conj(Tn_1_0)*op_1_0))*(Tn_1_1*(conj(Tn_1_1)*op_1_1)))))))
  // cpu_cost= 9.96154e+06  memory= 295168
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op_0_0,
      tensordot(
          Tn_0_0,
          tensordot(
              conj(Tn_0_0),
              tensordot(
                  tensordot(Tn_0_1,
                            tensordot(conj(Tn_0_1), op_0_1, Axes(4), Axes(1)),
                            Axes(1, 2, 4), Axes(1, 2, 4)),
                  tensordot(tensordot(Tn_1_0,
                                      tensordot(conj(Tn_1_0), op_1_0, Axes(4),
                                                Axes(1)),
                                      Axes(0, 3, 4), Axes(0, 3, 4)),
                            tensordot(Tn_1_1,
                                      tensordot(conj(Tn_1_1), op_1_1, Axes(4),
                                                Axes(1)),
                                      Axes(2, 3, 4), Axes(2, 3, 4)),
                            Axes(1, 3), Axes(0, 2)),
                  Axes(1, 3), Axes(2, 3)),
              Axes(2, 3), Axes(1, 3)),
          Axes(0, 1, 2, 3), Axes(0, 1, 3, 4)),
      Axes(0, 1), Axes(0, 1));
}

}  // namespace core
}  // namespace itps
}  // namespace tenes

#undef TENSOR_TYPE
#define TENSOR_TYPE real_tensor
#include "instantiate_mf.hpp"
#undef TENSOR_TYPE
#define TENSOR_TYPE complex_tensor
#include "instantiate_mf.hpp"
#undef TENSOR_TYPE
