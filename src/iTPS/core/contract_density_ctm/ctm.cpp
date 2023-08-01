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

#include "../contract_density_ctm.hpp"

#include "../../../tensor.hpp"

namespace tenes {
namespace itps {
namespace core {

using mptensor::Axes;

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
typename tensor::value_type Contract_density_CTM(
    const std::vector<const tensor *> &C,
    const std::vector<const tensor *> &eTt,
    const std::vector<const tensor *> &eTr,
    const std::vector<const tensor *> &eTb,
    const std::vector<const tensor *> &eTl,
    const std::vector<std::vector<const tensor *>> &Tn,
    const std::vector<std::vector<const tensor *>> &op) {
  const size_t nrow = Tn.size();
  const size_t ncol = Tn[0].size();

#define CALL_CONTRACT(NROW, NCOL)                                              \
  do {                                                                         \
    if (nrow == NROW && ncol == NCOL) {                                        \
      return Contract_##NROW##x##NCOL##_density_CTM(C, eTt, eTr, eTb, eTl, Tn, \
                                                    op);                       \
    }                                                                          \
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
  ss << "Contract_density_" << nrow << "_" << ncol << " is not implemented";
  throw std::runtime_error(ss.str());
}

template <class tensor>
typename tensor::value_type Contract_one_site_density_CTM(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &Tn1, const tensor &op1) {
  ////////////////////////////////////////////////////////////
  // contract_one_site_density.dat
  ////////////////////////////////////////////////////////////
  // (op1*(Tn1*((eT1*(C1*eT4))*((C2*eT2)*(C3*(C4*eT3))))))
  // cpu_cost= 5.96e+06  memory= 130004
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op1,
      tensordot(
          Tn1,
          tensordot(
              tensordot(eT1, tensordot(C1, eT4, Axes(0), Axes(1)), Axes(0),
                        Axes(0)),
              tensordot(tensordot(C2, eT2, Axes(1), Axes(0)),
                        tensordot(C3, tensordot(C4, eT3, Axes(0), Axes(1)),
                                  Axes(1), Axes(1)),
                        Axes(1), Axes(0)),
              Axes(0, 2), Axes(0, 2)),
          Axes(0, 1, 2, 3), Axes(1, 0, 2, 3)),
      Axes(0, 1), Axes(0, 1));
}

template <class tensor>
typename tensor::value_type Contract_two_sites_horizontal_density_CTM(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op1, const tensor &op2) {
  ////////////////////////////////////////////////////////////
  // contract_two_site_horizontal_density.dat
  ////////////////////////////////////////////////////////////
  // (op1*(Tn1*((eT1*(C1*eT6))*((C4*eT5)*(eT2*((Tn2*op2)*(eT4*(C2*(C3*eT3)))))))))
  // cpu_cost= 1.16e+07  memory= 178408
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op1,
      tensordot(
          Tn1,
          tensordot(
              tensordot(eT1, tensordot(C1, eT6, Axes(0), Axes(1)), Axes(0),
                        Axes(0)),
              tensordot(
                  tensordot(C4, eT5, Axes(0), Axes(1)),
                  tensordot(
                      eT2,
                      tensordot(tensordot(Tn2, op2, Axes(4, 5), Axes(0, 1)),
                                tensordot(eT4,
                                          tensordot(C2,
                                                    tensordot(C3, eT3, Axes(0),
                                                              Axes(1)),
                                                    Axes(1), Axes(1)),
                                          Axes(0), Axes(1)),
                                Axes(2, 3), Axes(3, 1)),
                      Axes(1, 2), Axes(3, 1)),
                  Axes(1), Axes(2)),
              Axes(0, 2), Axes(2, 0)),
          Axes(0, 1, 2, 3), Axes(1, 0, 3, 2)),
      Axes(0, 1), Axes(0, 1));
}

template <class tensor>
typename tensor::value_type Contract_two_sites_vertical_density_CTM(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op1, const tensor &op2) {
  ////////////////////////////////////////////////////////////
  // contract_two_site_vertical_density.dat
  ////////////////////////////////////////////////////////////
  // (op1*(Tn1*((eT1*(C1*eT6))*((C2*eT2)*(eT3*((Tn2*op2)*(eT5*(C3*(C4*eT4)))))))))
  // cpu_cost= 1.16e+07  memory= 178408
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op1,
      tensordot(
          Tn1,
          tensordot(
              tensordot(eT1, tensordot(C1, eT6, Axes(0), Axes(1)), Axes(0),
                        Axes(0)),
              tensordot(
                  tensordot(C2, eT2, Axes(1), Axes(0)),
                  tensordot(
                      eT3,
                      tensordot(tensordot(Tn2, op2, Axes(4, 5), Axes(0, 1)),
                                tensordot(eT5,
                                          tensordot(C3,
                                                    tensordot(C4, eT4, Axes(0),
                                                              Axes(1)),
                                                    Axes(1), Axes(1)),
                                          Axes(0), Axes(1)),
                                Axes(0, 3), Axes(1, 3)),
                      Axes(1, 2), Axes(3, 1)),
                  Axes(1), Axes(0)),
              Axes(0, 2), Axes(0, 3)),
          Axes(0, 1, 2, 3), Axes(1, 0, 2, 3)),
      Axes(0, 1), Axes(0, 1));
}

template <class tensor>
typename tensor::value_type Contract_two_sites_horizontal_op12_density_CTM(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op12) {
  ////////////////////////////////////////////////////////////
  // contract_two_site_horizontal_op12_density.dat
  ////////////////////////////////////////////////////////////
  // (op12*((eT1*(Tn1*(eT5*(C1*(C4*eT6)))))*(eT2*(Tn2*(eT4*(C2*(C3*eT3)))))))
  // cpu_cost= 4.0384e+07  memory= 296816
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op12,
      tensordot(
          tensordot(
              eT1,
              tensordot(
                  Tn1,
                  tensordot(eT5,
                            tensordot(C1, tensordot(C4, eT6, Axes(1), Axes(0)),
                                      Axes(0), Axes(1)),
                            Axes(1), Axes(1)),
                  Axes(0, 3), Axes(3, 1)),
              Axes(0, 2), Axes(5, 0)),
          tensordot(
              eT2,
              tensordot(
                  Tn2,
                  tensordot(eT4,
                            tensordot(C2, tensordot(C3, eT3, Axes(0), Axes(1)),
                                      Axes(1), Axes(1)),
                            Axes(0), Axes(1)),
                  Axes(2, 3), Axes(3, 1)),
              Axes(1, 2), Axes(5, 1)),
          Axes(0, 1, 4), Axes(0, 1, 4)),
      Axes(0, 1, 2, 3), Axes(0, 2, 1, 3));
}

template <class tensor>
typename tensor::value_type Contract_two_sites_vertical_op12_density_CTM(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op12) {
  ////////////////////////////////////////////////////////////
  // contract_two_site_vertical_op12_density.dat
  ////////////////////////////////////////////////////////////
  // (op12*((eT2*(Tn1*(eT6*(C1*(C2*eT1)))))*(eT3*(Tn2*(eT5*(C3*(C4*eT4)))))))
  // cpu_cost= 4.0384e+07  memory= 296816
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op12,
      tensordot(
          tensordot(
              eT2,
              tensordot(
                  Tn1,
                  tensordot(eT6,
                            tensordot(C1, tensordot(C2, eT1, Axes(0), Axes(1)),
                                      Axes(1), Axes(1)),
                            Axes(1), Axes(0)),
                  Axes(0, 1), Axes(1, 3)),
              Axes(0, 2), Axes(5, 0)),
          tensordot(
              eT3,
              tensordot(
                  Tn2,
                  tensordot(eT5,
                            tensordot(C3, tensordot(C4, eT4, Axes(0), Axes(1)),
                                      Axes(1), Axes(1)),
                            Axes(0), Axes(1)),
                  Axes(0, 3), Axes(1, 3)),
              Axes(1, 2), Axes(5, 1)),
          Axes(0, 1, 4), Axes(0, 1, 4)),
      Axes(0, 1, 2, 3), Axes(0, 2, 1, 3));
}

}  // namespace core
}  // namespace itps
}  // namespace tenes

#undef TENSOR_TYPE
#define TENSOR_TYPE real_tensor
#include "./instantiate_ctm.hpp"
#undef TENSOR_TYPE
#define TENSOR_TYPE complex_tensor
#include "./instantiate_ctm.hpp"
#undef TENSOR_TYPE
