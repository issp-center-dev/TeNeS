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

#include "../contract_ctm.hpp"

#include "../../../tensor.hpp"


namespace tenes {

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
typename tensor::value_type Contract(
    const std::vector<const tensor *> &C,
    const std::vector<const tensor *> &eTt,
    const std::vector<const tensor *> &eTr,
    const std::vector<const tensor *> &eTb,
    const std::vector<const tensor *> &eTl,
    const std::vector<std::vector<const tensor *>> &Tn,
    const std::vector<std::vector<const tensor *>> &op) {
  const size_t nrow = Tn.size();
  const size_t ncol = Tn[0].size();

#define CALL_CONTRACT(NROW, NCOL)                                     \
  do {                                                                \
    if (nrow == NROW && ncol == NCOL) {                               \
      return Contract_##NROW##x##NCOL(C, eTt, eTr, eTb, eTl, Tn, op); \
    }                                                                 \
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
  ss << "Contract_" << nrow << "_" << ncol << " is not implemented";
  throw std::runtime_error(ss.str());
}

template <class tensor>
typename tensor::value_type Contract_one_site(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &Tn1, const tensor &op1) {
  /*
    ##############################
    # (((((C2*(C3*eT2))*(C1*eT1))*(((C4*eT4)*eT3)*Tn1c))*Tn1)*op1)
    # cpu_cost= 6.04e+10  memory= 3.0207e+08
    # final_bond_order  ()
    ##############################
  */
  return trace(
      tensordot(
          tensordot(
              tensordot(tensordot(C2, tensordot(C3, eT2, Axes(0), Axes(1)),
                                  Axes(1), Axes(1)),
                        tensordot(C1, eT1, Axes(1), Axes(0)), Axes(0), Axes(1)),
              tensordot(tensordot(tensordot(C4, eT4, Axes(1), Axes(0)), eT3,
                                  Axes(0), Axes(1)),
                        conj(Tn1), Axes(2, 5), Axes(0, 3)),
              Axes(0, 2, 3, 5), Axes(2, 5, 0, 4)),
          Tn1, Axes(0, 1, 2, 3), Axes(2, 1, 0, 3)),
      op1, Axes(0, 1), Axes(1, 0));
}

template <class tensor>
typename tensor::value_type Contract_two_sites_horizontal(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op1, const tensor &op2) {
  /*
    ##   |  |
    ##  -T1-T2-
    ##   |  |
    INFO:19234428 (6,5) Finish 10615427/19234428 script=[10, 0, 4, -1, 3, 9, -1,
    -1, 11, 1, 5, -1, 12, 2, 7, -1, 6, -1, 13, 15, -1, -1, -1, -1, 8, -1, -1,
    -1, -1, 14, -1]
    ##############################
    #
    ((Tn1*(((C1*eT1)*(C4*eT6))*(Tn1c*(((C2*eT2)*(Tn2*(((C3*eT4)*eT3)*(Tn2c*op2))))*eT5))))*op1)
    # cpu_cost= 1.204e+11  memory= 3.041e+08
    # final_bond_order  ()
    ##############################
  */
  return trace(
      tensordot(
          Tn1,
          tensordot(
              tensordot(tensordot(C1, eT1, Axes(1), Axes(0)),
                        tensordot(C4, eT6, Axes(1), Axes(0)), Axes(0), Axes(1)),
              tensordot(
                  conj(Tn1),
                  tensordot(
                      tensordot(
                          tensordot(C2, eT2, Axes(0), Axes(1)),
                          tensordot(
                              Tn2,
                              tensordot(
                                  tensordot(
                                      tensordot(C3, eT4, Axes(1), Axes(0)), eT3,
                                      Axes(0), Axes(1)),
                                  tensordot(conj(Tn2), op2, Axes(4), Axes(1)),
                                  Axes(2, 5), Axes(3, 2)),
                              Axes(2, 3, 4), Axes(3, 1, 6)),
                          Axes(0, 2, 3), Axes(3, 1, 5)),
                      eT5, Axes(2), Axes(0)),
                  Axes(2, 3), Axes(2, 5)),
              Axes(0, 2, 3, 5), Axes(3, 1, 5, 0)),
          Axes(0, 1, 2, 3), Axes(1, 0, 3, 4)),
      op1, Axes(0, 1), Axes(0, 1));
}

template <class tensor>
typename tensor::value_type Contract_two_sites_vertical(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op1, const tensor &op2) {
  /*
    ##   |
    ##  -T1-
    ##   |
    ##  -T2-
    ##   |
  INFO:15775477 (3,3) Finish 5129244/15775477 script=[10, 1, 5, -1, 0, 4, -1,
  -1, 11, 2, 6, -1, 12, 3, 8, -1, 7, -1, 13, 15, -1, -1, -1, -1, 9, -1, -1, -1,
  -1, 14, -1]
  ##############################
  #
  ((Tn1*(((C2*eT2)*(C1*eT1))*(Tn1c*(((C3*eT3)*(Tn2*(((C4*eT5)*eT4)*(Tn2c*op2))))*eT6))))*op1)
  # cpu_cost= 1.204e+11  memory= 3.0411e+08
  # final_bond_order  ()
  ##############################
  */
  return trace(
      tensordot(
          Tn1,
          tensordot(
              tensordot(tensordot(C2, eT2, Axes(1), Axes(0)),
                        tensordot(C1, eT1, Axes(1), Axes(0)), Axes(0), Axes(1)),
              tensordot(
                  conj(Tn1),
                  tensordot(
                      tensordot(
                          tensordot(C3, eT3, Axes(0), Axes(1)),
                          tensordot(
                              Tn2,
                              tensordot(
                                  tensordot(
                                      tensordot(C4, eT5, Axes(1), Axes(0)), eT4,
                                      Axes(0), Axes(1)),
                                  tensordot(conj(Tn2), op2, Axes(4), Axes(1)),
                                  Axes(2, 5), Axes(0, 3)),
                              Axes(0, 3, 4), Axes(1, 3, 6)),
                          Axes(0, 2, 3), Axes(3, 1, 5)),
                      eT6, Axes(2), Axes(0)),
                  Axes(0, 3), Axes(5, 2)),
              Axes(0, 2, 3, 5), Axes(3, 1, 5, 0)),
          Axes(0, 1, 2, 3), Axes(4, 1, 0, 3)),
      op1, Axes(0, 1), Axes(0, 1));
}

template <class tensor>
typename tensor::value_type Contract_two_sites_horizontal_op12(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op12) {
  /*
    ##   |  |
    ##  -T1-T2-
    ##   |  |
  */
  ////////////////////////////////////////////////////////////
  // two_sites_horizontal_op12.dat
  ////////////////////////////////////////////////////////////
  // (op12*((eT1*(Tn1*(Tn1c*(eT5*(C1*(C4*eT6))))))*(eT2*(Tn2*(Tn2c*(eT4*(C2*(C3*eT3))))))))
  // cpu_cost= 2.20416e+11  memory= 6.0502e+08
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op12,
      tensordot(
          tensordot(
              eT1,
              tensordot(
                  Tn1,
                  tensordot(
                      conj(Tn1),
                      tensordot(
                          eT5,
                          tensordot(C1, tensordot(C4, eT6, Axes(1), Axes(0)),
                                    Axes(0), Axes(1)),
                          Axes(1), Axes(1)),
                      Axes(0, 3), Axes(5, 2)),
                  Axes(0, 3), Axes(6, 4)),
              Axes(0, 2, 3), Axes(7, 0, 3)),
          tensordot(
              eT2,
              tensordot(
                  Tn2,
                  tensordot(
                      conj(Tn2),
                      tensordot(
                          eT4,
                          tensordot(C2, tensordot(C3, eT3, Axes(0), Axes(1)),
                                    Axes(1), Axes(1)),
                          Axes(0), Axes(1)),
                      Axes(2, 3), Axes(5, 2)),
                  Axes(2, 3), Axes(6, 4)),
              Axes(1, 2, 3), Axes(7, 1, 4)),
          Axes(0, 1, 3, 5), Axes(0, 1, 3, 5)),
      Axes(0, 1, 2, 3), Axes(0, 2, 1, 3));
}

template <class tensor>
typename tensor::value_type Contract_two_sites_vertical_op12(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op12) {
  /*
    ##   |
    ##  -T1-
    ##   |
    ##  -T2-
    ##   |
  */
  ////////////////////////////////////////////////////////////
  // two_sites_vertical_op12.dat
  ////////////////////////////////////////////////////////////
  // (op12*((eT2*(Tn1*(Tn1c*(eT6*(C1*(C2*eT1))))))*(eT3*(Tn2*(Tn2c*(eT5*(C3*(C4*eT4))))))))
  // cpu_cost= 2.20416e+11  memory= 6.0502e+08
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op12,
      tensordot(
          tensordot(
              eT2,
              tensordot(
                  Tn1,
                  tensordot(
                      conj(Tn1),
                      tensordot(
                          eT6,
                          tensordot(C1, tensordot(C2, eT1, Axes(0), Axes(1)),
                                    Axes(1), Axes(1)),
                          Axes(1), Axes(0)),
                      Axes(0, 1), Axes(2, 5)),
                  Axes(0, 1), Axes(4, 6)),
              Axes(0, 2, 3), Axes(7, 0, 3)),
          tensordot(
              eT3,
              tensordot(
                  Tn2,
                  tensordot(
                      conj(Tn2),
                      tensordot(
                          eT5,
                          tensordot(C3, tensordot(C4, eT4, Axes(0), Axes(1)),
                                    Axes(1), Axes(1)),
                          Axes(0), Axes(1)),
                      Axes(0, 3), Axes(2, 5)),
                  Axes(0, 3), Axes(4, 6)),
              Axes(1, 2, 3), Axes(7, 1, 4)),
          Axes(0, 1, 3, 5), Axes(0, 1, 3, 5)),
      Axes(0, 1, 2, 3), Axes(0, 2, 1, 3));
}

template <class tensor>
typename tensor::value_type Contract_four_sites(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &eT7, const tensor &eT8,
    const tensor &Tn1, const tensor &Tn2, const tensor &Tn3, const tensor &Tn4,
    const tensor &op1, const tensor &op2, const tensor &op3,
    const tensor &op4) {
  /*
    INFO:535 (2,1) Finish 317/535 script=[1, 2, -1, 3, -1, 5, -1, 0, 4, -1, -1]
    ##############################
    # ((((C1*eT1)*eT8)*Tn1c)*(op1*Tn1))
    # cpu_cost= 5.01e+10  memory= 3.0004e+08
    # final_bond_order  (e12, e78, tc12, tc41, t12, t41)
    ##############################
  */
  tensor LT = tensordot(
      tensordot(tensordot(tensordot(C1, eT1, Axes(1), Axes(0)), eT8, Axes(0),
                          Axes(1)),
                conj(Tn1), Axes(2, 5), Axes(1, 0)),
      tensordot(op1, Tn1, Axes(0), Axes(4)), Axes(1, 3, 6), Axes(2, 1, 0));
  /*
    INFO:584 (2,1) Finish 368/584 script=[1, 3, -1, 2, -1, 5, -1, 0, 4, -1, -1]
    ##############################
    # ((((C2*eT3)*eT2)*Tn2c)*(op2*Tn2))
    # cpu_cost= 5.01e+10  memory= 3.0004e+08
    # final_bond_order  (e34, e12, tc12, tc23, t12, t23)
    ##############################
  */
  tensor RT = tensordot(
      tensordot(tensordot(tensordot(C2, eT3, Axes(1), Axes(0)), eT2, Axes(0),
                          Axes(1)),
                conj(Tn2), Axes(2, 5), Axes(2, 1)),
      tensordot(op2, Tn2, Axes(0), Axes(4)), Axes(1, 3, 6), Axes(3, 2, 0));
  /*
    INFO:584 (2,1) Finish 368/584 script=[1, 3, -1, 2, -1, 5, -1, 0, 4, -1, -1]
    ##############################
    # ((((C3*eT5)*eT4)*Tn3c)*(op3*Tn3))
    # cpu_cost= 5.01e+10  memory= 3.0004e+08
    # final_bond_order  (e56, e34, tc34, tc23, t34, t23)
    ##############################
  */
  tensor RB = tensordot(
      tensordot(tensordot(tensordot(C3, eT5, Axes(1), Axes(0)), eT4, Axes(0),
                          Axes(1)),
                conj(Tn3), Axes(2, 5), Axes(3, 2)),
      tensordot(op3, Tn3, Axes(0), Axes(4)), Axes(1, 3, 6), Axes(4, 3, 0));
  /*
    INFO:584 (2,1) Finish 368/584 script=[1, 3, -1, 2, -1, 5, -1, 0, 4, -1, -1]
    ##############################
    # ((((C4*eT7)*eT6)*Tn4c)*(op4*Tn4))
    # cpu_cost= 5.01e+10  memory= 3.0004e+08
    # final_bond_order  (e78, e56, tc41, tc34, t41, t34)
    ##############################
  */
  tensor LB = tensordot(
      tensordot(tensordot(tensordot(C4, eT7, Axes(1), Axes(0)), eT6, Axes(0),
                          Axes(1)),
                conj(Tn4), Axes(2, 5), Axes(0, 3)),
      tensordot(op4, Tn4, Axes(0), Axes(4)), Axes(1, 3, 6), Axes(1, 4, 0));
  /*
    INFO:101 (0,0) Finish 2/101 script=[0, 1, 2, 3, -1, -1, -1]
    ##############################
    # (LT*(RT*(RB*LB)))
    # cpu_cost= 2.0001e+12  memory= 5e+08
    # final_bond_order  ()
    ##############################
  */

  return trace(LT,
               tensordot(RT, tensordot(RB, LB, Axes(0, 2, 4), Axes(1, 3, 5)),
                         Axes(0, 3, 5), Axes(0, 1, 2)),
               Axes(0, 1, 2, 3, 4, 5), Axes(0, 3, 1, 4, 2, 5));
}

}  // end of namespace tenes

#undef TENSOR_TYPE
#define TENSOR_TYPE real_tensor
#include "./instantiate_ctm.hpp"
#undef TENSOR_TYPE
#define TENSOR_TYPE complex_tensor
#include "./instantiate_ctm.hpp"
#undef TENSOR_TYPE
