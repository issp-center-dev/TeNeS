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

#ifndef _PEPS_BASICS_HPP_
#define _PEPS_BASICS_HPP_


#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <mptensor/complex.hpp>
#include <mptensor/rsvd.hpp>
#include <mptensor/tensor.hpp>
#include <numeric>
#include <vector>

#include "exception.hpp"

#include "printlevel.hpp"
#include "PEPS_Parameters.hpp"
#include "mpi.hpp"

#include "PEPS_Contract_impl.hpp"
#include "PEPS_Contract_MF_impl.hpp"

namespace tenes {

using namespace mptensor;
// Contractions

template <class tensor>
typename tensor::value_type
Contract(const std::vector<const tensor*> &C,
         const std::vector<const tensor*> &eTt,
         const std::vector<const tensor*> &eTr,
         const std::vector<const tensor*> &eTb,
         const std::vector<const tensor*> &eTl,
         const std::vector<std::vector<const tensor*>> &Tn,
         const std::vector<std::vector<const tensor*>> &op
         ){
  const size_t nrow = Tn.size();
  const size_t ncol = Tn[0].size();

#define CALL_CONTRACT(NROW, NCOL) \
  do{\
    if(nrow == NROW && ncol == NCOL){\
      return Contract_ ## NROW ## x ## NCOL (C, eTt, eTr, eTb, eTl, Tn, op);\
    }\
  }while(false)

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
typename tensor::value_type
Contract_MF(const std::vector<std::vector<const tensor*>> &Tn,
            const std::vector<std::vector<const tensor*>> &op
            ){
  const size_t nrow = Tn.size();
  const size_t ncol = Tn[0].size();

#define CALL_CONTRACT_MF(NROW, NCOL) \
  do{\
    if(nrow == NROW && ncol == NCOL){\
      return Contract_MF_ ## NROW ## x ## NCOL (Tn, op);\
    }\
  }while(false)

  CALL_CONTRACT_MF(1, 1);
  CALL_CONTRACT_MF(2, 1);
  CALL_CONTRACT_MF(1, 2);
  CALL_CONTRACT_MF(2, 2);
  CALL_CONTRACT_MF(3, 1);
  CALL_CONTRACT_MF(1, 3);
  CALL_CONTRACT_MF(3, 2);
  CALL_CONTRACT_MF(2, 3);
  CALL_CONTRACT_MF(3, 3);
  CALL_CONTRACT_MF(1, 4);
  CALL_CONTRACT_MF(4, 1);
  CALL_CONTRACT_MF(2, 4);
  CALL_CONTRACT_MF(4, 2);
  CALL_CONTRACT_MF(4, 3);
  CALL_CONTRACT_MF(3, 4);
  CALL_CONTRACT_MF(4, 4);

#undef CALL_CONTRACT_MF

  std::stringstream ss;
  ss << "Contract_MF_" << nrow << "_" << ncol << " is not implemented";
  throw std::runtime_error(ss.str());
}


template <template <typename> class Matrix, typename C>
C Contract_one_site(const Tensor<Matrix, C> &C1, const Tensor<Matrix, C> &C2,
                    const Tensor<Matrix, C> &C3, const Tensor<Matrix, C> &C4,
                    const Tensor<Matrix, C> &eT1, const Tensor<Matrix, C> &eT2,
                    const Tensor<Matrix, C> &eT3, const Tensor<Matrix, C> &eT4,
                    const Tensor<Matrix, C> &Tn1,
                    const Tensor<Matrix, C> &op1) {
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
template <template <typename> class Matrix, typename C>
C Contract_two_sites_horizontal(
    const Tensor<Matrix, C> &C1, const Tensor<Matrix, C> &C2,
    const Tensor<Matrix, C> &C3, const Tensor<Matrix, C> &C4,
    const Tensor<Matrix, C> &eT1, const Tensor<Matrix, C> &eT2,
    const Tensor<Matrix, C> &eT3, const Tensor<Matrix, C> &eT4,
    const Tensor<Matrix, C> &eT5, const Tensor<Matrix, C> &eT6,
    const Tensor<Matrix, C> &Tn1, const Tensor<Matrix, C> &Tn2,
    const Tensor<Matrix, C> &op1, const Tensor<Matrix, C> &op2) {
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

template <template <typename> class Matrix, typename C>
C Contract_two_sites_vertical(
    const Tensor<Matrix, C> &C1, const Tensor<Matrix, C> &C2,
    const Tensor<Matrix, C> &C3, const Tensor<Matrix, C> &C4,
    const Tensor<Matrix, C> &eT1, const Tensor<Matrix, C> &eT2,
    const Tensor<Matrix, C> &eT3, const Tensor<Matrix, C> &eT4,
    const Tensor<Matrix, C> &eT5, const Tensor<Matrix, C> &eT6,
    const Tensor<Matrix, C> &Tn1, const Tensor<Matrix, C> &Tn2,
    const Tensor<Matrix, C> &op1, const Tensor<Matrix, C> &op2) {
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
template <template <typename> class Matrix, typename C>
C Contract_two_sites_horizontal_op12(
    const Tensor<Matrix, C> &C1, const Tensor<Matrix, C> &C2,
    const Tensor<Matrix, C> &C3, const Tensor<Matrix, C> &C4,
    const Tensor<Matrix, C> &eT1, const Tensor<Matrix, C> &eT2,
    const Tensor<Matrix, C> &eT3, const Tensor<Matrix, C> &eT4,
    const Tensor<Matrix, C> &eT5, const Tensor<Matrix, C> &eT6,
    const Tensor<Matrix, C> &Tn1, const Tensor<Matrix, C> &Tn2,
    const Tensor<Matrix, C> &op12) {
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

template <template <typename> class Matrix, typename C>
C Contract_two_sites_vertical_op12(
    const Tensor<Matrix, C> &C1, const Tensor<Matrix, C> &C2,
    const Tensor<Matrix, C> &C3, const Tensor<Matrix, C> &C4,
    const Tensor<Matrix, C> &eT1, const Tensor<Matrix, C> &eT2,
    const Tensor<Matrix, C> &eT3, const Tensor<Matrix, C> &eT4,
    const Tensor<Matrix, C> &eT5, const Tensor<Matrix, C> &eT6,
    const Tensor<Matrix, C> &Tn1, const Tensor<Matrix, C> &Tn2,
    const Tensor<Matrix, C> &op12) {
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
template <template <typename> class Matrix, typename C>
C Contract_four_sites(
    const Tensor<Matrix, C> &C1, const Tensor<Matrix, C> &C2,
    const Tensor<Matrix, C> &C3, const Tensor<Matrix, C> &C4,
    const Tensor<Matrix, C> &eT1, const Tensor<Matrix, C> &eT2,
    const Tensor<Matrix, C> &eT3, const Tensor<Matrix, C> &eT4,
    const Tensor<Matrix, C> &eT5, const Tensor<Matrix, C> &eT6,
    const Tensor<Matrix, C> &eT7, const Tensor<Matrix, C> &eT8,
    const Tensor<Matrix, C> &Tn1, const Tensor<Matrix, C> &Tn2,
    const Tensor<Matrix, C> &Tn3, const Tensor<Matrix, C> &Tn4,
    const Tensor<Matrix, C> &op1, const Tensor<Matrix, C> &op2,
    const Tensor<Matrix, C> &op3, const Tensor<Matrix, C> &op4) {
  /*
    INFO:535 (2,1) Finish 317/535 script=[1, 2, -1, 3, -1, 5, -1, 0, 4, -1, -1]
    ##############################
    # ((((C1*eT1)*eT8)*Tn1c)*(op1*Tn1))
    # cpu_cost= 5.01e+10  memory= 3.0004e+08
    # final_bond_order  (e12, e78, tc12, tc41, t12, t41)
    ##############################
  */
  Tensor<Matrix, C> LT = tensordot(
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
  Tensor<Matrix, C> RT = tensordot(
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
  Tensor<Matrix, C> RB = tensordot(
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
  Tensor<Matrix, C> LB = tensordot(
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

template <class tensor>
typename tensor::value_type
Contract_one_site_MF(
  const tensor &Tn_0_0,
  const tensor &op_0_0
)
{
  ////////////////////////////////////////////////////////////
  // (op_0_0*(Tn_0_0*conj(Tn_0_0)))
  // cpu_cost= 262208  memory= 65664
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return
  trace(
    op_0_0, tensordot(
      Tn_0_0, conj(Tn_0_0), Axes(0, 1, 2, 3), Axes(0, 1, 2, 3)
    ), Axes(0, 1), Axes(0, 1)
  )
  ;
}

template <class tensor>
typename tensor::value_type
Contract_two_sites_vertical_MF(
  const tensor &Tn_0_0,
  const tensor &Tn_1_0,
  const tensor &op_0_0,
  const tensor &op_1_0
)
{
  ////////////////////////////////////////////////////////////
  // (op_0_0*(Tn_0_0*(conj(Tn_0_0)*(Tn_1_0*(conj(Tn_1_0)*op_1_0)))))
  // cpu_cost= 1.04864e+06  memory= 163968
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return
  trace(
    op_0_0, tensordot(
      Tn_0_0, tensordot(
        conj(Tn_0_0), tensordot(
          Tn_1_0, tensordot(
            conj(Tn_1_0), op_1_0, Axes(4), Axes(1)
          ), Axes(0, 2, 3, 4), Axes(0, 2, 3, 4)
        ), Axes(3), Axes(1)
      ), Axes(0, 1, 2, 3), Axes(0, 1, 2, 4)
    ), Axes(0, 1), Axes(0, 1)
  )
  ;
}

template <class tensor>
typename tensor::value_type
Contract_two_sites_horizontal_MF(
  const tensor &Tn_0_0,
  const tensor &Tn_0_1,
  const tensor &op_0_0,
  const tensor &op_0_1
)
{
  ////////////////////////////////////////////////////////////
  // (op_0_0*(Tn_0_0*(conj(Tn_0_0)*(Tn_0_1*(conj(Tn_0_1)*op_0_1)))))
  // cpu_cost= 1.04864e+06  memory= 163968
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return
  trace(
    op_0_0, tensordot(
      Tn_0_0, tensordot(
        conj(Tn_0_0), tensordot(
          Tn_0_1, tensordot(
            conj(Tn_0_1), op_0_1, Axes(4), Axes(1)
          ), Axes(1, 2, 3, 4), Axes(1, 2, 3, 4)
        ), Axes(2), Axes(1)
      ), Axes(0, 1, 2, 3), Axes(0, 1, 4, 2)
    ), Axes(0, 1), Axes(0, 1)
  )
  ;
}

template <class tensor>
typename tensor::value_type
Contract_two_sites_vertical_op12_MF(
  const tensor &Tn_0_0,
  const tensor &Tn_0_1,
  const tensor &op12
)
{
  ////////////////////////////////////////////////////////////
  // (op12*((Tn_0_0*conj(Tn_0_0))*(Tn_0_1*conj(Tn_0_1))))
  // cpu_cost= 4.46054e+06  memory= 139264
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return
  trace(
    op12, tensordot(
      tensordot(
        Tn_0_0, conj(Tn_0_0), Axes(0, 1, 3), Axes(0, 1, 3)
      ), tensordot(
        Tn_0_1, conj(Tn_0_1), Axes(1, 2, 3), Axes(1, 2, 3)
      ), Axes(0, 2), Axes(0, 2)
    ), Axes(0, 1, 2, 3), Axes(0, 2, 1, 3)
  );
}

template <class tensor>
typename tensor::value_type
Contract_two_sites_horizontal_op12_MF(
  const tensor &Tn_0_0,
  const tensor &Tn_1_0,
  const tensor &op12
)
{
  ////////////////////////////////////////////////////////////
  // hoge.dat
  ////////////////////////////////////////////////////////////
  // (op12*((Tn_0_0*conj(Tn_0_0))*(Tn_1_0*conj(Tn_1_0))))
  // cpu_cost= 4.46054e+06  memory= 139264
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return
  trace(
    op12, tensordot(
      tensordot(
        Tn_0_0, conj(Tn_0_0), Axes(0, 1, 2), Axes(0, 1, 2)
      ), tensordot(
        Tn_1_0, conj(Tn_1_0), Axes(0, 2, 3), Axes(0, 2, 3)
      ), Axes(0, 2), Axes(0, 2)
    ), Axes(0, 1, 2, 3), Axes(0, 2, 1, 3)
  )
  ;
}



template <class tensor>
typename tensor::value_type
Contract_four_sites_MF(
  const tensor &Tn_0_0,
  const tensor &Tn_0_1,
  const tensor &Tn_1_0,
  const tensor &Tn_1_1,
  const tensor &op_0_0,
  const tensor &op_0_1,
  const tensor &op_1_0,
  const tensor &op_1_1
)
{
  ////////////////////////////////////////////////////////////
  // hoge.dat
  ////////////////////////////////////////////////////////////
  // (op_0_0*(Tn_0_0*(conj(Tn_0_0)*((Tn_0_1*(conj(Tn_0_1)*op_0_1))*((Tn_1_0*(conj(Tn_1_0)*op_1_0))*(Tn_1_1*(conj(Tn_1_1)*op_1_1)))))))
  // cpu_cost= 9.96154e+06  memory= 295168
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return
  trace(
    op_0_0, tensordot(
      Tn_0_0, tensordot(
        conj(Tn_0_0), tensordot(
          tensordot(
            Tn_0_1, tensordot(
              conj(Tn_0_1), op_0_1, Axes(4), Axes(1)
            ), Axes(1, 2, 4), Axes(1, 2, 4)
          ), tensordot(
            tensordot(
              Tn_1_0, tensordot(
                conj(Tn_1_0), op_1_0, Axes(4), Axes(1)
              ), Axes(0, 3, 4), Axes(0, 3, 4)
            ), tensordot(
              Tn_1_1, tensordot(
                conj(Tn_1_1), op_1_1, Axes(4), Axes(1)
              ), Axes(2, 3, 4), Axes(2, 3, 4)
            ), Axes(1, 3), Axes(0, 2)
          ), Axes(1, 3), Axes(2, 3)
        ), Axes(2, 3), Axes(1, 3)
      ), Axes(0, 1, 2, 3), Axes(0, 1, 3, 4)
    ), Axes(0, 1), Axes(0, 1)
  )
  ;
}


// environment

template <template <typename> class Matrix, typename C>
class Mult_col {
 public:
  Mult_col(const Tensor<Matrix, C> &LT, const Tensor<Matrix, C> &LB)
      : LT_(LT), LB_(LB){};
  Tensor<Matrix, C> operator()(const Tensor<Matrix, C> &T_in) {
    return tensordot(LT_, tensordot(LB_, T_in, Axes(1, 3, 5), Axes(0, 1, 2)),
                     Axes(1, 3, 5), Axes(0, 1, 2));
  };

 private:
  const Tensor<Matrix, C> &LT_;
  const Tensor<Matrix, C> &LB_;
};

template <template <typename> class Matrix, typename C>
class Mult_row {
 public:
  Mult_row(const Tensor<Matrix, C> &LT, const Tensor<Matrix, C> &LB)
      : LT_(LT), LB_(LB){};
  Tensor<Matrix, C> operator()(const Tensor<Matrix, C> &T_in) {
    return tensordot(tensordot(T_in, LT_, Axes(0, 1, 2), Axes(0, 2, 4)), LB_,
                     Axes(1, 2, 3), Axes(0, 2, 4));
  };

 private:
  const Tensor<Matrix, C> &LT_;
  const Tensor<Matrix, C> &LB_;
};
template <template <typename> class Matrix, typename C>
class Mult_col_ud {
 public:
  Mult_col_ud(const Tensor<Matrix, C> &LT, const Tensor<Matrix, C> &RT,
              const Tensor<Matrix, C> &RB, const Tensor<Matrix, C> &LB)
      : LT_(LT), RT_(RT), RB_(RB), LB_(LB){};
  Tensor<Matrix, C> operator()(const Tensor<Matrix, C> &T_in) {
    return tensordot(
        RT_,
        tensordot(
            LT_,
            tensordot(LB_, tensordot(RB_, T_in, Axes(1, 3, 5), Axes(0, 1, 2)),
                      Axes(1, 3, 5), Axes(0, 1, 2)),
            Axes(1, 3, 5), Axes(0, 1, 2)),
        Axes(0, 3, 4), Axes(0, 1, 2));
  };

 private:
  const Tensor<Matrix, C> &LT_;
  const Tensor<Matrix, C> &RT_;
  const Tensor<Matrix, C> &RB_;
  const Tensor<Matrix, C> &LB_;
};

template <template <typename> class Matrix, typename C>
class Mult_row_ud {
 public:
  Mult_row_ud(const Tensor<Matrix, C> &LT, const Tensor<Matrix, C> &RT,
              const Tensor<Matrix, C> &RB, const Tensor<Matrix, C> &LB)
      : LT_(LT), RT_(RT), RB_(RB), LB_(LB){};
  Tensor<Matrix, C> operator()(const Tensor<Matrix, C> &T_in) {
    return tensordot(
        tensordot(tensordot(tensordot(T_in, RT_, Axes(0, 1, 2), Axes(1, 2, 5)),
                            LT_, Axes(1, 2, 3), Axes(0, 2, 4)),
                  LB_, Axes(1, 2, 3), Axes(0, 2, 4)),
        RB_, Axes(1, 2, 3), Axes(0, 2, 4));
  };

 private:
  const Tensor<Matrix, C> &LT_;
  const Tensor<Matrix, C> &RT_;
  const Tensor<Matrix, C> &RB_;
  const Tensor<Matrix, C> &LB_;
};

template <template <typename> class Matrix, typename C>
void Calc_projector_left_block(
    const Tensor<Matrix, C> &C1, const Tensor<Matrix, C> &C4,
    const Tensor<Matrix, C> &eT1, const Tensor<Matrix, C> &eT6,
    const Tensor<Matrix, C> &eT7, const Tensor<Matrix, C> &eT8,
    const Tensor<Matrix, C> &Tn1, const Tensor<Matrix, C> &Tn4,
    const PEPS_Parameters peps_parameters, Tensor<Matrix, C> &PU,
    Tensor<Matrix, C> &PL) {
  // Original (cheaper version of P. Corboz, T.M.Rice and M. Troyer, PRL 113,
  // 046402(2014))

  int e12 = eT1.shape()[1];
  int e56 = eT6.shape()[0];
  int e78 = eT8.shape()[0];
  int t12 = Tn1.shape()[2];
  int t41 = Tn1.shape()[3];
  int t34 = Tn4.shape()[2];

  if (t41 != 1) {
    /*
      INFO:104 (2,3) Finish 72/104 script=[3, 0, 1, -1, 2, -1, 4, -1, -1]
      ##############################
      # (Tn1*(((C1*eT1)*eT8)*Tn1c))
      # cpu_cost= 5.01e+10  memory= 3.0004e+08
      # final_bond_order  (t12, t41, e12, e78, tc12, tc41)
      ##############################
    */
    Tensor<Matrix, C> LT =
        tensordot(Tn1,
                  tensordot(tensordot(tensordot(C1, eT1, Axes(1), Axes(0)), eT8,
                                      Axes(0), Axes(1)),
                            conj(Tn1), Axes(2, 5), Axes(1, 0)),
                  Axes(0, 1, 4), Axes(3, 1, 6));

    /*
      INFO:104 (2,3) Finish 74/104 script=[3, 0, 2, -1, 1, -1, 4, -1, -1]
      ##############################
      # (Tn4*(((C4*eT7)*eT6)*Tn4c))
      # cpu_cost= 5.01e+10  memory= 3.0004e+08
      # final_bond_order  (t41, t34, e78, e56, tc41, tc34)
      ##############################
    */
    Tensor<Matrix, C> LB =
        tensordot(Tn4,
                  tensordot(tensordot(tensordot(C4, eT7, Axes(1), Axes(0)), eT6,
                                      Axes(0), Axes(1)),
                            conj(Tn4), Axes(2, 5), Axes(0, 3)),
                  Axes(0, 3, 4), Axes(1, 3, 6));

    // Tensor<Matrix, C> R1 = LT;
    // Tensor<Matrix, C> R2 = LB;
    if ((t12 != 1 && t34 != 1) && peps_parameters.Use_RSVD) {
      Mult_col<Matrix, C> m_col(LT, LB);
      Mult_row<Matrix, C> m_row(LT, LB);

      Tensor<Matrix, C> U;
      Tensor<Matrix, C> VT;
      std::vector<double> s;

      Shape shape_row(t12, e12, t12);
      Shape shape_col(t34, e56, t34);

      /*int info ()= */
      rsvd(m_row, m_col, shape_row, shape_col, U, s, VT, e78,
           static_cast<size_t>(peps_parameters.RSVD_Oversampling_factor * e78));
      double denom = s[0];

      for (int i = 0; i < e78; ++i) {
        if (s[i] / denom > peps_parameters.Inverse_projector_cut) {
          s[i] = 1.0 / sqrt(s[i]);
        } else {
          s[i] = 0.0;
        }
      }

      U.multiply_vector(s, 3);
      VT.multiply_vector(s, 0);

      PU = tensordot(LB, conj(VT), Axes(1, 3, 5), Axes(1, 2, 3))
               .transpose(Axes(1, 0, 2, 3));
      PL = tensordot(LT, conj(U), Axes(0, 2, 4), Axes(0, 1, 2))
               .transpose(Axes(1, 0, 2, 3));
    } else {
      // full svd //
      Tensor<Matrix, C> U;
      Tensor<Matrix, C> VT;
      std::vector<double> s;

      /* int info = */
      svd(tensordot(LT, LB, Axes(1, 3, 5), Axes(0, 2, 4)), Axes(0, 1, 2),
          Axes(3, 4, 5), U, s, VT);
      double denom = s[0];
      std::vector<double> s_c;
      s_c.resize(e78);

      for (int i = 0; i < e78; ++i) {
        if (s[i] / denom > peps_parameters.Inverse_projector_cut) {
          s_c[i] = 1.0 / sqrt(s[i]);
        } else {
          s_c[i] = 0.0;
        }
      }
      // O(D^{10})
      Tensor<Matrix, C> U_c = slice(U, 3, 0, e78);
      Tensor<Matrix, C> VT_c = slice(VT, 0, 0, e78);

      U_c.multiply_vector(s_c, 3);
      VT_c.multiply_vector(s_c, 0);

      PU = tensordot(LB, conj(VT_c), Axes(1, 3, 5), Axes(1, 2, 3))
               .transpose(Axes(1, 0, 2, 3));
      PL = tensordot(LT, conj(U_c), Axes(0, 2, 4), Axes(0, 1, 2))
               .transpose(Axes(1, 0, 2, 3));
    }
  } else {
    Tensor<Matrix, C> identity_matrix(Shape(e78, e78));
    Index index;
    for (int i = 0; i < identity_matrix.local_size(); i++) {
      index = identity_matrix.global_index(i);
      if (index[0] == index[1]) {
        identity_matrix.set_value(index, 1.0);
      } else {
        identity_matrix.set_value(index, 0.0);
      }
    }
    PU = reshape(identity_matrix, Shape(e78, t41, t41, e78));
    PL = reshape(identity_matrix, Shape(e78, t41, t41, e78));
  }
}

template <template <typename> class Matrix, typename C>
void Calc_projector_updown_blocks(
    const Tensor<Matrix, C> &C1, const Tensor<Matrix, C> &C2,
    const Tensor<Matrix, C> &C3, const Tensor<Matrix, C> &C4,
    const Tensor<Matrix, C> &eT1, const Tensor<Matrix, C> &eT2,
    const Tensor<Matrix, C> &eT3, const Tensor<Matrix, C> &eT4,
    const Tensor<Matrix, C> &eT5, const Tensor<Matrix, C> &eT6,
    const Tensor<Matrix, C> &eT7, const Tensor<Matrix, C> &eT8,
    const Tensor<Matrix, C> &Tn1, const Tensor<Matrix, C> &Tn2,
    const Tensor<Matrix, C> &Tn3, const Tensor<Matrix, C> &Tn4,
    const PEPS_Parameters peps_parameters, Tensor<Matrix, C> &PU,
    Tensor<Matrix, C> &PL) {
  // based on P. Corboz, T.M.Rice and M. Troyer, PRL 113, 046402(2014)

  // comment out for unused variables
  // int e12 = eT1.shape()[1];
  int e34 = eT3.shape()[1];
  // int e56 = eT5.shape()[1];
  int e78 = eT7.shape()[1];
  // int t12 = Tn1.shape()[2];
  int t23 = Tn2.shape()[3];
  // int t34 = Tn3.shape()[0];
  int t41 = Tn4.shape()[1];

  if (t41 != 1) {
    /*
      INFO:104 (2,3) Finish 72/104 script=[3, 0, 1, -1, 2, -1, 4, -1, -1]
      ##############################
      # (Tn1*(((C1*eT1)*eT8)*Tn1c))
      # cpu_cost= 5.01e+10  memory= 3.0004e+08
      # final_bond_order  (t12, t41, e12, e78, tc12, tc41)
      ##############################
    */
    Tensor<Matrix, C> LT =
        tensordot(Tn1,
                  tensordot(tensordot(tensordot(C1, eT1, Axes(1), Axes(0)), eT8,
                                      Axes(0), Axes(1)),
                            conj(Tn1), Axes(2, 5), Axes(1, 0)),
                  Axes(0, 1, 4), Axes(3, 1, 6));

    /*
      INFO:104 (2,3) Finish 74/104 script=[3, 0, 2, -1, 1, -1, 4, -1, -1]
      ##############################
      # (Tn2*(((C2*eT3)*eT2)*Tn2c))
      # cpu_cost= 5.01e+10  memory= 3.0004e+08
      # final_bond_order  (t12, t23, e34, e12, tc12, tc23)
      ##############################
    */

    Tensor<Matrix, C> RT =
        tensordot(Tn2,
                  tensordot(tensordot(tensordot(C2, eT3, Axes(1), Axes(0)), eT2,
                                      Axes(0), Axes(1)),
                            conj(Tn2), Axes(2, 5), Axes(2, 1)),
                  Axes(1, 2, 4), Axes(3, 1, 6));

    /* INFO:104 (2,3) Finish 74/104 script=[3, 0, 2, -1, 1, -1, 4, -1, -1]
       ##############################
       # (Tn3*(((C3*eT5)*eT4)*Tn3c))
       # cpu_cost= 5.01e+10  memory= 3.0004e+08
       # final_bond_order  (t34, t23, e56, e34, tc34, tc23)
       ##############################
    */
    Tensor<Matrix, C> RB =
        tensordot(Tn3,
                  tensordot(tensordot(tensordot(C3, eT5, Axes(1), Axes(0)), eT4,
                                      Axes(0), Axes(1)),
                            conj(Tn3), Axes(2, 5), Axes(3, 2)),
                  Axes(2, 3, 4), Axes(3, 1, 6));

    /*
      INFO:104 (2,3) Finish 74/104 script=[3, 0, 2, -1, 1, -1, 4, -1, -1]
      ##############################
      # (Tn4*(((C4*eT7)*eT6)*Tn4c))
      # cpu_cost= 5.01e+10  memory= 3.0004e+08
      # final_bond_order  (t41, t34, e78, e56, tc41, tc34)
      ##############################
    */
    Tensor<Matrix, C> LB =
        tensordot(Tn4,
                  tensordot(tensordot(tensordot(C4, eT7, Axes(1), Axes(0)), eT6,
                                      Axes(0), Axes(1)),
                            conj(Tn4), Axes(2, 5), Axes(0, 3)),
                  Axes(0, 3, 4), Axes(1, 3, 6));
    if (t23 != 1 && peps_parameters.Use_RSVD) {
      Mult_col_ud<Matrix, C> m_col(LT, RT, RB, LB);
      Mult_row_ud<Matrix, C> m_row(LT, RT, RB, LB);

      Tensor<Matrix, C> U;
      Tensor<Matrix, C> VT;
      std::vector<double> s;

      Shape shape_row(t23, e34, t23);
      Shape shape_col(t23, e34, t23);

      /* int info = */
      rsvd(m_row, m_col, shape_row, shape_col, U, s, VT, e78,
           static_cast<size_t>(peps_parameters.RSVD_Oversampling_factor * e78));
      double denom = s[0];

      for (int i = 0; i < e78; ++i) {
        if (s[i] / denom > peps_parameters.Inverse_projector_cut) {
          s[i] = 1.0 / sqrt(s[i]);
        } else {
          s[i] = 0.0;
        }
      }

      U.multiply_vector(s, 3);
      VT.multiply_vector(s, 0);

      PU = tensordot(LB, tensordot(RB, conj(VT), Axes(1, 3, 5), Axes(1, 2, 3)),
                     Axes(1, 3, 5), Axes(0, 1, 2))
               .transpose(Axes(1, 0, 2, 3));
      PL = tensordot(LT, tensordot(RT, conj(U), Axes(1, 2, 5), Axes(0, 1, 2)),
                     Axes(0, 2, 4), Axes(0, 1, 2))
               .transpose(Axes(1, 0, 2, 3));
    } else {
      // full svd
      Tensor<Matrix, C> R1 = tensordot(RT, LT, Axes(0, 3, 4), Axes(0, 2, 4));
      Tensor<Matrix, C> R2 = tensordot(RB, LB, Axes(0, 2, 4), Axes(1, 3, 5));

      Tensor<Matrix, C> U;
      Tensor<Matrix, C> VT;
      std::vector<double> s;

      /* int info = */
      svd(tensordot(R1, R2, Axes(3, 4, 5), Axes(3, 4, 5)), Axes(0, 1, 2),
          Axes(3, 4, 5), U, s, VT);
      double denom = s[0];
      std::vector<double> s_c;
      s_c.resize(e78);

      for (int i = 0; i < e78; ++i) {
        if (s[i] / denom > peps_parameters.Inverse_projector_cut) {
          s_c[i] = 1.0 / sqrt(s[i]);
        } else {
          s_c[i] = 0.0;
        }
      }
      Tensor<Matrix, C> U_c = slice(U, 3, 0, e78);
      Tensor<Matrix, C> VT_c = slice(VT, 0, 0, e78);

      U_c.multiply_vector(s_c, 3);
      VT_c.multiply_vector(s_c, 0);

      PU = tensordot(R2, conj(VT_c), Axes(0, 1, 2), Axes(1, 2, 3))
               .transpose(Axes(1, 0, 2, 3));
      PL = tensordot(R1, conj(U_c), Axes(0, 1, 2), Axes(0, 1, 2))
               .transpose(Axes(1, 0, 2, 3));
    }
  } else {
    Tensor<Matrix, C> identity_matrix(Shape(e78, e78));
    Index index;
    for (int i = 0; i < identity_matrix.local_size(); i++) {
      index = identity_matrix.global_index(i);
      if (index[0] == index[1]) {
        identity_matrix.set_value(index, 1.0);
      } else {
        identity_matrix.set_value(index, 0.0);
      }
    }
    PU = reshape(identity_matrix, Shape(e78, t41, t41, e78));
    PL = reshape(identity_matrix, Shape(e78, t41, t41, e78));
  }
}

template <template <typename> class Matrix, typename C>
void Calc_Next_CTM(const Tensor<Matrix, C> &C1, const Tensor<Matrix, C> &C4,
                   const Tensor<Matrix, C> &eT1, const Tensor<Matrix, C> &eT6,
                   const Tensor<Matrix, C> &PU, const Tensor<Matrix, C> &PL,
                   Tensor<Matrix, C> &C1_out, Tensor<Matrix, C> &C4_out) {
  C1_out = tensordot(PU, tensordot(C1, eT1, Axes(1), Axes(0)), Axes(0, 1, 2),
                     Axes(0, 2, 3));
  C4_out = tensordot(tensordot(eT6, C4, Axes(1), Axes(0)), PL, Axes(3, 1, 2),
                     Axes(0, 1, 2));

  // normalization
  /*
  std::vector<double> w;
  svd(C1_out,w);
  double norm = std::accumulate(w.begin(),w.end(),0.0);
  C1_out /= norm;

  svd(C4_out,w);
  norm = std::accumulate(w.begin(),w.end(),0.0);
  C4_out /= norm;
  */

  double max_all = max_abs(C1_out);
  C1_out /= max_all;
  max_all = max_abs(C4_out);
  C4_out /= max_all;
}

template <template <typename> class Matrix, typename C>
void Calc_Next_eT(const Tensor<Matrix, C> &eT8, const Tensor<Matrix, C> &Tn1,
                  const Tensor<Matrix, C> &PU, const Tensor<Matrix, C> &PL,
                  Tensor<Matrix, C> &eT_out) {
  /*
    INFO:299 (1,1) Finish 258/299 script=[1, 2, 0, 4, -1, -1, -1, 3, -1]
    ##############################
    # ((Tn1*(Tn1c*(eT8*PL)))*PU)
    # cpu_cost= 6e+10  memory= 3.0104e+08
    # final_bond_order  (k2, l2, n1, n0)
    ##############################
  */

  eT_out = tensordot(tensordot(Tn1,
                               tensordot(conj(Tn1),
                                         tensordot(eT8, PL, Axes(1), Axes(0)),
                                         Axes(0, 1), Axes(2, 4)),
                               Axes(0, 1, 4), Axes(4, 5, 2)),
                     PU, Axes(1, 3, 4), Axes(1, 2, 0))
               .transpose(Axes(3, 2, 0, 1));

  // normalization
  /*
  std::vector<double> w;
  Tensor<Matrix, C> eT_temp;
  eT_temp = contract(eT_out,Axes(2),Axes(3));
  svd(eT_temp,w);
  double norm = std::accumulate(w.begin(),w.end(),0.0);
  eT_out /= norm;
  */

  double max_all = max_abs(eT_out);
  eT_out /= max_all;
}

template <template <typename> class Matrix, typename C>
void Simple_update_bond(const Tensor<Matrix, C> &Tn1,
                        const Tensor<Matrix, C> &Tn2,
                        const std::vector<std::vector<double>> &lambda1,
                        const std::vector<std::vector<double>> &lambda2,
                        const Tensor<Matrix, C> &op12, const int connect1,
                        const PEPS_Parameters peps_parameters,
                        Tensor<Matrix, C> &Tn1_new, Tensor<Matrix, C> &Tn2_new,
                        std::vector<double> &lambda_c) {
  using ptensor = Tensor<Matrix, C>;
  int connect2 = (connect1 + 2) % 4;

  std::vector<std::vector<double>> lambda1_inv(4);
  std::vector<std::vector<double>> lambda2_inv(4);

  for (int i = 0; i < 4; i++) {
    lambda1_inv[i] = std::vector<double>(lambda1[i].size());
    for (int j = 0; j < lambda1_inv[i].size(); ++j) {
      if (lambda1[i][j] > peps_parameters.Inverse_lambda_cut) {
        lambda1_inv[i][j] = 1.0 / lambda1[i][j];
      } else {
        lambda1_inv[i][j] = 0.0;
      }
    }
  };
  for (int i = 0; i < 4; i++) {
    lambda2_inv[i] = std::vector<double>(lambda2[i].size());
    for (int j = 0; j < lambda2_inv[i].size(); ++j) {
      if (lambda2[i][j] > peps_parameters.Inverse_lambda_cut) {
        lambda2_inv[i][j] = 1.0 / lambda2[i][j];
      } else {
        lambda2_inv[i][j] = 0.0;
      }
    }
  };

  int dc = Tn1.shape()[connect1];
  ptensor Tn1_lambda = Tn1;
  ptensor Tn2_lambda = Tn2;

  if (connect1 == 0) {
    Tn1_lambda.multiply_vector(lambda1[1], 1, lambda1[2], 2, lambda1[3], 3);
    Tn1_lambda.transpose(Axes(1, 2, 3, 0, 4));
  } else if (connect1 == 1) {
    Tn1_lambda.multiply_vector(lambda1[0], 0, lambda1[2], 2, lambda1[3], 3);
    Tn1_lambda.transpose(Axes(0, 2, 3, 1, 4));
  } else if (connect1 == 2) {
    Tn1_lambda.multiply_vector(lambda1[0], 0, lambda1[1], 1, lambda1[3], 3);
    Tn1_lambda.transpose(Axes(0, 1, 3, 2, 4));
  } else {
    Tn1_lambda.multiply_vector(lambda1[0], 0, lambda1[1], 1, lambda1[2], 2);
  }

  if (connect2 == 0) {
    Tn2_lambda.multiply_vector(lambda2[1], 1, lambda2[2], 2, lambda2[3], 3);
    Tn2_lambda.transpose(Axes(1, 2, 3, 0, 4));
  } else if (connect2 == 1) {
    Tn2_lambda.multiply_vector(lambda2[0], 0, lambda2[2], 2, lambda2[3], 3);
    Tn2_lambda.transpose(Axes(0, 2, 3, 1, 4));
  } else if (connect2 == 2) {
    Tn2_lambda.multiply_vector(lambda2[0], 0, lambda2[1], 1, lambda2[3], 3);
    Tn2_lambda.transpose(Axes(0, 1, 3, 2, 4));
  } else {
    Tn2_lambda.multiply_vector(lambda2[0], 0, lambda2[1], 1, lambda2[2], 2);
  };

  // QR
  ptensor Q1, R1, Q2, R2;
  int info = qr(Tn1_lambda, Axes(0, 1, 2), Axes(3, 4), Q1, R1);

  info = qr(Tn2_lambda, Axes(0, 1, 2), Axes(3, 4), Q2, R2);

  // connect R1, R2, op
  /*
    INFO:8 (1,2) Finish 7/8 script=[0, 1, -1, 2, -1]
    ##############################
    # ((R1*R2)*op12)
    # cpu_cost= 22400  memory= 3216
    # final_bond_order  (c1, c2, m1o, m2o)
    ##############################
  */
  ptensor Theta = tensordot(tensordot(R1, R2, Axes(1), Axes(1)), op12,
                                      Axes(1, 3), Axes(0, 1));

  // svd
  ptensor U, VT;
  std::vector<double> s;
  info = svd(Theta, Axes(0, 2), Axes(1, 3), U, s, VT);

  lambda_c = std::vector<double>(s.begin(), s.begin() + dc);
  ptensor Uc = slice(U, 2, 0, dc);
  ptensor VTc = slice(VT, 0, 0, dc);

  //  norm =
  //  std::inner_product(lambda_c.begin(),lambda_c.end(),lambda_c.begin(),0.0);

  double norm = 0.0;
  for (int i = 0; i < dc; ++i) {
    norm += lambda_c[i] * lambda_c[i];
  };
  norm = sqrt(norm);
  for (int i = 0; i < dc; ++i) {
    lambda_c[i] = sqrt(lambda_c[i] / norm);
  };

  /*for (int i=0; i < VTc.local_size();++i){
    Index index = VTc.global_index(i);
    std::cout<<"VTC[i,j]="<<index<<", "<<VTc[i]<<std::endl;
    }*/

  if(peps_parameters.Simple_Gauge_Fix){
    ////////////////////////////////////////////////////////////
    // (Uc*(conj(Uc)*(Q1*conj(Q1))))
    // cpu_cost= 1280  memory= 592
    // final_bond_order (UcD, UcD')
    ////////////////////////////////////////////////////////////
    ptensor M1 = tensordot(
      Uc, tensordot(
        conj(Uc), tensordot(
          Q1, conj(Q1), Axes(0, 1, 2), Axes(0, 1, 2)
        ), Axes(0), Axes(1)
      ), Axes(0, 1), Axes(2, 0)
    );

    ////////////////////////////////////////////////////////////
    // (VTc*(conj(VTc)*(Q2*conj(Q2))))
    // cpu_cost= 1280  memory= 592
    // final_bond_order (VTcD, VTcD')
    ////////////////////////////////////////////////////////////
    ptensor M2 = tensordot(
      VTc, tensordot(
        conj(VTc), tensordot(
          Q2, conj(Q2), Axes(0, 1, 2), Axes(0, 1, 2)
        ), Axes(1), Axes(1)
      ), Axes(1, 2), Axes(2, 1)
    );

    ptensor U1, U2;
    std::vector<double> D1, D2;
    eigh(M1, Axes(0), Axes(1), D1, U1);
    eigh(M2, Axes(0), Axes(1), D2, U2);

    std::vector<double> D1inv(dc), D2inv(dc), lcinv(dc);
    for(int d=0; d<dc; ++d){
      D1[d] = sqrt(D1[d]);
      D1inv[d] = 1.0/D1[d];

      D2[d] = sqrt(D2[d]);
      D2inv[d] = 1.0/D2[d];

      lcinv[d] = 1.0/lambda_c[d];
    }

    U1.multiply_vector(lambda_c, 0, D1, 1);
    U2.multiply_vector(lambda_c, 0, D2, 1);
    ptensor L = tensordot(U1, U2, Axes(0), Axes(0));

    // revert U1, U2
    U1.multiply_vector(lcinv, 0, D1inv, 1);
    U2.multiply_vector(lcinv, 0, D2inv, 1);

    ptensor W1, W2;
    info = svd(L, Axes(0), Axes(1), W1, lambda_c, W2);

    norm = 0.0;
    for(int d=0; d<dc; ++d){
      norm += lambda_c[d] * lambda_c[d];
    }
    norm = sqrt(norm);

    for(int d=0; d<dc; ++d){
      lambda_c[d] = sqrt(lambda_c[d] / norm);
      if (lambda_c[d] > peps_parameters.Inverse_lambda_cut){
        lcinv[d] = 1.0/lambda_c[d];
      }else{
        lcinv[d] = 0.0;
      }
    }

    W1.multiply_vector(D1inv, 0);
    W2.multiply_vector(D2inv, 1);
    
    ptensor X1inv = tensordot(conj(U1), W1, Axes(1), Axes(0));
    ptensor X2inv = tensordot(conj(U2), W2, Axes(1), Axes(1));
    Uc = tensordot(Uc, X1inv, Axes(2), Axes(0));
    VTc = tensordot(X2inv, VTc, Axes(0), Axes(0));
  } // end of local gauge fixing

  Uc.multiply_vector(lambda_c, 2);
  VTc.multiply_vector(lambda_c, 0);

  // Remove lambda effects from Qs
  // and create new tensors
  if (connect1 == 0) {
    Q1.multiply_vector(lambda1_inv[1], 0, lambda1_inv[2], 1, lambda1_inv[3], 2);
    Tn1_new =
        tensordot(Q1, Uc, Axes(3), Axes(0)).transpose(Axes(4, 0, 1, 2, 3));
  } else if (connect1 == 1) {
    Q1.multiply_vector(lambda1_inv[0], 0, lambda1_inv[2], 1, lambda1_inv[3], 2);
    Tn1_new =
        tensordot(Q1, Uc, Axes(3), Axes(0)).transpose(Axes(0, 4, 1, 2, 3));
  } else if (connect1 == 2) {
    Q1.multiply_vector(lambda1_inv[0], 0, lambda1_inv[1], 1, lambda1_inv[3], 2);
    Tn1_new =
        tensordot(Q1, Uc, Axes(3), Axes(0)).transpose(Axes(0, 1, 4, 2, 3));
  } else {
    Q1.multiply_vector(lambda1_inv[0], 0, lambda1_inv[1], 1, lambda1_inv[2], 2);
    Tn1_new =
        tensordot(Q1, Uc, Axes(3), Axes(0)).transpose(Axes(0, 1, 2, 4, 3));
  };

  if (connect2 == 0) {
    Q2.multiply_vector(lambda2_inv[1], 0, lambda2_inv[2], 1, lambda2_inv[3], 2);
    Tn2_new =
        tensordot(Q2, VTc, Axes(3), Axes(1)).transpose(Axes(3, 0, 1, 2, 4));
  } else if (connect2 == 1) {
    Q2.multiply_vector(lambda2_inv[0], 0, lambda2_inv[2], 1, lambda2_inv[3], 2);
    Tn2_new =
        tensordot(Q2, VTc, Axes(3), Axes(1)).transpose(Axes(0, 3, 1, 2, 4));
  } else if (connect2 == 2) {
    Q2.multiply_vector(lambda2_inv[0], 0, lambda2_inv[1], 1, lambda2_inv[3], 2);
    Tn2_new =
        tensordot(Q2, VTc, Axes(3), Axes(1)).transpose(Axes(0, 1, 3, 2, 4));
  } else {
    Q2.multiply_vector(lambda2_inv[0], 0, lambda2_inv[1], 1, lambda2_inv[2], 2);
    Tn2_new =
        tensordot(Q2, VTc, Axes(3), Axes(1)).transpose(Axes(0, 1, 2, 3, 4));
  };
}

// for full update
template <template <typename> class Matrix, typename C>
Tensor<Matrix, C> Create_Environment_two_sites(
    const Tensor<Matrix, C> &C1, const Tensor<Matrix, C> &C2,
    const Tensor<Matrix, C> &C3, const Tensor<Matrix, C> &C4,
    const Tensor<Matrix, C> &eT1, const Tensor<Matrix, C> &eT2,
    const Tensor<Matrix, C> &eT3, const Tensor<Matrix, C> &eT4,
    const Tensor<Matrix, C> &eT5, const Tensor<Matrix, C> &eT6,
    const Tensor<Matrix, C> &Q1, const Tensor<Matrix, C> &Q2) {
  /* C1 - eT1 - eT2 - C2
    ## |    |     |     |
    ##eT6 - Q1-  -Q2  - eT3
    ## |    |     |     |
    ## C4 - eT5 - eT4 - C3
    ##
    ## Q1,Q2 are double layered
  */
  return tensordot(
             tensordot(tensordot(C2, eT2, Axes(0), Axes(1)),
                       tensordot(tensordot(tensordot(tensordot(C3, eT4, Axes(1),
                                                               Axes(0)),
                                                     eT3, Axes(0), Axes(1)),
                                           conj(Q2), Axes(2, 5), Axes(3, 2)),
                                 Q2, Axes(1, 3), Axes(3, 2)),
                       Axes(0, 2, 3), Axes(1, 5, 3)),
             tensordot(tensordot(C1, eT1, Axes(1), Axes(0)),
                       tensordot(tensordot(tensordot(tensordot(C4, eT6, Axes(1),
                                                               Axes(0)),
                                                     eT5, Axes(0), Axes(1)),
                                           conj(Q1), Axes(2, 5), Axes(0, 2)),
                                 Q1, Axes(1, 3), Axes(0, 2)),
                       Axes(0, 2, 3), Axes(0, 4, 2)),
             Axes(0, 1), Axes(0, 1))
      .transpose(Axes(3, 1, 2, 0));
}

template <template <typename> class Matrix, typename C>
void Full_update_bond_horizontal(
    const Tensor<Matrix, C> &C1, const Tensor<Matrix, C> &C2,
    const Tensor<Matrix, C> &C3, const Tensor<Matrix, C> &C4,
    const Tensor<Matrix, C> &eT1, const Tensor<Matrix, C> &eT2,
    const Tensor<Matrix, C> &eT3, const Tensor<Matrix, C> &eT4,
    const Tensor<Matrix, C> &eT5, const Tensor<Matrix, C> &eT6,
    const Tensor<Matrix, C> &Tn1, const Tensor<Matrix, C> &Tn2,
    const Tensor<Matrix, C> &op12, const PEPS_Parameters peps_parameters,
    Tensor<Matrix, C> &Tn1_new, Tensor<Matrix, C> &Tn2_new) {
  Shape Tn1_shape = Tn1.shape();
  Shape Tn2_shape = Tn2.shape();

  int D_connect = Tn1_shape[2];

  // Connecting [2] bond of Tn1 and [0] bond of Tn2
  // QR decomposition
  Tensor<Matrix, C> Q1, R1, Q2, R2;

  int info = qr(Tn1, Axes(0, 1, 3), Axes(2, 4), Q1, R1);
  info = qr(Tn2, Axes(1, 2, 3), Axes(0, 4), Q2, R2);

  int envR1 = R1.shape()[0];
  int envR2 = R2.shape()[0];

  /*
    ## apply time evolution
    INFO:8 (1,2) Finish 7/8 script=[0, 1, -1, 2, -1]
    ##############################
    # ((R1*R2)*op12)
    # cpu_cost= 22400  memory= 3216
    # final_bond_order  (c1, c2, m1o, m2o)
    ##############################
  */

  Tensor<Matrix, C> Theta = tensordot(tensordot(R1, R2, Axes(1), Axes(1)), op12,
                                      Axes(1, 3), Axes(0, 1));
  // Environment
  // bond order (t1, t2, tc1, tc2)

  Tensor<Matrix, C> Environment =
      Create_Environment_two_sites(C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6,
                                   Q1, transpose(Q2, Axes(3, 0, 1, 2)));

  // Hermite
  Environment =
      0.5 * (Environment + conj(transpose(Environment, Axes(2, 3, 0, 1))));

  // diagonalization
  Tensor<Matrix, C> Z;
  std::vector<double> w;

  info = eigh(Environment, Axes(0, 1), Axes(2, 3), w, Z);

  // positive
  double w_max = fabs(w[envR1 * envR2 - 1]);
  for (int i = 0; i < envR1 * envR2; ++i) {
    if (w[i] / w_max > peps_parameters.Inverse_Env_cut) {
      w[i] = sqrt(w[i]);
    } else {
      w[i] = 0.0;
    }
  };

  Z.multiply_vector(w, 2);

  Tensor<Matrix, C> LR1, LR2, LR1_inv, LR2_inv;
  if (peps_parameters.Full_Gauge_Fix) {
    // gauge fix
    Tensor<Matrix, C> Q_temp;
    info = qr(Z, Axes(2, 1), Axes(0), Q_temp, LR1);
    info = qr(Z, Axes(2, 0), Axes(1), Q_temp, LR2);
    // for theta
    Theta = tensordot(LR1, tensordot(LR2, Theta, Axes(1), Axes(1)), Axes(1),
                      Axes(1));

    // for environment
    Tensor<Matrix, C> u, vt;
    std::vector<double> s;

    info = svd(LR1, u, s, vt);
    double s_max = s[0];
    for (int i = 0; i < envR1; ++i) {
      if (s[i] / s_max > peps_parameters.Inverse_Env_cut) {
        s[i] = 1.0 / s[i];
      } else {
        s[i] = 0.0;
      }
    }
    u.multiply_vector(s, 1);

    LR1_inv = tensordot(conj(u), conj(vt), Axes(1), Axes(0));

    info = svd(LR2, u, s, vt);
    s_max = s[0];
    for (int i = 0; i < envR2; ++i) {
      if (s[i] / s_max > peps_parameters.Inverse_Env_cut) {
        s[i] = 1.0 / s[i];
      } else {
        s[i] = 0.0;
      }
    }
    u.multiply_vector(s, 1);

    LR2_inv = tensordot(conj(u), conj(vt), Axes(1), Axes(0));

    Z = tensordot(tensordot(Z, LR1_inv, Axes(0), Axes(1)), LR2_inv, Axes(0),
                  Axes(1));

    Environment = tensordot(Z, conj(Z), Axes(0), Axes(0));

  } else {
    Environment = tensordot(Z, conj(Z), Axes(2), Axes(2));
  }

  /*
  // test//
  for (int i=0; i < Environment.local_size();++i){
    Index index = Environment.global_index(i);
    if (index[2]==0 && index[3]==0){
      std::cout<<"Environment[i,j,0,0]= "<<index<<",
  "<<Environment[i]<<std::endl;
    }
    }*/

  bool convergence = false;

  // create initial guess
  // SVD of Theta

  Tensor<Matrix, C> U, VT;
  std::vector<double> s;

  svd(Theta, Axes(0, 2), Axes(1, 3), U, s, VT);

  // truncation
  std::vector<double> lambda_c = s;
  lambda_c.resize(D_connect);
  U = slice(U, 2, 0, D_connect);
  VT = slice(VT, 0, 0, D_connect);

  double norm = 0.0;
  for (int i = 0; i < D_connect; ++i) {
    norm += lambda_c[i] * lambda_c[i];
  };

  norm = sqrt(norm);

  for (int i = 0; i < D_connect; ++i) {
    lambda_c[i] = sqrt(lambda_c[i] / norm);
  };

  U.multiply_vector(lambda_c, 2);
  VT.multiply_vector(lambda_c, 0);
  R1 = transpose(U, Axes(0, 2, 1));   // envR1 , D_connect, m1
  R2 = transpose(VT, Axes(1, 0, 2));  // envR2 , D_connect, m2

  int count = 0;
  C delta = 0;
  C C_phi = trace(tensordot(Environment, Theta, Axes(0, 1), Axes(0, 1)),
                conj(Theta), Axes(0, 1, 2, 3), Axes(0, 1, 2, 3));
  C Old_delta =
      (-2.0 * trace(tensordot(
                        tensordot(Environment, Theta, Axes(0, 1), Axes(0, 1)),
                        conj(R2), Axes(1, 3), Axes(0, 2)),
                    conj(R1), Axes(0, 1, 2), Axes(0, 2, 1)) +
       trace(
           R1,
           tensordot(R2,
                     tensordot(Environment,
                               tensordot(conj(R1), conj(R2), Axes(1), Axes(1)),
                               Axes(2, 3), Axes(0, 2)),
                     Axes(0, 2), Axes(1, 3)),
           Axes(0, 1, 2), Axes(1, 0, 2)));

  Tensor<Matrix, C> W_vec, N_mat, N_mat_inv;
  double denom;
  while (!convergence && (count < peps_parameters.Full_max_iteration)) {
    /*
      ## for R1
      ## create W
      ## ((envR1, m1, D_connect,m1)*)
    */

    W_vec =
        tensordot(tensordot(Environment, Theta, Axes(0, 1), Axes(0, 1)),
                  conj(R2), Axes(1, 3),
                  Axes(0, 2));  // transpose(0,2,1)).reshape(envR1*D_connect,m1)
    /*
      ## create N
      ## (envR1, envR1*, D_connect,D_connect*)
    */
    N_mat = tensordot(tensordot(Environment, R2, Axes(1), Axes(0)), conj(R2),
                      Axes(2, 4), Axes(0, 2));
    // transpose(1,3,0,2).reshape(envR1*D_connect,envR1*D_connect)

    // Moore-Penrose Psude Inverse (for Hermitian matrix)
    info = svd(N_mat, Axes(1, 3), Axes(0, 2), U, s, VT);
    denom = s[0];
    for (int i = 0; i < envR1 * D_connect; ++i) {
      if (s[i] / denom > peps_parameters.Full_Inverse_precision) {
        s[i] = 1.0 / s[i];
      } else {
        s[i] = 0.0;
      }
    };

    U.multiply_vector(s, 2);
    N_mat_inv = tensordot(conj(U), conj(VT), Axes(2), Axes(0));
    R1 = tensordot(N_mat_inv, W_vec, Axes(0, 1), Axes(0, 2));

    /*
      ## for R2
      ## create W
      ## ((envR2,m2, D_connect)*)
    */
    W_vec =
        tensordot(tensordot(Environment, Theta, Axes(0, 1), Axes(0, 1)),
                  conj(R1), Axes(0, 2),
                  Axes(0, 2));  //).transpose(0,2,1).reshape(envR2*D_connect,m2)

    /*
      ## create N
      ## (envR2, envR2*,D_connect, D_connect*)
    */
    N_mat = tensordot(
        tensordot(Environment, R1, Axes(0), Axes(0)), conj(R1), Axes(1, 4),
        Axes(
            0,
            2));  //.transpose(1,3,0,2).reshape(envR2*D_connect,envR2*D_connect)

    // Moore-Penrose Psude Inverse (for Hermitian matrix)
    info = svd(N_mat, Axes(1, 3), Axes(0, 2), U, s, VT);

    denom = s[0];
    for (int i = 0; i < envR2 * D_connect; ++i) {
      if (s[i] / denom > peps_parameters.Full_Inverse_precision) {
        s[i] = 1.0 / s[i];
      } else {
        s[i] = 0.0;
      }
    };

    U.multiply_vector(s, 2);
    N_mat_inv = tensordot(conj(U), conj(VT), Axes(2), Axes(0));

    R2 = tensordot(N_mat_inv, W_vec, Axes(0, 1), Axes(0, 2));

    delta =
        -2.0 * trace(tensordot(
                         tensordot(Environment, Theta, Axes(0, 1), Axes(0, 1)),
                         conj(R2), Axes(1, 3), Axes(0, 2)),
                     conj(R1), Axes(0, 1, 2), Axes(0, 2, 1)) +
        trace(
            R1,
            tensordot(R2,
                      tensordot(Environment,
                                tensordot(conj(R1), conj(R2), Axes(1), Axes(1)),
                                Axes(2, 3), Axes(0, 2)),
                      Axes(0, 2), Axes(1, 3)),
            Axes(0, 1, 2), Axes(1, 0, 2));

    // std::cout<<"delta "<<delta<<std::endl;

    /*
    if (fabs(Old_delta - delta)/C_phi < Full_Convergence_Epsilon){
      convergence = true;
    };
    */
    if (std::abs(Old_delta - delta) / std::abs(C_phi) <
        peps_parameters.Full_Convergence_Epsilon) {
      convergence = true;
    };
    Old_delta = delta;
    count += 1;
  }
  // Post processing
  if (!convergence &&
      peps_parameters.print_level >= PrintLevel::warn) {
    std::cout << "warning: Full update iteration was not conveged! count= "
              << count << std::endl;
  }
  if (peps_parameters.print_level >= PrintLevel::debug) {
    std::cout << "Full Update: count, delta,original_norm = " << count << " "
              << delta + C_phi << " " << C_phi << std::endl;
  }
  if (peps_parameters.Full_Gauge_Fix) {
    // remove gauge

    R1 = tensordot(LR1_inv, R1, Axes(0), Axes(0));
    R2 = tensordot(LR2_inv, R2, Axes(0), Axes(0));
  }
  // balancing and normalization

  Tensor<Matrix, C> q1, r1, q2, r2;

  info = qr(R1, Axes(0, 2), Axes(1), q1, r1);
  info = qr(R2, Axes(0, 2), Axes(1), q2, r2);

  info = svd(tensordot(r1, r2, Axes(1), Axes(1)), U, s, VT);

  norm = 0.0;
  for (int i = 0; i < D_connect; ++i) {
    norm += s[i] * s[i];
  }
  norm = sqrt(norm);
  for (int i = 0; i < D_connect; ++i) {
    s[i] = sqrt(s[i] / norm);
  }

  U.multiply_vector(s, 1);
  VT.multiply_vector(s, 0);
  R1 = tensordot(q1, U, Axes(2), Axes(0));
  R2 = tensordot(q2, VT, Axes(2), Axes(1));

  Tn1_new = tensordot(Q1, R1, Axes(3), Axes(0)).transpose(Axes(0, 1, 4, 2, 3));
  Tn2_new = tensordot(Q2, R2, Axes(3), Axes(0)).transpose(Axes(4, 0, 1, 2, 3));
}

template <template <typename> class Matrix, typename C>
void Full_update_bond(
    const Tensor<Matrix, C> &C1, const Tensor<Matrix, C> &C2,
    const Tensor<Matrix, C> &C3, const Tensor<Matrix, C> &C4,
    const Tensor<Matrix, C> &eT1, const Tensor<Matrix, C> &eT2,
    const Tensor<Matrix, C> &eT3, const Tensor<Matrix, C> &eT4,
    const Tensor<Matrix, C> &eT5, const Tensor<Matrix, C> &eT6,
    const Tensor<Matrix, C> &Tn1, const Tensor<Matrix, C> &Tn2,
    const Tensor<Matrix, C> &op12, const int connect1,
    const PEPS_Parameters peps_parameters, Tensor<Matrix, C> &Tn1_new,
    Tensor<Matrix, C> &Tn2_new) {
  Tensor<Matrix, C> Tn1_rot, Tn2_rot;
  if (connect1 == 0) {
    // Tn1_rot = Tn1;
    // Tn2_rot = Tn2;
    Tn1_rot = transpose(Tn1, Axes(2, 3, 0, 1, 4));
    Tn2_rot = transpose(Tn2, Axes(2, 3, 0, 1, 4));
  } else if (connect1 == 1) {
    Tn1_rot = transpose(Tn1, Axes(3, 0, 1, 2, 4));
    Tn2_rot = transpose(Tn2, Axes(3, 0, 1, 2, 4));
  } else if (connect1 == 2) {
    // Tn1_rot = transpose(Tn1,Axes(2,3,0,1,4));
    // Tn2_rot = transpose(Tn2,Axes(2,3,0,1,4));
    Tn1_rot = Tn1;
    Tn2_rot = Tn2;
  } else {
    Tn1_rot = transpose(Tn1, Axes(1, 2, 3, 0, 4));
    Tn2_rot = transpose(Tn2, Axes(1, 2, 3, 0, 4));
  }

  Full_update_bond_horizontal(C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6,
                              Tn1_rot, Tn2_rot, op12, peps_parameters, Tn1_new,
                              Tn2_new);

  if (connect1 == 0) {
    // Tn1_new = Tn1_new_rot;
    // Tn2_new = Tn2_new_rot;
    Tn1_new.transpose(Axes(2, 3, 0, 1, 4));
    Tn2_new.transpose(Axes(2, 3, 0, 1, 4));
  } else if (connect1 == 1) {
    Tn1_new.transpose(Axes(1, 2, 3, 0, 4));
    Tn2_new.transpose(Axes(1, 2, 3, 0, 4));
  } else if (connect1 == 2) {
    // Tn1_new = transpose(Tn1_new_rot,Axes(2,3,0,1,4));
    // Tn2_new = transpose(Tn2_new_rot,Axes(2,3,0,1,4));
  } else {
    Tn1_new.transpose(Axes(3, 0, 1, 2, 4));
    Tn2_new.transpose(Axes(3, 0, 1, 2, 4));
  }
}

template <template <typename> class Matrix, typename C>
void EvolutionaryTensor(Tensor<Matrix, C> &U, const Tensor<Matrix, C> &H,
                        double tau) {
  const Shape shape = H.shape();
  const size_t N = shape[0];

  Tensor<Matrix, C> V;
  std::vector<double> s;
  /* int info = */ eigh(H, s, V);
  for (size_t i = 0; i < N; ++i) {
    s[i] = std::exp(-tau * s[i]);
  }
  Tensor<Matrix, C> Vd = conj(V);
  V.multiply_vector(s, 1);
  U = tensordot(V, Vd, Axes(1), Axes(1));
}

template <template <typename> class Matrix, typename C>
Tensor<Matrix, C> EvolutionaryTensor(const Tensor<Matrix, C> &H, double tau) {
  Tensor<Matrix, C> U;
  EvolutionaryTensor(U, H, tau);
  return U;
}

template <template <typename> class Matrix, typename C>
void StartCorrelation(Tensor<Matrix, C> &A, const Tensor<Matrix, C> &C1,
                      const Tensor<Matrix, C> &C4, const Tensor<Matrix, C> &eT1,
                      const Tensor<Matrix, C> &eT3,
                      const Tensor<Matrix, C> &eT4,
                      const Tensor<Matrix, C> &Tn1,
                      const Tensor<Matrix, C> &op) {
  ////////////////////////////////////////////////////////////
  // (eT1*(Tn1*((Tn2*op)*(eT3*(C1*(C4*eT4))))))
  // cpu_cost= 7.5525e+06  memory= 192500
  // final_bond_order (e1r, e3r, n1r, n2r)
  ////////////////////////////////////////////////////////////
  A = transpose(
      tensordot(
          eT1,
          tensordot(
              Tn1,
              tensordot(
                  tensordot(conj(Tn1), op, Axes(4), Axes(1)  // Tn2 = conj(Tn1)
                            ),
                  tensordot(eT3,
                            tensordot(C1, tensordot(C4, eT4, Axes(1), Axes(0)),
                                      Axes(0), Axes(1)),
                            Axes(1), Axes(1)),
                  Axes(0, 3), Axes(5, 2)),
              Axes(0, 3, 4), Axes(6, 4, 2)),
          Axes(0, 2, 3), Axes(5, 0, 2)),
      Axes(0, 3, 1, 2));
}

template <template <typename> class Matrix, typename C>
void Transfer(Tensor<Matrix, C> &A, const Tensor<Matrix, C> &eT1,
              const Tensor<Matrix, C> &eT3, const Tensor<Matrix, C> &Tn1) {
  ////////////////////////////////////////////////////////////
  // (eT1*(Tn1*(Tn2*(A*eT3))))
  // cpu_cost= 7.5e+06  memory= 192500
  // final_bond_order (e1r, e3r, n1r, n2r)
  ////////////////////////////////////////////////////////////
  A = transpose(tensordot(eT1,
                          tensordot(Tn1,
                                    tensordot(conj(Tn1),
                                              tensordot(  // Tn2 = conj(Tn1)
                                                  A, eT3, Axes(1), Axes(1)),
                                              Axes(0, 3), Axes(2, 5)),
                                    Axes(0, 3, 4), Axes(4, 6, 2)),
                          Axes(0, 2, 3), Axes(4, 0, 2)),
                Axes(0, 3, 1, 2));
}

template <template <typename> class Matrix, typename C>
C FinishCorrelation(
    const Tensor<Matrix, C> &A, const Tensor<Matrix, C> &C2,
    const Tensor<Matrix, C> &C3, const Tensor<Matrix, C> &eT1,
    const Tensor<Matrix, C> &eT2, const Tensor<Matrix, C> &eT3,
    const Tensor<Matrix, C> &Tn1, const Tensor<Matrix, C> &op) {
  ////////////////////////////////////////////////////////////
  // (op*(Tn1*((A*eT1)*(Tn2*(eT3*(C2*(C3*eT2)))))))
  // cpu_cost= 7.5525e+06  memory= 252504
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op,
      tensordot(
          Tn1,
          tensordot(
              tensordot(A, eT1, Axes(0), Axes(0)),
              tensordot(conj(Tn1),
                        tensordot(  // Tn2 = conj(Tn1)
                            eT3,
                            tensordot(C2, tensordot(C3, eT2, Axes(0), Axes(1)),
                                      Axes(1), Axes(1)),
                            Axes(0), Axes(1)),
                        Axes(2, 3), Axes(5, 2)),
              Axes(0, 2, 3, 5), Axes(3, 0, 5, 1)),
          Axes(0, 1, 2, 3), Axes(0, 1, 4, 3)),
      Axes(0, 1), Axes(0, 1));
}

}  // end of namespace tenes

#endif  // _PEPS_BASICS_HPP_
