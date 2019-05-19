#ifndef _PEPS_BASICS_HPP_
#define _PEPS_BASICS_HPP_

/*
 Basic routines independent on unit cell structures.
 Using mptensor libraries
 (Test version)
 2015 Dec.  Tsuyoshi Okubo
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <mptensor/complex.hpp>
#include <mptensor/rsvd.hpp>
#include <mptensor/tensor.hpp>
#include <numeric>
#include <vector>

#include "PEPS_Parameters.hpp"
using namespace mptensor;
// Contractions
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
C Contract_two_sites_holizontal(
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
    # ((Tn1*(((C1*eT1)*(C4*eT6))*(Tn1c*(((C2*eT2)*(Tn2*(((C3*eT4)*eT3)*(Tn2c*op2))))*eT5))))*op1)
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
  # ((Tn1*(((C2*eT2)*(C1*eT1))*(Tn1c*(((C3*eT3)*(Tn2*(((C4*eT5)*eT4)*(Tn2c*op2))))*eT6))))*op1)
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
C Contract_two_sites_holizontal_op12(
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
  // two_sites_holizontal_op12.dat
  ////////////////////////////////////////////////////////////
  // (op12*((eT1*(Tn1*(Tn1c*(eT5*(C1*(C4*eT6))))))*(eT2*(Tn2*(Tn2c*(eT4*(C2*(C3*eT3))))))))
  // cpu_cost= 2.20416e+11  memory= 6.0502e+08
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op12,
      tensordot(
          tensordot(eT1,
                    tensordot(Tn1,
                              tensordot(conj(Tn1),
                                        tensordot(eT5,
                                                  tensordot(C1,
                                                            tensordot(C4, eT6,
                                                                      Axes(1),
                                                                      Axes(0)),
                                                            Axes(0), Axes(1)),
                                                  Axes(1), Axes(1)),
                                        Axes(0, 3), Axes(5, 2)),
                              Axes(0, 3), Axes(6, 4)),
                    Axes(0, 2, 3), Axes(7, 0, 3)),
          tensordot(eT2,
                    tensordot(Tn2,
                              tensordot(conj(Tn2),
                                        tensordot(eT4,
                                                  tensordot(C2,
                                                            tensordot(C3, eT3,
                                                                      Axes(0),
                                                                      Axes(1)),
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
          tensordot(eT2,
                    tensordot(Tn1,
                              tensordot(conj(Tn1),
                                        tensordot(eT6,
                                                  tensordot(C1,
                                                            tensordot(C2, eT1,
                                                                      Axes(0),
                                                                      Axes(1)),
                                                            Axes(1), Axes(1)),
                                                  Axes(1), Axes(0)),
                                        Axes(0, 1), Axes(2, 5)),
                              Axes(0, 1), Axes(4, 6)),
                    Axes(0, 2, 3), Axes(7, 0, 3)),
          tensordot(eT3,
                    tensordot(Tn2,
                              tensordot(conj(Tn2),
                                        tensordot(eT5,
                                                  tensordot(C3,
                                                            tensordot(C4, eT4,
                                                                      Axes(0),
                                                                      Axes(1)),
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
// environment

template <template <typename> class Matrix, typename C> class Mult_col {
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

template <template <typename> class Matrix, typename C> class Mult_row {
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
template <template <typename> class Matrix, typename C> class Mult_col_ud {
public:
  Mult_col_ud(const Tensor<Matrix, C> &LT, const Tensor<Matrix, C> &RT,
              const Tensor<Matrix, C> &RB, const Tensor<Matrix, C> &LB)
      : LT_(LT), RT_(RT), RB_(RB), LB_(LB){};
  Tensor<Matrix, C> operator()(const Tensor<Matrix, C> &T_in) {
    return tensordot(
        RT_,
        tensordot(LT_,
                  tensordot(LB_,
                            tensordot(RB_, T_in, Axes(1, 3, 5), Axes(0, 1, 2)),
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

template <template <typename> class Matrix, typename C> class Mult_row_ud {
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

      int info = rsvd(m_row, m_col, shape_row, shape_col, U, s, VT, e78,
                      peps_parameters.RSVD_Oversampling_factor * e78);
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

      int info = svd(tensordot(LT, LB, Axes(1, 3, 5), Axes(0, 2, 4)),
                     Axes(0, 1, 2), Axes(3, 4, 5), U, s, VT);
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

  int e12 = eT1.shape()[1];
  int e34 = eT3.shape()[1];
  int e56 = eT5.shape()[1];
  int e78 = eT7.shape()[1];
  int t12 = Tn1.shape()[2];
  int t23 = Tn2.shape()[3];
  int t34 = Tn3.shape()[0];
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

      int info = rsvd(m_row, m_col, shape_row, shape_col, U, s, VT, e78,
                      peps_parameters.RSVD_Oversampling_factor * e78);
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

      int info = svd(tensordot(R1, R2, Axes(3, 4, 5), Axes(3, 4, 5)),
                     Axes(0, 1, 2), Axes(3, 4, 5), U, s, VT);
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
                        const std::vector<std::vector<double> > &lambda1,
                        const std::vector<std::vector<double> > &lambda2,
                        const Tensor<Matrix, C> &op12, const int connect1,
                        const PEPS_Parameters peps_parameters,
                        Tensor<Matrix, C> &Tn1_new, Tensor<Matrix, C> &Tn2_new,
                        std::vector<double> &lambda_c) {
  // --前処理--
  int connect2;
  if (connect1 == 0) {
    connect2 = 2;
  } else if (connect1 == 1) {
    connect2 = 3;
  } else if (connect1 == 2) {
    connect2 = 0;
  } else {
    connect2 = 1;
  }

  std::vector<std::vector<double> > lambda1_inv(4);
  std::vector<std::vector<double> > lambda2_inv(4);

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

  int m1 = Tn1.shape()[4];
  int m2 = Tn2.shape()[4];
  int dc = Tn1.shape()[connect1];
  Tensor<Matrix, C> Tn1_lambda, Tn2_lambda;
  Tn1_lambda = Tn1;
  Tn2_lambda = Tn2;

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
  Tensor<Matrix, C> Q1, R1, Q2, R2;
  int info = qr(Tn1_lambda, Axes(0, 1, 2), Axes(3, 4), Q1, R1);

  info = qr(Tn2_lambda, Axes(0, 1, 2), Axes(3, 4), Q2, R2);
  // connetc R1, R2, op

  /*
    INFO:8 (1,2) Finish 7/8 script=[0, 1, -1, 2, -1]
    ##############################
    # ((R1*R2)*op12)
    # cpu_cost= 22400  memory= 3216
    # final_bond_order  (c1, c2, m1o, m2o)
    ##############################
  */

  Tensor<Matrix, C> Theta = tensordot(tensordot(R1, R2, Axes(1), Axes(1)), op12,
                                      Axes(1, 3), Axes(0, 1));

  // svd
  Tensor<Matrix, C> U, VT;
  std::vector<double> s;
  info = svd(Theta, Axes(0, 2), Axes(1, 3), U, s, VT);

  lambda_c = std::vector<double>(&s[0], &s[dc]);
  Tensor<Matrix, C> Uc = slice(U, 2, 0, dc);
  Tensor<Matrix, C> VTc = slice(VT, 0, 0, dc);

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

  Uc.multiply_vector(lambda_c, 2);
  VTc.multiply_vector(lambda_c, 0);

  // Create new tensors
  // remove lambda effets
  if (connect1 == 0) {
    Q1.multiply_vector(lambda1_inv[1], 0, lambda1_inv[2], 1, lambda1_inv[3], 2);
  } else if (connect1 == 1) {
    Q1.multiply_vector(lambda1_inv[0], 0, lambda1_inv[2], 1, lambda1_inv[3], 2);
  } else if (connect1 == 2) {
    Q1.multiply_vector(lambda1_inv[0], 0, lambda1_inv[1], 1, lambda1_inv[3], 2);
  } else {
    Q1.multiply_vector(lambda1_inv[0], 0, lambda1_inv[1], 1, lambda1_inv[2], 2);
  };

  if (connect2 == 0) {
    Q2.multiply_vector(lambda2_inv[1], 0, lambda2_inv[2], 1, lambda2_inv[3], 2);
  } else if (connect2 == 1) {
    Q2.multiply_vector(lambda2_inv[0], 0, lambda2_inv[2], 1, lambda2_inv[3], 2);
  } else if (connect2 == 2) {
    Q2.multiply_vector(lambda2_inv[0], 0, lambda2_inv[1], 1, lambda2_inv[3], 2);
  } else {
    Q2.multiply_vector(lambda2_inv[0], 0, lambda2_inv[1], 1, lambda2_inv[2], 2);
  };

  if (connect1 == 0) {
    Tn1_new =
        tensordot(Q1, Uc, Axes(3), Axes(0)).transpose(Axes(4, 0, 1, 2, 3));
  } else if (connect1 == 1) {
    Tn1_new =
        tensordot(Q1, Uc, Axes(3), Axes(0)).transpose(Axes(0, 4, 1, 2, 3));
  } else if (connect1 == 2) {
    Tn1_new =
        tensordot(Q1, Uc, Axes(3), Axes(0)).transpose(Axes(0, 1, 4, 2, 3));
  } else {
    Tn1_new =
        tensordot(Q1, Uc, Axes(3), Axes(0)).transpose(Axes(0, 1, 2, 4, 3));
  };

  if (connect2 == 0) {
    Tn2_new =
        tensordot(Q2, VTc, Axes(3), Axes(1)).transpose(Axes(3, 0, 1, 2, 4));
  } else if (connect2 == 1) {
    Tn2_new =
        tensordot(Q2, VTc, Axes(3), Axes(1)).transpose(Axes(0, 3, 1, 2, 4));
  } else if (connect2 == 2) {
    Tn2_new =
        tensordot(Q2, VTc, Axes(3), Axes(1)).transpose(Axes(0, 1, 3, 2, 4));
  } else {
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
  int m1 = Tn1_shape[4];
  int m2 = Tn2_shape[4];

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
  R1 = transpose(U, Axes(0, 2, 1));  // envR1 , D_connect, m1
  R2 = transpose(VT, Axes(1, 0, 2)); // envR2 , D_connect, m2

  int count = 0;
  C C_phi, Old_delta, delta;
  C_phi = trace(tensordot(Environment, Theta, Axes(0, 1), Axes(0, 1)),
                conj(Theta), Axes(0, 1, 2, 3), Axes(0, 1, 2, 3));
  Old_delta =
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

    W_vec = tensordot(
        tensordot(Environment, Theta, Axes(0, 1), Axes(0, 1)), conj(R2),
        Axes(1, 3), Axes(0, 2)); //transpose(0,2,1)).reshape(envR1*D_connect,m1)
    /*
      ## create N
      ## (envR1, envR1*, D_connect,D_connect*)
    */
    N_mat = tensordot(tensordot(Environment, R2, Axes(1), Axes(0)), conj(R2),
                      Axes(2, 4), Axes(0, 2));
    //transpose(1,3,0,2).reshape(envR1*D_connect,envR1*D_connect)

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
                  Axes(0, 2)); //).transpose(0,2,1).reshape(envR2*D_connect,m2)

    /*
      ## create N
      ## (envR2, envR2*,D_connect, D_connect*)
    */
    N_mat = tensordot(
        tensordot(Environment, R1, Axes(0), Axes(0)), conj(R1), Axes(1, 4),
        Axes(0,
             2)); //.transpose(1,3,0,2).reshape(envR2*D_connect,envR2*D_connect)

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

    delta = -2.0 * trace(tensordot(tensordot(Environment, Theta, Axes(0, 1),
                                             Axes(0, 1)),
                                   conj(R2), Axes(1, 3), Axes(0, 2)),
                         conj(R1), Axes(0, 1, 2), Axes(0, 2, 1)) +
            trace(R1,
                  tensordot(
                      R2,
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
  if (!convergence && peps_parameters.Warning_flag) {
    std::cout << "warning: Full update iteration was not conveged! count= "
              << count << std::endl;
  }
  if (peps_parameters.Debug_flag) {
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
#endif // _PEPS_BASICS_HPP_
