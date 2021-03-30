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

#include "ctm.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>
#include <mptensor/complex.hpp>
#include <mptensor/tensor.hpp>
#include <mptensor/rsvd.hpp>

#include "../../SquareLattice.hpp"
#include "../../PEPS_Parameters.hpp"
#include "../../tensor.hpp"

namespace tenes {

using mptensor::Axes;
using mptensor::Index;
using mptensor::Shape;

template <class tensor>
class Mult_col {
 public:
  Mult_col(const tensor &LT, const tensor &LB) : LT_(LT), LB_(LB){};
  tensor operator()(const tensor &T_in) {
    return tensordot(LT_, tensordot(LB_, T_in, Axes(1, 3, 5), Axes(0, 1, 2)),
                     Axes(1, 3, 5), Axes(0, 1, 2));
  }

 private:
  const tensor &LT_;
  const tensor &LB_;
};

template <class tensor>
class Mult_row {
 public:
  Mult_row(const tensor &LT, const tensor &LB) : LT_(LT), LB_(LB){};
  tensor operator()(const tensor &T_in) {
    return tensordot(tensordot(T_in, LT_, Axes(0, 1, 2), Axes(0, 2, 4)), LB_,
                     Axes(1, 2, 3), Axes(0, 2, 4));
  }

 private:
  const tensor &LT_;
  const tensor &LB_;
};

template <class tensor>
class Mult_col_ud {
 public:
  Mult_col_ud(const tensor &LT, const tensor &RT, const tensor &RB,
              const tensor &LB)
      : LT_(LT), RT_(RT), RB_(RB), LB_(LB){};
  tensor operator()(const tensor &T_in) {
    return tensordot(
        RT_,
        tensordot(
            LT_,
            tensordot(LB_, tensordot(RB_, T_in, Axes(1, 3, 5), Axes(0, 1, 2)),
                      Axes(1, 3, 5), Axes(0, 1, 2)),
            Axes(1, 3, 5), Axes(0, 1, 2)),
        Axes(0, 3, 4), Axes(0, 1, 2));
  }

 private:
  const tensor &LT_;
  const tensor &RT_;
  const tensor &RB_;
  const tensor &LB_;
};

template <class tensor>
class Mult_row_ud {
 public:
  Mult_row_ud(const tensor &LT, const tensor &RT, const tensor &RB,
              const tensor &LB)
      : LT_(LT), RT_(RT), RB_(RB), LB_(LB){};
  tensor operator()(const tensor &T_in) {
    return tensordot(
        tensordot(tensordot(tensordot(T_in, RT_, Axes(0, 1, 2), Axes(1, 2, 5)),
                            LT_, Axes(1, 2, 3), Axes(0, 2, 4)),
                  LB_, Axes(1, 2, 3), Axes(0, 2, 4)),
        RB_, Axes(1, 2, 3), Axes(0, 2, 4));
  }

 private:
  const tensor &LT_;
  const tensor &RT_;
  const tensor &RB_;
  const tensor &LB_;
};

template <class tensor>
void Calc_projector_left_block(const tensor &C1, const tensor &C4,
                               const tensor &eT1, const tensor &eT6,
                               const tensor &eT7, const tensor &eT8,
                               const tensor &Tn1, const tensor &Tn4,
                               const PEPS_Parameters peps_parameters,
                               tensor &PU, tensor &PL) {
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
    tensor LT =
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
    tensor LB =
        tensordot(Tn4,
                  tensordot(tensordot(tensordot(C4, eT7, Axes(1), Axes(0)), eT6,
                                      Axes(0), Axes(1)),
                            conj(Tn4), Axes(2, 5), Axes(0, 3)),
                  Axes(0, 3, 4), Axes(1, 3, 6));

    // tensor R1 = LT;
    // tensor R2 = LB;
    if ((t12 != 1 && t34 != 1) && peps_parameters.Use_RSVD) {
      Mult_col<tensor> m_col(LT, LB);
      Mult_row<tensor> m_row(LT, LB);

      tensor U;
      tensor VT;
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
      tensor U;
      tensor VT;
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
      tensor U_c = slice(U, 3, 0, e78);
      tensor VT_c = slice(VT, 0, 0, e78);

      U_c.multiply_vector(s_c, 3);
      VT_c.multiply_vector(s_c, 0);

      PU = tensordot(LB, conj(VT_c), Axes(1, 3, 5), Axes(1, 2, 3))
               .transpose(Axes(1, 0, 2, 3));
      PL = tensordot(LT, conj(U_c), Axes(0, 2, 4), Axes(0, 1, 2))
               .transpose(Axes(1, 0, 2, 3));
    }
  } else {
    tensor identity_matrix(Shape(e78, e78));
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

template <class tensor>
void Calc_projector_updown_blocks(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &eT7, const tensor &eT8,
    const tensor &Tn1, const tensor &Tn2, const tensor &Tn3, const tensor &Tn4,
    const PEPS_Parameters peps_parameters, tensor &PU, tensor &PL) {
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
    tensor LT =
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

    tensor RT =
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
    tensor RB =
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
    tensor LB =
        tensordot(Tn4,
                  tensordot(tensordot(tensordot(C4, eT7, Axes(1), Axes(0)), eT6,
                                      Axes(0), Axes(1)),
                            conj(Tn4), Axes(2, 5), Axes(0, 3)),
                  Axes(0, 3, 4), Axes(1, 3, 6));
    if (t23 != 1 && peps_parameters.Use_RSVD) {
      Mult_col_ud<tensor> m_col(LT, RT, RB, LB);
      Mult_row_ud<tensor> m_row(LT, RT, RB, LB);

      tensor U;
      tensor VT;
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
      tensor R1 = tensordot(RT, LT, Axes(0, 3, 4), Axes(0, 2, 4));
      tensor R2 = tensordot(RB, LB, Axes(0, 2, 4), Axes(1, 3, 5));

      tensor U;
      tensor VT;
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
      tensor U_c = slice(U, 3, 0, e78);
      tensor VT_c = slice(VT, 0, 0, e78);

      U_c.multiply_vector(s_c, 3);
      VT_c.multiply_vector(s_c, 0);

      PU = tensordot(R2, conj(VT_c), Axes(0, 1, 2), Axes(1, 2, 3))
               .transpose(Axes(1, 0, 2, 3));
      PL = tensordot(R1, conj(U_c), Axes(0, 1, 2), Axes(0, 1, 2))
               .transpose(Axes(1, 0, 2, 3));
    }
  } else {
    tensor identity_matrix(Shape(e78, e78));
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

template <class tensor>
void Calc_Next_CTM(const tensor &C1, const tensor &C4, const tensor &eT1,
                   const tensor &eT6, const tensor &PU, const tensor &PL,
                   tensor &C1_out, tensor &C4_out) {
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

template <class tensor>
void Calc_Next_eT(const tensor &eT8, const tensor &Tn1, const tensor &PU,
                  const tensor &PL, tensor &eT_out) {
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
  tensor eT_temp;
  eT_temp = contract(eT_out,Axes(2),Axes(3));
  svd(eT_temp,w);
  double norm = std::accumulate(w.begin(),w.end(),0.0);
  eT_out /= norm;
  */

  double max_all = max_abs(eT_out);
  eT_out /= max_all;
}

/*
 * corner and edge tensor
 *
 * C1 t C2
 * l  .  r
 * C4 b C3
 */

template <class tensor>
void Left_move(std::vector<tensor> &C1, const std::vector<tensor> &C2,
               const std::vector<tensor> &C3, std::vector<tensor> &C4,
               const std::vector<tensor> &eTt, const std::vector<tensor> &eTr,
               const std::vector<tensor> &eTb, std::vector<tensor> &eTl,
               const std::vector<tensor> &Tn, const int ix,
               const PEPS_Parameters peps_parameters, const SquareLattice lattice) {
  /* Do one step left move absoving X=ix column
     part of C1, C4, eTl will be modified */

  std::vector<tensor> PUs, PLs;
  PUs.resize(lattice.LY);
  PLs.resize(lattice.LY);
  int i, j, k, l;
  for (int iy = 0; iy < lattice.LY; ++iy) {
    i = lattice.index_fast(ix, iy);
    j = lattice.right(i);
    k = lattice.bottom(j);
    l = lattice.left(k);

    if (peps_parameters.CTM_Projector_corner) {
      Calc_projector_left_block(C1[i], C4[l], eTt[i], eTb[l], eTl[l], eTl[i],
                                Tn[i], Tn[l], peps_parameters, PUs[iy],
                                PLs[iy]);
    } else {
      Calc_projector_updown_blocks(C1[i], C2[j], C3[k], C4[l], eTt[i], eTt[j],
                                   eTr[j], eTr[k], eTb[k], eTb[l], eTl[l],
                                   eTl[i], Tn[i], Tn[j], Tn[k], Tn[l],
                                   peps_parameters, PUs[iy], PLs[iy]);
    }
  }
  // update
  std::vector<tensor> C1_bak(lattice.N_UNIT), C4_bak(lattice.N_UNIT),
      eTl_bak(lattice.N_UNIT);
  for (int num = 0; num < lattice.N_UNIT; num++) {
    C1_bak[num] = C1[num];
    C4_bak[num] = C4[num];
    eTl_bak[num] = eTl[num];
  }
  int iy_up, iy_down;
  for (int iy = 0; iy < lattice.LY; ++iy) {
    i = lattice.index_fast(ix, iy);
    j = lattice.right(i);
    k = lattice.bottom(j);
    l = lattice.left(k);
    iy_up = (iy + 1) % lattice.LY;
    iy_down = (iy - 1 + lattice.LY) % lattice.LY;

    Calc_Next_CTM(C1_bak[i], C4_bak[l], eTt[i], eTb[l], PUs[iy_up],
                  PLs[iy_down], C1[j], C4[k]);
    Calc_Next_eT(eTl_bak[i], Tn[i], PUs[iy], PLs[iy_up], eTl[j]);
    Calc_Next_eT(eTl_bak[l], Tn[l], PUs[iy_down], PLs[iy], eTl[k]);
  }
}

template <class tensor>
void Right_move(const std::vector<tensor> &C1, std::vector<tensor> &C2,
                std::vector<tensor> &C3, const std::vector<tensor> &C4,
                const std::vector<tensor> &eTt, std::vector<tensor> &eTr,
                const std::vector<tensor> &eTb, const std::vector<tensor> &eTl,
                const std::vector<tensor> &Tn, const int ix,
                const PEPS_Parameters peps_parameters, const SquareLattice lattice) {
  /*
    Do one step right move absorbing X=ix column
    part of C2, C3, eTr will be modified
  */
  std::vector<tensor> PUs, PLs;
  PUs.resize(lattice.LY);
  PLs.resize(lattice.LY);
  int i, j, k, l;
  for (int iy = 0; iy < lattice.LY; ++iy) {
    k = lattice.index_fast(ix, iy);
    l = lattice.left(k);
    i = lattice.top(l);
    j = lattice.right(i);

    if (peps_parameters.CTM_Projector_corner) {
      Calc_projector_left_block(C3[k], C2[j], eTb[k], eTt[j], eTr[j], eTr[k],
                                transpose(Tn[k], Axes(2, 3, 0, 1, 4)),
                                transpose(Tn[j], Axes(2, 3, 0, 1, 4)),
                                peps_parameters, PUs[iy], PLs[iy]);
    } else {
      Calc_projector_updown_blocks(
          C3[k], C4[l], C1[i], C2[j], eTb[k], eTb[l], eTl[l], eTl[i], eTt[i],
          eTt[j], eTr[j], eTr[k], transpose(Tn[k], Axes(2, 3, 0, 1, 4)),
          transpose(Tn[l], Axes(2, 3, 0, 1, 4)),
          transpose(Tn[i], Axes(2, 3, 0, 1, 4)),
          transpose(Tn[j], Axes(2, 3, 0, 1, 4)), peps_parameters, PUs[iy],
          PLs[iy]);
    }
  }
  // update
  std::vector<tensor> C2_bak(lattice.N_UNIT), C3_bak(lattice.N_UNIT),
      eTr_bak(lattice.N_UNIT);
  for (int num = 0; num < lattice.N_UNIT; num++) {
    C2_bak[num] = C2[num];
    C3_bak[num] = C3[num];
    eTr_bak[num] = eTr[num];
  }
  int iy_up, iy_down;
  for (int iy = 0; iy < lattice.LY; ++iy) {
    k = lattice.index_fast(ix, iy);
    l = lattice.left(k);
    i = lattice.top(l);
    j = lattice.right(i);

    iy_up = (iy + 1) % lattice.LY;
    iy_down = (iy - 1 + lattice.LY) % lattice.LY;

    Calc_Next_CTM(C3_bak[k], C2_bak[j], eTb[k], eTt[j], PUs[iy_down],
                  PLs[iy_up], C3[l], C2[i]);

    Calc_Next_eT(eTr_bak[k], transpose(Tn[k], Axes(2, 3, 0, 1, 4)), PUs[iy],
                 PLs[iy_down], eTr[l]);
    Calc_Next_eT(eTr_bak[j], transpose(Tn[j], Axes(2, 3, 0, 1, 4)), PUs[iy_up],
                 PLs[iy], eTr[i]);
  }
}

template <class tensor>
void Top_move(std::vector<tensor> &C1, std::vector<tensor> &C2,
              const std::vector<tensor> &C3, const std::vector<tensor> &C4,
              std::vector<tensor> &eTt, const std::vector<tensor> &eTr,
              const std::vector<tensor> &eTb, const std::vector<tensor> &eTl,
              const std::vector<tensor> &Tn, const int iy,
              const PEPS_Parameters peps_parameters, const SquareLattice lattice) {
  /*
    ## Do one step top move absorbing Y=iy row
    ## part of C1, C2, eTt will be modified
  */
  std::vector<tensor> PUs, PLs;
  PUs.resize(lattice.LX);
  PLs.resize(lattice.LX);
  int i, j, k, l;
  for (int ix = 0; ix < lattice.LX; ++ix) {
    j = lattice.index_fast(ix, iy);
    k = lattice.bottom(j);
    l = lattice.left(k);
    i = lattice.top(l);

    if (peps_parameters.CTM_Projector_corner) {
      Calc_projector_left_block(C2[j], C1[i], eTr[j], eTl[i], eTt[i], eTt[j],
                                transpose(Tn[j], Axes(1, 2, 3, 0, 4)),
                                transpose(Tn[i], Axes(1, 2, 3, 0, 4)),
                                peps_parameters, PUs[ix], PLs[ix]);
    } else {
      Calc_projector_updown_blocks(
          C2[j], C3[k], C4[l], C1[i], eTr[j], eTr[k], eTb[k], eTb[l], eTl[l],
          eTl[i], eTt[i], eTt[j], transpose(Tn[j], Axes(1, 2, 3, 0, 4)),
          transpose(Tn[k], Axes(1, 2, 3, 0, 4)),
          transpose(Tn[l], Axes(1, 2, 3, 0, 4)),
          transpose(Tn[i], Axes(1, 2, 3, 0, 4)), peps_parameters, PUs[ix],
          PLs[ix]);
    }
  }
  // update
  std::vector<tensor> C1_bak(lattice.N_UNIT), C2_bak(lattice.N_UNIT),
      eTt_bak(lattice.N_UNIT);
  for (int num = 0; num < lattice.N_UNIT; num++) {
    C1_bak[num] = C1[num];
    C2_bak[num] = C2[num];
    eTt_bak[num] = eTt[num];
  }
  int ix_right, ix_left;
  for (int ix = 0; ix < lattice.LX; ++ix) {
    j = lattice.index_fast(ix, iy);
    k = lattice.bottom(j);
    l = lattice.left(k);
    i = lattice.top(l);

    ix_right = (ix + 1) % lattice.LX;
    ix_left = (ix - 1 + lattice.LX) % lattice.LX;

    Calc_Next_CTM(C2_bak[j], C1_bak[i], eTr[j], eTl[i], PUs[ix_right],
                  PLs[ix_left], C2[k], C1[l]);

    Calc_Next_eT(eTt_bak[j], transpose(Tn[j], Axes(1, 2, 3, 0, 4)), PUs[ix],
                 PLs[ix_right], eTt[k]);
    Calc_Next_eT(eTt_bak[i], transpose(Tn[i], Axes(1, 2, 3, 0, 4)),
                 PUs[ix_left], PLs[ix], eTt[l]);
  }
}

template <class tensor>
void Bottom_move(const std::vector<tensor> &C1, const std::vector<tensor> &C2,
                 std::vector<tensor> &C3, std::vector<tensor> &C4,
                 const std::vector<tensor> &eTt, const std::vector<tensor> &eTr,
                 std::vector<tensor> &eTb, const std::vector<tensor> &eTl,
                 const std::vector<tensor> &Tn, const int iy,
                 const PEPS_Parameters peps_parameters, const SquareLattice lattice) {
  /*
    ## Do one step bottom move absorbing Y=iy row
    ## part of C3, C4, eTb will be modified
  */

  std::vector<tensor> PUs, PLs;
  PUs.resize(lattice.LX);
  PLs.resize(lattice.LX);
  int i, j, k, l;
  for (int ix = 0; ix < lattice.LX; ++ix) {
    l = lattice.index_fast(ix, iy);
    i = lattice.top(l);
    j = lattice.right(i);
    k = lattice.bottom(j);

    if (peps_parameters.CTM_Projector_corner) {
      Calc_projector_left_block(C4[l], C3[k], eTl[l], eTr[k], eTb[k], eTb[l],
                                transpose(Tn[l], Axes(3, 0, 1, 2, 4)),
                                transpose(Tn[k], Axes(3, 0, 1, 2, 4)),
                                peps_parameters, PUs[ix], PLs[ix]);
    } else {
      Calc_projector_updown_blocks(
          C4[l], C1[i], C2[j], C3[k], eTl[l], eTl[i], eTt[i], eTt[j], eTr[j],
          eTr[k], eTb[k], eTb[l], transpose(Tn[l], Axes(3, 0, 1, 2, 4)),
          transpose(Tn[i], Axes(3, 0, 1, 2, 4)),
          transpose(Tn[j], Axes(3, 0, 1, 2, 4)),
          transpose(Tn[k], Axes(3, 0, 1, 2, 4)), peps_parameters, PUs[ix],
          PLs[ix]);
    }
  }

  // update
  std::vector<tensor> C3_bak(lattice.N_UNIT), C4_bak(lattice.N_UNIT),
      eTb_bak(lattice.N_UNIT);
  for (int num = 0; num < lattice.N_UNIT; num++) {
    C3_bak[num] = C3[num];
    C4_bak[num] = C4[num];
    eTb_bak[num] = eTb[num];
  }
  int ix_left, ix_right;
  for (int ix = 0; ix < lattice.LX; ++ix) {
    l = lattice.index_fast(ix, iy);
    i = lattice.top(l);
    j = lattice.right(i);
    k = lattice.bottom(j);

    ix_right = (ix + 1) % lattice.LX;
    ix_left = (ix - 1 + lattice.LX) % lattice.LX;

    Calc_Next_CTM(C4_bak[l], C3_bak[k], eTl[l], eTr[k], PUs[ix_left],
                  PLs[ix_right], C4[i], C3[j]);

    Calc_Next_eT(eTb_bak[l], transpose(Tn[l], Axes(3, 0, 1, 2, 4)), PUs[ix],
                 PLs[ix_left], eTb[i]);
    Calc_Next_eT(eTb_bak[k], transpose(Tn[k], Axes(3, 0, 1, 2, 4)),
                 PUs[ix_right], PLs[ix], eTb[j]);
  }
}

template <class tensor>
bool Check_Convergence_CTM(
    const std::vector<tensor> &C1, const std::vector<tensor> &C2,
    const std::vector<tensor> &C3, const std::vector<tensor> &C4,
    const std::vector<tensor> &C1_old, const std::vector<tensor> &C2_old,
    const std::vector<tensor> &C3_old, const std::vector<tensor> &C4_old,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice,
    double &sig_max) {
  sig_max = 0.0;
  bool convergence = true;
  double sig, norm;
  std::vector<double> lam_new, lam_old;

  for (int i = 0; i < lattice.N_UNIT; ++i) {
    // C1
    svd(C1[i], lam_new);
    norm = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      norm += lam_new[j] * lam_new[j];
    }
    norm = sqrt(norm);
    for (int k = 0; k < lam_new.size(); ++k) {
      lam_new[k] /= norm;
    }
    svd(C1_old[i], lam_old);
    norm = 0.0;
    for (int j = 0; j < lam_old.size(); ++j) {
      norm += lam_old[j] * lam_old[j];
    }
    norm = sqrt(norm);
    for (int k = 0; k < lam_old.size(); ++k) {
      lam_old[k] /= norm;
    }

    sig = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      sig += (lam_new[j] - lam_old[j]) * (lam_new[j] - lam_old[j]);
    }
    sig = sqrt(sig);

    if (sig > peps_parameters.CTM_Convergence_Epsilon) {
      sig_max = sig;
      convergence = false;
      break;
    } else if (sig > sig_max) {
      sig_max = sig;
    }

    // C2
    svd(C2[i], lam_new);
    norm = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      norm += lam_new[j] * lam_new[j];
    }
    norm = sqrt(norm);
    for (int k = 0; k < lam_new.size(); ++k) {
      lam_new[k] /= norm;
    }

    svd(C2_old[i], lam_old);
    norm = 0.0;
    for (int j = 0; j < lam_old.size(); ++j) {
      norm += lam_old[j] * lam_old[j];
    }

    norm = sqrt(norm);
    for (int k = 0; k < lam_old.size(); ++k) {
      lam_old[k] /= norm;
    }

    sig = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      sig += (lam_new[j] - lam_old[j]) * (lam_new[j] - lam_old[j]);
    }
    sig = sqrt(sig);

    if (sig > peps_parameters.CTM_Convergence_Epsilon) {
      sig_max = sig;
      convergence = false;
      break;
    } else if (sig > sig_max) {
      sig_max = sig;
    }

    // C3
    svd(C3[i], lam_new);
    norm = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      norm += lam_new[j] * lam_new[j];
    }
    norm = sqrt(norm);
    for (int k = 0; k < lam_new.size(); ++k) {
      lam_new[k] /= norm;
    }
    svd(C3_old[i], lam_old);
    norm = 0.0;
    for (int j = 0; j < lam_old.size(); ++j) {
      norm += lam_old[j] * lam_old[j];
    }

    norm = sqrt(norm);
    for (int k = 0; k < lam_old.size(); ++k) {
      lam_old[k] /= norm;
    }

    sig = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      sig += (lam_new[j] - lam_old[j]) * (lam_new[j] - lam_old[j]);
    }
    sig = sqrt(sig);

    if (sig > peps_parameters.CTM_Convergence_Epsilon) {
      sig_max = sig;
      convergence = false;
      break;
    } else if (sig > sig_max) {
      sig_max = sig;
    }

    // C4
    svd(C4[i], lam_new);
    norm = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      norm += lam_new[j] * lam_new[j];
    }
    norm = sqrt(norm);
    for (int k = 0; k < lam_new.size(); ++k) {
      lam_new[k] /= norm;
    }
    svd(C4_old[i], lam_old);
    norm = 0.0;
    for (int j = 0; j < lam_old.size(); ++j) {
      norm += lam_old[j] * lam_old[j];
    }
    norm = sqrt(norm);
    for (int k = 0; k < lam_old.size(); ++k) {
      lam_old[k] /= norm;
    }

    sig = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      sig += (lam_new[j] - lam_old[j]) * (lam_new[j] - lam_old[j]);
    }
    sig = sqrt(sig);

    if (sig > peps_parameters.CTM_Convergence_Epsilon) {
      sig_max = sig;
      convergence = false;
      break;
    } else if (sig > sig_max) {
      sig_max = sig;
    }
  }
  return convergence;
}

template <class tensor>
int Calc_CTM_Environment(std::vector<tensor> &C1, std::vector<tensor> &C2,
                         std::vector<tensor> &C3, std::vector<tensor> &C4,
                         std::vector<tensor> &eTt, std::vector<tensor> &eTr,
                         std::vector<tensor> &eTb, std::vector<tensor> &eTl,
                         const std::vector<tensor> &Tn,
                         const PEPS_Parameters peps_parameters,
                         const SquareLattice lattice, bool initialize) {
  /*
    ## Calc environment tensors
    ## C1,C2,C3,C4 and eTt,eTl,eTr,eTb will be modified
  */
  // Initialize
  if (initialize) {
    int num, d1, d2, d34;
    Index index;
    tensor Projector;

    for (int i = 0; i < lattice.N_UNIT; ++i) {
      num = lattice.top(lattice.left(i));
      d1 = Tn[num].shape()[3] * Tn[num].shape()[3];
      d2 = Tn[num].shape()[2] * Tn[num].shape()[2];
      C1[i] = reshape(
          tensordot(Tn[num], conj(Tn[num]), Axes(0, 1, 4), Axes(0, 1, 4))
              .transpose(Axes(1, 3, 0, 2)),
          (Shape(d1, d2)));

      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        C1[i] = extend(C1[i], Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        C1[i] = slice(slice(C1[i], 0, 0, peps_parameters.CHI), 1, 0,
                      peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        C1[i] = extend(slice(C1[i], 1, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else {
        // d1 >= CHI
        // d2 < CHI
        C1[i] = extend(slice(C1[i], 0, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      }

      num = lattice.top(lattice.right(i));
      d1 = Tn[num].shape()[0] * Tn[num].shape()[0];
      d2 = Tn[num].shape()[3] * Tn[num].shape()[3];
      C2[i] = reshape(
          tensordot(Tn[num], conj(Tn[num]), Axes(1, 2, 4), Axes(1, 2, 4))
              .transpose(Axes(0, 2, 1, 3)),
          (Shape(d1, d2)));

      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        C2[i] = extend(C2[i], Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        C2[i] = slice(slice(C2[i], 0, 0, peps_parameters.CHI), 1, 0,
                      peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        C2[i] = extend(slice(C2[i], 1, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else {
        // d1 >= CHI
        // d2 < CHI
        C2[i] = extend(slice(C2[i], 0, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      }

      num = lattice.bottom(lattice.right(i));
      d1 = Tn[num].shape()[1] * Tn[num].shape()[1];
      d2 = Tn[num].shape()[0] * Tn[num].shape()[0];
      C3[i] = reshape(
          tensordot(Tn[num], conj(Tn[num]), Axes(2, 3, 4), Axes(2, 3, 4))
              .transpose(Axes(1, 3, 0, 2)),
          (Shape(d1, d2)));

      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        C3[i] = extend(C3[i], Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        C3[i] = slice(slice(C3[i], 0, 0, peps_parameters.CHI), 1, 0,
                      peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        C3[i] = extend(slice(C3[i], 1, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else {
        // d1 >= CHI
        // d2 < CHI
        C3[i] = extend(slice(C3[i], 0, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      }

      num = lattice.bottom(lattice.left(i));
      d1 = Tn[num].shape()[2] * Tn[num].shape()[2];
      d2 = Tn[num].shape()[1] * Tn[num].shape()[1];
      C4[i] = reshape(
          tensordot(Tn[num], conj(Tn[num]), Axes(0, 3, 4), Axes(0, 3, 4))
              .transpose(Axes(1, 3, 0, 2)),
          (Shape(d1, d2)));

      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        C4[i] = extend(C4[i], Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        C4[i] = slice(slice(C4[i], 0, 0, peps_parameters.CHI), 1, 0,
                      peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        C4[i] = extend(slice(C4[i], 1, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else {
        // d1 >= CHI
        // d2 < CHI
        C4[i] = extend(slice(C4[i], 0, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      }

      num = lattice.top(i);
      d1 = Tn[num].shape()[0] * Tn[num].shape()[0];
      d2 = Tn[num].shape()[2] * Tn[num].shape()[2];
      d34 = Tn[num].shape()[3];
      eTt[i] = reshape(tensordot(Tn[num], conj(Tn[num]), Axes(1, 4), Axes(1, 4))
                           .transpose(Axes(0, 3, 1, 4, 2, 5)),
                       (Shape(d1, d2, d34, d34)));

      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        eTt[i] = extend(
            eTt[i], Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        eTt[i] = slice(slice(eTt[i], 0, 0, peps_parameters.CHI), 1, 0,
                       peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        eTt[i] =
            extend(slice(eTt[i], 1, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else {
        // d1 >= CHI
        // d2 < CHI
        eTt[i] =
            extend(slice(eTt[i], 0, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      }

      num = lattice.right(i);
      d1 = Tn[num].shape()[1] * Tn[num].shape()[1];
      d2 = Tn[num].shape()[3] * Tn[num].shape()[3];
      d34 = Tn[num].shape()[0];
      eTr[i] = reshape(tensordot(Tn[num], conj(Tn[num]), Axes(2, 4), Axes(2, 4))
                           .transpose(Axes(1, 4, 2, 5, 0, 3)),
                       (Shape(d1, d2, d34, d34)));

      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        eTr[i] = extend(
            eTr[i], Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        eTr[i] = slice(slice(eTr[i], 0, 0, peps_parameters.CHI), 1, 0,
                       peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        eTr[i] =
            extend(slice(eTr[i], 1, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else {
        // d1 >= CHI
        // d2 < CHI
        eTr[i] =
            extend(slice(eTr[i], 0, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      }

      num = lattice.bottom(i);
      d1 = Tn[num].shape()[2] * Tn[num].shape()[2];
      d2 = Tn[num].shape()[0] * Tn[num].shape()[0];
      d34 = Tn[num].shape()[1];
      eTb[i] = reshape(tensordot(Tn[num], conj(Tn[num]), Axes(3, 4), Axes(3, 4))
                           .transpose(Axes(2, 5, 0, 3, 1, 4)),
                       (Shape(d1, d2, d34, d34)));

      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        eTb[i] = extend(
            eTb[i], Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        eTb[i] = slice(slice(eTb[i], 0, 0, peps_parameters.CHI), 1, 0,
                       peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        eTb[i] =
            extend(slice(eTb[i], 1, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else {
        // d1 >= CHI
        // d2 < CHI
        eTb[i] =
            extend(slice(eTb[i], 0, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      }

      num = lattice.left(i);
      d1 = Tn[num].shape()[3] * Tn[num].shape()[3];
      d2 = Tn[num].shape()[1] * Tn[num].shape()[1];
      d34 = Tn[num].shape()[2];
      eTl[i] = reshape(tensordot(Tn[num], conj(Tn[num]), Axes(0, 4), Axes(0, 4))
                           .transpose(Axes(2, 5, 0, 3, 1, 4)),
                       (Shape(d1, d2, d34, d34)));
      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        eTl[i] = extend(
            eTl[i], Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        eTl[i] = slice(slice(eTl[i], 0, 0, peps_parameters.CHI), 1, 0,
                       peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        eTl[i] =
            extend(slice(eTl[i], 1, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else {
        // d1 >= CHI
        // d2 < CHI
        eTl[i] =
            extend(slice(eTl[i], 0, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      }
    }
  }
  // Initialize done

  bool convergence = false;
  int count = 0;
  std::vector<tensor> C1_old = C1;
  std::vector<tensor> C2_old = C2;
  std::vector<tensor> C3_old = C3;
  std::vector<tensor> C4_old = C4;

  double sig_max = 0.0;
  while ((!convergence) && (count < peps_parameters.Max_CTM_Iteration)) {
    // left move
    for (int ix = 0; ix < lattice.LX; ++ix) {
      Left_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, ix, peps_parameters,
                lattice);
    }

    // right move
    for (int ix = 0; ix > -lattice.LX; --ix) {
      Right_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
                 (ix + 1 + lattice.LX) % lattice.LX, peps_parameters, lattice);
    }

    // top move
    for (int iy = 0; iy > -lattice.LY; --iy) {
      Top_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
               (iy + 1 + lattice.LY) % lattice.LY, peps_parameters, lattice);
    }

    // bottom move

    for (int iy = 0; iy < lattice.LY; ++iy) {
      Bottom_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, iy, peps_parameters,
                  lattice);
    }

    convergence =
        Check_Convergence_CTM(C1, C2, C3, C4, C1_old, C2_old, C3_old, C4_old,
                              peps_parameters, lattice, sig_max);
    count += 1;

    C1_old = C1;
    C2_old = C2;
    C3_old = C3;
    C4_old = C4;
    if (peps_parameters.print_level >= PrintLevel::debug) {
      std::cout << "CTM: count, sig_max " << count << " " << sig_max
                << std::endl;
    }
  }

  if (!convergence && peps_parameters.print_level >= PrintLevel::warn) {
    std::cout << "Warning: CTM did not converge! count, sig_max = " << count
              << " " << sig_max << std::endl;
  }
  if (peps_parameters.print_level >= PrintLevel::debug) {
    std::cout << "CTM: count to convergence= " << count << std::endl;
  }
  return count;
}

// template instantiate

template void Calc_projector_left_block(
    const real_tensor &C1, const real_tensor &C4, const real_tensor &eT1,
    const real_tensor &eT6, const real_tensor &eT7, const real_tensor &eT8,
    const real_tensor &Tn1, const real_tensor &Tn4,
    const PEPS_Parameters peps_parameters, real_tensor &PU, real_tensor &PL);
template void Calc_projector_left_block(
    const complex_tensor &C1, const complex_tensor &C4,
    const complex_tensor &eT1, const complex_tensor &eT6,
    const complex_tensor &eT7, const complex_tensor &eT8,
    const complex_tensor &Tn1, const complex_tensor &Tn4,
    const PEPS_Parameters peps_parameters, complex_tensor &PU,
    complex_tensor &PL);

template void Calc_projector_updown_blocks(
    const real_tensor &C1, const real_tensor &C2, const real_tensor &C3,
    const real_tensor &C4, const real_tensor &eT1, const real_tensor &eT2,
    const real_tensor &eT3, const real_tensor &eT4, const real_tensor &eT5,
    const real_tensor &eT6, const real_tensor &eT7, const real_tensor &eT8,
    const real_tensor &Tn1, const real_tensor &Tn2, const real_tensor &Tn3,
    const real_tensor &Tn4, const PEPS_Parameters peps_parameters,
    real_tensor &PU, real_tensor &PL);
template void Calc_projector_updown_blocks(
    const complex_tensor &C1, const complex_tensor &C2,
    const complex_tensor &C3, const complex_tensor &C4,
    const complex_tensor &eT1, const complex_tensor &eT2,
    const complex_tensor &eT3, const complex_tensor &eT4,
    const complex_tensor &eT5, const complex_tensor &eT6,
    const complex_tensor &eT7, const complex_tensor &eT8,
    const complex_tensor &Tn1, const complex_tensor &Tn2,
    const complex_tensor &Tn3, const complex_tensor &Tn4,
    const PEPS_Parameters peps_parameters, complex_tensor &PU,
    complex_tensor &PL);

template void Calc_Next_CTM(const real_tensor &C1, const real_tensor &C4,
                            const real_tensor &eT1, const real_tensor &eT6,
                            const real_tensor &PU, const real_tensor &PL,
                            real_tensor &C1_out, real_tensor &C4_out);
template void Calc_Next_CTM(const complex_tensor &C1, const complex_tensor &C4,
                            const complex_tensor &eT1,
                            const complex_tensor &eT6, const complex_tensor &PU,
                            const complex_tensor &PL, complex_tensor &C1_out,
                            complex_tensor &C4_out);

template void Calc_Next_eT(const real_tensor &eT8, const real_tensor &Tn1,
                           const real_tensor &PU, const real_tensor &PL,
                           real_tensor &eT_out);
template void Calc_Next_eT(const complex_tensor &eT8, const complex_tensor &Tn1,
                           const complex_tensor &PU, const complex_tensor &PL,
                           complex_tensor &eT_out);

template void Left_move(
    std::vector<real_tensor> &C1, const std::vector<real_tensor> &C2,
    const std::vector<real_tensor> &C3, std::vector<real_tensor> &C4,
    const std::vector<real_tensor> &eTt, const std::vector<real_tensor> &eTr,
    const std::vector<real_tensor> &eTb, std::vector<real_tensor> &eTl,
    const std::vector<real_tensor> &Tn, const int ix,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice);
template void Left_move(
    std::vector<complex_tensor> &C1, const std::vector<complex_tensor> &C2,
    const std::vector<complex_tensor> &C3, std::vector<complex_tensor> &C4,
    const std::vector<complex_tensor> &eTt,
    const std::vector<complex_tensor> &eTr,
    const std::vector<complex_tensor> &eTb, std::vector<complex_tensor> &eTl,
    const std::vector<complex_tensor> &Tn, const int ix,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice);

template void Right_move(
    const std::vector<real_tensor> &C1, std::vector<real_tensor> &C2,
    std::vector<real_tensor> &C3, const std::vector<real_tensor> &C4,
    const std::vector<real_tensor> &eTt, std::vector<real_tensor> &eTr,
    const std::vector<real_tensor> &eTb, const std::vector<real_tensor> &eTl,
    const std::vector<real_tensor> &Tn, const int ix,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice);
template void Right_move(
    const std::vector<complex_tensor> &C1, std::vector<complex_tensor> &C2,
    std::vector<complex_tensor> &C3, const std::vector<complex_tensor> &C4,
    const std::vector<complex_tensor> &eTt, std::vector<complex_tensor> &eTr,
    const std::vector<complex_tensor> &eTb,
    const std::vector<complex_tensor> &eTl,
    const std::vector<complex_tensor> &Tn, const int ix,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice);

template void Top_move(
    std::vector<real_tensor> &C1, std::vector<real_tensor> &C2,
    const std::vector<real_tensor> &C3, const std::vector<real_tensor> &C4,
    std::vector<real_tensor> &eTt, const std::vector<real_tensor> &eTr,
    const std::vector<real_tensor> &eTb, const std::vector<real_tensor> &eTl,
    const std::vector<real_tensor> &Tn, const int iy,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice);
template void Top_move(std::vector<complex_tensor> &C1,
                       std::vector<complex_tensor> &C2,
                       const std::vector<complex_tensor> &C3,
                       const std::vector<complex_tensor> &C4,
                       std::vector<complex_tensor> &eTt,
                       const std::vector<complex_tensor> &eTr,
                       const std::vector<complex_tensor> &eTb,
                       const std::vector<complex_tensor> &eTl,
                       const std::vector<complex_tensor> &Tn, const int iy,
                       const PEPS_Parameters peps_parameters,
                       const SquareLattice lattice);

template void Bottom_move(
    const std::vector<real_tensor> &C1, const std::vector<real_tensor> &C2,
    std::vector<real_tensor> &C3, std::vector<real_tensor> &C4,
    const std::vector<real_tensor> &eTt, const std::vector<real_tensor> &eTr,
    std::vector<real_tensor> &eTb, const std::vector<real_tensor> &eTl,
    const std::vector<real_tensor> &Tn, const int iy,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice);
template void Bottom_move(
    const std::vector<complex_tensor> &C1,
    const std::vector<complex_tensor> &C2, std::vector<complex_tensor> &C3,
    std::vector<complex_tensor> &C4, const std::vector<complex_tensor> &eTt,
    const std::vector<complex_tensor> &eTr, std::vector<complex_tensor> &eTb,
    const std::vector<complex_tensor> &eTl,
    const std::vector<complex_tensor> &Tn, const int iy,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice);

template bool Check_Convergence_CTM(const std::vector<real_tensor> &C1,
                                    const std::vector<real_tensor> &C2,
                                    const std::vector<real_tensor> &C3,
                                    const std::vector<real_tensor> &C4,
                                    const std::vector<real_tensor> &C1_old,
                                    const std::vector<real_tensor> &C2_old,
                                    const std::vector<real_tensor> &C3_old,
                                    const std::vector<real_tensor> &C4_old,
                                    const PEPS_Parameters peps_parameters,
                                    const SquareLattice lattice, double &sig_max);
template bool Check_Convergence_CTM(const std::vector<complex_tensor> &C1,
                                    const std::vector<complex_tensor> &C2,
                                    const std::vector<complex_tensor> &C3,
                                    const std::vector<complex_tensor> &C4,
                                    const std::vector<complex_tensor> &C1_old,
                                    const std::vector<complex_tensor> &C2_old,
                                    const std::vector<complex_tensor> &C3_old,
                                    const std::vector<complex_tensor> &C4_old,
                                    const PEPS_Parameters peps_parameters,
                                    const SquareLattice lattice, double &sig_max);

template int Calc_CTM_Environment(
    std::vector<real_tensor> &C1, std::vector<real_tensor> &C2,
    std::vector<real_tensor> &C3, std::vector<real_tensor> &C4,
    std::vector<real_tensor> &eTt, std::vector<real_tensor> &eTr,
    std::vector<real_tensor> &eTb, std::vector<real_tensor> &eTl,
    const std::vector<real_tensor> &Tn, const PEPS_Parameters peps_parameters,
    const SquareLattice lattice, bool initialize);

template int Calc_CTM_Environment(
    std::vector<complex_tensor> &C1, std::vector<complex_tensor> &C2,
    std::vector<complex_tensor> &C3, std::vector<complex_tensor> &C4,
    std::vector<complex_tensor> &eTt, std::vector<complex_tensor> &eTr,
    std::vector<complex_tensor> &eTb, std::vector<complex_tensor> &eTl,
    const std::vector<complex_tensor> &Tn,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice,
    bool initialize);

}  // end of namespace tenes
