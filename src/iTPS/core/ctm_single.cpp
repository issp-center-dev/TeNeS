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
#include "../../tensor.hpp"
#include "../PEPS_Parameters.hpp"

namespace tenes {
namespace itps {
namespace core {

using mptensor::Axes;
using mptensor::Index;
using mptensor::Shape;

template <class tensor>
std::vector<tensor> Make_single_tensor_density(const std::vector<tensor> &Tn){
  std::vector<tensor> Tn_single;
  for (int num = 0; num < Tn.size(); ++num) {
    Tn_single.push_back(mptensor::contract(Tn[num], Axes(4), Axes(5)));
  }
  return Tn_single;
}
  
template <class tensor>
class Mult_col_single {
 public:
  Mult_col_single(const tensor &LT, const tensor &LB) : LT_(LT), LB_(LB){};
  tensor operator()(const tensor &T_in) {
    return tensordot(LT_, tensordot(LB_, T_in, Axes(1), Axes(0)),
                     Axes(0), Axes(0));
  }

 private:
  const tensor &LT_;
  const tensor &LB_;
};

template <class tensor>
class Mult_row_single {
 public:
  Mult_row_single(const tensor &LT, const tensor &LB) : LT_(LT), LB_(LB){};
  tensor operator()(const tensor &T_in) {
    return tensordot(tensordot(T_in, LT_, Axes(0), Axes(1)), LB_,
                     Axes(1), Axes(0));
  }

 private:
  const tensor &LT_;
  const tensor &LB_;
};

template <class tensor>
class Mult_col_ud_single {
 public:
  Mult_col_ud_single(const tensor &LT, const tensor &RT, const tensor &RB,
              const tensor &LB)
      : LT_(LT), RT_(RT), RB_(RB), LB_(LB){};
  tensor operator()(const tensor &T_in) {
    return tensordot(
        RT_,
        tensordot(
            LT_,
            tensordot(LB_, tensordot(RB_, T_in, Axes(1), Axes(0)),
                      Axes(1), Axes(0)),
            Axes(0), Axes(0)),
        Axes(0), Axes(0));
  }

 private:
  const tensor &LT_;
  const tensor &RT_;
  const tensor &RB_;
  const tensor &LB_;
};

template <class tensor>
class Mult_row_ud_single {
 public:
  Mult_row_ud_single(const tensor &LT, const tensor &RT, const tensor &RB,
              const tensor &LB)
      : LT_(LT), RT_(RT), RB_(RB), LB_(LB){};
  tensor operator()(const tensor &T_in) {
    return tensordot(
        tensordot(tensordot(tensordot(T_in, RT_, Axes(0), Axes(1)),
                            LT_, Axes(1), Axes(1)),
                  LB_, Axes(1), Axes(0)),
        RB_, Axes(1), Axes(0));
  }

 private:
  const tensor &LT_;
  const tensor &RT_;
  const tensor &RB_;
  const tensor &LB_;
};

template <class tensor>
void Calc_projector_left_block_single(const tensor &C1, const tensor &C4,
                               const tensor &eT1, const tensor &eT6,
                               const tensor &eT7, const tensor &eT8,
                               const tensor &Tn1_s, const tensor &Tn4_s,
                               const PEPS_Parameters peps_parameters,
                               tensor &PU, tensor &PL) {

  int e12 = eT1.shape()[1];
  int e56 = eT6.shape()[0];
  int e78 = eT8.shape()[0];
  int t12 = Tn1_s.shape()[2];
  int t41 = Tn1_s.shape()[3];
  int t34 = Tn4_s.shape()[2];

  if (t41 != 1) {
    /*
        ##############################
        # (((C1*eT1)*eT8)*Tn1_s)
        # cpu_cost= 2.1e+08  memory= 2.01e+06
        # final_bond_order  (e78, t41, e12, t12)
        ##############################
    */
    tensor LT =
      reshape(tensordot(tensordot(tensordot(C1, eT1, Axes(1), Axes(0)), eT8,
                                      Axes(0), Axes(1)),
		  Tn1_s, Axes(1, 3), Axes(1, 0)).transpose(Axes(1, 3, 0, 2)),
	      Shape(e78 * t41, e12 * t12));

    /*
        ##############################
        # (((C4*eT7)*eT6)*Tn4_s)
        # cpu_cost= 2.1e+08  memory= 2.01e+06
        # final_bond_order  (e78, t41, e56, t34)
        ##############################
    */
    tensor LB =
      reshape(tensordot(tensordot(tensordot(C4, eT7, Axes(1), Axes(0)), eT6,
                                      Axes(0), Axes(1)),
			Tn4_s, Axes(1, 3), Axes(0, 3)).transpose(Axes(0, 2, 1, 3)),
	      Shape(e78 * t41, e56 * t34));

    // tensor R1 = LT;
    // tensor R2 = LB;
    if ((t12 != 1 && t34 != 1) && peps_parameters.Use_RSVD) {
      Mult_col_single<tensor> m_col(LT, LB);
      Mult_row_single<tensor> m_row(LT, LB);

      tensor U;
      tensor VT;
      std::vector<double> s;

      Shape shape_row(e12 * t12);
      Shape shape_col(e56 * t34);

      /*int info ()= */
      int cut = std::min(std::min(std::min(peps_parameters.CHI, e78 * t41), e12 * t12), e56 * t34);
      rsvd(m_row, m_col, shape_row, shape_col, U, s, VT, cut,
           static_cast<size_t>(peps_parameters.RSVD_Oversampling_factor * cut));
      double denom = s[0];
      std::vector<double> s_c;
      s_c.resize(cut);

      for (int i = 0; i < s_c.size(); ++i) {
        if (s[i] / denom > peps_parameters.Inverse_projector_cut) {
          s_c[i] = 1.0 / sqrt(s[i]);
        } else {
          s_c[i] = 0.0;
	  cut = i;
	  break;
        }
      }

      tensor U_c = mptensor::slice(U, 1, 0, cut);
      tensor VT_c = mptensor::slice(VT, 0, 0, cut);
      s_c.resize(cut);

      U_c.multiply_vector(s_c, 1);
      VT_c.multiply_vector(s_c, 0);

      PU = reshape(tensordot(LB, conj(VT_c), Axes(1), Axes(1)), Shape(e78, t41, cut));
      PL = reshape(tensordot(LT, conj(U_c), Axes(1), Axes(0)), Shape(e78, t41, cut));
      
    } else {
      // full svd //
      tensor U;
      tensor VT;
      std::vector<double> s;

      /* int info = */
      svd(tensordot(LT, LB, Axes(0), Axes(0)), Axes(0),
          Axes(1), U, s, VT);
      double denom = s[0];
      std::vector<double> s_c;
      //s_c.resize(e78);
      int cut = std::min(std::min(std::min(peps_parameters.CHI, e78 * t41), e12 * t12), e56 * t34);
      s_c.resize(cut);
      for (int i = 0; i < cut; ++i) {
        if (s[i] / denom > peps_parameters.Inverse_projector_cut) {
          s_c[i] = 1.0 / sqrt(s[i]);
        } else {
          s_c[i] = 0.0;
	  cut = i;
	  break;	    
        }
      }
      // O(D^{10})
      //tensor U_c = mptensor::slice(U, 1, 0, e78);
      //tensor VT_c = mptensor::slice(VT, 0, 0, e78);
      tensor U_c = mptensor::slice(U, 1, 0, cut);
      tensor VT_c = mptensor::slice(VT, 0, 0, cut);
      s_c.resize(cut);

      
      U_c.multiply_vector(s_c, 1);
      VT_c.multiply_vector(s_c, 0);

      
      //PU = reshape(tensordot(LB, conj(VT), Axes(1), Axes(1)), Shape(e78, t41, e78));
      //PL = reshape(tensordot(LT, conj(U), Axes(1), Axes(0)), Shape(e78, t41, e78));
      PU = reshape(tensordot(LB, conj(VT_c), Axes(1), Axes(1)), Shape(e78, t41, cut));
      PL = reshape(tensordot(LT, conj(U_c), Axes(1), Axes(0)), Shape(e78, t41, cut));

      
    }
  } else {
    const auto comm = C1.get_comm();
    tensor identity_matrix(comm, Shape(e78, e78));
    Index index;
    for (int i = 0; i < identity_matrix.local_size(); i++) {
      index = identity_matrix.global_index(i);
      if (index[0] == index[1]) {
        identity_matrix.set_value(index, 1.0);
      } else {
        identity_matrix.set_value(index, 0.0);
      }
    }
    PU = reshape(identity_matrix, Shape(e78, t41, e78));
    PL = reshape(identity_matrix, Shape(e78, t41, e78));
  }
}

template <class tensor>
void Calc_projector_updown_blocks_single(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &eT7, const tensor &eT8,
    const tensor &Tn1_s, const tensor &Tn2_s, const tensor &Tn3_s, const tensor &Tn4_s,
    const PEPS_Parameters peps_parameters, tensor &PU, tensor &PL) {
  // based on P. Corboz, T.M.Rice and M. Troyer, PRL 113, 046402(2014)

  // comment out for unused variables
  int e12 = eT1.shape()[1];
  int e34 = eT3.shape()[1];
  int e56 = eT5.shape()[1];
  int e78 = eT7.shape()[1];
  int t12 = Tn1_s.shape()[2];
  int t23 = Tn2_s.shape()[3];
  int t34 = Tn3_s.shape()[0];
  int t41 = Tn4_s.shape()[1];

  if (t41 != 1) {
    /* 
        ##############################
        # (((C1*eT1)*eT8)*Tn1_s)
        # cpu_cost= 2.1e+08  memory= 2.01e+06
        # final_bond_order  (e78, t41, e12, t12)
        ##############################
    */
    tensor LT =
      reshape(tensordot(tensordot(tensordot(C1, eT1, Axes(1), Axes(0)), eT8,
                                      Axes(0), Axes(1)),
			Tn1_s, Axes(1, 3), Axes(1, 0)).transpose(Axes(1, 3, 0, 2)),
	      Shape(e78 * t41, e12 * t12));

    /*
        ##############################
        # (((C2*eT3)*eT2)*Tn2_s)
        # cpu_cost= 2.1e+08  memory= 2.01e+06
        # final_bond_order  (e12, t12, e34, t23)
        ##############################
    */

    tensor RT =
      reshape(tensordot(tensordot(tensordot(C2, eT3, Axes(1), Axes(0)), eT2,
                                      Axes(0), Axes(1)),
		  Tn2_s, Axes(1, 3), Axes(2, 1)).transpose(Axes(1, 2, 0, 3)),
	      Shape(e12 * t12, e34 * t23));

    /* 
        ##############################
        # (((C3*eT5)*eT4)*Tn3_s)
        # cpu_cost= 2.1e+08  memory= 2.01e+06
        # final_bond_order  (e56, t34, e34, t23)
        ##############################
    */
    tensor RB =
      reshape(tensordot(tensordot(tensordot(C3, eT5, Axes(1), Axes(0)), eT4,
                                      Axes(0), Axes(1)),
		  Tn3_s, Axes(1, 3), Axes(3, 2)).transpose(Axes(0, 2, 1, 3)),
	      Shape(e56 * t34, e34 * t23));

    /*
        ##############################
        # (((C4*eT7)*eT6)*Tn4_s)
        # cpu_cost= 2.1e+08  memory= 2.01e+06
        # final_bond_order  (e78, t41, e56, t34)
        ##############################
    */
    tensor LB =
      reshape(tensordot(tensordot(tensordot(C4, eT7, Axes(1), Axes(0)), eT6,
                                      Axes(0), Axes(1)),
		  Tn4_s, Axes(1, 3), Axes(0, 3)).transpose(Axes(0, 2, 1, 3)),
	      Shape(e78 * t41, e56 * t34));
    
    if (t23 != 1 && peps_parameters.Use_RSVD) {
      Mult_col_ud_single<tensor> m_col(LT, RT, RB, LB);
      Mult_row_ud_single<tensor> m_row(LT, RT, RB, LB);

      tensor U;
      tensor VT;
      std::vector<double> s;

      Shape shape_row(e34 * t23);
      Shape shape_col(e34 * t23);

      /* int info = */
      int cut = std::min(peps_parameters.CHI, e34 * t23);
      rsvd(m_row, m_col, shape_row, shape_col, U, s, VT, cut,
           static_cast<size_t>(peps_parameters.RSVD_Oversampling_factor * cut));
      double denom = s[0];
      std::vector<double> s_c;
      s_c.resize(cut);

      for (int i = 0; i < s.size(); ++i) {
        if (s[i] / denom > peps_parameters.Inverse_projector_cut) {
          s_c[i] = 1.0 / sqrt(s[i]);
        } else {
          s_c[i] = 0.0;
	  cut = i;
	  break;
        }
      }

      tensor U_c = mptensor::slice(U, 1, 0, cut);
      tensor VT_c = mptensor::slice(VT, 0, 0, cut);
      s_c.resize(cut);

      U_c.multiply_vector(s_c, 1);
      VT_c.multiply_vector(s_c, 0);

      PU = reshape(tensordot(LB, tensordot(RB, conj(VT_c), Axes(1), Axes(1)),
			     Axes(1), Axes(0)),
		   Shape(e78, t41, cut));
      PL = reshape(tensordot(LT, tensordot(RT, conj(U_c), Axes(1), Axes(0)),
			     Axes(1), Axes(0)),
		   Shape(e78, t41, cut));
               
    } else {
      // full svd
      tensor R1 = tensordot(RT, LT, Axes(0), Axes(1));
      tensor R2 = tensordot(RB, LB, Axes(0), Axes(1));

      tensor U;
      tensor VT;
      std::vector<double> s;

      /* int info = */
      svd(tensordot(R1, R2, Axes(1), Axes(1)), Axes(0),
          Axes(0), U, s, VT);
      double denom = s[0];
      std::vector<double> s_c;
      int cut = std::min(peps_parameters.CHI, e34 * t23);
      s_c.resize(cut);

      for (int i = 0; i < e78; ++i) {
        if (s[i] / denom > peps_parameters.Inverse_projector_cut) {
          s_c[i] = 1.0 / sqrt(s[i]);
        } else {
          s_c[i] = 0.0;
	  cut = i;
	  break;
        }
      }
      //tensor U_c = mptensor::slice(U, 1, 0, e78);
      //tensor VT_c = mptensor::slice(VT, 0, 0, e78);
      tensor U_c = mptensor::slice(U, 1, 0, cut);
      tensor VT_c = mptensor::slice(VT, 0, 0, cut);
      s_c.resize(cut);

      U_c.multiply_vector(s_c, 1);
      VT_c.multiply_vector(s_c, 0);

      PU = reshape(tensordot(R2, conj(VT_c), Axes(0), Axes(1)), Shape(e78, t41, cut));
      PL = reshape(tensordot(R1, conj(U_c), Axes(0), Axes(0)), Shape(e78, t41, cut));
    }
  } else {
    const MPI_Comm comm = C1.get_comm();
    tensor identity_matrix(comm, Shape(e78, e78));
    Index index;
    for (int i = 0; i < identity_matrix.local_size(); i++) {
      index = identity_matrix.global_index(i);
      if (index[0] == index[1]) {
        identity_matrix.set_value(index, 1.0);
      } else {
        identity_matrix.set_value(index, 0.0);
      }
    }
    PU = reshape(identity_matrix, Shape(e78, t41, e78));
    PL = reshape(identity_matrix, Shape(e78, t41, e78));
  }
}

template <class tensor>
void Calc_Next_CTM_single(const tensor &C1, const tensor &C4, const tensor &eT1,
                   const tensor &eT6, const tensor &PU, const tensor &PL,
                   tensor &C1_out, tensor &C4_out) {
  C1_out = tensordot(PU, tensordot(C1, eT1, Axes(1), Axes(0)), Axes(0, 1),
                     Axes(0, 2));
  C4_out = tensordot(tensordot(eT6, C4, Axes(1), Axes(0)), PL, Axes(2, 1),
                     Axes(0, 1));

  double max_all = max_abs(C1_out);
  C1_out /= max_all;
  max_all = max_abs(C4_out);
  C4_out /= max_all;
}

template <class tensor>
void Calc_Next_eT_single(const tensor &eT8, const tensor &Tn1_s, const tensor &PU,
                  const tensor &PL, tensor &eT_out) {
  /*
    ##############################
    # ((Tn1_s*(eT8*PL))*PU)
    # cpu_cost= 3e+08  memory= 2.11e+06
    # final_bond_order  (n0, n1, k2)
    ##############################
  */

  eT_out = tensordot(tensordot(Tn1_s,
                               tensordot(eT8, PL, Axes(1), Axes(0)),
                                         Axes(0, 1), Axes(1, 2)),
                     PU, Axes(1, 2), Axes(1, 0))
               .transpose(Axes(2, 1, 0));


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
void Left_move_single(std::vector<tensor> &C1, const std::vector<tensor> &C2,
               const std::vector<tensor> &C3, std::vector<tensor> &C4,
               const std::vector<tensor> &eTt, const std::vector<tensor> &eTr,
               const std::vector<tensor> &eTb, std::vector<tensor> &eTl,
               const std::vector<tensor> &Tn_single, const int ix,
               const PEPS_Parameters peps_parameters,
               const SquareLattice lattice) {
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
      Calc_projector_left_block_single(C1[i], C4[l], eTt[i], eTb[l], eTl[l], eTl[i],
                                Tn_single[i], Tn_single[l], peps_parameters, PUs[iy],
                                PLs[iy]);
    } else {
      Calc_projector_updown_blocks_single(C1[i], C2[j], C3[k], C4[l], eTt[i], eTt[j],
                                   eTr[j], eTr[k], eTb[k], eTb[l], eTl[l],
                                   eTl[i], Tn_single[i], Tn_single[j], Tn_single[k], Tn_single[l],
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

    Calc_Next_CTM_single(C1_bak[i], C4_bak[l], eTt[i], eTb[l], PUs[iy_up],
                  PLs[iy_down], C1[j], C4[k]);
    Calc_Next_eT_single(eTl_bak[i], Tn_single[i], PUs[iy], PLs[iy_up], eTl[j]);
    Calc_Next_eT_single(eTl_bak[l], Tn_single[l], PUs[iy_down], PLs[iy], eTl[k]);
  }
}

template <class tensor>
void Right_move_single(const std::vector<tensor> &C1, std::vector<tensor> &C2,
                std::vector<tensor> &C3, const std::vector<tensor> &C4,
                const std::vector<tensor> &eTt, std::vector<tensor> &eTr,
                const std::vector<tensor> &eTb, const std::vector<tensor> &eTl,
                const std::vector<tensor> &Tn_single, const int ix,
                const PEPS_Parameters peps_parameters,
                const SquareLattice lattice) {
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
      Calc_projector_left_block_single(C3[k], C2[j], eTb[k], eTt[j], eTr[j], eTr[k],
                                transpose(Tn_single[k], Axes(2, 3, 0, 1)),
                                transpose(Tn_single[j], Axes(2, 3, 0, 1)),
                                peps_parameters, PUs[iy], PLs[iy]);
    } else {
      Calc_projector_updown_blocks_single(
          C3[k], C4[l], C1[i], C2[j], eTb[k], eTb[l], eTl[l], eTl[i], eTt[i],
          eTt[j], eTr[j], eTr[k], transpose(Tn_single[k], Axes(2, 3, 0, 1)),
          transpose(Tn_single[l], Axes(2, 3, 0, 1)),
          transpose(Tn_single[i], Axes(2, 3, 0, 1)),
          transpose(Tn_single[j], Axes(2, 3, 0, 1)), peps_parameters, PUs[iy],
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

    Calc_Next_CTM_single(C3_bak[k], C2_bak[j], eTb[k], eTt[j], PUs[iy_down],
                  PLs[iy_up], C3[l], C2[i]);

    Calc_Next_eT_single(eTr_bak[k], transpose(Tn_single[k], Axes(2, 3, 0, 1)), PUs[iy],
                 PLs[iy_down], eTr[l]);
    Calc_Next_eT_single(eTr_bak[j], transpose(Tn_single[j], Axes(2, 3, 0, 1)), PUs[iy_up],
                 PLs[iy], eTr[i]);
  }
}

template <class tensor>
void Top_move_single(std::vector<tensor> &C1, std::vector<tensor> &C2,
              const std::vector<tensor> &C3, const std::vector<tensor> &C4,
              std::vector<tensor> &eTt, const std::vector<tensor> &eTr,
              const std::vector<tensor> &eTb, const std::vector<tensor> &eTl,
              const std::vector<tensor> &Tn_single, const int iy,
              const PEPS_Parameters peps_parameters,
              const SquareLattice lattice) {
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
      Calc_projector_left_block_single(C2[j], C1[i], eTr[j], eTl[i], eTt[i], eTt[j],
                                transpose(Tn_single[j], Axes(1, 2, 3, 0)),
                                transpose(Tn_single[i], Axes(1, 2, 3, 0)),
                                peps_parameters, PUs[ix], PLs[ix]);
    } else {
      Calc_projector_updown_blocks_single(
          C2[j], C3[k], C4[l], C1[i], eTr[j], eTr[k], eTb[k], eTb[l], eTl[l],
          eTl[i], eTt[i], eTt[j], transpose(Tn_single[j], Axes(1, 2, 3, 0)),
          transpose(Tn_single[k], Axes(1, 2, 3, 0)),
          transpose(Tn_single[l], Axes(1, 2, 3, 0)),
          transpose(Tn_single[i], Axes(1, 2, 3, 0)), peps_parameters, PUs[ix],
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

    Calc_Next_CTM_single(C2_bak[j], C1_bak[i], eTr[j], eTl[i], PUs[ix_right],
                  PLs[ix_left], C2[k], C1[l]);

    Calc_Next_eT_single(eTt_bak[j], transpose(Tn_single[j], Axes(1, 2, 3, 0)), PUs[ix],
                 PLs[ix_right], eTt[k]);
    Calc_Next_eT_single(eTt_bak[i], transpose(Tn_single[i], Axes(1, 2, 3, 0)),
                 PUs[ix_left], PLs[ix], eTt[l]);
  }
}

template <class tensor>
void Bottom_move_single(const std::vector<tensor> &C1, const std::vector<tensor> &C2,
                 std::vector<tensor> &C3, std::vector<tensor> &C4,
                 const std::vector<tensor> &eTt, const std::vector<tensor> &eTr,
                 std::vector<tensor> &eTb, const std::vector<tensor> &eTl,
                 const std::vector<tensor> &Tn_single, const int iy,
                 const PEPS_Parameters peps_parameters,
                 const SquareLattice lattice) {
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
      Calc_projector_left_block_single(C4[l], C3[k], eTl[l], eTr[k], eTb[k], eTb[l],
                                transpose(Tn_single[l], Axes(3, 0, 1, 2)),
                                transpose(Tn_single[k], Axes(3, 0, 1, 2)),
                                peps_parameters, PUs[ix], PLs[ix]);
    } else {
      Calc_projector_updown_blocks_single(
          C4[l], C1[i], C2[j], C3[k], eTl[l], eTl[i], eTt[i], eTt[j], eTr[j],
          eTr[k], eTb[k], eTb[l], transpose(Tn_single[l], Axes(3, 0, 1, 2)),
          transpose(Tn_single[i], Axes(3, 0, 1, 2)),
          transpose(Tn_single[j], Axes(3, 0, 1, 2)),
          transpose(Tn_single[k], Axes(3, 0, 1, 2)), peps_parameters, PUs[ix],
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

    Calc_Next_CTM_single(C4_bak[l], C3_bak[k], eTl[l], eTr[k], PUs[ix_left],
                  PLs[ix_right], C4[i], C3[j]);

    Calc_Next_eT_single(eTb_bak[l], transpose(Tn_single[l], Axes(3, 0, 1, 2)), PUs[ix],
                 PLs[ix_left], eTb[i]);
    Calc_Next_eT_single(eTb_bak[k], transpose(Tn_single[k], Axes(3, 0, 1, 2)),
                 PUs[ix_right], PLs[ix], eTb[j]);
  }
}

template <class tensor>
int Calc_CTM_Environment_density(std::vector<tensor> &C1, std::vector<tensor> &C2,
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

  std::vector<tensor> Tn_single = Make_single_tensor_density(Tn);
  if (initialize) {
    int num, d1, d2, d34;
    Index index;
    tensor Projector;

    for (int i = 0; i < lattice.N_UNIT; ++i) {
      num = lattice.top(lattice.left(i));
      d1 = Tn[num].shape()[3];
      d2 = Tn[num].shape()[2];
      C1[i] = reshape(mptensor::slice(mptensor::slice(Tn_single[num], 0, 0, 1), 1, 0, 1), Shape(d2, d1)).transpose(Axes(1, 0));
      /*
      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
	C1[i] =	extend(C1[i], Shape(peps_parameters.CHI, peps_parameters.CHI));

      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        C1[i] = mptensor::slice(mptensor::slice(C1[i], 0, 0, peps_parameters.CHI), 1, 0,
                      peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        C1[i] = extend(mptensor::slice(C1[i], 1, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else {
        // d1 >= CHI
        // d2 < CHI
	
        C1[i] = extend(mptensor::slice(C1[i], 0, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      }
      */

      num = lattice.top(lattice.right(i));
      d1 = Tn[num].shape()[0];
      d2 = Tn[num].shape()[3];
      C2[i] = reshape(mptensor::slice(mptensor::slice(Tn_single[num], 1, 0, 1), 2, 0, 1), Shape(d1,d2));
      /*
      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        C2[i] = extend(C2[i], Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        C2[i] = mptensor::slice(mptensor::slice(C2[i], 0, 0, peps_parameters.CHI), 1, 0,
                      peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        C2[i] = extend(mptensor::slice(C2[i], 1, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else {
        // d1 >= CHI
        // d2 < CHI
        C2[i] = extend(mptensor::slice(C2[i], 0, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      }
      */
      num = lattice.bottom(lattice.right(i));
      d1 = Tn[num].shape()[1];
      d2 = Tn[num].shape()[0];
      C3[i] = reshape(mptensor::slice(mptensor::slice(Tn_single[num], 2, 0, 1), 3, 0, 1), Shape(d2, d1)).transpose(Axes(1, 0));

      /*
      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        C3[i] = extend(C3[i], Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        C3[i] = mptensor::slice(mptensor::slice(C3[i], 0, 0, peps_parameters.CHI), 1, 0,
                      peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        C3[i] = extend(mptensor::slice(C3[i], 1, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else {
        // d1 >= CHI
        // d2 < CHI
        C3[i] = extend(mptensor::slice(C3[i], 0, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      }
      */
      
      num = lattice.bottom(lattice.left(i));
      d1 = Tn[num].shape()[2];
      d2 = Tn[num].shape()[1];
      C4[i] = reshape(mptensor::slice(mptensor::slice(Tn_single[num], 0, 0, 1), 3, 0, 1), Shape(d2, d1)).transpose(Axes(1, 0));

      /*
      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        C4[i] = extend(C4[i], Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        C4[i] = mptensor::slice(mptensor::slice(C4[i], 0, 0, peps_parameters.CHI), 1, 0,
                      peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        C4[i] = extend(mptensor::slice(C4[i], 1, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else {
        // d1 >= CHI
        // d2 < CHI
        C4[i] = extend(mptensor::slice(C4[i], 0, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      }
      */

      num = lattice.top(i);
      d1 = Tn[num].shape()[0];
      d2 = Tn[num].shape()[2];
      d34 = Tn[num].shape()[3];
      eTt[i] = reshape(mptensor::slice(Tn_single[num], 1, 0, 1), Shape(d1, d2, d34));

      /*
      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        eTt[i] = extend(
            eTt[i], Shape(peps_parameters.CHI, peps_parameters.CHI, d34));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        eTt[i] = mptensor::slice(mptensor::slice(eTt[i], 0, 0, peps_parameters.CHI), 1, 0,
                       peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        eTt[i] =
            extend(mptensor::slice(eTt[i], 1, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34));
      } else {
        // d1 >= CHI
        // d2 < CHI
        eTt[i] =
            extend(mptensor::slice(eTt[i], 0, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34));
      }
      */
      
      num = lattice.right(i);
      d1 = Tn[num].shape()[1];
      d2 = Tn[num].shape()[3];
      d34 = Tn[num].shape()[0];
      eTr[i] = reshape(mptensor::slice(Tn_single[num], 2, 0, 1), Shape(d34, d1, d2)).transpose(Axes(1, 2, 0));

      /*
      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        eTr[i] = extend(
            eTr[i], Shape(peps_parameters.CHI, peps_parameters.CHI, d34));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        eTr[i] = mptensor::slice(mptensor::slice(eTr[i], 0, 0, peps_parameters.CHI), 1, 0,
                       peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        eTr[i] =
            extend(mptensor::slice(eTr[i], 1, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34));
      } else {
        // d1 >= CHI
        // d2 < CHI
        eTr[i] =
            extend(mptensor::slice(eTr[i], 0, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34));
      }
      */
      
      num = lattice.bottom(i);
      d1 = Tn[num].shape()[2];
      d2 = Tn[num].shape()[0];
      d34 = Tn[num].shape()[1];
      eTb[i] = reshape(mptensor::slice(Tn_single[num], 3, 0, 1), Shape(d2, d34, d1)).transpose(Axes(2, 0, 1));

      /*
      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        eTb[i] = extend(
            eTb[i], Shape(peps_parameters.CHI, peps_parameters.CHI, d34));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        eTb[i] = mptensor::slice(mptensor::slice(eTb[i], 0, 0, peps_parameters.CHI), 1, 0,
                       peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        eTb[i] =
            extend(mptensor::slice(eTb[i], 1, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34));
      } else {
        // d1 >= CHI
        // d2 < CHI
        eTb[i] =
            extend(mptensor::slice(eTb[i], 0, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34));
      }
      */
      
      num = lattice.left(i);
      d1 = Tn[num].shape()[3];
      d2 = Tn[num].shape()[1];
      d34 = Tn[num].shape()[2];
      eTl[i] = reshape(mptensor::slice(Tn_single[num], 0, 0, 1), Shape(d2, d34, d1)).transpose(Axes(2, 0, 1));

      /*
      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        eTl[i] = extend(
            eTl[i], Shape(peps_parameters.CHI, peps_parameters.CHI, d34));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        eTl[i] = mptensor::slice(mptensor::slice(eTl[i], 0, 0, peps_parameters.CHI), 1, 0,
                       peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        eTl[i] =
            extend(mptensor::slice(eTl[i], 1, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34));
      } else {
        // d1 >= CHI
        // d2 < CHI
        eTl[i] =
            extend(mptensor::slice(eTl[i], 0, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34));
      }
      */
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
      Left_move_single(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn_single, ix, peps_parameters,
                lattice);
    }

    // right move
    for (int ix = 0; ix > -lattice.LX; --ix) {
      Right_move_single(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn_single,
                 (ix + 1 + lattice.LX) % lattice.LX, peps_parameters, lattice);
    }

    // top move
    for (int iy = 0; iy > -lattice.LY; --iy) {
      Top_move_single(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn_single,
               (iy + 1 + lattice.LY) % lattice.LY, peps_parameters, lattice);
    }

    // bottom move

    for (int iy = 0; iy < lattice.LY; ++iy) {
      Bottom_move_single(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn_single, iy, peps_parameters,
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

template void Calc_projector_left_block_single(
    const real_tensor &C1, const real_tensor &C4, const real_tensor &eT1,
    const real_tensor &eT6, const real_tensor &eT7, const real_tensor &eT8,
    const real_tensor &Tn1, const real_tensor &Tn4,
    const PEPS_Parameters peps_parameters, real_tensor &PU, real_tensor &PL);
template void Calc_projector_left_block_single(
    const complex_tensor &C1, const complex_tensor &C4,
    const complex_tensor &eT1, const complex_tensor &eT6,
    const complex_tensor &eT7, const complex_tensor &eT8,
    const complex_tensor &Tn1, const complex_tensor &Tn4,
    const PEPS_Parameters peps_parameters, complex_tensor &PU,
    complex_tensor &PL);

template void Calc_projector_updown_blocks_single(
    const real_tensor &C1, const real_tensor &C2, const real_tensor &C3,
    const real_tensor &C4, const real_tensor &eT1, const real_tensor &eT2,
    const real_tensor &eT3, const real_tensor &eT4, const real_tensor &eT5,
    const real_tensor &eT6, const real_tensor &eT7, const real_tensor &eT8,
    const real_tensor &Tn1, const real_tensor &Tn2, const real_tensor &Tn3,
    const real_tensor &Tn4, const PEPS_Parameters peps_parameters,
    real_tensor &PU, real_tensor &PL);
template void Calc_projector_updown_blocks_single(
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

template void Calc_Next_CTM_single(const real_tensor &C1, const real_tensor &C4,
                            const real_tensor &eT1, const real_tensor &eT6,
                            const real_tensor &PU, const real_tensor &PL,
                            real_tensor &C1_out, real_tensor &C4_out);
template void Calc_Next_CTM_single(const complex_tensor &C1, const complex_tensor &C4,
                            const complex_tensor &eT1,
                            const complex_tensor &eT6, const complex_tensor &PU,
                            const complex_tensor &PL, complex_tensor &C1_out,
                            complex_tensor &C4_out);

template void Calc_Next_eT_single(const real_tensor &eT8, const real_tensor &Tn1,
                           const real_tensor &PU, const real_tensor &PL,
                           real_tensor &eT_out);
template void Calc_Next_eT_single(const complex_tensor &eT8, const complex_tensor &Tn1,
                           const complex_tensor &PU, const complex_tensor &PL,
                           complex_tensor &eT_out);

template void Left_move_single(
    std::vector<real_tensor> &C1, const std::vector<real_tensor> &C2,
    const std::vector<real_tensor> &C3, std::vector<real_tensor> &C4,
    const std::vector<real_tensor> &eTt, const std::vector<real_tensor> &eTr,
    const std::vector<real_tensor> &eTb, std::vector<real_tensor> &eTl,
    const std::vector<real_tensor> &Tn, const int ix,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice);
template void Left_move_single(
    std::vector<complex_tensor> &C1, const std::vector<complex_tensor> &C2,
    const std::vector<complex_tensor> &C3, std::vector<complex_tensor> &C4,
    const std::vector<complex_tensor> &eTt,
    const std::vector<complex_tensor> &eTr,
    const std::vector<complex_tensor> &eTb, std::vector<complex_tensor> &eTl,
    const std::vector<complex_tensor> &Tn, const int ix,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice);

template void Right_move_single(
    const std::vector<real_tensor> &C1, std::vector<real_tensor> &C2,
    std::vector<real_tensor> &C3, const std::vector<real_tensor> &C4,
    const std::vector<real_tensor> &eTt, std::vector<real_tensor> &eTr,
    const std::vector<real_tensor> &eTb, const std::vector<real_tensor> &eTl,
    const std::vector<real_tensor> &Tn, const int ix,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice);
template void Right_move_single(
    const std::vector<complex_tensor> &C1, std::vector<complex_tensor> &C2,
    std::vector<complex_tensor> &C3, const std::vector<complex_tensor> &C4,
    const std::vector<complex_tensor> &eTt, std::vector<complex_tensor> &eTr,
    const std::vector<complex_tensor> &eTb,
    const std::vector<complex_tensor> &eTl,
    const std::vector<complex_tensor> &Tn, const int ix,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice);

template void Top_move_single(
    std::vector<real_tensor> &C1, std::vector<real_tensor> &C2,
    const std::vector<real_tensor> &C3, const std::vector<real_tensor> &C4,
    std::vector<real_tensor> &eTt, const std::vector<real_tensor> &eTr,
    const std::vector<real_tensor> &eTb, const std::vector<real_tensor> &eTl,
    const std::vector<real_tensor> &Tn, const int iy,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice);
template void Top_move_single(std::vector<complex_tensor> &C1,
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

template void Bottom_move_single(
    const std::vector<real_tensor> &C1, const std::vector<real_tensor> &C2,
    std::vector<real_tensor> &C3, std::vector<real_tensor> &C4,
    const std::vector<real_tensor> &eTt, const std::vector<real_tensor> &eTr,
    std::vector<real_tensor> &eTb, const std::vector<real_tensor> &eTl,
    const std::vector<real_tensor> &Tn, const int iy,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice);
template void Bottom_move_single(
    const std::vector<complex_tensor> &C1,
    const std::vector<complex_tensor> &C2, std::vector<complex_tensor> &C3,
    std::vector<complex_tensor> &C4, const std::vector<complex_tensor> &eTt,
    const std::vector<complex_tensor> &eTr, std::vector<complex_tensor> &eTb,
    const std::vector<complex_tensor> &eTl,
    const std::vector<complex_tensor> &Tn, const int iy,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice);


template int Calc_CTM_Environment_density(
    std::vector<real_tensor> &C1, std::vector<real_tensor> &C2,
    std::vector<real_tensor> &C3, std::vector<real_tensor> &C4,
    std::vector<real_tensor> &eTt, std::vector<real_tensor> &eTr,
    std::vector<real_tensor> &eTb, std::vector<real_tensor> &eTl,
    const std::vector<real_tensor> &Tn, const PEPS_Parameters peps_parameters,
    const SquareLattice lattice, bool initialize);

template int Calc_CTM_Environment_density(
    std::vector<complex_tensor> &C1, std::vector<complex_tensor> &C2,
    std::vector<complex_tensor> &C3, std::vector<complex_tensor> &C4,
    std::vector<complex_tensor> &eTt, std::vector<complex_tensor> &eTr,
    std::vector<complex_tensor> &eTb, std::vector<complex_tensor> &eTl,
    const std::vector<complex_tensor> &Tn,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice,
    bool initialize);

}  // end of namespace core
}  // namespace itps
}  // namespace tenes
