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

#include "full_update.hpp"

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

#include "../../mpi.hpp"
#include "../../tensor.hpp"
#include "../PEPS_Parameters.hpp"

namespace tenes {

using mptensor::Axes;
using mptensor::Shape;

template <class tensor>
tensor Create_Environment_two_sites(const tensor &C1, const tensor &C2,
                                    const tensor &C3, const tensor &C4,
                                    const tensor &eT1, const tensor &eT2,
                                    const tensor &eT3, const tensor &eT4,
                                    const tensor &eT5, const tensor &eT6,
                                    const tensor &Q1, const tensor &Q2) {
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

template <class tensor>
void Full_update_bond_horizontal(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op12, const PEPS_Parameters peps_parameters, tensor &Tn1_new,
    tensor &Tn2_new) {
  Shape Tn1_shape = Tn1.shape();
  Shape Tn2_shape = Tn2.shape();

  int D_connect = Tn1_shape[2];

  // Connecting [2] bond of Tn1 and [0] bond of Tn2
  // QR decomposition
  tensor Q1, R1, Q2, R2;

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

  tensor Theta = tensordot(tensordot(R1, R2, Axes(1), Axes(1)), op12,
                           Axes(1, 3), Axes(0, 1));
  // Environment
  // bond order (t1, t2, tc1, tc2)

  tensor Environment =
      Create_Environment_two_sites(C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6,
                                   Q1, transpose(Q2, Axes(3, 0, 1, 2)));

  // Hermite
  Environment =
      0.5 * (Environment + conj(transpose(Environment, Axes(2, 3, 0, 1))));

  // diagonalization
  tensor Z;
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
  }

  Z.multiply_vector(w, 2);

  tensor LR1, LR2, LR1_inv, LR2_inv;
  if (peps_parameters.Full_Gauge_Fix) {
    // gauge fix
    tensor Q_temp;
    info = qr(Z, Axes(2, 1), Axes(0), Q_temp, LR1);
    info = qr(Z, Axes(2, 0), Axes(1), Q_temp, LR2);
    // for theta
    Theta = tensordot(LR1, tensordot(LR2, Theta, Axes(1), Axes(1)), Axes(1),
                      Axes(1));

    // for environment
    tensor u, vt;
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

  tensor U, VT;
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
  }

  norm = sqrt(norm);

  for (int i = 0; i < D_connect; ++i) {
    lambda_c[i] = sqrt(lambda_c[i] / norm);
  }

  U.multiply_vector(lambda_c, 2);
  VT.multiply_vector(lambda_c, 0);
  R1 = transpose(U, Axes(0, 2, 1));   // envR1 , D_connect, m1
  R2 = transpose(VT, Axes(1, 0, 2));  // envR2 , D_connect, m2

  int count = 0;
  typename tensor::value_type delta = 0;
  auto C_phi = trace(tensordot(Environment, Theta, Axes(0, 1), Axes(0, 1)),
                     conj(Theta), Axes(0, 1, 2, 3), Axes(0, 1, 2, 3));
  auto Old_delta =
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

  tensor W_vec, N_mat, N_mat_inv;
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

    // Moore-Penrose Pseudo Inverse (for Hermitian matrix)
    info = svd(N_mat, Axes(1, 3), Axes(0, 2), U, s, VT);
    denom = s[0];
    for (int i = 0; i < envR1 * D_connect; ++i) {
      if (s[i] / denom > peps_parameters.Full_Inverse_precision) {
        s[i] = 1.0 / s[i];
      } else {
        s[i] = 0.0;
      }
    }

    U.multiply_vector(s, 2);
    N_mat_inv = tensordot(conj(U), conj(VT), Axes(2), Axes(0));
    R1 = tensordot(N_mat_inv, W_vec, Axes(0, 1), Axes(0, 2));

    /*
      ## for R2
      ## create W
      ## ((envR2,m2, D_connect)*)
    */
    W_vec = tensordot(tensordot(Environment, Theta, Axes(0, 1), Axes(0, 1)),
                      conj(R1), Axes(0, 2), Axes(0, 2));

    /*
      ## create N
      ## (envR2, envR2*,D_connect, D_connect*)
    */
    N_mat = tensordot(tensordot(Environment, R1, Axes(0), Axes(0)), conj(R1),
                      Axes(1, 4), Axes(0, 2));

    // Moore-Penrose Pseudo Inverse (for Hermitian matrix)
    info = svd(N_mat, Axes(1, 3), Axes(0, 2), U, s, VT);

    denom = s[0];
    for (int i = 0; i < envR2 * D_connect; ++i) {
      if (s[i] / denom > peps_parameters.Full_Inverse_precision) {
        s[i] = 1.0 / s[i];
      } else {
        s[i] = 0.0;
      }
    }

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
    }
    Old_delta = delta;
    count += 1;
  }
  // Post processing
  if (!convergence && peps_parameters.print_level >= PrintLevel::warn) {
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

  tensor q1, r1, q2, r2;

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

template <class tensor>
void Full_update_bond(const tensor &C1, const tensor &C2, const tensor &C3,
                      const tensor &C4, const tensor &eT1, const tensor &eT2,
                      const tensor &eT3, const tensor &eT4, const tensor &eT5,
                      const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
                      const tensor &op12, const int connect1,
                      const PEPS_Parameters peps_parameters, tensor &Tn1_new,
                      tensor &Tn2_new) {
  tensor Tn1_rot, Tn2_rot;
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

// template instantiations

template void Full_update_bond_horizontal(
    const real_tensor &C1, const real_tensor &C2, const real_tensor &C3,
    const real_tensor &C4, const real_tensor &eT1, const real_tensor &eT2,
    const real_tensor &eT3, const real_tensor &eT4, const real_tensor &eT5,
    const real_tensor &eT6, const real_tensor &Tn1, const real_tensor &Tn2,
    const real_tensor &op12, const PEPS_Parameters peps_parameters,
    real_tensor &Tn1_new, real_tensor &Tn2_new);

template void Full_update_bond_horizontal(
    const complex_tensor &C1, const complex_tensor &C2,
    const complex_tensor &C3, const complex_tensor &C4,
    const complex_tensor &eT1, const complex_tensor &eT2,
    const complex_tensor &eT3, const complex_tensor &eT4,
    const complex_tensor &eT5, const complex_tensor &eT6,
    const complex_tensor &Tn1, const complex_tensor &Tn2,
    const complex_tensor &op12, const PEPS_Parameters peps_parameters,
    complex_tensor &Tn1_new, complex_tensor &Tn2_new);

template void Full_update_bond(const real_tensor &C1, const real_tensor &C2,
                               const real_tensor &C3, const real_tensor &C4,
                               const real_tensor &eT1, const real_tensor &eT2,
                               const real_tensor &eT3, const real_tensor &eT4,
                               const real_tensor &eT5, const real_tensor &eT6,
                               const real_tensor &Tn1, const real_tensor &Tn2,
                               const real_tensor &op12, const int connect1,
                               const PEPS_Parameters peps_parameters,
                               real_tensor &Tn1_new, real_tensor &Tn2_new);

template void Full_update_bond(
    const complex_tensor &C1, const complex_tensor &C2,
    const complex_tensor &C3, const complex_tensor &C4,
    const complex_tensor &eT1, const complex_tensor &eT2,
    const complex_tensor &eT3, const complex_tensor &eT4,
    const complex_tensor &eT5, const complex_tensor &eT6,
    const complex_tensor &Tn1, const complex_tensor &Tn2,
    const complex_tensor &op12, const int connect1,
    const PEPS_Parameters peps_parameters, complex_tensor &Tn1_new,
    complex_tensor &Tn2_new);

}  // end of namespace tenes
