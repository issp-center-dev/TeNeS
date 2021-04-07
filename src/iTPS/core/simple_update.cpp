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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <mptensor/complex.hpp>
#include <mptensor/rsvd.hpp>
#include <mptensor/tensor.hpp>

#include "../../tensor.hpp"
#include "../PEPS_Parameters.hpp"

namespace tenes {

using mptensor::Axes;
using mptensor::Shape;

// environment

template <class tensor>
void Simple_update_bond(const tensor &Tn1, const tensor &Tn2,
                        const std::vector<std::vector<double>> &lambda1,
                        const std::vector<std::vector<double>> &lambda2,
                        const tensor &op12, const int connect1,
                        const PEPS_Parameters peps_parameters, tensor &Tn1_new,
                        tensor &Tn2_new, std::vector<double> &lambda_c) {
  using ptensor = tensor;
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
  }
  for (int i = 0; i < 4; i++) {
    lambda2_inv[i] = std::vector<double>(lambda2[i].size());
    for (int j = 0; j < lambda2_inv[i].size(); ++j) {
      if (lambda2[i][j] > peps_parameters.Inverse_lambda_cut) {
        lambda2_inv[i][j] = 1.0 / lambda2[i][j];
      } else {
        lambda2_inv[i][j] = 0.0;
      }
    }
  }

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
  }

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
  }
  norm = sqrt(norm);
  for (int i = 0; i < dc; ++i) {
    lambda_c[i] = sqrt(lambda_c[i] / norm);
  }

  /*for (int i=0; i < VTc.local_size();++i){
    Index index = VTc.global_index(i);
    std::cout<<"VTC[i,j]="<<index<<", "<<VTc[i]<<std::endl;
    }*/

  /*
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
  */

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
  }

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
  }
}

// template instantiations

template void Simple_update_bond(
    const real_tensor &Tn1, const real_tensor &Tn2,
    const std::vector<std::vector<double>> &lambda1,
    const std::vector<std::vector<double>> &lambda2, const real_tensor &op12,
    const int connect1, const PEPS_Parameters peps_parameters,
    real_tensor &Tn1_new, real_tensor &Tn2_new, std::vector<double> &lambda_c);

template void Simple_update_bond(
    const complex_tensor &Tn1, const complex_tensor &Tn2,
    const std::vector<std::vector<double>> &lambda1,
    const std::vector<std::vector<double>> &lambda2, const complex_tensor &op12,
    const int connect1, const PEPS_Parameters peps_parameters,
    complex_tensor &Tn1_new, complex_tensor &Tn2_new,
    std::vector<double> &lambda_c);

}  // end of namespace tenes
