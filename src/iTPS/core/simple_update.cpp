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

#include "../../fermion/fops.hpp"
#include "../../fermion/ftensor.hpp"
#include "../../tensor.hpp"
#include "../PEPS_Parameters.hpp"

namespace tenes::itps::core {

using mptensor::Axes;
using mptensor::Shape;

namespace {
// The plain-tensor overloads of enforce_even_parity / log_theta_blocks are
// no-ops so the shared kernel below compiles for both the bosonic path
// (plain tensors) and the fermionic one (ftensor).
template <class tensor>
void enforce_even_parity(tensor &) {}

// Diagnostic knob: TENES_FERMION_THETA_LOG_LIMIT caps how many theta
// tensors have their per-parity-sector norms dumped to stderr by
// log_theta_blocks (default 0: none).
int fermion_theta_log_limit() {
  const char *raw = std::getenv("TENES_FERMION_THETA_LOG_LIMIT");
  if (raw == nullptr) {
    return 0;
  }
  return std::max(0, std::atoi(raw));
}

template <class tensor>
void log_theta_blocks(const tensor &, const char *, const Axes &,
                      const Axes &) {}

template <class tensor>
void log_theta_blocks(const tenes::fermion::ftensor<tensor> &theta,
                      const char *label, const Axes &rows, const Axes &cols) {
  static int theta_log_count = 0;
  const int limit = fermion_theta_log_limit();
  if (theta_log_count >= limit) {
    return;
  }
  double even_norm2 = 0.0;
  double odd_norm2 = 0.0;
  double even_max = 0.0;
  double odd_max = 0.0;
  mptensor::Index idx;
  idx.resize(theta.t.shape().size());
  for (std::size_t n = 0; n < theta.t.local_size(); ++n) {
    theta.t.global_index_fast(n, idx);
    bool row_odd = false;
    for (std::size_t i = 0; i < rows.size(); ++i) {
      row_odd = row_odd != theta.parity[rows[i]][idx[rows[i]]];
    }
    bool col_odd = false;
    for (std::size_t i = 0; i < cols.size(); ++i) {
      col_odd = col_odd != theta.parity[cols[i]][idx[cols[i]]];
    }
    if (row_odd != col_odd) {
      continue;
    }
    const double abs_value = std::abs(theta.t[n]);
    if (row_odd) {
      odd_norm2 += abs_value * abs_value;
      odd_max = std::max(odd_max, abs_value);
    } else {
      even_norm2 += abs_value * abs_value;
      even_max = std::max(even_max, abs_value);
    }
  }
  std::cerr << "TENES_FERMION_THETA call=" << theta_log_count
            << " label=" << label << " even_norm=" << std::sqrt(even_norm2)
            << " even_max=" << even_max << " odd_norm=" << std::sqrt(odd_norm2)
            << " odd_max=" << odd_max << "\n";
  ++theta_log_count;
}

// Guard: the graded update must keep the state parity even; elements in the
// odd sector above the tolerance indicate a sign-bookkeeping bug upstream.
template <class tensor>
void enforce_even_parity(tenes::fermion::ftensor<tensor> &a) {
  const double v = tenes::fermion::parity_violation(a);
  const double scale = std::max(1.0, tenes::fermion::max_abs(a));
  const double threshold = 1.0e-10 * scale;
  if (v > threshold) {
    std::stringstream ss;
    ss << "fermion Simple_update_bond produced odd-parity elements: max_abs="
       << v << " threshold=" << threshold;
    throw std::runtime_error(ss.str());
  }
  mptensor::Index index;
  index.resize(a.t.shape().size());
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    a.t.global_index_fast(n, index);
    if (tenes::fermion::count_odd(a.parity, index) % 2 == 1) {
      a.t[n] = typename tensor::value_type{};
    }
  }
}
}  // namespace

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
    for (size_t j = 0; j < lambda1_inv[i].size(); ++j) {
      if (lambda1[i][j] > peps_parameters.Inverse_lambda_cut) {
        lambda1_inv[i][j] = 1.0 / lambda1[i][j];
      } else {
        lambda1_inv[i][j] = 0.0;
      }
    }
  }
  for (int i = 0; i < 4; i++) {
    lambda2_inv[i] = std::vector<double>(lambda2[i].size());
    for (size_t j = 0; j < lambda2_inv[i].size(); ++j) {
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
  ptensor Theta_before = tensordot(R1, R2, Axes(1), Axes(1));
  log_theta_blocks(Theta_before, "before_gate", Axes(0, 1), Axes(2, 3));
  // Theta legs: (aux1, aux2, out1, out2); the bipartition is site 1 =
  // (aux1, out1) against site 2 = (aux2, out2). svd_trunc performs the
  // graded transpose implied by these axes itself, so no regrouping here.
  ptensor Theta = tensordot(Theta_before, op12, Axes(1, 3), Axes(0, 1));
  log_theta_blocks(Theta, "after_gate", Axes(0, 2), Axes(1, 3));

  // svd
  ptensor U, VT;
  std::vector<double> s;
  info = svd_trunc(Theta, Axes(0, 2), Axes(1, 3), U, s, VT, dc);

  lambda_c = std::vector<double>(s.begin(), s.end());
  ptensor Uc = U;
  ptensor VTc = VT;

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

  enforce_even_parity(Tn1_new);
  enforce_even_parity(Tn2_new);
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

template void Simple_update_bond(
    const tenes::fermion::ftensor<real_tensor> &Tn1,
    const tenes::fermion::ftensor<real_tensor> &Tn2,
    const std::vector<std::vector<double>> &lambda1,
    const std::vector<std::vector<double>> &lambda2,
    const tenes::fermion::ftensor<real_tensor> &op12, const int connect1,
    const PEPS_Parameters peps_parameters,
    tenes::fermion::ftensor<real_tensor> &Tn1_new,
    tenes::fermion::ftensor<real_tensor> &Tn2_new,
    std::vector<double> &lambda_c);

template void Simple_update_bond(
    const tenes::fermion::ftensor<complex_tensor> &Tn1,
    const tenes::fermion::ftensor<complex_tensor> &Tn2,
    const std::vector<std::vector<double>> &lambda1,
    const std::vector<std::vector<double>> &lambda2,
    const tenes::fermion::ftensor<complex_tensor> &op12, const int connect1,
    const PEPS_Parameters peps_parameters,
    tenes::fermion::ftensor<complex_tensor> &Tn1_new,
    tenes::fermion::ftensor<complex_tensor> &Tn2_new,
    std::vector<double> &lambda_c);

}  // namespace tenes::itps::core
