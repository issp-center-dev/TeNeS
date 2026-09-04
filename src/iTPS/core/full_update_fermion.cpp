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

#include "full_update_fermion.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../../fermion/fops.hpp"
#include "../../fermion/full_update_env.hpp"
#include "../../tensor.hpp"
#include "../PEPS_Parameters.hpp"
#include "full_update.hpp"

namespace tenes::itps::core {

using mptensor::Axes;

namespace {

bool fermion_sector_log_enabled() {
  return std::getenv("TENES_FERMION_SECTOR_LOG") != nullptr;
}

bool fermion_full_update_log_enabled() {
  return std::getenv("TENES_FERMION_FULL_UPDATE_LOG") != nullptr;
}

const char* direction_name(tenes::fermion::reduced_pair_direction direction) {
  return direction == tenes::fermion::reduced_pair_direction::horizontal
             ? "horizontal"
             : "vertical";
}

int connect_leg_a(tenes::fermion::reduced_pair_direction direction) {
  return direction == tenes::fermion::reduced_pair_direction::horizontal ? 2 : 3;
}

int connect_leg_b(tenes::fermion::reduced_pair_direction direction) {
  return direction == tenes::fermion::reduced_pair_direction::horizontal ? 0 : 1;
}

template <class tensor>
void apply_pair_mask(tensor& a, const tenes::fermion::parity_vector& p1,
                     const tenes::fermion::parity_vector& p2, int ax1,
                     int ax2) {
  mptensor::Index idx;
  idx.resize(a.shape().size());
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    a.global_index_fast(n, idx);
    if (p1[idx[ax1]] && p2[idx[ax2]]) {
      a[n] = -a[n];
    }
  }
}

template <class tensor>
void project_even(tenes::fermion::ftensor<tensor>& a) {
  mptensor::Index idx;
  idx.resize(a.t.shape().size());
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    a.t.global_index_fast(n, idx);
    if (tenes::fermion::count_odd(a.parity, idx) % 2 == 1) {
      a.t[n] = typename tensor::value_type{};
    }
  }
}

template <class tensor>
double parity_ratio(const tenes::fermion::ftensor<tensor>& a,
                    double& forbidden_abs, double& scale) {
  forbidden_abs = tenes::fermion::parity_violation(a);
  scale = std::max(1.0, tenes::fermion::max_abs(a));
  return forbidden_abs / scale;
}

template <class tensor>
double project_or_throw(tenes::fermion::ftensor<tensor>& a,
                        const std::string& label, bool throw_on_large = true) {
  double forbidden_abs = 0.0;
  double scale = 1.0;
  const double ratio = parity_ratio(a, forbidden_abs, scale);
  if (fermion_full_update_log_enabled() && a.t.get_comm_rank() == 0) {
    std::cerr << "TENES_FERMION_FULL_UPDATE_PARITY label=" << label
              << " ratio=" << ratio << " forbidden_max=" << forbidden_abs
              << " max_abs=" << scale
              << " projected=" << (forbidden_abs > 0.0 ? 1 : 0) << "\n";
  }
  if (throw_on_large && ratio > 1.0e-8) {
    std::stringstream ss;
    ss << "fermion Full_update_bond parity projection rejected " << label
       << ": forbidden parity block ratio " << ratio << " exceeds 1e-8";
    throw std::runtime_error(ss.str());
  }
  project_even(a);
  return ratio;
}

std::pair<int, int> sector_dimensions(
    const tenes::fermion::parity_vector& parity) {
  int even = 0;
  int odd = 0;
  for (bool p : parity) {
    if (p) {
      ++odd;
    } else {
      ++even;
    }
  }
  return std::make_pair(even, odd);
}

template <class tensor>
void log_sector_dimensions(
    const tenes::fermion::ftensor<tensor>& U,
    tenes::fermion::reduced_pair_direction direction,
    const PEPS_Parameters& peps_parameters) {
  const auto dims = sector_dimensions(U.parity[2]);
  if (fermion_sector_log_enabled() && U.t.get_comm_rank() == 0) {
    std::cerr << "TENES_FERMION_SECTOR step=-1 source=-1 leg="
              << connect_leg_a(direction) << " target=-1 Deven=" << dims.first
              << " Dodd=" << dims.second
              << " direction=" << direction_name(direction)
              << " context=full_update\n";
  }
  if ((dims.first == 0 || dims.second == 0) &&
      peps_parameters.print_level >= PrintLevel::warn &&
      U.t.get_comm_rank() == 0) {
    std::cerr << "warning: fermion full update kept an empty parity sector: "
              << "direction=" << direction_name(direction)
              << " Deven=" << dims.first << " Dodd=" << dims.second << "\n";
  }
}

template <class tensor>
std::vector<double> normalized_square_root_weights(
    const std::vector<double>& s) {
  double norm = 0.0;
  for (std::size_t i = 0; i < s.size(); ++i) {
    norm += s[i] * s[i];
  }
  norm = std::sqrt(norm);
  if (norm == 0.0) {
    throw std::runtime_error(
        "fermion Full_update_bond: zero singular-value norm");
  }
  std::vector<double> weights(s.size());
  for (std::size_t i = 0; i < s.size(); ++i) {
    weights[i] = std::sqrt(s[i] / norm);
  }
  return weights;
}

}  // namespace

template <class tensor>
void Full_update_bond_fermion(
    const tensor& C1, const tensor& C2, const tensor& C3, const tensor& C4,
    const tensor& eT1, const tensor& eT2, const tensor& eT3, const tensor& eT4,
    const tensor& eT5, const tensor& eT6,
    const tenes::fermion::ftensor<tensor>& Tn1,
    const tenes::fermion::ftensor<tensor>& Tn2,
    const tenes::fermion::ftensor<tensor>& wrapped_gate,
    tenes::fermion::reduced_pair_direction direction,
    const PEPS_Parameters& peps_parameters,
    tenes::fermion::ftensor<tensor>& Tn1_new,
    tenes::fermion::ftensor<tensor>& Tn2_new) {
  if (direction != tenes::fermion::reduced_pair_direction::horizontal &&
      direction != tenes::fermion::reduced_pair_direction::vertical) {
    throw std::runtime_error("Full_update_bond_fermion: invalid direction");
  }

  tenes::fermion::ftensor<tensor> QA, RA, QB, RB;
  int info = 0;
  if (direction == tenes::fermion::reduced_pair_direction::horizontal) {
    info = tenes::fermion::qr(Tn1, Axes(0, 1, 3), Axes(2, 4), QA, RA);
    info = tenes::fermion::qr(Tn2, Axes(1, 2, 3), Axes(0, 4), QB, RB);
  } else {
    info = tenes::fermion::qr(Tn1, Axes(0, 1, 2), Axes(3, 4), QA, RA);
    info = tenes::fermion::qr(Tn2, Axes(0, 2, 3), Axes(1, 4), QB, RB);
  }
  static_cast<void>(info);

  const auto env = tenes::fermion::build_full_update_environment(
      C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6, QA, QB, direction);
  if (fermion_full_update_log_enabled() && Tn1.t.get_comm_rank() == 0) {
    std::cerr << "TENES_FERMION_FULL_UPDATE_ENV direction="
              << direction_name(direction)
              << " N_forbidden_ratio=" << env.forbidden_ratio << "\n";
  }

  const auto X =
      tenes::fermion::transpose(tenes::fermion::tensordot(RA, RB, Axes(1),
                                                          Axes(1)),
                                Axes(0, 2, 1, 3));
  const auto Theta = tenes::fermion::tensordot(X, wrapped_gate, Axes(2, 3),
                                               Axes(0, 1));
  auto Theta_input_checked = Theta;
  project_or_throw(Theta_input_checked, "Theta_before_prepare");

  tensor Theta_tilde = Theta_input_checked.t;
  apply_pair_mask(Theta_tilde, Theta.parity[0], Theta.parity[1], 0, 1);

  tensor Env_out, Theta_out, LR1_inv, LR2_inv;
  prepare_environment(env.N_plain, Theta_tilde, peps_parameters, Env_out,
                      Theta_out, LR1_inv, LR2_inv);

  const tenes::fermion::parity_vector& pA = Theta.parity[0];
  const tenes::fermion::parity_vector& pB = Theta.parity[1];
  const tenes::fermion::parity_vector& ps1 = Theta.parity[2];
  const tenes::fermion::parity_vector& ps2 = Theta.parity[3];

  tenes::fermion::ftensor<tensor> Env_checked{
      Env_out, {pA, pB, pA, pB}};
  project_or_throw(Env_checked, "Env_out", false);
  Env_out = Env_checked.t;

  tenes::fermion::ftensor<tensor> Theta_checked{
      Theta_out, {pA, pB, ps1, ps2}};
  apply_pair_mask(Theta_checked.t, pA, pB, 0, 1);
  project_or_throw(Theta_checked, "Theta_out", false);
  apply_pair_mask(Theta_checked.t, pA, pB, 0, 1);
  Theta_out = Theta_checked.t;

  if (peps_parameters.Full_Gauge_Fix) {
    tenes::fermion::ftensor<tensor> LR1_checked{LR1_inv, {pA, pA}};
    tenes::fermion::ftensor<tensor> LR2_checked{LR2_inv, {pB, pB}};
    project_or_throw(LR1_checked, "LR1_inv", false);
    project_or_throw(LR2_checked, "LR2_inv", false);
    LR1_inv = LR1_checked.t;
    LR2_inv = LR2_checked.t;
  }

  tenes::fermion::ftensor<tensor> Theta_graded{
      Theta_out, {pA, pB, ps1, ps2}};
  apply_pair_mask(Theta_graded.t, pA, pB, 0, 1);

  tenes::fermion::ftensor<tensor> U, VT;
  std::vector<double> s;
  const int D_connect =
      static_cast<int>(Tn1.parity[connect_leg_a(direction)].size());
  info = tenes::fermion::svd_trunc(Theta_graded, Axes(0, 2), Axes(1, 3), U, s,
                                   VT, D_connect);
  static_cast<void>(info);
  log_sector_dimensions(U, direction, peps_parameters);

  std::vector<double> lambda_c =
      normalized_square_root_weights<tensor>(s);
  U.multiply_vector(lambda_c, 2);
  VT.multiply_vector(lambda_c, 0);
  auto R1 = tenes::fermion::transpose(U, Axes(0, 2, 1));
  auto R2 = tenes::fermion::transpose(VT, Axes(1, 0, 2));
  const tenes::fermion::parity_vector pm = U.parity[2];

  tensor R1_plain = R1.t;
  apply_pair_mask(R1_plain, pm, ps1, 1, 2);
  tensor R2_plain = R2.t;
  als_iterate(Env_out, Theta_out, peps_parameters, R1_plain, R2_plain);

  if (peps_parameters.Full_Gauge_Fix) {
    R1_plain = mptensor::tensordot(LR1_inv, R1_plain, Axes(0), Axes(0));
    R2_plain = mptensor::tensordot(LR2_inv, R2_plain, Axes(0), Axes(0));
  }

  R1 = tenes::fermion::ftensor<tensor>{R1_plain, {pA, pm, ps1}};
  apply_pair_mask(R1.t, pm, ps1, 1, 2);
  R2 = tenes::fermion::ftensor<tensor>{R2_plain, {pB, pm, ps2}};
  project_or_throw(R1, "R1_after_ALS");
  project_or_throw(R2, "R2_after_ALS");

  tenes::fermion::ftensor<tensor> q1, r1, q2, r2;
  info = tenes::fermion::qr(R1, Axes(0, 2), Axes(1), q1, r1);
  info = tenes::fermion::qr(R2, Axes(0, 2), Axes(1), q2, r2);

  tenes::fermion::ftensor<tensor> U2, VT2;
  std::vector<double> s2;
  info = tenes::fermion::svd(
      tenes::fermion::tensordot(r1, r2, Axes(1), Axes(1)), Axes(0), Axes(1),
      U2, s2, VT2);
  static_cast<void>(info);

  std::vector<double> lambda2 = normalized_square_root_weights<tensor>(s2);
  U2.multiply_vector(lambda2, 1);
  VT2.multiply_vector(lambda2, 0);
  R1 = tenes::fermion::tensordot(q1, U2, Axes(2), Axes(0));
  R2 = tenes::fermion::tensordot(q2, VT2, Axes(2), Axes(1));

  if (direction == tenes::fermion::reduced_pair_direction::horizontal) {
    Tn1_new =
        tenes::fermion::transpose(tenes::fermion::tensordot(QA, R1, Axes(3),
                                                            Axes(0)),
                                  Axes(0, 1, 4, 2, 3));
    Tn2_new =
        tenes::fermion::transpose(tenes::fermion::tensordot(QB, R2, Axes(3),
                                                            Axes(0)),
                                  Axes(4, 0, 1, 2, 3));
  } else {
    Tn1_new =
        tenes::fermion::transpose(tenes::fermion::tensordot(QA, R1, Axes(3),
                                                            Axes(0)),
                                  Axes(0, 1, 2, 4, 3));
    Tn2_new =
        tenes::fermion::transpose(tenes::fermion::tensordot(QB, R2, Axes(3),
                                                            Axes(0)),
                                  Axes(0, 4, 1, 2, 3));
  }

  project_or_throw(Tn1_new, "Tn1_new");
  project_or_throw(Tn2_new, "Tn2_new");
  if (Tn1_new.parity[connect_leg_a(direction)] !=
      Tn2_new.parity[connect_leg_b(direction)]) {
    throw std::runtime_error(
        "Full_update_bond_fermion: output bond parity ledgers differ");
  }
}

template void Full_update_bond_fermion(
    const real_tensor& C1, const real_tensor& C2, const real_tensor& C3,
    const real_tensor& C4, const real_tensor& eT1, const real_tensor& eT2,
    const real_tensor& eT3, const real_tensor& eT4, const real_tensor& eT5,
    const real_tensor& eT6,
    const tenes::fermion::ftensor<real_tensor>& Tn1,
    const tenes::fermion::ftensor<real_tensor>& Tn2,
    const tenes::fermion::ftensor<real_tensor>& wrapped_gate,
    tenes::fermion::reduced_pair_direction direction,
    const PEPS_Parameters& peps_parameters,
    tenes::fermion::ftensor<real_tensor>& Tn1_new,
    tenes::fermion::ftensor<real_tensor>& Tn2_new);

template void Full_update_bond_fermion(
    const complex_tensor& C1, const complex_tensor& C2, const complex_tensor& C3,
    const complex_tensor& C4, const complex_tensor& eT1,
    const complex_tensor& eT2, const complex_tensor& eT3,
    const complex_tensor& eT4, const complex_tensor& eT5,
    const complex_tensor& eT6,
    const tenes::fermion::ftensor<complex_tensor>& Tn1,
    const tenes::fermion::ftensor<complex_tensor>& Tn2,
    const tenes::fermion::ftensor<complex_tensor>& wrapped_gate,
    tenes::fermion::reduced_pair_direction direction,
    const PEPS_Parameters& peps_parameters,
    tenes::fermion::ftensor<complex_tensor>& Tn1_new,
    tenes::fermion::ftensor<complex_tensor>& Tn2_new);

}  // namespace tenes::itps::core
