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

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "test_workdir.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <vector>

#include "../src/fermion/fops.hpp"
#include "../src/fermion/ftensor.hpp"
#include "../src/mpi.hpp"
#include "../src/tensor.hpp"
#include "../src/iTPS/PEPS_Parameters.hpp"
#include "../src/iTPS/core/simple_update.hpp"

TEST_CASE("testing simple update") {
  using tensor = tenes::real_tensor;

  using mptensor::Index;
  using mptensor::Shape;

  const int ldof = 2;
  const int D = 2;
  const int nleg = 4;

  const double tol = 1.0e-8;

  // make inputs

  std::vector<tensor> T(2, tensor(Shape(D, D, D, D, ldof)));

  std::vector<std::vector<std::vector<double>>> lambda;
  for (int i = 0; i < 2; ++i) {
    lambda.push_back(
        std::vector<std::vector<double>>(nleg, std::vector<double>(D, 1.0)));
  }
  std::vector<std::vector<double>>& lambda_1 = lambda[0];
  std::vector<std::vector<double>>& lambda_2 = lambda[1];

  for (int i = 0; i < D; ++i)
    for (int j = 0; j < D; ++j)
      for (int k = 0; k < D; ++k)
        for (int l = 0; l < D; ++l)
          for (int m = 0; m < ldof; ++m) {
            for (int a = 0; a < 2; ++a)
              T[a].set_value(Index(i, j, k, l, m), 1.0);
          }

  tensor op(Shape(ldof, ldof, ldof, ldof));

  for (int i = 0; i < ldof; ++i)
    for (int j = 0; j < ldof; ++j)
      for (int k = 0; k < ldof; ++k)
        for (int l = 0; l < ldof; ++l) {
          op.set_value(Index(i, j, k, l), 0.0);
        }
  for (int i = 0; i < ldof; ++i) {
    op.set_value(Index(i, i, i, i), 1.0);
  }

  // load answer

  std::ifstream ifs("data/simple_update.dat");

  std::vector<tensor> ans_T(2, tensor(Shape(D, D, D, D, ldof)));
  std::vector<double> ans_lambda(D, 0.0);

  for (int a = 0; a < 2; ++a)
    for (int i = 0; i < D; ++i)
      for (int j = 0; j < D; ++j)
        for (int k = 0; k < D; ++k)
          for (int l = 0; l < D; ++l)
            for (int m = 0; m < ldof; ++m) {
              double val;
              ifs >> val;
              ans_T[a].set_value(Index(i, j, k, l, m), val);
            }
  for (int i = 0; i < D; ++i) {
    double val;
    ifs >> val;
    ans_lambda[i] = val;
  }

  // calculation

  tenes::itps::PEPS_Parameters peps_parameters;
  int connect = 2;
  std::vector<tensor> new_T(2);
  std::vector<double> new_lambda;

  tenes::itps::core::Simple_update_bond(T[0], T[1], lambda_1, lambda_2, op,
                                        connect, peps_parameters, new_T[0],
                                        new_T[1], new_lambda);

  // check results

  std::ofstream ofs("res_simple_update.dat");
  ofs << std::setprecision(std::numeric_limits<double>::digits10);

  int sign = 0;
  for (int a = 0; a < 2; ++a) {
    for (int i = 0; i < D; ++i) {
      for (int j = 0; j < D; ++j) {
        for (int k = 0; k < D; ++k) {
          for (int l = 0; l < D; ++l) {
            sign = 0;
            for (int m = 0; m < ldof; ++m) {
              double result, answer;
              new_T[a].get_value(Index(i, j, k, l, m), result);
              ans_T[a].get_value(Index(i, j, k, l, m), answer);
              if (sign == 0) {
                if (result != 0.0) {
                  if (answer * result > 0.0) {
                    sign = 1;
                  } else {
                    sign = -1;
                  }
                }
              }
              CHECK(result * sign == doctest::Approx(answer).epsilon(tol));
              ofs << result << " ";
            }
          }
        }
      }
    }
    ofs << std::endl;
  }

  sign = 0;

  for (int i = 0; i < D; ++i) {
    double result = new_lambda[i];
    double answer = ans_lambda[i];
    if (sign == 0) {
      if (result != 0.0) {
        if (answer * result > 0.0) {
          sign = 1;
        } else {
          sign = -1;
        }
      }
    }
    CHECK(result == doctest::Approx(answer).epsilon(tol));
    ofs << new_lambda[i] << " ";
  }
  ofs << std::endl;
}

namespace {
using real_ftensor = tenes::fermion::ftensor<tenes::real_tensor>;

tenes::real_tensor make_even_Tn(int seed) {
  tenes::real_tensor t(mptensor::Shape(2, 2, 2, 2, 2));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const auto idx = t.global_index(n);
    const int odd = (idx[0] + idx[1] + idx[2] + idx[3] + idx[4]) % 2;
    if (odd == 0) {
      const double x = static_cast<double>((seed + 3) * (n + 1));
      t.set_value(idx, 0.37 * std::sin(x) + 0.19 * std::cos(0.7 * x));
    }
  }
  return t;
}

tenes::real_tensor make_hopping_gate() {
  tenes::real_tensor op(mptensor::Shape(2, 2, 2, 2));
  op.set_value(mptensor::Index(0, 0, 0, 0), 1.0);
  op.set_value(mptensor::Index(0, 1, 0, 1), 1.0);
  op.set_value(mptensor::Index(1, 0, 1, 0), 1.0);
  op.set_value(mptensor::Index(1, 1, 1, 1), 1.0);
  op.set_value(mptensor::Index(0, 1, 1, 0), 0.1);
  op.set_value(mptensor::Index(1, 0, 0, 1), 0.1);
  return op;
}

tenes::fermion::leg_parities Tn_parities(bool all_even = false) {
  const tenes::fermion::parity_vector p =
      all_even ? tenes::fermion::parity_vector{false, false}
               : tenes::fermion::parity_vector{false, true};
  return {p, p, p, p, p};
}

tenes::real_tensor make_even_Tn_with_parity(
    int seed, const tenes::fermion::leg_parities& parity) {
  tenes::real_tensor t(mptensor::Shape(2, 2, 2, 2, 2));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const auto idx = t.global_index(n);
    if (tenes::fermion::count_odd(parity, idx) % 2 == 0) {
      const double x = static_cast<double>((seed + 3) * (n + 1));
      t.set_value(idx, 0.37 * std::sin(x) + 0.19 * std::cos(0.7 * x));
    }
  }
  return t;
}

tenes::fermion::leg_parities op_parities(bool all_even = false) {
  const tenes::fermion::parity_vector p =
      all_even ? tenes::fermion::parity_vector{false, false}
               : tenes::fermion::parity_vector{false, true};
  return {p, p, p, p};
}

tenes::real_tensor make_hopping_gate(double hopping) {
  auto op = make_hopping_gate();
  op.set_value(mptensor::Index(0, 1, 1, 0), hopping);
  op.set_value(mptensor::Index(1, 0, 0, 1), hopping);
  return op;
}

std::vector<double> sorted_spectrum(std::vector<double> values) {
  std::sort(values.begin(), values.end());
  return values;
}

std::vector<double> run_boson_lambda_spectrum(double ty) {
  const auto virtual_parity = (ty == 0.0)
                                  ? tenes::fermion::parity_vector{false, false}
                                  : tenes::fermion::parity_vector{false, true};
  const auto physical_parity = (ty == 0.0)
                                   ? tenes::fermion::parity_vector{false, false}
                                   : tenes::fermion::parity_vector{false, true};
  const tenes::fermion::leg_parities init_parity{virtual_parity, virtual_parity,
                                                 virtual_parity, virtual_parity,
                                                 physical_parity};
  std::vector<tenes::real_tensor> T{make_even_Tn_with_parity(1, init_parity),
                                    make_even_Tn_with_parity(3, init_parity),
                                    make_even_Tn_with_parity(5, init_parity),
                                    make_even_Tn_with_parity(7, init_parity)};
  std::vector<std::vector<std::vector<double>>> lambda(
      4, std::vector<std::vector<double>>(4, std::vector<double>{1.0, 1.0}));
  tenes::itps::PEPS_Parameters peps_parameters;
  peps_parameters.Inverse_lambda_cut = 1.0e-12;
  const auto hop_x = make_hopping_gate(0.12);
  const auto hop_y = make_hopping_gate(ty);

  auto update = [&](int source, int target, int leg,
                    const tenes::real_tensor& op) {
    tenes::real_tensor out0, out1;
    std::vector<double> bond_lambda;
    tenes::itps::core::Simple_update_bond(
        T[source], T[target], lambda[source], lambda[target], op, leg,
        peps_parameters, out0, out1, bond_lambda);
    T[source] = out0;
    T[target] = out1;
    lambda[source][leg] = bond_lambda;
    lambda[target][(leg + 2) % 4] = bond_lambda;
  };

  if (ty != 0.0) {
    update(0, 2, 1, hop_y);
    update(1, 3, 1, hop_y);
  }
  update(0, 1, 2, hop_x);
  return sorted_spectrum(lambda[0][2]);
}

std::vector<double> run_fermion_lambda_spectrum(double ty) {
  const auto virtual_parity = (ty == 0.0)
                                  ? tenes::fermion::parity_vector{false, false}
                                  : tenes::fermion::parity_vector{false, true};
  const auto physical_parity = (ty == 0.0)
                                   ? tenes::fermion::parity_vector{false, false}
                                   : tenes::fermion::parity_vector{false, true};
  const tenes::fermion::leg_parities init_parity{virtual_parity, virtual_parity,
                                                 virtual_parity, virtual_parity,
                                                 physical_parity};
  std::vector<tenes::real_tensor> T{make_even_Tn_with_parity(1, init_parity),
                                    make_even_Tn_with_parity(3, init_parity),
                                    make_even_Tn_with_parity(5, init_parity),
                                    make_even_Tn_with_parity(7, init_parity)};
  std::vector<std::array<tenes::fermion::parity_vector, 4>> virt(
      4, {virtual_parity, virtual_parity, virtual_parity, virtual_parity});
  const std::vector<tenes::fermion::parity_vector> phys(4, physical_parity);
  std::vector<std::vector<std::vector<double>>> lambda(
      4, std::vector<std::vector<double>>(4, std::vector<double>{1.0, 1.0}));
  tenes::itps::PEPS_Parameters peps_parameters;
  peps_parameters.Inverse_lambda_cut = 1.0e-12;
  const auto hop_x = make_hopping_gate(0.12);
  const auto hop_y = make_hopping_gate(ty);

  auto wrap_T = [&](int site) {
    return real_ftensor{T[site],
                        {virt[site][0], virt[site][1], virt[site][2],
                         virt[site][3], phys[site]}};
  };
  auto update = [&](int source, int target, int leg,
                    const tenes::real_tensor& op) {
    real_ftensor out0, out1;
    std::vector<double> bond_lambda;
    real_ftensor fop{op,
                     {phys[source], phys[target], phys[source], phys[target]}};
    tenes::itps::core::Simple_update_bond(
        wrap_T(source), wrap_T(target), lambda[source], lambda[target], fop,
        leg, peps_parameters, out0, out1, bond_lambda);
    T[source] = out0.t;
    T[target] = out1.t;
    virt[source][leg] = out0.parity[leg];
    virt[target][(leg + 2) % 4] = out1.parity[(leg + 2) % 4];
    lambda[source][leg] = bond_lambda;
    lambda[target][(leg + 2) % 4] = bond_lambda;
  };

  if (ty != 0.0) {
    update(0, 2, 1, hop_y);
    update(1, 3, 1, hop_y);
  }
  update(0, 1, 2, hop_x);
  return sorted_spectrum(lambda[0][2]);
}
}  // namespace

TEST_CASE("fermion Simple_update_bond preserves even parity") {
  const auto T0 = make_even_Tn(1);
  const auto T1 = make_even_Tn(7);
  const auto op = make_hopping_gate();
  const std::vector<std::vector<double>> lambda(4,
                                                std::vector<double>{1.0, 1.0});
  tenes::itps::PEPS_Parameters peps_parameters;
  peps_parameters.Inverse_lambda_cut = 1.0e-12;

  real_ftensor fT0{T0, Tn_parities()};
  real_ftensor fT1{T1, Tn_parities()};
  real_ftensor fop{op, op_parities()};
  real_ftensor out0, out1;
  std::vector<double> lambda_c;
  tenes::itps::core::Simple_update_bond(fT0, fT1, lambda, lambda, fop, 2,
                                        peps_parameters, out0, out1, lambda_c);

  CHECK(tenes::fermion::parity_violation(out0) == doctest::Approx(0.0));
  CHECK(tenes::fermion::parity_violation(out1) == doctest::Approx(0.0));
  REQUIRE(lambda_c.size() == 2);
  for (double x : lambda_c) {
    CHECK(std::isfinite(x));
    CHECK(x >= 0.0);
  }
}

TEST_CASE("parity-trivial ftensor Simple_update_bond matches bosonic kernel") {
  const auto T0 = make_even_Tn(3);
  const auto T1 = make_even_Tn(11);
  const auto op = make_hopping_gate();
  const std::vector<std::vector<double>> lambda(4,
                                                std::vector<double>{1.0, 1.0});
  tenes::itps::PEPS_Parameters peps_parameters;
  peps_parameters.Inverse_lambda_cut = 1.0e-12;

  tenes::real_tensor boson_out0, boson_out1;
  std::vector<double> boson_lambda;
  tenes::itps::core::Simple_update_bond(T0, T1, lambda, lambda, op, 2,
                                        peps_parameters, boson_out0, boson_out1,
                                        boson_lambda);

  real_ftensor fT0{T0, Tn_parities(true)};
  real_ftensor fT1{T1, Tn_parities(true)};
  real_ftensor fop{op, op_parities(true)};
  real_ftensor fermion_out0, fermion_out1;
  std::vector<double> fermion_lambda;
  tenes::itps::core::Simple_update_bond(fT0, fT1, lambda, lambda, fop, 2,
                                        peps_parameters, fermion_out0,
                                        fermion_out1, fermion_lambda);

  REQUIRE(fermion_lambda.size() == boson_lambda.size());
  for (std::size_t i = 0; i < boson_lambda.size(); ++i) {
    CHECK(fermion_lambda[i] == doctest::Approx(boson_lambda[i]).epsilon(1e-12));
  }
  for (std::size_t n = 0; n < boson_out0.local_size(); ++n) {
    const auto idx = boson_out0.global_index(n);
    double b0, f0, b1, f1;
    boson_out0.get_value(idx, b0);
    fermion_out0.t.get_value(idx, f0);
    boson_out1.get_value(idx, b1);
    fermion_out1.t.get_value(idx, f1);
    CHECK(f0 == doctest::Approx(b0).epsilon(1e-12));
    CHECK(f1 == doctest::Approx(b1).epsilon(1e-12));
  }
}

TEST_CASE(
    "fermion simple update lambda spectrum matches the 1D hard-core "
    "boson limit only") {
  const auto boson_1d = run_boson_lambda_spectrum(0.0);
  const auto fermion_1d = run_fermion_lambda_spectrum(0.0);
  REQUIRE(boson_1d.size() == fermion_1d.size());
  for (std::size_t i = 0; i < boson_1d.size(); ++i) {
    CHECK(fermion_1d[i] == doctest::Approx(boson_1d[i]).epsilon(1.0e-8));
  }

  const auto boson_2d = run_boson_lambda_spectrum(0.12);
  const auto fermion_2d = run_fermion_lambda_spectrum(0.12);
  REQUIRE(boson_2d.size() == fermion_2d.size());
  double max_diff = 0.0;
  for (std::size_t i = 0; i < boson_2d.size(); ++i) {
    max_diff = std::max(max_diff, std::abs(fermion_2d[i] - boson_2d[i]));
  }
  CHECK(max_diff > 1.0e-8);
}
