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

#include <fstream>
#include <vector>

#include "../src/tensor.hpp"
#include "../src/PEPS_Parameters.hpp"
#include "../src/mpi.hpp"
#include "../src/iTPS/core/full_update.hpp"

#include "doctest.h"

TEST_CASE("testing full update") {
  using tensor = tenes::real_tensor;

  using mptensor::Index;
  using mptensor::Shape;

  const int ldof = 2;
  const int D = 2;
  const int chi = 4;
  const int nleg = 4;

  const double tol = 1.0e-8;

  // make inputs

  std::vector<tensor> T(2, tensor(Shape(D, D, D, D, ldof)));
  std::vector<tensor> C(4, tensor(Shape(chi, chi)));
  std::vector<tensor> E(6, tensor(Shape(chi, chi, D, D)));
  tensor op(Shape(ldof, ldof, ldof, ldof));

  for (int a = 0; a < 2; ++a)
    for (int i = 0; i < D; ++i)
      for (int j = 0; j < D; ++j)
        for (int k = 0; k < D; ++k)
          for (int l = 0; l < D; ++l)
            for (int m = 0; m < ldof; ++m) {
              T[a].set_value(Index(i, j, k, l, m), 1.0);
            }

  for (int a = 0; a < 4; ++a)
    for (int i = 0; i < D; ++i)
      for (int j = 0; j < D; ++j) {
        C[a].set_value(Index(i, j), 1.0);
      }

  for (int a = 0; a < 6; ++a)
    for (int i = 0; i < chi; ++i)
      for (int j = 0; j < chi; ++j)
        for (int k = 0; k < D; ++k)
          for (int l = 0; l < D; ++l) {
            E[a].set_value(Index(i, j, k, l), 1.0);
          }

  for (int i = 0; i < ldof; ++i)
    for (int j = 0; j < ldof; ++j)
      for (int k = 0; k < ldof; ++k)
        for (int l = 0; l < ldof; ++l) {
          op.set_value(Index(i, j, k, l), 0.0);
        }
  for (int i = 0; i < ldof; ++i) {
    op.set_value(Index(i, i, i, i), 1.0);
  }

  // calculation

  tenes::PEPS_Parameters peps_parameters;
  int connect = 2;
  std::vector<tensor> new_T(2);
  std::vector<double> new_lambda;

  tenes::Full_update_bond(C[0], C[1], C[2], C[3], E[0], E[1], E[2], E[3], E[4],
                          E[5], T[0], T[1], op, 2, peps_parameters, new_T[0],
                          new_T[1]);

  // load answer

  std::ifstream ifs("data/full_update.dat");

  std::vector<tensor> ans_T(2, tensor(Shape(D, D, D, D, ldof)));

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

  // check results

  std::ofstream ofs("res_full_update.dat");
  ofs << std::setprecision(std::numeric_limits<double>::digits10);

  std::vector<int> sign(ldof);
  for (int a = 0; a < 2; ++a) {
    std::fill(sign.begin(), sign.end(), 0.0);
    for (int i = 0; i < D; ++i)
      for (int j = 0; j < D; ++j)
        for (int k = 0; k < D; ++k)
          for (int l = 0; l < D; ++l)
            for (int m = 0; m < ldof; ++m) {
              double result, answer;
              new_T[a].get_value(Index(i, j, k, l, m), result);
              ans_T[a].get_value(Index(i, j, k, l, m), answer);
              if (sign[m] == 0) {
                if (result != 0.0) {
                  if (answer * result > 0.0) {
                    sign[m] = 1;
                  } else {
                    sign[m] = -1;
                  }
                }
              }
              CHECK(result * sign[m] == doctest::Approx(answer).epsilon(tol));
            }
    ofs << std::endl;
  }

  ofs << std::endl;
}
