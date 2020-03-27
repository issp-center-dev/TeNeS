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

#include <fstream>
#include <vector>

#include <PEPS_Basics.hpp>
#include <PEPS_Parameters.cpp>
#include <mpi.cpp>

TEST_CASE("testing simple update") {
#ifdef _NO_MPI
  using tensor = mptensor::Tensor<mptensor::lapack::Matrix, double>;
#else
  using tensor = mptensor::Tensor<mptensor::scalapack::Matrix, double>;
#endif

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
  std::vector<std::vector<double>> &lambda_1 = lambda[0];
  std::vector<std::vector<double>> &lambda_2 = lambda[1];

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

  tenes::PEPS_Parameters peps_parameters;
  int connect = 2;
  std::vector<tensor> new_T(2);
  std::vector<double> new_lambda;

  Simple_update_bond(T[0], T[1], lambda_1, lambda_2, op, connect,
                     peps_parameters, new_T[0], new_T[1], new_lambda);

  // check results

  std::ofstream ofs("res_simple_update.dat");
  ofs << std::setprecision(std::numeric_limits<double>::digits10);

  int sign = 0;
  for (int a = 0; a < 2; ++a) {
    for (int i = 0; i < D; ++i){
      for (int j = 0; j < D; ++j){
        for (int k = 0; k < D; ++k){
          for (int l = 0; l < D; ++l){
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
            }}}}}
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
