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

#include "local_gauge.hpp"

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
#include <mptensor/tensor.hpp>

#include "../../tensor.hpp"
#include "../PEPS_Parameters.hpp"

#undef TENES_DEBUG_LOCAL_GAUGE_FIXING
// #define TENES_DEBUG_LOCAL_GAUGE_FIXING

namespace tenes {
namespace itps {
namespace core {

using mptensor::Axes;
using mptensor::Shape;

// environment

template <class tensor>
void fix_local_gauge(const tensor &Tn1, const tensor &Tn2,
                       const std::vector<std::vector<double>> &lambda1,
                       const std::vector<std::vector<double>> &lambda2,
                       const int connect1,
                       const PEPS_Parameters peps_parameters, tensor &Tn1_new,
                       tensor &Tn2_new, std::vector<double> &lambda_c) {
  using ptensor = tensor;
  constexpr int NLEG = 4;
  int connect2 = (connect1 + 2) % NLEG;


#ifdef TENES_DEBUG_LOCAL_GAUGE_FIXING
  std::cout << "Before fixing" << std::endl;
  {
    auto lambda = lambda1;
    lambda[connect1] = lambda_c;
    auto M = boundary_tensor(Tn1, lambda, connect1, peps_parameters);
    std::vector<double> D1;
    ptensor U1;
    eigh(M, Axes(0), Axes(1), D1, U1);
    std::cout << "D1 = ";
    for (auto d : D1) {
      std::cout << d << " ";
    }
    std::cout << std::endl;
  }
  {
    auto lambda = lambda2;
    lambda[connect2] = lambda_c;
    auto M = boundary_tensor(Tn2, lambda, connect2, peps_parameters);
    std::vector<double> D2;
    ptensor U2;
    eigh(M, Axes(0), Axes(1), D2, U2);
    std::cout << "D2 = ";
    for (auto d : D2) {
      std::cout << d << " ";
    }
    std::cout << std::endl;
  }
#endif

  std::vector<std::vector<double>> lambda1_inv(NLEG);
  std::vector<std::vector<double>> lambda2_inv(NLEG);

  for (int i = 0; i < NLEG; i++) {
    lambda1_inv[i] = std::vector<double>(lambda1[i].size());
    for (int j = 0; j < lambda1_inv[i].size(); ++j) {
      if (lambda1[i][j] > peps_parameters.Inverse_lambda_cut) {
        lambda1_inv[i][j] = 1.0 / lambda1[i][j];
      } else {
        lambda1_inv[i][j] = 0.0;
      }
    }
  }
  for (int i = 0; i < NLEG; i++) {
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
    Tn1_lambda.multiply_vector(lambda1_inv[0], 0, lambda1[1], 1, lambda1[2], 2,
                               lambda1[3], 3);
    Tn1_lambda.transpose(Axes(1, 2, 3, 0, 4));
  } else if (connect1 == 1) {
    Tn1_lambda.multiply_vector(lambda1[0], 0, lambda1_inv[1], 1, lambda1[2], 2,
                               lambda1[3], 3);
    Tn1_lambda.transpose(Axes(0, 2, 3, 1, 4));
  } else if (connect1 == 2) {
    Tn1_lambda.multiply_vector(lambda1[0], 0, lambda1[1], 1, lambda1_inv[2], 2,
                               lambda1[3], 3);
    Tn1_lambda.transpose(Axes(0, 1, 3, 2, 4));
  } else {
    Tn1_lambda.multiply_vector(lambda1[0], 0, lambda1[1], 1, lambda1[2], 2,
                               lambda1_inv[3], 3);
  }

  if (connect2 == 0) {
    Tn2_lambda.multiply_vector(lambda2_inv[0], 0, lambda2[1], 1, lambda2[2], 2,
                               lambda2[3], 3);
    Tn2_lambda.transpose(Axes(1, 2, 3, 0, 4));
  } else if (connect2 == 1) {
    Tn2_lambda.multiply_vector(lambda2[0], 0, lambda2_inv[1], 1, lambda2[2], 2,
                               lambda2[3], 3);
    Tn2_lambda.transpose(Axes(0, 2, 3, 1, 4));
  } else if (connect2 == 2) {
    Tn2_lambda.multiply_vector(lambda2[0], 0, lambda2[1], 1, lambda2_inv[2], 2,
                               lambda2[3], 3);
    Tn2_lambda.transpose(Axes(0, 1, 3, 2, 4));
  } else {
    Tn2_lambda.multiply_vector(lambda2[0], 0, lambda2[1], 1, lambda2[2], 2,
                               lambda2_inv[3], 3);
  }

  ptensor M1 = tensordot(conj(Tn1_lambda), Tn1_lambda, Axes(0, 1, 2, 4),
                         Axes(0, 1, 2, 4));
  ptensor M2 = tensordot(conj(Tn2_lambda), Tn2_lambda, Axes(0, 1, 2, 4),
                         Axes(0, 1, 2, 4));

  ptensor U1, U2;
  std::vector<double> D1, D2;
  eigh(M1, Axes(0), Axes(1), D1, U1);
  eigh(M2, Axes(0), Axes(1), D2, U2);

  std::vector<double> D1inv(dc), D2inv(dc);
  for (int d = 0; d < dc; ++d) {
    D1[d] = sqrt(D1[d]);
    D1inv[d] = 1.0 / D1[d];

    D2[d] = sqrt(D2[d]);
    D2inv[d] = 1.0 / D2[d];
  }

  ptensor U1_ = U1;
  ptensor U2_ = U2;
  U1_.multiply_vector(lambda1[connect1], 0, D1, 1);
  U2_.multiply_vector(lambda2[connect2], 0, D2, 1);
  ptensor L = tensordot(conj(U1_), conj(U2_), Axes(0), Axes(0));

  ptensor W1, W2;
  svd(L, Axes(0), Axes(1), W1, lambda_c, W2);

  double norm = 0.0;
  for (int d = dc - 1; d >= 0; --d) {
    norm += lambda_c[d] * lambda_c[d];
  }
  norm = sqrt(norm);

  for (int d = 0; d < dc; ++d) {
    lambda_c[d] = sqrt(lambda_c[d] / norm);
  }

  W1.multiply_vector(D1inv, 0);
  W2.multiply_vector(D2inv, 1);

  ptensor X1inv = tensordot(U1, W1, Axes(1), Axes(0));
  Tn1_new = tensordot(Tn1_lambda, X1inv, Axes(3), Axes(0));

  ptensor X2inv = tensordot(U2, W2, Axes(1), Axes(1));
  Tn2_new = tensordot(Tn2_lambda, X2inv, Axes(3), Axes(0));

  if (connect1 == 0) {
    Tn1_new.transpose(Axes(4, 0, 1, 2, 3));
    Tn1_new.multiply_vector(lambda_c, 0, lambda1_inv[1], 1, lambda1_inv[2], 2,
                            lambda1_inv[3], 3);
  } else if (connect1 == 1) {
    Tn1_new.transpose(Axes(0, 4, 1, 2, 3));
    Tn1_new.multiply_vector(lambda1_inv[0], 0, lambda_c, 1, lambda1_inv[2], 2,
                            lambda1_inv[3], 3);
  } else if (connect1 == 2) {
    Tn1_new.transpose(Axes(0, 1, 4, 2, 3));
    Tn1_new.multiply_vector(lambda1_inv[0], 0, lambda1_inv[1], 1, lambda_c, 2,
                            lambda1_inv[3], 3);
  } else {
    Tn1_new.transpose(Axes(0, 1, 2, 4, 3));
    Tn1_new.multiply_vector(lambda1_inv[0], 0, lambda1_inv[1], 1,
                            lambda1_inv[2], 2, lambda_c, 3);
  }

  if (connect2 == 0) {
    Tn2_new.transpose(Axes(4, 0, 1, 2, 3));
    Tn2_new.multiply_vector(lambda_c, 0, lambda2_inv[1], 1, lambda2_inv[2], 2,
                            lambda2_inv[3], 3);
  } else if (connect2 == 1) {
    Tn2_new.transpose(Axes(0, 4, 1, 2, 3));
    Tn2_new.multiply_vector(lambda2_inv[0], 0, lambda_c, 1, lambda2_inv[2], 2,
                            lambda2_inv[3], 3);
  } else if (connect2 == 2) {
    Tn2_new.transpose(Axes(0, 1, 4, 2, 3));
    Tn2_new.multiply_vector(lambda2_inv[0], 0, lambda2_inv[1], 1, lambda_c, 2,
                            lambda2_inv[3], 3);
  } else {
    Tn2_new.transpose(Axes(0, 1, 2, 4, 3));
    Tn2_new.multiply_vector(lambda2_inv[0], 0, lambda2_inv[1], 1,
                            lambda2_inv[2], 2, lambda_c, 3);
  }

#ifdef TENES_DEBUG_LOCAL_GAUGE_FIXING
  std::cout << "After fixing" << std::endl;
  {
    auto lambda = lambda1;
    lambda[connect1] = lambda_c;
    auto M = boundary_tensor(Tn1_new, lambda, connect1, peps_parameters);
    eigh(M, Axes(0), Axes(1), D1, U1);
    std::cout << "D1 = ";
    for (auto d : D1) {
      std::cout << d << " ";
    }
    std::cout << std::endl;
  }
  {
    auto lambda = lambda2;
    lambda[connect2] = lambda_c;
    auto M = boundary_tensor(Tn2_new, lambda, connect2, peps_parameters);
    eigh(M, Axes(0), Axes(1), D2, U2);
    std::cout << "D2 = ";
    for (auto d : D2) {
      std::cout << d << " ";
    }
    std::cout << std::endl;
  }
#endif
}

template <class tensor>
tensor boundary_tensor(const tensor &Tn,
                       const std::vector<std::vector<double>> &lambda,
                       const int connect,
                       const PEPS_Parameters peps_parameters) {
  tensor Tn_lambda = Tn;
  const size_t dc = Tn.shape()[connect];
  std::vector<double> lambda_inv(dc);
  for (size_t i = 0; i < dc; ++i) {
    if (lambda[connect][i] > peps_parameters.Inverse_lambda_cut) {
      lambda_inv[i] = 1.0 / lambda[connect][i];
    } else {
      lambda_inv[i] = 0.0;
    }
  }

  if (connect == 0) {
    Tn_lambda.multiply_vector(lambda_inv, 0, lambda[1], 1, lambda[2], 2,
                              lambda[3], 3);
    Tn_lambda.transpose(Axes(1, 2, 3, 0, 4));
  } else if (connect == 1) {
    Tn_lambda.multiply_vector(lambda[0], 0, lambda_inv, 1, lambda[2], 2,
                              lambda[3], 3);
    Tn_lambda.transpose(Axes(0, 2, 3, 1, 4));
  } else if (connect == 2) {
    Tn_lambda.multiply_vector(lambda[0], 0, lambda[1], 1, lambda_inv, 2,
                              lambda[3], 3);
    Tn_lambda.transpose(Axes(0, 1, 3, 2, 4));
  } else {
    Tn_lambda.multiply_vector(lambda[0], 0, lambda[1], 1, lambda[2], 2,
                              lambda_inv, 3);
  }

  return tensordot(conj(Tn_lambda), Tn_lambda, Axes(0, 1, 2, 4),
                   Axes(0, 1, 2, 4));
}

// template instantiations

template void fix_local_gauge(const real_tensor &Tn1, const real_tensor &Tn2,
                              const std::vector<std::vector<double>> &lambda1,
                              const std::vector<std::vector<double>> &lambda2,
                              const int connect1,
                              const PEPS_Parameters peps_parameters,
                              real_tensor &Tn1_new, real_tensor &Tn2_new,
                              std::vector<double> &lambda_c);

template void fix_local_gauge(const complex_tensor &Tn1,
                              const complex_tensor &Tn2,
                              const std::vector<std::vector<double>> &lambda1,
                              const std::vector<std::vector<double>> &lambda2,
                              const int connect1,
                              const PEPS_Parameters peps_parameters,
                              complex_tensor &Tn1_new, complex_tensor &Tn2_new,
                              std::vector<double> &lambda_c);

template real_tensor boundary_tensor(
    const real_tensor &Tn, const std::vector<std::vector<double>> &lambda,
    const int connect, const PEPS_Parameters peps_parameters);

template complex_tensor boundary_tensor(
    const complex_tensor &Tn, const std::vector<std::vector<double>> &lambda,
    const int connect, const PEPS_Parameters peps_parameters);

}  // end of namespace core
}  // namespace itps
}  // namespace tenes
