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

#ifndef TENES_SRC_ARNOLDI_HPP_
#define TENES_SRC_ARNOLDI_HPP_

#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <vector>
#include <complex>
#include <functional>
#include <random>
#include <tuple>

#include "util/abs.hpp"
#include "tensor.hpp"

namespace tenes {

template <class ptensor>
double norm(ptensor const &A) {
  return std::sqrt(std::real(trace(tensordot(conj(A), A, {0}, {0}), {}, {})));
}

template <class ptensor>
struct Arnoldi {
  using value_type = typename ptensor::value_type;
  using Eigenmatrix = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>;

  size_t N;
  size_t maxvec;
  size_t nev;

  std::vector<ptensor> Q;
  Eigenmatrix H;

  Arnoldi(size_t N, size_t maxvec)
      : N(N), maxvec(maxvec), Q(maxvec + 1), H(maxvec + 1, maxvec) {
    for (size_t i = 0; i < maxvec; ++i) {
      Q[i] = ptensor(mptensor::Shape(N));
    }
  }

  void initialize(ptensor const &initial) {
    Q[0] = initial;
    orthonormalize(0);
  }

  void run(std::function<void(ptensor &, ptensor const &)> A, size_t nev,
           int maxiter = 10) {
    this->nev = nev;
    for (size_t k = 1; k <= nev; ++k) {
      A(Q[k], Q[k - 1]);
      orthonormalize(k);
    }
    for (size_t iter = 1; iter <= maxiter; ++iter) {
      for (size_t k = nev + 1; k <= maxvec; ++k) {
        A(Q[k], Q[k - 1]);
        orthonormalize(k);
      }
      if (iter < maxiter) {
        restart(nev);
      }
    }
  }

  std::vector<std::complex<double>> eigenvalues() const {
    auto evs = H.topLeftCorner(maxvec, maxvec).eigenvalues();
    std::vector<std::complex<double>> ret(maxvec);
    for (size_t i = 0; i < maxvec; ++i) {
      ret[i] = evs[i];
    }
    std::partial_sort(ret.begin(), ret.begin() + nev, ret.end(),
                      [](std::complex<double> a, std::complex<double> b) {
                        return util::abs2(a) > util::abs2(b);
                      });
    ret.erase(ret.begin() + nev, ret.end());
    return ret;
  }

 private:
  void orthonormalize(size_t k) {
    for (size_t j = 0; j < k; ++j) {
      auto v = tensordot(conj(Q[j]), Q[k], {0}, {0});
      auto x = trace(v, mptensor::Axes(), mptensor::Axes());
      H(j, k - 1) = x;
      Q[k] -= Q[j] * x;
    }
    auto x = norm(Q[k]);
    Q[k] *= 1.0 / x;
    if (k > 0) {
      H(k, k - 1) = x;
    }
  }

  void restart(size_t minvec, double cutoff = 1.0e-13) {
    size_t k = maxvec;
    auto h = H.topLeftCorner(k, k);
    auto eigvals = h.eigenvalues();
    std::vector<std::complex<double>> evals(k);
    for (size_t i = 0; i < k; ++i) {
      evals[i] = eigvals[i];
    }
    std::sort(evals.begin(), evals.end(),
              [](std::complex<double> a, std::complex<double> b) {
                return util::abs2(a) < util::abs2(b);
              });

    for (size_t i = 0; i < k - minvec; ++i) {
      auto sigma = std::real(evals[i]);
      Eigenmatrix A = h - sigma * Eigenmatrix::Identity(k, k);
      auto qr = A.householderQr();
      auto q = qr.householderQ();
      h = (q.adjoint() * h * q).eval();

      size_t n = 0;
      const size_t max_block_size = 1;
      while (n < N) {
        size_t block_size = std::min(max_block_size, N - n);
        Eigenmatrix newq(block_size, k);
        for (size_t j = 0; j < k; ++j) {
          for (size_t b = 0; b < block_size; ++b) {
            value_type temp;
            Q[j].get_value({n + b}, temp);
            newq(b, j) = temp;
          }
        }
        newq = (newq * q).eval();
        for (size_t j = 0; j < k; ++j) {
          for (size_t b = 0; b < block_size; ++b) {
            value_type temp = newq(b, j);
            Q[j].set_value({n + b}, temp);
          }
        }
        n += block_size;
      }
    }
    H.topLeftCorner(k, k) = h;
  }
};

}  // end of namespace tenes

#endif  // TENES_SRC_ARNOLDI_HPP_
