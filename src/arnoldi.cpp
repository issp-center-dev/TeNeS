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

#include "arnoldi.hpp"

#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <vector>
#include <complex>
#include <functional>
#include <type_traits>
#include <numeric>

#include "util/abs.hpp"
#include "tensor.hpp"

namespace tenes {

template <class ptensor>
Arnoldi<ptensor>::Arnoldi(size_t N, size_t maxvec)
    : N(N), maxvec(maxvec), Q(maxvec + 1), H(maxvec + 1, maxvec) {
  for (size_t i = 0; i < maxvec; ++i) {
    Q[i] = ptensor(mptensor::Shape(N));
  }
}

template <class ptensor>
void Arnoldi<ptensor>::initialize(ptensor const &initial) {
  Q[0] = initial;
  orthonormalize(0);
}

template <class ptensor>
double norm(ptensor const &A) {
  return std::sqrt(std::real(trace(tensordot(conj(A), A, {0}, {0}), {}, {})));
}

template <class ptensor>
void Arnoldi<ptensor>::run(std::function<void(ptensor &, ptensor const &)> A,
                           size_t nev, int mindim, int maxiter, double rtol) {
  if (mindim < nev) {
    mindim = nev;
  }
  this->nev = nev;
  for (size_t k = 1; k <= nev; ++k) {
    A(Q[k], Q[k - 1]);
    orthonormalize(k);
  }
  double res_nev = 1.0;
  for (size_t iter = 1; iter <= maxiter; ++iter) {
    for (size_t k = nev + 1; k <= maxvec; ++k) {
      A(Q[k], Q[k - 1]);
      orthonormalize(k);
      auto res = residue(k);
      res_nev = res[nev - 1];
      if (res_nev <= rtol) {
        break;
      }
    }
    if (res_nev <= rtol || iter == maxiter) {
      break;
    }
    restart(mindim);
  }
}

template <class ptensor>
std::vector<std::complex<double>> Arnoldi<ptensor>::eigenvalues() const {
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

template <class ptensor>
void Arnoldi<ptensor>::orthonormalize(size_t k) {
  for (size_t j = 0; j < k; ++j) {
    auto x = trace(tensordot(conj(Q[j]), Q[k], {0}, {0}), {}, {});
    H(j, k - 1) = x;
    Q[k] -= Q[j] * x;
  }
  auto x = norm(Q[k]);
  Q[k] *= 1.0 / x;
  if (k > 0) {
    H(k, k - 1) = x;
  }
}

double real_if_real(std::complex<double> v, double) { return std::real(v); }
std::complex<double> real_if_real(std::complex<double> v,
                                  std::complex<double>) {
  return v;
}

template <class ptensor>
void Arnoldi<ptensor>::restart(size_t minvec, double cutoff) {
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
    value_type sigma = real_if_real(evals[i], value_type());
    Eigenmatrix A = h - sigma * Eigenmatrix::Identity(k, k);
    auto qr = A.householderQr();
    auto q = qr.householderQ();
    h = (q.adjoint() * h * q).eval();

    size_t n = 0;
    constexpr size_t max_block_size = 10;
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

  cutoff *= cutoff;
  for (size_t col = 0; col < k; ++col) {
    for (size_t row = 0; row < k; ++row) {
      if (util::abs2(H(row, col)) < cutoff) {
        H(row, col) = 0.0;
      }
    }
  }
}

template <class ptensor>
std::vector<double> Arnoldi<ptensor>::residue(size_t k) const {
  using EigSolver =
      typename std::conditional<std::is_floating_point<value_type>::value,
                                Eigen::EigenSolver<Eigenmatrix>,
                                Eigen::ComplexEigenSolver<Eigenmatrix>>::type;
  auto h = H.topLeftCorner(k, k);
  EigSolver eigensolver(h);
  auto eigvals = eigensolver.eigenvalues();
  auto eigvecs = eigensolver.eigenvectors();

  std::vector<int> sorted_indices(k);
  std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
  std::partial_sort(sorted_indices.begin(), sorted_indices.begin() + nev,
                    sorted_indices.end(), [&eigvals](int a, int b) {
                      return util::abs2(eigvals[a]) > util::abs2(eigvals[b]);
                    });

  std::vector<double> res(nev);
  for (size_t i = 0; i < nev; ++i) {
    int j = sorted_indices[i];
    res[i] = std::abs(H(k, k - 1) * eigvecs(k - 1, j)) / std::abs(eigvals[j]);
  }
  return res;
}

template class Arnoldi<real_tensor>;
template class Arnoldi<complex_tensor>;

}  // end of namespace tenes
