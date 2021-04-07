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

#include <algorithm>   // for min
#include <cmath>       // for sqrt
#include <complex>     // for complex, ope...
#include <functional>  // for function
#include <vector>      // for __vector_bas...

#include "util/abs.hpp"
#include "tensor.hpp"

using mptensor::Shape;
using dcomplex = std::complex<double>;

namespace tenes {

template <class ptensor>
Arnoldi<ptensor>::Arnoldi(size_t N, size_t maxvec)
    : N(N), maxvec(maxvec), Q(maxvec + 1), H(Shape(maxvec + 1, maxvec)) {
  for (size_t i = 0; i < maxvec; ++i) {
    Q[i] = ptensor(mptensor::Shape(N));
  }
}

template <class ptensor>
void Arnoldi<ptensor>::initialize(ptensor const& initial) {
  Q[0] = initial;
  orthonormalize(0);
}

template <class ptensor>
double norm(ptensor const& A) {
  auto x = tensordot(conj(A), A, {0}, {0});
  return std::sqrt(std::real(trace(x, {}, {})));
}

template <class ptensor>
void Arnoldi<ptensor>::run(std::function<void(ptensor&, ptensor const&)> A,
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
  auto h = slice(H, 0, 0, maxvec);
  std::vector<dcomplex> evs(maxvec);
  std::vector<dcomplex> ret(maxvec);
  eigen(h, ret, evs, nev);
  return ret;
}

template <class ptensor>
void Arnoldi<ptensor>::orthonormalize(size_t k) {
  for (size_t j = 0; j < k; ++j) {
    auto xx = tensordot(conj(Q[j]), Q[k], {0}, {0});
    // std::cout << xx << std::endl;
    auto x = trace(xx, {}, {});
    // auto x = trace(tensordot(conj(Q[j]), Q[k], {0}, {0}), {}, {});
    H.set_value({j, k - 1}, x);
    Q[k] -= Q[j] * x;
  }
  auto x = norm(Q[k]);
  Q[k] *= 1.0 / x;
  if (k > 0) {
    H.set_value({k, k - 1}, x);
  }
}

double real_if_real(dcomplex v, double) { return std::real(v); }
dcomplex real_if_real(dcomplex v, dcomplex) { return v; }

template <class ptensor>
void Arnoldi<ptensor>::restart(size_t minvec, double cutoff) {
  size_t k = maxvec;
  auto h = slice(H, 0, 0, maxvec);
  std::vector<dcomplex> eigvals, vecs_last;
  eigen(h, eigvals, vecs_last, maxvec);

  for (size_t i = 0; i < k - minvec; ++i) {
    value_type sigma = real_if_real(eigvals[i], value_type());
    small_tensor<value_type> A = h - identity(k, sigma);
    small_tensor<value_type> q, r;
    qr(A, q, r);

    h = tensordot(tensordot(conj(q), h, {0}, {0}), q, {1}, {0});

    size_t n = 0;
    constexpr size_t max_block_size = 10;
    while (n < N) {
      size_t block_size = std::min(max_block_size, N - n);
      small_tensor<value_type> newq{Shape(block_size, k)};
      for (size_t j = 0; j < k; ++j) {
        for (size_t b = 0; b < block_size; ++b) {
          value_type temp;
          Q[j].get_value({n + b}, temp);
          newq.set_value({b, j}, temp);
        }
      }
      newq = tensordot(newq, q, {1}, {0});
      for (size_t j = 0; j < k; ++j) {
        for (size_t b = 0; b < block_size; ++b) {
          value_type temp;
          newq.get_value({b, j}, temp);
          Q[j].set_value({n + b}, temp);
        }
      }
      n += block_size;
    }
  }
  H.set_slice(h, 0, 0, maxvec);

  cutoff *= cutoff;
  for (size_t col = 0; col < k; ++col) {
    for (size_t row = 0; row < k; ++row) {
      value_type temp;
      H.get_value({row, col}, temp);
      if (util::abs2(temp) < cutoff) {
        H.set_value({row, col}, 0.0);
      }
    }
  }
}

template <class ptensor>
std::vector<double> Arnoldi<ptensor>::residue(size_t k) const {
  auto h = slice(H, 0, 0, k);
  std::vector<dcomplex> eigvals, eigvecs;
  eigen(h, eigvals, eigvecs, nev);
  std::vector<double> res(nev);
  value_type hk;
  H.get_value({k, k - 1}, hk);
  for (size_t i = 0; i < nev; ++i) {
    res[i] = std::abs(hk * eigvecs[i]) / std::abs(eigvals[i]);
  }
  return res;
}

template class Arnoldi<real_tensor>;
template class Arnoldi<complex_tensor>;

}  // end of namespace tenes
