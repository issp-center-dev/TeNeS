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

#include <cstddef>
#include <complex>
#include <vector>      // IWYU pragma: export
#include <functional>  // IWYU pragma: export
#include "tensor.hpp"  // IWYU pragma: export

namespace tenes {

template <class ptensor>
class Arnoldi {
 public:
  using value_type = typename ptensor::value_type;

  Arnoldi(std::size_t N, std::size_t maxvec);
  void initialize(ptensor const &initial);
  void run(std::function<void(ptensor &, ptensor const &)> A, std::size_t nev,
           int mindim = 20, int maxiter = 10, double rtol = 1.0e-8);
  std::vector<std::complex<double>> eigenvalues() const;

 private:
  void orthonormalize(std::size_t k);
  void restart(std::size_t minvec, double cutoff = 1.0e-12);
  std::vector<double> residue(std::size_t k) const;

  std::size_t N;
  std::size_t maxvec;
  std::size_t nev;

  std::vector<ptensor> Q;
  small_tensor<value_type> H;
};

}  // end of namespace tenes

#endif  // TENES_SRC_ARNOLDI_HPP_
