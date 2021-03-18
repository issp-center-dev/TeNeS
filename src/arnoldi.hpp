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

#include <vector>
#include <functional>

namespace tenes {

template <class ptensor>
class Arnoldi {
 public:
  using value_type = typename ptensor::value_type;
  using Eigenmatrix = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>;

  Arnoldi(size_t N, size_t maxvec);
  void initialize(ptensor const &initial);
  void run(std::function<void(ptensor &, ptensor const &)> A, size_t nev,
           int mindim = 20, int maxiter = 10, double rtol = 1.0e-8);
  std::vector<std::complex<double>> eigenvalues() const;

 private:
  void orthonormalize(size_t k);
  void restart(size_t minvec, double cutoff = 1.0e-12);
  std::vector<double> residue(size_t k) const;

  size_t N;
  size_t maxvec;
  size_t nev;

  std::vector<ptensor> Q;
  Eigenmatrix H;
};

}  // end of namespace tenes

#endif  // TENES_SRC_ARNOLDI_HPP_
