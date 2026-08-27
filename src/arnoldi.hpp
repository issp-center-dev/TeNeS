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

/*! @file
 *  @brief The builtin implicit-restart Arnoldi eigensolver
 *         (::tenes::Arnoldi), used for transfer-matrix spectra when
 *         ARPACK-NG is not available or not selected.
 */

#ifndef TENES_SRC_ARNOLDI_HPP_
#define TENES_SRC_ARNOLDI_HPP_

#include <cstddef>
#include <complex>
#include <vector>      // IWYU pragma: export
#include <functional>  // IWYU pragma: export
#include "tensor.hpp"  // IWYU pragma: export

namespace tenes {

/*! @brief Eigenvalue calculator for mptensor::Tensor
 *
 * by using (implicit restart) Arnoldi method
 */
template <class ptensor>
class Arnoldi {
 public:
  using value_type = typename ptensor::value_type;  //!< scalar type

  //! Prepare for a problem of dimension N with a Krylov subspace of at
  //! most maxvec vectors.
  Arnoldi(std::size_t N, std::size_t maxvec);
  //! Set (and normalize) the start vector of the iteration.
  void initialize(ptensor const &initial);

  /*! @brief perform Arnoldi method
   *
   * @param[in] A "Matrix" as a function taking a "vector" `x` and returning
   * another, `Ax`. The first [out] argument of `A` is for `Ax` and the second
   * [in] argument is for `x`.
   * @param[in] nev Number of eigenvalues to be calculated
   * @param[in] mindim Number of vectors re-generated in restart
   * @param[in] maxiter Maximum number of iterations
   * @param[in] rtol Relative torelance for convergence check
   *
   */
  void run(std::function<void(ptensor &, ptensor const &)> A, std::size_t nev,
           int mindim = 20, int maxiter = 10, double rtol = 1.0e-8);

  /*! @brief return eigenvalues as a vector
   */
  std::vector<std::complex<double>> eigenvalues() const;

 private:
  //! Orthonormalize the k-th Krylov vector against the previous ones.
  void orthonormalize(std::size_t k);
  //! Implicit restart: compress the subspace back to minvec vectors.
  void restart(std::size_t minvec, double cutoff = 1.0e-12);
  //! Residual estimate of each of the k current Ritz pairs.
  std::vector<double> residue(std::size_t k) const;

  std::size_t N;       //!< dimension of the problem
  std::size_t maxvec;  //!< maximum Krylov subspace size
  std::size_t nev;     //!< number of requested eigenvalues

  std::vector<ptensor> Q;      //!< Krylov basis vectors
  small_tensor<value_type> H;  //!< projected (Hessenberg) matrix

  int mpisize;  //!< size of the communicator the vectors live on
  int mpirank;  //!< rank within that communicator
};

}  // end of namespace tenes

#endif  // TENES_SRC_ARNOLDI_HPP_
