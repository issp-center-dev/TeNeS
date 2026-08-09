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

#ifndef TENES_SRC_ARPACK_SOLVER_HPP_
#define TENES_SRC_ARPACK_SOLVER_HPP_

#include <complex>
#include <cstddef>
#include <functional>
#include <vector>      // IWYU pragma: export
#include "tensor.hpp"  // IWYU pragma: export

namespace tenes {

//! @return true if TeNeS was built with ARPACK-NG
constexpr bool arpack_available() {
#ifdef TENES_USE_ARPACK
  return true;
#else
  return false;
#endif
}

/*! @brief Largest-magnitude eigenvalues via ARPACK-NG (dnaupd/znaupd)
 *
 * Every rank redundantly runs the same serial ARPACK iteration; only the
 * matrix-vector product runs on the distributed tensor.
 *
 * @param[in] A "Matrix" as a function taking a "vector" `x` and returning
 *              another, `Ax`. The first [out] argument of `A` is for `Ax`
 *              and the second [in] argument is for `x` (same as Arnoldi).
 * @param[in] initial Initial (distributed) vector; also fixes the problem
 *                    size N via its shape.
 * @param[in] nev Number of eigenvalues to be calculated
 * @param[in] ncv Dimension of the Krylov subspace (clamped into
 *                [nev+2, N] internally)
 * @param[in] maxiter Maximum number of implicit restarts
 * @param[in] tol Relative tolerance on the Ritz values (<= 0 means
 *                machine epsilon)
 * @param[in] grow_ncv If true and fewer than nev eigenvalues converge,
 *                     double ncv (capped at N, where convergence is exact)
 *                     and retry until all nev converge or ncv reaches N.
 * @return Eigenvalues sorted by decreasing magnitude, size nev. If fewer
 *         than nev eigenvalues converge, the tail is NaN and a warning is
 *         printed on rank 0.
 *
 * Throws tenes::runtime_error on ARPACK failures and std::logic_error if
 * TeNeS was built without ARPACK-NG.
 */
template <class ptensor>
std::vector<std::complex<double>> arpack_eigenvalues(
    std::function<void(ptensor &, ptensor const &)> A, ptensor const &initial,
    std::size_t nev, int ncv, int maxiter, double tol, bool grow_ncv = false);

}  // end of namespace tenes

#endif  // TENES_SRC_ARPACK_SOLVER_HPP_
