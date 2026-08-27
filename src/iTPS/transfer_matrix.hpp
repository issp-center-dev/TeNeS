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
 *  @brief Transfer matrices of the iTPS, whose leading eigenvalues give
 *         the correlation lengths.
 *
 *  A row (or column) of the contracted double-layer network defines a
 *  transfer matrix; the ratio of its subleading to leading eigenvalue is
 *  @f$e^{-1/\xi}@f$ with @f$\xi@f$ the correlation length along that
 *  direction. ::tenes::itps::TransferMatrix provides the eigensolver
 *  driver; the two subclasses realize the matrix-vector product with the
 *  CTM environment or the mean-field environment.
 */

#ifndef TENES_SRC_ITPS_TRANSFER_MATRIX_HPP_
#define TENES_SRC_ITPS_TRANSFER_MATRIX_HPP_

#include <algorithm>
#include <vector>
#include <random>

#include "../tensor.hpp"
#include "../SquareLattice.hpp"

namespace tenes::itps {

//! Resolve the automatic (non-positive) arnoldi_maxdim.
//!
//! The subspace dimension scales with the number of requested eigenvalues,
//! and the floor depends on the solver: ARPACK-NG converges reliably from a
//! small subspace by restarting, while the builtin solver has a weak restart
//! and needs a subspace large enough to converge in a single sweep.
inline int effective_arnoldi_maxdim(int arnoldi_maxdim, int num_eigvals,
                                    bool use_arpack) {
  if (arnoldi_maxdim > 0) {
    return arnoldi_maxdim;
  }
  return std::max(2 * num_eigvals + 1, use_arpack ? 25 : 50);
}

//! Resolve the automatic (non-positive) arnoldi_maxiterations: ARPACK-NG
//! is given room to restart; the builtin solver is not (see above).
inline int effective_arnoldi_maxiter(int arnoldi_maxiter, bool use_arpack) {
  if (arnoldi_maxiter > 0) {
    return arnoldi_maxiter;
  }
  return use_arpack ? 10 : 1;
}

//! Resolve the automatic (non-positive) arnoldi_restartdim: keep half of
//! the Krylov subspace on restart (at least num_eigvals + 1 vectors).
inline int effective_arnoldi_restartdim(int arnoldi_restartdim, int num_eigvals,
                                        int maxdim) {
  if (arnoldi_restartdim > 0) {
    return arnoldi_restartdim;
  }
  return std::max(num_eigvals + 1, maxdim / 2);
}

//! How to solve the transfer-matrix eigenproblem (when larger than
//! maxdim_dense_eigensolver)
enum class TransferMatrixEigensolver : int {
  automatic = 0,  //!< ARPACK-NG if built in, otherwise builtin
  arpack = 1,     //!< ARPACK-NG (input error if not built in)
  builtin = 2,    //!< builtin implicit-restart Arnoldi
};

/*! @brief Settings of the correlation-length measurement, read from the
 *         [correlation_length] table of input.toml.
 */
struct TransferMatrix_Parameters {
  //! Whether to measure the correlation length at all (measure).
  bool to_calculate;
  //! Number of eigenvalues to compute per transfer matrix
  //! (num_eigvals).
  int num_eigvals;
  //! Transfer matrices up to this dimension are diagonalized densely;
  //! larger ones go to an Arnoldi solver
  //! (maxdim_dense_eigensolver).
  int maxdim_dense_eigensolver;
  //! Krylov subspace dimension (arnoldi_maxdim); non-positive means
  //! automatic, see effective_arnoldi_maxdim().
  int arnoldi_maxdim;
  //! Vectors kept on restart (arnoldi_restartdim); non-positive means
  //! automatic, see effective_arnoldi_restartdim().
  int arnoldi_restartdim;
  //! Maximum number of Arnoldi restarts (arnoldi_maxiter); non-positive
  //! means automatic, see effective_arnoldi_maxiter().
  int arnoldi_maxiter;
  //! Relative tolerance of the Arnoldi eigensolver (arnoldi_rtol).
  double arnoldi_rtol;
  //! Which Arnoldi implementation to use (solver).
  TransferMatrixEigensolver eigensolver;

  //! Set the documented default of every setting.
  TransferMatrix_Parameters()
      : to_calculate(true),
        num_eigvals(4),
        maxdim_dense_eigensolver(200),
        arnoldi_maxdim(0),      // 0 means automatic
        arnoldi_restartdim(0),  // 0 means automatic
        arnoldi_maxiter(0),     // 0 means automatic
        arnoldi_rtol(1.0e-10),
        eigensolver(TransferMatrixEigensolver::automatic) {}
  //! Broadcast all settings from root to every rank.
  void Bcast(MPI_Comm comm, int root = 0);
};

/*! @brief Abstract transfer matrix of one row or column of the network.
 *
 *  Concrete environments (TransferMatrix_ctm, TransferMatrix_mf) supply
 *  the matrix-vector product; eigenvalues() picks the dense or Arnoldi
 *  eigensolver by the matrix dimension and returns the leading spectrum.
 *
 *  @tparam ptensor ::tenes::real_tensor or ::tenes::complex_tensor.
 */
template <class ptensor>
class TransferMatrix {
 public:
  //! Hold references to the lattice and the center tensors.
  TransferMatrix(SquareLattice const &lattice, std::vector<ptensor> const &Tn)
      : lattice(lattice), Tn(Tn) {}
  virtual ~TransferMatrix() {}
  /*!
   * @brief Leading eigenvalues of one transfer matrix.
   *
   * @param[in] dir 0 for horizontal (a row), 1 for vertical (a column).
   * @param[in] fixed_coord The y (dir = 0) or x (dir = 1) coordinate of
   *            the row/column.
   * @param[in] param Eigensolver settings.
   * @param[in] rng Source of the random initial vector.
   * @return Eigenvalues sorted by descending magnitude.
   */
  std::vector<std::complex<double>> eigenvalues(
      int dir, int fixed_coord, TransferMatrix_Parameters const &param,
      std::mt19937 &rng) const;

 private:
  //! Apply the row transfer matrix at height y to invec.
  virtual void matvec_horizontal(ptensor &outvec, ptensor const &invec,
                                 int y) const = 0;
  //! Apply the column transfer matrix at position x to invec.
  virtual void matvec_vertical(ptensor &outvec, ptensor const &invec,
                               int x) const = 0;

  //! Dense matrix of the row transfer matrix (for small dimensions).
  virtual ptensor matrix_horizontal(int y) const = 0;
  //! Dense matrix of the column transfer matrix.
  virtual ptensor matrix_vertical(int y) const = 0;

  //! Random start vector of the iterative eigensolver.
  virtual ptensor initial_vector(int dir, int fixed_coord,
                                 std::mt19937 &rng) const = 0;

  //! Dimension of the transfer matrix as a linear map.
  virtual size_t dim(int dir, int fixed_coord) const = 0;

 protected:
  SquareLattice lattice;           //!< unit-cell geometry
  const std::vector<ptensor> &Tn;  //!< center tensors, per site
};

/*! @brief Transfer matrix built from the CTM edge tensors: one
 *         application contracts the pair of opposing edge tensors
 *         (already carrying the renormalized double layer) across the
 *         whole row/column, acting on a (chi x chi) boundary vector.
 */
template <class ptensor>
class TransferMatrix_ctm : public TransferMatrix<ptensor> {
 public:
  //! Hold references to the state and its full CTM environment.
  TransferMatrix_ctm(
      const SquareLattice &lattice, const std::vector<ptensor> &Tn,
      const std::vector<ptensor> &C1, const std::vector<ptensor> &C2,
      const std::vector<ptensor> &C3, const std::vector<ptensor> &C4,
      const std::vector<ptensor> &eTl, const std::vector<ptensor> &eTt,
      const std::vector<ptensor> &eTr, const std::vector<ptensor> &eTb)
      : TransferMatrix<ptensor>(lattice, Tn),
        C1(C1),
        C2(C2),
        C3(C3),
        C4(C4),
        eTl(eTl),
        eTt(eTt),
        eTr(eTr),
        eTb(eTb) {}

 private:
  const std::vector<ptensor> &C1;  //!< left-top corners, per site
  const std::vector<ptensor> &C2;  //!< right-top corners, per site
  const std::vector<ptensor> &C3;  //!< right-bottom corners, per site
  const std::vector<ptensor> &C4;  //!< left-bottom corners, per site

  const std::vector<ptensor> &eTl;  //!< left edge tensors, per site
  const std::vector<ptensor> &eTt;  //!< top edge tensors, per site
  const std::vector<ptensor> &eTr;  //!< right edge tensors, per site
  const std::vector<ptensor> &eTb;  //!< bottom edge tensors, per site

  void matvec_horizontal(ptensor &outvec, ptensor const &invec,
                         int y) const override;
  void matvec_vertical(ptensor &outvec, ptensor const &invec,
                       int x) const override;
  ptensor matrix_horizontal(int y) const override;
  ptensor matrix_vertical(int x) const override;
  ptensor initial_vector(int dir, int fixed_coord,
                         std::mt19937 &gen) const override;
  size_t dim(int dir, int fixed_coord) const override;
};

/*! @brief Transfer matrix with the mean-field environment: the
 *         lambda-dressed center tensors alone form the row/column and the
 *         vector lives on the (D x D) double-layer bond.
 */
template <class ptensor>
class TransferMatrix_mf : public TransferMatrix<ptensor> {
 public:
  //! Hold references to the state and its mean-field weights.
  TransferMatrix_mf(
      const SquareLattice &lattice, const std::vector<ptensor> &Tn,
      const std::vector<std::vector<std::vector<double>>> lambda_tensor)
      : TransferMatrix<ptensor>(lattice, Tn), lambda_tensor(lambda_tensor) {}

 private:
  //! Mean-field weights, indexed by [site][leg].
  const std::vector<std::vector<std::vector<double>>> lambda_tensor;

  void matvec_horizontal(ptensor &outvec, ptensor const &invec,
                         int y) const override;
  void matvec_vertical(ptensor &outvec, ptensor const &invec,
                       int x) const override;
  ptensor matrix_horizontal(int y) const override;
  ptensor matrix_vertical(int x) const override;
  ptensor initial_vector(int dir, int fixed_coord,
                         std::mt19937 &rng) const override;
  size_t dim(int dir, int fixed_coord) const override;
};

}  // namespace tenes::itps

#endif  // TENES_SRC_ITPS_TRANSFER_MATRIX_HPP_
