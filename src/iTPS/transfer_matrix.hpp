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

#ifndef TENES_SRC_ITPS_TRANSFER_MATRIX_HPP_
#define TENES_SRC_ITPS_TRANSFER_MATRIX_HPP_

#include <vector>
#include <random>

#include "../tensor.hpp"
#include "../SquareLattice.hpp"

namespace tenes {

struct TransferMatrix_Parameters {
  bool to_calculate;
  int num_eigvals;
  int maxdim_dense_eigensolver;
  int arnoldi_maxdim;
  int arnoldi_restartdim;
  int arnoldi_maxiter;
  double arnoldi_rtol;

  TransferMatrix_Parameters()
      : to_calculate(true),
        num_eigvals(4),
        maxdim_dense_eigensolver(200),
        arnoldi_maxdim(50),
        arnoldi_restartdim(20),
        arnoldi_maxiter(1),
        arnoldi_rtol(1.0e-10) {}
  void Bcast(MPI_Comm comm, int root = 0);
};

template <class ptensor>
class TransferMatrix {
 public:
  TransferMatrix(SquareLattice const &lattice, std::vector<ptensor> const &Tn)
      : lattice(lattice), Tn(Tn) {}
  virtual ~TransferMatrix() {}
  std::vector<std::complex<double>> eigenvalues(
      int dir, int fixed_coord, TransferMatrix_Parameters const &param,
      std::mt19937 &rng) const;

 private:
  virtual void matvec_horizontal(ptensor &outvec, ptensor const &invec,
                                 int y) const = 0;
  virtual void matvec_vertical(ptensor &outvec, ptensor const &invec,
                               int x) const = 0;

  virtual ptensor matrix_horizontal(int y) const = 0;
  virtual ptensor matrix_vertical(int y) const = 0;

  virtual ptensor initial_vector(int dir, int fixed_coord,
                                 std::mt19937 &rng) const = 0;

  virtual size_t dim(int dir, int fixed_coord) const = 0;

 protected:
  SquareLattice lattice;
  const std::vector<ptensor> &Tn;
};

template <class ptensor>
class TransferMatrix_ctm : public TransferMatrix<ptensor> {
 public:
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
  const std::vector<ptensor> &C1;
  const std::vector<ptensor> &C2;
  const std::vector<ptensor> &C3;
  const std::vector<ptensor> &C4;

  const std::vector<ptensor> &eTl;
  const std::vector<ptensor> &eTt;
  const std::vector<ptensor> &eTr;
  const std::vector<ptensor> &eTb;

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

template <class ptensor>
class TransferMatrix_mf : public TransferMatrix<ptensor> {
 public:
  TransferMatrix_mf(
      const SquareLattice &lattice, const std::vector<ptensor> &Tn,
      const std::vector<std::vector<std::vector<double>>> lambda_tensor)
      : TransferMatrix<ptensor>(lattice, Tn),
        lambda_tensor(lambda_tensor) {}

 private:
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

}  // end of namespace tenes

#endif  // TENES_SRC_ITPS_TRANSFER_MATRIX_HPP_
