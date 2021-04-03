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

#include "transfer_matrix.hpp"

#include <algorithm>

#include "../tensor.hpp"
#include "../arnoldi.hpp"

#include "../util/abs.hpp"

namespace tenes {

#define SAVE_PARAM(name, type) params_##type[I_##name] = static_cast<type>(name)
#define LOAD_PARAM(name, type) \
  name = static_cast<decltype(name)>(params_##type[I_##name])

void TransferMatrix_Parameters::Bcast(MPI_Comm comm, int root) {
  enum PARAMS_INT_INDEX {
    I_to_calculate,
    I_num_eigvals,
    I_maxdim_dense_eigensolver,
    I_arnoldi_maxdim,
    I_arnoldi_restartdim,
    I_arnoldi_maxiter,

    N_PARAMS_INT_INDEX,
  };
  enum PARAMS_DOUBLE_INDEX {
    I_arnoldi_rtol,

    N_PARAMS_DOUBLE_INDEX,
  };

  int irank;
  MPI_Comm_rank(comm, &irank);

  std::vector<int> params_int(N_PARAMS_INT_INDEX);
  std::vector<double> params_double(N_PARAMS_DOUBLE_INDEX);

  if (irank == root) {
    SAVE_PARAM(to_calculate, int);
    SAVE_PARAM(num_eigvals, int);
    SAVE_PARAM(maxdim_dense_eigensolver, int);
    SAVE_PARAM(arnoldi_maxdim, int);
    SAVE_PARAM(arnoldi_restartdim, int);
    SAVE_PARAM(arnoldi_maxiter, int);
    SAVE_PARAM(arnoldi_rtol, double);

    bcast(params_int, 0, comm);
    bcast(params_double, 0, comm);
  } else {
    bcast(params_int, 0, comm);
    bcast(params_double, 0, comm);

    LOAD_PARAM(to_calculate, int);
    LOAD_PARAM(num_eigvals, int);
    LOAD_PARAM(maxdim_dense_eigensolver, int);
    LOAD_PARAM(arnoldi_maxdim, int);
    LOAD_PARAM(arnoldi_restartdim, int);
    LOAD_PARAM(arnoldi_maxiter, int);
    LOAD_PARAM(arnoldi_rtol, double);
  }
}

#undef SAVE_PARAM
#undef LOAD_PARAM

template <class ptensor>
std::vector<std::complex<double>> TransferMatrix<ptensor>::eigenvalues(
    int dir, int fixed_coord, TransferMatrix_Parameters const &params,
    std::mt19937 &rng) const {
  using value_type = typename ptensor::value_type;
  const size_t N = dim(dir, fixed_coord);
  const size_t nev = std::min(static_cast<size_t>(params.num_eigvals), N);

  std::vector<std::complex<double>> eigvals;

  if (N == 1) {
    eigvals.resize(nev);
    eigvals[0] = 1.0;
    return eigvals;
  }

  if (N <= params.maxdim_dense_eigensolver) {
    ptensor matrix = dir == 0 ? matrix_horizontal(fixed_coord)
                              : matrix_vertical(fixed_coord);

    small_tensor<value_type> matrix_2{mptensor::Shape(N, N)};
    for (size_t row = 0; row < N; ++row) {
      for (size_t col = 0; col < N; ++col) {
        typename ptensor::value_type v;
        matrix.get_value({row, col}, v);
        matrix_2.set_value({row, col}, v);
      }
    }
    std::vector<std::complex<double>> evecs;
    eigen(matrix_2, eigvals, evecs, nev);
  } else {  // use Arnoldi
    auto maxvec = params.arnoldi_maxdim;
    auto maxiter = params.arnoldi_maxiter;
    if (N < maxvec) {
      maxvec = N;
      maxiter = 1;
    }

    ptensor initial_vec = initial_vector(dir, fixed_coord, rng);
    Arnoldi<ptensor> arnoldi(N, maxvec);
    arnoldi.initialize(initial_vec);
    if (dir == 0) {
      arnoldi.run(
          [&](ptensor &out, ptensor const &in) {
            matvec_horizontal(out, in, fixed_coord);
          },
          nev, params.arnoldi_restartdim, maxiter, params.arnoldi_rtol);
    } else {
      arnoldi.run(
          [&](ptensor &out, ptensor const &in) {
            matvec_vertical(out, in, fixed_coord);
          },
          nev, params.arnoldi_restartdim, maxiter, params.arnoldi_rtol);
    }
    eigvals = arnoldi.eigenvalues();
  }
  return eigvals;
}

template <class ptensor>
ptensor TransferMatrix_mf<ptensor>::initial_vector(int dir, int fixed_coord,
                                                   std::mt19937 &rng) const {
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  auto D = dir == 0
               ? this->Tn[this->lattice.index_fast(0, fixed_coord)].shape()[0]
               : this->Tn[this->lattice.index_fast(fixed_coord, 0)].shape()[3];
  auto N = D * D;
  mptensor::Shape sp(N);
  ptensor initial_vec = ptensor(sp);
  for (size_t i = 0; i < N; ++i) {
    double v = dist(rng);
    initial_vec.set_value({i}, v);
  }
  return initial_vec;
}

template <class ptensor>
ptensor TransferMatrix_ctm<ptensor>::initial_vector(int dir, int fixed_coord,
                                                    std::mt19937 &rng) const {
  const auto &lattice = this->lattice;
  const size_t N = C1[0].shape()[0];

  ptensor initial_vec;
  if (dir == 0) {
    int site = lattice.index(0, fixed_coord);
    ptensor left_top = C1[site];
    ptensor left_bottom = C4[lattice.top(site)];
    initial_vec = reshape(tensordot(left_top, left_bottom, {0}, {1}), {N * N});
  } else {
    int site = lattice.index(fixed_coord, 0);
    auto left_bottom = C4[site];
    auto right_bottom = C3[lattice.left(site)];
    initial_vec =
        reshape(tensordot(left_bottom, right_bottom, {0}, {1}), {N * N});
  }
  return initial_vec;
}

template <class ptensor>
void TransferMatrix_ctm<ptensor>::matvec_horizontal(ptensor &outvec,
                                                    ptensor const &invec,
                                                    int y) const {
  using mptensor::Axes;
  const auto &lattice = this->lattice;
  const size_t CHI = C1[0].shape()[0];
  outvec = reshape(invec, {CHI, CHI});

  for (int x = 0; x < lattice.LX; ++x) {
    int site = lattice.index(x, y);
    ptensor top = eTt[site];
    ptensor bottom = eTb[lattice.top(site)];

    ////////////////////////////////////////////////////////////
    // transfer_matvec_horizontal_ctm.dat
    ////////////////////////////////////////////////////////////
    // (top*(bottom*vec))
    // cpu_cost= 1.35e+06  memory= 68400
    // final_bond_order (top_right, bottom_right)
    ////////////////////////////////////////////////////////////
    outvec = tensordot(top, tensordot(bottom, outvec, Axes(1), Axes(1)),
                       Axes(0, 2, 3), Axes(3, 1, 2));
  }
  outvec = reshape(outvec, {CHI * CHI});
}

template <class ptensor>
void TransferMatrix_ctm<ptensor>::matvec_vertical(ptensor &outvec,
                                                  ptensor const &invec,
                                                  int x) const {
  using mptensor::Axes;
  const auto &lattice = this->lattice;
  const size_t CHI = C1[0].shape()[0];
  const auto x_orig = x;
  outvec = reshape(invec, {CHI, CHI});

  do {
    for (int y = 0; y < lattice.LY; ++y) {
      const int site = lattice.index(x, y);
      auto left = eTl[site];
      auto right = eTr[lattice.left(site)];

      ////////////////////////////////////////////////////////////
      // transfer_matvec_vertical_ctm.dat
      ////////////////////////////////////////////////////////////
      // (left*(right*vec))
      // cpu_cost= 1.35e+06  memory= 68400
      // final_bond_order (left_top, right_top)
      ////////////////////////////////////////////////////////////
      outvec = tensordot(left, tensordot(right, outvec, Axes(1), Axes(1)),
                         Axes(0, 2, 3), Axes(3, 1, 2));
    }

    x = (x + lattice.skew + lattice.LX) % lattice.LX;
  } while (x != x_orig);
  outvec = reshape(outvec, {CHI * CHI});
}

template <class ptensor>
void TransferMatrix_mf<ptensor>::matvec_horizontal(ptensor &outvec,
                                                   ptensor const &invec,
                                                   int y) const {
  using mptensor::Axes;
  const auto &lattice = this->lattice;
  auto dim = static_cast<size_t>(std::sqrt(invec.shape()[0]));
  outvec = reshape(invec, {dim, dim});

  for (int x = 0; x < lattice.LX; ++x) {
    const int site = lattice.index(x, y);
    ptensor center = this->Tn[site];
    center.multiply_vector(lambda_tensor[site][1], 1);
    center.multiply_vector(lambda_tensor[site][3], 3);

    ////////////////////////////////////////////////////////////
    // transfer_matvec_horizontal_mf.dat
    ////////////////////////////////////////////////////////////
    // (center*(vec*conj(center)))
    // cpu_cost= 400000  memory= 60100
    // final_bond_order (cr, ccr)
    ////////////////////////////////////////////////////////////
    outvec =
        tensordot(center, tensordot(outvec, conj(center), Axes(1), Axes(0)),
                  Axes(0, 1, 3, 4), Axes(0, 1, 3, 4));
  }
  dim = outvec.shape()[0];
  outvec = reshape(outvec, {dim * dim});
}

template <class ptensor>
void TransferMatrix_mf<ptensor>::matvec_vertical(ptensor &outvec,
                                                 ptensor const &invec,
                                                 int x) const {
  using mptensor::Axes;
  const auto &lattice = this->lattice;
  const auto x_orig = x;
  auto dim = static_cast<size_t>(std::sqrt(invec.shape()[0]));
  outvec = reshape(invec, {dim, dim});

  do {
    for (int y = 0; y < lattice.LY; ++y) {
      const int site = lattice.index(x, y);
      auto center = this->Tn[site];
      center.multiply_vector(lambda_tensor[site][0], 0);
      center.multiply_vector(lambda_tensor[site][2], 2);

      ////////////////////////////////////////////////////////////
      // transfer_matvec_vertical_mf.dat
      ////////////////////////////////////////////////////////////
      // (center*(vec*conj(center)))
      // cpu_cost= 400000  memory= 60100
      // final_bond_order (ct, cct)
      ////////////////////////////////////////////////////////////
      outvec =
          tensordot(center, tensordot(outvec, conj(center), Axes(1), Axes(3)),
                    Axes(0, 2, 3, 4), Axes(1, 3, 0, 4));
    }

    x -= lattice.skew;
    x %= lattice.LX;
    if (x < 0) {
      x += lattice.LX;
    }
  } while (x != x_orig);
  dim = outvec.shape()[0];
  outvec = reshape(outvec, {dim * dim});
}

template <class ptensor>
ptensor TransferMatrix_ctm<ptensor>::matrix_horizontal(int y) const {
  using mptensor::Axes;
  const auto &lattice = this->lattice;
  const size_t CHI = C1[0].shape()[0];
  int site = lattice.index(0, y);
  ptensor top = eTt[site];
  ptensor bottom = eTb[lattice.top(site)];
  ptensor res = transpose(tensordot(top, bottom, Axes(2, 3), Axes(2, 3)),
                          Axes(0, 3, 1, 2));
  for (int x = 1; x < lattice.LX; ++x) {
    site = lattice.index(x, y);
    top = eTt[site];
    bottom = eTb[lattice.top(site)];
    res = tensordot(res, tensordot(top, bottom, Axes(2, 3), Axes(2, 3)),
                    Axes(2, 3), Axes(0, 3));
  }
  res = reshape(res, {CHI * CHI, CHI * CHI});
  return res;
}

template <class ptensor>
ptensor TransferMatrix_ctm<ptensor>::matrix_vertical(int x) const {
  using mptensor::Axes;
  using mptensor::Shape;
  const auto &lattice = this->lattice;
  const size_t CHI = C1[0].shape()[0];
  const auto x_orig = x;
  ptensor res{Shape(CHI * CHI, CHI * CHI)};
#pragma omp parallel for shared(res)
  for (size_t i = 0; i < CHI * CHI; ++i) {
    typename ptensor::value_type v = 1.0;
    res.set_value({i, i}, v);
  }
  res = reshape(res, {CHI, CHI, CHI, CHI});

  do {
    for (int y = 0; y < lattice.LY; ++y) {
      const int site = lattice.index(x, y);
      auto left = eTl[site];
      auto right = eTr[lattice.left(site)];
      res = tensordot(res, tensordot(left, right, Axes(2, 3), Axes(2, 3)),
                      Axes(2, 3), Axes(0, 3));
    }
    x = (x + lattice.skew + lattice.LX) % lattice.LX;
  } while (x != x_orig);
  res = reshape(res, {CHI * CHI, CHI * CHI});
  return res;
}

template <class ptensor>
ptensor TransferMatrix_mf<ptensor>::matrix_horizontal(int y) const {
  using mptensor::Axes;
  const auto &lattice = this->lattice;
  int site = lattice.index(0, y);
  ptensor T = this->Tn[site];
  T.multiply_vector(lambda_tensor[site][1], 1);
  T.multiply_vector(lambda_tensor[site][3], 3);
  size_t dim = T.shape()[0];
  ptensor res = transpose(tensordot(T, conj(T), Axes(1, 3, 4), Axes(1, 3, 4)),
                          Axes(0, 2, 1, 3));

  for (int x = 1; x < lattice.LX; ++x) {
    site = lattice.index(x, y);
    T = this->Tn[site];
    T.multiply_vector(lambda_tensor[site][1], 1);
    T.multiply_vector(lambda_tensor[site][3], 3);
    res = tensordot(res, tensordot(T, conj(T), Axes(1, 3, 4), Axes(1, 3, 4)),
                    Axes(2, 3), Axes(0, 2));
  }
  res = reshape(res, {dim * dim, dim * dim});
  return res;
}

template <class ptensor>
ptensor TransferMatrix_mf<ptensor>::matrix_vertical(int x) const {
  using mptensor::Axes;
  using mptensor::Shape;
  const auto &lattice = this->lattice;
  const size_t N = dim(1, x);
  const size_t Nsqrt = static_cast<size_t>(std::sqrt(N));
  const auto x_orig = x;
  ptensor res{Shape(N, N)};
  for (size_t i = 0; i < N; ++i) {
    typename ptensor::value_type v = 1.0;
    res.set_value({i, i}, v);
  }
  res = reshape(res, {Nsqrt, Nsqrt, Nsqrt, Nsqrt});

  do {
    for (int y = 0; y < lattice.LY; ++y) {
      const int site = lattice.index(x, y);
      auto center = this->Tn[site];
      center.multiply_vector(lambda_tensor[site][0], 0);
      center.multiply_vector(lambda_tensor[site][2], 2);
      res = tensordot(
          res, tensordot(center, conj(center), Axes(0, 2, 4), Axes(0, 2, 4)),
          Axes(2, 3), Axes(1, 3));
    }
    x = (x + lattice.skew + lattice.LX) % lattice.LX;
  } while (x != x_orig);
  res = reshape(res, {N, N});
  return res;
}

template <class ptensor>
size_t TransferMatrix_ctm<ptensor>::dim(int dir, int fixed_coord) const {
  auto ret = C1[0].shape()[0];
  return ret * ret;
}

template <class ptensor>
size_t TransferMatrix_mf<ptensor>::dim(int dir, int fixed_coord) const {
  if (dir == 0) {
    auto ret = this->Tn[this->lattice.index(0, fixed_coord)].shape()[0];
    return ret * ret;
  } else {
    auto ret = this->Tn[this->lattice.index(fixed_coord, 0)].shape()[3];
    return ret * ret;
  }
}

template class TransferMatrix<real_tensor>;
template class TransferMatrix<complex_tensor>;

template class TransferMatrix_ctm<real_tensor>;
template class TransferMatrix_ctm<complex_tensor>;

template class TransferMatrix_mf<real_tensor>;
template class TransferMatrix_mf<complex_tensor>;

}  // end of namespace tenes
