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

#include "arpack_solver.hpp"

#include <algorithm>
#include <complex>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "exception.hpp"
#include "mpi.hpp"
#include "util/abs.hpp"

namespace tenes {

#ifdef TENES_USE_ARPACK

// Direct declarations of the four ARPACK-NG routines (arpack++ style).
// The trailing std::size_t arguments are the hidden lengths of the
// CHARACTER dummy arguments appended by Fortran compilers (gfortran etc.).
// All INTEGER arguments (n, nev, ncv, iparam, ipntr, info, ...) are assumed
// to be 32-bit (standard/"LP64" ARPACK-NG, not ILP64), consistent with the
// 32-bit LAPACK/BLAS integers TeNeS already assumes elsewhere.
extern "C" {
void dnaupd_(int *ido, const char *bmat, const int *n, const char *which,
             const int *nev, const double *tol, double *resid, const int *ncv,
             double *v, const int *ldv, int *iparam, int *ipntr, double *workd,
             double *workl, const int *lworkl, int *info, std::size_t bmat_len,
             std::size_t which_len);
void dneupd_(const int *rvec, const char *howmny, const int *select, double *dr,
             double *di, double *z, const int *ldz, const double *sigmar,
             const double *sigmai, double *workev, const char *bmat,
             const int *n, const char *which, const int *nev, const double *tol,
             double *resid, const int *ncv, double *v, const int *ldv,
             int *iparam, int *ipntr, double *workd, double *workl,
             const int *lworkl, int *info, std::size_t howmny_len,
             std::size_t bmat_len, std::size_t which_len);
void znaupd_(int *ido, const char *bmat, const int *n, const char *which,
             const int *nev, const double *tol, std::complex<double> *resid,
             const int *ncv, std::complex<double> *v, const int *ldv,
             int *iparam, int *ipntr, std::complex<double> *workd,
             std::complex<double> *workl, const int *lworkl, double *rwork,
             int *info, std::size_t bmat_len, std::size_t which_len);
void zneupd_(const int *rvec, const char *howmny, const int *select,
             std::complex<double> *d, std::complex<double> *z, const int *ldz,
             const std::complex<double> *sigma, std::complex<double> *workev,
             const char *bmat, const int *n, const char *which, const int *nev,
             const double *tol, std::complex<double> *resid, const int *ncv,
             std::complex<double> *v, const int *ldv, int *iparam, int *ipntr,
             std::complex<double> *workd, std::complex<double> *workl,
             const int *lworkl, double *rwork, int *info,
             std::size_t howmny_len, std::size_t bmat_len,
             std::size_t which_len);
}

namespace {

using dcomplex = std::complex<double>;

[[noreturn]] void throw_arpack_error(const char *routine, int info) {
  std::stringstream ss;
  ss << "ARPACK routine " << routine << " failed with info = " << info
     << " (see the ARPACK-NG documentation of " << routine << ")";
  if (info == 3 || info == -3) {
    ss << "; hint: increase correlation_length.arnoldi_maxdim";
  }
  throw tenes::runtime_error(ss.str());
}

// gather a distributed rank-1 tensor into a full local copy on every rank
template <class ptensor>
std::vector<typename ptensor::value_type> gather_vector(ptensor const &t) {
  using value_type = typename ptensor::value_type;
  const std::size_t N = t.shape()[0];
  std::vector<value_type> buf(N, 0.0);
  for (std::size_t i = 0; i < t.local_size(); ++i) {
    const auto index = t.global_index(i);
    value_type v;
    if (t.get_value(index, v)) {
      buf[index[0]] = v;
    }
  }
  // every element has exactly one owner, so the sum only adds zeros and
  // the result is bitwise identical on all ranks
  allreduce_sum(buf, t.get_comm());
  return buf;
}

// scatter a full local copy into the distributed rank-1 tensor
template <class ptensor>
void scatter_vector(std::vector<typename ptensor::value_type> const &buf,
                    ptensor &t) {
  for (std::size_t i = 0; i < t.local_size(); ++i) {
    const auto index = t.global_index(i);
    t.set_value(index, buf[index[0]]);
  }
}

using serial_matvec_real =
    std::function<void(std::vector<double> &, std::vector<double> const &)>;
using serial_matvec_complex =
    std::function<void(std::vector<dcomplex> &, std::vector<dcomplex> const &)>;

std::vector<dcomplex> run_arpack(serial_matvec_real const &av,
                                 std::vector<double> &resid, int nev, int ncv,
                                 int maxiter, double tol, bool print_warn) {
  const int n = static_cast<int>(resid.size());
  int ido = 0;
  int info = 1;  // resid contains the initial vector
  int iparam[11] = {};
  int ipntr[14] = {};
  iparam[0] = 1;        // exact shifts
  iparam[2] = maxiter;  // max implicit restarts
  iparam[6] = 1;        // mode 1: standard eigenvalue problem
  const int lworkl = 3 * ncv * ncv + 6 * ncv;
  std::vector<double> v(static_cast<std::size_t>(n) * ncv), workd(3 * n),
      workl(lworkl);
  std::vector<double> x(n), y(n);

  while (true) {
    dnaupd_(&ido, "I", &n, "LM", &nev, &tol, resid.data(), &ncv, v.data(), &n,
            iparam, ipntr, workd.data(), workl.data(), &lworkl, &info, 1, 2);
    if (ido != 1 && ido != -1) {
      break;
    }
    std::copy_n(workd.data() + ipntr[0] - 1, n, x.begin());
    av(y, x);
    std::copy_n(y.data(), n, workd.data() + ipntr[1] - 1);
  }
  if (info < 0 || info == 3) {
    throw_arpack_error("dnaupd", info);
  }
  if (info == 1 && print_warn) {
    std::cerr << "WARNING: ARPACK (dnaupd) reached the maximum number of "
                 "restarts ("
              << maxiter << ") before all eigenvalues converged" << std::endl;
  }

  const int rvec = 0;  // eigenvalues only
  std::vector<int> select(ncv);
  std::vector<double> dr(nev + 1), di(nev + 1), workev(3 * ncv);
  const double sigmar = 0.0, sigmai = 0.0;
  int ierr = 0;
  dneupd_(&rvec, "A", select.data(), dr.data(), di.data(), v.data(), &n,
          &sigmar, &sigmai, workev.data(), "I", &n, "LM", &nev, &tol,
          resid.data(), &ncv, v.data(), &n, iparam, ipntr, workd.data(),
          workl.data(), &lworkl, &ierr, 1, 1, 2);
  if (ierr != 0) {
    throw_arpack_error("dneupd", ierr);
  }

  // ARPACK does not return dr/di sorted by decreasing |lambda|, and when
  // the nev-th and (nev+1)-th Ritz values form a complex-conjugate pair,
  // dngets bumps the internal nev up by one to keep the pair together, so
  // iparam[4] ("nconv") can be nev + 1. Read every value ARPACK actually
  // wrote (dr/di are sized nev + 1) rather than blindly slicing the first
  // nev, or the true largest-magnitude eigenvalue can be dropped; the
  // caller sorts by |lambda| and truncates to nev afterwards.
  const int nconv = std::min(iparam[4], nev + 1);
  std::vector<dcomplex> ev;
  ev.reserve(nconv);
  for (int i = 0; i < nconv; ++i) {
    ev.emplace_back(dr[i], di[i]);
  }
  return ev;
}

std::vector<dcomplex> run_arpack(serial_matvec_complex const &av,
                                 std::vector<dcomplex> &resid, int nev, int ncv,
                                 int maxiter, double tol, bool print_warn) {
  const int n = static_cast<int>(resid.size());
  int ido = 0;
  int info = 1;  // resid contains the initial vector
  int iparam[11] = {};
  int ipntr[14] = {};
  iparam[0] = 1;        // exact shifts
  iparam[2] = maxiter;  // max implicit restarts
  iparam[6] = 1;        // mode 1: standard eigenvalue problem
  const int lworkl = 3 * ncv * ncv + 5 * ncv;
  std::vector<dcomplex> v(static_cast<std::size_t>(n) * ncv), workd(3 * n),
      workl(lworkl);
  std::vector<double> rwork(ncv);
  std::vector<dcomplex> x(n), y(n);

  while (true) {
    znaupd_(&ido, "I", &n, "LM", &nev, &tol, resid.data(), &ncv, v.data(), &n,
            iparam, ipntr, workd.data(), workl.data(), &lworkl, rwork.data(),
            &info, 1, 2);
    if (ido != 1 && ido != -1) {
      break;
    }
    std::copy_n(workd.data() + ipntr[0] - 1, n, x.begin());
    av(y, x);
    std::copy_n(y.data(), n, workd.data() + ipntr[1] - 1);
  }
  if (info < 0 || info == 3) {
    throw_arpack_error("znaupd", info);
  }
  if (info == 1 && print_warn) {
    std::cerr << "WARNING: ARPACK (znaupd) reached the maximum number of "
                 "restarts ("
              << maxiter << ") before all eigenvalues converged" << std::endl;
  }

  const int rvec = 0;  // eigenvalues only
  std::vector<int> select(ncv);
  std::vector<dcomplex> d(nev + 1), workev(2 * ncv);
  const dcomplex sigma(0.0, 0.0);
  int ierr = 0;
  zneupd_(&rvec, "A", select.data(), d.data(), v.data(), &n, &sigma,
          workev.data(), "I", &n, "LM", &nev, &tol, resid.data(), &ncv,
          v.data(), &n, iparam, ipntr, workd.data(), workl.data(), &lworkl,
          rwork.data(), &ierr, 1, 1, 2);
  if (ierr != 0) {
    throw_arpack_error("zneupd", ierr);
  }

  // See the real overload above: nconv (iparam[4]) can be nev + 1 when a
  // conjugate pair straddles the requested nev, and d is sized nev + 1 to
  // allow for that.
  const int nconv = std::min(iparam[4], nev + 1);
  return std::vector<dcomplex>(d.begin(), d.begin() + nconv);
}

}  // namespace

template <class ptensor>
std::vector<std::complex<double>> arpack_eigenvalues(
    std::function<void(ptensor &, ptensor const &)> A, ptensor const &initial,
    std::size_t nev, int ncv, int maxiter, double tol) {
  using value_type = typename ptensor::value_type;
  const int N = static_cast<int>(initial.shape()[0]);
  const int nev_ = static_cast<int>(nev);
  const bool print_warn = initial.get_comm_rank() == 0;

  int ncv_ = std::min(ncv, N);
  if (ncv_ < nev_ + 2) {
    ncv_ = std::min(nev_ + 2, N);
  }
  if (nev_ >= ncv_ - 1) {
    std::stringstream ss;
    ss << "ARPACK requires num_eigvals + 2 <= min(arnoldi_maxdim, N) but "
          "num_eigvals = "
       << nev_ << ", N = " << N
       << "; hint: increase correlation_length.maxdim_dense_eigensolver "
          "above N to use the dense eigensolver instead";
    throw tenes::input_error(ss.str());
  }

  const MPI_Comm comm = initial.get_comm();
  std::vector<value_type> resid = gather_vector(initial);
  auto serial_matvec = [&](std::vector<value_type> &out,
                           std::vector<value_type> const &in) {
    ptensor x(comm, mptensor::Shape(N));
    scatter_vector(in, x);
    ptensor y;
    A(y, x);
    out = gather_vector(y);
  };

  std::vector<dcomplex> ev =
      run_arpack(serial_matvec, resid, nev_, ncv_, maxiter, tol, print_warn);

  std::sort(ev.begin(), ev.end(), [](dcomplex const &a, dcomplex const &b) {
    return util::abs2(a) > util::abs2(b);
  });
  if (ev.size() < nev) {
    if (print_warn) {
      std::cerr << "WARNING: ARPACK converged only " << ev.size() << " of "
                << nev << " eigenvalues; the rest are reported as NaN"
                << std::endl;
    }
  }
  // Unconditional: run_arpack() may return up to nev + 1 entries (see the
  // comment at its nconv computation), so this both pads (fewer than nev
  // converged) and truncates (a conjugate pair pushed nconv to nev + 1) to
  // exactly nev, keeping the nev largest-|lambda| values now that ev is
  // sorted.
  ev.resize(nev, dcomplex(std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN()));
  return ev;
}

#else  // TENES_USE_ARPACK

template <class ptensor>
std::vector<std::complex<double>> arpack_eigenvalues(
    std::function<void(ptensor &, ptensor const &)> /* A */,
    ptensor const & /* initial */, std::size_t /* nev */, int /* ncv */,
    int /* maxiter */, double /* tol */) {
  throw std::logic_error(
      "internal error: arpack_eigenvalues is called, but TeNeS was built "
      "without ARPACK-NG");
}

#endif  // TENES_USE_ARPACK

template std::vector<std::complex<double>> arpack_eigenvalues<real_tensor>(
    std::function<void(real_tensor &, real_tensor const &)>,
    real_tensor const &, std::size_t, int, int, double);
template std::vector<std::complex<double>> arpack_eigenvalues<complex_tensor>(
    std::function<void(complex_tensor &, complex_tensor const &)>,
    complex_tensor const &, std::size_t, int, int, double);

}  // end of namespace tenes
