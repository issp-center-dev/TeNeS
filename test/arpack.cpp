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

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "../src/arpack_solver.hpp"
#include "../src/arnoldi.hpp"
#include "../src/tensor.hpp"
#include "../src/mpi.hpp"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  doctest::Context context(argc, argv);
  const int res = context.run();
  MPI_Finalize();
  return res;
}

namespace {

using rtensor = tenes::real_tensor;
using ctensor = tenes::complex_tensor;

// The matvec functions below are written to be safe to run with more than
// one rank (even though this build/test only ever exercises a single rank:
// ENABLE_MPI is off on this platform). They only touch entries of `out`
// that `out.local_size()`/`out.global_index()` say this rank owns, and
// initialize scratch values before a possibly-failing `get_value()`
// (get_value returns false, leaving its output argument untouched, for an
// index this rank does not own) instead of reading uninitialized memory.

// A x for the diagonal matrix A = diag(2^0, 2^-1, ..., 2^-(N-1)), scaled by
// 16. Diagonal, so `out`'s and `in`'s local entries line up index-for-index
// (both share the same shape and comm, hence the same distribution).
void matvec_real(rtensor &out, rtensor const &in, std::size_t N) {
  out = rtensor(in.get_comm(), mptensor::Shape(N));
  for (std::size_t n = 0; n < out.local_size(); ++n) {
    const auto index = out.global_index(n);
    double v = 0.0;
    in.get_value(index, v);
    out.set_value(index,
                  16.0 * std::pow(2.0, -static_cast<double>(index[0])) * v);
  }
}

// the same spectrum rotated to the imaginary axis: A = diag(16i * 2^-i)
void matvec_complex(ctensor &out, ctensor const &in, std::size_t N) {
  out = ctensor(in.get_comm(), mptensor::Shape(N));
  const std::complex<double> I(0.0, 1.0);
  for (std::size_t n = 0; n < out.local_size(); ++n) {
    const auto index = out.global_index(n);
    std::complex<double> v(0.0, 0.0);
    in.get_value(index, v);
    out.set_value(index,
                  16.0 * I * std::pow(2.0, -static_cast<double>(index[0])) * v);
  }
}

// gather a distributed rank-1 real tensor into a full copy on every rank.
// Needed (unlike the diagonal matvecs above) because the conjugate-pair
// matvec below couples entries 1 and 2, which is not guaranteed to be a
// purely-local operation once run with more than one rank.
std::vector<double> gather_real(rtensor const &t, std::size_t N) {
  std::vector<double> buf(N, 0.0);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const auto index = t.global_index(n);
    double v = 0.0;
    if (t.get_value(index, v)) {
      buf[index[0]] = v;
    }
  }
  tenes::allreduce_sum(buf, t.get_comm());
  return buf;
}

// A with an isolated eigenvalue 10 and a complex-conjugate pair 8 +/- 1i
// coming from the 2x2 block [[8, 1], [-1, 8]] (char. poly (8-l)^2 + 1 = 0),
// with all remaining eigenvalues much smaller in magnitude. Regression case
// for ARPACK converging nconv = nev + 1: requesting nev = 2 largest-|l|
// eigenvalues would otherwise split the conjugate pair, so ARPACK reports
// all three of {10, 8+1i, 8-1i} and a naive "keep the first nev" slice
// (before sorting by magnitude) can drop the dominant eigenvalue 10.
void matvec_conjugate_pair(rtensor &out, rtensor const &in, std::size_t N) {
  out = rtensor(in.get_comm(), mptensor::Shape(N));
  std::vector<double> x = gather_real(in, N);
  std::vector<double> y(N, 0.0);
  y[0] = 10.0 * x[0];
  y[1] = 8.0 * x[1] + 1.0 * x[2];
  y[2] = -1.0 * x[1] + 8.0 * x[2];
  for (std::size_t i = 3; i < N; ++i) {
    y[i] = std::pow(2.0, -static_cast<double>(i - 3)) * x[i];
  }
  for (std::size_t n = 0; n < out.local_size(); ++n) {
    const auto index = out.global_index(n);
    out.set_value(index, y[index[0]]);
  }
}

// A single N x N Jordan block for eigenvalue 0: A e_i = e_{i-1} (the
// superdiagonal is 1, everywhere else 0). This operator is nilpotent
// (A^N = 0, all eigenvalues are exactly 0) and defective (only one
// linearly independent eigenvector regardless of N), which makes it
// adversarial for Krylov methods: resolving more than one Ritz pair
// requires restarts. Used below to force ARPACK to report zero converged
// Ritz values (dneupd ierr = -14) deterministically, by combining it with
// maxiter = 1 (no restarts allowed) and a very tight tolerance.
void matvec_jordan_block(rtensor &out, rtensor const &in, std::size_t N) {
  out = rtensor(in.get_comm(), mptensor::Shape(N));
  std::vector<double> x = gather_real(in, N);
  for (std::size_t n = 0; n < out.local_size(); ++n) {
    const auto index = out.global_index(n);
    const std::size_t i = index[0];
    out.set_value(index, (i + 1 < N) ? x[i + 1] : 0.0);
  }
}

rtensor ones_real(std::size_t N) {
  rtensor v{mptensor::Shape(N)};
  for (std::size_t i = 0; i < N; ++i) {
    v.set_value({i}, 1.0);
  }
  return v;
}

ctensor ones_complex(std::size_t N) {
  ctensor v{mptensor::Shape(N)};
  for (std::size_t i = 0; i < N; ++i) {
    v.set_value({i}, std::complex<double>(1.0, 0.0));
  }
  return v;
}

}  // namespace

TEST_CASE("arpack_eigenvalues, real") {
  const std::size_t N = 20;
  const std::size_t nev = 2;
  auto A = [](rtensor &out, rtensor const &in) { matvec_real(out, in, N); };

  auto ev =
      tenes::arpack_eigenvalues<rtensor>(A, ones_real(N), nev, 10, 10, 1.0e-10);
  REQUIRE(ev.size() == nev);
  CHECK(std::abs(ev[0]) == doctest::Approx(16.0).epsilon(1.0e-8));
  CHECK(std::abs(ev[1]) == doctest::Approx(8.0).epsilon(1.0e-8));
  CHECK(ev[0].imag() == doctest::Approx(0.0).epsilon(1.0e-8));
}

TEST_CASE("arpack_eigenvalues, complex") {
  const std::size_t N = 20;
  const std::size_t nev = 2;
  auto A = [](ctensor &out, ctensor const &in) { matvec_complex(out, in, N); };

  auto ev = tenes::arpack_eigenvalues<ctensor>(A, ones_complex(N), nev, 10, 10,
                                               1.0e-10);
  REQUIRE(ev.size() == nev);
  // eigenvalues are purely imaginary: 16i, 8i, ...
  CHECK(ev[0].real() == doctest::Approx(0.0).epsilon(1.0e-8));
  CHECK(ev[0].imag() == doctest::Approx(16.0).epsilon(1.0e-8));
  CHECK(std::abs(ev[1]) == doctest::Approx(8.0).epsilon(1.0e-8));
}

TEST_CASE("arpack agrees with the builtin Arnoldi") {
  const std::size_t N = 20;
  const std::size_t nev = 2;
  auto A = [](rtensor &out, rtensor const &in) { matvec_real(out, in, N); };

  auto ev_arpack =
      tenes::arpack_eigenvalues<rtensor>(A, ones_real(N), nev, 10, 10, 1.0e-10);

  tenes::Arnoldi<rtensor> arnoldi(N, 10);
  arnoldi.initialize(ones_real(N));
  arnoldi.run(A, nev, 5, 10, 1.0e-8);
  auto ev_builtin = arnoldi.eigenvalues();

  for (std::size_t i = 0; i < nev; ++i) {
    CHECK(std::abs(ev_arpack[i]) ==
          doctest::Approx(std::abs(ev_builtin[i])).epsilon(1.0e-6));
  }
}

TEST_CASE(
    "arpack_eigenvalues keeps the dominant eigenvalue when a conjugate pair "
    "makes ARPACK converge nev + 1 Ritz values") {
  const std::size_t N = 20;
  const std::size_t nev = 2;
  auto A = [](rtensor &out, rtensor const &in) {
    matvec_conjugate_pair(out, in, N);
  };

  auto ev =
      tenes::arpack_eigenvalues<rtensor>(A, ones_real(N), nev, 10, 50, 1.0e-10);
  REQUIRE(ev.size() == nev);
  // dominant, isolated eigenvalue must survive the nev+1 -> nev truncation
  CHECK(std::abs(ev[0]) == doctest::Approx(10.0).epsilon(1.0e-6));
  // one of the 8 +/- 1i conjugate pair
  CHECK(std::abs(ev[1]) == doctest::Approx(std::sqrt(65.0)).epsilon(1.0e-6));
}

TEST_CASE(
    "arpack_eigenvalues returns a NaN-padded vector (not an exception) when "
    "ARPACK converges zero Ritz values") {
  // A Jordan block with maxiter = 1 (no restarts) and a very tight tolerance
  // reliably makes ARPACK-NG's dneupd report ierr = -14 ("no eigenvalues of
  // sufficient accuracy") on the first call, since resolving nev > 1 Ritz
  // pairs of a defective operator needs at least one restart. Verified to
  // reproduce deterministically (same nconv = 0 every time, no randomness
  // involved) via a standalone probe against libarpack directly before
  // being written up as this test.
  const std::size_t N = 12;
  const std::size_t nev = 4;
  auto A = [](rtensor &out, rtensor const &in) {
    matvec_jordan_block(out, in, N);
  };

  auto ev =
      tenes::arpack_eigenvalues<rtensor>(A, ones_real(N), nev, 6, 1, 1.0e-14);
  REQUIRE(ev.size() == nev);
  for (std::size_t i = 0; i < nev; ++i) {
    CHECK(std::isnan(ev[i].real()));
    CHECK(std::isnan(ev[i].imag()));
  }
}
