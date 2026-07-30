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

// A x for the diagonal matrix A = diag(2^0, 2^-1, ..., 2^-(N-1)), scaled by 16
void matvec_real(rtensor &out, rtensor const &in, std::size_t N) {
  out = rtensor(mptensor::Shape(N));
  for (std::size_t i = 0; i < N; ++i) {
    double v;
    in.get_value({i}, v);
    out.set_value({i}, 16.0 * std::pow(2.0, -static_cast<double>(i)) * v);
  }
}

// the same spectrum rotated to the imaginary axis: A = diag(16i * 2^-i)
void matvec_complex(ctensor &out, ctensor const &in, std::size_t N) {
  out = ctensor(mptensor::Shape(N));
  const std::complex<double> I(0.0, 1.0);
  for (std::size_t i = 0; i < N; ++i) {
    std::complex<double> v;
    in.get_value({i}, v);
    out.set_value({i},
                  16.0 * I * std::pow(2.0, -static_cast<double>(i)) * v);
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

  auto ev = tenes::arpack_eigenvalues<rtensor>(A, ones_real(N), nev, 10, 10,
                                               1.0e-10);
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

  auto ev_arpack = tenes::arpack_eigenvalues<rtensor>(A, ones_real(N), nev, 10,
                                                      10, 1.0e-10);

  tenes::Arnoldi<rtensor> arnoldi(N, 10);
  arnoldi.initialize(ones_real(N));
  arnoldi.run(A, nev, 5, 10, 1.0e-8);
  auto ev_builtin = arnoldi.eigenvalues();

  for (std::size_t i = 0; i < nev; ++i) {
    CHECK(std::abs(ev_arpack[i]) ==
          doctest::Approx(std::abs(ev_builtin[i])).epsilon(1.0e-6));
  }
}
