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

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <cmath>
#include <complex>
#include <cstddef>

#include "../src/arnoldi.hpp"
#include "../src/tensor.hpp"
#include "../src/mpi.hpp"

namespace {

using tensor = tenes::real_tensor;

// A x for the diagonal matrix A = diag(2^0, 2^-1, ..., 2^-(N-1)), scaled by 16
void matvec(tensor &out, tensor const &in, std::size_t N) {
  out = tensor(mptensor::Shape(N));
  for (std::size_t i = 0; i < N; ++i) {
    double v;
    in.get_value({i}, v);
    out.set_value({i}, 16.0 * std::pow(2.0, -static_cast<double>(i)) * v);
  }
}

tensor ones(std::size_t N) {
  tensor v{mptensor::Shape(N)};
  for (std::size_t i = 0; i < N; ++i) {
    v.set_value({i}, 1.0);
  }
  return v;
}

}  // namespace

TEST_CASE("Arnoldi") {
  const std::size_t N = 20;
  const std::size_t maxvec = 10;
  const std::size_t nev = 2;

  auto A = [](tensor &out, tensor const &in) { matvec(out, in, N); };

  SUBCASE("eigenvalues of a diagonal matrix") {
    tenes::Arnoldi<tensor> arnoldi(N, maxvec);
    arnoldi.initialize(ones(N));
    arnoldi.run(A, nev, 5, 10, 1.0e-8);
    auto ev = arnoldi.eigenvalues();
    REQUIRE(ev.size() >= nev);
    CHECK(std::abs(ev[0]) == doctest::Approx(16.0).epsilon(1.0e-6));
    CHECK(std::abs(ev[1]) == doctest::Approx(8.0).epsilon(1.0e-6));
  }

  SUBCASE("mindim larger than maxvec must not hang") {
    tenes::Arnoldi<tensor> arnoldi(N, maxvec);
    arnoldi.initialize(ones(N));
    // an unsatisfiable tolerance forces restarts
    arnoldi.run(A, nev, 100, 3, 1.0e-300);
    auto ev = arnoldi.eigenvalues();
    REQUIRE(ev.size() >= nev);
    CHECK(std::abs(ev[0]) == doctest::Approx(16.0).epsilon(1.0e-6));
  }
}
