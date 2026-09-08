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

// ===== is_finite_value() under finite-math-only ============================
//
// This translation unit is compiled with -ffast-math on purpose. Under it
// (and under icpx's default -fp-model=fast, which CMakeLists.txt has to
// override for exactly this reason) the compiler is allowed to assume no
// operand is Inf or NaN, and std::isfinite folds to a constant true.
//
// The predicate guarded here is the one that tells a user their state had
// already diverged before a decomposition was reached. A diagnostic that
// the optimizer is allowed to delete is worse than none, because it reports
// "all elements finite" about a matrix full of NaN - which is what a real
// oneAPI 2022.2.1 run did. So it must not be written in terms the
// finite-math assumption can reason about.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../doctest.h"

#include <complex>
#include <cstdint>
#include <cstring>
#include <limits>

#include "../../src/fermion/finite_check.hpp"

namespace {

//! Build a value from its bits, so that no compile-time constant the
//! optimizer could fold reaches the predicate.
double from_bits(std::uint64_t bits) {
  double v = 0.0;
  std::memcpy(&v, &bits, sizeof(v));
  return v;
}

}  // namespace

TEST_CASE("is_finite_value accepts finite values under -ffast-math") {
  namespace fd = tenes::fermion::detail;
  CHECK(fd::is_finite_value(from_bits(0x3ff0000000000000ULL)));  // 1.0
  CHECK(fd::is_finite_value(from_bits(0x0000000000000000ULL)));  // +0.0
  CHECK(fd::is_finite_value(from_bits(0x8000000000000000ULL)));  // -0.0
  CHECK(fd::is_finite_value(from_bits(0x7fefffffffffffffULL)));  // DBL_MAX
  CHECK(fd::is_finite_value(from_bits(0x0000000000000001ULL)));  // denormal
}

TEST_CASE("is_finite_value still rejects Inf and NaN under -ffast-math") {
  namespace fd = tenes::fermion::detail;
  CHECK_FALSE(fd::is_finite_value(from_bits(0x7ff0000000000000ULL)));  // +Inf
  CHECK_FALSE(fd::is_finite_value(from_bits(0xfff0000000000000ULL)));  // -Inf
  CHECK_FALSE(fd::is_finite_value(from_bits(0x7ff8000000000000ULL)));  // qNaN
  CHECK_FALSE(fd::is_finite_value(from_bits(0x7ff0000000000001ULL)));  // sNaN
}

TEST_CASE("is_finite_value rejects a complex value with one bad part") {
  namespace fd = tenes::fermion::detail;
  const double nan = from_bits(0x7ff8000000000000ULL);
  const double inf = from_bits(0x7ff0000000000000ULL);
  CHECK(fd::is_finite_value(std::complex<double>(1.0, -2.0)));
  CHECK_FALSE(fd::is_finite_value(std::complex<double>(nan, 0.0)));
  CHECK_FALSE(fd::is_finite_value(std::complex<double>(0.0, inf)));
}
