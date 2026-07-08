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

#include "../src/iTPS/correlation_length.hpp"

TEST_CASE("calc_correlation_length") {
  using tenes::itps::calc_correlation_length;

  SUBCASE("ordinary values") {
    // e1/e0 = exp(-2)  ->  xi = L / 2
    CHECK(calc_correlation_length(1.0, std::exp(-2.0), 1) ==
          doctest::Approx(0.5));
    CHECK(calc_correlation_length(1.0, std::exp(-2.0), 3) ==
          doctest::Approx(1.5));
    // overall scale of the eigenvalues does not matter
    CHECK(calc_correlation_length(2.0, 2.0 * std::exp(-2.0), 1) ==
          doctest::Approx(0.5));
  }

  SUBCASE("degenerate leading eigenvalues give +infinity") {
    const double xi = calc_correlation_length(1.0, 1.0, 4);
    CHECK(std::isinf(xi));
    CHECK(xi > 0.0);
  }

  SUBCASE("e1 larger than e0 also gives +infinity") {
    const double xi = calc_correlation_length(1.0, 1.5, 4);
    CHECK(std::isinf(xi));
    CHECK(xi > 0.0);
  }

  SUBCASE("vanishing leading eigenvalue gives 0") {
    CHECK(calc_correlation_length(0.0, 0.0, 4) == 0.0);
  }

  SUBCASE("vanishing subleading eigenvalue gives 0") {
    CHECK(calc_correlation_length(1.0, 0.0, 4) == 0.0);
  }

  SUBCASE("result is never NaN") {
    CHECK_FALSE(std::isnan(calc_correlation_length(0.0, 0.0, 4)));
    CHECK_FALSE(std::isnan(calc_correlation_length(1.0, 1.0, 4)));
    CHECK_FALSE(std::isnan(calc_correlation_length(1.0, 0.0, 4)));
  }
}
