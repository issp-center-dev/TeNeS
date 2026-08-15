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

#include "../src/fermion/parity.hpp"

using namespace tenes::fermion;

TEST_CASE("fuse combines parities row-major with XOR") {
  parity_vector a{false, true};
  parity_vector b{false, true, true};
  parity_vector f = fuse(a, b);
  REQUIRE(f.size() == 6);
  CHECK(f == parity_vector{false, true, true, true, false, false});
}

TEST_CASE("count_odd counts odd legs at an element") {
  leg_parities p{{false, true}, {false, true}, {false}};
  CHECK(count_odd(p, mptensor::Index(1, 1, 0)) == 2);
  CHECK(count_odd(p, mptensor::Index(0, 1, 0)) == 1);
}

TEST_CASE("parity_sort_perm is stable even-first") {
  parity_vector p{true, false, true, false};
  CHECK(parity_sort_perm(p) == std::vector<size_t>{1, 3, 0, 2});
}
