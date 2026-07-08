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

#include "../src/tensor.hpp"
#include "../src/exception.hpp"
#include "../src/mpi.hpp"

TEST_CASE("identity") {
  using mptensor::Index;

  const std::size_t k = 3;
  const double coeff = 2.0;
  auto I = tenes::identity(k, coeff);
  REQUIRE(I.shape() == mptensor::Shape(k, k));
  double v;
  for (std::size_t i = 0; i < k; ++i) {
    for (std::size_t j = 0; j < k; ++j) {
      I.get_value(Index(i, j), v);
      CHECK(v == (i == j ? coeff : 0.0));
    }
  }
}

TEST_CASE("resize_tensor") {
  using tensor = tenes::real_tensor;
  using mptensor::Index;
  using mptensor::Shape;

  tensor src(Shape(2, 2));
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      src.set_value(Index(i, j), 10.0 * i + j);
    }
  }

  SUBCASE("ndim mismatch throws") {
    CHECK_THROWS_AS(tenes::resize_tensor(src, Shape(2, 2, 2)),
                    tenes::logic_error);
  }

  SUBCASE("extend pads with zero") {
    tensor A = tenes::resize_tensor(src, Shape(3, 3));
    REQUIRE(A.shape() == Shape(3, 3));
    double v;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        A.get_value(Index(i, j), v);
        const double expected = (i < 2 && j < 2) ? 10.0 * i + j : 0.0;
        CHECK(v == expected);
      }
    }
  }

  SUBCASE("shrink keeps the leading block") {
    tensor A = tenes::resize_tensor(src, Shape(1, 2));
    REQUIRE(A.shape() == Shape(1, 2));
    double v;
    for (int j = 0; j < 2; ++j) {
      A.get_value(Index(0, j), v);
      CHECK(v == static_cast<double>(j));
    }
  }
}
