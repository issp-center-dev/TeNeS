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

#include <string>
#include <vector>

#include "../src/util/string.hpp"
#include "../src/util/file.hpp"
#include "../src/util/read_tensor.hpp"
#include "../src/tensor.hpp"
#include "../src/mpi.hpp"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  doctest::Context context(argc, argv);
  const int res = context.run();
  MPI_Finalize();
  return res;
}

TEST_CASE("split") {
  using tenes::util::split;

  SUBCASE("words do not contain delimiters") {
    auto words = split("a b c");
    REQUIRE(words.size() == 3);
    CHECK(words[0] == "a");
    CHECK(words[1] == "b");
    CHECK(words[2] == "c");
  }

  SUBCASE("consecutive, leading, and trailing delimiters are skipped") {
    auto words = split("  1\t\t22  333 ");
    REQUIRE(words.size() == 3);
    CHECK(words[0] == "1");
    CHECK(words[1] == "22");
    CHECK(words[2] == "333");
  }

  SUBCASE("empty or delimiter-only input yields no words") {
    CHECK(split("").empty());
    CHECK(split(" \t ").empty());
  }

  SUBCASE("single word") {
    auto words = split("word");
    REQUIRE(words.size() == 1);
    CHECK(words[0] == "word");
  }
}

TEST_CASE("strip and drop_comment") {
  using namespace tenes::util;

  CHECK(strip("  x  ") == "x");
  CHECK(lstrip("  x") == "x");
  CHECK(rstrip("x  ") == "x");
  CHECK(drop_comment("1 2 # comment") == "1 2 ");
  CHECK(drop_comment("no comment") == "no comment");
  CHECK(drop_comment("# all comment") == "");
}

TEST_CASE("basename") {
  using tenes::util::basename;

  CHECK(basename("foo/bar.dat") == "bar.dat");
  CHECK(basename("bar.dat") == "bar.dat");
  CHECK(basename("") == "");
}

TEST_CASE("read_tensor") {
  using tensor = tenes::real_tensor;
  using mptensor::Index;
  using mptensor::Shape;

  SUBCASE("blank and comment lines are skipped") {
    const std::string str =
        "0 0 1.0 0.0\n"
        "\n"
        "   \n"
        "# comment only\n"
        "1 1 2.0 0.0  # trailing comment\n";
    tensor A = tenes::util::read_tensor<tensor>(str, Shape(2, 2));
    double v;
    A.get_value(Index(0, 0), v);
    CHECK(v == 1.0);
    A.get_value(Index(1, 1), v);
    CHECK(v == 2.0);
    A.get_value(Index(0, 1), v);
    CHECK(v == 0.0);
  }

  SUBCASE("wrong number of columns throws") {
    CHECK_THROWS_AS(
        tenes::util::read_tensor<tensor>("0 0 1.0\n", mptensor::Shape(2, 2)),
        tenes::input_error);
  }
}
