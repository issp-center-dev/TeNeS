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

#include <random>

#include "../src/fermion/fops.hpp"
#include "../src/fermion/ftensor.hpp"
#include "../src/fermion/parity.hpp"
#include "../src/tensor.hpp"

using namespace tenes::fermion;

using ft = tenes::fermion::ftensor<tenes::real_tensor>;

static ft make_random_ft(const mptensor::Shape& sh,
                         const tenes::fermion::leg_parities& p, unsigned seed) {
  tenes::real_tensor t(sh);
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    t.set_value(t.global_index(n), dist(gen));
  }
  return ft{t, p};
}

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

TEST_CASE("apply_swap negates odd-odd elements only") {
  tenes::fermion::leg_parities p{{false, true}, {false, true}, {false, true}};
  ft a = make_random_ft(mptensor::Shape(2, 2, 2), p, 42);
  ft b = a;
  tenes::fermion::apply_swap(b, 0, 2);
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    auto idx = a.t.global_index(n);
    double va, vb;
    a.t.get_value(idx, va);
    b.t.get_value(idx, vb);
    double expected = (p[0][idx[0]] && p[2][idx[2]]) ? -va : va;
    CHECK(vb == doctest::Approx(expected));
  }
}

TEST_CASE("parity_violation detects odd-total elements") {
  tenes::fermion::leg_parities p{{false, true}, {false, true}};
  ft a{tenes::real_tensor(mptensor::Shape(2, 2)), p};
  a.t.set_value(mptensor::Index(0, 0), 1.0);
  CHECK(tenes::fermion::parity_violation(a) == doctest::Approx(0.0));
  a.t.set_value(mptensor::Index(0, 1), 0.5);
  CHECK(tenes::fermion::parity_violation(a) == doctest::Approx(0.5));
}

static int inversion_sign(const tenes::fermion::leg_parities& p,
                          const mptensor::Index& idx,
                          const mptensor::Axes& axes) {
  int k = 0;
  for (std::size_t x = 0; x < axes.size(); ++x) {
    for (std::size_t y = x + 1; y < axes.size(); ++y) {
      if (axes[x] > axes[y] && p[axes[y]][idx[axes[y]]] &&
          p[axes[x]][idx[axes[x]]]) {
        ++k;
      }
    }
  }
  return (k % 2 == 0) ? 1 : -1;
}

TEST_CASE("graded transpose applies Koszul signs") {
  tenes::fermion::leg_parities p{
      {false, true}, {false, true, true}, {false, true}};
  ft a = make_random_ft(mptensor::Shape(2, 3, 2), p, 7);
  mptensor::Axes ax(2, 0, 1);
  ft b = tenes::fermion::transpose(a, ax);
  REQUIRE(b.parity[0] == p[2]);
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    auto idx = a.t.global_index(n);
    double va, vb;
    a.t.get_value(idx, va);
    b.t.get_value(mptensor::Index(idx[2], idx[0], idx[1]), vb);
    CHECK(vb == doctest::Approx(inversion_sign(p, idx, ax) * va));
  }
}

TEST_CASE("transpose twice returns original") {
  tenes::fermion::leg_parities p{{false, true}, {true, false}, {false, true}};
  ft a = make_random_ft(mptensor::Shape(2, 2, 2), p, 8);
  ft b = tenes::fermion::transpose(
      tenes::fermion::transpose(a, mptensor::Axes(1, 2, 0)),
      mptensor::Axes(2, 0, 1));
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    auto idx = a.t.global_index(n);
    double va, vb;
    a.t.get_value(idx, va);
    b.t.get_value(idx, vb);
    CHECK(vb == doctest::Approx(va));
  }
}

static ft make_even_ft(const mptensor::Shape& sh,
                       const tenes::fermion::leg_parities& p, unsigned seed) {
  ft a = make_random_ft(sh, p, seed);
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    auto idx = a.t.global_index(n);
    if (tenes::fermion::count_odd(p, idx) % 2 == 1) {
      a.t.set_value(idx, 0.0);
    }
  }
  return a;
}

TEST_CASE("tensordot rejects parity mismatch") {
  ft a =
      make_random_ft(mptensor::Shape(2, 2), {{false, true}, {false, true}}, 1);
  ft b =
      make_random_ft(mptensor::Shape(2, 2), {{false, false}, {false, true}}, 2);
  CHECK_THROWS_AS(
      tenes::fermion::tensordot(a, b, mptensor::Axes(1), mptensor::Axes(0)),
      std::runtime_error);
}

TEST_CASE("tensordot is contraction-order independent") {
  tenes::fermion::leg_parities pab{{false, true}, {false, true}};
  ft a = make_even_ft(mptensor::Shape(2, 2), pab, 11);
  ft b = make_even_ft(mptensor::Shape(2, 2, 2),
                      {{false, true}, {false, true}, {false, true}}, 12);
  ft c = make_even_ft(mptensor::Shape(2, 2), pab, 13);
  ft ab = tenes::fermion::tensordot(a, b, mptensor::Axes(1), mptensor::Axes(0));
  ft abc1 =
      tenes::fermion::tensordot(ab, c, mptensor::Axes(2), mptensor::Axes(0));
  ft bc = tenes::fermion::tensordot(b, c, mptensor::Axes(2), mptensor::Axes(0));
  ft abc2 =
      tenes::fermion::tensordot(a, bc, mptensor::Axes(1), mptensor::Axes(0));
  for (std::size_t n = 0; n < abc1.t.local_size(); ++n) {
    auto idx = abc1.t.global_index(n);
    double v1, v2;
    abc1.t.get_value(idx, v1);
    abc2.t.get_value(idx, v2);
    CHECK(v1 == doctest::Approx(v2));
  }
}

TEST_CASE("tensordot matches manual swap for a crossing contraction") {
  tenes::fermion::leg_parities p{{false, true}, {false, true}};
  ft a = make_random_ft(mptensor::Shape(2, 2), p, 31);
  ft b = make_random_ft(mptensor::Shape(2, 2), p, 32);
  ft got =
      tenes::fermion::tensordot(a, b, mptensor::Axes(1), mptensor::Axes(1));
  ft b_swapped = b;
  tenes::fermion::apply_swap(b_swapped, 0, 1);
  tenes::real_tensor expected = mptensor::tensordot(
      a.t, b_swapped.t, mptensor::Axes(1), mptensor::Axes(1));
  for (std::size_t n = 0; n < got.t.local_size(); ++n) {
    auto idx = got.t.global_index(n);
    double vg, ve;
    got.t.get_value(idx, vg);
    expected.get_value(idx, ve);
    CHECK(vg == doctest::Approx(ve));
  }
}

TEST_CASE("trace matches manual swap for reversed contracted ordering") {
  tenes::fermion::leg_parities p{{false, true}, {false, true}};
  ft a = make_random_ft(mptensor::Shape(2, 2), p, 41);
  ft b = make_random_ft(mptensor::Shape(2, 2), p, 42);
  double got =
      tenes::fermion::trace(a, b, mptensor::Axes(0, 1), mptensor::Axes(0, 1));
  ft b_swapped = b;
  tenes::fermion::apply_swap(b_swapped, 0, 1);
  double expected = mptensor::trace(a.t, b_swapped.t, mptensor::Axes(0, 1),
                                    mptensor::Axes(0, 1));
  CHECK(got == doctest::Approx(expected));
}
