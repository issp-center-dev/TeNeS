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

#include <iostream>
#include <random>

#include "../src/fermion/fermion_info.hpp"
#include "../src/fermion/fops.hpp"
#include "../src/fermion/ftensor.hpp"
#include "../src/fermion/parity.hpp"
#include "../src/SquareLattice.hpp"
#include "../src/tensor.hpp"
#include "../src/iTPS/PEPS_Parameters.hpp"
#include "../src/iTPS/core/ctm.hpp"
#include "../src/iTPS/core/contract_itps_ctm.hpp"

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

TEST_CASE("fuse combines parities column-major with XOR") {
  parity_vector a{false, true};
  parity_vector b{false, true, true};
  parity_vector f = fuse(a, b);
  REQUIRE(f.size() == 6);
  // mptensor flattening: index = i + a.size() * j (first leg fastest)
  CHECK(f == parity_vector{false, true, true, false, true, false});
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

TEST_CASE("norm of even vector tensor is positive via conj") {
  tenes::fermion::leg_parities p{{false, false, true, true}};
  ft v = make_random_ft(mptensor::Shape(4), p, 21);
  ft cv = tenes::fermion::conj(v);
  double nrm =
      tenes::fermion::trace(cv, v, mptensor::Axes(0), mptensor::Axes(0));
  CHECK(nrm > 0.0);
}

TEST_CASE("conj applies candidate two-odd-leg twist") {
  tenes::fermion::leg_parities p{{false, true}, {false, true}};
  ft a{tenes::real_tensor(mptensor::Shape(2, 2)), p};
  a.t.set_value(mptensor::Index(1, 1), 2.0);
  ft ca = tenes::fermion::conj(a);
  double v;
  ca.t.get_value(mptensor::Index(1, 1), v);
  CHECK(v == doctest::Approx(-2.0));
}

TEST_CASE("conj is an involution up to elementwise conj") {
  tenes::fermion::leg_parities p{{false, true}, {false, true}};
  ft a = make_random_ft(mptensor::Shape(2, 2), p, 22);
  ft cca = tenes::fermion::conj(tenes::fermion::conj(a));
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    auto idx = a.t.global_index(n);
    double va, vb;
    a.t.get_value(idx, va);
    cca.t.get_value(idx, vb);
    CHECK(vb == doctest::Approx(va));
  }
}

TEST_CASE("slice carries the sliced leg parity interval") {
  tenes::fermion::leg_parities p{{false, true, false}, {false, true}};
  ft a = make_random_ft(mptensor::Shape(3, 2), p, 51);
  ft s = tenes::fermion::slice(a, 0, 1, 3);
  REQUIRE(s.shape() == mptensor::Shape(2, 2));
  CHECK(s.parity[0] == tenes::fermion::parity_vector{true, false});
  for (std::size_t n = 0; n < s.t.local_size(); ++n) {
    auto idx = s.t.global_index(n);
    double vs, va;
    s.t.get_value(idx, vs);
    a.t.get_value(mptensor::Index(idx[0] + 1, idx[1]), va);
    CHECK(vs == doctest::Approx(va));
  }
}

TEST_CASE("extend pads new parity entries as even") {
  tenes::fermion::leg_parities p{{false, true}, {true}};
  ft a = make_random_ft(mptensor::Shape(2, 1), p, 52);
  ft e = tenes::fermion::extend(a, mptensor::Shape(3, 2));
  CHECK(e.shape() == mptensor::Shape(3, 2));
  CHECK(e.parity[0] == tenes::fermion::parity_vector{false, true, false});
  CHECK(e.parity[1] == tenes::fermion::parity_vector{true, false});
}

TEST_CASE("FermionInfo wraps and unwraps Tn roundtrip") {
  tenes::fermion::FermionInfo fi;
  fi.enabled = true;
  fi.phys = {{false, true}};
  fi.virt = {{{tenes::fermion::parity_vector{false, true},
               tenes::fermion::parity_vector{false, true},
               tenes::fermion::parity_vector{false, true},
               tenes::fermion::parity_vector{false, true}}}};

  tenes::real_tensor raw(mptensor::Shape(2, 2, 2, 2, 2));
  for (std::size_t n = 0; n < raw.local_size(); ++n) {
    raw.set_value(raw.global_index(n), static_cast<double>(n + 1));
  }

  ft wrapped = tenes::fermion::wrap_Tn(raw, fi, 0);
  REQUIRE(wrapped.parity.size() == 5);
  CHECK(wrapped.parity[0] == fi.virt[0][0]);
  CHECK(wrapped.parity[4] == fi.phys[0]);

  tenes::real_tensor unwrapped(raw.shape());
  tenes::fermion::FermionInfo fi_out = fi;
  tenes::fermion::unwrap_Tn(wrapped, unwrapped, fi_out, 0);
  for (std::size_t n = 0; n < raw.local_size(); ++n) {
    auto idx = raw.global_index(n);
    double before, after;
    raw.get_value(idx, before);
    unwrapped.get_value(idx, after);
    CHECK(after == doctest::Approx(before));
  }
  CHECK(fi_out.virt[0][2] == fi.virt[0][2]);
  CHECK(fi_out.phys[0] == fi.phys[0]);
}

TEST_CASE("reshape fuses adjacent leg parities") {
  tenes::fermion::leg_parities p{
      {false, true}, {false, true, true}, {true, false}};
  ft a = make_random_ft(mptensor::Shape(2, 3, 2), p, 53);
  ft r = tenes::fermion::reshape(a, mptensor::Shape(6, 2));
  REQUIRE(r.shape() == mptensor::Shape(6, 2));
  REQUIRE(r.parity.size() == 2);
  CHECK(r.parity[0] == tenes::fermion::fuse(p[0], p[1]));
  CHECK(r.parity[1] == p[2]);
  // column-major fusion: fused index 4 = i + 2*j -> (i, j) = (0, 2)
  CHECK(tenes::fermion::count_odd(r.parity, mptensor::Index(4, 0)) ==
        (static_cast<int>(p[0][0] != p[1][2]) + static_cast<int>(p[2][0])));
}

TEST_CASE("reshape fusion preserves contraction values") {
  tenes::fermion::leg_parities p{{false, true}, {false, true}, {false, true}};
  ft a = make_even_ft(mptensor::Shape(2, 2, 2), p, 54);
  ft b = make_even_ft(mptensor::Shape(4, 2),
                      {tenes::fermion::fuse(p[0], p[1]), {false, true}}, 55);
  ft ar = tenes::fermion::reshape(a, mptensor::Shape(4, 2));
  ft got =
      tenes::fermion::tensordot(ar, b, mptensor::Axes(0), mptensor::Axes(0));
  for (std::size_t n = 0; n < got.t.local_size(); ++n) {
    auto idx = got.t.global_index(n);
    double vg;
    got.t.get_value(idx, vg);
    double vr = 0.0;
    for (std::size_t fused = 0; fused < 4; ++fused) {
      double va, vb;
      ar.t.get_value(mptensor::Index(fused, idx[0]), va);
      b.t.get_value(mptensor::Index(fused, idx[1]), vb);
      const int sign = (ar.parity[0][fused] && ar.parity[1][idx[0]]) ? -1 : 1;
      vr += sign * va * vb;
    }
    CHECK(vg == doctest::Approx(vr));
  }
}

TEST_CASE("max_abs and multiply_vector forward to the dense tensor") {
  tenes::fermion::leg_parities p{{false, true}, {false, true}, {false, true}};
  ft a = make_random_ft(mptensor::Shape(2, 2, 2), p, 56);
  ft b = a;
  std::vector<double> v0{2.0, 3.0};
  std::vector<double> v1{5.0, 7.0};
  std::vector<double> v2{11.0, 13.0};
  b.multiply_vector(v0, 0, v1, 1, v2, 2);
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    auto idx = a.t.global_index(n);
    double va, vb;
    a.t.get_value(idx, va);
    b.t.get_value(idx, vb);
    CHECK(vb == doctest::Approx(va * v0[idx[0]] * v1[idx[1]] * v2[idx[2]]));
  }
  CHECK(tenes::fermion::max_abs(b) == doctest::Approx(mptensor::max_abs(b.t)));
}

TEST_CASE("make_perm_matrix uses new-position to old-position convention") {
  std::vector<std::size_t> perm{1, 3, 0, 2};
  tenes::real_tensor pmat =
      tenes::fermion::make_perm_matrix<tenes::real_tensor>(perm);
  for (std::size_t i = 0; i < 4; ++i) {
    for (std::size_t j = 0; j < 4; ++j) {
      double v;
      pmat.get_value(mptensor::Index(i, j), v);
      CHECK(v == doctest::Approx(perm[i] == j ? 1.0 : 0.0));
    }
  }
}

TEST_CASE("graded QR reconstructs an even tensor with block parity") {
  tenes::fermion::leg_parities p{{false, true, false}, {false, true, true}};
  ft a = make_even_ft(mptensor::Shape(3, 3), p, 61);
  ft q, r;
  int info = tenes::fermion::qr(a, mptensor::Axes(0), mptensor::Axes(1), q, r);
  CHECK(info == 0);
  CHECK(tenes::fermion::parity_violation(q) == doctest::Approx(0.0));
  CHECK(tenes::fermion::parity_violation(r) == doctest::Approx(0.0));
  ft recon =
      tenes::fermion::tensordot(q, r, mptensor::Axes(1), mptensor::Axes(0));
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    auto idx = a.t.global_index(n);
    double va, vr;
    a.t.get_value(idx, va);
    recon.t.get_value(idx, vr);
    CHECK(vr == doctest::Approx(va).epsilon(1.0e-10));
  }
  ft qtq = tenes::fermion::tensordot(tenes::fermion::conj(q), q,
                                     mptensor::Axes(0), mptensor::Axes(0));
  for (std::size_t n = 0; n < qtq.t.local_size(); ++n) {
    auto idx = qtq.t.global_index(n);
    double v;
    qtq.t.get_value(idx, v);
    CHECK(v == doctest::Approx(idx[0] == idx[1] ? 1.0 : 0.0).epsilon(1.0e-10));
  }
}

TEST_CASE("graded SVD reconstructs an even tensor with block parity") {
  tenes::fermion::leg_parities p{{false, true, false}, {false, true, true}};
  ft a = make_even_ft(mptensor::Shape(3, 3), p, 71);
  ft u, vt;
  std::vector<double> s;
  int info =
      tenes::fermion::svd(a, mptensor::Axes(0), mptensor::Axes(1), u, s, vt);
  CHECK(info == 0);
  CHECK(tenes::fermion::parity_violation(u) == doctest::Approx(0.0));
  CHECK(tenes::fermion::parity_violation(vt) == doctest::Approx(0.0));
  ft us = u;
  us.multiply_vector(s, 1);
  ft recon =
      tenes::fermion::tensordot(us, vt, mptensor::Axes(1), mptensor::Axes(0));
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    auto idx = a.t.global_index(n);
    double va, vr;
    a.t.get_value(idx, va);
    recon.t.get_value(idx, vr);
    CHECK(vr == doctest::Approx(va).epsilon(1.0e-10));
  }
}

TEST_CASE("svd_trunc selects across sectors and returns even-first metadata") {
  tenes::fermion::leg_parities p{{false, false, true}, {false, false, true}};
  ft a{tenes::real_tensor(mptensor::Shape(3, 3)), p};
  a.t.set_value(mptensor::Index(0, 0), 3.0);
  a.t.set_value(mptensor::Index(1, 1), 1.0);
  a.t.set_value(mptensor::Index(2, 2), 2.0);
  ft u, vt;
  std::vector<double> s;
  int info = tenes::fermion::svd_trunc(a, mptensor::Axes(0), mptensor::Axes(1),
                                       u, s, vt, 2);
  CHECK(info == 0);
  REQUIRE(s.size() == 2);
  CHECK(s[0] == doctest::Approx(3.0));
  CHECK(s[1] == doctest::Approx(2.0));
  REQUIRE(u.parity.size() == 2);
  CHECK(u.parity[1] == tenes::fermion::parity_vector{false, true});
  CHECK(vt.parity[0] == tenes::fermion::parity_vector{false, true});
  CHECK(tenes::fermion::parity_violation(u) == doctest::Approx(0.0));
  CHECK(tenes::fermion::parity_violation(vt) == doctest::Approx(0.0));
}

TEST_CASE("svd_trunc breaks cross-sector ties by preferring even parity") {
  tenes::fermion::leg_parities p{{false, true}, {false, true}};
  ft a{tenes::real_tensor(mptensor::Shape(2, 2)), p};
  a.t.set_value(mptensor::Index(0, 0), 2.0);
  a.t.set_value(mptensor::Index(1, 1), 2.0);
  ft u, vt;
  std::vector<double> s;
  int info = tenes::fermion::svd_trunc(a, mptensor::Axes(0), mptensor::Axes(1),
                                       u, s, vt, 1);
  CHECK(info == 0);
  REQUIRE(s.size() == 1);
  CHECK(s[0] == doctest::Approx(2.0));
  CHECK(u.parity[1] == tenes::fermion::parity_vector{false});
  CHECK(vt.parity[0] == tenes::fermion::parity_vector{false});
}

static ft make_reference_mps_tensor(int site) {
  const double weights[4][2][2] = {{{1.10, 0.70}, {0.0, 0.0}},
                                   {{0.90, -0.40}, {1.30, 0.80}},
                                   {{1.20, 0.50}, {-0.60, 1.10}},
                                   {{0.75, -0.35}, {1.40, 0.65}}};
  const std::size_t left_dim = (site == 0) ? 1 : 2;
  const std::size_t right_dim = (site == 3) ? 1 : 2;
  ft a{tenes::real_tensor(mptensor::Shape(left_dim, 2, right_dim)),
       {parity_vector(left_dim, false),
        {false, true},
        parity_vector(right_dim, false)}};
  if (site != 0) {
    a.parity[0] = {false, true};
  }
  if (site != 3) {
    a.parity[2] = {false, true};
  }
  for (std::size_t left = 0; left < left_dim; ++left) {
    for (std::size_t occ = 0; occ < 2; ++occ) {
      const std::size_t right = left ^ occ;
      if (right < right_dim) {
        a.t.set_value(mptensor::Index(left, occ, right),
                      weights[site][left][occ]);
      }
    }
  }
  return a;
}

static ft reference_chain_left_to_right() {
  ft a0 = make_reference_mps_tensor(0);
  ft a1 = make_reference_mps_tensor(1);
  ft a2 = make_reference_mps_tensor(2);
  ft a3 = make_reference_mps_tensor(3);
  ft psi =
      tenes::fermion::tensordot(a0, a1, mptensor::Axes(2), mptensor::Axes(0));
  psi =
      tenes::fermion::tensordot(psi, a2, mptensor::Axes(3), mptensor::Axes(0));
  psi =
      tenes::fermion::tensordot(psi, a3, mptensor::Axes(4), mptensor::Axes(0));
  return psi;
}

static ft reference_chain_pairwise() {
  ft a0 = make_reference_mps_tensor(0);
  ft a1 = make_reference_mps_tensor(1);
  ft a2 = make_reference_mps_tensor(2);
  ft a3 = make_reference_mps_tensor(3);
  ft left =
      tenes::fermion::tensordot(a0, a1, mptensor::Axes(2), mptensor::Axes(0));
  ft right =
      tenes::fermion::tensordot(a2, a3, mptensor::Axes(2), mptensor::Axes(0));
  return tenes::fermion::tensordot(left, right, mptensor::Axes(3),
                                   mptensor::Axes(0));
}

static double full_norm(const ft& psi) {
  return tenes::fermion::trace(tenes::fermion::conj(psi), psi,
                               mptensor::Axes(0, 1, 2, 3, 4, 5),
                               mptensor::Axes(0, 1, 2, 3, 4, 5));
}

static ft apply_one_site_op(const ft& psi, int site, const ft& op) {
  const int axis = site + 1;
  ft applied = tenes::fermion::tensordot(psi, op, mptensor::Axes(axis),
                                         mptensor::Axes(0));
  mptensor::Axes perm;
  for (int i = 0; i < axis; ++i) {
    perm.push(i);
  }
  perm.push(5);
  for (int i = axis; i < 5; ++i) {
    perm.push(i);
  }
  return tenes::fermion::transpose(applied, perm);
}

static ft apply_two_site_op_01(const ft& psi, const ft& op) {
  ft applied = tenes::fermion::tensordot(psi, op, mptensor::Axes(1, 2),
                                         mptensor::Axes(0, 1));
  return tenes::fermion::transpose(applied, mptensor::Axes(0, 4, 5, 1, 2, 3));
}

TEST_CASE("JW four-site reference matches f-primitive contractions") {
  constexpr double ref_norm = 2.0353391950000002;
  constexpr double ref_n0 = 0.47988049112374104;
  constexpr double ref_n1 = 0.2202854252015719;
  constexpr double ref_n2 = 0.40879752477817338;
  constexpr double ref_n3 = 0.14348843338616096;
  constexpr double ref_hop01 = -0.32763393032383481;

  ft psi_ltr = reference_chain_left_to_right();
  ft psi_pair = reference_chain_pairwise();
  for (std::size_t n = 0; n < psi_ltr.t.local_size(); ++n) {
    auto idx = psi_ltr.t.global_index(n);
    double vl, vp;
    psi_ltr.t.get_value(idx, vl);
    psi_pair.t.get_value(idx, vp);
    CHECK(vl == doctest::Approx(vp).epsilon(1.0e-12));
  }

  const double norm = full_norm(psi_ltr);
  CHECK(norm == doctest::Approx(ref_norm).epsilon(1.0e-12));

  ft n_op{tenes::real_tensor(mptensor::Shape(2, 2)),
          {{false, true}, {false, true}}};
  n_op.t.set_value(mptensor::Index(1, 1), 1.0);
  const double refs[4] = {ref_n0, ref_n1, ref_n2, ref_n3};
  for (int site = 0; site < 4; ++site) {
    ft npsi = apply_one_site_op(psi_ltr, site, n_op);
    const double value =
        tenes::fermion::trace(tenes::fermion::conj(psi_ltr), npsi,
                              mptensor::Axes(0, 1, 2, 3, 4, 5),
                              mptensor::Axes(0, 1, 2, 3, 4, 5)) /
        norm;
    CHECK(value == doctest::Approx(refs[site]).epsilon(1.0e-12));
  }

  ft hop{tenes::real_tensor(mptensor::Shape(2, 2, 2, 2)),
         {{false, true}, {false, true}, {false, true}, {false, true}}};
  hop.t.set_value(mptensor::Index(0, 1, 1, 0), 1.0);
  hop.t.set_value(mptensor::Index(1, 0, 0, 1), 1.0);
  ft hpsi = apply_two_site_op_01(psi_ltr, hop);
  const double hop_value =
      tenes::fermion::trace(tenes::fermion::conj(psi_ltr), hpsi,
                            mptensor::Axes(0, 1, 2, 3, 4, 5),
                            mptensor::Axes(0, 1, 2, 3, 4, 5)) /
      norm;
  CHECK(hop_value == doctest::Approx(ref_hop01).epsilon(1.0e-12));
}

TEST_CASE("manual four-swap reduced tensor matches f-primitive contraction") {
  tenes::fermion::leg_parities p{{false, true},
                                 {false, true},
                                 {false, true},
                                 {false, true},
                                 {false, true}};
  ft ket = make_even_ft(mptensor::Shape(2, 2, 2, 2, 2), p, 91);
  ft f_reduced = tenes::fermion::tensordot(
      tenes::fermion::conj(ket), ket, mptensor::Axes(4), mptensor::Axes(4));
  ft manual_ket = ket;
  tenes::fermion::apply_swap(manual_ket, 4, 0);
  tenes::fermion::apply_swap(manual_ket, 4, 1);
  tenes::fermion::apply_swap(manual_ket, 4, 2);
  tenes::fermion::apply_swap(manual_ket, 4, 3);
  tenes::real_tensor manual =
      mptensor::tensordot(tenes::fermion::conj(ket).t, manual_ket.t,
                          mptensor::Axes(4), mptensor::Axes(4));
  for (std::size_t n = 0; n < f_reduced.t.local_size(); ++n) {
    auto idx = f_reduced.t.global_index(n);
    double vf, vm;
    f_reduced.t.get_value(idx, vf);
    manual.get_value(idx, vm);
    CHECK(vf == doctest::Approx(vm).epsilon(1.0e-12));
  }
}

TEST_CASE("graded svd preserves parity for permuting axes") {
  // regression: svd with axes (0,2),(1,3) requires a non-identity internal
  // transpose; lazy-transpose + mptensor::reshape produced the storage-order
  // (wrong) unfolding. Legs deliberately carry DIFFERENT parity patterns so a
  // wrong unfolding cannot masquerade as a consistent graded SVD.
  namespace f = tenes::fermion;
  using FT = f::ftensor<tenes::real_tensor>;
  f::parity_vector p2{false, true};               // dim 2: {e,o}
  f::parity_vector p4{false, false, true, true};  // dim 4: {e,e,o,o}
  // leg layout matches the simple-update Theta: (aux1, aux2, out1, out2),
  // so rows (0,2) and cols (1,3) fuse legs with DIFFERENT parity patterns
  f::leg_parities lp{p4, p4, p2, p2};
  tenes::real_tensor t(mptensor::Shape(4, 4, 2, 2));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    auto idx = t.global_index(n);
    if (f::count_odd(lp, idx) % 2 == 0) {
      const double x = static_cast<double>(3 * (n + 2));
      t.set_value(idx, 0.31 * std::sin(x) + 0.23 * std::cos(0.6 * x));
    }
  }
  FT a{t, lp};
  FT u, vt;
  std::vector<double> s;
  f::svd(a, mptensor::Axes(0, 2), mptensor::Axes(1, 3), u, s, vt);
  CHECK(f::parity_violation(u) == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(f::parity_violation(vt) == doctest::Approx(0.0).epsilon(1e-12));
  FT us = u;
  us.multiply_vector(s, 2);
  FT rec = f::tensordot(us, vt, mptensor::Axes(2), mptensor::Axes(0));
  FT rec2 = f::transpose(rec, mptensor::Axes(0, 2, 1, 3));
  double maxdiff = 0.0;
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    auto idx = a.t.global_index(n);
    double va, vb;
    a.t.get_value(idx, va);
    rec2.t.get_value(idx, vb);
    maxdiff = std::max(maxdiff, std::abs(va - vb));
  }
  CHECK(maxdiff == doctest::Approx(0.0).epsilon(1e-10));
}

TEST_CASE("ftensor CTM environment preserves even parity on asymmetric legs") {
  namespace f = tenes::fermion;
  using FT = f::ftensor<tenes::real_tensor>;

  tenes::SquareLattice lattice(2, 2);
  tenes::itps::PEPS_Parameters params;
  params.CHI = 8;
  params.Max_CTM_Iteration = 2;
  params.CTM_Convergence_Epsilon = 1.0e-12;
  params.Use_RSVD = false;

  const f::parity_vector phys{false, true};
  const f::parity_vector p_h{false, true};
  const f::parity_vector p_v{false, false, true, true};
  std::vector<FT> Tn;
  Tn.reserve(lattice.N_UNIT);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    f::leg_parities parity{p_h, p_v, p_h, p_v, phys};
    FT t{tenes::real_tensor(mptensor::Shape(2, 4, 2, 4, 2)), parity};
    for (std::size_t n = 0; n < t.t.local_size(); ++n) {
      auto idx = t.t.global_index(n);
      if (f::count_odd(parity, idx) % 2 == 0) {
        const double x = static_cast<double>((site + 1) * (n + 3));
        t.t.set_value(idx, 0.15 * std::sin(x) + 0.07 * std::cos(0.5 * x));
      }
    }
    Tn.push_back(t);
  }

  std::vector<FT> C1(lattice.N_UNIT), C2(lattice.N_UNIT), C3(lattice.N_UNIT),
      C4(lattice.N_UNIT), eTt(lattice.N_UNIT), eTr(lattice.N_UNIT),
      eTb(lattice.N_UNIT), eTl(lattice.N_UNIT);
  tenes::fermion::parity_cleanup_observations().clear();
  int iterations = tenes::itps::core::Calc_CTM_Environment(
      C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, params, lattice);
  double max_cleanup = 0.0;
  for (double v : tenes::fermion::parity_cleanup_observations()) {
    max_cleanup = std::max(max_cleanup, v);
  }
  std::cout << "ctm parity cleanup max=" << max_cleanup << std::endl;

  CHECK(iterations > 0);
  CHECK(max_cleanup < 1.0e-12);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    CHECK(f::parity_violation(C1[site]) == doctest::Approx(0.0));
    CHECK(f::parity_violation(C2[site]) == doctest::Approx(0.0));
    CHECK(f::parity_violation(C3[site]) == doctest::Approx(0.0));
    CHECK(f::parity_violation(C4[site]) == doctest::Approx(0.0));
    CHECK(f::parity_violation(eTt[site]) == doctest::Approx(0.0));
    CHECK(f::parity_violation(eTr[site]) == doctest::Approx(0.0));
    CHECK(f::parity_violation(eTb[site]) == doctest::Approx(0.0));
    CHECK(f::parity_violation(eTl[site]) == doctest::Approx(0.0));
  }
}

TEST_CASE("ftensor CTM measurement kernels produce finite positive norm") {
  namespace f = tenes::fermion;
  using FT = f::ftensor<tenes::real_tensor>;

  tenes::SquareLattice lattice(2, 2);
  tenes::itps::PEPS_Parameters params;
  params.CHI = 8;
  params.Max_CTM_Iteration = 2;
  params.CTM_Convergence_Epsilon = 1.0e-12;
  params.Use_RSVD = false;

  const f::parity_vector phys{false, true};
  const f::parity_vector p_h{false, true};
  const f::parity_vector p_v{false, false, true, true};
  std::vector<FT> Tn;
  Tn.reserve(lattice.N_UNIT);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    f::leg_parities parity{p_h, p_v, p_h, p_v, phys};
    FT t{tenes::real_tensor(mptensor::Shape(2, 4, 2, 4, 2)), parity};
    for (std::size_t n = 0; n < t.t.local_size(); ++n) {
      auto idx = t.t.global_index(n);
      if (f::count_odd(parity, idx) % 2 == 0) {
        const double x = static_cast<double>((site + 2) * (n + 5));
        t.t.set_value(idx, 0.11 * std::sin(x) + 0.09 * std::cos(0.4 * x));
      }
    }
    Tn.push_back(t);
  }

  std::vector<FT> C1(lattice.N_UNIT), C2(lattice.N_UNIT), C3(lattice.N_UNIT),
      C4(lattice.N_UNIT), eTt(lattice.N_UNIT), eTr(lattice.N_UNIT),
      eTb(lattice.N_UNIT), eTl(lattice.N_UNIT);
  tenes::fermion::parity_cleanup_observations().clear();
  tenes::itps::core::Calc_CTM_Environment(C1, C2, C3, C4, eTt, eTr, eTb, eTl,
                                          Tn, params, lattice);
  double max_cleanup = 0.0;
  for (double v : tenes::fermion::parity_cleanup_observations()) {
    max_cleanup = std::max(max_cleanup, v);
  }
  std::cout << "ctm measurement parity cleanup max=" << max_cleanup
            << std::endl;
  CHECK(max_cleanup < 1.0e-12);

  FT id{tenes::real_tensor(mptensor::Shape(2, 2)), {phys, phys}};
  id.t.set_value(mptensor::Index(0, 0), 1.0);
  id.t.set_value(mptensor::Index(1, 1), 1.0);
  const double norm = tenes::itps::core::Contract_one_site_iTPS_CTM(
      C1[0], C2[0], C3[0], C4[0], eTt[0], eTr[0], eTb[0], eTl[0], Tn[0], id);
  CHECK(std::isfinite(norm));
  CHECK(norm > 0.0);

  FT hop{tenes::real_tensor(mptensor::Shape(2, 2, 2, 2)),
         {phys, phys, phys, phys}};
  hop.t.set_value(mptensor::Index(0, 1, 1, 0), 1.0);
  hop.t.set_value(mptensor::Index(1, 0, 0, 1), 1.0);
  const int left = 0;
  const int right = lattice.right(left);
  const double hop_value =
      tenes::itps::core::Contract_two_sites_horizontal_op12_iTPS_CTM(
          C1[left], C2[right], C3[right], C4[left], eTt[left], eTt[right],
          eTr[right], eTb[right], eTb[left], eTl[left], Tn[left], Tn[right],
          hop);
  CHECK(std::isfinite(hop_value));
}

#include "fermion/r2_convention.cpp"
