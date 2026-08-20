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

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>

#include "../src/fermion/fermion_info.hpp"
#include "../src/fermion/fops.hpp"
#include "../src/fermion/ftensor.hpp"
#include "../src/fermion/parity.hpp"
#include "../src/fermion/reduced.hpp"
#include "../src/fermion/reduced_measure.hpp"
#include "../src/SquareLattice.hpp"
#include "../src/tensor.hpp"
#include "../src/iTPS/PEPS_Parameters.hpp"
#include "../src/iTPS/iTPS.hpp"
#include "../src/iTPS/core/ctm.hpp"
#include "../src/iTPS/core/contract.hpp"
#include "../src/iTPS/core/contract_itps_ctm.hpp"
#include "../src/iTPS/core/simple_update.hpp"

using namespace tenes::fermion;

using ft = tenes::fermion::ftensor<tenes::real_tensor>;

namespace tenes::itps {
struct iTPSTestAccessor {
  template <class tensor>
  static std::vector<tensor>& Tn(iTPS<tensor>& state) {
    return state.Tn;
  }

  template <class tensor>
  static std::vector<std::vector<std::vector<double>>>& lambda_tensor(
      iTPS<tensor>& state) {
    return state.lambda_tensor;
  }

  template <class tensor>
  static tenes::fermion::FermionInfo& finfo(iTPS<tensor>& state) {
    return state.finfo;
  }

  template <class tensor>
  static std::vector<std::string>& twosite_operator_names(iTPS<tensor>& state) {
    return state.twosite_operator_names;
  }

  template <class tensor>
  static void update_reduced_density_environment(iTPS<tensor>& state) {
    // Bare Tn, matching measure.cpp: the CTM provides the environment.
    const std::vector<tensor> reduced =
        tenes::fermion::build_reduced_density_tensors(state.Tn, state.finfo);
    core::Calc_CTM_Environment_density(
        state.C1, state.C2, state.C3, state.C4, state.eTt, state.eTr, state.eTb,
        state.eTl, reduced, state.peps_parameters, state.lattice);
  }
};
}  // namespace tenes::itps

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

TEST_CASE(
    "two-site operator doubly-odd input channel matches one-site composition") {
  // The hopping operator pinned in the JW reference test has no matrix
  // element whose two INPUT legs are both odd, so that test cannot fix the
  // sign convention of e.g. the |11><11| channel present in every Trotter
  // gate exp(tau h). Pin it here with n0*n1, whose only nonzero element is
  // exactly that channel: the two-site graded application must agree with
  // the (already validated) composition of one-site applications.
  ft psi = reference_chain_left_to_right();
  const double norm = full_norm(psi);

  ft n_op{tenes::real_tensor(mptensor::Shape(2, 2)),
          {{false, true}, {false, true}}};
  n_op.t.set_value(mptensor::Index(1, 1), 1.0);
  ft composed = apply_one_site_op(apply_one_site_op(psi, 1, n_op), 0, n_op);
  const double ref = tenes::fermion::trace(tenes::fermion::conj(psi), composed,
                                           mptensor::Axes(0, 1, 2, 3, 4, 5),
                                           mptensor::Axes(0, 1, 2, 3, 4, 5)) /
                     norm;
  CHECK(ref > 0.0);

  tenes::real_tensor nn_plain(mptensor::Shape(2, 2, 2, 2));
  nn_plain.set_value(mptensor::Index(1, 1, 1, 1), 1.0);
  ft nn = tenes::fermion::wrap_twosite_op(nn_plain, parity_vector{false, true},
                                          parity_vector{false, true});
  double channel = 0.0;
  nn.t.get_value(mptensor::Index(1, 1, 1, 1), channel);
  CHECK(channel == doctest::Approx(-1.0));
  ft nnpsi = apply_two_site_op_01(psi, nn);
  const double value = tenes::fermion::trace(tenes::fermion::conj(psi), nnpsi,
                                             mptensor::Axes(0, 1, 2, 3, 4, 5),
                                             mptensor::Axes(0, 1, 2, 3, 4, 5)) /
                       norm;
  CHECK(value == doctest::Approx(ref).epsilon(1.0e-12));
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

TEST_CASE("graded svd preserves parity for asymmetric non-crossing axes") {
  // Legs deliberately carry DIFFERENT parity patterns so a wrong fused-index
  // convention cannot masquerade as a consistent graded SVD.
  namespace f = tenes::fermion;
  using FT = f::ftensor<tenes::real_tensor>;
  f::parity_vector p2{false, true};               // dim 2: {e,o}
  f::parity_vector p4{false, false, true, true};  // dim 4: {e,e,o,o}
  // leg layout matches the fixed simple-update Theta:
  // (aux1, out1, out2, aux2), rows (0,1) and cols (2,3).
  f::leg_parities lp{p4, p2, p2, p4};
  tenes::real_tensor t(mptensor::Shape(4, 2, 2, 4));
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
  f::svd(a, mptensor::Axes(0, 1), mptensor::Axes(2, 3), u, s, vt);
  CHECK(f::parity_violation(u) == doctest::Approx(0.0).epsilon(1e-12));
  CHECK(f::parity_violation(vt) == doctest::Approx(0.0).epsilon(1e-12));
  FT us = u;
  us.multiply_vector(s, 2);
  FT rec = f::tensordot(us, vt, mptensor::Axes(2), mptensor::Axes(0));
  double maxdiff = 0.0;
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    auto idx = a.t.global_index(n);
    double va, vb;
    a.t.get_value(idx, va);
    rec.t.get_value(idx, vb);
    maxdiff = std::max(maxdiff, std::abs(va - vb));
  }
  CHECK(maxdiff == doctest::Approx(0.0).epsilon(1e-10));
}

TEST_CASE("graded svd is invariant under regrouping to a contiguous split") {
  // svd(a, rows, cols) starts by transposing a into (rows..., cols...) order,
  // so requesting an interleaved split is exactly the same computation as
  // first applying that graded transpose and then splitting contiguously.
  // The kernel therefore hands svd_trunc the interleaved axes directly. The
  // third variant below - a metadata-only regroup that drops the Koszul mask
  // - is what a past version of the kernel did, and it does NOT agree: that
  // dropped sign (cancelling against a missing gate sign) was the actual bug
  // once misattributed to interleaved splits.
  namespace f = tenes::fermion;
  f::leg_parities p{
      {false, true}, {false, true, true}, {false, true}, {false, true, true}};
  ft a = make_even_ft(mptensor::Shape(2, 3, 2, 3), p, 4321);

  ft u1, vt1, u2, vt2;
  std::vector<double> s1, s2;
  f::svd(f::transpose(a, mptensor::Axes(0, 2, 1, 3)), mptensor::Axes(0, 1),
         mptensor::Axes(2, 3), u1, s1, vt1);
  f::svd(a, mptensor::Axes(0, 2), mptensor::Axes(1, 3), u2, s2, vt2);
  REQUIRE(s1.size() == s2.size());
  for (std::size_t i = 0; i < s1.size(); ++i) {
    CHECK(s1[i] == doctest::Approx(s2[i]).epsilon(1.0e-12));
  }

  ft maskless;
  maskless.t = a.t;
  maskless.t.transpose(mptensor::Axes(0, 2, 1, 3));
  maskless.parity = {a.parity[0], a.parity[2], a.parity[1], a.parity[3]};
  ft u3, vt3;
  std::vector<double> s3;
  f::svd(maskless, mptensor::Axes(0, 1), mptensor::Axes(2, 3), u3, s3, vt3);
  double maskless_diff = 0.0;
  for (std::size_t i = 0; i < s1.size(); ++i) {
    maskless_diff = std::max(maskless_diff, std::abs(s1[i] - s3[i]));
  }
  CHECK(maskless_diff > 1.0e-3);
}

TEST_CASE("reduced density CTM measurement gives positive fermion norms") {
  namespace f = tenes::fermion;

  tenes::SquareLattice lattice(2, 2);
  tenes::itps::PEPS_Parameters params;
  params.CHI = 8;
  params.Max_CTM_Iteration = 2;
  params.CTM_Convergence_Epsilon = 1.0e-12;
  params.Use_RSVD = false;

  const f::parity_vector phys{false, true};
  const f::parity_vector p_h{false, true};
  const f::parity_vector p_v{false, true};

  f::FermionInfo finfo;
  finfo.enabled = true;
  finfo.phys.assign(lattice.N_UNIT, phys);
  finfo.virt.assign(lattice.N_UNIT,
                    std::array<f::parity_vector, 4>{p_h, p_v, p_h, p_v});

  std::vector<tenes::real_tensor> Tn;
  Tn.reserve(lattice.N_UNIT);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    f::leg_parities parity = f::Tn_parity(finfo, site);
    tenes::real_tensor t(mptensor::Shape(2, 2, 2, 2, 2));
    for (std::size_t n = 0; n < t.local_size(); ++n) {
      auto idx = t.global_index(n);
      if (f::count_odd(parity, idx) % 2 == 0) {
        const double x = static_cast<double>((site + 3) * (n + 7));
        t.set_value(idx, 0.12 * std::sin(x) + 0.08 * std::cos(0.7 * x));
      }
    }
    Tn.push_back(t);
  }

  std::vector<tenes::real_tensor> reduced =
      f::build_reduced_density_tensors(Tn, finfo);
  std::vector<tenes::real_tensor> C1(lattice.N_UNIT), C2(lattice.N_UNIT),
      C3(lattice.N_UNIT), C4(lattice.N_UNIT), eTt(lattice.N_UNIT),
      eTr(lattice.N_UNIT), eTb(lattice.N_UNIT), eTl(lattice.N_UNIT);
  tenes::itps::core::Calc_CTM_Environment_density(
      C1, C2, C3, C4, eTt, eTr, eTb, eTl, reduced, params, lattice);

  tenes::real_tensor id(mptensor::Shape(phys.size(), phys.size()));
  tenes::real_tensor number(mptensor::Shape(phys.size(), phys.size()));
  id.set_value(mptensor::Index(0, 0), 1.0);
  id.set_value(mptensor::Index(1, 1), 1.0);
  number.set_value(mptensor::Index(1, 1), 1.0);

  std::vector<double> density(lattice.N_UNIT);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    const double norm = tenes::itps::core::Contract_one_site_density_CTM(
        C1[site], C2[site], C3[site], C4[site], eTt[site], eTr[site], eTb[site],
        eTl[site], reduced[site], id);
    const double n =
        tenes::itps::core::Contract_one_site_density_CTM(
            C1[site], C2[site], C3[site], C4[site], eTt[site], eTr[site],
            eTb[site], eTl[site], reduced[site], number) /
        norm;
    density[site] = n;
    CHECK(std::isfinite(norm));
    CHECK(norm > 0.0);
    CHECK(n >= doctest::Approx(0.0).epsilon(1.0e-12));
    CHECK(n <= doctest::Approx(1.0).epsilon(1.0e-12));
  }

  tenes::real_tensor id2(
      mptensor::Shape(phys.size(), phys.size(), phys.size(), phys.size()));
  tenes::real_tensor density2(
      mptensor::Shape(phys.size(), phys.size(), phys.size(), phys.size()));
  for (std::size_t in_a = 0; in_a < phys.size(); ++in_a) {
    for (std::size_t in_b = 0; in_b < phys.size(); ++in_b) {
      id2.set_value(mptensor::Index(in_a, in_b, in_a, in_b), 1.0);
      density2.set_value(mptensor::Index(in_a, in_b, in_a, in_b),
                         static_cast<double>(in_a + in_b));
    }
  }
  auto check_pair_density = [&](int first, int second,
                                f::reduced_pair_direction direction) {
    f::ftensor<tenes::real_tensor> fid{id2, {phys, phys, phys, phys}};
    f::ftensor<tenes::real_tensor> fdensity{density2, {phys, phys, phys, phys}};
    const tenes::real_tensor norm_blob = f::build_reduced_identity_pair(
        f::wrap_Tn(Tn[first], finfo, first),
        f::wrap_Tn(Tn[second], finfo, second), direction);
    const tenes::real_tensor density_blob = f::build_reduced_pair(
        f::wrap_Tn(Tn[first], finfo, first),
        f::wrap_Tn(Tn[second], finfo, second), fdensity, direction);
    double norm = 0.0;
    double val = 0.0;
    double ref_norm = 0.0;
    if (direction == f::reduced_pair_direction::horizontal) {
      norm = f::contract_reduced_pair_horizontal_density_CTM(
          C1[first], C2[second], C3[second], C4[first], eTt[first], eTt[second],
          eTr[second], eTb[second], eTb[first], eTl[first], norm_blob);
      val = f::contract_reduced_pair_horizontal_density_CTM(
          C1[first], C2[second], C3[second], C4[first], eTt[first], eTt[second],
          eTr[second], eTb[second], eTb[first], eTl[first], density_blob);
      ref_norm =
          tenes::itps::core::Contract_two_sites_horizontal_op12_density_CTM(
              C1[first], C2[second], C3[second], C4[first], eTt[first],
              eTt[second], eTr[second], eTb[second], eTb[first], eTl[first],
              reduced[first], reduced[second], id2);
    } else {
      norm = f::contract_reduced_pair_vertical_density_CTM(
          C1[first], C2[first], C3[second], C4[second], eTt[first], eTr[first],
          eTr[second], eTb[second], eTl[second], eTl[first], norm_blob);
      val = f::contract_reduced_pair_vertical_density_CTM(
          C1[first], C2[first], C3[second], C4[second], eTt[first], eTr[first],
          eTr[second], eTb[second], eTl[second], eTl[first], density_blob);
      ref_norm =
          tenes::itps::core::Contract_two_sites_vertical_op12_density_CTM(
              C1[first], C2[first], C3[second], C4[second], eTt[first],
              eTr[first], eTr[second], eTb[second], eTl[second], eTl[first],
              reduced[first], reduced[second], id2);
    }
    CHECK(norm > 0.0);
    CHECK(norm == doctest::Approx(ref_norm).epsilon(1.0e-12));
    CHECK(std::isfinite(val / norm));
    CHECK(val / norm >= doctest::Approx(0.0).epsilon(1.0e-12));
    CHECK(val / norm <= doctest::Approx(2.0).epsilon(1.0e-12));
  };
  check_pair_density(0, lattice.right(0),
                     f::reduced_pair_direction::horizontal);
  check_pair_density(0, lattice.bottom(0), f::reduced_pair_direction::vertical);
}

static tenes::real_tensor make_free_fermion_gate(double tau) {
  tenes::real_tensor gate(mptensor::Shape(2, 2, 2, 2));
  gate.set_value(mptensor::Index(0, 0, 0, 0), 1.0);
  gate.set_value(mptensor::Index(1, 1, 1, 1), 1.0);
  gate.set_value(mptensor::Index(0, 1, 0, 1), std::cosh(tau));
  gate.set_value(mptensor::Index(1, 0, 1, 0), std::cosh(tau));
  gate.set_value(mptensor::Index(0, 1, 1, 0), std::sinh(tau));
  gate.set_value(mptensor::Index(1, 0, 0, 1), std::sinh(tau));
  return gate;
}

struct TensorDiffReport {
  double max_abs = 0.0;
  int site = -1;
  mptensor::Index index;
  double fermion_value = 0.0;
  double boson_value = 0.0;
};

static TensorDiffReport max_tn_diff(
    const std::vector<tenes::real_tensor>& fermion_tensors,
    const std::vector<tenes::real_tensor>& boson_tensors) {
  TensorDiffReport report;
  for (std::size_t site = 0; site < fermion_tensors.size(); ++site) {
    const auto& fT = fermion_tensors[site];
    const auto& bT = boson_tensors[site];
    for (std::size_t n = 0; n < fT.local_size(); ++n) {
      const auto idx = fT.global_index(n);
      double fv = 0.0;
      double bv = 0.0;
      fT.get_value(idx, fv);
      bT.get_value(idx, bv);
      const double diff = std::abs(fv - bv);
      if (diff > report.max_abs) {
        report.max_abs = diff;
        report.site = static_cast<int>(site);
        report.index = idx;
        report.fermion_value = fv;
        report.boson_value = bv;
      }
    }
  }
  return report;
}

static std::string index_to_string(const mptensor::Index& idx) {
  std::ostringstream os;
  os << "(";
  for (std::size_t i = 0; i < idx.size(); ++i) {
    if (i != 0) {
      os << ",";
    }
    os << idx[i];
  }
  os << ")";
  return os.str();
}

static std::vector<double> sorted_desc(std::vector<double> values) {
  std::sort(values.begin(), values.end(), std::greater<double>());
  return values;
}

static double lambda_relative_diff(const std::vector<double>& a,
                                   const std::vector<double>& b) {
  double diff = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    const double scale = std::max({1.0e-300, std::abs(a[i]), std::abs(b[i])});
    diff = std::max(diff, std::abs(a[i] - b[i]) / scale);
  }
  return diff;
}

static std::string vector_to_string(const std::vector<double>& values) {
  std::ostringstream os;
  os << "[";
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i != 0) {
      os << ",";
    }
    os << std::setprecision(17) << values[i];
  }
  os << "]";
  return os.str();
}

// Reorder theta from the kernel's (aux1, aux2, out1, out2) layout to the
// site-contiguous (aux1, out1, aux2, out2) one, so the diagnostics below can
// report row/col sector data with a plain (0,1) x (2,3) split. Equivalent to
// the interleaved split the kernel hands to svd_trunc.
template <class Tensor>
static Tensor diagnostic_site_ordered_theta(const Tensor& theta) {
  Tensor ret = theta;
  ret.transpose(mptensor::Axes(0, 2, 1, 3));
  return ret;
}

static tenes::fermion::ftensor<tenes::real_tensor>
diagnostic_site_ordered_theta(
    const tenes::fermion::ftensor<tenes::real_tensor>& theta) {
  return tenes::fermion::transpose(theta, mptensor::Axes(0, 2, 1, 3));
}

template <class Tensor>
static std::pair<Tensor, Tensor> diagnostic_update_thetas(
    const Tensor& Tn1, const Tensor& Tn2,
    const std::vector<std::vector<double>>& lambda1,
    const std::vector<std::vector<double>>& lambda2, const Tensor& op12,
    int connect1) {
  const int connect2 = (connect1 + 2) % 4;
  Tensor Tn1_lambda = Tn1;
  Tensor Tn2_lambda = Tn2;

  if (connect1 == 0) {
    Tn1_lambda.multiply_vector(lambda1[1], 1, lambda1[2], 2, lambda1[3], 3);
    Tn1_lambda.transpose(mptensor::Axes(1, 2, 3, 0, 4));
  } else if (connect1 == 1) {
    Tn1_lambda.multiply_vector(lambda1[0], 0, lambda1[2], 2, lambda1[3], 3);
    Tn1_lambda.transpose(mptensor::Axes(0, 2, 3, 1, 4));
  } else if (connect1 == 2) {
    Tn1_lambda.multiply_vector(lambda1[0], 0, lambda1[1], 1, lambda1[3], 3);
    Tn1_lambda.transpose(mptensor::Axes(0, 1, 3, 2, 4));
  } else {
    Tn1_lambda.multiply_vector(lambda1[0], 0, lambda1[1], 1, lambda1[2], 2);
  }

  if (connect2 == 0) {
    Tn2_lambda.multiply_vector(lambda2[1], 1, lambda2[2], 2, lambda2[3], 3);
    Tn2_lambda.transpose(mptensor::Axes(1, 2, 3, 0, 4));
  } else if (connect2 == 1) {
    Tn2_lambda.multiply_vector(lambda2[0], 0, lambda2[2], 2, lambda2[3], 3);
    Tn2_lambda.transpose(mptensor::Axes(0, 2, 3, 1, 4));
  } else if (connect2 == 2) {
    Tn2_lambda.multiply_vector(lambda2[0], 0, lambda2[1], 1, lambda2[3], 3);
    Tn2_lambda.transpose(mptensor::Axes(0, 1, 3, 2, 4));
  } else {
    Tn2_lambda.multiply_vector(lambda2[0], 0, lambda2[1], 1, lambda2[2], 2);
  }

  Tensor Q1, R1, Q2, R2;
  qr(Tn1_lambda, mptensor::Axes(0, 1, 2), mptensor::Axes(3, 4), Q1, R1);
  qr(Tn2_lambda, mptensor::Axes(0, 1, 2), mptensor::Axes(3, 4), Q2, R2);
  Tensor theta_before = tensordot(R1, R2, mptensor::Axes(1), mptensor::Axes(1));
  Tensor theta = diagnostic_site_ordered_theta(tensordot(
      theta_before, op12, mptensor::Axes(1, 3), mptensor::Axes(0, 1)));
  return {theta_before, theta};
}

template <class Tensor>
struct DiagnosticGateData {
  Tensor R1;
  Tensor R2;
  Tensor theta_before;
  Tensor theta_after;
};

template <class Tensor>
static DiagnosticGateData<Tensor> diagnostic_gate_data(
    const Tensor& Tn1, const Tensor& Tn2,
    const std::vector<std::vector<double>>& lambda1,
    const std::vector<std::vector<double>>& lambda2, const Tensor& op12,
    int connect1) {
  const int connect2 = (connect1 + 2) % 4;
  Tensor Tn1_lambda = Tn1;
  Tensor Tn2_lambda = Tn2;

  if (connect1 == 0) {
    Tn1_lambda.multiply_vector(lambda1[1], 1, lambda1[2], 2, lambda1[3], 3);
    Tn1_lambda.transpose(mptensor::Axes(1, 2, 3, 0, 4));
  } else if (connect1 == 1) {
    Tn1_lambda.multiply_vector(lambda1[0], 0, lambda1[2], 2, lambda1[3], 3);
    Tn1_lambda.transpose(mptensor::Axes(0, 2, 3, 1, 4));
  } else if (connect1 == 2) {
    Tn1_lambda.multiply_vector(lambda1[0], 0, lambda1[1], 1, lambda1[3], 3);
    Tn1_lambda.transpose(mptensor::Axes(0, 1, 3, 2, 4));
  } else {
    Tn1_lambda.multiply_vector(lambda1[0], 0, lambda1[1], 1, lambda1[2], 2);
  }

  if (connect2 == 0) {
    Tn2_lambda.multiply_vector(lambda2[1], 1, lambda2[2], 2, lambda2[3], 3);
    Tn2_lambda.transpose(mptensor::Axes(1, 2, 3, 0, 4));
  } else if (connect2 == 1) {
    Tn2_lambda.multiply_vector(lambda2[0], 0, lambda2[2], 2, lambda2[3], 3);
    Tn2_lambda.transpose(mptensor::Axes(0, 2, 3, 1, 4));
  } else if (connect2 == 2) {
    Tn2_lambda.multiply_vector(lambda2[0], 0, lambda2[1], 1, lambda2[3], 3);
    Tn2_lambda.transpose(mptensor::Axes(0, 1, 3, 2, 4));
  } else {
    Tn2_lambda.multiply_vector(lambda2[0], 0, lambda2[1], 1, lambda2[2], 2);
  }

  DiagnosticGateData<Tensor> data;
  Tensor Q1, Q2;
  qr(Tn1_lambda, mptensor::Axes(0, 1, 2), mptensor::Axes(3, 4), Q1, data.R1);
  qr(Tn2_lambda, mptensor::Axes(0, 1, 2), mptensor::Axes(3, 4), Q2, data.R2);
  data.theta_before =
      tensordot(data.R1, data.R2, mptensor::Axes(1), mptensor::Axes(1));
  data.theta_after = diagnostic_site_ordered_theta(tensordot(
      data.theta_before, op12, mptensor::Axes(1, 3), mptensor::Axes(0, 1)));
  return data;
}

static std::string parity_vector_to_string(
    const tenes::fermion::parity_vector& parity) {
  std::ostringstream os;
  os << "[";
  for (std::size_t i = 0; i < parity.size(); ++i) {
    if (i != 0) {
      os << ",";
    }
    os << (parity[i] ? 1 : 0);
  }
  os << "]";
  return os.str();
}

static const char* parity_name(bool odd) { return odd ? "odd" : "even"; }

static bool combined_parity(const tenes::fermion::leg_parities& parity,
                            const mptensor::Index& idx,
                            const mptensor::Axes& axes) {
  bool odd = false;
  for (std::size_t i = 0; i < axes.size(); ++i) {
    odd = odd != parity[axes[i]][idx[axes[i]]];
  }
  return odd;
}

static void dump_parity_metadata(
    const std::string& label,
    const tenes::fermion::ftensor<tenes::real_tensor>& r1,
    const tenes::fermion::ftensor<tenes::real_tensor>& r2,
    const tenes::fermion::ftensor<tenes::real_tensor>& theta_before) {
  for (std::size_t ax = 0; ax < r1.parity.size(); ++ax) {
    std::cout << label << " R1.parity[" << ax
              << "]=" << parity_vector_to_string(r1.parity[ax]) << std::endl;
  }
  for (std::size_t ax = 0; ax < r2.parity.size(); ++ax) {
    std::cout << label << " R2.parity[" << ax
              << "]=" << parity_vector_to_string(r2.parity[ax]) << std::endl;
  }
  for (std::size_t ax = 0; ax < theta_before.parity.size(); ++ax) {
    std::cout << label << " theta_before.parity[" << ax
              << "]=" << parity_vector_to_string(theta_before.parity[ax])
              << std::endl;
  }
}

static void dump_tensordot_mask_pattern(
    const std::string& label,
    const tenes::fermion::ftensor<tenes::real_tensor>& theta_before,
    const tenes::fermion::ftensor<tenes::real_tensor>& op12) {
  const mptensor::Axes axes_a(1, 3);
  const mptensor::Axes axes_b(0, 1);
  const mptensor::Axes left_perm = tenes::fermion::detail::tensordot_left_perm(
      theta_before.parity.size(), axes_a);
  const mptensor::Axes right_perm =
      tenes::fermion::detail::tensordot_right_perm(op12.parity.size(), axes_b);

  std::cout << label << " left_perm=" << index_to_string(left_perm)
            << " right_perm=" << index_to_string(right_perm) << std::endl;

  std::size_t left_count = 0;
  for (std::size_t n = 0; n < theta_before.t.local_size(); ++n) {
    const auto idx = theta_before.t.global_index(n);
    if (tenes::fermion::detail::transpose_sign(theta_before.parity, idx,
                                               left_perm) < 0) {
      ++left_count;
      const bool row =
          combined_parity(theta_before.parity, idx, mptensor::Axes(0, 1));
      const bool col =
          combined_parity(theta_before.parity, idx, mptensor::Axes(2, 3));
      std::cout << label << " left_mask_minus idx=" << index_to_string(idx)
                << " sector=" << parity_name(row) << "/" << parity_name(col)
                << std::endl;
    }
  }
  std::size_t right_count = 0;
  for (std::size_t n = 0; n < op12.t.local_size(); ++n) {
    const auto idx = op12.t.global_index(n);
    if (tenes::fermion::detail::transpose_sign(op12.parity, idx, right_perm) <
        0) {
      ++right_count;
      std::cout << label << " right_mask_minus idx=" << index_to_string(idx)
                << " in_sector="
                << parity_name(
                       combined_parity(op12.parity, idx, mptensor::Axes(0, 1)))
                << " out_sector="
                << parity_name(
                       combined_parity(op12.parity, idx, mptensor::Axes(2, 3)))
                << std::endl;
    }
  }
  std::cout << label << " mask_minus_counts left=" << left_count
            << " right=" << right_count << std::endl;
}

static void dump_post_gate_diff(
    const std::string& label,
    const tenes::fermion::ftensor<tenes::real_tensor>& ftheta,
    const tenes::real_tensor& btheta) {
  std::size_t diff_count = 0;
  double max_abs_diff = 0.0;
  std::array<std::size_t, 4> sector_counts{0, 0, 0, 0};
  for (std::size_t n = 0; n < ftheta.t.local_size(); ++n) {
    const auto idx = ftheta.t.global_index(n);
    double fv = 0.0;
    double bv = 0.0;
    ftheta.t.get_value(idx, fv);
    btheta.get_value(idx, bv);
    const double scale = std::max({1.0, std::abs(fv), std::abs(bv)});
    const double diff = std::abs(fv - bv);
    if (diff <= 1.0e-12 * scale) {
      continue;
    }
    const bool row = combined_parity(ftheta.parity, idx, mptensor::Axes(0, 1));
    const bool col = combined_parity(ftheta.parity, idx, mptensor::Axes(2, 3));
    const std::size_t sector = (row ? 2 : 0) + (col ? 1 : 0);
    ++sector_counts[sector];
    ++diff_count;
    max_abs_diff = std::max(max_abs_diff, diff);
    std::cout << std::setprecision(17) << label
              << " post_gate_diff idx=" << index_to_string(idx)
              << " sector=" << parity_name(row) << "/" << parity_name(col)
              << " fermion=" << fv << " boson=" << bv << " diff=" << diff
              << " ratio=" << (std::abs(bv) > 0.0 ? fv / bv : 0.0) << std::endl;
  }
  std::cout << label << " post_gate_diff_summary count=" << diff_count
            << " max_abs_diff=" << max_abs_diff
            << " even/even=" << sector_counts[0]
            << " even/odd=" << sector_counts[1]
            << " odd/even=" << sector_counts[2]
            << " odd/odd=" << sector_counts[3] << std::endl;
}

struct SectorSpectrum {
  std::vector<double> even;
  std::vector<double> odd;
};

static SectorSpectrum sector_spectrum(
    const tenes::fermion::ftensor<tenes::real_tensor>& theta,
    const mptensor::Axes& rows, const mptensor::Axes& cols) {
  tenes::fermion::ftensor<tenes::real_tensor> u, vt;
  std::vector<double> s;
  tenes::fermion::svd(theta, rows, cols, u, s, vt);
  SectorSpectrum out;
  const auto& internal = u.parity.back();
  for (std::size_t i = 0; i < s.size(); ++i) {
    (internal[i] ? out.odd : out.even).push_back(s[i]);
  }
  return out;
}

static double off_block_max(
    const tenes::fermion::ftensor<tenes::real_tensor>& theta,
    const mptensor::Axes& rows, const mptensor::Axes& cols) {
  double max_value = 0.0;
  for (std::size_t n = 0; n < theta.t.local_size(); ++n) {
    const auto idx = theta.t.global_index(n);
    bool row_odd = false;
    for (std::size_t i = 0; i < rows.size(); ++i) {
      row_odd = row_odd != theta.parity[rows[i]][idx[rows[i]]];
    }
    bool col_odd = false;
    for (std::size_t i = 0; i < cols.size(); ++i) {
      col_odd = col_odd != theta.parity[cols[i]][idx[cols[i]]];
    }
    if (row_odd == col_odd) {
      continue;
    }
    double value = 0.0;
    theta.t.get_value(idx, value);
    max_value = std::max(max_value, std::abs(value));
  }
  return max_value;
}

static std::vector<double> normalized_lambda_from_s(std::vector<double> s) {
  double norm2 = 0.0;
  for (const double value : s) {
    norm2 += value * value;
  }
  const double norm = std::sqrt(norm2);
  for (double& value : s) {
    value = std::sqrt(value / norm);
  }
  return s;
}

TEST_CASE("diagnostic horizontal-chain fermion/boson trajectory divergence") {
  if (std::getenv("TENES_RUN_HORIZONTAL_TRAJECTORY_DIAG") == nullptr) {
    return;
  }

  namespace f = tenes::fermion;
  using tensor = tenes::real_tensor;

  auto make_lattice = []() {
    tenes::SquareLattice lattice(2, 2);
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      lattice.physical_dims[site] = 2;
      lattice.virtual_dims[site] = {2, 1, 2, 1};
      lattice.initial_dirs[site] = {0.0};
      lattice.noises[site] = 1.0;
    }
    return lattice;
  };
  auto make_params = [](bool fermion, const std::string& outdir) {
    tenes::itps::PEPS_Parameters params;
    params.fermion = fermion;
    if (fermion) {
      params.phys_parity.assign(4, f::parity_vector{false, true});
    }
    params.print_level = tenes::PrintLevel::none;
    params.outdir = outdir;
    params.CHI = 8;
    params.Max_CTM_Iteration = 10;
    params.CTM_Convergence_Epsilon = 1.0e-8;
    params.Use_RSVD = false;
    params.seed = 11;
    return params;
  };
  auto make_updates = []() {
    std::vector<tenes::EvolutionOperator<tensor>> updates;
    const tensor gate = make_free_fermion_gate(0.01);
    for (int site = 0; site < 4; ++site) {
      updates.push_back(
          tenes::make_twosite_EvolutionOperator(site, 2, 0, gate));
    }
    return updates;
  };

  auto lattice_f = make_lattice();
  auto lattice_b = make_lattice();
  const auto updates = make_updates();
  auto fparams = make_params(true, "output_test_horizontal_diag_fermion");
  auto bparams = make_params(false, "output_test_horizontal_diag_boson");
  tenes::itps::iTPS<tensor> fstate(
      MPI_COMM_WORLD, fparams, lattice_f, updates,
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  tenes::itps::iTPS<tensor> bstate(
      MPI_COMM_WORLD, bparams, lattice_b, updates,
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});

  auto& fTn = tenes::itps::iTPSTestAccessor::Tn(fstate);
  auto& bTn = tenes::itps::iTPSTestAccessor::Tn(bstate);
  auto diff = max_tn_diff(fTn, bTn);
  std::cout << std::setprecision(17)
            << "horizontal trajectory raw-initial maxdiff=" << diff.max_abs
            << " site=" << diff.site << " index=" << index_to_string(diff.index)
            << " fermion=" << diff.fermion_value
            << " boson=" << diff.boson_value << std::endl;

  auto fparams2 =
      make_params(true, "output_test_horizontal_diag_fermion_forced");
  auto bparams2 =
      make_params(false, "output_test_horizontal_diag_boson_forced");
  auto lattice_f2 = make_lattice();
  auto lattice_b2 = make_lattice();
  tenes::itps::iTPS<tensor> fstate2(
      MPI_COMM_WORLD, fparams2, lattice_f2, updates,
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  tenes::itps::iTPS<tensor> bstate2(
      MPI_COMM_WORLD, bparams2, lattice_b2, updates,
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  auto& fTn2 = tenes::itps::iTPSTestAccessor::Tn(fstate2);
  auto& bTn2 = tenes::itps::iTPSTestAccessor::Tn(bstate2);
  auto& flambda2 = tenes::itps::iTPSTestAccessor::lambda_tensor(fstate2);
  auto& blambda2 = tenes::itps::iTPSTestAccessor::lambda_tensor(bstate2);
  bTn2 = fTn2;
  blambda2 = flambda2;
  diff = max_tn_diff(fTn2, bTn2);
  std::cout << std::setprecision(17)
            << "horizontal trajectory forced-initial maxdiff=" << diff.max_abs
            << std::endl;

  bool found = false;
  for (int step = 0; step < 3 && !found; ++step) {
    for (std::size_t update_index = 0; update_index < updates.size();
         ++update_index) {
      fstate2.simple_update(updates[update_index]);
      bstate2.simple_update(updates[update_index]);
      diff = max_tn_diff(fTn2, bTn2);
      std::cout << std::setprecision(17)
                << "horizontal trajectory after step=" << step
                << " update=" << update_index
                << " source_site=" << updates[update_index].source_site
                << " source_leg=" << updates[update_index].source_leg
                << " maxdiff=" << diff.max_abs << " site=" << diff.site
                << " index=" << index_to_string(diff.index)
                << " fermion=" << diff.fermion_value
                << " boson=" << diff.boson_value << std::endl;
      if (diff.max_abs > 1.0e-10) {
        found = true;
        break;
      }
    }
  }
}

TEST_CASE("diagnostic update1 gate application component diff") {
  if (std::getenv("TENES_RUN_UPDATE1_GATE_DIAG") == nullptr) {
    return;
  }

  namespace f = tenes::fermion;
  using tensor = tenes::real_tensor;

  auto make_lattice = []() {
    tenes::SquareLattice lattice(2, 2);
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      lattice.physical_dims[site] = 2;
      lattice.virtual_dims[site] = {2, 1, 2, 1};
      lattice.initial_dirs[site] = {0.0};
      lattice.noises[site] = 1.0;
    }
    return lattice;
  };
  auto make_params = [](bool fermion, const std::string& outdir) {
    tenes::itps::PEPS_Parameters params;
    params.fermion = fermion;
    if (fermion) {
      params.phys_parity.assign(4, f::parity_vector{false, true});
    }
    params.print_level = tenes::PrintLevel::none;
    params.outdir = outdir;
    params.CHI = 8;
    params.Max_CTM_Iteration = 10;
    params.CTM_Convergence_Epsilon = 1.0e-8;
    params.Use_RSVD = false;
    params.seed = 11;
    return params;
  };
  auto make_updates = []() {
    std::vector<tenes::EvolutionOperator<tensor>> updates;
    const tensor gate = make_free_fermion_gate(0.01);
    for (int site = 0; site < 4; ++site) {
      updates.push_back(
          tenes::make_twosite_EvolutionOperator(site, 2, 0, gate));
    }
    return updates;
  };

  auto lattice_f = make_lattice();
  auto lattice_b = make_lattice();
  const auto updates = make_updates();
  auto fparams = make_params(true, "output_test_update1_gate_diag_fermion");
  auto bparams = make_params(false, "output_test_update1_gate_diag_boson");
  tenes::itps::iTPS<tensor> fstate(
      MPI_COMM_WORLD, fparams, lattice_f, updates,
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  tenes::itps::iTPS<tensor> bstate(
      MPI_COMM_WORLD, bparams, lattice_b, updates,
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});

  auto& fTn = tenes::itps::iTPSTestAccessor::Tn(fstate);
  auto& bTn = tenes::itps::iTPSTestAccessor::Tn(bstate);
  auto& flambda = tenes::itps::iTPSTestAccessor::lambda_tensor(fstate);
  auto& blambda = tenes::itps::iTPSTestAccessor::lambda_tensor(bstate);
  bTn = fTn;
  blambda = flambda;

  auto make_fermion_gate = [&](const tenes::EvolutionOperator<tensor>& update,
                               int source, int target) {
    const auto& finfo = tenes::itps::iTPSTestAccessor::finfo(fstate);
    return f::wrap_twosite_op(update.op, finfo.phys[source],
                              finfo.phys[target]);
  };

  auto dump_update = [&](int update_index, const std::string& label) {
    const auto& update = updates[update_index];
    const int source = update.source_site;
    const int source_leg = update.source_leg;
    const int target = lattice_f.neighbor(source, source_leg);
    const auto& finfo = tenes::itps::iTPSTestAccessor::finfo(fstate);
    const auto fTn1 = f::wrap_Tn(fTn[source], finfo, source);
    const auto fTn2 = f::wrap_Tn(fTn[target], finfo, target);
    const auto fop = make_fermion_gate(update, source, target);
    const auto fdata = diagnostic_gate_data(fTn1, fTn2, flambda[source],
                                            flambda[target], fop, source_leg);
    std::cout << label << " source=" << source << " target=" << target
              << " leg=" << source_leg << std::endl;
    dump_parity_metadata(label, fdata.R1, fdata.R2, fdata.theta_before);
    dump_tensordot_mask_pattern(label, fdata.theta_before, fop);
    return fdata;
  };

  const auto update0 = dump_update(0, "update0");

  fstate.simple_update(updates[0]);
  bstate.simple_update(updates[0]);

  const auto update1 = dump_update(1, "update1");
  const int source = updates[1].source_site;
  const int source_leg = updates[1].source_leg;
  const int target = lattice_b.neighbor(source, source_leg);
  const auto bdata =
      diagnostic_gate_data(bTn[source], bTn[target], blambda[source],
                           blambda[target], updates[1].op, source_leg);
  const f::ftensor<tensor> btheta_before_as_f{bdata.theta_before,
                                              update1.theta_before.parity};
  std::cout << "update1 pre_gate_maxdiff="
            << max_tn_diff(std::vector<tensor>{update1.theta_before.t},
                           std::vector<tensor>{bdata.theta_before})
                   .max_abs
            << " pre_gate_parity_violation="
            << f::parity_violation(btheta_before_as_f) << std::endl;
  dump_post_gate_diff("update1_raw", update1.theta_after, bdata.theta_after);

  const auto fop = make_fermion_gate(updates[1], source, target);
  const auto same_input_fermion_after =
      diagnostic_site_ordered_theta(f::tensordot(
          btheta_before_as_f, fop, mptensor::Axes(1, 3), mptensor::Axes(0, 1)));
  dump_post_gate_diff("update1_same_input_gate_only", same_input_fermion_after,
                      bdata.theta_after);

  static_cast<void>(update0);
}

TEST_CASE("diagnostic horizontal-chain sorted lambda trajectory") {
  if (std::getenv("TENES_RUN_LAMBDA_TRAJECTORY_DIAG") == nullptr) {
    return;
  }

  namespace f = tenes::fermion;
  using tensor = tenes::real_tensor;

  const char* d_env = std::getenv("TENES_LAMBDA_DIAG_D");
  const int D = (d_env != nullptr) ? std::atoi(d_env) : 2;
  auto make_lattice = [D]() {
    tenes::SquareLattice lattice(2, 2);
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      lattice.physical_dims[site] = 2;
      lattice.virtual_dims[site] = {D, 1, D, 1};
      lattice.initial_dirs[site] = {0.0};
      lattice.noises[site] = 1.0;
    }
    return lattice;
  };
  auto make_params = [](bool fermion, const std::string& outdir) {
    tenes::itps::PEPS_Parameters params;
    params.fermion = fermion;
    if (fermion) {
      params.phys_parity.assign(4, f::parity_vector{false, true});
    }
    params.print_level = tenes::PrintLevel::none;
    params.outdir = outdir;
    params.CHI = 8;
    params.Max_CTM_Iteration = 10;
    params.CTM_Convergence_Epsilon = 1.0e-8;
    params.Use_RSVD = false;
    params.seed = 11;
    return params;
  };
  auto make_updates = []() {
    std::vector<tenes::EvolutionOperator<tensor>> updates;
    const char* mu_env = std::getenv("TENES_LAMBDA_DIAG_MU");
    const double mu = (mu_env != nullptr) ? std::atof(mu_env) : 0.0;
    tensor gate = make_free_fermion_gate(0.01);
    for (std::size_t n = 0; n < gate.local_size(); ++n) {
      const auto idx = gate.global_index(n);
      double v = 0.0;
      gate.get_value(idx, v);
      const int ntot = static_cast<int>(idx[2]) + static_cast<int>(idx[3]);
      gate.set_value(idx, v * std::exp(0.25 * 0.01 * mu * ntot));
    }
    for (int site = 0; site < 4; ++site) {
      updates.push_back(
          tenes::make_twosite_EvolutionOperator(site, 2, 0, gate));
    }
    return updates;
  };

  auto lattice_f = make_lattice();
  auto lattice_b = make_lattice();
  const auto updates = make_updates();
  auto fparams = make_params(true, "output_test_lambda_diag_fermion");
  auto bparams = make_params(false, "output_test_lambda_diag_boson");
  tenes::itps::iTPS<tensor> fstate(
      MPI_COMM_WORLD, fparams, lattice_f, updates,
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  tenes::itps::iTPS<tensor> bstate(
      MPI_COMM_WORLD, bparams, lattice_b, updates,
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});

  auto& fTn = tenes::itps::iTPSTestAccessor::Tn(fstate);
  auto& bTn = tenes::itps::iTPSTestAccessor::Tn(bstate);
  auto& flambda = tenes::itps::iTPSTestAccessor::lambda_tensor(fstate);
  auto& blambda = tenes::itps::iTPSTestAccessor::lambda_tensor(bstate);
  bTn = fTn;
  blambda = flambda;

  int first_step = -1;
  int first_bond = -1;
  double first_rel = 0.0;
  int first_call_step = -1;
  int first_call_update = -1;
  int first_call_bond = -1;
  double first_call_rel = 0.0;
  double max_rel = 0.0;
  int max_rel_step = 0;
  int max_rel_bond = 0;
  std::vector<double> first_fermion;
  std::vector<double> first_boson;

  for (int step = 0; step <= 300; ++step) {
    for (int site = 0; site < lattice_f.N_UNIT; ++site) {
      const auto fvals = sorted_desc(flambda[site][2]);
      const auto bvals = sorted_desc(blambda[site][2]);
      const double rel = lambda_relative_diff(fvals, bvals);
      if (rel > max_rel) {
        max_rel = rel;
        max_rel_step = step;
        max_rel_bond = site;
      }
      if (first_step < 0 && rel > 1.0e-8) {
        first_step = step;
        first_bond = site;
        first_rel = rel;
        first_fermion = fvals;
        first_boson = bvals;
      }
    }
    if (step == 300) {
      break;
    }
    for (std::size_t update_index = 0; update_index < updates.size();
         ++update_index) {
      const auto& update = updates[update_index];
      fstate.simple_update(update);
      bstate.simple_update(update);
      if (first_call_step < 0) {
        for (int site = 0; site < lattice_f.N_UNIT; ++site) {
          const auto fvals = sorted_desc(flambda[site][2]);
          const auto bvals = sorted_desc(blambda[site][2]);
          const double rel = lambda_relative_diff(fvals, bvals);
          if (rel > 1.0e-8) {
            first_call_step = step;
            first_call_update = static_cast<int>(update_index);
            first_call_bond = site;
            first_call_rel = rel;
            break;
          }
        }
      }
    }
  }

  std::cout << std::setprecision(17) << "lambda trajectory max_rel=" << max_rel
            << " max_step=" << max_rel_step << " max_bond=" << max_rel_bond
            << std::endl;
  if (first_step < 0) {
    std::cout << "lambda trajectory case=A all horizontal bond spectra match "
                 "rtol=1e-8 through 300 steps"
              << std::endl;
    const auto& finfo = tenes::itps::iTPSTestAccessor::finfo(fstate);
    for (int site = 0; site < lattice_f.N_UNIT; ++site) {
      for (int leg = 0; leg < 4; ++leg) {
        std::cout << "lambda trajectory final_finfo site=" << site
                  << " leg=" << leg << " parity=[";
        for (std::size_t i = 0; i < finfo.virt[site][leg].size(); ++i) {
          if (i != 0) {
            std::cout << ",";
          }
          std::cout << (finfo.virt[site][leg][i] ? 1 : 0);
        }
        std::cout << "]" << std::endl;
      }
    }
  } else {
    std::cout << "lambda trajectory case=B first_step=" << first_step
              << " first_bond=" << first_bond << " rel=" << first_rel
              << " fermion=[";
    for (std::size_t i = 0; i < first_fermion.size(); ++i) {
      if (i != 0) {
        std::cout << ",";
      }
      std::cout << first_fermion[i];
    }
    std::cout << "] boson=[";
    for (std::size_t i = 0; i < first_boson.size(); ++i) {
      if (i != 0) {
        std::cout << ",";
      }
      std::cout << first_boson[i];
    }
    std::cout << "]" << std::endl;
    std::cout << "lambda trajectory first_divergent_call step="
              << first_call_step << " update=" << first_call_update
              << " source_site=" << updates[first_call_update].source_site
              << " source_leg=" << updates[first_call_update].source_leg
              << " changed_bond=" << first_call_bond
              << " rel=" << first_call_rel << std::endl;

    auto lattice_fd = make_lattice();
    auto lattice_bd = make_lattice();
    auto fparamsd =
        make_params(true, "output_test_lambda_diag_fermion_first_call");
    auto bparamsd =
        make_params(false, "output_test_lambda_diag_boson_first_call");
    tenes::itps::iTPS<tensor> fdiag(
        MPI_COMM_WORLD, fparamsd, lattice_fd, updates,
        tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
        tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
        tenes::itps::CorrelationParameter{},
        tenes::itps::TransferMatrix_Parameters{});
    tenes::itps::iTPS<tensor> bdiag(
        MPI_COMM_WORLD, bparamsd, lattice_bd, updates,
        tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
        tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
        tenes::itps::CorrelationParameter{},
        tenes::itps::TransferMatrix_Parameters{});
    auto& fdiag_tn = tenes::itps::iTPSTestAccessor::Tn(fdiag);
    auto& bdiag_tn = tenes::itps::iTPSTestAccessor::Tn(bdiag);
    auto& fdiag_lambda = tenes::itps::iTPSTestAccessor::lambda_tensor(fdiag);
    auto& bdiag_lambda = tenes::itps::iTPSTestAccessor::lambda_tensor(bdiag);
    bdiag_tn = fdiag_tn;
    bdiag_lambda = fdiag_lambda;

    for (int step = 0; step < first_call_step; ++step) {
      for (const auto& update : updates) {
        fdiag.simple_update(update);
        bdiag.simple_update(update);
      }
    }
    for (int update_index = 0; update_index < first_call_update;
         ++update_index) {
      fdiag.simple_update(updates[update_index]);
      bdiag.simple_update(updates[update_index]);
    }

    const auto& first_update = updates[first_call_update];
    const int source = first_update.source_site;
    const int source_leg = first_update.source_leg;
    const int target = lattice_fd.neighbor(source, source_leg);
    const auto& finfo = tenes::itps::iTPSTestAccessor::finfo(fdiag);
    const auto fTn1 = f::wrap_Tn(fdiag_tn[source], finfo, source);
    const auto fTn2 = f::wrap_Tn(fdiag_tn[target], finfo, target);
    const auto fop = f::wrap_twosite_op(first_update.op, finfo.phys[source],
                                        finfo.phys[target]);
    const auto ftheta =
        diagnostic_update_thetas(fTn1, fTn2, fdiag_lambda[source],
                                 fdiag_lambda[target], fop, source_leg);
    const auto btheta = diagnostic_update_thetas(
        bdiag_tn[source], bdiag_tn[target], bdiag_lambda[source],
        bdiag_lambda[target], first_update.op, source_leg);
    const f::ftensor<tensor> btheta_before_as_f{btheta.first,
                                                ftheta.first.parity};
    const f::ftensor<tensor> btheta_as_f{btheta.second, ftheta.second.parity};

    const auto f_before_spec = sector_spectrum(
        ftheta.first, mptensor::Axes(0, 1), mptensor::Axes(2, 3));
    const auto b_before_spec = sector_spectrum(
        btheta_before_as_f, mptensor::Axes(0, 1), mptensor::Axes(2, 3));
    const auto f_after_spec = sector_spectrum(
        ftheta.second, mptensor::Axes(0, 1), mptensor::Axes(2, 3));
    const auto b_after_spec = sector_spectrum(btheta_as_f, mptensor::Axes(0, 1),
                                              mptensor::Axes(2, 3));
    std::cout << "lambda trajectory theta first_call before_gate "
              << "fermion_even=" << vector_to_string(f_before_spec.even)
              << " fermion_odd=" << vector_to_string(f_before_spec.odd)
              << " boson_even=" << vector_to_string(b_before_spec.even)
              << " boson_odd=" << vector_to_string(b_before_spec.odd)
              << std::endl;
    std::cout << "lambda trajectory theta first_call after_gate "
              << "fermion_even=" << vector_to_string(f_after_spec.even)
              << " fermion_odd=" << vector_to_string(f_after_spec.odd)
              << " boson_even=" << vector_to_string(b_after_spec.even)
              << " boson_odd=" << vector_to_string(b_after_spec.odd)
              << std::endl;

    tensor fU_plain, fVT_plain;
    std::vector<double> f_plain_s;
    mptensor::svd(ftheta.second.t, mptensor::Axes(0, 1), mptensor::Axes(2, 3),
                  fU_plain, f_plain_s, fVT_plain);
    std::cout << "lambda trajectory theta_f_invariants first_call "
              << "parity_violation=" << f::parity_violation(ftheta.second)
              << " off_block_max="
              << off_block_max(ftheta.second, mptensor::Axes(0, 1),
                               mptensor::Axes(2, 3))
              << " raw_full_svd=" << vector_to_string(f_plain_s) << std::endl;

    f::ftensor<tensor> fU, fVT, bU, bVT;
    tensor bU_plain, bVT_plain;
    std::vector<double> f_s, b_s;
    std::vector<double> b_plain_s;
    f::svd_trunc(ftheta.second, mptensor::Axes(0, 1), mptensor::Axes(2, 3), fU,
                 f_s, fVT, 2);
    f::svd_trunc(btheta_as_f, mptensor::Axes(0, 1), mptensor::Axes(2, 3), bU,
                 b_s, bVT, 2);
    mptensor::svd_trunc(btheta.second, mptensor::Axes(0, 1),
                        mptensor::Axes(2, 3), bU_plain, b_plain_s, bVT_plain,
                        2);
    std::cout << "lambda trajectory fsvd_trunc_selected first_call "
              << "fermion_s=" << vector_to_string(f_s)
              << " boson_manual_sector_s=" << vector_to_string(b_s)
              << " fermion_lambda="
              << vector_to_string(normalized_lambda_from_s(f_s))
              << " boson_manual_sector_lambda="
              << vector_to_string(normalized_lambda_from_s(b_s)) << std::endl;
    std::cout << "lambda trajectory plain_boson_svd_trunc first_call "
              << "boson_plain_s=" << vector_to_string(b_plain_s)
              << " boson_plain_lambda="
              << vector_to_string(normalized_lambda_from_s(b_plain_s))
              << std::endl;
  }
}

TEST_CASE("diagnostic weak-2d coupled chains trajectory") {
  if (std::getenv("TENES_RUN_WEAK2D_DIAG") == nullptr) {
    return;
  }

  namespace f = tenes::fermion;
  using tensor = tenes::real_tensor;

  const char* ty_env = std::getenv("TENES_WEAK2D_TY");
  const double ty = (ty_env != nullptr) ? std::atof(ty_env) : 0.1;
  const char* mu_env = std::getenv("TENES_WEAK2D_MU");
  const double mu = (mu_env != nullptr) ? std::atof(mu_env) : 0.0;
  const char* d_env = std::getenv("TENES_WEAK2D_D");
  const int D = (d_env != nullptr) ? std::atoi(d_env) : 2;
  const char* steps_env = std::getenv("TENES_WEAK2D_STEPS");
  const int nsteps = (steps_env != nullptr) ? std::atoi(steps_env) : 1000;

  tenes::SquareLattice lattice(2, 2);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    lattice.physical_dims[site] = 2;
    lattice.virtual_dims[site] = {D, D, D, D};
    lattice.initial_dirs[site] = {0.0};
    lattice.noises[site] = 1.0;
  }

  tenes::itps::PEPS_Parameters params;
  params.fermion = true;
  params.phys_parity.assign(4, f::parity_vector{false, true});
  params.print_level = tenes::PrintLevel::none;
  params.outdir = "output_test_weak2d";
  params.CHI = 8;
  params.Max_CTM_Iteration = 100;
  params.CTM_Convergence_Epsilon = 1.0e-10;
  params.Use_RSVD = false;
  params.seed = 11;

  // Gate exp(-tau h) with h = -t (c+c + h.c.) - 0.25 mu (n1+n2); the two
  // terms commute, so multiply the hopping gate by the diagonal mu factor.
  auto make_gate = [&](double t, double tau) {
    tensor gate = make_free_fermion_gate(tau * t);
    for (std::size_t n = 0; n < gate.local_size(); ++n) {
      const auto idx = gate.global_index(n);
      double v = 0.0;
      gate.get_value(idx, v);
      const int ntot = static_cast<int>(idx[2]) + static_cast<int>(idx[3]);
      gate.set_value(idx, v * std::exp(0.25 * tau * mu * ntot));
    }
    return gate;
  };
  auto make_updates_tau = [&](double tau) {
    std::vector<tenes::EvolutionOperator<tensor>> ups;
    const tensor hgate = make_gate(1.0, tau);
    const tensor vgate = make_gate(ty, tau);
    for (int site = 0; site < 4; ++site) {
      ups.push_back(tenes::make_twosite_EvolutionOperator(site, 2, 0, hgate));
    }
    for (int site = 0; site < 4; ++site) {
      ups.push_back(tenes::make_twosite_EvolutionOperator(site, 1, 0, vgate));
    }
    return ups;
  };
  const bool anneal = std::getenv("TENES_WEAK2D_ANNEAL") != nullptr;
  std::vector<tenes::EvolutionOperator<tensor>> updates =
      make_updates_tau(0.01);

  tensor hopping(mptensor::Shape(2, 2, 2, 2));
  hopping.set_value(mptensor::Index(0, 1, 1, 0), 1.0);
  hopping.set_value(mptensor::Index(1, 0, 0, 1), 1.0);
  tenes::Operators<tensor> twosite_ops;
  for (int site = 0; site < 4; ++site) {
    twosite_ops.emplace_back("hop", 0, site, std::vector<int>{1},
                             std::vector<int>{0}, hopping);
    twosite_ops.emplace_back("hop", 0, site, std::vector<int>{0},
                             std::vector<int>{1}, hopping);
  }

  tensor n_op(mptensor::Shape(2, 2));
  n_op.set_value(mptensor::Index(1, 1), 1.0);
  tenes::Operators<tensor> onesite_ops;
  for (int site = 0; site < 4; ++site) {
    onesite_ops.emplace_back("n", 0, site, n_op);
  }

  tenes::itps::iTPS<tensor> state(
      MPI_COMM_WORLD, params, lattice, updates,
      tenes::EvolutionOperators<tensor>{}, onesite_ops, twosite_ops,
      tenes::Operators<tensor>{}, tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  auto& lambda = tenes::itps::iTPSTestAccessor::lambda_tensor(state);
  auto& finfo = tenes::itps::iTPSTestAccessor::finfo(state);

  for (int step = 0; step < nsteps; ++step) {
    if (anneal) {
      // Four equal phases with decreasing Trotter step.
      const double taus[4] = {0.1, 0.05, 0.02, 0.01};
      const int phase = std::min(3, 4 * step / std::max(1, nsteps));
      if (step == 0 || phase != std::min(3, 4 * (step - 1) / nsteps)) {
        updates = make_updates_tau(taus[phase]);
        std::cout << "weak2d anneal step=" << step << " tau=" << taus[phase]
                  << std::endl;
      }
    }
    for (const auto& update : updates) {
      state.simple_update(update);
    }
    if ((step + 1) % (nsteps / 10) == 0) {
      std::cout << std::setprecision(10) << "weak2d step=" << step + 1;
      for (int site = 0; site < 4; ++site) {
        const auto h = sorted_desc(lambda[site][2]);
        std::cout << " lamH" << site << "=[" << h[0] << "," << h[1] << "]";
      }
      for (int site = 0; site < 2; ++site) {
        const auto v = sorted_desc(lambda[site][1]);
        std::cout << " lamV" << site << "=[" << v[0] << "," << v[1] << "]";
      }
      std::cout << std::endl;
    }
  }
  for (int site = 0; site < 4; ++site) {
    for (int leg : {2, 1}) {
      std::cout << "weak2d final parity site=" << site << " leg=" << leg
                << " =[";
      for (std::size_t i = 0; i < finfo.virt[site][leg].size(); ++i) {
        std::cout << (finfo.virt[site][leg][i] ? 1 : 0);
      }
      std::cout << "] lambda=[";
      for (const double v : lambda[site][leg]) {
        std::cout << std::setprecision(10) << v << ",";
      }
      std::cout << "]" << std::endl;
    }
  }

  tenes::itps::iTPSTestAccessor::update_reduced_density_environment(state);
  const auto measured = state.measure_twosite();
  REQUIRE(measured.size() >= 1);
  for (int site = 0; site < 4; ++site) {
    std::cout << std::setprecision(10) << "weak2d measured hop site=" << site
              << " dx=" << measured[0].at(tenes::itps::Bond{site, 1, 0})
              << " dy=" << measured[0].at(tenes::itps::Bond{site, 0, 1})
              << std::endl;
  }

  // CTM-independent mean-field (lambda-gauge) estimate per bond via the
  // f-primitive open network: <theta|h|theta>/<theta|theta>.
  auto& Tn = tenes::itps::iTPSTestAccessor::Tn(state);
  const auto fop_h = f::wrap_twosite_op(hopping, f::parity_vector{false, true},
                                        f::parity_vector{false, true});
  auto mf_bond = [&](int source, int source_leg) {
    // Normalize orientation exactly like the driver.
    int s1 = source;
    int s2 = lattice.neighbor(source, source_leg);
    int s1_leg = source_leg;
    if (source_leg == 0 || source_leg == 1) {
      std::swap(s1, s2);
      s1_leg = (source_leg + 2) % 4;
    }
    const auto fTn1 = f::wrap_Tn(Tn[s1], finfo, s1);
    const auto fTn2 = f::wrap_Tn(Tn[s2], finfo, s2);
    const auto data =
        diagnostic_gate_data(fTn1, fTn2, lambda[s1], lambda[s2], fop_h, s1_leg);
    // theta_before: (aux1, p1, aux2, p2); numerator via graded application.
    auto applied = f::tensordot(data.theta_before, fop_h, mptensor::Axes(1, 3),
                                mptensor::Axes(0, 1));
    applied = f::transpose(applied, mptensor::Axes(0, 2, 1, 3));
    const double num =
        f::trace(f::conj(data.theta_before), applied,
                 mptensor::Axes(0, 1, 2, 3), mptensor::Axes(0, 1, 2, 3));
    const double den =
        f::trace(f::conj(data.theta_before), data.theta_before,
                 mptensor::Axes(0, 1, 2, 3), mptensor::Axes(0, 1, 2, 3));
    return num / den;
  };
  for (int site = 0; site < 4; ++site) {
    std::cout << std::setprecision(10) << "weak2d MF hop site=" << site
              << " dx=" << mf_bond(site, 2) << " dy=" << mf_bond(site, 1)
              << std::endl;
  }

  const auto n_measured = state.measure_onesite();
  const auto fop_n = f::ftensor<tensor>{
      n_op, {f::parity_vector{false, true}, f::parity_vector{false, true}}};
  auto mf_n = [&](int site) {
    const auto fTn1 = f::wrap_Tn(Tn[site], finfo, site);
    const int target = lattice.neighbor(site, 2);
    const auto fTn2 = f::wrap_Tn(Tn[target], finfo, target);
    const auto data = diagnostic_gate_data(fTn1, fTn2, lambda[site],
                                           lambda[target], fop_h, 2);
    auto applied = f::tensordot(data.theta_before, fop_n, mptensor::Axes(1),
                                mptensor::Axes(0));
    applied = f::transpose(applied, mptensor::Axes(0, 3, 1, 2));
    const double num =
        f::trace(f::conj(data.theta_before), applied,
                 mptensor::Axes(0, 1, 2, 3), mptensor::Axes(0, 1, 2, 3));
    const double den =
        f::trace(f::conj(data.theta_before), data.theta_before,
                 mptensor::Axes(0, 1, 2, 3), mptensor::Axes(0, 1, 2, 3));
    return num / den;
  };
  for (int site = 0; site < 4; ++site) {
    std::cout << std::setprecision(10) << "weak2d n site=" << site
              << " measured=" << n_measured[0][site] << " MF=" << mf_n(site)
              << std::endl;
  }
}

static tenes::real_tensor electron_gate(double t, double u, double mu,
                                        double tau);

TEST_CASE("diagnostic plaquette kernel vs exact trotter") {
  if (std::getenv("TENES_RUN_PLAQUETTE_TROTTER_DIAG") == nullptr) {
    return;
  }

  namespace f = tenes::fermion;
  using tensor = tenes::real_tensor;
  using fts = f::ftensor<tensor>;

  // Bond dimension 8 with initial support on {index 0 (even), index 1 (odd)}
  // only: theta rank stays <= 8 = dc through the first TWO sweeps, so the
  // kernel performs NO truncation there and must match the exact reference
  // to numerical precision. Sweep 3 onward truncates (reported only).
  constexpr int D = 8;
  f::parity_vector bond_par(D, false);
  bond_par[1] = true;
  const f::parity_vector triv{false};
  // Physical dimension 2 (spinless, parity {e,o}) by default; 4 (electron,
  // parity {e,o,o,e}) with TENES_PLAQUETTE_PHYS_DIM=4. At d = 4 only the
  // first sweep is truncation-free (theta rank <= 2 * 4 = 8 = dc).
  const char* pd_env = std::getenv("TENES_PLAQUETTE_PHYS_DIM");
  const std::size_t PD = (pd_env != nullptr) ? std::atoi(pd_env) : 2;
  const f::parity_vector phys =
      PD == 2 ? f::parity_vector{false, true}
              : f::parity_vector{false, true, true, false};
  std::cout << "plaquette trotter phys_dim=" << PD << std::endl;

  // Open 2x2 patch, sites in raster order: 0=TL, 1=TR, 2=BL, 3=BR.
  // Leg order (l, t, r, b, p); open legs have dimension 1.
  auto make_site = [&](std::size_t l, std::size_t t, std::size_t r,
                       std::size_t b, unsigned seed) {
    f::leg_parities p{l == 1 ? triv : bond_par, t == 1 ? triv : bond_par,
                      r == 1 ? triv : bond_par, b == 1 ? triv : bond_par, phys};
    tensor a(mptensor::Shape(l, t, r, b, PD));
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (std::size_t n = 0; n < a.local_size(); ++n) {
      const auto idx = a.global_index(n);
      const double v = dist(gen);
      bool in_support = true;
      for (int leg = 0; leg < 4; ++leg) {
        if (a.shape()[leg] > 1 && idx[leg] > 1) {
          in_support = false;
        }
      }
      if (in_support && f::count_odd(p, idx) % 2 == 0) {
        a.set_value(idx, v);
      }
    }
    return fts{a, p};
  };

  std::array<fts, 4> Tn = {
      make_site(1, 1, D, D, 101), make_site(D, 1, 1, D, 102),
      make_site(1, D, D, 1, 103), make_site(D, D, 1, 1, 104)};
  std::array<std::vector<std::vector<double>>, 4> lambda;
  for (int s = 0; s < 4; ++s) {
    lambda[s].resize(4);
    for (int leg = 0; leg < 4; ++leg) {
      lambda[s][leg].assign(Tn[s].shape()[leg], 1.0);
    }
  }

  tenes::itps::PEPS_Parameters params;
  params.Inverse_lambda_cut = 1.0e-12;

  const double tau = 0.05;
  const tensor gate =
      PD == 2 ? make_free_fermion_gate(tau) : electron_gate(1.0, 0.0, 0.0, tau);
  const auto fop = [&](int s1, int s2) {
    static_cast<void>(s2);
    static_cast<void>(s1);
    return f::wrap_twosite_op(gate, phys, phys);
  };

  // Contract the patch to the four-site wavefunction psi[p0,p1,p2,p3]
  // (raster order), squeezing the trivial open legs.
  auto contract_patch = [&](const std::array<fts, 4>& T) {
    // TL(l0,t0,r0,b0,p0) x TR over r0-l1
    fts ab = f::tensordot(T[0], T[1], mptensor::Axes(2), mptensor::Axes(0));
    // ab: (l0,t0,b0,p0, t1,r1,b1,p1)
    fts abc = f::tensordot(ab, T[2], mptensor::Axes(2), mptensor::Axes(1));
    // contracted b0 with BL top. abc: (l0,t0,p0,t1,r1,b1,p1, l2,r2,b2,p2)
    fts abcd =
        f::tensordot(abc, T[3], mptensor::Axes(5, 8), mptensor::Axes(1, 0));
    // contracted b1(TR bottom)-t3 and r2(BL right)-l3.
    // abcd: (l0,t0,p0,t1,r1,p1, l2,b2,p2, r3,b3,p3)
    // Move physical legs to raster order, trivial legs first.
    fts ordered = f::transpose(
        abcd, mptensor::Axes(0, 1, 3, 4, 6, 7, 9, 10, 2, 5, 8, 11));
    return f::reshape(ordered, mptensor::Shape(PD, PD, PD, PD));
  };

  // Alternative contraction order (validated R2 plaquette pattern: rows
  // first, then join both vertical bonds) to confirm order independence of
  // the extraction itself.
  auto contract_patch_alt = [&](const std::array<fts, 4>& T) {
    fts top = f::tensordot(T[0], T[1], mptensor::Axes(2), mptensor::Axes(0));
    fts bottom = f::tensordot(T[2], T[3], mptensor::Axes(2), mptensor::Axes(0));
    fts joined =
        f::tensordot(top, bottom, mptensor::Axes(2, 6), mptensor::Axes(1, 4));
    fts ordered = f::transpose(
        joined, mptensor::Axes(0, 1, 3, 4, 6, 7, 9, 10, 2, 5, 8, 11));
    return f::reshape(ordered, mptensor::Shape(PD, PD, PD, PD));
  };

  auto normalized_amplitudes = [&](const fts& psi) {
    std::vector<double> amps(PD * PD * PD * PD, 0.0);
    double norm2 = 0.0;
    for (std::size_t n = 0; n < psi.t.local_size(); ++n) {
      const auto idx = psi.t.global_index(n);
      double v = 0.0;
      psi.t.get_value(idx, v);
      amps[((idx[0] * PD + idx[1]) * PD + idx[2]) * PD + idx[3]] = v;
      norm2 += v * v;
    }
    const double norm = std::sqrt(norm2);
    for (double& a : amps) {
      a /= norm;
    }
    return amps;
  };

  auto compare = [&](const std::vector<double>& a,
                     const std::vector<double>& b) {
    double best = std::numeric_limits<double>::max();
    for (const double sign : {1.0, -1.0}) {
      double diff = 0.0;
      for (std::size_t i = 0; i < a.size(); ++i) {
        diff = std::max(diff, std::abs(a[i] - sign * b[i]));
      }
      best = std::min(best, diff);
    }
    return best;
  };

  // Exact reference: graded two-site gate application on the rank-4 psi.
  auto apply_exact = [&](fts psi, int axis0, int axis1) {
    fts applied =
        f::tensordot(psi, fop(axis0, axis1), mptensor::Axes(axis0, axis1),
                     mptensor::Axes(0, 1));
    // applied: (free..., out0, out1); restore raster order.
    mptensor::Axes perm;
    int free_axis = 0;
    for (int ax = 0; ax < 4; ++ax) {
      if (ax == axis0) {
        perm.push(2);
      } else if (ax == axis1) {
        perm.push(3);
      } else {
        perm.push(free_axis++);
      }
    }
    return f::transpose(applied, perm);
  };

  fts psi_ref = contract_patch(Tn);
  std::cout << std::setprecision(17)
            << "plaquette trotter extraction order check maxdiff="
            << compare(normalized_amplitudes(psi_ref),
                       normalized_amplitudes(contract_patch_alt(Tn)))
            << std::endl;

  // --- third referee: exact Fock-space evolution with Jordan-Wigner ------
  // Modes in raster-site order; d = 4 has (up, dn) per site, so hopping on
  // a vertical bond automatically carries the string over the intervening
  // site's modes. The gates here are pure hopping (t = 1, U = mu = 0).
  const int spins = (PD == 4) ? 2 : 1;
  const int nmodes = 4 * spins;
  const int fock_dim = 1 << nmodes;
  auto fock_mode_op = [&](int g, int mode, bool dagger, double& sign) -> int {
    const int bit = nmodes - 1 - mode;
    const int occ = (g >> bit) & 1;
    if (dagger == (occ == 1)) {
      sign = 0.0;
      return 0;
    }
    int cnt = 0;
    for (int m = 0; m < mode; ++m) {
      cnt += (g >> (nmodes - 1 - m)) & 1;
    }
    sign = (cnt % 2 == 0) ? 1.0 : -1.0;
    return g ^ (1 << bit);
  };
  auto fock_apply_h = [&](const std::vector<double>& v, int sa, int sb) {
    std::vector<double> hv(fock_dim, 0.0);
    for (int g = 0; g < fock_dim; ++g) {
      if (v[g] == 0.0) {
        continue;
      }
      for (int sp = 0; sp < spins; ++sp) {
        const int mode_a = spins * sa + sp;
        const int mode_b = spins * sb + sp;
        for (int dir = 0; dir < 2; ++dir) {
          const int m_annih = dir == 0 ? mode_b : mode_a;
          const int m_creat = dir == 0 ? mode_a : mode_b;
          double s1 = 0.0, s2 = 0.0;
          const int g1 = fock_mode_op(g, m_annih, false, s1);
          if (s1 == 0.0) {
            continue;
          }
          const int g2 = fock_mode_op(g1, m_creat, true, s2);
          if (s2 == 0.0) {
            continue;
          }
          hv[g2] += -1.0 * s1 * s2 * v[g];
        }
      }
    }
    return hv;
  };
  auto fock_apply_gate = [&](std::vector<double>& psi, int sa, int sb) {
    std::vector<double> term = psi;
    for (int order = 1; order <= 24; ++order) {
      term = fock_apply_h(term, sa, sb);
      for (int g = 0; g < fock_dim; ++g) {
        term[g] *= -tau / order;
      }
      for (int g = 0; g < fock_dim; ++g) {
        psi[g] += term[g];
      }
    }
  };
  // local index i -> per-spin occupations; amplitude array index -> Fock g
  auto amps_to_fock = [&](const std::vector<double>& amps) {
    std::vector<double> psi(fock_dim, 0.0);
    const int npd = static_cast<int>(PD);
    for (int i0 = 0; i0 < npd; ++i0) {
      for (int i1 = 0; i1 < npd; ++i1) {
        for (int i2 = 0; i2 < npd; ++i2) {
          for (int i3 = 0; i3 < npd; ++i3) {
            const int idx = ((i0 * npd + i1) * npd + i2) * npd + i3;
            const int loc[4] = {i0, i1, i2, i3};
            int g = 0;
            for (int site = 0; site < 4; ++site) {
              for (int sp = 0; sp < spins; ++sp) {
                // local index i = n_up + 2 n_dn (d=4) or i = n (d=2)
                const int occ = (loc[site] >> sp) & 1;
                const int mode = spins * site + sp;
                g |= occ << (nmodes - 1 - mode);
              }
            }
            psi[g] = amps[idx];
          }
        }
      }
    }
    return psi;
  };
  auto fock_to_amps = [&](const std::vector<double>& psi) {
    std::vector<double> amps(PD * PD * PD * PD, 0.0);
    const int npd = static_cast<int>(PD);
    for (int i0 = 0; i0 < npd; ++i0) {
      for (int i1 = 0; i1 < npd; ++i1) {
        for (int i2 = 0; i2 < npd; ++i2) {
          for (int i3 = 0; i3 < npd; ++i3) {
            const int idx = ((i0 * npd + i1) * npd + i2) * npd + i3;
            const int loc[4] = {i0, i1, i2, i3};
            int g = 0;
            for (int site = 0; site < 4; ++site) {
              for (int sp = 0; sp < spins; ++sp) {
                const int occ = (loc[site] >> sp) & 1;
                const int mode = spins * site + sp;
                g |= occ << (nmodes - 1 - mode);
              }
            }
            amps[idx] = psi[g];
          }
        }
      }
    }
    double norm2 = 0.0;
    for (const double a : amps) {
      norm2 += a * a;
    }
    for (double& a : amps) {
      a /= std::sqrt(norm2);
    }
    return amps;
  };
  std::vector<double> psi_fock = amps_to_fock(normalized_amplitudes(psi_ref));

  struct BondSpec {
    int s1;
    int s1_leg;
    int s2;
    const char* name;
    int axis0;
    int axis1;
  };
  // Kernel bonds in normalized (raster) orientation; exact axes in raster
  // order of (s1, s2).
  std::array<BondSpec, 4> bonds = {
      BondSpec{0, 2, 1, "TL-TR", 0, 1}, BondSpec{2, 2, 3, "BL-BR", 2, 3},
      BondSpec{0, 3, 2, "TL-BL", 0, 2}, BondSpec{1, 3, 3, "TR-BR", 1, 3}};
  if (std::getenv("TENES_PLAQUETTE_VERTICAL_FIRST") != nullptr) {
    std::swap(bonds[0], bonds[2]);
    std::swap(bonds[1], bonds[3]);
    std::cout << "plaquette trotter vertical bonds first" << std::endl;
  }

  for (int step = 0; step < 6; ++step) {
    for (const auto& bond : bonds) {
      fts out1, out2;
      std::vector<double> lambda_work;
      tenes::itps::core::Simple_update_bond(
          Tn[bond.s1], Tn[bond.s2], lambda[bond.s1], lambda[bond.s2],
          fop(bond.s1, bond.s2), bond.s1_leg, params, out1, out2, lambda_work);
      const int s2_leg = (bond.s1_leg + 2) % 4;
      Tn[bond.s1] = out1;
      Tn[bond.s2] = out2;
      lambda[bond.s1][bond.s1_leg] = lambda_work;
      lambda[bond.s2][s2_leg] = lambda_work;

      psi_ref = apply_exact(psi_ref, bond.axis0, bond.axis1);
      fock_apply_gate(psi_fock, bond.axis0, bond.axis1);

      const auto kernel_amps = normalized_amplitudes(contract_patch(Tn));
      const auto ref_amps = normalized_amplitudes(psi_ref);
      const auto fock_amps = fock_to_amps(psi_fock);
      const double diff = compare(kernel_amps, ref_amps);
      std::cout << std::setprecision(17) << "plaquette trotter step=" << step
                << " bond=" << bond.name << " maxdiff=" << diff
                << " kernel_vs_fock=" << compare(kernel_amps, fock_amps)
                << " fprim_vs_fock=" << compare(ref_amps, fock_amps)
                << std::endl;
      if (diff > 1.0e-8 && step == 0) {
        std::cout << "plaquette trotter first divergence amplitudes:"
                  << std::endl;
        for (int i = 0; i < 16; ++i) {
          std::cout << "  amp[" << ((i >> 3) & 1) << ((i >> 2) & 1)
                    << ((i >> 1) & 1) << (i & 1)
                    << "] kernel=" << kernel_amps[i] << " ref=" << ref_amps[i]
                    << std::endl;
        }
      }
    }
  }
}

// --- two-site, four-mode Fock helpers for the electron (d = 4) diagnostics --
//
// Local basis per site: 0:|0>, 1:|up>, 2:|dn>, 3:|up dn> with the internal
// convention |up dn> = c^dag_up c^dag_dn |0>, i.e. local index i = n_up +
// 2 n_dn. Global modes are ordered (1up, 1dn, 2up, 2dn) and a global Fock
// state is encoded in bits g = a<<3 | b<<2 | c<<1 | d for occupations
// (a, b, c, d) of those modes; the creation-string order matches the mode
// order, so Jordan-Wigner signs are (-1)^(number of occupied earlier modes).

// Apply c_mode (dagger=false) or c^dag_mode to basis state g.
// Returns {g', sign}; sign = 0 means the result vanishes.
static std::pair<int, double> electron_mode_op(int g, int mode, bool dagger) {
  const int bit = 3 - mode;  // mode 0 (1up) is the highest bit
  const int occ = (g >> bit) & 1;
  if (dagger == (occ == 1)) {
    return {0, 0.0};
  }
  int string_count = 0;
  for (int m = 0; m < mode; ++m) {
    string_count += (g >> (3 - m)) & 1;
  }
  const double sign = (string_count % 2 == 0) ? 1.0 : -1.0;
  return {g ^ (1 << bit), sign};
}

// 16x16 matrix of the electron bond Hamiltonian
//   h = -t sum_sigma (c^dag_{1 sigma} c_{2 sigma} + h.c.)
//       + (U/4) (n_{1up} n_{1dn} + n_{2up} n_{2dn})   [U split over 4 bonds]
//       - (mu/4) (n_1 + n_2)
// in the ordered Fock basis, indexed by (i1, i2) as i1 * 4 + i2.
static std::array<std::array<double, 16>, 16> electron_bond_hamiltonian(
    double t, double u, double mu) {
  auto local_to_bits = [](int i) {  // i = n_up + 2 n_dn -> (n_up, n_dn)
    return std::pair<int, int>{i & 1, (i >> 1) & 1};
  };
  auto pair_to_g = [&](int i1, int i2) {
    const auto [a, b] = local_to_bits(i1);
    const auto [c, d] = local_to_bits(i2);
    return a << 3 | b << 2 | c << 1 | d;
  };
  std::array<int, 16> g_of_index{};
  std::array<int, 16> index_of_g{};
  for (int i1 = 0; i1 < 4; ++i1) {
    for (int i2 = 0; i2 < 4; ++i2) {
      const int idx = i1 * 4 + i2;
      g_of_index[idx] = pair_to_g(i1, i2);
      index_of_g[g_of_index[idx]] = idx;
    }
  }

  std::array<std::array<double, 16>, 16> h{};
  for (int idx = 0; idx < 16; ++idx) {
    const int g = g_of_index[idx];
    // diagonal: U and mu terms
    const int a = (g >> 3) & 1, b = (g >> 2) & 1, c = (g >> 1) & 1, d = g & 1;
    h[idx][idx] += 0.25 * u * (a * b + c * d) - 0.25 * mu * (a + b + c + d);
    // hopping: -t (c^dag_{1s} c_{2s} + c^dag_{2s} c_{1s});
    // modes: 1up=0, 1dn=1, 2up=2, 2dn=3
    const int hop_pairs[4][2] = {{0, 2}, {2, 0}, {1, 3}, {3, 1}};
    for (const auto& mp : hop_pairs) {
      const auto [g1, s1] = electron_mode_op(g, mp[1], false);
      if (s1 == 0.0) {
        continue;
      }
      const auto [g2, s2] = electron_mode_op(g1, mp[0], true);
      if (s2 == 0.0) {
        continue;
      }
      h[index_of_g[g2]][idx] += -t * s1 * s2;
    }
  }
  return h;
}

// Gate tensor exp(-tau h) as op[in1][in2][out1][out2] via a Taylor series
// (norm(tau h) << 1 for the parameters used here).
static tenes::real_tensor electron_gate(double t, double u, double mu,
                                        double tau) {
  const auto h = electron_bond_hamiltonian(t, u, mu);
  std::array<std::array<double, 16>, 16> gate{};
  std::array<std::array<double, 16>, 16> term{};
  for (int i = 0; i < 16; ++i) {
    gate[i][i] = 1.0;
    term[i][i] = 1.0;
  }
  for (int order = 1; order <= 20; ++order) {
    std::array<std::array<double, 16>, 16> next{};
    for (int i = 0; i < 16; ++i) {
      for (int k = 0; k < 16; ++k) {
        if (term[i][k] == 0.0) {
          continue;
        }
        const double w = term[i][k] * (-tau) / order;
        for (int j = 0; j < 16; ++j) {
          next[i][j] += w * h[k][j];
        }
      }
    }
    term = next;
    for (int i = 0; i < 16; ++i) {
      for (int j = 0; j < 16; ++j) {
        gate[i][j] += term[i][j];
      }
    }
  }
  tenes::real_tensor op(mptensor::Shape(4, 4, 4, 4));
  for (int i1 = 0; i1 < 4; ++i1) {
    for (int i2 = 0; i2 < 4; ++i2) {
      for (int o1 = 0; o1 < 4; ++o1) {
        for (int o2 = 0; o2 < 4; ++o2) {
          const double v = gate[o1 * 4 + o2][i1 * 4 + i2];
          if (std::abs(v) > 1.0e-16) {
            op.set_value(mptensor::Index(i1, i2, o1, o2), v);
          }
        }
      }
    }
  }
  return op;
}

TEST_CASE("diagnostic electron chain sorted lambda trajectory") {
  if (std::getenv("TENES_RUN_ELECTRON_CHAIN_DIAG") == nullptr) {
    return;
  }

  namespace f = tenes::fermion;
  using tensor = tenes::real_tensor;

  const char* u_env = std::getenv("TENES_ELECTRON_U");
  const double u = (u_env != nullptr) ? std::atof(u_env) : 0.0;
  const char* mu_env = std::getenv("TENES_ELECTRON_MU");
  const double mu = (mu_env != nullptr) ? std::atof(mu_env) : 0.0;
  const char* leg_env = std::getenv("TENES_ELECTRON_LEG");
  const int source_leg = (leg_env != nullptr) ? std::atoi(leg_env) : 2;
  const bool horizontal = (source_leg == 0 || source_leg == 2);
  // Optional parity-sorted local basis |0>, |up dn>, |up>, |dn| to test
  // whether an unsorted physical parity pattern is handled correctly.
  const bool sorted_basis = std::getenv("TENES_ELECTRON_SORTED") != nullptr;

  std::cout << "electron chain diag U=" << u << " mu=" << mu
            << " source_leg=" << source_leg
            << " sorted_basis=" << (sorted_basis ? 1 : 0) << std::endl;

  const char* lx_env = std::getenv("TENES_ELECTRON_LX");
  const int LX = (lx_env != nullptr) ? std::atoi(lx_env) : 2;
  auto make_lattice = [&]() {
    tenes::SquareLattice lattice(LX, 2);
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      lattice.physical_dims[site] = 4;
      lattice.virtual_dims[site] = horizontal ? std::array<int, 4>{2, 1, 2, 1}
                                              : std::array<int, 4>{1, 2, 1, 2};
      lattice.initial_dirs[site] = {0.0};
      lattice.noises[site] = 1.0;
    }
    return lattice;
  };
  const f::parity_vector electron_parity =
      sorted_basis ? f::parity_vector{false, false, true, true}
                   : f::parity_vector{false, true, true, false};
  auto make_params = [&](bool fermion, const std::string& outdir) {
    tenes::itps::PEPS_Parameters params;
    params.fermion = fermion;
    if (fermion) {
      params.phys_parity.assign(2 * LX, electron_parity);
    }
    params.print_level = tenes::PrintLevel::none;
    params.outdir = outdir;
    params.CHI = 8;
    params.Max_CTM_Iteration = 10;
    params.CTM_Convergence_Epsilon = 1.0e-8;
    params.Use_RSVD = false;
    params.seed = 11;
    return params;
  };
  tensor gate = electron_gate(1.0, u, mu, 0.01);
  if (sorted_basis) {
    // basis permutation: new index -> old index {0, 3, 1, 2}
    const int perm[4] = {0, 3, 1, 2};
    tensor permuted(mptensor::Shape(4, 4, 4, 4));
    for (int i1 = 0; i1 < 4; ++i1) {
      for (int i2 = 0; i2 < 4; ++i2) {
        for (int o1 = 0; o1 < 4; ++o1) {
          for (int o2 = 0; o2 < 4; ++o2) {
            double v = 0.0;
            gate.get_value(
                mptensor::Index(perm[i1], perm[i2], perm[o1], perm[o2]), v);
            if (std::abs(v) > 1.0e-16) {
              permuted.set_value(mptensor::Index(i1, i2, o1, o2), v);
            }
          }
        }
      }
    }
    gate = permuted;
  }
  std::vector<tenes::EvolutionOperator<tensor>> updates;
  for (int site = 0; site < 2 * LX; ++site) {
    updates.push_back(
        tenes::make_twosite_EvolutionOperator(site, source_leg, 0, gate));
  }

  auto lattice_f = make_lattice();
  auto lattice_b = make_lattice();
  auto fparams = make_params(true, "output_test_electron_chain_fermion");
  auto bparams = make_params(false, "output_test_electron_chain_boson");
  tenes::itps::iTPS<tensor> fstate(
      MPI_COMM_WORLD, fparams, lattice_f, updates,
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  tenes::itps::iTPS<tensor> bstate(
      MPI_COMM_WORLD, bparams, lattice_b, updates,
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});

  auto& fTn = tenes::itps::iTPSTestAccessor::Tn(fstate);
  auto& bTn = tenes::itps::iTPSTestAccessor::Tn(bstate);
  auto& flambda = tenes::itps::iTPSTestAccessor::lambda_tensor(fstate);
  auto& blambda = tenes::itps::iTPSTestAccessor::lambda_tensor(bstate);
  bTn = fTn;
  blambda = flambda;

  const int bond_leg = source_leg;
  double max_rel = 0.0;
  int first_step = -1;
  int first_bond = -1;
  double first_rel = 0.0;
  int first_call_update = -1;
  for (int step = 0; step <= 300; ++step) {
    for (int site = 0; site < lattice_f.N_UNIT; ++site) {
      const auto fvals = sorted_desc(flambda[site][bond_leg]);
      const auto bvals = sorted_desc(blambda[site][bond_leg]);
      const double rel = lambda_relative_diff(fvals, bvals);
      max_rel = std::max(max_rel, rel);
      if (first_step < 0 && rel > 1.0e-8) {
        first_step = step;
        first_bond = site;
        first_rel = rel;
      }
    }
    if (step == 300) {
      break;
    }
    for (std::size_t iu = 0; iu < updates.size(); ++iu) {
      const auto& update = updates[iu];
      if (first_call_update < 0) {
        // Compare theta spectra for this call BEFORE applying it.
        const int src = update.source_site;
        const int tgt = lattice_f.neighbor(src, update.source_leg);
        int s1 = src, s2 = tgt, s1_leg = update.source_leg;
        tensor op12 = update.op;
        if (s1_leg == 0 || s1_leg == 1) {
          std::swap(s1, s2);
          s1_leg = (s1_leg + 2) % 4;
          op12 = mptensor::transpose(update.op, mptensor::Axes(1, 0, 3, 2));
        }
        const auto& finfo = tenes::itps::iTPSTestAccessor::finfo(fstate);
        const auto fTn1 = f::wrap_Tn(fTn[s1], finfo, s1);
        const auto fTn2 = f::wrap_Tn(fTn[s2], finfo, s2);
        const auto fop =
            f::wrap_twosite_op(op12, finfo.phys[s1], finfo.phys[s2]);
        const auto fthetas = diagnostic_update_thetas(fTn1, fTn2, flambda[s1],
                                                      flambda[s2], fop, s1_leg);
        const auto bthetas = diagnostic_update_thetas(
            bTn[s1], bTn[s2], blambda[s1], blambda[s2], op12, s1_leg);
        auto raw_svals = [](const tensor& th) {
          tensor uu, vv;
          std::vector<double> ss;
          mptensor::svd(th, mptensor::Axes(0, 1), mptensor::Axes(2, 3), uu, ss,
                        vv);
          return ss;
        };
        const auto f_before = raw_svals(fthetas.first.t);
        const auto b_before = raw_svals(bthetas.first);
        const auto f_after = raw_svals(fthetas.second.t);
        const auto b_after = raw_svals(bthetas.second);
        const double rel_before =
            lambda_relative_diff(sorted_desc(f_before), sorted_desc(b_before));
        const double rel_after =
            lambda_relative_diff(sorted_desc(f_after), sorted_desc(b_after));
        if (rel_before > 1.0e-9 || rel_after > 1.0e-9) {
          first_call_update = static_cast<int>(iu);
          std::cout << std::setprecision(17)
                    << "electron first divergent call step=" << step
                    << " update=" << iu << " source=" << src
                    << " leg=" << static_cast<int>(update.source_leg)
                    << " theta_rel_before=" << rel_before
                    << " theta_rel_after=" << rel_after << std::endl;
          std::cout << "  f_after=" << vector_to_string(f_after) << std::endl;
          std::cout << "  b_after=" << vector_to_string(b_after) << std::endl;
          if (std::getenv("TENES_ELECTRON_DUMP_THETA") != nullptr) {
            auto dump = [&](const char* label, const tensor& th) {
              for (std::size_t n = 0; n < th.local_size(); ++n) {
                const auto tidx = th.global_index(n);
                double v = 0.0;
                th.get_value(tidx, v);
                std::cout << "THETA " << label << " " << tidx[0] << " "
                          << tidx[1] << " " << tidx[2] << " " << tidx[3] << " "
                          << std::setprecision(17) << v << std::endl;
              }
            };
            dump("f_before", fthetas.first.t);
            dump("f_after", fthetas.second.t);
            dump("b_after", bthetas.second);
          }
          // Same-input check: graded gate application on the BOSON theta
          // components (with fermion parity metadata). Separates the QR
          // sign-gauge question from the gate-application question.
          const f::ftensor<tensor> btheta_as_f{bthetas.first,
                                               fthetas.first.parity};
          const auto same_after = diagnostic_site_ordered_theta(f::tensordot(
              btheta_as_f, fop, mptensor::Axes(1, 3), mptensor::Axes(0, 1)));
          const auto s_same = raw_svals(same_after.t);
          std::cout << "  same_input_rel_after="
                    << lambda_relative_diff(sorted_desc(s_same),
                                            sorted_desc(b_after))
                    << " s_same=" << vector_to_string(s_same) << std::endl;
        }
      }
      fstate.simple_update(update);
      bstate.simple_update(update);
    }
  }

  std::cout << std::setprecision(17)
            << "electron chain lambda max_rel=" << max_rel << std::endl;
  if (first_step < 0) {
    std::cout << "electron chain lambda case=A all bond spectra match "
                 "rtol=1e-8 through 300 steps"
              << std::endl;
  } else {
    std::cout << "electron chain lambda case=B first_step=" << first_step
              << " first_bond=" << first_bond << " rel=" << first_rel
              << std::endl;
  }
}

TEST_CASE("diagnostic vertical-chain sorted lambda trajectory") {
  if (std::getenv("TENES_RUN_VERTICAL_LAMBDA_DIAG") == nullptr) {
    return;
  }

  namespace f = tenes::fermion;
  using tensor = tenes::real_tensor;

  auto make_lattice = []() {
    tenes::SquareLattice lattice(2, 2);
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      lattice.physical_dims[site] = 2;
      lattice.virtual_dims[site] = {1, 2, 1, 2};
      lattice.initial_dirs[site] = {0.0};
      lattice.noises[site] = 1.0;
    }
    return lattice;
  };
  auto make_params = [](bool fermion, const std::string& outdir) {
    tenes::itps::PEPS_Parameters params;
    params.fermion = fermion;
    if (fermion) {
      params.phys_parity.assign(4, f::parity_vector{false, true});
    }
    params.print_level = tenes::PrintLevel::none;
    params.outdir = outdir;
    params.CHI = 8;
    params.Max_CTM_Iteration = 10;
    params.CTM_Convergence_Epsilon = 1.0e-8;
    params.Use_RSVD = false;
    params.seed = 11;
    return params;
  };
  const char* leg_env = std::getenv("TENES_VERTICAL_DIAG_LEG");
  const int source_leg = (leg_env != nullptr) ? std::atoi(leg_env) : 1;
  const int bond_leg = source_leg;
  auto make_updates = [source_leg]() {
    std::vector<tenes::EvolutionOperator<tensor>> updates;
    const char* mu_env = std::getenv("TENES_LAMBDA_DIAG_MU");
    const double mu = (mu_env != nullptr) ? std::atof(mu_env) : 0.0;
    tensor gate = make_free_fermion_gate(0.01);
    for (std::size_t n = 0; n < gate.local_size(); ++n) {
      const auto idx = gate.global_index(n);
      double v = 0.0;
      gate.get_value(idx, v);
      const int ntot = static_cast<int>(idx[2]) + static_cast<int>(idx[3]);
      gate.set_value(idx, v * std::exp(0.25 * 0.01 * mu * ntot));
    }
    for (int site = 0; site < 4; ++site) {
      updates.push_back(
          tenes::make_twosite_EvolutionOperator(site, source_leg, 0, gate));
    }
    return updates;
  };

  auto lattice_f = make_lattice();
  auto lattice_b = make_lattice();
  const auto updates = make_updates();
  std::cout << "vertical lambda trajectory source_leg=" << source_leg
            << std::endl;
  auto fparams = make_params(true, "output_test_vertical_lambda_fermion");
  auto bparams = make_params(false, "output_test_vertical_lambda_boson");
  tenes::itps::iTPS<tensor> fstate(
      MPI_COMM_WORLD, fparams, lattice_f, updates,
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  tenes::itps::iTPS<tensor> bstate(
      MPI_COMM_WORLD, bparams, lattice_b, updates,
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});

  auto& fTn = tenes::itps::iTPSTestAccessor::Tn(fstate);
  auto& bTn = tenes::itps::iTPSTestAccessor::Tn(bstate);
  auto& flambda = tenes::itps::iTPSTestAccessor::lambda_tensor(fstate);
  auto& blambda = tenes::itps::iTPSTestAccessor::lambda_tensor(bstate);
  bTn = fTn;
  blambda = flambda;

  double max_rel = 0.0;
  int max_rel_step = -1;
  int max_rel_bond = -1;
  int first_step = -1;
  int first_bond = -1;
  double first_rel = 0.0;
  for (int step = 0; step <= 300; ++step) {
    for (int site = 0; site < lattice_f.N_UNIT; ++site) {
      const auto fvals = sorted_desc(flambda[site][bond_leg]);
      const auto bvals = sorted_desc(blambda[site][bond_leg]);
      const double rel = lambda_relative_diff(fvals, bvals);
      if (rel > max_rel) {
        max_rel = rel;
        max_rel_step = step;
        max_rel_bond = site;
      }
      if (first_step < 0 && rel > 1.0e-8) {
        first_step = step;
        first_bond = site;
        first_rel = rel;
      }
    }
    if (step == 300) {
      break;
    }
    for (const auto& update : updates) {
      fstate.simple_update(update);
      bstate.simple_update(update);
    }
  }

  std::cout << std::setprecision(17)
            << "vertical lambda trajectory max_rel=" << max_rel
            << " max_step=" << max_rel_step << " max_bond=" << max_rel_bond
            << std::endl;
  if (first_step < 0) {
    std::cout << "vertical lambda trajectory case=A all vertical bond spectra "
                 "match rtol=1e-8 through 300 steps"
              << std::endl;
  } else {
    std::cout << "vertical lambda trajectory case=B first_step=" << first_step
              << " first_bond=" << first_bond << " rel=" << first_rel
              << std::endl;
  }
}

TEST_CASE("diagnostic horizontal-chain measured energy fermion vs boson") {
  if (std::getenv("TENES_RUN_CHAIN_MEASURE_DIAG") == nullptr) {
    return;
  }

  namespace f = tenes::fermion;
  using tensor = tenes::real_tensor;

  auto make_lattice = []() {
    tenes::SquareLattice lattice(2, 2);
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      lattice.physical_dims[site] = 2;
      lattice.virtual_dims[site] = {2, 1, 2, 1};
      lattice.initial_dirs[site] = {0.0};
      lattice.noises[site] = 1.0;
    }
    return lattice;
  };
  auto make_params = [](bool fermion, const std::string& outdir) {
    tenes::itps::PEPS_Parameters params;
    params.fermion = fermion;
    if (fermion) {
      params.phys_parity.assign(4, f::parity_vector{false, true});
    }
    params.print_level = tenes::PrintLevel::none;
    params.outdir = outdir;
    params.CHI = 8;
    params.Max_CTM_Iteration = 100;
    params.CTM_Convergence_Epsilon = 1.0e-10;
    params.Use_RSVD = false;
    params.seed = 11;
    return params;
  };
  auto make_updates = []() {
    std::vector<tenes::EvolutionOperator<tensor>> updates;
    const tensor gate = make_free_fermion_gate(0.01);
    for (int site = 0; site < 4; ++site) {
      updates.push_back(
          tenes::make_twosite_EvolutionOperator(site, 2, 0, gate));
    }
    return updates;
  };

  tensor hopping(mptensor::Shape(2, 2, 2, 2));
  hopping.set_value(mptensor::Index(0, 1, 1, 0), 1.0);
  hopping.set_value(mptensor::Index(1, 0, 0, 1), 1.0);
  tenes::Operators<tensor> twosite_ops;
  for (int site = 0; site < 4; ++site) {
    twosite_ops.emplace_back("hop", 0, site, std::vector<int>{1},
                             std::vector<int>{0}, hopping);
  }

  auto lattice_f = make_lattice();
  auto lattice_b = make_lattice();
  const auto updates = make_updates();
  auto fparams = make_params(true, "output_test_chain_measure_fermion");
  auto bparams = make_params(false, "output_test_chain_measure_boson");
  tenes::itps::iTPS<tensor> fstate(MPI_COMM_WORLD, fparams, lattice_f, updates,
                                   tenes::EvolutionOperators<tensor>{},
                                   tenes::Operators<tensor>{}, twosite_ops,
                                   tenes::Operators<tensor>{},
                                   tenes::itps::CorrelationParameter{},
                                   tenes::itps::TransferMatrix_Parameters{});
  tenes::itps::iTPS<tensor> bstate(MPI_COMM_WORLD, bparams, lattice_b, updates,
                                   tenes::EvolutionOperators<tensor>{},
                                   tenes::Operators<tensor>{}, twosite_ops,
                                   tenes::Operators<tensor>{},
                                   tenes::itps::CorrelationParameter{},
                                   tenes::itps::TransferMatrix_Parameters{});

  auto& fTn = tenes::itps::iTPSTestAccessor::Tn(fstate);
  auto& bTn = tenes::itps::iTPSTestAccessor::Tn(bstate);
  auto& flambda = tenes::itps::iTPSTestAccessor::lambda_tensor(fstate);
  auto& blambda = tenes::itps::iTPSTestAccessor::lambda_tensor(bstate);
  bTn = fTn;
  blambda = flambda;

  for (int step = 0; step < 300; ++step) {
    for (const auto& update : updates) {
      fstate.simple_update(update);
      bstate.simple_update(update);
    }
  }

  const auto diff = max_tn_diff(fTn, bTn);
  std::cout << std::setprecision(17)
            << "chain measure diag post-evolution component maxdiff="
            << diff.max_abs << std::endl;
  for (int site = 0; site < 4; ++site) {
    std::cout << "chain measure diag lambda site=" << site << " fermion=["
              << flambda[site][2][0] << "," << flambda[site][2][1]
              << "] boson=[" << blambda[site][2][0] << ","
              << blambda[site][2][1] << "]" << std::endl;
  }

  tenes::itps::iTPSTestAccessor::update_reduced_density_environment(fstate);
  const auto fmeasured = fstate.measure_twosite();
  bstate.update_CTM();
  const auto bmeasured = bstate.measure_twosite();
  REQUIRE(fmeasured.size() >= 1);
  REQUIRE(bmeasured.size() >= 1);
  for (int site = 0; site < 4; ++site) {
    const auto key = tenes::itps::Bond{site, 1, 0};
    std::cout << "chain measure diag hop site=" << site
              << " fermion=" << fmeasured[0].at(key)
              << " boson=" << bmeasured[0].at(key) << std::endl;
  }
}

TEST_CASE(
    "reduced density measurement remains positive after fermion simple "
    "update") {
  namespace f = tenes::fermion;
  using tensor = tenes::real_tensor;

  tenes::SquareLattice lattice(2, 2);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    lattice.physical_dims[site] = 2;
    lattice.virtual_dims[site] = {2, 2, 2, 2};
    lattice.initial_dirs[site] = {0.0};
    lattice.noises[site] = 1.0;
  }

  tenes::itps::PEPS_Parameters params;
  params.fermion = true;
  params.phys_parity.assign(lattice.N_UNIT, f::parity_vector{false, true});
  params.print_level = tenes::PrintLevel::none;
  params.outdir = "output_test_fermion_post_update_measurement";
  params.CHI = 8;
  params.Max_CTM_Iteration = 10;
  params.CTM_Convergence_Epsilon = 1.0e-8;
  params.Use_RSVD = false;
  params.seed = 11;

  std::vector<tenes::EvolutionOperator<tensor>> updates;
  const tensor gate = make_free_fermion_gate(0.01);
  for (int leg : {2, 1}) {
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      updates.push_back(
          tenes::make_twosite_EvolutionOperator(site, leg, 0, gate));
    }
  }

  tenes::itps::iTPS<tensor> state(
      MPI_COMM_WORLD, params, lattice, updates,
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  for (int step = 0; step < 50; ++step) {
    for (const auto& update : updates) {
      state.simple_update(update);
    }
  }

  auto& finfo = tenes::itps::iTPSTestAccessor::finfo(state);
  f::validate_neighbor_consistency(finfo, lattice);
  auto& Tn = tenes::itps::iTPSTestAccessor::Tn(state);
  auto& lambda = tenes::itps::iTPSTestAccessor::lambda_tensor(state);
  static_cast<void>(lambda);
  std::vector<tensor> reduced = f::build_reduced_density_tensors(Tn, finfo);
  std::vector<tensor> C1(lattice.N_UNIT), C2(lattice.N_UNIT),
      C3(lattice.N_UNIT), C4(lattice.N_UNIT), eTt(lattice.N_UNIT),
      eTr(lattice.N_UNIT), eTb(lattice.N_UNIT), eTl(lattice.N_UNIT);
  tenes::itps::core::Calc_CTM_Environment_density(
      C1, C2, C3, C4, eTt, eTr, eTb, eTl, reduced, params, lattice);

  tensor id(mptensor::Shape(2, 2));
  tensor number(mptensor::Shape(2, 2));
  id.set_value(mptensor::Index(0, 0), 1.0);
  id.set_value(mptensor::Index(1, 1), 1.0);
  number.set_value(mptensor::Index(1, 1), 1.0);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    const double norm = tenes::itps::core::Contract_one_site_density_CTM(
        C1[site], C2[site], C3[site], C4[site], eTt[site], eTr[site], eTb[site],
        eTl[site], reduced[site], id);
    const double n =
        tenes::itps::core::Contract_one_site_density_CTM(
            C1[site], C2[site], C3[site], C4[site], eTt[site], eTr[site],
            eTb[site], eTl[site], reduced[site], number) /
        norm;
    CHECK(norm > 0.0);
    CHECK(n >= doctest::Approx(0.0).epsilon(1.0e-12));
    CHECK(n <= doctest::Approx(1.0).epsilon(1.0e-12));
  }
}

TEST_CASE("fermion twosite measurement is translation invariant across wraps") {
  namespace f = tenes::fermion;
  using tensor = tenes::real_tensor;

  tenes::SquareLattice lattice(2, 2);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    lattice.physical_dims[site] = 2;
    lattice.virtual_dims[site] = {2, 2, 2, 2};
    lattice.initial_dirs[site] = {0.0};
    lattice.noises[site] = 1.0;
  }

  tenes::itps::PEPS_Parameters params;
  params.fermion = true;
  params.phys_parity.assign(lattice.N_UNIT, f::parity_vector{false, true});
  params.print_level = tenes::PrintLevel::none;
  params.outdir = "output_test_fermion_wrapped_bond_translation";
  params.CHI = 8;
  params.Max_CTM_Iteration = 100;
  params.CTM_Convergence_Epsilon = 1.0e-12;
  params.Use_RSVD = false;
  params.seed = 17;

  std::vector<tenes::EvolutionOperator<tensor>> updates;
  const tensor gate = make_free_fermion_gate(0.01);
  for (int leg : {2, 1}) {
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      updates.push_back(
          tenes::make_twosite_EvolutionOperator(site, leg, 0, gate));
    }
  }

  tensor hopping(mptensor::Shape(2, 2, 2, 2));
  hopping.set_value(mptensor::Index(0, 1, 0, 1), -0.175);
  hopping.set_value(mptensor::Index(1, 0, 1, 0), -0.175);
  hopping.set_value(mptensor::Index(1, 1, 1, 1), -0.35);
  hopping.set_value(mptensor::Index(0, 1, 1, 0), -1.0);
  hopping.set_value(mptensor::Index(1, 0, 0, 1), -1.0);
  tenes::Operators<tensor> twosite_ops;
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    twosite_ops.emplace_back("hop", 0, site, std::vector<int>{1},
                             std::vector<int>{0}, hopping);
    twosite_ops.emplace_back("hop", 0, site, std::vector<int>{0},
                             std::vector<int>{1}, hopping);
  }

  tenes::itps::iTPS<tensor> state(MPI_COMM_WORLD, params, lattice, updates,
                                  tenes::EvolutionOperators<tensor>{},
                                  tenes::Operators<tensor>{}, twosite_ops,
                                  tenes::Operators<tensor>{},
                                  tenes::itps::CorrelationParameter{},
                                  tenes::itps::TransferMatrix_Parameters{});
  auto& finfo = tenes::itps::iTPSTestAccessor::finfo(state);
  auto& Tn = tenes::itps::iTPSTestAccessor::Tn(state);
  tensor uniform(mptensor::Shape(2, 2, 2, 2, 2));
  const auto parity = f::Tn_parity(finfo, 0);
  for (std::size_t n = 0; n < uniform.local_size(); ++n) {
    const auto idx = uniform.global_index(n);
    if (f::count_odd(parity, idx) % 2 == 0) {
      uniform.set_value(idx, 0.125 * static_cast<double>(1 + n));
    }
  }
  for (auto& site_tensor : Tn) {
    site_tensor = uniform;
  }
  auto& lambda = tenes::itps::iTPSTestAccessor::lambda_tensor(state);
  for (auto& site_lambda : lambda) {
    for (auto& leg_lambda : site_lambda) {
      leg_lambda = {1.0, 0.75};
    }
  }
  tenes::itps::iTPSTestAccessor::update_reduced_density_environment(state);
  const auto measured = state.measure_twosite();
  REQUIRE(measured.size() >= 1);

  const auto& hop = measured[0];
  const std::array<double, 4> hx = {
      hop.at(tenes::itps::Bond{0, 1, 0}), hop.at(tenes::itps::Bond{1, 1, 0}),
      hop.at(tenes::itps::Bond{2, 1, 0}), hop.at(tenes::itps::Bond{3, 1, 0})};
  const std::array<double, 4> hy = {
      hop.at(tenes::itps::Bond{0, 0, 1}), hop.at(tenes::itps::Bond{1, 0, 1}),
      hop.at(tenes::itps::Bond{2, 0, 1}), hop.at(tenes::itps::Bond{3, 0, 1})};

  std::cout << std::setprecision(17) << "fermion wrapped-bond translation hx=["
            << hx[0] << ", " << hx[1] << ", " << hx[2] << ", " << hx[3]
            << "] hy=[" << hy[0] << ", " << hy[1] << ", " << hy[2] << ", "
            << hy[3] << "]" << std::endl;

  const double horizontal_scale =
      std::max({1.0e-12, std::abs(hx[0]), std::abs(hx[2])});
  const double vertical_scale =
      std::max({1.0e-12, std::abs(hy[0]), std::abs(hy[1])});
  CHECK(std::abs(hx[1]) < 10.0 * horizontal_scale);
  CHECK(std::abs(hx[3]) < 10.0 * horizontal_scale);
  CHECK(std::abs(hy[2]) < 10.0 * vertical_scale);
  CHECK(std::abs(hy[3]) < 10.0 * vertical_scale);
}

#include "fermion/r2_convention.cpp"
