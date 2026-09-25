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

//! @file
//! Task T2 of docs/superpowers/plans/2026-09-25-fermion-longrange-measure.md:
//! long-range two-site observables, the ops form and parity-odd one-site
//! operators in fermion mode with the CTM environment (design
//! docs/superpowers/specs/2026-09-25-fermion-longrange-measure-design.md,
//! sections 3 to 5). Behaviour contract: work/fermion-longrange/t2/
//! contract.md, items 1 to 11 (item 12 lives in the rewritten guard tests of
//! test/input.cpp and test/fermion/fold_geometry.cpp).
//!
//! Conventions fixed for every case in this file:
//!   - One-site operators use TeNeS' layout op[in, out] = <out|A|in>; two-site
//!     operators op[in_s, in_t, out_s, out_t], source first.
//!   - d = 2 basis {|0>, |1>}, ledger [e, o]. d = 4 basis {|0>, |up>, |dn>,
//!     |updn> = c+_up c+_dn |0>}, ledger [e, o, o, e].
//!   - The product A_s B_t is written out by the test itself (lr_product_el),
//!     from the contract's formula op4[i_s, i_t, o_s, o_t] = (-1)^{p_B p(i_s)}
//!     A[i_s, o_s] B[i_t, o_t]; that table is pinned to Jordan-Wigner
//!     matrices in a [truth] case, never through product_twosite_op().
//!   - A window of the solver follows twosite_obs.cpp's bosonic CTM branch:
//!     ncol = |dx| + 1, nrow = |dy| + 1, the source at column 0 (dx >= 0) or
//!     ncol - 1, at row nrow - 1 (dy >= 0) or 0, the target in the opposite
//!     corner; window cell (row, col) holds lattice.other(source, col -
//!     source_col, source_row - row); the corners C1..C4 and the edges eTt,
//!     eTr, eTb, eTl are those of the window's corner / edge sites.
//!
//! Truth sources (never the code under test):
//!   - contract item 3: the exact contraction of an open patch (every
//!     perimeter leg dimension 1 and even), built with T1's relay builders,
//!     against literals from test/fermion/fock_oracle.py;
//!   - contract items 5, 8, 10, 11: the "direct" window value, assembled here
//!     from the solver's own environment tensors with T1's
//!     build_relay_window() and core::Contract_density_CTM();
//!   - contract item 6 anchors that direct value on the existing
//!     nearest-neighbour path of measure_twosite() (bundled-k), which was
//!     verified end to end before this task;
//!   - contract item 4: the one-site reduced tensor closed with
//!     core::Contract_one_site_density_CTM(), as onesite_obs.cpp does today.
//!
//! Cases whose name carries [truth] check the reference side only, and cases
//! whose name carries [kept] pin a rejection that the stub already performs
//! (a guard that must survive the change); both are expected to pass against
//! the stub. Every other case must fail until T2 is implemented.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../test_fermion_common.hpp"

#include <algorithm>
#include <complex>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "../../src/fermion/relay.hpp"
#include "../../src/iTPS/load_toml.hpp"

namespace {

namespace rf = tenes::fermion;
using lr_acc = tenes::itps::iTPSTestAccessor;
using tenes::complex_tensor;
using tenes::real_tensor;
using tenes::itps::Bond;

template <class tensor>
using lr_ft = rf::ftensor<tensor>;

template <class tensor>
using lr_state = tenes::itps::iTPS<tensor>;

// ---- small helpers ---------------------------------------------------------

template <class tensor>
const char* lr_type_name() {
  return std::is_same<tensor, complex_tensor>::value ? "complex" : "real";
}

template <class tensor>
double lr_max_abs_entry(const tensor& a) {
  double m = 0.0;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    m = std::max(m, std::abs(a[n]));
  }
  return m;
}

// |got - want| <= rtol * max(|want|, scale). scale is a magnitude that does
// not cancel (the sum of the absolute channel contributions, or the largest
// operator element), so that a small result is not judged by a tolerance
// proportional to itself.
template <class V>
bool lr_check_close(const std::string& label, V got, V want, double scale,
                    double rtol) {
  const double tol = rtol * std::max(std::abs(want), scale);
  const double diff = std::abs(got - want);
  INFO(label << ": got=" << got << " want=" << want << " |diff|=" << diff
             << " tol=" << tol);
  CHECK(diff <= tol);
  return diff <= tol;
}

template <class tensor>
void lr_check_allclose(const tensor& got, const tensor& want,
                       const std::string& label, double rtol) {
  {
    INFO(label << ": shape mismatch");
    REQUIRE(got.shape() == want.shape());
  }
  const double scale =
      std::max(1.0, std::max(lr_max_abs_entry(got), lr_max_abs_entry(want)));
  double max_dev = 0.0;
  for (std::size_t n = 0; n < want.local_size(); ++n) {
    const mptensor::Index idx = want.global_index(n);
    typename tensor::value_type w, g;
    want.get_value(idx, w);
    got.get_value(idx, g);
    max_dev = std::max(max_dev, std::abs(g - w));
  }
  INFO(label << ": max |got-want| = " << max_dev << " (tol " << rtol * scale
             << ")");
  CHECK(max_dev <= rtol * scale);
}

inline rf::parity_vector lr_phys(int d) {
  if (d == 2) {
    return {false, true};
  }
  if (d == 4) {
    return {false, true, true, false};
  }
  throw std::runtime_error("lr_phys: unsupported physical dimension");
}

// ---- one-site operators as element tables -----------------------------------

using lr_cplx = std::complex<double>;

// m[in * d + out] = <out|A|in>.
struct lr_onesite {
  int d = 2;
  std::vector<lr_cplx> m;
  bool odd = false;
};

inline lr_onesite lr_onesite_op(const std::string& name, int d) {
  lr_onesite o;
  o.d = d;
  o.m.assign(d * d, 0.0);
  auto set = [&](int in, int out, double v) { o.m[in * d + out] = v; };
  if (name == "I") {
    for (int i = 0; i < d; ++i) {
      set(i, i, 1.0);
    }
    return o;
  }
  if (d == 2) {
    if (name == "c") {
      set(1, 0, 1.0);
      o.odd = true;
    } else if (name == "cdag") {
      set(0, 1, 1.0);
      o.odd = true;
    } else if (name == "n") {
      set(1, 1, 1.0);
    } else {
      throw std::runtime_error("lr_onesite_op: unknown d=2 operator " + name);
    }
    return o;
  }
  if (d != 4) {
    throw std::runtime_error("lr_onesite_op: unsupported dimension");
  }
  if (name == "c_up") {
    set(1, 0, 1.0);  // c_up |up> = |0>
    set(3, 2, 1.0);  // c_up |updn> = |dn>
    o.odd = true;
  } else if (name == "c_dn") {
    set(2, 0, 1.0);   // c_dn |dn> = |0>
    set(3, 1, -1.0);  // c_dn |updn> = -|up>
    o.odd = true;
  } else if (name == "cdag_up") {
    set(0, 1, 1.0);
    set(2, 3, 1.0);
    o.odd = true;
  } else if (name == "cdag_dn") {
    set(0, 2, 1.0);
    set(1, 3, -1.0);
    o.odd = true;
  } else if (name == "n") {
    set(1, 1, 1.0);
    set(2, 2, 1.0);
    set(3, 3, 2.0);
  } else if (name == "pair") {
    // |updn><0| + |0><updn|: even, off-diagonal inside the even sector.
    set(0, 3, 1.0);
    set(3, 0, 1.0);
  } else if (name == "flip") {
    // c+_up c_dn: |dn> -> |up>, even, off-diagonal inside the odd sector.
    set(2, 1, 1.0);
  } else {
    throw std::runtime_error("lr_onesite_op: unknown d=4 operator " + name);
  }
  return o;
}

inline void lr_set_scalar(real_tensor& t, const mptensor::Index& idx,
                          lr_cplx v) {
  REQUIRE(v.imag() == 0.0);
  t.set_value(idx, v.real());
}

inline void lr_set_scalar(complex_tensor& t, const mptensor::Index& idx,
                          lr_cplx v) {
  t.set_value(idx, v);
}

template <class tensor>
tensor lr_onesite_tensor(const lr_onesite& o, lr_cplx coef = 1.0) {
  tensor t(mptensor::Shape(o.d, o.d));
  for (int in = 0; in < o.d; ++in) {
    for (int out = 0; out < o.d; ++out) {
      const lr_cplx v = coef * o.m[in * o.d + out];
      if (v != 0.0) {
        lr_set_scalar(t, mptensor::Index(in, out), v);
      }
    }
  }
  return t;
}

// ---- two-site operators written out by the test -----------------------------

// coef * A_s B_t.
struct lr_term {
  lr_cplx coef;
  std::string a;
  std::string b;
};

// op4[i_s, i_t, o_s, o_t] = sum coef (-1)^{p_B p(i_s)} A[i_s, o_s]
// B[i_t, o_t], flattened ((i_s * dt + i_t) * ds + o_s) * dt + o_t.
inline std::vector<lr_cplx> lr_product_el(const std::vector<lr_term>& terms,
                                          int ds, int dt) {
  const rf::parity_vector ps = lr_phys(ds);
  std::vector<lr_cplx> el(ds * dt * ds * dt, 0.0);
  for (const lr_term& term : terms) {
    const lr_onesite A = lr_onesite_op(term.a, ds);
    const lr_onesite B = lr_onesite_op(term.b, dt);
    for (int is = 0; is < ds; ++is) {
      for (int it = 0; it < dt; ++it) {
        for (int os = 0; os < ds; ++os) {
          for (int ot = 0; ot < dt; ++ot) {
            const double sign = (B.odd && ps[is]) ? -1.0 : 1.0;
            el[((is * dt + it) * ds + os) * dt + ot] +=
                term.coef * sign * A.m[is * ds + os] * B.m[it * dt + ot];
          }
        }
      }
    }
  }
  return el;
}

template <class tensor>
tensor lr_twosite_tensor(const std::vector<lr_term>& terms, int ds, int dt) {
  const auto el = lr_product_el(terms, ds, dt);
  tensor op(mptensor::Shape(ds, dt, ds, dt));
  for (int is = 0; is < ds; ++is) {
    for (int it = 0; it < dt; ++it) {
      for (int os = 0; os < ds; ++os) {
        for (int ot = 0; ot < dt; ++ot) {
          const lr_cplx v = el[((is * dt + it) * ds + os) * dt + ot];
          if (v != 0.0) {
            lr_set_scalar(op, mptensor::Index(is, it, os, ot), v);
          }
        }
      }
    }
  }
  return op;
}

// Complex coefficients of the non-Hermitian operators: a conj() on the wrong
// layer or a source / target swap changes their value.
const lr_cplx lr_alpha(0.8, 0.6);
const lr_cplx lr_beta(-0.3, 0.5);
const lr_cplx lr_gamma(0.4, -0.7);
const lr_cplx lr_delta(-0.25, 0.35);

// Operator kinds of the solver cases:
//   cdagc  c+_s c_t (d = 2), odd x odd only
//   hop    sum_sigma (c+_s c_t + c+_t c_s), c+_t c_s = -c_s c+_t
//   nn     n_s n_t
//   hopnn  hop + nn
//   asym   c+_s c_t - 0.6 c_s c+_t + 0.8 n_s 1_t + 0.35 n_s n_t (d = 2):
//          even and odd channels, and not symmetric under s <-> t
//   asym4  c+_up,s c_up,t - 0.7 c_dn,s c+_dn,t + 0.4 n_s 1_t (d = 4)
//   n1mn   n_s (1 - n_t) = n_s 1_t - n_s n_t
//   cplx   sum_sigma (alpha c+_s c_t - beta c_s c+_t) + gamma n_s n_t
//          + delta n_s 1_t
inline std::vector<lr_term> lr_kind_terms(const std::string& kind, int d) {
  std::vector<std::pair<std::string, std::string>> flavours;
  if (d == 2) {
    flavours.push_back({"cdag", "c"});
  } else {
    flavours.push_back({"cdag_up", "c_up"});
    flavours.push_back({"cdag_dn", "c_dn"});
  }
  std::vector<lr_term> t;
  if (kind == "cdagc") {
    t.push_back({1.0, flavours[0].first, flavours[0].second});
  } else if (kind == "hop" || kind == "hopnn") {
    for (const auto& f : flavours) {
      t.push_back({1.0, f.first, f.second});
      t.push_back({-1.0, f.second, f.first});
    }
    if (kind == "hopnn") {
      t.push_back({1.0, "n", "n"});
    }
  } else if (kind == "nn") {
    t.push_back({1.0, "n", "n"});
  } else if (kind == "asym") {
    t.push_back({1.0, "cdag", "c"});
    t.push_back({-0.6, "c", "cdag"});
    t.push_back({0.8, "n", "I"});
    t.push_back({0.35, "n", "n"});
  } else if (kind == "asym4") {
    t.push_back({1.0, "cdag_up", "c_up"});
    t.push_back({-0.7, "c_dn", "cdag_dn"});
    t.push_back({0.4, "n", "I"});
  } else if (kind == "n1mn") {
    t.push_back({1.0, "n", "I"});
    t.push_back({-1.0, "n", "n"});
  } else if (kind == "cplx") {
    for (const auto& f : flavours) {
      t.push_back({lr_alpha, f.first, f.second});
      t.push_back({-lr_beta, f.second, f.first});
    }
    t.push_back({lr_gamma, "n", "n"});
    t.push_back({lr_delta, "n", "I"});
  } else {
    throw std::runtime_error("lr_kind_terms: unknown kind " + kind);
  }
  return t;
}

template <class tensor>
tensor lr_kind_tensor(const std::string& kind, int d) {
  return lr_twosite_tensor<tensor>(lr_kind_terms(kind, d), d, d);
}

}  // namespace

// ============================================================================
// Contract item 2: operator_parity
// ============================================================================

namespace {

template <class tensor>
void lr_run_operator_parity() {
  using P = rf::op_parity;
  INFO("tensor type " << lr_type_name<tensor>());
  const auto op = [](int d, const std::vector<lr_term>& sum_of_ones) {
    // A one-site operator as a sum of coef * named operator (b ignored).
    lr_onesite acc;
    acc.d = d;
    acc.m.assign(d * d, 0.0);
    for (const lr_term& t : sum_of_ones) {
      const lr_onesite o = lr_onesite_op(t.a, d);
      for (int k = 0; k < d * d; ++k) {
        acc.m[k] += t.coef * o.m[k];
      }
    }
    return lr_onesite_tensor<tensor>(acc);
  };
  const rf::parity_vector p2 = lr_phys(2);
  const rf::parity_vector p4 = lr_phys(4);

  // d = 2.
  CHECK(rf::operator_parity(op(2, {{1.0, "n", ""}}), p2) == P::even);
  CHECK(rf::operator_parity(op(2, {{1.0, "I", ""}, {-0.5, "n", ""}}), p2) ==
        P::even);
  CHECK(rf::operator_parity(op(2, {{1.0, "c", ""}}), p2) == P::odd);
  CHECK(rf::operator_parity(op(2, {{1.0, "c", ""}, {0.3, "cdag", ""}}), p2) ==
        P::odd);
  CHECK(rf::operator_parity(op(2, {{1.0, "c", ""}, {0.3, "n", ""}}), p2) ==
        P::mixed);
  CHECK(rf::operator_parity(tensor(mptensor::Shape(2, 2)), p2) == P::even);
  // The parity is read from the ledger, not from the index: with the
  // odd-first ledger [o, e] the diagonal is still even and c still odd,
  // and an element (0, 0) + (0, 1) is still mixed.
  const rf::parity_vector odd_first{true, false};
  CHECK(rf::operator_parity(op(2, {{1.0, "n", ""}}), odd_first) == P::even);
  CHECK(rf::operator_parity(op(2, {{1.0, "cdag", ""}}), odd_first) == P::odd);
  CHECK(rf::operator_parity(op(2, {{1.0, "I", ""}, {1.0, "cdag", ""}}),
                            odd_first) == P::mixed);
  // A ledger [e, e] makes every operator even (no index is odd).
  CHECK(rf::operator_parity(op(2, {{1.0, "c", ""}}),
                            rf::parity_vector{false, false}) == P::even);

  // d = 4, ledger [e, o, o, e].
  CHECK(rf::operator_parity(op(4, {{1.0, "n", ""}}), p4) == P::even);
  CHECK(rf::operator_parity(op(4, {{1.0, "pair", ""}}), p4) == P::even);
  CHECK(rf::operator_parity(op(4, {{1.0, "flip", ""}}), p4) == P::even);
  CHECK(rf::operator_parity(op(4, {{1.0, "c_up", ""}}), p4) == P::odd);
  CHECK(rf::operator_parity(op(4, {{1.0, "cdag_dn", ""}, {0.7, "c_up", ""}}),
                            p4) == P::odd);
  CHECK(rf::operator_parity(op(4, {{1.0, "c_dn", ""}, {0.2, "n", ""}}), p4) ==
        P::mixed);
  CHECK(rf::operator_parity(op(4, {{1.0, "flip", ""}, {1.0, "c_up", ""}}),
                            p4) == P::mixed);
  CHECK(rf::operator_parity(tensor(mptensor::Shape(4, 4)), p4) == P::even);
  // A single mixing element among many even ones is enough.
  {
    lr_onesite o = lr_onesite_op("n", 4);
    o.m[3 * 4 + 1] = 1.0e-3;  // |updn> -> |up>: odd
    CHECK(rf::operator_parity(lr_onesite_tensor<tensor>(o), p4) == P::mixed);
  }
}

}  // namespace

TEST_CASE(
    "longrange T2-2: operator_parity classifies even, odd, mixed and zero "
    "one-site operators (d = 2 and d = 4)") {
  lr_run_operator_parity<real_tensor>();
  lr_run_operator_parity<complex_tensor>();
}

// ============================================================================
// Contract item 3: product_twosite_op
// ============================================================================

TEST_CASE(
    "longrange [truth] the test's product table reproduces the "
    "Jordan-Wigner hopping tables of fold_geometry.cpp") {
  // d = 2: c+_1 c_2 + c+_2 c_1 has exactly (0,1,1,0) = (1,0,0,1) = 1.
  {
    const auto hop = lr_kind_tensor<real_tensor>("hop", 2);
    real_tensor want(mptensor::Shape(2, 2, 2, 2));
    want.set_value(mptensor::Index(0, 1, 1, 0), 1.0);
    want.set_value(mptensor::Index(1, 0, 0, 1), 1.0);
    lr_check_allclose(hop, want, "d=2 hop", 1.0e-15);
  }
  // d = 4: c_A = c x 1, c_B = P x c with P = (-1)^n on site A.
  {
    real_tensor want(mptensor::Shape(4, 4, 4, 4));
    double c[2][4][4] = {};  // c[sigma][out][in]
    c[0][0][1] = 1.0;
    c[0][2][3] = 1.0;
    c[1][0][2] = 1.0;
    c[1][1][3] = -1.0;
    const double ps[4] = {1.0, -1.0, -1.0, 1.0};
    for (int i1 = 0; i1 < 4; ++i1) {
      for (int i2 = 0; i2 < 4; ++i2) {
        for (int o1 = 0; o1 < 4; ++o1) {
          for (int o2 = 0; o2 < 4; ++o2) {
            double v = 0.0;
            for (int s = 0; s < 2; ++s) {
              v += c[s][i1][o1] * ps[i1] * c[s][o2][i2];
              v += ps[o1] * c[s][o1][i1] * c[s][i2][o2];
            }
            if (v != 0.0) {
              want.set_value(mptensor::Index(i1, i2, o1, o2), v);
            }
          }
        }
      }
    }
    lr_check_allclose(lr_kind_tensor<real_tensor>("hop", 4), want, "d=4 hop",
                      1.0e-15);
  }
}

namespace {

template <class tensor>
void lr_run_product_elements(const std::string& a, const std::string& b, int ds,
                             int dt, lr_cplx coef_a) {
  const std::string label = std::string("product_twosite_op(") + a + ", " + b +
                            ") ds=" + std::to_string(ds) +
                            " dt=" + std::to_string(dt) + " " +
                            lr_type_name<tensor>();
  INFO(label);
  lr_onesite A = lr_onesite_op(a, ds);
  const lr_onesite B = lr_onesite_op(b, dt);
  const tensor At = lr_onesite_tensor<tensor>(A, coef_a);
  const tensor Bt = lr_onesite_tensor<tensor>(B);
  const tensor got =
      rf::product_twosite_op(At, Bt, lr_phys(ds), lr_phys(dt), B.odd);
  const tensor want = lr_twosite_tensor<tensor>({{coef_a, a, b}}, ds, dt);
  lr_check_allclose(got, want, label, 1.0e-15);
}

}  // namespace

TEST_CASE(
    "longrange T2-3a: product_twosite_op is the contract's table "
    "(-1)^{p_B p(i_s)} A[i_s,o_s] B[i_t,o_t], elementwise") {
  // d = 2, all four odd x odd products and the even ones.
  for (const char* a : {"cdag", "c"}) {
    for (const char* b : {"cdag", "c"}) {
      lr_run_product_elements<real_tensor>(a, b, 2, 2, 1.0);
    }
  }
  lr_run_product_elements<real_tensor>("n", "n", 2, 2, 1.0);
  lr_run_product_elements<real_tensor>("n", "c", 2, 2, 1.0);
  lr_run_product_elements<real_tensor>("c", "n", 2, 2, 1.0);
  // d = 4: the sign is taken from the source ledger [e, o, o, e].
  lr_run_product_elements<real_tensor>("cdag_up", "c_up", 4, 4, 1.0);
  lr_run_product_elements<real_tensor>("c_dn", "cdag_dn", 4, 4, 1.0);
  lr_run_product_elements<real_tensor>("cdag_dn", "c_up", 4, 4, 1.0);
  lr_run_product_elements<real_tensor>("n", "n", 4, 4, 1.0);
  lr_run_product_elements<real_tensor>("pair", "flip", 4, 4, 1.0);
  // Different physical dimensions on the two ends: phys_s and phys_t must
  // each be used for its own leg.
  lr_run_product_elements<real_tensor>("c", "cdag_up", 2, 4, 1.0);
  lr_run_product_elements<real_tensor>("c_dn", "cdag", 4, 2, 1.0);
  // Complex entries are copied, not conjugated.
  lr_run_product_elements<complex_tensor>("c", "cdag", 2, 2, lr_alpha);
  lr_run_product_elements<complex_tensor>("c_dn", "cdag_up", 4, 4, lr_beta);
}

// ---- contract item 3b: open patch against the Fock oracle -------------------
//
// Same framework as test/fermion/relay_window.cpp (T1): a patch "NROW x
// NCOL" in the Contract_* orientation (row 0 on top), site s = col + ncol *
// row = fock_oracle's x + lx * y, every perimeter leg of dimension 1 and even,
// every internal bond [e, o], physical [e, o], and the deterministic
// parity-even site tensors of fock_oracle.deterministic_tensor (seed 0). The
// relay value is the exact contraction of build_relay_window() summed over
// the channels, divided by the exact contraction of the identity window.
//
// Literals generated from test/fermion/fock_oracle.py (unchanged), run in a
// copy of test/fermion/ with
//
//   from fock_oracle import Oracle, make_case
//   for lx, ly, pairs, dd in [
//           (3, 1, [(0, 1), (1, 0), (0, 2), (2, 0), (2, 1)], [(0, 1), (0, 2)]),
//           (2, 2, [(0, 1), (1, 0), (2, 1), (1, 2)], [(2, 1)])]:
//       patch, tensors, lp = make_case(lx, ly, [False, True], 0)
//       o = Oracle(patch, tensors, lp); n = o.norm()
//       print(n, [o.one_body(s, t) / n for s, t in pairs],
//             [o.density_density(s, t) / n for s, t in dd])
//
// one_body(i, j) = <c+_i c_j>, density_density(i, j) = <n_i n_j> (both times
// the norm). The states are real, so one_body is symmetric; the product
// c_s c+_t is then -<c+_t c_s> = -one_body(s, t). (A 1x2 patch is useless
// here: its state is a|00> + b|11>, whose hopping vanishes identically.)

namespace {

struct lp_patch {
  int nrow;
  int ncol;
  int nsite() const { return nrow * ncol; }
  int site(int row, int col) const { return col + ncol * row; }
  int row_of(int s) const { return s / ncol; }
  int col_of(int s) const { return s % ncol; }
  // Neighbour across virtual leg 0 = l, 1 = t, 2 = r, 3 = b; -1 if outside.
  int neighbour(int s, int leg) const {
    const int r = row_of(s);
    const int c = col_of(s);
    switch (leg) {
      case 0:
        return c > 0 ? s - 1 : -1;
      case 1:
        return r > 0 ? s - ncol : -1;
      case 2:
        return c + 1 < ncol ? s + 1 : -1;
      case 3:
        return r + 1 < nrow ? s + ncol : -1;
      default:
        return -1;
    }
  }
  // Leg labels site * 5 + leg; partner label across an internal bond.
  int partner(int label) const {
    const int s = label / 5;
    const int leg = label % 5;
    if (leg == 4) {
      return -1;
    }
    const int n = neighbour(s, leg);
    return n < 0 ? -1 : n * 5 + (leg + 2) % 4;
  }
};

// fock_oracle.deterministic_tensor(site, parities, seed).
inline double lp_det_entry(int site, int seed, const mptensor::Index& idx) {
  double x = 1.0 + seed;
  for (int ax = 0; ax < 5; ++ax) {
    x += static_cast<double>((ax + 3 + seed % 5) * idx[ax]);
  }
  x *= site + 2;
  return 0.19 * std::sin(x) + 0.13 * std::cos(0.37 * x);
}

inline lr_ft<real_tensor> lp_make_site(const lp_patch& p, int s) {
  const rf::parity_vector edge{false};
  const rf::parity_vector eo{false, true};
  rf::leg_parities lp;
  for (int leg = 0; leg < 4; ++leg) {
    lp.push_back(p.neighbour(s, leg) >= 0 ? eo : edge);
  }
  lp.push_back(eo);
  real_tensor t(mptensor::Shape(lp[0].size(), lp[1].size(), lp[2].size(),
                                lp[3].size(), lp[4].size()));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    if (rf::count_odd(lp, idx) % 2 == 0) {
      t.set_value(idx, lp_det_entry(s, 0, idx));
    }
  }
  return lr_ft<real_tensor>{t, lp};
}

// Exact contraction of a folded window: physical pair of each rank-6 cell
// traced (the identity operator the kernels are given), fused bonds
// contracted plainly.
inline double lp_contract_window(const lp_patch& p,
                                 const std::vector<std::vector<real_tensor>>& w,
                                 const std::string& label) {
  INFO(label);
  REQUIRE(w.size() == static_cast<std::size_t>(p.nrow));
  for (int r = 0; r < p.nrow; ++r) {
    REQUIRE(w[r].size() == static_cast<std::size_t>(p.ncol));
    for (int c = 0; c < p.ncol; ++c) {
      REQUIRE(w[r][c].shape().size() == 6);
    }
  }
  for (int s = 0; s < p.nsite(); ++s) {
    const auto& t = w[p.row_of(s)][p.col_of(s)];
    for (int leg = 0; leg < 4; ++leg) {
      const int n = p.neighbour(s, leg);
      if (n < 0) {
        REQUIRE(t.shape()[leg] == 1);
      } else {
        REQUIRE(t.shape()[leg] ==
                w[p.row_of(n)][p.col_of(n)].shape()[(leg + 2) % 4]);
      }
    }
  }
  real_tensor acc;
  std::vector<int> labels;
  for (int s = 0; s < p.nsite(); ++s) {
    const real_tensor traced = mptensor::contract(
        w[p.row_of(s)][p.col_of(s)], mptensor::Axes(4), mptensor::Axes(5));
    const std::vector<int> node{s * 5 + 0, s * 5 + 1, s * 5 + 2, s * 5 + 3};
    if (s == 0) {
      acc = traced;
      labels = node;
      continue;
    }
    mptensor::Axes axes_a, axes_b;
    std::vector<bool> used_a(labels.size(), false);
    std::vector<bool> used_b(node.size(), false);
    for (std::size_t j = 0; j < node.size(); ++j) {
      const int partner = p.partner(node[j]);
      const auto it = std::find(labels.begin(), labels.end(), partner);
      if (partner < 0 || it == labels.end()) {
        continue;
      }
      const std::size_t k = static_cast<std::size_t>(it - labels.begin());
      axes_a.push(static_cast<int>(k));
      axes_b.push(static_cast<int>(j));
      used_a[k] = true;
      used_b[j] = true;
    }
    std::vector<int> next;
    for (std::size_t k = 0; k < labels.size(); ++k) {
      if (!used_a[k]) {
        next.push_back(labels[k]);
      }
    }
    for (std::size_t j = 0; j < node.size(); ++j) {
      if (!used_b[j]) {
        next.push_back(node[j]);
      }
    }
    acc = mptensor::tensordot(acc, traced, axes_a, axes_b);
    labels = next;
  }
  mptensor::Index idx;
  idx.resize(acc.shape().size());
  for (std::size_t ax = 0; ax < acc.shape().size(); ++ax) {
    REQUIRE(acc.shape()[ax] == 1);
    idx[ax] = 0;
  }
  double v = 0.0;
  acc.get_value(idx, v);
  return v;
}

struct lp_state {
  lp_patch patch;
  std::vector<std::vector<lr_ft<real_tensor>>> grid;
  double norm_window = 0.0;
};

inline lp_state lp_make_state(int nrow, int ncol) {
  lp_state st;
  st.patch = lp_patch{nrow, ncol};
  st.grid.resize(nrow);
  std::vector<std::vector<real_tensor>> w(nrow);
  for (int r = 0; r < nrow; ++r) {
    for (int c = 0; c < ncol; ++c) {
      st.grid[r].push_back(lp_make_site(st.patch, st.patch.site(r, c)));
      w[r].push_back(rf::build_reduced_op(st.grid[r][c]));
    }
  }
  st.norm_window = lp_contract_window(st.patch, w, "identity window");
  return st;
}

// sum_k [exact window k] / [exact identity window] for a plain two-site
// operator loaded with wrap_twosite_gate().
inline double lp_relay_value(const lp_state& st, rf::window_cell src,
                             rf::window_cell tgt, const real_tensor& op4,
                             const std::string& label) {
  const rf::parity_vector eo{false, true};
  const auto channels = rf::relay_channels(rf::wrap_twosite_gate(op4, eo, eo));
  REQUIRE(!channels.empty());
  double sum = 0.0;
  for (std::size_t k = 0; k < channels.size(); ++k) {
    sum += lp_contract_window(
        st.patch, rf::build_relay_window(st.grid, src, tgt, channels[k]),
        label + " channel " + std::to_string(k));
  }
  return sum / st.norm_window;
}

struct lp_anchor {
  const char* name;
  int nrow;
  int ncol;
  rf::window_cell src;
  rf::window_cell tgt;
  const char* a;  // source operator
  const char* b;  // target operator
  double value;   // Fock-oracle <A_s B_t>
};

const double lp_norm_1x3 = 2.67265746667580693e-05;
const double lp_norm_2x2 = 9.97335335163249894e-06;

// make_case(3, 1): one_body(0,1) = one_body(1,0) = 2.05141392064624936e-02,
// one_body(0,2) = one_body(2,0) = -1.10806860734770890e-01, one_body(2,1) =
// 3.53504037997332812e-03, density_density(0,1) = 1.90944754314671984e-02,
// density_density(0,2) = 6.54456862818448941e-04. make_case(2, 2):
// one_body(2,1) = one_body(1,2) = 1.42176220745016541e-01.
const lp_anchor lp_anchors[] = {
    {"(1,0) 1x3 <c+_0 c_1> = one_body(0,1)",
     1,
     3,
     {0, 0},
     {0, 1},
     "cdag",
     "c",
     2.05141392064624936e-02},
    {"(1,0) 1x3 <c_0 c+_1> = -one_body(1,0)",
     1,
     3,
     {0, 0},
     {0, 1},
     "c",
     "cdag",
     -2.05141392064624936e-02},
    {"(1,0) 1x3 <n_0 n_1> = density_density(0,1)",
     1,
     3,
     {0, 0},
     {0, 1},
     "n",
     "n",
     1.90944754314671984e-02},
    {"(2,0) 1x3 <c+_0 c_2> = one_body(0,2)",
     1,
     3,
     {0, 0},
     {0, 2},
     "cdag",
     "c",
     -1.10806860734770890e-01},
    {"(2,0) 1x3 <c_0 c+_2> = -one_body(2,0)",
     1,
     3,
     {0, 0},
     {0, 2},
     "c",
     "cdag",
     1.10806860734770890e-01},
    {"(2,0) 1x3 <n_0 n_2> = density_density(0,2)",
     1,
     3,
     {0, 0},
     {0, 2},
     "n",
     "n",
     6.54456862818448941e-04},
    {"extra (-1,0) 1x3 <c+_2 c_1> = one_body(2,1)",
     1,
     3,
     {0, 2},
     {0, 1},
     "cdag",
     "c",
     3.53504037997332812e-03},
    {"extra (-1,0) 1x3 <c_2 c+_1> = -one_body(1,2)",
     1,
     3,
     {0, 2},
     {0, 1},
     "c",
     "cdag",
     -3.53504037997332812e-03},
    {"extra (1,1) 2x2 <c+_2 c_1> = one_body(2,1)",
     2,
     2,
     {1, 0},
     {0, 1},
     "cdag",
     "c",
     1.42176220745016541e-01},
    {"extra (1,1) 2x2 <c_2 c+_1> = -one_body(1,2)",
     2,
     2,
     {1, 0},
     {0, 1},
     "c",
     "cdag",
     -1.42176220745016541e-01},
};

}  // namespace

TEST_CASE(
    "longrange [truth] the open-patch identity windows reproduce the Fock "
    "oracle norms") {
  const lp_state s13 = lp_make_state(1, 3);
  lr_check_close("1x3 identity window vs oracle norm", s13.norm_window,
                 lp_norm_1x3, 0.0, 1.0e-12);
  const lp_state s22 = lp_make_state(2, 2);
  lr_check_close("2x2 identity window vs oracle norm", s22.norm_window,
                 lp_norm_2x2, 0.0, 1.0e-12);
  // The test's own product table fed through the same relay gives the
  // oracle values, so the anchors below judge product_twosite_op() alone.
  for (const lp_anchor& a : lp_anchors) {
    const lp_state st = lp_make_state(a.nrow, a.ncol);
    const real_tensor op4 =
        lr_twosite_tensor<real_tensor>({{1.0, a.a, a.b}}, 2, 2);
    lr_check_close(std::string(a.name) + " [test table through the relay]",
                   lp_relay_value(st, a.src, a.tgt, op4, a.name), a.value, 1.0,
                   1.0e-12);
  }
}

TEST_CASE(
    "longrange T2-3b: product_twosite_op through wrap_twosite_gate, "
    "relay_channels and build_relay_window equals the Fock oracle "
    "(<c+_s c_t>, -<c+_t c_s>, <n_s n_t>; (1,0) and (2,0))") {
  const rf::parity_vector eo{false, true};
  for (const lp_anchor& a : lp_anchors) {
    INFO(a.name);
    const lp_state st = lp_make_state(a.nrow, a.ncol);
    const lr_onesite A = lr_onesite_op(a.a, 2);
    const lr_onesite B = lr_onesite_op(a.b, 2);
    const real_tensor op4 = rf::product_twosite_op(
        lr_onesite_tensor<real_tensor>(A), lr_onesite_tensor<real_tensor>(B),
        eo, eo, B.odd);
    lr_check_close(std::string(a.name) + " [product_twosite_op vs oracle]",
                   lp_relay_value(st, a.src, a.tgt, op4, a.name), a.value, 1.0,
                   1.0e-12);
  }
}

// ============================================================================
// The solver fixture (contract items 4 to 11)
// ============================================================================

namespace {

constexpr int lr_D = 2;
constexpr int lr_chi = 4;
constexpr int lr_ctm_iteration_max = 30;
constexpr double lr_ctm_epsilon = 1.0e-10;
//! Odd virtual components are scaled by this per odd index: the CTM of a
//! fully random state converges poorly, and a vanishing odd component would
//! leave no signal in the odd channels.
constexpr double lr_odd_scale = 0.6;
constexpr unsigned lr_seed = 2027;

//! Contract items 5, 8, 10: relative tolerance of the wiring comparison.
constexpr double lr_rtol_window = 1.0e-12;
//! Contract item 6.
constexpr double lr_rtol_nn = 1.0e-10;

tenes::SquareLattice lr_lattice(int lx, int ly, int skew, int d, int D) {
  tenes::SquareLattice lattice(lx, ly, skew);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    lattice.physical_dims[site] = d;
    lattice.virtual_dims[site] = {D, D, D, D};
    lattice.initial_dirs[site] = {0.0};
    lattice.noises[site] = 1.0;
  }
  return lattice;
}

template <class tensor>
tenes::itps::PEPS_Parameters lr_params(int n_unit, int d, bool meanfield,
                                       int chi, int iteration_max,
                                       double epsilon) {
  tenes::itps::PEPS_Parameters p;
  p.fermion = true;
  p.is_real = std::is_same<tensor, real_tensor>::value;
  p.phys_parity.assign(n_unit, lr_phys(d));
  p.print_level = tenes::PrintLevel::none;
  p.outdir = "output_test_fermion_longrange";
  p.CHI = chi;
  p.Max_CTM_Iteration = iteration_max;
  p.CTM_Convergence_Epsilon = epsilon;
  p.Use_RSVD = false;
  p.MeanField_Env = meanfield;
  return p;
}

inline void lr_set_random(real_tensor& t, const mptensor::Index& idx,
                          std::mt19937& gen, double scale) {
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  t.set_value(idx, scale * dist(gen));
}

inline void lr_set_random(complex_tensor& t, const mptensor::Index& idx,
                          std::mt19937& gen, double scale) {
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  const double re = dist(gen);
  const double im = dist(gen);
  t.set_value(idx, scale * lr_cplx(re, im));
}

//! Random, parity-even, site-distinct Tn (site s draws from seed + 97 s),
//! every element multiplied by odd_scale per odd virtual index.
template <class tensor>
void lr_seed_Tn(lr_state<tensor>& state, double odd_scale, unsigned seed) {
  auto& Tn = lr_acc::Tn(state);
  const auto& fi = lr_acc::finfo(state);
  for (std::size_t s = 0; s < Tn.size(); ++s) {
    const rf::leg_parities parity = rf::Tn_parity(fi, static_cast<int>(s));
    mptensor::Shape sh;
    for (const auto& leg : parity) {
      sh.push(leg.size());
    }
    tensor t(sh);
    std::mt19937 gen(seed + 97u * static_cast<unsigned>(s));
    for (std::size_t n = 0; n < t.local_size(); ++n) {
      const mptensor::Index idx = t.global_index(n);
      double scale = 1.0;
      for (int leg = 0; leg < 4; ++leg) {
        if (parity[leg][idx[leg]]) {
          scale *= odd_scale;
        }
      }
      if (rf::count_odd(parity, idx) % 2 == 0) {
        lr_set_random(t, idx, gen, scale);
      }
    }
    REQUIRE(rf::parity_violation(lr_ft<tensor>{t, parity}) == 0.0);
    Tn[s] = t;
  }
}

//! Additive random noise on every environment tensor, of relative size eps:
//! t += eps * max|t| * u. Additive rather than multiplicative, so that
//! elements the converged CTM left exactly zero (a parity block, say) become
//! non-zero too.
template <class tensor>
void lr_perturb_env(lr_state<tensor>& state, double eps, unsigned seed) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  auto perturb = [&](std::vector<tensor>& ts) {
    for (tensor& t : ts) {
      const double scale = eps * lr_max_abs_entry(t);
      for (std::size_t n = 0; n < t.local_size(); ++n) {
        if constexpr (std::is_same<tensor, complex_tensor>::value) {
          const double re = dist(gen);
          const double im = dist(gen);
          t[n] += scale * lr_cplx(re, im);
        } else {
          t[n] += scale * dist(gen);
        }
      }
    }
  };
  perturb(lr_acc::C1(state));
  perturb(lr_acc::C2(state));
  perturb(lr_acc::C3(state));
  perturb(lr_acc::C4(state));
  perturb(lr_acc::eTt(state));
  perturb(lr_acc::eTr(state));
  perturb(lr_acc::eTb(state));
  perturb(lr_acc::eTl(state));
}

template <class tensor>
tensor lr_identity(int d) {
  return lr_onesite_tensor<tensor>(lr_onesite_op("I", d));
}

//! The reference ("direct") value of one window measurement.
template <class tensor>
struct lr_ref {
  typename tensor::value_type value = 0.0;  //!< coeff * sum_k c_k / norm
  double scale = 0.0;  //!< |coeff| * sum_k |c_k| / |norm| (no cancellation)
  typename tensor::value_type norm = 0.0;
  std::size_t nchannel = 0;
};

//! Deliberate miswirings of the direct value, used only by the [truth]
//! discrimination case (they show that the item-5 comparison can fail):
//!   y_first      the other relay path (must NOT change the value)
//!   mirror_edges the left / right (or top / bottom) edges of the opposite
//!                column (row)
//!   swap_ends    the relay started at the target and ended at the source
//!   drop_odd     odd channels left out
//!   drop_even    even channels left out
//!   flip_odd     odd channels with the opposite sign (a lost string sign)
enum class lr_mut {
  none,
  y_first,
  mirror_edges,
  swap_ends,
  drop_odd,
  drop_even,
  flip_odd
};

//! Contract item 5's direct computation: the window of twosite_obs.cpp's
//! bosonic CTM branch, the solver's environment tensors, T1's relay window
//! per channel and core::Contract_density_CTM. The norm is the same window
//! of plain reduced tensors. Never calls measure_twosite().
template <class tensor>
lr_ref<tensor> lr_direct(lr_state<tensor>& state, int source, int dx, int dy,
                         const tensor& op4, typename tensor::value_type coeff,
                         lr_mut mut = lr_mut::none) {
  const tenes::SquareLattice& lat = lr_acc::lattice(state);
  const auto& fi = lr_acc::finfo(state);
  const auto& Tn = lr_acc::Tn(state);
  const int ncol = std::abs(dx) + 1;
  const int nrow = std::abs(dy) + 1;
  REQUIRE(ncol <= 4);
  REQUIRE(nrow <= 4);
  const int scol = dx >= 0 ? 0 : ncol - 1;
  const int tcol = ncol - 1 - scol;
  const int srow = dy >= 0 ? nrow - 1 : 0;
  const int trow = nrow - 1 - srow;

  std::vector<std::vector<int>> idx(nrow, std::vector<int>(ncol));
  std::vector<std::vector<lr_ft<tensor>>> grid(nrow);
  std::vector<std::vector<tensor>> reduced(nrow);
  std::vector<std::vector<tensor>> ident(nrow);
  for (int r = 0; r < nrow; ++r) {
    for (int c = 0; c < ncol; ++c) {
      idx[r][c] = lat.other(source, c - scol, srow - r);
      grid[r].push_back(rf::wrap_Tn(Tn[idx[r][c]], fi, idx[r][c]));
      reduced[r].push_back(rf::build_reduced_op(grid[r][c]));
      ident[r].push_back(
          lr_identity<tensor>(static_cast<int>(fi.phys[idx[r][c]].size())));
    }
  }
  const int target = idx[trow][tcol];
  REQUIRE(target == lat.other(source, dx, dy));

  std::vector<const tensor*> C{&lr_acc::C1(state)[idx[0][0]],
                               &lr_acc::C2(state)[idx[0][ncol - 1]],
                               &lr_acc::C3(state)[idx[nrow - 1][ncol - 1]],
                               &lr_acc::C4(state)[idx[nrow - 1][0]]};
  std::vector<const tensor*> eTt(ncol), eTb(ncol), eTl(nrow), eTr(nrow);
  const bool mirror_h = mut == lr_mut::mirror_edges && ncol > 1;
  const bool mirror_v = mut == lr_mut::mirror_edges && ncol == 1;
  for (int c = 0; c < ncol; ++c) {
    eTt[c] = &lr_acc::eTt(state)[idx[mirror_v ? nrow - 1 : 0][c]];
    eTb[c] = &lr_acc::eTb(state)[idx[mirror_v ? 0 : nrow - 1][c]];
  }
  for (int r = 0; r < nrow; ++r) {
    eTl[r] = &lr_acc::eTl(state)[idx[r][mirror_h ? ncol - 1 : 0]];
    eTr[r] = &lr_acc::eTr(state)[idx[r][mirror_h ? 0 : ncol - 1]];
  }
  const auto pointers = [](const std::vector<std::vector<tensor>>& w) {
    std::vector<std::vector<const tensor*>> p(w.size());
    for (std::size_t r = 0; r < w.size(); ++r) {
      for (const tensor& t : w[r]) {
        p[r].push_back(&t);
      }
    }
    return p;
  };
  const auto op_ptr = pointers(ident);

  lr_ref<tensor> ref;
  ref.norm = tenes::itps::core::Contract_density_CTM(C, eTt, eTr, eTb, eTl,
                                                     pointers(reduced), op_ptr);
  REQUIRE(std::abs(ref.norm) > 0.0);

  const lr_ft<tensor> op12 =
      rf::wrap_twosite_gate(op4, fi.phys[source], fi.phys[target]);
  const auto channels = rf::relay_channels(op12);
  REQUIRE(!channels.empty());
  typename tensor::value_type sum = 0.0;
  double abs_sum = 0.0;
  rf::window_cell from{srow, scol};
  rf::window_cell to{trow, tcol};
  if (mut == lr_mut::swap_ends) {
    std::swap(from, to);
  }
  const rf::relay_order order = mut == lr_mut::y_first
                                    ? rf::relay_order::y_first
                                    : rf::relay_order::x_first;
  for (const auto& ch : channels) {
    const bool odd = ch.u.parity[2][0];
    if ((odd && mut == lr_mut::drop_odd) ||
        (!odd && mut == lr_mut::drop_even)) {
      continue;
    }
    const auto w = rf::build_relay_window(grid, from, to, ch, order);
    auto ck = tenes::itps::core::Contract_density_CTM(C, eTt, eTr, eTb, eTl,
                                                      pointers(w), op_ptr);
    if (odd && mut == lr_mut::flip_odd) {
      ck = -ck;
    }
    sum += ck;
    abs_sum += std::abs(ck);
  }
  ref.value = coeff * sum / ref.norm;
  ref.scale = std::abs(coeff) * abs_sum / std::abs(ref.norm);
  ref.nchannel = channels.size();
  return ref;
}

struct lr_disp {
  int dx;
  int dy;
};

inline std::string lr_bond_label(int source, int dx, int dy) {
  return "source " + std::to_string(source) + " (dx,dy)=(" +
         std::to_string(dx) + "," + std::to_string(dy) + ")";
}

//! One two-site observable group of a solver case.
struct lr_group {
  std::string kind;
  lr_cplx coeff;
};

template <class tensor>
typename tensor::value_type lr_value(lr_cplx v) {
  if constexpr (std::is_same<tensor, real_tensor>::value) {
    return v.real();
  } else {
    return v;
  }
}

//! A fermion iTPS on the given cell with one two-site group per entry of
//! `groups`, each measured at every (source, displacement) pair.
template <class tensor>
std::unique_ptr<lr_state<tensor>> lr_make_state(
    const tenes::SquareLattice& lattice, int d, bool meanfield,
    const tenes::Operators<tensor>& onesite,
    const tenes::Operators<tensor>& twosite, int chi = lr_chi,
    int iteration_max = lr_ctm_iteration_max, double epsilon = lr_ctm_epsilon) {
  return std::make_unique<lr_state<tensor>>(
      MPI_COMM_WORLD,
      lr_params<tensor>(lattice.N_UNIT, d, meanfield, chi, iteration_max,
                        epsilon),
      lattice, tenes::EvolutionOperators<tensor>{},
      tenes::EvolutionOperators<tensor>{}, onesite, twosite,
      tenes::Operators<tensor>{}, tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
}

template <class tensor>
tenes::Operators<tensor> lr_twosite_ops(const std::vector<lr_group>& groups,
                                        int d, const std::vector<int>& sources,
                                        const std::vector<lr_disp>& disps) {
  tenes::Operators<tensor> ops;
  for (std::size_t g = 0; g < groups.size(); ++g) {
    const tensor op = lr_kind_tensor<tensor>(groups[g].kind, d);
    for (const int s : sources) {
      for (const lr_disp& dd : disps) {
        ops.emplace_back(groups[g].kind, static_cast<int>(g), s, dd.dx, dd.dy,
                         op, lr_value<tensor>(groups[g].coeff));
      }
    }
  }
  return ops;
}

//! Contract items 5, 8, 10: measure_twosite() against the direct value, for
//! every (group, source, displacement). Returns the direct values.
template <class tensor>
std::map<std::tuple<int, int, int, int>, lr_ref<tensor>> lr_run_window_case(
    const std::string& label, int lx, int ly, int d,
    const std::vector<lr_group>& groups, const std::vector<int>& sources,
    const std::vector<lr_disp>& disps, double env_noise) {
  INFO(label);
  const tenes::SquareLattice lattice = lr_lattice(lx, ly, 0, d, lr_D);
  auto state =
      lr_make_state<tensor>(lattice, d, false, tenes::Operators<tensor>{},
                            lr_twosite_ops<tensor>(groups, d, sources, disps));
  lr_seed_Tn(*state, lr_odd_scale, lr_seed);
  state->update_CTM();
  if (env_noise > 0.0) {
    lr_perturb_env(*state, env_noise, lr_seed + 1);
  }
  std::map<std::tuple<int, int, int, int>, lr_ref<tensor>> refs;
  for (std::size_t g = 0; g < groups.size(); ++g) {
    const tensor op = lr_kind_tensor<tensor>(groups[g].kind, d);
    for (const int s : sources) {
      for (const lr_disp& dd : disps) {
        refs[{static_cast<int>(g), s, dd.dx, dd.dy}] = lr_direct(
            *state, s, dd.dx, dd.dy, op, lr_value<tensor>(groups[g].coeff));
      }
    }
  }
  std::vector<std::map<Bond, typename tensor::value_type>> measured;
  REQUIRE_NOTHROW(measured = state->measure_twosite());
  REQUIRE(measured.size() >= groups.size());
  for (const auto& [key, ref] : refs) {
    const auto [g, s, dx, dy] = key;
    const std::string what = label + " group " + groups[g].kind + " " +
                             lr_bond_label(s, dx, dy) + " [" +
                             std::to_string(ref.nchannel) + " channels]";
    INFO(what);
    REQUIRE(measured[g].count(Bond{s, dx, dy}) == 1);
    lr_check_close(what + " [measure_twosite vs direct]",
                   measured[g].at(Bond{s, dx, dy}), ref.value, ref.scale,
                   lr_rtol_window);
  }
  return refs;
}

const std::vector<lr_disp> lr_contract5_disps = {{2, 0},  {0, 2}, {1, 1},
                                                 {-1, 2}, {3, 0}, {2, -1}};

}  // namespace

// ============================================================================
// Contract item 4: one-site measurement
// ============================================================================

namespace {

struct lr_onesite_group {
  std::string name;
  lr_cplx coeff;
};

template <class tensor>
void lr_run_onesite(int d, const std::vector<lr_onesite_group>& groups,
                    double env_noise) {
  using value_type = typename tensor::value_type;
  const std::string label =
      std::string("one-site d=") + std::to_string(d) + " " +
      lr_type_name<tensor>() +
      (env_noise > 0.0 ? " perturbed environment" : " converged environment");
  INFO(label);
  const tenes::SquareLattice lattice = lr_lattice(2, 2, 0, d, lr_D);
  tenes::Operators<tensor> onesite;
  std::vector<rf::op_parity> expected_parity;
  for (std::size_t g = 0; g < groups.size(); ++g) {
    const lr_onesite o = lr_onesite_op(groups[g].name, d);
    for (int s = 0; s < lattice.N_UNIT; ++s) {
      onesite.emplace_back(groups[g].name, static_cast<int>(g), s,
                           lr_onesite_tensor<tensor>(o),
                           lr_value<tensor>(groups[g].coeff));
      expected_parity.push_back(o.odd ? rf::op_parity::odd
                                      : rf::op_parity::even);
    }
  }
  auto state = lr_make_state<tensor>(lattice, d, false, onesite,
                                     tenes::Operators<tensor>{});
  lr_seed_Tn(*state, lr_odd_scale, lr_seed + 5);
  state->update_CTM();
  if (env_noise > 0.0) {
    lr_perturb_env(*state, env_noise, lr_seed + 6);
  }

  // Direct values: build_reduced_op + Contract_one_site_density_CTM, as
  // onesite_obs.cpp does before this task.
  const auto& fi = lr_acc::finfo(*state);
  std::vector<std::vector<value_type>> want(groups.size());
  std::vector<std::vector<double>> want_scale(groups.size());
  double min_odd_raw = std::numeric_limits<double>::infinity();
  for (std::size_t g = 0; g < groups.size(); ++g) {
    const lr_onesite o = lr_onesite_op(groups[g].name, d);
    for (int s = 0; s < lattice.N_UNIT; ++s) {
      const tensor reduced =
          rf::build_reduced_op(rf::wrap_Tn(lr_acc::Tn(*state)[s], fi, s));
      const auto one_site = [&](const tensor& op) {
        return tenes::itps::core::Contract_one_site_density_CTM(
            lr_acc::C1(*state)[s], lr_acc::C2(*state)[s], lr_acc::C3(*state)[s],
            lr_acc::C4(*state)[s], lr_acc::eTt(*state)[s],
            lr_acc::eTr(*state)[s], lr_acc::eTb(*state)[s],
            lr_acc::eTl(*state)[s], reduced, op);
      };
      const value_type norm = one_site(lr_identity<tensor>(d));
      const value_type raw = one_site(lr_onesite_tensor<tensor>(o));
      want[g].push_back(lr_value<tensor>(groups[g].coeff) * raw / norm);
      want_scale[g].push_back(std::abs(groups[g].coeff) *
                              lr_max_abs_entry(lr_onesite_tensor<tensor>(o)));
      if (o.odd) {
        min_odd_raw = std::min(min_odd_raw, std::abs(raw / norm));
      }
    }
  }
  // Premise (perturbed environment only): the plain contraction of an odd
  // operator is not zero here, so an implementation that simply contracts
  // it, as the stub does, cannot pass by accident.
  if (env_noise > 0.0) {
    INFO("premise: smallest |plain contraction| of an odd operator "
         << min_odd_raw);
    REQUIRE(min_odd_raw > 1.0e-6);
  }

  const auto& parity = lr_acc::onesite_parity(*state);
  {
    INFO(
        "iTPS::onesite_parity is filled at construction, one entry per "
        "one-site operator");
    CHECK(parity.size() == expected_parity.size());
    for (std::size_t i = 0; i < std::min(parity.size(), expected_parity.size());
         ++i) {
      INFO("one-site operator " << i);
      CHECK(parity[i] == expected_parity[i]);
    }
  }

  const auto measured = state->measure_onesite();
  REQUIRE(measured.size() >= groups.size());
  for (std::size_t g = 0; g < groups.size(); ++g) {
    const lr_onesite o = lr_onesite_op(groups[g].name, d);
    REQUIRE(measured[g].size() == static_cast<std::size_t>(lattice.N_UNIT));
    for (int s = 0; s < lattice.N_UNIT; ++s) {
      const value_type v = measured[g][s];
      const std::string what =
          label + " operator " + groups[g].name + " site " + std::to_string(s);
      INFO(what << ": measured " << v << ", plain contraction " << want[g][s]);
      if (o.odd) {
        // Exactly zero, real and imaginary parts (no contraction).
        CHECK(std::real(v) == 0.0);
        CHECK(std::imag(v) == 0.0);
      } else {
        lr_check_close(what + " [even operator unchanged]", v, want[g][s],
                       want_scale[g][s], 1.0e-13);
      }
    }
  }
}

}  // namespace

TEST_CASE(
    "longrange T2-4a: an odd one-site operator measures exactly 0, even ones "
    "are unchanged (d = 2, converged and perturbed environment)") {
  const std::vector<lr_onesite_group> groups = {
      {"n", 1.0}, {"cdag", 1.0}, {"c", 1.3}, {"n", 0.7}};
  lr_run_onesite<real_tensor>(2, groups, 0.0);
  lr_run_onesite<real_tensor>(2, groups, 0.3);
  const std::vector<lr_onesite_group> cgroups = {
      {"n", lr_alpha}, {"cdag", 1.0}, {"c", lr_beta}};
  lr_run_onesite<complex_tensor>(2, cgroups, 0.3);
}

TEST_CASE(
    "longrange T2-4b: an odd one-site operator measures exactly 0, even ones "
    "are unchanged (d = 4)") {
  const std::vector<lr_onesite_group> groups = {{"n", 1.0},
                                                {"cdag_up", 1.0},
                                                {"c_dn", 0.9},
                                                {"pair", 1.0},
                                                {"flip", 1.0}};
  lr_run_onesite<real_tensor>(4, groups, 0.3);
}

// ============================================================================
// Contract items 5 and 10: the window measurement is wired as specified
// ============================================================================

TEST_CASE(
    "longrange T2-5a: measure_twosite equals the direct relay window, "
    "d = 2 real, 2x2 cell, converged CTM") {
  // cdagc carries a coefficient, asym mixes even and odd channels and is not
  // symmetric under exchanging its ends.
  lr_run_window_case<real_tensor>("T2-5a", 2, 2, 2,
                                  {{"cdagc", 1.7}, {"asym", 1.0}}, {0, 1, 2, 3},
                                  lr_contract5_disps, 0.0);
}

TEST_CASE(
    "longrange T2-5b: measure_twosite equals the direct relay window, "
    "d = 2 real, 2x2 cell, perturbed environment") {
  // An environment that is no CTM fixed point: every corner and edge differs
  // from every other, so a window closed with a wrong corner or edge (or a
  // mirrored window) cannot agree by symmetry.
  lr_run_window_case<real_tensor>("T2-5b", 2, 2, 2,
                                  {{"cdagc", 1.0}, {"asym", -0.8}},
                                  {0, 1, 2, 3}, lr_contract5_disps, 0.3);
}

TEST_CASE(
    "longrange T2-5c: measure_twosite equals the direct relay window, "
    "d = 4 real, 2x2 cell") {
  lr_run_window_case<real_tensor>("T2-5c converged", 2, 2, 4,
                                  {{"hopnn", 1.0}, {"asym4", 1.0}}, {0},
                                  lr_contract5_disps, 0.0);
  lr_run_window_case<real_tensor>("T2-5c perturbed", 2, 2, 4, {{"asym4", 1.0}},
                                  {3}, lr_contract5_disps, 0.3);
}

// The comparison of items 5, 8 and 10 must be able to fail: on the states
// those cases use, each miswiring of lr_mut moves the direct value by far
// more than the tolerance, and the other relay path (y_first) does not move
// it (so an implementation is free to choose either path).
TEST_CASE(
    "longrange [truth] T2-5 discrimination: miswired direct values differ "
    "from the direct value by far more than the item-5 tolerance") {
  const int d = 2;
  const tenes::SquareLattice lattice = lr_lattice(2, 2, 0, d, lr_D);
  for (const double noise : {0.0, 0.3}) {
    auto state = lr_make_state<real_tensor>(lattice, d, false,
                                            tenes::Operators<real_tensor>{},
                                            tenes::Operators<real_tensor>{});
    lr_seed_Tn(*state, lr_odd_scale, lr_seed);
    state->update_CTM();
    if (noise > 0.0) {
      lr_perturb_env(*state, noise, lr_seed + 1);
    }
    const real_tensor op = lr_kind_tensor<real_tensor>("asym", d);
    double min_sep[7];
    std::fill(std::begin(min_sep), std::end(min_sep),
              std::numeric_limits<double>::infinity());
    double max_path = 0.0;
    for (const int s : {0, 3}) {
      for (const lr_disp& dd : lr_contract5_disps) {
        const auto ref = lr_direct(*state, s, dd.dx, dd.dy, op, 1.0);
        const double tol =
            lr_rtol_window * std::max(std::abs(ref.value), ref.scale);
        for (const lr_mut m :
             {lr_mut::y_first, lr_mut::mirror_edges, lr_mut::swap_ends,
              lr_mut::drop_odd, lr_mut::drop_even, lr_mut::flip_odd}) {
          // In a 2x2 cell the mirrored column (row) holds the same sites
          // when |dx| (|dy|) is even: the miswiring is invisible there.
          const bool mirror_moves = std::abs(dd.dx) > 0
                                        ? std::abs(dd.dx) % 2 == 1
                                        : std::abs(dd.dy) % 2 == 1;
          if (m == lr_mut::mirror_edges && !mirror_moves) {
            continue;
          }
          const auto alt = lr_direct(*state, s, dd.dx, dd.dy, op, 1.0, m);
          const double sep = std::abs(alt.value - ref.value) / tol;
          if (m == lr_mut::y_first) {
            max_path = std::max(max_path, sep);
          } else {
            min_sep[static_cast<int>(m)] =
                std::min(min_sep[static_cast<int>(m)], sep);
          }
        }
      }
    }
    INFO("environment noise " << noise);
    std::cout << std::setprecision(3) << "longrange discrimination (noise "
              << noise << "): |y_first - x_first| / tol max " << max_path
              << "; min |miswired - direct| / tol: mirror_edges " << min_sep[2]
              << ", swap_ends " << min_sep[3] << ", drop_odd " << min_sep[4]
              << ", drop_even " << min_sep[5] << ", flip_odd " << min_sep[6]
              << std::endl;
    CHECK(max_path <= 1.0);
    for (int m = 2; m <= 6; ++m) {
      INFO("miswiring " << m);
      CHECK(min_sep[m] > 1.0e3);
    }
  }
}

TEST_CASE(
    "longrange T2-10: measure_twosite equals the direct relay window, "
    "complex_tensor (Review Focus 3)") {
  lr_run_window_case<complex_tensor>("T2-10 d=2 converged", 2, 2, 2,
                                     {{"cplx", 1.0}, {"cdagc", lr_gamma}},
                                     {0, 1, 2, 3}, lr_contract5_disps, 0.0);
  lr_run_window_case<complex_tensor>("T2-10 d=2 perturbed", 2, 2, 2,
                                     {{"cplx", 1.0}}, {0, 3},
                                     lr_contract5_disps, 0.3);
  lr_run_window_case<complex_tensor>("T2-10 d=4 perturbed", 2, 2, 4,
                                     {{"cplx", lr_delta}}, {0},
                                     {{1, 1}, {2, -1}, {-1, 2}}, 0.3);
}

// ============================================================================
// Contract item 8: windows larger than the unit cell (Review Focus 1)
// ============================================================================

TEST_CASE(
    "longrange T2-8: windows larger than the 2x2 unit cell, (2,0) onto the "
    "same sublattice and (3,1) (Review Focus 1)") {
  {
    const tenes::SquareLattice lattice = lr_lattice(2, 2, 0, 2, lr_D);
    for (int s = 0; s < lattice.N_UNIT; ++s) {
      INFO("premise: (2,0) lands on the source's own sublattice, site " << s);
      REQUIRE(lattice.other(s, 2, 0) == s);
    }
  }
  const std::vector<lr_disp> disps = {{2, 0}, {3, 1}};
  lr_run_window_case<real_tensor>("T2-8 d=2 converged", 2, 2, 2,
                                  {{"cdagc", 1.0}, {"asym", 1.0}}, {0, 1, 2, 3},
                                  disps, 0.0);
  lr_run_window_case<real_tensor>("T2-8 d=2 perturbed", 2, 2, 2,
                                  {{"asym", 1.0}}, {0, 1, 2, 3}, disps, 0.3);
  lr_run_window_case<real_tensor>("T2-8 d=4 perturbed", 2, 2, 4,
                                  {{"asym4", 1.0}}, {0, 3}, disps, 0.3);
}

// ============================================================================
// Contract item 11: mixed-parity channels (Review Focus 5)
// ============================================================================

TEST_CASE(
    "longrange T2-11: d = 4 hopping + nn at (2,0) equals hopping plus nn "
    "(Review Focus 5)") {
  const std::vector<lr_group> groups = {
      {"hopnn", 1.0}, {"hop", 1.0}, {"nn", 1.0}};
  const std::vector<int> sources = {0, 1, 2, 3};
  const tenes::SquareLattice lattice = lr_lattice(2, 2, 0, 4, lr_D);
  auto state = lr_make_state<real_tensor>(
      lattice, 4, false, tenes::Operators<real_tensor>{},
      lr_twosite_ops<real_tensor>(groups, 4, sources, {{2, 0}}));
  lr_seed_Tn(*state, lr_odd_scale, lr_seed + 11);
  state->update_CTM();
  // Premise: both parts carry signal far above the tolerance, so dropping
  // either the even or the odd channels of the mixed operator is visible.
  for (const int s : sources) {
    const auto hop =
        lr_direct(*state, s, 2, 0, lr_kind_tensor<real_tensor>("hop", 4), 1.0);
    const auto nn =
        lr_direct(*state, s, 2, 0, lr_kind_tensor<real_tensor>("nn", 4), 1.0);
    INFO("premise " << lr_bond_label(s, 2, 0) << ": direct <hop> = "
                    << hop.value << ", <nn> = " << nn.value);
    REQUIRE(std::abs(hop.value) > 1.0e-6);
    REQUIRE(std::abs(nn.value) > 1.0e-6);
  }
  std::vector<std::map<Bond, double>> measured;
  REQUIRE_NOTHROW(measured = state->measure_twosite());
  REQUIRE(measured.size() >= 3);
  for (const int s : sources) {
    const Bond b{s, 2, 0};
    REQUIRE(measured[0].count(b) == 1);
    REQUIRE(measured[1].count(b) == 1);
    REQUIRE(measured[2].count(b) == 1);
    const double sum = measured[1].at(b) + measured[2].at(b);
    lr_check_close(lr_bond_label(s, 2, 0) + " [<hop + nn> vs <hop> + <nn>]",
                   measured[0].at(b), sum,
                   std::abs(measured[1].at(b)) + std::abs(measured[2].at(b)),
                   1.0e-12);
  }
}

// ============================================================================
// Contract item 6: nearest neighbours
// ============================================================================

namespace {

const std::vector<lr_disp> lr_nn_disps = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};

template <class tensor>
void lr_run_nn_anchor(const std::string& label, int d,
                      const std::vector<lr_group>& groups, double env_noise,
                      const std::vector<int>& sources = {0, 1, 2, 3}) {
  INFO(label);
  const tenes::SquareLattice lattice = lr_lattice(2, 2, 0, d, lr_D);
  auto state = lr_make_state<tensor>(
      lattice, d, false, tenes::Operators<tensor>{},
      lr_twosite_ops<tensor>(groups, d, sources, lr_nn_disps));
  lr_seed_Tn(*state, lr_odd_scale, lr_seed + 21);
  state->update_CTM();
  if (env_noise > 0.0) {
    lr_perturb_env(*state, env_noise, lr_seed + 22);
  }
  const auto measured = state->measure_twosite();
  double max_rel = 0.0;
  for (std::size_t g = 0; g < groups.size(); ++g) {
    const tensor op = lr_kind_tensor<tensor>(groups[g].kind, d);
    for (const int s : sources) {
      for (const lr_disp& dd : lr_nn_disps) {
        const auto ref = lr_direct(*state, s, dd.dx, dd.dy, op,
                                   lr_value<tensor>(groups[g].coeff));
        const std::string what = label + " group " + groups[g].kind + " " +
                                 lr_bond_label(s, dd.dx, dd.dy);
        REQUIRE(measured[g].count(Bond{s, dd.dx, dd.dy}) == 1);
        const auto got = measured[g].at(Bond{s, dd.dx, dd.dy});
        max_rel =
            std::max(max_rel, std::abs(got - ref.value) /
                                  std::max(std::abs(ref.value), ref.scale));
        lr_check_close(what + " [existing nearest-neighbour path vs relay]",
                       got, ref.value, ref.scale, lr_rtol_nn);
      }
    }
  }
  std::cout << std::setprecision(3) << "longrange " << label
            << ": max relative |bundled-k - relay| " << max_rel << std::endl;
}

}  // namespace

// Contract item 6 compares T1's relay window, closed with the solver's
// environment exactly as the direct value of items 5, 8, 10 and 11 is, with
// the EXISTING nearest-neighbour path of measure_twosite() (bundled-k), which
// the task leaves unchanged. Neither side is T2 code, so this case passes
// against the stub by construction: it is the anchor that makes the direct
// value of the other cases trustworthy (window orientation, source / target
// cells, corner and edge assignment, norm and coefficient), on an
// independently verified path. The T2 check of item 6 proper is T2-6b.
TEST_CASE(
    "longrange [truth] T2-6a: at nearest neighbours the direct relay window "
    "equals the existing measure_twosite path") {
  lr_run_nn_anchor<real_tensor>("T2-6a d=2 converged", 2,
                                {{"hop", 1.0}, {"asym", 1.3}}, 0.0);
  lr_run_nn_anchor<real_tensor>("T2-6a d=2 perturbed", 2,
                                {{"hop", 1.0}, {"asym", 1.0}}, 0.3);
  lr_run_nn_anchor<real_tensor>("T2-6a d=4 converged", 4,
                                {{"hopnn", 1.0}, {"asym4", 1.0}}, 0.0, {0, 3});
  lr_run_nn_anchor<complex_tensor>("T2-6a d=2 complex perturbed", 2,
                                   {{"cplx", 1.0}}, 0.3);
}

TEST_CASE(
    "longrange T2-6b: nearest-neighbour values are unchanged when long-range "
    "observables are measured alongside") {
  // Two solvers with the same state and environment: one measures the
  // nearest-neighbour bonds alone (the existing path), the other the same
  // bonds together with long-range ones. The shared norm cache and the
  // group bookkeeping must not let the long-range windows leak into the
  // nearest-neighbour values.
  const int d = 2;
  const std::vector<lr_group> groups = {{"hop", 1.0}, {"asym", 1.0}};
  const std::vector<int> sources = {0, 1, 2, 3};
  std::vector<lr_disp> mixed = lr_nn_disps;
  for (const lr_disp& dd : lr_contract5_disps) {
    mixed.push_back(dd);
  }
  const tenes::SquareLattice lattice = lr_lattice(2, 2, 0, d, lr_D);
  auto nn_only = lr_make_state<real_tensor>(
      lattice, d, false, tenes::Operators<real_tensor>{},
      lr_twosite_ops<real_tensor>(groups, d, sources, lr_nn_disps));
  auto with_long = lr_make_state<real_tensor>(
      lattice, d, false, tenes::Operators<real_tensor>{},
      lr_twosite_ops<real_tensor>(groups, d, sources, mixed));
  lr_seed_Tn(*nn_only, lr_odd_scale, lr_seed + 31);
  lr_seed_Tn(*with_long, lr_odd_scale, lr_seed + 31);
  nn_only->update_CTM();
  // The same environment, copied rather than recomputed, so that the two
  // solvers differ in nothing but the observables they measure.
  lr_acc::C1(*with_long) = lr_acc::C1(*nn_only);
  lr_acc::C2(*with_long) = lr_acc::C2(*nn_only);
  lr_acc::C3(*with_long) = lr_acc::C3(*nn_only);
  lr_acc::C4(*with_long) = lr_acc::C4(*nn_only);
  lr_acc::eTt(*with_long) = lr_acc::eTt(*nn_only);
  lr_acc::eTr(*with_long) = lr_acc::eTr(*nn_only);
  lr_acc::eTb(*with_long) = lr_acc::eTb(*nn_only);
  lr_acc::eTl(*with_long) = lr_acc::eTl(*nn_only);
  const auto a = nn_only->measure_twosite();
  std::vector<std::map<Bond, double>> b;
  REQUIRE_NOTHROW(b = with_long->measure_twosite());
  REQUIRE(b.size() >= groups.size());
  for (std::size_t g = 0; g < groups.size(); ++g) {
    for (const auto& [bond, v] : a[g]) {
      INFO("group " << groups[g].kind << " "
                    << lr_bond_label(bond.source_site, bond.dx, bond.dy));
      REQUIRE(b[g].count(bond) == 1);
      lr_check_close("nearest neighbour alone vs with long range",
                     b[g].at(bond), v, 0.0, 1.0e-12);
    }
  }
}

// ============================================================================
// Contract item 7: the ops form
// ============================================================================

namespace {

template <class tensor>
void lr_run_ops_form(const std::string& label, double env_noise) {
  INFO(label);
  const int d = 2;
  const tenes::SquareLattice lattice = lr_lattice(2, 2, 0, d, lr_D);
  // One-site groups: 0 = n, 1 = c+, 2 = c.
  const char* names[3] = {"n", "cdag", "c"};
  tenes::Operators<tensor> onesite;
  for (int g = 0; g < 3; ++g) {
    for (int s = 0; s < lattice.N_UNIT; ++s) {
      onesite.emplace_back(
          names[g], g, s,
          lr_onesite_tensor<tensor>(lr_onesite_op(names[g], d)));
    }
  }
  // Pairs (ops form, explicit form of the same product, coefficient):
  // group 2k is ops = [i, j], group 2k + 1 the explicit A_s B_t table.
  struct pair_case {
    int i;
    int j;
    lr_cplx coeff;
  };
  const pair_case pairs[] = {{1, 2, 1.0}, {2, 1, 0.9}, {0, 0, 1.0}};
  // The contract's (1,0) and (2,1), and extra ones whose source is the
  // right / lower end of the window.
  const std::vector<lr_disp> disps = {
      {1, 0}, {2, 1}, {-1, 0}, {0, -1}, {-2, -1}};
  const std::vector<int> sources = {0, 1, 2, 3};
  tenes::Operators<tensor> twosite;
  for (int k = 0; k < 3; ++k) {
    const pair_case& pc = pairs[k];
    const tensor op =
        lr_twosite_tensor<tensor>({{1.0, names[pc.i], names[pc.j]}}, d, d);
    for (const int s : sources) {
      for (const lr_disp& dd : disps) {
        twosite.emplace_back("ops", 2 * k, s, dd.dx, dd.dy,
                             std::vector<int>{pc.i, pc.j},
                             lr_value<tensor>(pc.coeff));
        twosite.emplace_back("explicit", 2 * k + 1, s, dd.dx, dd.dy, op,
                             lr_value<tensor>(pc.coeff));
      }
    }
  }
  auto state = lr_make_state<tensor>(lattice, d, false, onesite, twosite);
  lr_seed_Tn(*state, lr_odd_scale, lr_seed + 41);
  state->update_CTM();
  if (env_noise > 0.0) {
    lr_perturb_env(*state, env_noise, lr_seed + 42);
  }
  // Premise: every compared value carries signal (direct reference, whose
  // non-cancelling scale also sets the tolerance below).
  std::map<std::tuple<int, int, int, int>, double> scales;
  for (int k = 0; k < 3; ++k) {
    const pair_case& pc = pairs[k];
    const tensor op =
        lr_twosite_tensor<tensor>({{1.0, names[pc.i], names[pc.j]}}, d, d);
    for (const int s : sources) {
      for (const lr_disp& dd : disps) {
        const auto ref =
            lr_direct(*state, s, dd.dx, dd.dy, op, lr_value<tensor>(pc.coeff));
        INFO("premise ops [" << pc.i << ", " << pc.j << "] "
                             << lr_bond_label(s, dd.dx, dd.dy) << ": "
                             << ref.value);
        REQUIRE(std::abs(ref.value) > 1.0e-6);
        scales[{k, s, dd.dx, dd.dy}] = ref.scale;
      }
    }
  }
  std::vector<std::map<Bond, typename tensor::value_type>> measured;
  REQUIRE_NOTHROW(measured = state->measure_twosite());
  REQUIRE(measured.size() >= 6);
  for (int k = 0; k < 3; ++k) {
    for (const int s : sources) {
      for (const lr_disp& dd : disps) {
        const Bond b{s, dd.dx, dd.dy};
        const std::string what = label + " ops [" + std::to_string(pairs[k].i) +
                                 ", " + std::to_string(pairs[k].j) + "] " +
                                 lr_bond_label(s, dd.dx, dd.dy);
        REQUIRE(measured[2 * k].count(b) == 1);
        REQUIRE(measured[2 * k + 1].count(b) == 1);
        const auto want = measured[2 * k + 1].at(b);
        lr_check_close(what + " [ops form vs explicit form]",
                       measured[2 * k].at(b), want,
                       scales.at({k, s, dd.dx, dd.dy}), 1.0e-12);
      }
    }
  }
}

}  // namespace

TEST_CASE(
    "longrange T2-7: the ops form ops = [i, j] equals the explicit product "
    "A_s B_t at (1,0) and (2,1)") {
  lr_run_ops_form<real_tensor>("T2-7 real converged", 0.0);
  lr_run_ops_form<real_tensor>("T2-7 real perturbed", 0.3);
  lr_run_ops_form<complex_tensor>("T2-7 complex perturbed", 0.3);
}

// ============================================================================
// Contract item 9: skew (Review Focus 2)
// ============================================================================
//
// The cell [2, 1] with skew 1 tiles the plane like [2, 2] with skew 0: with
// T(x, y + k LY) = T(x - k s, y) (src/SquareLattice.hpp), unfolded site
// (x, y) holds skew-cell site ((x - floor(y / LY)) mod 2, 0), i.e. unfolded
// sites 0, 1, 2, 3 hold skew sites 0, 1, 1, 0 (test/fermion/skew_unfold.cpp,
// computed here by arithmetic). The two cells run different CTM sweeps, so
// their environments agree only up to the finite-chi residual; the state and
// CTM settings are skew_unfold.cpp's short-ranged ones (odd virtual
// components scaled by 0.3, D = 2, chi = 8, epsilon 1e-12), where that
// residual is ~1e-11 for nearest neighbours.

namespace {

constexpr int lr_skew_chi = 8;
constexpr int lr_skew_ctm_iteration_max = 100;
constexpr double lr_skew_ctm_epsilon = 1.0e-12;
constexpr double lr_skew_odd_scale = 0.3;
constexpr unsigned lr_skew_seed = 1009;
//! Tolerance of item 9, anchored on the direct values of the two cells (the
//! [truth] case below prints them). Measured 2026-09-25 (Debug, g++-16,
//! macOS arm64, one thread): both CTMs converge in 8 sweeps, the largest
//! difference between the cells' direct values is 6.3e-12 and the smallest
//! |value| 7.9e-3. The tolerance is skew_unfold.cpp's sku_tol_ctm: 160 times
//! the residual, and 8e6 times below the smallest value, so a geometry or
//! sign error (a change of the order of the value) cannot hide in it.
constexpr double lr_skew_tol = 1.0e-9;

int lr_skew_preimage(int u) {
  const int x = u % 2;
  const int y = u / 2;
  return ((x - y) % 2 + 2) % 2;
}

struct lr_skew_pair {
  tenes::SquareLattice skew_lattice = lr_lattice(2, 1, 1, 2, lr_D);
  tenes::SquareLattice unfolded_lattice = lr_lattice(2, 2, 0, 2, lr_D);
  std::unique_ptr<lr_state<real_tensor>> skew;
  std::unique_ptr<lr_state<real_tensor>> unfolded;
  int count_skew = 0;
  int count_unfolded = 0;
};

const std::vector<lr_group> lr_skew_groups = {
    {"hop", -1.0}, {"n1mn", 1.0}, {"cdagc", 1.0}};
const std::vector<lr_disp> lr_skew_disps = {{2, 0}, {1, 1}};

int lr_converge_ctm(lr_state<real_tensor>& state) {
  const std::vector<real_tensor> reduced = rf::build_reduced_density_tensors(
      lr_acc::Tn(state), lr_acc::finfo(state));
  return tenes::itps::core::Calc_CTM_Environment_density(
      lr_acc::C1(state), lr_acc::C2(state), lr_acc::C3(state),
      lr_acc::C4(state), lr_acc::eTt(state), lr_acc::eTr(state),
      lr_acc::eTb(state), lr_acc::eTl(state), reduced,
      lr_acc::peps_parameters(state), lr_acc::lattice(state), true, true);
}

std::unique_ptr<lr_skew_pair> lr_make_skew_pair() {
  auto p = std::make_unique<lr_skew_pair>();
  for (int s = 0; s < 2; ++s) {
    REQUIRE(lr_skew_preimage(s) == s);
  }
  // Premise on the geometry this relies on (tested in skew_unfold.cpp).
  for (int u = 0; u < 4; ++u) {
    for (int leg = 0; leg < 4; ++leg) {
      REQUIRE(lr_skew_preimage(p->unfolded_lattice.neighbor(u, leg)) ==
              p->skew_lattice.neighbor(lr_skew_preimage(u), leg));
    }
  }
  const auto make = [](const tenes::SquareLattice& lat) {
    return lr_make_state<real_tensor>(
        lat, 2, false, tenes::Operators<real_tensor>{},
        lr_twosite_ops<real_tensor>(
            lr_skew_groups, 2,
            [&] {
              std::vector<int> all;
              for (int s = 0; s < lat.N_UNIT; ++s) {
                all.push_back(s);
              }
              return all;
            }(),
            lr_skew_disps),
        lr_skew_chi, lr_skew_ctm_iteration_max, lr_skew_ctm_epsilon);
  };
  p->skew = make(p->skew_lattice);
  p->unfolded = make(p->unfolded_lattice);
  lr_seed_Tn(*p->skew, lr_skew_odd_scale, lr_skew_seed);
  auto& Ts = lr_acc::Tn(*p->skew);
  auto& Tu = lr_acc::Tn(*p->unfolded);
  auto& fs = lr_acc::finfo(*p->skew);
  auto& fu = lr_acc::finfo(*p->unfolded);
  for (int u = 0; u < 4; ++u) {
    Tu[u] = Ts[lr_skew_preimage(u)];
    fu.virt[u] = fs.virt[lr_skew_preimage(u)];
  }
  REQUIRE_NOTHROW(rf::validate_neighbor_consistency(fu, p->unfolded_lattice));
  REQUIRE_NOTHROW(rf::validate_neighbor_consistency(fs, p->skew_lattice));
  p->count_skew = lr_converge_ctm(*p->skew);
  p->count_unfolded = lr_converge_ctm(*p->unfolded);
  INFO("premise: CTM sweeps " << p->count_skew << " (skew) and "
                              << p->count_unfolded << " (unfolded)");
  CHECK(p->count_skew < lr_skew_ctm_iteration_max);
  CHECK(p->count_unfolded < lr_skew_ctm_iteration_max);
  return p;
}

}  // namespace

TEST_CASE(
    "longrange [truth] T2-9 anchor: the direct values of the skew cell and "
    "of its unfolded cell agree within the item-9 tolerance and carry "
    "signal") {
  auto p = lr_make_skew_pair();
  double max_diff = 0.0;
  double min_abs = std::numeric_limits<double>::infinity();
  for (std::size_t g = 0; g < lr_skew_groups.size(); ++g) {
    const real_tensor op =
        lr_kind_tensor<real_tensor>(lr_skew_groups[g].kind, 2);
    const double coeff = lr_skew_groups[g].coeff.real();
    for (int u = 0; u < 4; ++u) {
      for (const lr_disp& dd : lr_skew_disps) {
        const auto vu = lr_direct(*p->unfolded, u, dd.dx, dd.dy, op, coeff);
        const auto vs =
            lr_direct(*p->skew, lr_skew_preimage(u), dd.dx, dd.dy, op, coeff);
        max_diff = std::max(max_diff, std::abs(vu.value - vs.value));
        min_abs = std::min(min_abs, std::abs(vu.value));
      }
    }
  }
  std::cout << std::setprecision(3) << "longrange skew anchor: CTM sweeps "
            << p->count_skew << " / " << p->count_unfolded
            << ", max |direct skew - direct unfolded| " << max_diff
            << ", min |value| " << min_abs << " (tolerance " << lr_skew_tol
            << ")" << std::endl;
  CHECK(max_diff <= lr_skew_tol / 30.0);
  CHECK(min_abs > 1.0e3 * lr_skew_tol);
}

TEST_CASE(
    "longrange T2-9: skew = 1 [2,1] cell equals its unfolded skew-0 [2,2] "
    "cell at (2,0) and (1,1) (Review Focus 2)") {
  auto p = lr_make_skew_pair();
  std::vector<std::map<Bond, double>> ms, mu;
  REQUIRE_NOTHROW(ms = p->skew->measure_twosite());
  REQUIRE_NOTHROW(mu = p->unfolded->measure_twosite());
  REQUIRE(ms.size() >= lr_skew_groups.size());
  REQUIRE(mu.size() >= lr_skew_groups.size());
  double max_diff = 0.0;
  for (std::size_t g = 0; g < lr_skew_groups.size(); ++g) {
    for (int u = 0; u < 4; ++u) {
      for (const lr_disp& dd : lr_skew_disps) {
        const Bond bu{u, dd.dx, dd.dy};
        const Bond bs{lr_skew_preimage(u), dd.dx, dd.dy};
        INFO("group " << lr_skew_groups[g].kind << " unfolded "
                      << lr_bond_label(u, dd.dx, dd.dy) << " -> skew site "
                      << bs.source_site);
        REQUIRE(mu[g].count(bu) == 1);
        REQUIRE(ms[g].count(bs) == 1);
        const double diff = std::abs(mu[g].at(bu) - ms[g].at(bs));
        max_diff = std::max(max_diff, diff);
        CHECK(diff <= lr_skew_tol);
      }
    }
  }
  std::cout << std::setprecision(3) << "longrange skew: max |skew - unfolded| "
            << max_diff << std::endl;
}

// ============================================================================
// Contract item 1: input acceptance and rejection
// ============================================================================

namespace {

using lr_gtensor = real_tensor;

struct lr_guard_input {
  tenes::itps::PEPS_Parameters params;
  tenes::SquareLattice lattice = lr_lattice(2, 2, 0, 2, 2);
  tenes::EvolutionOperators<lr_gtensor> simple_updates;
  tenes::Operators<lr_gtensor> onesite;
  tenes::Operators<lr_gtensor> twosite;
  tenes::Operators<lr_gtensor> multisite;
  tenes::itps::CorrelationParameter corparam;

  //! One-site groups on every site: 0 = n (even) and, with odd_onesite,
  //! 1 = c+ and 2 = c (odd). Cases about anything but odd one-site operators
  //! or the ops form leave the odd ones out, so that against the stub (which
  //! rejects odd one-site operators first) they fail or pass for their own
  //! reason only.
  explicit lr_guard_input(int d = 2, bool meanfield = false,
                          bool odd_onesite = true) {
    lattice = lr_lattice(2, 2, 0, d, 2);
    params = lr_params<lr_gtensor>(lattice.N_UNIT, d, meanfield, lr_chi,
                                   lr_ctm_iteration_max, lr_ctm_epsilon);
    const std::vector<std::string> names =
        d == 2 ? std::vector<std::string>{"n", "cdag", "c"}
               : std::vector<std::string>{"n", "cdag_up", "c_up"};
    for (int g = 0; g < (odd_onesite ? 3 : 1); ++g) {
      for (int s = 0; s < lattice.N_UNIT; ++s) {
        onesite.emplace_back(
            names[g], g, s,
            lr_onesite_tensor<lr_gtensor>(lr_onesite_op(names[g], d)));
      }
    }
  }

  void add_twosite(int dx, int dy, const std::string& kind = "hop") {
    const int d = lattice.physical_dims[0];
    twosite.emplace_back(kind, 0, 0, dx, dy,
                         lr_kind_tensor<lr_gtensor>(kind, d));
  }

  void add_ops(int dx, int dy, int i, int j) {
    twosite.emplace_back("ops", 0, 0, dx, dy, std::vector<int>{i, j});
  }

  void validate() const {
    tenes::itps::validate_fermion_constraints(
        params, lattice, simple_updates,
        tenes::EvolutionOperators<lr_gtensor>{}, onesite, twosite, multisite,
        corparam);
  }
};

void lr_check_accepts(const lr_guard_input& in, const std::string& what) {
  INFO(what);
  try {
    in.validate();
  } catch (const std::exception& e) {
    FAIL_CHECK(what << " was rejected: " << std::string(e.what()));
  }
}

void lr_check_rejects(const lr_guard_input& in, const std::string& what,
                      const std::string& needle = "") {
  INFO(what);
  try {
    in.validate();
    FAIL_CHECK(what << " was accepted");
  } catch (const tenes::input_error& e) {
    if (!needle.empty()) {
      INFO("message: " << std::string(e.what()));
      CHECK(std::string(e.what()).find(needle) != std::string::npos);
    }
  } catch (const std::exception& e) {
    FAIL_CHECK(what << " threw something other than tenes::input_error: "
                    << std::string(e.what()));
  }
}

std::string lr_disp_name(int dx, int dy) {
  return "(" + std::to_string(dx) + ", " + std::to_string(dy) + ")";
}

}  // namespace

TEST_CASE(
    "longrange T2-1a: two-site observables inside the 4x4 window are "
    "accepted") {
  for (const lr_disp dd :
       {lr_disp{2, 0}, lr_disp{0, 2}, lr_disp{1, 1}, lr_disp{-1, 2},
        lr_disp{3, 0}, lr_disp{2, -1}, lr_disp{3, 3}, lr_disp{-3, -3},
        lr_disp{0, -3}, lr_disp{-3, 1}}) {
    lr_guard_input in(2, false, false);
    in.add_twosite(dd.dx, dd.dy);
    lr_check_accepts(in, "hopping at " + lr_disp_name(dd.dx, dd.dy));
  }
  // d = 4 as well.
  lr_guard_input in4(4, false, false);
  in4.add_twosite(2, 1, "hopnn");
  lr_check_accepts(in4, "d = 4 hopping + nn at (2, 1)");
}

TEST_CASE(
    "longrange T2-1b: two-site observables beyond the 4x4 window are "
    "rejected, naming 4x4") {
  for (const lr_disp dd : {lr_disp{4, 0}, lr_disp{0, 4}, lr_disp{-4, 1},
                           lr_disp{1, -5}, lr_disp{4, 4}}) {
    lr_guard_input in(2, false, false);
    in.add_twosite(dd.dx, dd.dy);
    lr_check_rejects(in, "hopping at " + lr_disp_name(dd.dx, dd.dy), "4x4");
  }
}

TEST_CASE(
    "longrange [kept] T2-1c: a same-site (0, 0) two-site observable is "
    "rejected") {
  lr_guard_input in(2, false, false);
  in.add_twosite(0, 0);
  lr_check_rejects(in, "hopping at (0, 0)");
  lr_guard_input ops(2, false, false);
  ops.add_ops(0, 0, 0, 0);
  lr_check_rejects(ops, "ops = [0, 0] at (0, 0)");
}

TEST_CASE(
    "longrange T2-1d: parity-odd one-site operators are accepted, mixed "
    "ones are rejected naming mixed parity") {
  for (const int d : {2, 4}) {
    INFO("d = " << d);
    lr_guard_input odd(d);  // groups 1 and 2 are odd
    lr_check_accepts(odd, "odd one-site operators c+, c");
    lr_guard_input odd_mf(d, true);
    lr_check_accepts(odd_mf, "odd one-site operators c+, c (mean field)");

    lr_guard_input mixed(d, false, false);
    lr_onesite o = lr_onesite_op(d == 2 ? "c" : "c_dn", d);
    const lr_onesite n = lr_onesite_op("n", d);
    for (int k = 0; k < d * d; ++k) {
      o.m[k] += 0.5 * n.m[k];
    }
    mixed.onesite.emplace_back("mixed", 1, 1, lr_onesite_tensor<lr_gtensor>(o));
    lr_check_rejects(mixed, "a mixed-parity one-site operator", "mixed parity");
  }
}

TEST_CASE(
    "longrange T2-1e: the ops form is accepted when both one-site operators "
    "have the same parity") {
  for (const lr_disp dd : {lr_disp{1, 0}, lr_disp{0, -1}, lr_disp{2, 0},
                           lr_disp{2, 1}, lr_disp{-3, 3}}) {
    for (const auto& ij :
         std::vector<std::pair<int, int>>{{1, 2}, {2, 1}, {1, 1}, {0, 0}}) {
      // [0, 0] needs no odd one-site operator: against the stub it is then
      // refused for being an ops form, not for the odd operators.
      lr_guard_input in(2, false, ij.first != 0 || ij.second != 0);
      in.add_ops(dd.dx, dd.dy, ij.first, ij.second);
      lr_check_accepts(in, "ops = [" + std::to_string(ij.first) + ", " +
                               std::to_string(ij.second) + "] at " +
                               lr_disp_name(dd.dx, dd.dy));
    }
  }
  // An ops form beyond the window is rejected like an explicit one.
  lr_guard_input far;
  far.add_ops(0, 4, 1, 2);
  lr_check_rejects(far, "ops = [1, 2] at (0, 4)", "4x4");
}

// Green against the stub by construction (the stub rejects every ops form);
// it pins that the new acceptance does not extend to operators of different
// parity.
TEST_CASE(
    "longrange [kept] T2-1f: an ops form of two one-site operators of "
    "different parity is rejected") {
  for (const lr_disp dd : {lr_disp{1, 0}, lr_disp{2, 1}}) {
    for (const auto& ij : std::vector<std::pair<int, int>>{{0, 1}, {2, 0}}) {
      lr_guard_input in;
      in.add_ops(dd.dx, dd.dy, ij.first, ij.second);
      lr_check_rejects(in, "ops = [" + std::to_string(ij.first) + ", " +
                               std::to_string(ij.second) + "] at " +
                               lr_disp_name(dd.dx, dd.dy));
    }
  }
}

TEST_CASE(
    "longrange T2-1g: with meanfield_env the nearest-neighbour ops form is "
    "accepted") {
  for (const lr_disp dd : lr_nn_disps) {
    lr_guard_input in(2, true);
    in.add_ops(dd.dx, dd.dy, 1, 2);
    lr_check_accepts(
        in, "mean field, ops = [1, 2] at " + lr_disp_name(dd.dx, dd.dy));
  }
}

// Task T5 of the same plan (contract item 1) lifts the mean-field
// restriction: this case used to pin that the MF environment rejects
// long-range observables until T5 ("[kept] T2-1h ... stay rejected") and now
// requires that they be accepted. Their values are checked in
// test/fermion/longrange_mf.cpp. Red until T5 is implemented.
TEST_CASE(
    "longrange T2-1h: with meanfield_env long-range two-site observables and "
    "long-range ops forms are accepted (task T5)") {
  for (const lr_disp dd : {lr_disp{2, 0}, lr_disp{1, 1}, lr_disp{0, -3}}) {
    lr_guard_input in(2, true, false);
    in.add_twosite(dd.dx, dd.dy);
    lr_check_accepts(in,
                     "mean field, hopping at " + lr_disp_name(dd.dx, dd.dy));
    lr_guard_input ops(2, true);
    ops.add_ops(dd.dx, dd.dy, 1, 2);
    lr_check_accepts(
        ops, "mean field, ops = [1, 2] at " + lr_disp_name(dd.dx, dd.dy));
  }
}

TEST_CASE(
    "longrange [kept] T2-1i: odd one-site gates, parity-odd two-site "
    "observables and multisite observables stay rejected") {
  {
    lr_guard_input in(2, false, false);
    in.simple_updates.push_back(
        tenes::make_onesite_EvolutionOperator<lr_gtensor>(
            0, 0, lr_onesite_tensor<lr_gtensor>(lr_onesite_op("c", 2))));
    lr_check_rejects(in, "an odd one-site gate");
  }
  {
    lr_guard_input in(2, false, false);
    // c+_s 1_t is odd as a two-site operator, at nearest and longer range.
    in.twosite.emplace_back(
        "odd", 0, 0, 1, 0,
        lr_twosite_tensor<lr_gtensor>({{1.0, "cdag", "I"}}, 2, 2));
    lr_check_rejects(in, "a parity-odd two-site observable at (1, 0)");
    lr_guard_input far(2, false, false);
    far.twosite.emplace_back(
        "odd", 0, 0, 2, 0,
        lr_twosite_tensor<lr_gtensor>({{1.0, "cdag", "I"}}, 2, 2));
    lr_check_rejects(far, "a parity-odd two-site observable at (2, 0)");
  }
  {
    lr_guard_input in(2, false, false);
    lr_gtensor op3(mptensor::Shape(2, 2, 2, 2, 2, 2));
    op3.set_value(mptensor::Index(1, 1, 1, 1, 1, 1), 1.0);
    in.multisite.emplace_back("nnn", 0, 0, std::vector<int>{1, 2},
                              std::vector<int>{0, 0}, op3);
    lr_check_rejects(in, "a multisite observable");
  }
}

// ---- validate_fermion_ctm_measurement: the same conditions ----------------

namespace {

template <class Fn>
void lr_check_measure_rejects(Fn&& fn, const std::string& what) {
  INFO(what);
  try {
    fn();
    FAIL_CHECK(what << " was measured");
  } catch (const tenes::input_error&) {
  } catch (const std::exception& e) {
    FAIL_CHECK(what << " threw something other than tenes::input_error: "
                    << std::string(e.what()));
  }
}

}  // namespace

// Green against the stub (it rejects every non-nearest pair); pins that the
// bosonic fallback (a warning and a skipped pair beyond 4x4, or a
// raw-Tn contraction at (0, 0)) never runs in fermion mode.
TEST_CASE(
    "longrange [kept] T2-1j: measure_twosite rejects (0, 0) and pairs beyond "
    "the 4x4 window in fermion CTM mode") {
  for (const lr_disp dd : {lr_disp{0, 0}, lr_disp{4, 0}, lr_disp{0, -4}}) {
    const tenes::SquareLattice lattice = lr_lattice(2, 2, 0, 2, lr_D);
    auto state = lr_make_state<real_tensor>(
        lattice, 2, false, tenes::Operators<real_tensor>{},
        lr_twosite_ops<real_tensor>({{"hop", 1.0}}, 2, {0}, {dd}));
    lr_seed_Tn(*state, lr_odd_scale, lr_seed);
    state->update_CTM();
    lr_check_measure_rejects(
        [&] { state->measure_twosite(); },
        "measure_twosite at " + lr_disp_name(dd.dx, dd.dy));
  }
}
