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
//! Task T5 of docs/superpowers/plans/2026-09-25-fermion-longrange-measure.md:
//! long-range two-site observables, the ops form and the correlation
//! function in fermion mode with the mean-field environment
//! (meanfield_env = true; design
//! docs/superpowers/specs/2026-09-25-fermion-longrange-measure-design.md,
//! section 6). Behaviour contract: work/fermion-longrange/t5/contract.md,
//! items 1 to 5 (item 6 is test/fermion/free_fermion_longrange_mf.py.in).
//!
//! Conventions fixed for every case in this file:
//!   - One-site operators use TeNeS' layout op[in, out] = <out|A|in>; two-site
//!     operators op[in_s, in_t, out_s, out_t], source first.
//!   - d = 2 basis {|0>, |1>}, ledger [e, o]. d = 4 basis {|0>, |up>, |dn>,
//!     |updn> = c+_up c+_dn |0>}, ledger [e, o, o, e].
//!   - The product A_s B_t is written out by the test itself (lm_product_el),
//!     op4[i_s, i_t, o_s, o_t] = (-1)^{p_B p(i_s)} A[i_s, o_s] B[i_t, o_t]
//!     (design section 4.2), never through product_twosite_op(), except in
//!     contract item 5, whose window side is by definition "the same pair
//!     measured as a window" and uses the T2 function the solver's ops form
//!     uses.
//!   - A window follows twosite_obs.cpp: ncol = |dx| + 1, nrow = |dy| + 1,
//!     the source at column 0 (dx >= 0) or ncol - 1, at row nrow - 1
//!     (dy >= 0) or 0, the target in the opposite corner; window cell
//!     (row, col) holds lattice.other(source, col - source_col,
//!     source_row - row). Row 0 is the top row.
//!   - The mean-field weights (design section 6, the bosonic convention of
//!     twosite_obs.cpp and measure_correlation_mf()): a window cell carries
//!     lambda[site][leg] on each of its legs that lies on the window's
//!     perimeter (row 0: t, last row: b, column 0: l, last column: r), and
//!     nothing on the window's internal bonds. A correlation chain from the
//!     left (lower) site to the right (upper) site is the window of the same
//!     pair, so it carries the same weights.
//!   - Every solver state has distinct horizontal and vertical bond
//!     dimensions where the cost allows it (Dh != Dv), so that a delta edge
//!     of the wrong leg's dimension cannot go unnoticed, and a lambda that
//!     differs from bond to bond, so that a weight on the wrong leg or the
//!     wrong site changes the value.
//!
//! Truth sources (never the T5 code: not the delta environment, not the
//! mean-field branches of twosite_obs.cpp / correlation_function.cpp):
//!   - contract_pair_MF() of the lambda-dressed pair state, the existing
//!     single-layer graded mean-field path of the nearest-neighbour
//!     measurement (verified against the open-leg Fock oracle in
//!     test/fermion/mf_measure.cpp); this is contract item 3's reference;
//!   - its window generalization, the "graded window truth" (lm_truth_*):
//!     every cell of the window, lambda-dressed as above, contracted into
//!     one single-layer graded state psi whose perimeter legs stay open,
//!     and <O> = trace(conj(psi), O psi) / trace(conj(psi), psi) over every
//!     open leg, which is exactly what contract_pair_MF() does for a pair.
//!     A [truth] case pins it to contract_pair_MF() at every nearest-
//!     neighbour direction and, with every parity even, to the bosonic
//!     mean-field solver (Contract_iTPS_MF and the *Correlation_iTPS_MF
//!     kernels), which fixes its lambda placement independently of T5;
//!   - the bosonic mean-field solver itself for contract item 4 (every
//!     parity even);
//!   - for contract item 5, the window measurement of the same pair on a
//!     second solver without [correlation].
//!
//! Why contract item 3 (gate 1) is not measured through measure_twosite() at
//! the nearest-neighbour displacements: that path stays on contract_pair_MF()
//! by design (plan, Global Constraints: the nearest-neighbour measurement
//! path is not changed), so the comparison would be the existing path against
//! itself and would pass against the stub. Gate 1 is checked instead where
//! the T5 construction is actually used at nearest-neighbour distance:
//!   - T5-3a: the construction itself (lambda-dressed window, relay window of
//!     T1, the delta environment of make_delta_corner() / make_delta_edge(),
//!     core::Contract_density_CTM()), assembled here, at the four
//!     nearest-neighbour directions, against contract_pair_MF();
//!   - T5-3b: the correlation function at r = 1 (the solver's T5 path)
//!     against contract_pair_MF(), and at r = 2, 3 against the graded window
//!     truth;
//!   - T5-3c: long-range windows and ops forms of measure_twosite() against
//!     the graded window truth.
//!
//! Cases whose name carries [truth] check the reference side only, and cases
//! whose name carries [kept] pin a rejection that the stub already performs;
//! both are expected to pass against the stub. Every other case must fail
//! until T5 is implemented.

#include "../test_fermion_common.hpp"

#include <algorithm>
#include <complex>
#include <cstdint>
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
using lm_acc = tenes::itps::iTPSTestAccessor;
using tenes::complex_tensor;
using tenes::real_tensor;
using tenes::itps::Bond;
using tenes::itps::Correlation;
using tenes::itps::CorrelationParameter;
using lm_cplx = std::complex<double>;

template <class tensor>
using lm_ft = rf::ftensor<tensor>;

template <class tensor>
using lm_state = tenes::itps::iTPS<tensor>;

// ---- small helpers ---------------------------------------------------------

template <class tensor>
const char* lm_type_name() {
  return std::is_same<tensor, complex_tensor>::value ? "complex" : "real";
}

template <class tensor>
double lm_max_abs_entry(const tensor& a) {
  double m = 0.0;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    m = std::max(m, std::abs(a[n]));
  }
  return m;
}

// |got - want| <= rtol * max(|want|, scale). scale is a magnitude that does
// not cancel (the largest operator element times |coeff|: every value
// compared here is a convex combination of expectation values of that
// operator), so that a small result is not judged by a tolerance
// proportional to itself.
bool lm_check_close(const std::string& label, lm_cplx got, lm_cplx want,
                    double scale, double rtol) {
  const double tol = rtol * std::max(std::abs(want), scale);
  const double diff = std::abs(got - want);
  INFO(label << ": got=" << got << " want=" << want << " |diff|=" << diff
             << " tol=" << tol);
  CHECK(diff <= tol);
  return diff <= tol;
}

inline rf::parity_vector lm_phys(int d) {
  if (d == 2) {
    return {false, true};
  }
  if (d == 4) {
    return {false, true, true, false};
  }
  throw std::runtime_error("lm_phys: unsupported physical dimension");
}

// ---- deterministic entries -----------------------------------------------
//
// Every entry is a hash of (seed, salt, global index), not a draw from a
// stream over local elements, so a state means the same at any rank count.

inline std::uint64_t lm_mix(std::uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

//! Uniform in [-1, 1).
inline double lm_unit(std::uint64_t seed, std::uint64_t salt,
                      const mptensor::Index& idx) {
  std::uint64_t h = lm_mix(seed * 1000003ULL + salt);
  for (std::size_t k = 0; k < idx.size(); ++k) {
    h = lm_mix(h ^ (static_cast<std::uint64_t>(idx[k]) + 0x100ULL * (k + 1)));
  }
  return 2.0 * (static_cast<double>(h >> 11) * 0x1.0p-53) - 1.0;
}

inline void lm_set_hashed(real_tensor& t, const mptensor::Index& idx,
                          std::uint64_t seed, double scale) {
  t.set_value(idx, scale * lm_unit(seed, 0, idx));
}

inline void lm_set_hashed(complex_tensor& t, const mptensor::Index& idx,
                          std::uint64_t seed, double scale) {
  t.set_value(idx,
              scale * lm_cplx(lm_unit(seed, 0, idx), lm_unit(seed, 1, idx)));
}

//! Odd virtual components are scaled by this per odd index, so that the odd
//! sector is present but does not dominate.
constexpr double lm_odd_scale = 0.8;

//! A rank-5 site (l, t, r, b, s) with ledgers lp. Parity-even entries only,
//! unless dense (every entry, for the all-even comparison of item 4).
template <class tensor>
tensor lm_site_tensor(const rf::leg_parities& lp, std::uint64_t seed,
                      bool dense) {
  mptensor::Shape sh;
  for (const auto& leg : lp) {
    sh.push(leg.size());
  }
  tensor t(sh);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    if (!dense && rf::count_odd(lp, idx) % 2 != 0) {
      continue;
    }
    double scale = 1.0;
    for (int leg = 0; leg < 4; ++leg) {
      if (lp[leg][idx[leg]]) {
        scale *= lm_odd_scale;
      }
    }
    lm_set_hashed(t, idx, seed, scale);
  }
  return t;
}

//! A mean-field weight vector of dimension D: 1 first, then distinct values
//! in [0.25, 0.85).
inline std::vector<double> lm_lambda_vector(int D, std::uint64_t seed, int site,
                                            int leg) {
  std::vector<double> v(D, 1.0);
  for (int k = 1; k < D; ++k) {
    v[k] = 0.55 + 0.3 * lm_unit(seed, 7, mptensor::Index(site, leg, k));
  }
  return v;
}

// ---- one-site operators as element tables -----------------------------------

// m[in * d + out] = <out|A|in>.
struct lm_onesite {
  int d = 2;
  std::vector<lm_cplx> m;
  bool odd = false;
};

inline lm_onesite lm_onesite_op(const std::string& name, int d) {
  lm_onesite o;
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
      throw std::runtime_error("lm_onesite_op: unknown d=2 operator " + name);
    }
    return o;
  }
  if (d != 4) {
    throw std::runtime_error("lm_onesite_op: unsupported dimension");
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
  } else {
    throw std::runtime_error("lm_onesite_op: unknown d=4 operator " + name);
  }
  return o;
}

inline void lm_set_scalar(real_tensor& t, const mptensor::Index& idx,
                          lm_cplx v) {
  REQUIRE(v.imag() == 0.0);
  t.set_value(idx, v.real());
}

inline void lm_set_scalar(complex_tensor& t, const mptensor::Index& idx,
                          lm_cplx v) {
  t.set_value(idx, v);
}

template <class tensor>
tensor lm_onesite_tensor(const lm_onesite& o) {
  tensor t(mptensor::Shape(o.d, o.d));
  for (int in = 0; in < o.d; ++in) {
    for (int out = 0; out < o.d; ++out) {
      const lm_cplx v = o.m[in * o.d + out];
      if (v != 0.0) {
        lm_set_scalar(t, mptensor::Index(in, out), v);
      }
    }
  }
  return t;
}

template <class tensor>
tensor lm_identity(int d) {
  return lm_onesite_tensor<tensor>(lm_onesite_op("I", d));
}

// ---- two-site operators written out by the test -----------------------------

// coef * A_s B_t.
struct lm_term {
  lm_cplx coef;
  std::string a;
  std::string b;
};

// op4[i_s, i_t, o_s, o_t] = sum coef (-1)^{p_B p(i_s)} A[i_s, o_s]
// B[i_t, o_t], flattened ((i_s * d + i_t) * d + o_s) * d + o_t.
inline std::vector<lm_cplx> lm_product_el(const std::vector<lm_term>& terms,
                                          int d) {
  const rf::parity_vector ps = lm_phys(d);
  std::vector<lm_cplx> el(d * d * d * d, 0.0);
  for (const lm_term& term : terms) {
    const lm_onesite A = lm_onesite_op(term.a, d);
    const lm_onesite B = lm_onesite_op(term.b, d);
    for (int is = 0; is < d; ++is) {
      for (int it = 0; it < d; ++it) {
        for (int os = 0; os < d; ++os) {
          for (int ot = 0; ot < d; ++ot) {
            const double sign = (B.odd && ps[is]) ? -1.0 : 1.0;
            el[((is * d + it) * d + os) * d + ot] +=
                term.coef * sign * A.m[is * d + os] * B.m[it * d + ot];
          }
        }
      }
    }
  }
  return el;
}

template <class tensor>
tensor lm_twosite_tensor(const std::vector<lm_term>& terms, int d) {
  const auto el = lm_product_el(terms, d);
  tensor op(mptensor::Shape(d, d, d, d));
  for (int is = 0; is < d; ++is) {
    for (int it = 0; it < d; ++it) {
      for (int os = 0; os < d; ++os) {
        for (int ot = 0; ot < d; ++ot) {
          const lm_cplx v = el[((is * d + it) * d + os) * d + ot];
          if (v != 0.0) {
            lm_set_scalar(op, mptensor::Index(is, it, os, ot), v);
          }
        }
      }
    }
  }
  return op;
}

// Complex coefficients of the non-Hermitian operators: a conj() on the wrong
// layer or a source / target swap changes their value.
const lm_cplx lm_alpha(0.8, 0.6);
const lm_cplx lm_beta(-0.3, 0.5);
const lm_cplx lm_gamma(0.4, -0.7);
const lm_cplx lm_delta(-0.25, 0.35);

// Operator kinds:
//   cdagc  c+_s c_t (the first flavour), odd x odd only
//   hop    sum_sigma (c+_s c_t + c+_t c_s), c+_t c_s = -c_s c+_t
//   asym   c+_s c_t - 0.6 c_s c+_t + 0.8 n_s 1_t + 0.35 n_s n_t (d = 2):
//          even and odd channels, and not symmetric under s <-> t
//   asym4  c+_up,s c_up,t - 0.7 c_dn,s c+_dn,t + 0.4 n_s 1_t (d = 4)
//   cplx   sum_sigma (alpha c+_s c_t - beta c_s c+_t) + gamma n_s n_t
//          + delta n_s 1_t
inline std::vector<lm_term> lm_kind_terms(const std::string& kind, int d) {
  std::vector<std::pair<std::string, std::string>> flavours;
  if (d == 2) {
    flavours.push_back({"cdag", "c"});
  } else {
    flavours.push_back({"cdag_up", "c_up"});
    flavours.push_back({"cdag_dn", "c_dn"});
  }
  std::vector<lm_term> t;
  if (kind == "cdagc") {
    t.push_back({1.0, flavours[0].first, flavours[0].second});
  } else if (kind == "hop") {
    for (const auto& f : flavours) {
      t.push_back({1.0, f.first, f.second});
      t.push_back({-1.0, f.second, f.first});
    }
  } else if (kind == "asym") {
    t.push_back({1.0, "cdag", "c"});
    t.push_back({-0.6, "c", "cdag"});
    t.push_back({0.8, "n", "I"});
    t.push_back({0.35, "n", "n"});
  } else if (kind == "asym4") {
    t.push_back({1.0, "cdag_up", "c_up"});
    t.push_back({-0.7, "c_dn", "cdag_dn"});
    t.push_back({0.4, "n", "I"});
  } else if (kind == "cplx") {
    for (const auto& f : flavours) {
      t.push_back({lm_alpha, f.first, f.second});
      t.push_back({-lm_beta, f.second, f.first});
    }
    t.push_back({lm_gamma, "n", "n"});
    t.push_back({lm_delta, "n", "I"});
  } else {
    throw std::runtime_error("lm_kind_terms: unknown kind " + kind);
  }
  return t;
}

template <class tensor>
tensor lm_kind_tensor(const std::string& kind, int d) {
  return lm_twosite_tensor<tensor>(lm_kind_terms(kind, d), d);
}

//! The odd x odd part of a kind (its c / c+ terms only).
inline std::vector<lm_term> lm_odd_terms(const std::string& kind, int d) {
  std::vector<lm_term> out;
  for (const lm_term& t : lm_kind_terms(kind, d)) {
    if (lm_onesite_op(t.a, d).odd) {
      out.push_back(t);
    }
  }
  return out;
}

template <class tensor>
typename tensor::value_type lm_value(lm_cplx v) {
  if constexpr (std::is_same<tensor, real_tensor>::value) {
    return v.real();
  } else {
    return v;
  }
}

// ---- the solver fixture ----------------------------------------------------

//! A unit cell: lx x ly sites, physical dimension d, bond dimension Dh on
//! the horizontal bonds (legs l, r) and Dv on the vertical ones (t, b).
struct lm_geom {
  int lx;
  int ly;
  int d;
  int Dh;
  int Dv;
};

inline std::string lm_geom_label(const lm_geom& g) {
  return std::to_string(g.lx) + "x" + std::to_string(g.ly) +
         " cell d=" + std::to_string(g.d) + " Dh=" + std::to_string(g.Dh) +
         " Dv=" + std::to_string(g.Dv);
}

tenes::SquareLattice lm_lattice(const lm_geom& g) {
  tenes::SquareLattice lattice(g.lx, g.ly, 0);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    lattice.physical_dims[site] = g.d;
    lattice.virtual_dims[site] = {g.Dh, g.Dv, g.Dh, g.Dv};
    lattice.initial_dirs[site] = {0.0};
    lattice.noises[site] = 1.0;
  }
  return lattice;
}

template <class tensor>
tenes::itps::PEPS_Parameters lm_params(int n_unit, bool fermion,
                                       const rf::parity_vector& phys,
                                       bool meanfield) {
  tenes::itps::PEPS_Parameters p;
  p.fermion = fermion;
  p.is_real = std::is_same<tensor, real_tensor>::value;
  p.phys_parity.assign(n_unit, phys);
  p.print_level = tenes::PrintLevel::none;
  p.outdir = "output_test_fermion_longrange_mf";
  p.CHI = 4;
  p.Max_CTM_Iteration = 30;
  p.CTM_Convergence_Epsilon = 1.0e-10;
  p.Use_RSVD = false;
  p.MeanField_Env = meanfield;
  return p;
}

//! A solver of `g` with the mean-field environment.
template <class tensor>
std::unique_ptr<lm_state<tensor>> lm_make(
    const lm_geom& g, bool fermion, const rf::parity_vector& phys,
    const tenes::Operators<tensor>& onesite,
    const tenes::Operators<tensor>& twosite, const CorrelationParameter& cp) {
  const tenes::SquareLattice lattice = lm_lattice(g);
  return std::make_unique<lm_state<tensor>>(
      MPI_COMM_WORLD, lm_params<tensor>(lattice.N_UNIT, fermion, phys, true),
      lattice, tenes::EvolutionOperators<tensor>{},
      tenes::EvolutionOperators<tensor>{}, onesite, twosite,
      tenes::Operators<tensor>{}, cp, tenes::itps::TransferMatrix_Parameters{});
}

//! Seed Tn (site s from seed + 97 s) and lambda. A fermion solver draws
//! parity-even entries under its ledgers, unless `all_even`: then every
//! ledger (virtual and, through the parameters, physical) is even and every
//! entry is drawn, which a bosonic solver of the same seed reproduces.
template <class tensor>
void lm_seed(lm_state<tensor>& state, std::uint64_t seed, bool all_even) {
  const tenes::SquareLattice& lat = lm_acc::lattice(state);
  auto& fi = lm_acc::finfo(state);
  auto& Tn = lm_acc::Tn(state);
  REQUIRE(Tn.size() == static_cast<std::size_t>(lat.N_UNIT));
  for (int s = 0; s < lat.N_UNIT; ++s) {
    rf::leg_parities lp;
    for (int leg = 0; leg < 4; ++leg) {
      lp.push_back(rf::parity_vector(lat.virtual_dims[s][leg], false));
    }
    lp.push_back(rf::parity_vector(lat.physical_dims[s], false));
    if (fi.enabled) {
      if (all_even) {
        for (int leg = 0; leg < 4; ++leg) {
          fi.virt[s][leg] = lp[leg];
        }
        for (const bool p : fi.phys[s]) {
          REQUIRE(!p);
        }
      }
      lp = rf::Tn_parity(fi, s);
    }
    Tn[s] = lm_site_tensor<tensor>(lp, seed + 97u * s, all_even);
    if (fi.enabled) {
      REQUIRE(rf::parity_violation(lm_ft<tensor>{Tn[s], lp}) == 0.0);
    }
  }
  auto& lambda = lm_acc::lambda_tensor(state);
  lambda.assign(lat.N_UNIT, std::vector<std::vector<double>>(4));
  for (int s = 0; s < lat.N_UNIT; ++s) {
    const auto h = lm_lambda_vector(lat.virtual_dims[s][2], seed, s, 2);
    lambda[s][2] = h;
    lambda[lat.right(s)][0] = h;
    const auto v = lm_lambda_vector(lat.virtual_dims[s][1], seed, s, 1);
    lambda[s][1] = v;
    lambda[lat.top(s)][3] = v;
  }
}

// ---- windows ---------------------------------------------------------------

struct lm_window {
  int nrow = 0;
  int ncol = 0;
  int srow = 0;
  int scol = 0;
  int trow = 0;
  int tcol = 0;
  std::vector<std::vector<int>> idx;  //!< [row][col] -> site
  int source() const { return idx[srow][scol]; }
  int target() const { return idx[trow][tcol]; }
  bool outer(int row, int col, int leg) const {
    switch (leg) {
      case 0:
        return col == 0;
      case 1:
        return row == 0;
      case 2:
        return col == ncol - 1;
      case 3:
        return row == nrow - 1;
      default:
        return false;
    }
  }
};

inline lm_window lm_make_window(const tenes::SquareLattice& lat, int source,
                                int dx, int dy) {
  lm_window w;
  w.ncol = std::abs(dx) + 1;
  w.nrow = std::abs(dy) + 1;
  w.scol = dx >= 0 ? 0 : w.ncol - 1;
  w.tcol = w.ncol - 1 - w.scol;
  w.srow = dy >= 0 ? w.nrow - 1 : 0;
  w.trow = w.nrow - 1 - w.srow;
  w.idx.assign(w.nrow, std::vector<int>(w.ncol));
  for (int r = 0; r < w.nrow; ++r) {
    for (int c = 0; c < w.ncol; ++c) {
      w.idx[r][c] = lat.other(source, c - w.scol, w.srow - r);
    }
  }
  return w;
}

//! Deliberate misdressings of the truth, used only by the [truth]
//! discrimination case.
struct lm_dress {
  double power = 1.0;     //!< lambda^power (0.5: the weight on one layer only)
  bool internal = false;  //!< also dress the window's internal bonds
  bool none = false;      //!< no weight at all
};

//! The window's cells, wrapped, with the mean-field weights on the perimeter
//! legs.
template <class tensor>
std::vector<std::vector<lm_ft<tensor>>> lm_dressed_grid(lm_state<tensor>& state,
                                                        const lm_window& w,
                                                        lm_dress dress = {}) {
  const auto& fi = lm_acc::finfo(state);
  const auto& Tn = lm_acc::Tn(state);
  const auto& lambda = lm_acc::lambda_tensor(state);
  std::vector<std::vector<lm_ft<tensor>>> grid(w.nrow);
  for (int r = 0; r < w.nrow; ++r) {
    for (int c = 0; c < w.ncol; ++c) {
      const int site = w.idx[r][c];
      lm_ft<tensor> f = rf::wrap_Tn(Tn[site], fi, site);
      for (int leg = 0; leg < 4; ++leg) {
        if (dress.none || !(w.outer(r, c, leg) || dress.internal)) {
          continue;
        }
        std::vector<double> lam = lambda[site][leg];
        for (double& x : lam) {
          x = std::pow(x, dress.power);
        }
        f.t.multiply_vector(lam, leg);
      }
      grid[r].push_back(f);
    }
  }
  return grid;
}

// ---- the graded window truth
// -------------------------------------------------

//! Single-layer graded state of a dressed window: every cell contracted in
//! row-major order over the window's internal bonds (the left or upper
//! cell always the first operand, the orientation of the fold and of
//! build_pair_state()); the perimeter legs and the physical legs stay open.
template <class tensor>
struct lm_truth_state {
  lm_window w;
  lm_ft<tensor> psi;
  lm_ft<tensor> psi_conj;   //!< conj(psi), kept for the value of each operator
  std::vector<int> labels;  //!< cell * 5 + leg per axis, cell = row*ncol+col
  typename tensor::value_type norm = 0.0;
};

template <class tensor>
lm_truth_state<tensor> lm_truth_build(
    const lm_window& w, const std::vector<std::vector<lm_ft<tensor>>>& grid) {
  lm_truth_state<tensor> ts;
  ts.w = w;
  const int ncol = w.ncol;
  const auto partner = [&](int label) {
    const int cell = label / 5;
    const int leg = label % 5;
    const int r = cell / ncol;
    const int c = cell % ncol;
    int nr = r;
    int nc = c;
    switch (leg) {
      case 0:
        nc = c - 1;
        break;
      case 1:
        nr = r - 1;
        break;
      case 2:
        nc = c + 1;
        break;
      case 3:
        nr = r + 1;
        break;
      default:
        return -1;
    }
    if (nr < 0 || nr >= w.nrow || nc < 0 || nc >= ncol) {
      return -1;
    }
    return (nr * ncol + nc) * 5 + (leg + 2) % 4;
  };
  lm_ft<tensor> acc = grid[0][0];
  std::vector<int> labels{0, 1, 2, 3, 4};
  for (int cell = 1; cell < w.nrow * ncol; ++cell) {
    const lm_ft<tensor>& node = grid[cell / ncol][cell % ncol];
    mptensor::Axes axes_a;
    mptensor::Axes axes_b;
    std::vector<bool> used_a(labels.size(), false);
    std::vector<bool> used_b(5, false);
    for (int j = 0; j < 5; ++j) {
      const int p = partner(cell * 5 + j);
      if (p < 0) {
        continue;
      }
      const auto it = std::find(labels.begin(), labels.end(), p);
      if (it == labels.end()) {
        continue;
      }
      const std::size_t k = static_cast<std::size_t>(it - labels.begin());
      axes_a.push(static_cast<int>(k));
      axes_b.push(j);
      used_a[k] = true;
      used_b[j] = true;
    }
    std::vector<int> next;
    for (std::size_t k = 0; k < labels.size(); ++k) {
      if (!used_a[k]) {
        next.push_back(labels[k]);
      }
    }
    for (int j = 0; j < 5; ++j) {
      if (!used_b[j]) {
        next.push_back(cell * 5 + j);
      }
    }
    acc = rf::tensordot(acc, node, axes_a, axes_b);
    labels = next;
  }
  ts.psi = acc;
  ts.labels = labels;
  mptensor::Axes all;
  for (int ax = 0; ax < ts.psi.rank(); ++ax) {
    all.push(ax);
  }
  ts.psi_conj = rf::conj(ts.psi);
  ts.norm = rf::trace(ts.psi_conj, ts.psi, all, all);
  REQUIRE(std::abs(ts.norm) > 0.0);
  return ts;
}

template <class tensor>
lm_truth_state<tensor> lm_truth_build(lm_state<tensor>& state, int source,
                                      int dx, int dy, lm_dress dress = {}) {
  const lm_window w = lm_make_window(lm_acc::lattice(state), source, dx, dy);
  return lm_truth_build<tensor>(w, lm_dressed_grid(state, w, dress));
}

//! Normalized <op4> (source first, loaded with wrap_twosite_gate()). The
//! cells stay in row-major (= Jordan-Wigner) order; a source that comes
//! later is expressed on the gate by the graded transpose (1, 0, 3, 2), the
//! convention of relay_window.cpp's truth and of the existing
//! nearest-neighbour path.
template <class tensor>
typename tensor::value_type lm_truth_value(const lm_truth_state<tensor>& ts,
                                           const tensor& op4,
                                           const rf::parity_vector& phys_s,
                                           const rf::parity_vector& phys_t) {
  const lm_window& w = ts.w;
  const int s = w.srow * w.ncol + w.scol;
  const int t = w.trow * w.ncol + w.tcol;
  lm_ft<tensor> gate = rf::wrap_twosite_gate(op4, phys_s, phys_t);
  int a = s;
  int b = t;
  if (s > t) {
    a = t;
    b = s;
    gate = rf::transpose(gate, mptensor::Axes(1, 0, 3, 2));
  }
  const int rank = ts.psi.rank();
  int ax_a = -1;
  int ax_b = -1;
  for (int k = 0; k < rank; ++k) {
    if (ts.labels[k] == a * 5 + 4) {
      ax_a = k;
    }
    if (ts.labels[k] == b * 5 + 4) {
      ax_b = k;
    }
  }
  REQUIRE(ax_a >= 0);
  REQUIRE(ax_b >= 0);
  lm_ft<tensor> applied = rf::tensordot(
      ts.psi, gate, mptensor::Axes(ax_a, ax_b), mptensor::Axes(0, 1));
  mptensor::Axes perm;
  int run = 0;
  for (int ax = 0; ax < rank; ++ax) {
    if (ax == ax_a) {
      perm.push(rank - 2);
    } else if (ax == ax_b) {
      perm.push(rank - 1);
    } else {
      perm.push(run++);
    }
  }
  applied = rf::transpose(applied, perm);
  mptensor::Axes all;
  for (int ax = 0; ax < rank; ++ax) {
    all.push(ax);
  }
  return rf::trace(ts.psi_conj, applied, all, all) / ts.norm;
}

//! contract_pair_MF() of the dressed nearest-neighbour window (the existing
//! mean-field path of twosite_obs.cpp, assembled here), normalized.
template <class tensor>
typename tensor::value_type lm_pair_mf_value(
    const lm_window& w, const std::vector<std::vector<lm_ft<tensor>>>& grid,
    const tensor& op4, const rf::parity_vector& phys_s,
    const rf::parity_vector& phys_t) {
  REQUIRE(w.nrow * w.ncol == 2);
  const bool horizontal = w.nrow == 1;
  const lm_ft<tensor>& A = grid[0][0];
  const lm_ft<tensor>& B = horizontal ? grid[0][1] : grid[1][0];
  const bool source_first = w.srow == 0 && w.scol == 0;
  lm_ft<tensor> o = rf::wrap_twosite_gate(op4, phys_s, phys_t);
  if (!source_first) {
    o = rf::transpose(o, mptensor::Axes(1, 0, 3, 2));
  }
  const auto pair =
      rf::build_pair_state(A, B,
                           horizontal ? rf::reduced_pair_direction::horizontal
                                      : rf::reduced_pair_direction::vertical);
  return rf::contract_pair_MF(pair, o) / rf::contract_pair_MF(pair);
}

// ---- the delta environment (contract item 3, T5-3a) -------------------------

//! sum_k [window of channel k] / [identity window], each closed with the
//! CHI = 1 environment of make_delta_corner() / make_delta_edge() through
//! core::Contract_density_CTM(): the construction of design section 6,
//! assembled here from T1's relay window.
template <class tensor>
typename tensor::value_type lm_delta_value(
    const lm_window& w, const std::vector<std::vector<lm_ft<tensor>>>& grid,
    const tensor& op4, const rf::parity_vector& phys_s,
    const rf::parity_vector& phys_t) {
  const auto comm = grid[0][0].t.get_comm();
  tensor C;
  REQUIRE_NOTHROW(C = rf::make_delta_corner<tensor>(comm));
  REQUIRE(C.shape() == mptensor::Shape(1, 1));
  const auto edge = [&](int D) {
    tensor e;
    REQUIRE_NOTHROW(e = rf::make_delta_edge<tensor>(D, comm));
    REQUIRE(e.shape() == mptensor::Shape(1, 1, D * D));
    return e;
  };
  std::vector<tensor> et(w.ncol), eb(w.ncol), el(w.nrow), er(w.nrow);
  for (int c = 0; c < w.ncol; ++c) {
    et[c] = edge(static_cast<int>(grid[0][c].shape()[1]));
    eb[c] = edge(static_cast<int>(grid[w.nrow - 1][c].shape()[3]));
  }
  for (int r = 0; r < w.nrow; ++r) {
    el[r] = edge(static_cast<int>(grid[r][0].shape()[0]));
    er[r] = edge(static_cast<int>(grid[r][w.ncol - 1].shape()[2]));
  }
  const std::vector<const tensor*> Cs{&C, &C, &C, &C};
  std::vector<const tensor*> pt, pb, pl, pr;
  for (auto& e : et) {
    pt.push_back(&e);
  }
  for (auto& e : eb) {
    pb.push_back(&e);
  }
  for (auto& e : el) {
    pl.push_back(&e);
  }
  for (auto& e : er) {
    pr.push_back(&e);
  }
  std::vector<std::vector<tensor>> ident(w.nrow);
  std::vector<std::vector<tensor>> plain(w.nrow);
  for (int r = 0; r < w.nrow; ++r) {
    for (int c = 0; c < w.ncol; ++c) {
      ident[r].push_back(
          lm_identity<tensor>(static_cast<int>(grid[r][c].shape()[4])));
      plain[r].push_back(rf::detail::doubled_pipeline(grid[r][c], grid[r][c]));
    }
  }
  const auto pointers = [](const std::vector<std::vector<tensor>>& v) {
    std::vector<std::vector<const tensor*>> p(v.size());
    for (std::size_t r = 0; r < v.size(); ++r) {
      for (const tensor& t : v[r]) {
        p[r].push_back(&t);
      }
    }
    return p;
  };
  const auto op_ptr = pointers(ident);
  const auto norm = tenes::itps::core::Contract_density_CTM(
      Cs, pt, pr, pb, pl, pointers(plain), op_ptr);
  REQUIRE(std::abs(norm) > 0.0);
  const auto channels =
      rf::relay_channels(rf::wrap_twosite_gate(op4, phys_s, phys_t));
  REQUIRE(!channels.empty());
  typename tensor::value_type sum = 0.0;
  for (const auto& ch : channels) {
    const auto window =
        rf::build_relay_window(grid, rf::window_cell{w.srow, w.scol},
                               rf::window_cell{w.trow, w.tcol}, ch);
    sum += tenes::itps::core::Contract_density_CTM(Cs, pt, pr, pb, pl,
                                                   pointers(window), op_ptr);
  }
  return sum / norm;
}

struct lm_disp {
  int dx;
  int dy;
};

const std::vector<lm_disp> lm_nn_disps = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};

inline std::string lm_bond_label(int source, int dx, int dy) {
  return "source " + std::to_string(source) + " (dx,dy)=(" +
         std::to_string(dx) + "," + std::to_string(dy) + ")";
}

//! Contract items 3 and 4: relative tolerance of the gates.
constexpr double lm_rtol_gate = 1.0e-12;
//! Contract item 5.
constexpr double lm_rtol_window = 1.0e-10;
//! Premise: a compared value carries signal far above the tolerance.
constexpr double lm_min_signal = 1.0e-3;

}  // namespace

// ============================================================================
// Contract item 2 / Interfaces: the delta environment
// ============================================================================

namespace {

template <class tensor>
void lm_run_delta_corner() {
  INFO("tensor type " << lm_type_name<tensor>());
  tensor C;
  REQUIRE_NOTHROW(C = rf::make_delta_corner<tensor>(MPI_COMM_WORLD));
  REQUIRE(C.shape() == mptensor::Shape(1, 1));
  typename tensor::value_type v;
  C.get_value(mptensor::Index(0, 0), v);
  CHECK(v == typename tensor::value_type(1.0));
}

template <class tensor>
void lm_run_delta_edge(int D) {
  INFO("tensor type " << lm_type_name<tensor>() << ", D = " << D);
  tensor e;
  REQUIRE_NOTHROW(e = rf::make_delta_edge<tensor>(D, MPI_COMM_WORLD));
  REQUIRE(e.shape() == mptensor::Shape(1, 1, D * D));
  // Every element: the fused index f = x + D * xb (x fastest, as
  // detail::doubled_pipeline() fuses [x xb]) is 1 exactly when x == xb.
  int ones = 0;
  for (std::size_t n = 0; n < e.local_size(); ++n) {
    const mptensor::Index idx = e.global_index(n);
    const int f = static_cast<int>(idx[2]);
    const int x = f % D;
    const int xb = f / D;
    const typename tensor::value_type want = x == xb ? 1.0 : 0.0;
    INFO("element (0, 0, " << f << ") = [x=" << x << ", xb=" << xb << "]");
    CHECK(e[n] == want);
    if (x == xb) {
      ++ones;
    }
  }
  if (e.get_comm_size() == 1) {
    CHECK(ones == D);
  }
}

}  // namespace

TEST_CASE(
    "longrange-mf T5-2a: make_delta_corner is the CHI = 1 corner, shape "
    "(1, 1) and value 1") {
  lm_run_delta_corner<real_tensor>();
  lm_run_delta_corner<complex_tensor>();
}

TEST_CASE(
    "longrange-mf T5-2b: make_delta_edge closes the fused leg [x xb] with "
    "delta_{x, xb}, shape (1, 1, D * D)") {
  for (int D = 1; D <= 4; ++D) {
    lm_run_delta_edge<real_tensor>(D);
    lm_run_delta_edge<complex_tensor>(D);
  }
}

// ============================================================================
// [truth]: the graded window truth
// ============================================================================

namespace {

//! Dressed nearest-neighbour windows of a fermion solver state: the graded
//! window truth against contract_pair_MF(), every direction.
template <class tensor>
void lm_run_truth_nn(const lm_geom& g, const std::string& kind,
                     std::uint64_t seed) {
  const std::string label = lm_geom_label(g) + " " + kind + " " +
                            lm_type_name<tensor>() + " seed " +
                            std::to_string(seed);
  INFO(label);
  auto state =
      lm_make<tensor>(g, true, lm_phys(g.d), {}, {}, CorrelationParameter{});
  lm_seed(*state, seed, false);
  const tensor op4 = lm_kind_tensor<tensor>(kind, g.d);
  const double scale = lm_max_abs_entry(op4);
  const auto& lat = lm_acc::lattice(*state);
  for (int source = 0; source < lat.N_UNIT; ++source) {
    for (const lm_disp dd : lm_nn_disps) {
      const std::string what =
          label + " " + lm_bond_label(source, dd.dx, dd.dy);
      const lm_window w = lm_make_window(lat, source, dd.dx, dd.dy);
      const auto grid = lm_dressed_grid(*state, w);
      const auto phys = lm_phys(g.d);
      const auto truth =
          lm_truth_value(lm_truth_build<tensor>(w, grid), op4, phys, phys);
      const auto pair = lm_pair_mf_value(w, grid, op4, phys, phys);
      lm_check_close(what + " [graded window truth vs contract_pair_MF]", truth,
                     pair, scale, 1.0e-13);
    }
  }
}

}  // namespace

TEST_CASE(
    "longrange-mf [truth] the graded window truth reproduces "
    "contract_pair_MF at the four nearest-neighbour directions") {
  lm_run_truth_nn<real_tensor>({3, 2, 2, 2, 3}, "asym", 501);
  lm_run_truth_nn<complex_tensor>({3, 2, 2, 2, 3}, "cplx", 502);
  lm_run_truth_nn<real_tensor>({2, 2, 4, 3, 2}, "asym4", 503);
  lm_run_truth_nn<complex_tensor>({2, 2, 4, 2, 2}, "cplx", 504);
}

namespace {

//! Contract item 4's premise, and the lambda placement of the truth: with
//! every parity even the graded truth is the bosonic mean-field value, both
//! for windows (Contract_iTPS_MF) and for the correlation chain
//! (*Correlation_iTPS_MF), as measured by a bosonic solver.
template <class tensor>
void lm_run_truth_boson(std::uint64_t seed) {
  using value_type = typename tensor::value_type;
  const lm_geom g{3, 3, 2, 2, 2};
  const std::string label =
      lm_geom_label(g) + " all even " + lm_type_name<tensor>();
  INFO(label);
  const rf::parity_vector even{false, false};
  // A dense one-site pair and a dense two-site operator.
  auto dense_onesite = [&](std::uint64_t s) {
    tensor t(mptensor::Shape(2, 2));
    for (std::size_t n = 0; n < t.local_size(); ++n) {
      lm_set_hashed(t, t.global_index(n), s, 1.0);
    }
    return t;
  };
  const tensor A = dense_onesite(seed + 11);
  const tensor B = dense_onesite(seed + 12);
  tensor op4(mptensor::Shape(2, 2, 2, 2));
  for (std::size_t n = 0; n < op4.local_size(); ++n) {
    lm_set_hashed(op4, op4.global_index(n), seed + 13, 1.0);
  }
  // A_s B_t with every parity even: the plain outer product.
  tensor AB(mptensor::Shape(2, 2, 2, 2));
  for (int is = 0; is < 2; ++is) {
    for (int it = 0; it < 2; ++it) {
      for (int os = 0; os < 2; ++os) {
        for (int ot = 0; ot < 2; ++ot) {
          value_type a, b;
          A.get_value(mptensor::Index(is, os), a);
          B.get_value(mptensor::Index(it, ot), b);
          AB.set_value(mptensor::Index(is, it, os, ot), a * b);
        }
      }
    }
  }
  // Windows with at least two rows and two columns only: the bosonic
  // mean-field branch of measure_twosite() asserts
  // boundaries.size() == 2 * (ncol + nrow - 2) (twosite_obs.cpp), which a
  // straight window of three or more sites violates, and a Debug build
  // aborts there (a bug since 2021, reported with T5; Release builds skip
  // the assert and measure correctly). The straight chains are pinned below
  // through the correlation function, whose kernels do not have it.
  const std::vector<lm_disp> disps = {{1, 1}, {-1, 2}, {2, -1}, {-1, -1}};
  const std::vector<int> sources = {0, 4};
  tenes::Operators<tensor> onesite;
  for (int s = 0; s < 9; ++s) {
    onesite.emplace_back("a", 0, s, A);
    onesite.emplace_back("b", 1, s, B);
  }
  tenes::Operators<tensor> twosite;
  for (const int s : sources) {
    for (const lm_disp dd : disps) {
      twosite.emplace_back("op", 0, s, dd.dx, dd.dy, op4);
    }
  }
  auto boson = lm_make<tensor>(g, false, even, onesite, twosite,
                               CorrelationParameter(3, {{0, 1}}));
  lm_seed(*boson, seed, true);
  // The truth needs a fermion solver's ledgers (all even) for wrap_Tn().
  auto fermion =
      lm_make<tensor>(g, true, even, onesite, {}, CorrelationParameter{});
  lm_seed(*fermion, seed, true);
  for (int s = 0; s < 9; ++s) {
    REQUIRE(lm_acc::Tn(*fermion)[s].shape() == lm_acc::Tn(*boson)[s].shape());
  }

  const auto measured = boson->measure_twosite();
  REQUIRE(!measured.empty());
  const double scale = lm_max_abs_entry(op4);
  for (const int s : sources) {
    for (const lm_disp dd : disps) {
      const std::string what = label + " " + lm_bond_label(s, dd.dx, dd.dy);
      REQUIRE(measured[0].count(Bond{s, dd.dx, dd.dy}) == 1);
      const auto truth = lm_truth_value(
          lm_truth_build(*fermion, s, dd.dx, dd.dy), op4, even, even);
      lm_check_close(what + " [truth vs bosonic Contract_iTPS_MF]",
                     measured[0].at(Bond{s, dd.dx, dd.dy}), truth, scale,
                     lm_rtol_gate);
    }
  }
  const auto rows = boson->measure_correlation();
  REQUIRE(rows.size() == 9u * 3u * 2u);
  const double cscale = lm_max_abs_entry(A) * lm_max_abs_entry(B);
  for (const Correlation& row : rows) {
    if (std::find(sources.begin(), sources.end(), row.left_index) ==
        sources.end()) {
      continue;  // the truth costs a window contraction per row
    }
    const std::string what =
        label + " correlation " +
        lm_bond_label(row.left_index, row.right_dx, row.right_dy);
    const auto truth = lm_truth_value(
        lm_truth_build(*fermion, row.left_index, row.right_dx, row.right_dy),
        AB, even, even);
    lm_check_close(what + " [truth vs bosonic *Correlation_iTPS_MF]",
                   lm_cplx(row.real, row.imag), truth, cscale, lm_rtol_gate);
  }
}

}  // namespace

TEST_CASE(
    "longrange-mf [truth] with every parity even the graded window truth is "
    "the bosonic mean-field value, for windows and correlation chains") {
  lm_run_truth_boson<real_tensor>(601);
  lm_run_truth_boson<complex_tensor>(602);
}

namespace {

template <class tensor>
void lm_run_truth_discrimination(const lm_geom& g, const std::string& kind,
                                 std::uint64_t seed) {
  const std::string label =
      lm_geom_label(g) + " " + kind + " " + lm_type_name<tensor>();
  INFO(label);
  auto state =
      lm_make<tensor>(g, true, lm_phys(g.d), {}, {}, CorrelationParameter{});
  lm_seed(*state, seed, false);
  const tensor op4 = lm_kind_tensor<tensor>(kind, g.d);
  const tensor odd4 = lm_twosite_tensor<tensor>(lm_odd_terms(kind, g.d), g.d);
  const double scale = lm_max_abs_entry(op4);
  const double tol = lm_rtol_gate * scale;
  const auto phys = lm_phys(g.d);
  for (const lm_disp dd :
       std::vector<lm_disp>{{1, 0}, {0, -1}, {2, 0}, {1, 1}, {0, 2}}) {
    const std::string what = label + " " + lm_bond_label(0, dd.dx, dd.dy);
    INFO(what);
    const auto ts = lm_truth_build(*state, 0, dd.dx, dd.dy);
    const auto truth = lm_truth_value(ts, op4, phys, phys);
    // The odd x odd channels carry signal: dropping them or flipping their
    // sign moves the value by far more than the tolerance.
    const auto odd = lm_truth_value(ts, odd4, phys, phys);
    INFO("odd part " << odd << ", tolerance " << tol);
    CHECK(std::abs(odd) > lm_min_signal * scale);
    for (const auto& mut : std::vector<std::pair<std::string, lm_dress>>{
             {"lambda on one layer only", {0.5, false, false}},
             {"lambda on the internal bonds too", {1.0, true, false}},
             {"no lambda", {1.0, false, true}}}) {
      const auto bad = lm_truth_value(
          lm_truth_build(*state, 0, dd.dx, dd.dy, mut.second), op4, phys, phys);
      INFO(mut.first << ": " << bad << " against " << truth);
      CHECK(std::abs(bad - truth) > 1.0e6 * tol);
    }
  }
}

}  // namespace

TEST_CASE(
    "longrange-mf [truth] the states used make a misplaced lambda and a lost "
    "odd channel visible far above the tolerance") {
  lm_run_truth_discrimination<real_tensor>({3, 3, 2, 2, 3}, "asym", 701);
  lm_run_truth_discrimination<complex_tensor>({2, 2, 4, 2, 2}, "cplx", 702);
}

// ============================================================================
// Contract item 3 (gate 1): the delta closure against contract_pair_MF
// ============================================================================

namespace {

template <class tensor>
void lm_run_gate1_construction(const lm_geom& g, const std::string& kind,
                               std::uint64_t seed) {
  const std::string label = lm_geom_label(g) + " " + kind + " " +
                            lm_type_name<tensor>() + " seed " +
                            std::to_string(seed);
  INFO(label);
  auto state =
      lm_make<tensor>(g, true, lm_phys(g.d), {}, {}, CorrelationParameter{});
  lm_seed(*state, seed, false);
  const tensor op4 = lm_kind_tensor<tensor>(kind, g.d);
  const double scale = lm_max_abs_entry(op4);
  const auto phys = lm_phys(g.d);
  const auto& lat = lm_acc::lattice(*state);
  for (const int source : {0, lat.N_UNIT - 1}) {
    for (const lm_disp dd : lm_nn_disps) {
      const std::string what =
          label + " " + lm_bond_label(source, dd.dx, dd.dy);
      INFO(what);
      const lm_window w = lm_make_window(lat, source, dd.dx, dd.dy);
      const auto grid = lm_dressed_grid(*state, w);
      const auto want = lm_pair_mf_value(w, grid, op4, phys, phys);
      const auto got = lm_delta_value(w, grid, op4, phys, phys);
      lm_check_close(what + " [delta closure vs contract_pair_MF]", got, want,
                     scale, lm_rtol_gate);
    }
  }
}

}  // namespace

TEST_CASE(
    "longrange-mf T5-3a: gate 1, the lambda-dressed relay window closed by "
    "the delta environment in the density kernel equals contract_pair_MF at "
    "(1,0), (0,1), (-1,0), (0,-1) (d = 2 and 4, real and complex)") {
  lm_run_gate1_construction<real_tensor>({3, 2, 2, 2, 3}, "asym", 801);
  lm_run_gate1_construction<complex_tensor>({3, 2, 2, 3, 2}, "cplx", 802);
  lm_run_gate1_construction<real_tensor>({2, 2, 4, 2, 3}, "asym4", 803);
  lm_run_gate1_construction<complex_tensor>({2, 2, 4, 2, 2}, "cplx", 804);
}

// ============================================================================
// Contract item 3 (gate 1) through the solver: the correlation function
// ============================================================================

namespace {

struct lm_corr_case {
  lm_geom g;
  std::vector<std::string> names;
  std::vector<std::pair<int, int>> pairs;
  int r_max;
  std::uint64_t seed;
  //! Left sites whose rows are compared with the graded window truth (it
  //! costs a window contraction per row); every row at r = 1 is compared
  //! with contract_pair_MF.
  std::vector<int> truth_sites;
};

inline CorrelationParameter lm_corparam(const lm_corr_case& c) {
  std::vector<std::tuple<int, int>> ops;
  for (const auto& p : c.pairs) {
    ops.emplace_back(p.first, p.second);
  }
  return CorrelationParameter(c.r_max, ops);
}

template <class tensor>
tenes::Operators<tensor> lm_onesite_groups(
    const std::vector<std::string>& names, int d, int n_unit) {
  tenes::Operators<tensor> onesite;
  for (std::size_t g = 0; g < names.size(); ++g) {
    const tensor op = lm_onesite_tensor<tensor>(lm_onesite_op(names[g], d));
    for (int s = 0; s < n_unit; ++s) {
      onesite.emplace_back(names[g], static_cast<int>(g), s, op);
    }
  }
  return onesite;
}

inline bool lm_pair_even(const lm_corr_case& c, int lop, int rop) {
  return lm_onesite_op(c.names[lop], c.g.d).odd ==
         lm_onesite_op(c.names[rop], c.g.d).odd;
}

template <class tensor>
void lm_run_correlation_truth(const lm_corr_case& c) {
  const std::string label =
      lm_geom_label(c.g) + " " + lm_type_name<tensor>() + " correlation";
  INFO(label);
  const int n_unit = c.g.lx * c.g.ly;
  auto state = lm_make<tensor>(
      c.g, true, lm_phys(c.g.d),
      lm_onesite_groups<tensor>(c.names, c.g.d, n_unit), {}, lm_corparam(c));
  lm_seed(*state, c.seed, false);
  std::vector<Correlation> rows;
  REQUIRE_NOTHROW(rows = state->measure_correlation());
  CHECK(rows.size() == static_cast<std::size_t>(n_unit) * c.pairs.size() *
                           static_cast<std::size_t>(c.r_max) * 2u);
  const auto phys = lm_phys(c.g.d);
  const auto& lat = lm_acc::lattice(*state);
  std::map<std::tuple<int, int, int>, lm_truth_state<tensor>> cache;
  std::size_t checked = 0;
  double max_signal = 0.0;
  for (const Correlation& row : rows) {
    const int r = std::max(std::abs(row.right_dx), std::abs(row.right_dy));
    const std::string what =
        label + " [" + c.names[row.left_op] + ", " + c.names[row.right_op] +
        "] " + lm_bond_label(row.left_index, row.right_dx, row.right_dy);
    INFO(what);
    REQUIRE(row.left_op >= 0);
    REQUIRE(row.right_op >= 0);
    REQUIRE(row.left_op < static_cast<int>(c.names.size()));
    REQUIRE(row.right_op < static_cast<int>(c.names.size()));
    const lm_cplx got(row.real, row.imag);
    if (!lm_pair_even(c, row.left_op, row.right_op)) {
      // Different parity: exactly 0, not contracted.
      CHECK(row.real == 0.0);
      CHECK(row.imag == 0.0);
      continue;
    }
    const tensor op4 = lm_twosite_tensor<tensor>(
        {{1.0, c.names[row.left_op], c.names[row.right_op]}}, c.g.d);
    const double scale = lm_max_abs_entry(op4);
    if (std::find(c.truth_sites.begin(), c.truth_sites.end(), row.left_index) !=
        c.truth_sites.end()) {
      const auto key =
          std::make_tuple(row.left_index, row.right_dx, row.right_dy);
      if (cache.count(key) == 0) {
        cache.emplace(key, lm_truth_build(*state, row.left_index, row.right_dx,
                                          row.right_dy));
      }
      const auto truth = lm_truth_value(cache.at(key), op4, phys, phys);
      lm_check_close(what + " [correlation vs graded window truth]", got, truth,
                     scale, lm_rtol_gate);
      max_signal = std::max(max_signal, std::abs(truth) / scale);
      ++checked;
    }
    if (r == 1) {
      // Gate 1 as the contract states it: the existing contract_pair_MF.
      const lm_window w =
          lm_make_window(lat, row.left_index, row.right_dx, row.right_dy);
      const auto pair =
          lm_pair_mf_value(w, lm_dressed_grid(*state, w), op4, phys, phys);
      lm_check_close(what + " [r = 1 correlation vs contract_pair_MF]", got,
                     pair, scale, lm_rtol_gate);
      max_signal = std::max(max_signal, std::abs(pair) / scale);
      ++checked;
    }
  }
  CHECK(checked > 0);
  INFO("largest |value| / scale " << max_signal);
  CHECK(max_signal > lm_min_signal);
}

const std::vector<std::string> lm_names2 = {"n", "cdag", "c"};
const std::vector<std::string> lm_names4 = {"n", "cdag_up", "c_up", "cdag_dn",
                                            "c_dn"};

}  // namespace

TEST_CASE(
    "longrange-mf T5-3b: gate 1 through the solver, the mean-field "
    "correlation function at r = 1 equals contract_pair_MF and at r = 2, 3 "
    "the graded window truth (d = 2 and 4, real and complex)") {
  // [n,n], [cdag,c], [c,cdag], [c,c] and the mixed [n,c] (exactly 0).
  const lm_corr_case c2{{3, 2, 2, 3, 2},
                        lm_names2,
                        {{0, 0}, {1, 2}, {2, 2}, {0, 2}},
                        3,
                        901,
                        {0}};
  lm_run_correlation_truth<real_tensor>(c2);
  lm_corr_case c2c = c2;
  c2c.g = {2, 3, 2, 2, 2};
  c2c.pairs = {{2, 1}, {0, 0}, {1, 1}};
  c2c.seed = 902;
  c2c.truth_sites = {0, 5};
  lm_run_correlation_truth<complex_tensor>(c2c);
  const lm_corr_case c4{{2, 2, 4, 2, 2},
                        lm_names4,
                        {{1, 2}, {4, 3}, {0, 0}, {2, 0}},
                        3,
                        903,
                        {0}};
  lm_run_correlation_truth<real_tensor>(c4);
  lm_corr_case c4c = c4;
  c4c.seed = 904;
  c4c.truth_sites = {3};
  lm_run_correlation_truth<complex_tensor>(c4c);
}

// ============================================================================
// Contract items 2, 3: long-range windows and ops forms through the solver
// ============================================================================

namespace {

//! One two-site group of a window case: an explicit kind with a
//! coefficient, or (ops_a, ops_b non-empty) the ops form of two one-site
//! groups.
struct lm_group {
  std::string kind;
  lm_cplx coeff;
  int ops_a = -1;
  int ops_b = -1;
};

struct lm_window_case {
  lm_geom g;
  std::vector<std::string> names;
  std::vector<lm_group> groups;
  std::vector<int> sources;
  std::vector<lm_disp> disps;
  std::uint64_t seed;
};

template <class tensor>
void lm_run_window_truth(const lm_window_case& c) {
  const std::string label =
      lm_geom_label(c.g) + " " + lm_type_name<tensor>() + " windows";
  INFO(label);
  const int n_unit = c.g.lx * c.g.ly;
  const int d = c.g.d;
  tenes::Operators<tensor> twosite;
  std::vector<tensor> op4s;
  for (std::size_t gi = 0; gi < c.groups.size(); ++gi) {
    const lm_group& gr = c.groups[gi];
    if (gr.ops_a >= 0) {
      op4s.push_back(lm_twosite_tensor<tensor>(
          {{1.0, c.names[gr.ops_a], c.names[gr.ops_b]}}, d));
    } else {
      op4s.push_back(lm_kind_tensor<tensor>(gr.kind, d));
    }
    for (const int s : c.sources) {
      for (const lm_disp dd : c.disps) {
        if (gr.ops_a >= 0) {
          twosite.emplace_back("ops", static_cast<int>(gi), s, dd.dx, dd.dy,
                               std::vector<int>{gr.ops_a, gr.ops_b});
        } else {
          twosite.emplace_back(gr.kind, static_cast<int>(gi), s, dd.dx, dd.dy,
                               op4s.back(), lm_value<tensor>(gr.coeff));
        }
      }
    }
  }
  auto state = lm_make<tensor>(c.g, true, lm_phys(d),
                               lm_onesite_groups<tensor>(c.names, d, n_unit),
                               twosite, CorrelationParameter{});
  lm_seed(*state, c.seed, false);
  std::vector<std::map<Bond, typename tensor::value_type>> measured;
  REQUIRE_NOTHROW(measured = state->measure_twosite());
  REQUIRE(measured.size() >= c.groups.size());
  const auto phys = lm_phys(d);
  for (const int s : c.sources) {
    for (const lm_disp dd : c.disps) {
      const auto ts = lm_truth_build(*state, s, dd.dx, dd.dy);
      for (std::size_t gi = 0; gi < c.groups.size(); ++gi) {
        const lm_group& gr = c.groups[gi];
        const lm_cplx coeff = gr.ops_a >= 0 ? lm_cplx(1.0) : gr.coeff;
        const std::string what =
            label + " group " +
            (gr.ops_a >= 0
                 ? "ops [" + c.names[gr.ops_a] + ", " + c.names[gr.ops_b] + "]"
                 : gr.kind) +
            " " + lm_bond_label(s, dd.dx, dd.dy);
        INFO(what);
        REQUIRE(measured[gi].count(Bond{s, dd.dx, dd.dy}) == 1);
        const lm_cplx want =
            coeff * lm_cplx(lm_truth_value(ts, op4s[gi], phys, phys));
        lm_check_close(what + " [measure_twosite vs graded window truth]",
                       measured[gi].at(Bond{s, dd.dx, dd.dy}), want,
                       std::abs(coeff) * lm_max_abs_entry(op4s[gi]),
                       lm_rtol_gate);
      }
    }
  }
}

}  // namespace

TEST_CASE(
    "longrange-mf T5-3c: long-range two-site observables and ops forms with "
    "meanfield_env equal the graded window truth (d = 2, real)") {
  lm_run_window_truth<real_tensor>(
      {{3, 3, 2, 3, 2},
       lm_names2,
       {{"asym", 1.0}, {"hop", 0.7}, {"", 1.0, 1, 2}, {"", 1.0, 2, 1}},
       {0, 5},
       {{2, 0}, {0, -2}, {1, 1}, {-1, -1}, {3, 0}},
       1001});
}

TEST_CASE(
    "longrange-mf T5-3c: long-range two-site observables and ops forms with "
    "meanfield_env equal the graded window truth (d = 2, complex)") {
  lm_run_window_truth<complex_tensor>(
      {{3, 3, 2, 2, 2},
       lm_names2,
       {{"cplx", lm_alpha}, {"", 1.0, 1, 2}},
       {0, 4},
       {{2, 0}, {0, 2}, {1, 1}, {-1, 1}, {2, -1}, {-2, 0}, {0, -3}},
       1002});
}

TEST_CASE(
    "longrange-mf T5-3c: long-range two-site observables and ops forms with "
    "meanfield_env equal the graded window truth (d = 4, real and "
    "complex)") {
  lm_run_window_truth<real_tensor>(
      {{2, 2, 4, 2, 2},
       lm_names4,
       {{"asym4", 1.0}, {"", 1.0, 1, 2}, {"", 1.0, 4, 3}},
       {0, 3},
       {{2, 0}, {1, 1}, {0, -2}, {-1, 1}},
       1003});
  lm_run_window_truth<complex_tensor>({{2, 2, 4, 2, 2},
                                       lm_names4,
                                       {{"cplx", lm_beta}, {"", 1.0, 3, 4}},
                                       {1},
                                       {{-2, 0}, {1, -1}, {0, 3}},
                                       1004});
}

// ============================================================================
// Contract item 4 (gate 2): every parity even, against the bosonic solver
// ============================================================================

namespace {

template <class tensor>
void lm_run_gate2(std::uint64_t seed) {
  using value_type = typename tensor::value_type;
  const lm_geom g{3, 2, 2, 2, 3};
  const std::string label =
      lm_geom_label(g) + " all even " + lm_type_name<tensor>();
  INFO(label);
  const rf::parity_vector even{false, false};
  const int n_unit = g.lx * g.ly;
  auto dense = [&](const mptensor::Shape& sh, std::uint64_t s) {
    tensor t(sh);
    for (std::size_t n = 0; n < t.local_size(); ++n) {
      lm_set_hashed(t, t.global_index(n), s, 1.0);
    }
    return t;
  };
  const tensor A = dense(mptensor::Shape(2, 2), seed + 1);
  const tensor B = dense(mptensor::Shape(2, 2), seed + 2);
  const tensor op4 = dense(mptensor::Shape(2, 2, 2, 2), seed + 3);
  const value_type coeff = lm_value<tensor>(lm_cplx(0.9, -0.4));
  tenes::Operators<tensor> onesite;
  for (int s = 0; s < n_unit; ++s) {
    onesite.emplace_back("a", 0, s, A);
    onesite.emplace_back("b", 1, s, B);
  }
  const std::vector<lm_disp> disps = {{2, 0},  {1, 1},  {3, 0},
                                      {0, -2}, {-1, 2}, {-2, -1}};
  const std::vector<int> sources = {0, 4};
  tenes::Operators<tensor> twosite;
  for (const int s : sources) {
    for (const lm_disp dd : disps) {
      twosite.emplace_back("op", 0, s, dd.dx, dd.dy, op4, coeff);
      twosite.emplace_back("ops", 1, s, dd.dx, dd.dy, std::vector<int>{0, 1});
    }
  }
  const CorrelationParameter cp(3, {{0, 1}, {1, 0}, {1, 1}});
  auto boson = lm_make<tensor>(g, false, even, onesite, twosite, cp);
  auto fermion = lm_make<tensor>(g, true, even, onesite, twosite, cp);
  lm_seed(*boson, seed, true);
  lm_seed(*fermion, seed, true);

  // The fermion solver first: against the stub it stops here, and the
  // bosonic reference below would abort a Debug build on the straight
  // windows (the assert in twosite_obs.cpp's mean-field branch, which the
  // fermionic mean-field path shares and has to fix; see
  // lm_run_truth_boson()).
  std::vector<std::map<Bond, value_type>> got;
  REQUIRE_NOTHROW(got = fermion->measure_twosite());
  const auto want = boson->measure_twosite();
  REQUIRE(want.size() >= 2);
  REQUIRE(got.size() >= 2);
  const double scales[2] = {std::abs(coeff) * lm_max_abs_entry(op4),
                            lm_max_abs_entry(A) * lm_max_abs_entry(B)};
  for (int gi = 0; gi < 2; ++gi) {
    for (const int s : sources) {
      for (const lm_disp dd : disps) {
        const Bond b{s, dd.dx, dd.dy};
        const std::string what = label + " group " + std::to_string(gi) + " " +
                                 lm_bond_label(s, dd.dx, dd.dy);
        INFO(what);
        REQUIRE(want[gi].count(b) == 1);
        REQUIRE(got[gi].count(b) == 1);
        lm_check_close(what + " [fermion vs boson, mean field]", got[gi].at(b),
                       want[gi].at(b), scales[gi], lm_rtol_gate);
      }
    }
  }

  std::vector<Correlation> got_rows;
  REQUIRE_NOTHROW(got_rows = fermion->measure_correlation());
  const auto want_rows = boson->measure_correlation();
  REQUIRE(want_rows.size() == static_cast<std::size_t>(n_unit) * 3u * 3u * 2u);
  CHECK(got_rows.size() == want_rows.size());
  std::map<std::tuple<int, int, int, int, int>, lm_cplx> got_map;
  for (const Correlation& row : got_rows) {
    got_map[std::make_tuple(row.left_index, row.right_dx, row.right_dy,
                            row.left_op, row.right_op)] =
        lm_cplx(row.real, row.imag);
  }
  const double cscale = lm_max_abs_entry(A) * lm_max_abs_entry(B);
  for (const Correlation& row : want_rows) {
    const auto key = std::make_tuple(row.left_index, row.right_dx, row.right_dy,
                                     row.left_op, row.right_op);
    const std::string what =
        label + " correlation [" + std::to_string(row.left_op) + ", " +
        std::to_string(row.right_op) + "] " +
        lm_bond_label(row.left_index, row.right_dx, row.right_dy);
    INFO(what);
    REQUIRE(got_map.count(key) == 1);
    lm_check_close(what + " [fermion vs boson, mean field]", got_map.at(key),
                   lm_cplx(row.real, row.imag), cscale, lm_rtol_gate);
  }
}

}  // namespace

TEST_CASE(
    "longrange-mf T5-4: gate 2, with every parity even the long-range "
    "windows (2,0), (1,1), (3,0), ... and the correlation function r = 1..3 "
    "equal the bosonic mean-field path") {
  lm_run_gate2<real_tensor>(1101);
  lm_run_gate2<complex_tensor>(1102);
}

// ============================================================================
// Contract item 5: the correlation function against the window measurement
// ============================================================================

namespace {

template <class tensor>
void lm_run_window_vs_correlation(const lm_corr_case& c) {
  const std::string label =
      lm_geom_label(c.g) + " " + lm_type_name<tensor>() + " item 5";
  INFO(label);
  const int n_unit = c.g.lx * c.g.ly;
  const int d = c.g.d;
  const auto phys = lm_phys(d);
  const auto onesite = lm_onesite_groups<tensor>(c.names, d, n_unit);
  // The window side: product_twosite_op(A, B) of every same-parity pair at
  // (r, 0) and (0, r) on a solver without [correlation].
  tenes::Operators<tensor> twosite;
  std::map<std::pair<int, int>, int> group_of;
  int group = 0;
  for (const auto& p : c.pairs) {
    if (!lm_pair_even(c, p.first, p.second)) {
      continue;
    }
    const lm_onesite A = lm_onesite_op(c.names[p.first], d);
    const lm_onesite B = lm_onesite_op(c.names[p.second], d);
    const tensor op4 =
        rf::product_twosite_op(lm_onesite_tensor<tensor>(A),
                               lm_onesite_tensor<tensor>(B), phys, phys, B.odd);
    group_of[p] = group;
    for (int s = 0; s < n_unit; ++s) {
      for (int r = 1; r <= c.r_max; ++r) {
        twosite.emplace_back("pair", group, s, r, 0, op4);
        twosite.emplace_back("pair", group, s, 0, r, op4);
      }
    }
    ++group;
  }
  auto ref = lm_make<tensor>(c.g, true, phys, onesite, twosite,
                             CorrelationParameter{});
  auto corr = lm_make<tensor>(c.g, true, phys, onesite, {}, lm_corparam(c));
  lm_seed(*ref, c.seed, false);
  lm_seed(*corr, c.seed, false);
  std::vector<std::map<Bond, typename tensor::value_type>> windows;
  REQUIRE_NOTHROW(windows = ref->measure_twosite());
  std::vector<Correlation> rows;
  REQUIRE_NOTHROW(rows = corr->measure_correlation());
  CHECK(rows.size() == static_cast<std::size_t>(n_unit) * c.pairs.size() *
                           static_cast<std::size_t>(c.r_max) * 2u);
  double max_signal = 0.0;
  for (const Correlation& row : rows) {
    const std::string what =
        label + " [" + c.names[row.left_op] + ", " + c.names[row.right_op] +
        "] " + lm_bond_label(row.left_index, row.right_dx, row.right_dy);
    INFO(what);
    const std::pair<int, int> p{row.left_op, row.right_op};
    if (!lm_pair_even(c, p.first, p.second)) {
      CHECK(row.real == 0.0);
      CHECK(row.imag == 0.0);
      continue;
    }
    REQUIRE(group_of.count(p) == 1);
    const int gi = group_of.at(p);
    const Bond b{row.left_index, row.right_dx, row.right_dy};
    REQUIRE(windows[gi].count(b) == 1);
    const lm_cplx want = windows[gi].at(b);
    const double scale =
        lm_max_abs_entry(
            lm_onesite_tensor<tensor>(lm_onesite_op(c.names[p.first], d))) *
        lm_max_abs_entry(
            lm_onesite_tensor<tensor>(lm_onesite_op(c.names[p.second], d)));
    lm_check_close(what + " [correlation vs window]",
                   lm_cplx(row.real, row.imag), want, scale, lm_rtol_window);
    max_signal = std::max(max_signal, std::abs(want) / scale);
  }
  INFO("largest |value| / scale " << max_signal);
  CHECK(max_signal > lm_min_signal);
}

}  // namespace

TEST_CASE(
    "longrange-mf T5-5: with meanfield_env the correlation function r = "
    "1..3 equals the window measurement of the same pair (d = 2 real, d = 4 "
    "complex)") {
  const lm_corr_case c2{{3, 2, 2, 2, 3},
                        lm_names2,
                        {{1, 2}, {0, 0}, {1, 1}, {1, 0}},
                        3,
                        1201,
                        {}};
  lm_run_window_vs_correlation<real_tensor>(c2);
  const lm_corr_case c4{{2, 2, 4, 2, 2},
                        lm_names4,
                        {{1, 2}, {3, 4}, {0, 0}, {0, 1}},
                        3,
                        1204,
                        {}};
  lm_run_window_vs_correlation<complex_tensor>(c4);
}

// ============================================================================
// Contract item 1: input acceptance and rejection
// ============================================================================

namespace {

using lm_gtensor = real_tensor;

struct lm_guard_input {
  tenes::itps::PEPS_Parameters params;
  tenes::SquareLattice lattice = lm_lattice({2, 2, 2, 2, 2});
  tenes::Operators<lm_gtensor> onesite;
  tenes::Operators<lm_gtensor> twosite;
  tenes::Operators<lm_gtensor> multisite;
  CorrelationParameter corparam;

  //! One-site groups 0 = n, 1 = c+, 2 = c on every site, mean field.
  explicit lm_guard_input(int d = 2) {
    lattice = lm_lattice({2, 2, d, 2, 2});
    params = lm_params<lm_gtensor>(lattice.N_UNIT, true, lm_phys(d), true);
    const std::vector<std::string> names =
        d == 2 ? std::vector<std::string>{"n", "cdag", "c"}
               : std::vector<std::string>{"n", "cdag_up", "c_up"};
    onesite = lm_onesite_groups<lm_gtensor>(names, d, lattice.N_UNIT);
  }

  void add_twosite(int dx, int dy, const std::string& kind = "hop") {
    const int d = lattice.physical_dims[0];
    twosite.emplace_back(kind, 0, 0, dx, dy,
                         lm_kind_tensor<lm_gtensor>(kind, d));
  }

  void add_ops(int dx, int dy, int i, int j) {
    twosite.emplace_back("ops", 0, 0, dx, dy, std::vector<int>{i, j});
  }

  void validate() const {
    tenes::itps::validate_fermion_constraints(
        params, lattice, tenes::EvolutionOperators<lm_gtensor>{},
        tenes::EvolutionOperators<lm_gtensor>{}, onesite, twosite, multisite,
        corparam);
  }

  //! The measurement-time guard of a solver built from this input.
  void validate_measurement() const {
    lm_state<lm_gtensor> state(MPI_COMM_WORLD, params, lattice,
                               tenes::EvolutionOperators<lm_gtensor>{},
                               tenes::EvolutionOperators<lm_gtensor>{}, onesite,
                               twosite, multisite, corparam,
                               tenes::itps::TransferMatrix_Parameters{});
    lm_acc::validate_fermion_ctm_measurement(state);
  }
};

std::string lm_disp_name(int dx, int dy) {
  return "(" + std::to_string(dx) + ", " + std::to_string(dy) + ")";
}

void lm_check_accepts(const lm_guard_input& in, const std::string& what) {
  INFO(what);
  try {
    in.validate();
  } catch (const std::exception& e) {
    FAIL_CHECK(what << " was rejected at load: " << std::string(e.what()));
  }
  try {
    in.validate_measurement();
  } catch (const std::exception& e) {
    FAIL_CHECK(what << " was rejected at measurement: "
                    << std::string(e.what()));
  }
}

void lm_check_rejects(const lm_guard_input& in, const std::string& what) {
  INFO(what);
  try {
    in.validate();
    FAIL_CHECK(what << " was accepted");
  } catch (const tenes::input_error&) {
  } catch (const std::exception& e) {
    FAIL_CHECK(what << " threw something other than tenes::input_error: "
                    << std::string(e.what()));
  }
}

}  // namespace

TEST_CASE(
    "longrange-mf T5-1a: with meanfield_env, long-range two-site observables "
    "inside the 4x4 window are accepted, at load and at measurement") {
  for (const lm_disp dd : {lm_disp{2, 0}, lm_disp{0, 2}, lm_disp{1, 1},
                           lm_disp{-1, 2}, lm_disp{3, 0}, lm_disp{2, -1},
                           lm_disp{3, 3}, lm_disp{-3, -3}, lm_disp{0, -3}}) {
    lm_guard_input in;
    in.add_twosite(dd.dx, dd.dy);
    lm_check_accepts(in,
                     "mean field, hopping at " + lm_disp_name(dd.dx, dd.dy));
  }
  lm_guard_input in4(4);
  in4.add_twosite(2, 1, "asym4");
  lm_check_accepts(in4, "mean field, d = 4 asym4 at (2, 1)");
}

TEST_CASE(
    "longrange-mf T5-1b: with meanfield_env, long-range ops forms of "
    "same-parity operators are accepted, at load and at measurement") {
  for (const lm_disp dd : {lm_disp{2, 0}, lm_disp{1, 1}, lm_disp{0, -3},
                           lm_disp{-2, 1}, lm_disp{3, 3}}) {
    for (const auto& ij :
         std::vector<std::pair<int, int>>{{1, 2}, {2, 1}, {0, 0}}) {
      lm_guard_input in;
      in.add_ops(dd.dx, dd.dy, ij.first, ij.second);
      lm_check_accepts(in, "mean field, ops = [" + std::to_string(ij.first) +
                               ", " + std::to_string(ij.second) + "] at " +
                               lm_disp_name(dd.dx, dd.dy));
    }
  }
}

TEST_CASE(
    "longrange-mf T5-1c: with meanfield_env, r_max > 0 is accepted for odd, "
    "even and mixed-parity pairs, at load and at measurement") {
  {
    lm_guard_input in;
    in.corparam = CorrelationParameter(3, {{1, 2}, {2, 1}, {0, 0}, {0, 2}});
    lm_check_accepts(in,
                     "mean field, r_max = 3, pairs [1,2] [2,1] [0,0] [0,2]");
  }
  {
    lm_guard_input in;
    in.corparam = CorrelationParameter(1, {{0, 0}});
    lm_check_accepts(in, "mean field, r_max = 1, pair [0,0]");
  }
  {
    lm_guard_input in(4);
    in.corparam = CorrelationParameter(5, {{1, 2}, {0, 0}});
    lm_check_accepts(in, "mean field, d = 4, r_max = 5, pairs [1,2] [0,0]");
  }
  {
    // Everything at once.
    lm_guard_input in;
    in.add_twosite(2, 1);
    in.corparam = CorrelationParameter(2, {{1, 2}});
    lm_check_accepts(in, "mean field, hopping at (2, 1) and r_max = 2");
  }
}

// Green against the stub by construction; pins that the MF acceptance does
// not extend past what the CTM environment accepts either.
TEST_CASE(
    "longrange-mf [kept] T5-1d: with meanfield_env, pairs beyond 4x4, (0, 0), "
    "ops forms of different parity and multisite observables stay "
    "rejected") {
  for (const lm_disp dd : {lm_disp{4, 0}, lm_disp{0, -4}, lm_disp{0, 0}}) {
    lm_guard_input in;
    in.add_twosite(dd.dx, dd.dy);
    lm_check_rejects(in,
                     "mean field, hopping at " + lm_disp_name(dd.dx, dd.dy));
  }
  for (const lm_disp dd : {lm_disp{1, 0}, lm_disp{2, 1}}) {
    lm_guard_input in;
    in.add_ops(dd.dx, dd.dy, 0, 1);
    lm_check_rejects(
        in, "mean field, ops = [0, 1] at " + lm_disp_name(dd.dx, dd.dy));
  }
  {
    lm_guard_input in;
    lm_gtensor op3(mptensor::Shape(2, 2, 2, 2, 2, 2));
    op3.set_value(mptensor::Index(1, 1, 1, 1, 1, 1), 1.0);
    in.multisite.emplace_back("nnn", 0, 0, std::vector<int>{1, 2},
                              std::vector<int>{0, 0}, op3);
    lm_check_rejects(in, "mean field, a multisite observable");
  }
}
