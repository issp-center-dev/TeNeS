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

// ===== T3: one full-update bond in fermion mode ============================
//
// Contract: docs/superpowers/specs/2026-09-04-fermion-full-update-contract.md
// section 3.4. Pins tenes::itps::core::Full_update_bond_fermion().
//
// Included into the test_fermion_layer TU after fermion/fold_geometry.cpp and
// fermion/full_update_env.cpp, and it reuses their fixtures: the window
// geometry, the exactly contracted patch environment (fue_make_case) and the
// judgment helpers.
//
// The contract's own section 4 says what these tests can and cannot see:
// T3-i is an EXACTNESS test, and when the target is exactly representable the
// minimum of ||psi - Theta||^2_N is zero for ANY positive definite metric.
// So nothing here is asked to detect a sign error inside N - that is T2-iv's
// job. What T3-i does pin is the graded composition of Theta, the mask on the
// source physical leg, the initial guess and the final transposes.
//
// Reference states are never taken from the routine under test: the target is
// always apply_pair_op(build_pair_state(...)), the primitives fold_geometry
// already ties to the Fock oracle.

#include "../../src/iTPS/core/full_update_fermion.hpp"

namespace {

// ---- judgment: two pair states agree up to norm and global phase ----------

template <class tensor>
double fub_frobenius(const tensor& a) {
  double s = 0.0;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    typename tensor::value_type v;
    a.get_value(a.global_index(n), v);
    s += std::abs(v) * std::abs(v);
  }
  return std::sqrt(s);
}

template <class tensor>
void fub_check_proportional(const fg_ftensor<tensor>& got,
                            const fg_ftensor<tensor>& want, double tol,
                            const std::string& label) {
  REQUIRE(got.rank() == want.rank());
  REQUIRE(got.shape() == want.shape());
  REQUIRE(got.parity == want.parity);
  const double ng = fub_frobenius(got.t);
  const double nw = fub_frobenius(want.t);
  INFO(label << " [norms]: got=" << ng << " want=" << nw);
  REQUIRE(ng > 0.0);
  REQUIRE(nw > 0.0);

  // Global phase is fixed on the largest element of the reference.
  mptensor::Index anchor;
  double best = -1.0;
  for (std::size_t n = 0; n < want.t.local_size(); ++n) {
    const mptensor::Index idx = want.t.global_index(n);
    typename tensor::value_type v;
    want.t.get_value(idx, v);
    if (std::abs(v) > best) {
      best = std::abs(v);
      anchor = idx;
    }
  }
  REQUIRE(best > 0.0);
  typename tensor::value_type wa, ga;
  want.t.get_value(anchor, wa);
  got.t.get_value(anchor, ga);
  const typename tensor::value_type phase = (ga / ng) / (wa / nw);
  INFO(label << " [phase]: " << phase << " at " << ib_index_string(anchor));
  CHECK(std::abs(std::abs(phase) - 1.0) <= tol);

  double dev = 0.0;
  mptensor::Index worst;
  for (std::size_t n = 0; n < want.t.local_size(); ++n) {
    const mptensor::Index idx = want.t.global_index(n);
    typename tensor::value_type w, g;
    want.t.get_value(idx, w);
    got.t.get_value(idx, g);
    const double d = std::abs(g / ng - phase * (w / nw));
    if (d > dev) {
      dev = d;
      worst = idx;
    }
  }
  INFO(label << ": max normalized deviation = " << dev << " (tol " << tol
             << ") at " << ib_index_string(worst));
  CHECK(dev <= tol);
}

// ---- Schmidt spectrum across the connecting bond of a pair state ----------

struct fub_schmidt {
  std::size_t total = 0;
  std::size_t even = 0;
  std::size_t odd = 0;
  double smax = 0.0;
  double smin_kept = 0.0;
};

//! Count the nonzero singular values of a rank-8 pair state cut into
//! (site A legs) x (site B legs), split by the parity of the Schmidt index.
//! This is the Schmidt rank "across the connecting bond": the number of
//! product terms the state needs, sector by sector.
template <class tensor>
fub_schmidt fub_bond_schmidt(const fg_ftensor<tensor>& pair) {
  REQUIRE(pair.rank() == 8);
  fg_ftensor<tensor> u, vt;
  std::vector<double> s;
  REQUIRE(fgf::svd(pair, mptensor::Axes(0, 1, 2, 3), mptensor::Axes(4, 5, 6, 7),
                   u, s, vt) == 0);
  REQUIRE(!s.empty());
  const fgf::parity_vector& internal = u.parity.back();
  fub_schmidt r;
  for (const double v : s) {
    r.smax = std::max(r.smax, v);
  }
  REQUIRE(r.smax > 0.0);
  r.smin_kept = r.smax;
  for (std::size_t i = 0; i < s.size(); ++i) {
    if (s[i] <= 1.0e-11 * r.smax) {
      continue;
    }
    ++r.total;
    if (internal[i]) {
      ++r.odd;
    } else {
      ++r.even;
    }
    r.smin_kept = std::min(r.smin_kept, s[i]);
  }
  return r;
}

// ---- gates ----------------------------------------------------------------

/*! Plain matrix elements of  a * 1  +  b * c^dag_{A,up} c_{B,up}.
 *
 * Both terms are parity even. Their graded SVD across (in_A, out_A) |
 * (in_B, out_B) has exactly two channels, one in each parity sector: the
 * identity is a product of two even factors, the single hop is a product of
 * two odd ones. That is the largest gate whose exact action on a
 * Schmidt-rank-1 state still fits a D = 2 bond with ledger [even, odd] -
 * the full two-spin (or two-way) hopping would need two odd channels and
 * could not be represented, which would make T3-i's "exact" claim false.
 *
 * The d = 4 elements are the s = up part of ib_hop_plain(4)'s c^dag_A c_B
 * term, in the same Jordan-Wigner representation (c_A = c x 1,
 * c_B = P x c), so the two physical dimensions share one convention.
 */
template <class tensor>
tensor fub_gate_plain(int d, double a, double b) {
  tensor op(mptensor::Shape(d, d, d, d));
  for (int i = 0; i < d; ++i) {
    for (int j = 0; j < d; ++j) {
      op.set_value(mptensor::Index(i, j, i, j),
                   a * typename tensor::value_type(1.0));
    }
  }
  const typename tensor::value_type bv = b * typename tensor::value_type(1.0);
  if (d == 2) {
    // <1 0| c^dag_A c_B |0 1>
    typename tensor::value_type v;
    op.get_value(mptensor::Index(0, 1, 1, 0), v);
    op.set_value(mptensor::Index(0, 1, 1, 0), v + bv);
    return op;
  }
  REQUIRE(d == 4);
  const int i1[4] = {0, 0, 2, 2};
  const int i2[4] = {1, 3, 1, 3};
  const int o1[4] = {1, 1, 3, 3};
  const int o2[4] = {0, 2, 0, 2};
  const double sg[4] = {1.0, 1.0, -1.0, -1.0};
  for (int k = 0; k < 4; ++k) {
    typename tensor::value_type v;
    const mptensor::Index idx(i1[k], i2[k], o1[k], o2[k]);
    op.get_value(idx, v);
    op.set_value(idx, v + sg[k] * bv);
  }
  return op;
}

//! Zero every element whose index on `axis` is not 0. Applied to the two
//! connecting-bond legs it makes the pair state a product across the bond,
//! i.e. Schmidt rank 1 (checked, not assumed, by the tests).
template <class tensor>
fg_ftensor<tensor> fub_restrict_bond(const fg_ftensor<tensor>& T, int axis) {
  fg_ftensor<tensor> r = T;
  for (std::size_t n = 0; n < r.t.local_size(); ++n) {
    const mptensor::Index idx = r.t.global_index(n);
    if (idx[axis] != 0) {
      r.t.set_value(idx, typename tensor::value_type(0.0));
    }
  }
  return r;
}

// ---- parameter sets -------------------------------------------------------

//! Contract 3.4, premise 4: the cutoffs must not throw away directions the
//! exact solution needs.
inline tenes::itps::PEPS_Parameters fub_exact_params(bool gauge_fix) {
  tenes::itps::PEPS_Parameters p;
  p.print_level = tenes::PrintLevel::none;
  p.Inverse_Env_cut = 1.0e-14;
  p.Full_Inverse_precision = 1.0e-14;
  p.Full_Convergence_Epsilon = 1.0e-14;
  p.Full_max_iteration = 200;
  p.Full_Gauge_Fix = gauge_fix;
  p.Full_Use_FastFullUpdate = false;
  p.fermion = true;
  return p;
}

// ---- drivers --------------------------------------------------------------

template <class tensor>
void fub_run_fermion(const fue_case<tensor>& c, const fg_ftensor<tensor>& Tn1,
                     const fg_ftensor<tensor>& Tn2,
                     const fg_ftensor<tensor>& gate,
                     const tenes::itps::PEPS_Parameters& params,
                     fg_ftensor<tensor>& out1, fg_ftensor<tensor>& out2) {
  tenes::itps::core::Full_update_bond_fermion(
      c.env[0], c.env[1], c.env[2], c.env[3], c.env[4], c.env[5], c.env[6],
      c.env[7], c.env[8], c.env[9], Tn1, Tn2, gate, c.dir_e, params, out1,
      out2);
}

//! The bosonic routine on the same window. The bosonic builder only knows the
//! horizontal layout, so a vertical window enters rotated by one environment
//! slot - the wiring iTPS<tensor>::full_update() uses for source_leg == 3.
template <class tensor>
void fub_run_boson(const fue_case<tensor>& c, const fg_ftensor<tensor>& Tn1,
                   const fg_ftensor<tensor>& Tn2, const tensor& op_plain,
                   const tenes::itps::PEPS_Parameters& params, tensor& out1,
                   tensor& out2) {
  std::vector<tensor> s;
  for (int i = 4; i < 10; ++i) {
    s.push_back(fue_split_edge(c.env[i], c.D));
  }
  if (c.dir == 'h') {
    tenes::itps::core::Full_update_bond(
        c.env[0], c.env[1], c.env[2], c.env[3], s[0], s[1], s[2], s[3], s[4],
        s[5], Tn1.t, Tn2.t, op_plain, 2, params, out1, out2);
  } else {
    tenes::itps::core::Full_update_bond(
        c.env[1], c.env[2], c.env[3], c.env[0], s[1], s[2], s[3], s[4], s[5],
        s[0], Tn1.t, Tn2.t, op_plain, 3, params, out1, out2);
  }
}

// ---- shared output checks (T3-iii) ---------------------------------------

template <class tensor>
void fub_check_output_parity(const fue_case<tensor>& c,
                             const fg_ftensor<tensor>& Tn1_new,
                             const fg_ftensor<tensor>& Tn2_new,
                             const fg_ftensor<tensor>& Tn1,
                             const fg_ftensor<tensor>& Tn2,
                             const std::string& label) {
  INFO(label << " [T3-iii]");
  REQUIRE(Tn1_new.rank() == 5);
  REQUIRE(Tn2_new.rank() == 5);
  {
    // Every leg but the connecting bond keeps its dimension; the bond's own
    // dimension is only required to agree at the two ends (the contract does
    // not promise it stays at D).
    const int a_bond0 = c.dir == 'h' ? 2 : 3;
    const int b_bond0 = c.dir == 'h' ? 0 : 1;
    for (int ax = 0; ax < 5; ++ax) {
      if (ax != a_bond0) {
        CHECK(Tn1_new.shape()[ax] == Tn1.shape()[ax]);
      }
      if (ax != b_bond0) {
        CHECK(Tn2_new.shape()[ax] == Tn2.shape()[ax]);
      }
    }
    REQUIRE(Tn1_new.shape()[a_bond0] == Tn2_new.shape()[b_bond0]);
  }
  const double m1 = std::max(1.0, fgf::max_abs(Tn1_new));
  const double m2 = std::max(1.0, fgf::max_abs(Tn2_new));
  REQUIRE(fgf::max_abs(Tn1_new) > 0.0);
  REQUIRE(fgf::max_abs(Tn2_new) > 0.0);
  INFO(label << " [parity violation]: " << fgf::parity_violation(Tn1_new)
             << " / " << fgf::parity_violation(Tn2_new));
  CHECK(fgf::parity_violation(Tn1_new) <= 1.0e-12 * m1);
  CHECK(fgf::parity_violation(Tn2_new) <= 1.0e-12 * m2);

  // The connecting bond's ledger must agree at both ends, and every other
  // leg must have kept the ledger it came in with.
  const int a_bond = c.dir == 'h' ? 2 : 3;
  const int b_bond = c.dir == 'h' ? 0 : 1;
  CHECK(Tn1_new.parity[a_bond] == Tn2_new.parity[b_bond]);
  for (int ax = 0; ax < 5; ++ax) {
    if (ax != a_bond) {
      CHECK(Tn1_new.parity[ax] == Tn1.parity[ax]);
    }
    if (ax != b_bond) {
      CHECK(Tn2_new.parity[ax] == Tn2.parity[ax]);
    }
  }
}

//! Premise: the metric the optimizer works in is positive DEFINITE. Without
//! that, ||psi - Theta||_N = 0 no longer forces psi == Theta and T3-i would
//! be comparing against one of many equally optimal answers.
template <class tensor>
void fub_require_definite_environment(const fue_case<tensor>& c,
                                      const std::string& label) {
  const auto N = fue_build_N(c.env, c.QA, c.QB, c.dir_e);
  REQUIRE(N.forbidden_ratio <= 1.0e-10);
  std::vector<double> w;
  REQUIRE(mptensor::eigh(N.N_plain, mptensor::Axes(0, 1), mptensor::Axes(2, 3),
                         w) == 0);
  REQUIRE(!w.empty());
  double wmax = w[0];
  double wmin = w[0];
  for (const double v : w) {
    wmax = std::max(wmax, v);
    wmin = std::min(wmin, v);
  }
  INFO(label << " [premise: environment spectrum] min=" << wmin
             << " max=" << wmax);
  REQUIRE(wmax > 0.0);
  REQUIRE(wmin > 1.0e-10 * wmax);
}

// ---- T3-i / T3-v: exactness ----------------------------------------------

template <class tensor>
void fub_run_exact_case(char dir, int d, int seed, bool gauge_fix,
                        const char* type_name) {
  const std::string tag = std::string("T3-i ") + type_name +
                          (gauge_fix ? "" : " (T3-v, no gauge fix)");
  const fue_case<tensor> c = fue_make_case<tensor>(
      dir, fue_full_geom(dir), "eo", d, d, 2, seed, false, tag);
  INFO(c.label);
  fub_require_definite_environment(c, c.label);

  const int a_bond = dir == 'h' ? 2 : 3;
  const int b_bond = dir == 'h' ? 0 : 1;
  const fg_ftensor<tensor> Tn1 = fub_restrict_bond(c.sites[c.a], a_bond);
  const fg_ftensor<tensor> Tn2 = fub_restrict_bond(c.sites[c.b], b_bond);
  REQUIRE(fgf::max_abs(Tn1) > 0.0);
  REQUIRE(fgf::max_abs(Tn2) > 0.0);
  REQUIRE(Tn1.parity[a_bond] == Tn2.parity[b_bond]);
  REQUIRE(fg_ledger_mixed(Tn1.parity[a_bond]));

  const tensor plain = fub_gate_plain<tensor>(d, 1.0, 0.8);
  const fg_ftensor<tensor> gate = fgf::wrap_twosite_gate(
      plain, Tn1.parity[4], Tn2.parity[4]);

  const fg_ftensor<tensor> psi = fgf::build_pair_state(Tn1, Tn2, c.dir_e);
  const fg_ftensor<tensor> theta = fgf::apply_pair_op(psi, gate);

  // ---- the four premises of contract 3.4 (T3-i) ----------------------
  // 1. the initial state is a product across the connecting bond
  const fub_schmidt s0 = fub_bond_schmidt(psi);
  INFO(c.label << " [premise 1] initial Schmidt: total=" << s0.total
               << " even=" << s0.even << " odd=" << s0.odd);
  REQUIRE(s0.total == 1);
  // 2./3. the exact target fits the bond, and it occupies BOTH sectors
  const fub_schmidt s1 = fub_bond_schmidt(theta);
  const std::size_t D = Tn1.shape()[a_bond];
  INFO(c.label << " [premise 2/3] target Schmidt: total=" << s1.total
               << " even=" << s1.even << " odd=" << s1.odd << " D=" << D
               << " smin/smax=" << s1.smin_kept / s1.smax);
  REQUIRE(s1.total <= D);
  REQUIRE(s1.even >= 1);
  REQUIRE(s1.odd >= 1);
  // the odd channel must not be a rounding-level contribution
  REQUIRE(s1.smin_kept > 1.0e-4 * s1.smax);
  // 4. the metric is positive definite and the cutoffs are small - the
  // former is checked above, the latter is a property of the parameters
  const tenes::itps::PEPS_Parameters params = fub_exact_params(gauge_fix);
  REQUIRE(params.Inverse_Env_cut <= 1.0e-14);
  REQUIRE(params.Full_Inverse_precision <= 1.0e-14);

  fg_ftensor<tensor> Tn1_new, Tn2_new;
  fub_run_fermion(c, Tn1, Tn2, gate, params, Tn1_new, Tn2_new);
  fub_check_output_parity(c, Tn1_new, Tn2_new, Tn1, Tn2, c.label);

  const fg_ftensor<tensor> got =
      fgf::build_pair_state(Tn1_new, Tn2_new, c.dir_e);
  fub_check_proportional(got, theta, 1.0e-8,
                         c.label + " [T3-i exact pair state]");
}

}  // namespace

TEST_CASE(
    "full update T3-i: an exactly representable gate reproduces the exact "
    "pair state") {
  fub_run_exact_case<tenes::real_tensor>('h', 2, 3100, true, "real d=2");
  fub_run_exact_case<tenes::real_tensor>('v', 2, 3200, true, "real d=2");
  fub_run_exact_case<tenes::real_tensor>('h', 4, 3300, true, "real d=4");
  fub_run_exact_case<tenes::real_tensor>('v', 4, 3400, true, "real d=4");
  fub_run_exact_case<tenes::complex_tensor>('h', 2, 3500, true, "complex d=2");
  fub_run_exact_case<tenes::complex_tensor>('v', 4, 3600, true, "complex d=4");
}

TEST_CASE("full update T3-v: exactness also holds without gauge fixing") {
  fub_run_exact_case<tenes::real_tensor>('h', 2, 3100, false, "real d=2");
  fub_run_exact_case<tenes::real_tensor>('v', 2, 3200, false, "real d=2");
  fub_run_exact_case<tenes::real_tensor>('h', 4, 3300, false, "real d=4");
  fub_run_exact_case<tenes::complex_tensor>('v', 4, 3600, false,
                                            "complex d=4");
}

// ---- T3-ii: the identity gate ---------------------------------------------

namespace {

template <class tensor>
void fub_run_identity_case(char dir, int d, int seed, const char* type_name) {
  const fue_case<tensor> c =
      fue_make_case<tensor>(dir, fue_full_geom(dir), "eo", d, d, 2, seed, false,
                            std::string("T3-ii ") + type_name);
  INFO(c.label);
  fub_require_definite_environment(c, c.label);

  const fg_ftensor<tensor> Tn1 = c.sites[c.a];
  const fg_ftensor<tensor> Tn2 = c.sites[c.b];
  const tensor plain = fub_gate_plain<tensor>(d, 1.0, 0.0);
  const fg_ftensor<tensor> gate =
      fgf::wrap_twosite_gate(plain, Tn1.parity[4], Tn2.parity[4]);

  const fg_ftensor<tensor> psi = fgf::build_pair_state(Tn1, Tn2, c.dir_e);
  const fg_ftensor<tensor> theta = fgf::apply_pair_op(psi, gate);
  // The identity really is the identity on this state.
  fub_check_proportional(theta, psi, 1.0e-13,
                         c.label + " [premise: identity gate acts trivially]");
  // Anti-vacuity: the state is genuinely entangled across the bond and uses
  // both parity sectors, so "unchanged" is a nontrivial statement.
  const int a_bond = dir == 'h' ? 2 : 3;
  const fub_schmidt s0 = fub_bond_schmidt(psi);
  INFO(c.label << " [premise] Schmidt: total=" << s0.total
               << " even=" << s0.even << " odd=" << s0.odd);
  REQUIRE(s0.total == Tn1.shape()[a_bond]);
  REQUIRE(s0.even >= 1);
  REQUIRE(s0.odd >= 1);

  const tenes::itps::PEPS_Parameters params = fub_exact_params(true);
  fg_ftensor<tensor> Tn1_new, Tn2_new;
  fub_run_fermion(c, Tn1, Tn2, gate, params, Tn1_new, Tn2_new);
  fub_check_output_parity(c, Tn1_new, Tn2_new, Tn1, Tn2, c.label);

  const fg_ftensor<tensor> got =
      fgf::build_pair_state(Tn1_new, Tn2_new, c.dir_e);
  fub_check_proportional(got, psi, 1.0e-8,
                         c.label + " [T3-ii identity leaves the state]");
}

}  // namespace

TEST_CASE("full update T3-ii: the identity gate leaves the pair state alone") {
  fub_run_identity_case<tenes::real_tensor>('h', 2, 3700, "real d=2");
  fub_run_identity_case<tenes::real_tensor>('v', 2, 3800, "real d=2");
  fub_run_identity_case<tenes::real_tensor>('h', 4, 3900, "real d=4");
  fub_run_identity_case<tenes::complex_tensor>('v', 2, 4000, "complex d=2");
}

// ---- T3-iv: bosonic degeneracy -------------------------------------------

namespace {

//! `exact` selects the Schmidt-rank-1 start plus the two-channel gate (the
//! optimum is then the same exact state for both routines); otherwise a
//! generic state and a dense gate that really has to be truncated.
template <class tensor>
void fub_run_boson_case(char dir, int d, int seed, bool exact,
                        const char* type_name) {
  const fue_case<tensor> c = fue_make_case<tensor>(
      dir, fue_full_geom(dir), "ee", d, d, 2, seed, true,
      std::string("T3-iv ") + type_name + (exact ? " exact" : " generic"));
  INFO(c.label);
  // Premise: this really is the all-even control, so the fermionic routine
  // has no sign to apply anywhere.
  for (int s = 0; s < c.patch.nsite(); ++s) {
    for (const auto& p : c.sites[s].parity) {
      for (std::size_t i = 0; i < p.size(); ++i) {
        REQUIRE(p[i] == false);
      }
    }
  }
  fub_require_definite_environment(c, c.label);

  const int a_bond = dir == 'h' ? 2 : 3;
  const int b_bond = dir == 'h' ? 0 : 1;
  fg_ftensor<tensor> Tn1 = c.sites[c.a];
  fg_ftensor<tensor> Tn2 = c.sites[c.b];
  tensor plain;
  if (exact) {
    Tn1 = fub_restrict_bond(Tn1, a_bond);
    Tn2 = fub_restrict_bond(Tn2, b_bond);
    plain = fub_gate_plain<tensor>(d, 1.0, 0.8);
  } else {
    plain = fue_random_even_plain<tensor>(Tn1.parity[4], Tn2.parity[4],
                                          seed + 17);
  }
  REQUIRE(fg_max_abs_entry(plain) > 0.0);
  // With every ledger even, wrap_twosite_gate() is the identity, so the two
  // routines really do get the same numbers.
  const fg_ftensor<tensor> gate =
      fgf::wrap_twosite_gate(plain, Tn1.parity[4], Tn2.parity[4]);
  fg_check_allclose(gate.t, plain,
                    c.label + " [premise: all-even wrap is a no-op]");

  const tenes::itps::PEPS_Parameters params = fub_exact_params(true);
  fg_ftensor<tensor> Tn1_new, Tn2_new;
  fub_run_fermion(c, Tn1, Tn2, gate, params, Tn1_new, Tn2_new);
  fub_check_output_parity(c, Tn1_new, Tn2_new, Tn1, Tn2, c.label);

  tensor b1, b2;
  fub_run_boson(c, Tn1, Tn2, plain, params, b1, b2);
  const fg_ftensor<tensor> Bn1{b1, Tn1.parity};
  const fg_ftensor<tensor> Bn2{b2, Tn2.parity};

  // Compared as pair states: the two routines are free to differ by the
  // gauge on the connecting bond.
  const fg_ftensor<tensor> got =
      fgf::build_pair_state(Tn1_new, Tn2_new, c.dir_e);
  const fg_ftensor<tensor> want = fgf::build_pair_state(Bn1, Bn2, c.dir_e);
  fub_check_proportional(got, want, 1.0e-8,
                         c.label + " [T3-iv fermion vs bosonic pair state]");
}

}  // namespace

TEST_CASE(
    "full update T3-iv: with all-even ledgers the fermionic bond update "
    "reproduces the bosonic Full_update_bond") {
  fub_run_boson_case<tenes::real_tensor>('h', 2, 4100, true, "real d=2");
  fub_run_boson_case<tenes::real_tensor>('v', 2, 4200, true, "real d=2");
  fub_run_boson_case<tenes::real_tensor>('h', 2, 4300, false, "real d=2");
  fub_run_boson_case<tenes::real_tensor>('v', 2, 4400, false, "real d=2");
  fub_run_boson_case<tenes::complex_tensor>('h', 2, 4500, true, "complex d=2");
}

// ---- T3-vi: pathological environments -------------------------------------

namespace {

/*! Ten environment tensors that close every window leg with a delta between
 *  the ket and the bra layer (chi = 1 everywhere).
 *
 *  With QA and QB isometric - which the graded QR makes them, block by
 *  block - this is the environment of the isolated two-site cluster, so
 *  N_plain is the identity: every eigenvalue is degenerate, across the
 *  parity sectors included. That is the first pathological case of the
 *  contract.
 *
 *  `deficient` additionally collapses two of site A's three environment legs
 *  onto their (ket, bra) = (0, 0) component. A's Gram matrix then has rank
 *  at most the dimension of the single remaining leg, so N_plain acquires
 *  exactly-zero eigenvalues and the Inverse_Env_cut branch is taken - the
 *  second pathological case.
 */
template <class tensor>
std::vector<tensor> fub_delta_env(std::size_t D, bool deficient) {
  std::vector<tensor> e;
  for (int i = 0; i < 4; ++i) {
    e.push_back(fue_ones<tensor>(mptensor::Shape(1, 1)));
  }
  for (int i = 0; i < 6; ++i) {
    tensor t(mptensor::Shape(1, 1, D * D));
    // eT1 (i == 0) and eT6 (i == 5) are the two legs of site A that are not
    // the connecting bond in either direction.
    const std::size_t kmax = (deficient && (i == 0 || i == 5)) ? 1 : D;
    for (std::size_t k = 0; k < kmax; ++k) {
      t.set_value(mptensor::Index(0, 0, k + D * k),
                  typename tensor::value_type(1.0));
    }
    e.push_back(t);
  }
  return e;
}

template <class tensor>
void fub_run_pathological_case(char dir, int d, int seed, bool deficient,
                               const char* type_name) {
  const std::string tag = std::string("T3-vi ") + type_name +
                          (deficient ? " rank-deficient" : " degenerate");
  fue_case<tensor> c = fue_make_case<tensor>(dir, fue_full_geom(dir), "eo", d,
                                             d, 2, seed, false, tag);
  INFO(c.label);
  c.env = fub_delta_env<tensor>(c.D, deficient);

  const auto N = fue_build_N(c.env, c.QA, c.QB, c.dir_e);
  CHECK(N.forbidden_ratio <= 1.0e-10);
  std::vector<double> w;
  REQUIRE(mptensor::eigh(N.N_plain, mptensor::Axes(0, 1), mptensor::Axes(2, 3),
                         w) == 0);
  REQUIRE(!w.empty());
  double wmax = w[0];
  double wmin = w[0];
  for (const double v : w) {
    wmax = std::max(wmax, v);
    wmin = std::min(wmin, v);
  }
  INFO(c.label << " [premise: metric spectrum] min=" << wmin
               << " max=" << wmax << " size=" << w.size());
  REQUIRE(wmax > 0.0);
  if (deficient) {
    // Premise: the metric really is singular, so Inverse_Env_cut fires.
    REQUIRE(wmin <= 1.0e-12 * wmax);
    REQUIRE(fg_ledger_mixed(c.pa));
  } else {
    // Premise: the spectrum really is degenerate across both sectors, so the
    // even/odd tie-breaking really is exercised.
    REQUIRE(fg_ledger_mixed(c.pa));
    REQUIRE(fg_ledger_mixed(c.pb));
    REQUIRE(wmax - wmin <= 1.0e-10 * wmax);
  }

  const fg_ftensor<tensor> Tn1 = c.sites[c.a];
  const fg_ftensor<tensor> Tn2 = c.sites[c.b];
  const tensor plain = fub_gate_plain<tensor>(d, 1.0, 0.8);
  const fg_ftensor<tensor> gate =
      fgf::wrap_twosite_gate(plain, Tn1.parity[4], Tn2.parity[4]);

  tenes::itps::PEPS_Parameters params = fub_exact_params(true);
  // Back to the shipped cutoffs: the rank-deficient case is only meaningful
  // when the cutoff actually removes the null directions.
  params.Inverse_Env_cut = 1.0e-12;
  params.Full_Inverse_precision = 1.0e-12;
  params.Full_Convergence_Epsilon = 1.0e-10;

  fg_ftensor<tensor> Tn1_new, Tn2_new;
  REQUIRE_NOTHROW(
      fub_run_fermion(c, Tn1, Tn2, gate, params, Tn1_new, Tn2_new));
  fub_check_output_parity(c, Tn1_new, Tn2_new, Tn1, Tn2, c.label);
  // No silent NaN: the result must still be a finite tensor.
  CHECK(std::isfinite(fgf::max_abs(Tn1_new)));
  CHECK(std::isfinite(fgf::max_abs(Tn2_new)));
}

}  // namespace

TEST_CASE(
    "full update T3-vi: degenerate and rank-deficient metrics run through "
    "without an exception and keep the parity clean") {
  for (const char dir : {'h', 'v'}) {
    fub_run_pathological_case<tenes::real_tensor>(dir, 2, 4600, false,
                                                  "real d=2");
    fub_run_pathological_case<tenes::real_tensor>(dir, 2, 4700, true,
                                                  "real d=2");
  }
  fub_run_pathological_case<tenes::complex_tensor>('h', 2, 4800, false,
                                                   "complex d=2");
  fub_run_pathological_case<tenes::complex_tensor>('v', 2, 4900, true,
                                                   "complex d=2");
}
