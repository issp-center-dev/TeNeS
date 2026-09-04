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

// ===== T2: the open-channel two-site environment N =========================
//
// Contract: docs/superpowers/specs/2026-09-04-fermion-full-update-contract.md
// sections 3.2 and 3.3. Pins tenes::fermion::build_full_update_environment().
//
// This file is included into the test_fermion_layer TU AFTER
// fermion/fold_geometry.cpp and reuses its fixtures (fg_patch, fg_make_site,
// fg_patch_state, fg_graded_norm, fg_graded_twosite, fg_folded_norm,
// fg_blob_closure, fg_check_allclose, fg_max_abs_entry, fg_ledger_mixed, ...).
//
// Environments here are never random: every one of the ten CTM slots is the
// EXACT contraction of the corresponding block of an open patch, so
//   - the folded coefficient tensor has no forbidden parity block (a generic
//     random environment would trip the guard: the parity flowing through a
//     cut of a double-layer network is conserved, and a structureless
//     environment does not respect that), and
//   - the environment is a Gram matrix in the two-site boundary index, hence
//     Hermitian and positive semidefinite, which is what T2-iii asserts.
// The wiring of the ten slots is not asserted by inspection: every fixture
// checks it by closing an identity and an operator blob in it and comparing
// against the exact plain contraction of the same patch
// (fg_folded_norm / fg_blob_closure), which fold_geometry T1/T9 already tie
// to the single-layer graded truth.
//
// Truth sources (never the machinery under test):
//   - T2-i / T2-ii / T2-iii / T2-v: the already-pinned closed path
//     contract_reduced_pair_halves_density_CTM(build_reduced_pair_halves(...)),
//     the exact patch contraction, and the bosonic
//     Create_Environment_two_sites().
//   - T2-iv: the SINGLE-LAYER graded contraction of the same patch with the
//     two window sites replaced by their pseudo-site QR factors. This shares
//     no code with the fold, and is the only test here with power over the
//     sign of N (contract section 4).

#include "../../src/fermion/full_update_env.hpp"
#include "../../src/iTPS/core/full_update.hpp"

namespace {

namespace fue = tenes::fermion;

// ---- small scalar helpers -------------------------------------------------

template <class tensor>
typename tensor::value_type fue_one() {
  return typename tensor::value_type(1.0);
}

//! Deterministic real draw, independent of ib_random_value()'s weights.
inline double fue_draw(int seed, const mptensor::Index& idx, std::size_t rank) {
  const double w[6] = {2.3, 4.7, 6.1, 8.9, 10.3, 12.7};
  double x = 0.41 * (seed + 3) + 0.017 * seed * seed;
  for (std::size_t ax = 0; ax < rank; ++ax) {
    x += w[ax] * static_cast<double>(idx[ax]);
  }
  return 0.61 * std::sin(x) + 0.29 * std::cos(1.31 * x + 0.11 * seed);
}

inline void fue_set(tenes::real_tensor& t, const mptensor::Index& idx, int seed,
                    std::size_t rank) {
  t.set_value(idx, fue_draw(seed, idx, rank));
}

inline void fue_set(tenes::complex_tensor& t, const mptensor::Index& idx,
                    int seed, std::size_t rank) {
  t.set_value(idx, std::complex<double>(fue_draw(seed, idx, rank),
                                        fue_draw(seed + 700, idx, rank)));
}

//! Tensor of all ones (used for the trivial environment slots of a window
//! that touches the patch boundary; every leg then has dimension 1).
template <class tensor>
tensor fue_ones(const mptensor::Shape& sh) {
  tensor t(sh);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    t.set_value(t.global_index(n), fue_one<tensor>());
  }
  return t;
}

//! Dense parity-even rank-4 operator over two arbitrary ledgers, legs
//! (in1, in2, out1, out2). The plain (unwrapped) matrix elements.
template <class tensor>
tensor fue_random_even_plain(const fue::parity_vector& p1,
                             const fue::parity_vector& p2, int seed) {
  tensor op(mptensor::Shape(p1.size(), p2.size(), p1.size(), p2.size()));
  for (std::size_t n = 0; n < op.local_size(); ++n) {
    const mptensor::Index idx = op.global_index(n);
    const bool odd =
        ((p1[idx[0]] != p2[idx[1]]) != (p1[idx[2]] != p2[idx[3]]));
    if (!odd) {
      fue_set(op, idx, seed, 4);
    }
  }
  return op;
}

//! Plain elementwise sum over two identically shaped tensors,
//! sum_x a(x) * b(x). NOT fermion::trace: the contract (section 3.2)
//! defines <O> through the plain product, because fermion::trace would add
//! the four-odd-leg supertrace sign.
template <class tensor>
typename tensor::value_type fue_sum_product(const tensor& a, const tensor& b) {
  REQUIRE(a.shape() == b.shape());
  typename tensor::value_type v = typename tensor::value_type(0.0);
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    const mptensor::Index idx = a.global_index(n);
    typename tensor::value_type x, y;
    a.get_value(idx, x);
    b.get_value(idx, y);
    v += x * y;
  }
  return v;
}

//! fue_sum_product() with an extra (-1)^{p(leg 0)} on every term. Used only
//! as a detection-power premise: if this does not move the value, the probe
//! operator has no support on the odd sector of N's first leg and the
//! comparison would be blind to it.
template <class tensor>
typename tensor::value_type fue_sum_product_flipped(
    const tensor& a, const tensor& b, const fue::parity_vector& p0) {
  REQUIRE(a.shape() == b.shape());
  typename tensor::value_type v = typename tensor::value_type(0.0);
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    const mptensor::Index idx = a.global_index(n);
    typename tensor::value_type x, y;
    a.get_value(idx, x);
    b.get_value(idx, y);
    v += (p0[idx[0]] ? -1.0 : 1.0) * x * y;
  }
  return v;
}

//! Relative judgment, tol = rel * max(|got|, |ref|) (fold_geometry T15's
//! convention).
template <class V>
void fue_check_rel(const std::string& label, V got, V ref, double rel) {
  const double tol = rel * std::max(std::abs(got), std::abs(ref));
  INFO(label << ": got=" << got << " ref=" << ref
             << " |diff|=" << std::abs(got - ref) << " tol=" << tol);
  CHECK(std::abs(got - ref) <= tol);
}

// ---- window geometry inside an open patch ---------------------------------

// Slot order below is the ARGUMENT order of the closure calls:
//   0..3   C1, C2, C3, C4
//   4..9   eT1, eT2, eT3, eT4, eT5, eT6
// The offsets are relative to the raster-earlier window site A. `perm`
// reorders a site's reduced legs (l, t, r, b) into the environment tensor's
// own leg order; the trailing 4 - keep legs are patch-boundary legs of
// dimension 1 and get dropped.
//
// The orders are read off the argument wiring of
// contract_reduced_pair_{horizontal,vertical}_density_CTM():
//   horizontal   C1(->eT6, ->eT1)  C2(->eT2, ->eT3)  C3(->eT3, ->eT4)
//                C4(->eT5, ->eT6)
//                eT1(C1, eT2, T_A) eT2(eT1, C2, T_B) eT3(C2, C3, R_B)
//                eT4(C3, eT5, B_B) eT5(eT4, C4, B_A) eT6(C4, C1, L_A)
//   vertical     C1(->eT6, ->eT1)  C2(->eT1, ->eT2)  C3(->eT3, ->eT4)
//                C4(->eT4, ->eT5)
//                eT1(C1, C2, T_A)  eT2(C2, eT3, R_A) eT3(eT2, C3, R_B)
//                eT4(C3, C4, B_B)  eT5(C4, eT6, L_B) eT6(eT5, C1, L_A)
struct fue_slot {
  int dx;
  int dy;
  int perm[4];
  int keep;
};

const fue_slot fue_slots_h[10] = {
    {-1, -1, {3, 2, 0, 1}, 2},  // C1  = (b, r)
    {+2, -1, {0, 3, 1, 2}, 2},  // C2  = (l, b)
    {+2, +1, {1, 0, 2, 3}, 2},  // C3  = (t, l)
    {-1, +1, {2, 1, 0, 3}, 2},  // C4  = (r, t)
    {0, -1, {0, 2, 3, 1}, 3},   // eT1 = (l, r, b)
    {+1, -1, {0, 2, 3, 1}, 3},  // eT2 = (l, r, b)
    {+2, 0, {1, 3, 0, 2}, 3},   // eT3 = (t, b, l)
    {+1, +1, {2, 0, 1, 3}, 3},  // eT4 = (r, l, t)
    {0, +1, {2, 0, 1, 3}, 3},   // eT5 = (r, l, t)
    {-1, 0, {3, 1, 2, 0}, 3},   // eT6 = (b, t, r)
};

const fue_slot fue_slots_v[10] = {
    {-1, -1, {3, 2, 0, 1}, 2},  // C1  = (b, r)
    {+1, -1, {0, 3, 1, 2}, 2},  // C2  = (l, b)
    {+1, +2, {1, 0, 2, 3}, 2},  // C3  = (t, l)
    {-1, +2, {2, 1, 0, 3}, 2},  // C4  = (r, t)
    {0, -1, {0, 2, 3, 1}, 3},   // eT1 = (l, r, b)
    {+1, 0, {1, 3, 0, 2}, 3},   // eT2 = (t, b, l)
    {+1, +1, {1, 3, 0, 2}, 3},  // eT3 = (t, b, l)
    {0, +2, {2, 0, 1, 3}, 3},   // eT4 = (r, l, t)
    {-1, +1, {3, 1, 2, 0}, 3},  // eT5 = (b, t, r)
    {-1, 0, {3, 1, 2, 0}, 3},   // eT6 = (b, t, r)
};

//! One environment slot from a site's rank-4 reduced tensor.
template <class tensor>
tensor fue_block(const tensor& R, const mptensor::Axes& perm, int keep) {
  const tensor t = mptensor::transpose(R, perm);
  mptensor::Shape sh;
  for (int ax = 0; ax < keep; ++ax) {
    sh.push(t.shape()[ax]);
  }
  for (std::size_t ax = static_cast<std::size_t>(keep); ax < t.shape().size();
       ++ax) {
    REQUIRE(t.shape()[ax] == 1);
  }
  return mptensor::reshape(t, sh);
}

//! The ten exact environment tensors of the window whose first site is
//! (ax, ay). Slots outside the patch become all-ones tensors of dimension 1;
//! the corresponding window leg is then a patch-boundary leg of dimension 1
//! too, so the closure is still the exact patch contraction.
template <class tensor>
std::vector<tensor> fue_build_env(const fg_patch& p,
                                  const std::vector<tensor>& reduced, int ax,
                                  int ay, char dir) {
  const fue_slot* slots = (dir == 'h') ? fue_slots_h : fue_slots_v;
  std::vector<tensor> e;
  for (int i = 0; i < 10; ++i) {
    const int x = ax + slots[i].dx;
    const int y = ay + slots[i].dy;
    if (x < 0 || x >= p.lx || y < 0 || y >= p.ly) {
      mptensor::Shape sh;
      for (int k = 0; k < slots[i].keep; ++k) {
        sh.push(1);
      }
      e.push_back(fue_ones<tensor>(sh));
      continue;
    }
    mptensor::Axes perm;
    for (int k = 0; k < 4; ++k) {
      perm.push(slots[i].perm[k]);
    }
    e.push_back(fue_block(reduced[p.site(x, y)], perm, slots[i].keep));
  }
  return e;
}

template <class tensor>
typename tensor::value_type fue_close_blob(const std::vector<tensor>& e,
                                           char dir, const tensor& blob) {
  if (dir == 'h') {
    return fgf::contract_reduced_pair_horizontal_density_CTM(
        e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9], blob);
  }
  return fgf::contract_reduced_pair_vertical_density_CTM(
      e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9], blob);
}

template <class tensor>
typename tensor::value_type fue_close_halves(
    const std::vector<tensor>& e, const fgf::reduced_pair_halves<tensor>& h) {
  return fgf::contract_reduced_pair_halves_density_CTM(
      e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9], h);
}

template <class tensor>
fue::full_update_environment<tensor> fue_build_N(
    const std::vector<tensor>& e, const fg_ftensor<tensor>& QA,
    const fg_ftensor<tensor>& QB, fgf::reduced_pair_direction dir) {
  return fue::build_full_update_environment(e[0], e[1], e[2], e[3], e[4], e[5],
                                            e[6], e[7], e[8], e[9], QA, QB,
                                            dir);
}

// ---- the fixture ----------------------------------------------------------

//! Physical-leg ledger: the fermionic one from ib_phys_parity(), or the
//! all-even control used by T2-v.
inline fue::parity_vector fue_phys(int d, bool all_even) {
  if (all_even) {
    return fue::parity_vector(static_cast<std::size_t>(d), false);
  }
  return ib_phys_parity(d);
}

template <class tensor>
struct fue_case {
  fg_patch patch{0, 0};
  int a = 0;
  int b = 0;
  char dir = 'h';
  fgf::reduced_pair_direction dir_e = fgf::reduced_pair_direction::horizontal;
  std::size_t D = 0;
  std::vector<fg_ftensor<tensor>> sites;
  std::vector<tensor> reduced;
  std::vector<tensor> env;
  fg_ftensor<tensor> QA, QB;    // rank 4, environment-side QR factors
  fg_ftensor<tensor> QAp, QBp;  // rank 5 pseudo-sites (internal leg = phys)
  fue::parity_vector pa, pb;    // ledgers of the two internal legs
  std::string label;
};

//! Insert a dimension-1 even leg at position `slot` of a rank-4 factor,
//! producing the rank-5 pseudo-site of section 3.2 (T2-i). mptensor flattens
//! column-major, so inserting a dimension-1 axis is a pure reshape: the
//! element order is unchanged.
template <class tensor>
fg_ftensor<tensor> fue_pseudo_site(const fg_ftensor<tensor>& Q, int slot) {
  const mptensor::Shape qs = Q.shape();
  REQUIRE(qs.size() == 4);
  mptensor::Shape sh;
  fue::leg_parities lp;
  for (int ax = 0; ax < 5; ++ax) {
    if (ax == slot) {
      sh.push(1);
      lp.push_back(fue::parity_vector{false});
    } else {
      const int src = ax < slot ? ax : ax - 1;
      sh.push(qs[src]);
      lp.push_back(Q.parity[src]);
    }
  }
  return fg_ftensor<tensor>{mptensor::reshape(Q.t, sh), lp};
}

//! Patch size and the position of the raster-earlier window site inside it.
struct fue_geom {
  int lx;
  int ly;
  int ax;
  int ay;
};

//! The geometry in which EVERY one of the six window legs is a genuine
//! interior bond and no environment slot is trivial: 4 columns x 3 rows for
//! a horizontal window, 3 columns x 4 rows for a vertical one.
inline fue_geom fue_full_geom(char dir) {
  return dir == 'h' ? fue_geom{4, 3, 1, 1} : fue_geom{3, 4, 1, 1};
}

/*! Build one window fixture.
 *
 * The window sites are (ax, ay) and its right (horizontal) or lower
 * (vertical) neighbour. Environment slots that fall outside the patch become
 * trivial; the corresponding window leg is then a patch-boundary leg of
 * dimension 1, so the closure is still the exact patch contraction. Use
 * fue_full_geom() when the test needs every slot to carry weight.
 */
template <class tensor>
fue_case<tensor> fue_make_case(char dir, const fue_geom& g,
                               const std::string& vps, int dA, int dB,
                               int env_d, int seed, bool all_even,
                               const std::string& tag) {
  fue_case<tensor> c;
  c.dir = dir;
  c.dir_e = dir == 'h' ? fgf::reduced_pair_direction::horizontal
                       : fgf::reduced_pair_direction::vertical;
  c.patch = fg_patch{g.lx, g.ly};
  const int ax = g.ax;
  const int ay = g.ay;
  c.a = c.patch.site(ax, ay);
  c.b = dir == 'h' ? c.patch.site(ax + 1, ay) : c.patch.site(ax, ay + 1);
  c.label = tag + " dir=" + dir + " patch=" + std::to_string(g.lx) + "x" +
            std::to_string(g.ly) + " at(" + std::to_string(ax) + "," +
            std::to_string(ay) + ")" + " vp=" + vps +
            " dA=" + std::to_string(dA) + " dB=" + std::to_string(dB) +
            " env_d=" + std::to_string(env_d) +
            " seed=" + std::to_string(seed) + (all_even ? " all-even" : "");
  INFO(c.label);

  const fue::parity_vector vp =
      all_even ? fue::parity_vector(vps.size(), false) : fg_pv(vps);
  c.D = vp.size();
  const fue::parity_vector phys_env = fue_phys(env_d, all_even);
  const fue::parity_vector phys_a = fue_phys(dA, all_even);
  const fue::parity_vector phys_b = fue_phys(dB, all_even);

  for (int s = 0; s < c.patch.nsite(); ++s) {
    const fue::parity_vector& phys =
        (s == c.a) ? phys_a : ((s == c.b) ? phys_b : phys_env);
    c.sites.push_back(fg_make_site<tensor>(c.patch, s, vp, phys, seed + s));
  }
  for (int s = 0; s < c.patch.nsite(); ++s) {
    REQUIRE(fgf::parity_violation(c.sites[s]) == 0.0);
    c.reduced.push_back(fgf::build_reduced(c.sites[s]));
  }
  c.env = fue_build_env(c.patch, c.reduced, ax, ay, dir);

  // Premise: the ten slots are wired the way the closures expect. Judged
  // against the exact plain contraction of the same patch, which
  // fold_geometry T1/T9 already tie to the single-layer graded truth.
  const tensor id_blob =
      fgf::build_reduced_identity_pair(c.sites[c.a], c.sites[c.b], c.dir_e);
  const auto norm_patch = fg_folded_norm(c.patch, c.reduced);
  INFO(c.label << " [premise: exact patch norm = " << norm_patch << "]");
  REQUIRE(std::abs(norm_patch) > 0.0);
  fue_check_rel(c.label + " [premise: env wiring, identity blob]",
                fue_close_blob(c.env, dir, id_blob), norm_patch, 1.0e-12);
  {
    // ... and once more with a dense operator blob, so a slot swap that the
    // norm happens to be blind to is caught too.
    const tensor plain_op =
        fue_random_even_plain<tensor>(phys_a, phys_b, seed + 41);
    const fg_ftensor<tensor> gate =
        fgf::wrap_twosite_gate(plain_op, phys_a, phys_b);
    const tensor op_blob = fgf::build_reduced_pair_direct(
        c.sites[c.a], c.sites[c.b], gate, c.dir_e);
    const auto ref =
        fg_blob_closure(c.patch, c.reduced, c.a, c.b, dir, op_blob);
    REQUIRE(std::abs(ref) > 0.0);
    fue_check_rel(c.label + " [premise: env wiring, operator blob]",
                  fue_close_blob(c.env, dir, op_blob), ref, 1.0e-12);
  }

  // Environment-side QR factors, leg orders exactly as tabulated in
  // contract section 2.1.
  fg_ftensor<tensor> RA, RB;
  if (dir == 'h') {
    REQUIRE(fgf::qr(c.sites[c.a], mptensor::Axes(0, 1, 3), mptensor::Axes(2, 4),
                    c.QA, RA) == 0);
    REQUIRE(fgf::qr(c.sites[c.b], mptensor::Axes(1, 2, 3), mptensor::Axes(0, 4),
                    c.QB, RB) == 0);
    c.QAp = fue_pseudo_site(c.QA, 2);  // (l, t, [dummy r], b, a)
    c.QBp = fue_pseudo_site(c.QB, 0);  // ([dummy l], t, r, b, beta)
  } else {
    REQUIRE(fgf::qr(c.sites[c.a], mptensor::Axes(0, 1, 2), mptensor::Axes(3, 4),
                    c.QA, RA) == 0);
    REQUIRE(fgf::qr(c.sites[c.b], mptensor::Axes(0, 2, 3), mptensor::Axes(1, 4),
                    c.QB, RB) == 0);
    c.QAp = fue_pseudo_site(c.QA, 3);  // (l, t, r, [dummy b], a)
    c.QBp = fue_pseudo_site(c.QB, 1);  // (l, [dummy t], r, b, beta)
  }
  c.pa = c.QA.parity[3];
  c.pb = c.QB.parity[3];
  REQUIRE(fgf::parity_violation(c.QAp) <=
          1.0e-13 * std::max(1.0, fgf::max_abs(c.QAp)));
  REQUIRE(fgf::parity_violation(c.QBp) <=
          1.0e-13 * std::max(1.0, fgf::max_abs(c.QBp)));
  REQUIRE(fg_max_abs_entry(c.QA.t) > 0.0);
  REQUIRE(fg_max_abs_entry(c.QB.t) > 0.0);
  if (!all_even) {
    // Premise: the internal legs really carry both parities, otherwise no
    // odd channel exists and every sign check below is vacuous.
    REQUIRE(fg_ledger_mixed(c.pa));
    REQUIRE(fg_ledger_mixed(c.pb));
  }
  return c;
}

// ---- T2-i: the open coefficient tensor reproduces the closed value --------

//! The identity on the internal-leg pair, as a plain wrap-form operator.
//! Its closed value is the natural magnitude of this window and serves as
//! the signal floor for the random probes - a scale that does not come from
//! the tensor under test.
template <class tensor>
tensor fue_internal_identity(const fue::parity_vector& p1,
                             const fue::parity_vector& p2) {
  tensor op(mptensor::Shape(p1.size(), p2.size(), p1.size(), p2.size()));
  for (std::size_t i = 0; i < p1.size(); ++i) {
    for (std::size_t j = 0; j < p2.size(); ++j) {
      op.set_value(mptensor::Index(i, j, i, j), fue_one<tensor>());
    }
  }
  return op;
}

//! One probe operator: sum_x N(x) O(x) against the closed fold path.
template <class tensor>
void fue_run_probe(const fue_case<tensor>& c,
                   const fue::full_update_environment<tensor>& N, int op_seed,
                   double ref_scale, const std::string& what) {
  const std::string label = c.label + " " + what +
                            " op_seed=" + std::to_string(op_seed);
  INFO(label);
  const tensor plain_op = fue_random_even_plain<tensor>(c.pa, c.pb, op_seed);
  REQUIRE(fg_max_abs_entry(plain_op) > 0.0);
  const fg_ftensor<tensor> O = fgf::wrap_twosite_gate(plain_op, c.pa, c.pb);

  const auto halves =
      fgf::build_reduced_pair_halves(c.QAp, c.QBp, O, c.dir_e);
  const auto ref = fue_close_halves(c.env, halves);
  const auto got = fue_sum_product(N.N.t, O.t);

  // Anti-vacuity: the probe must produce signal. The floor is built from the
  // identity closure, i.e. from the reference path only.
  const double floor = 1.0e-6 * ref_scale * fg_max_abs_entry(O.t);
  INFO(label << " [signal]: ref=" << ref << " floor=" << floor);
  REQUIRE(std::abs(ref) > floor);

  fue_check_rel(label + " [T2-i open vs closed]", got, ref, 1.0e-12);

  // ... and it must reach the odd sector of N's first leg, otherwise a sign
  // error confined to that block would not show up in this scalar.
  const auto flipped = fue_sum_product_flipped(N.N.t, O.t, c.pa);
  INFO(label << " [odd-sector reach]: plain=" << got
             << " odd-flipped=" << flipped);
  CHECK(std::abs(flipped - got) >
        1.0e-6 * (std::abs(got) + std::abs(flipped)));
}

template <class tensor>
void fue_run_t2i_case(char dir, const std::string& vps, int dA, int dB,
                      int seed, const char* type_name) {
  const fue_case<tensor> c =
      fue_make_case<tensor>(dir, fue_full_geom(dir), vps, dA, dB, 2, seed,
                            false, std::string("T2-i ") + type_name);
  INFO(c.label);
  const auto N = fue_build_N(c.env, c.QA, c.QB, c.dir_e);

  // Shape / ledger contract of N: legs (in_A, in_B, out_A, out_B).
  const std::size_t nA = c.pa.size();
  const std::size_t nB = c.pb.size();
  REQUIRE(N.N.rank() == 4);
  REQUIRE(N.N.shape()[0] == nA);
  REQUIRE(N.N.shape()[1] == nB);
  REQUIRE(N.N.shape()[2] == nA);
  REQUIRE(N.N.shape()[3] == nB);
  CHECK(N.N.parity[0] == c.pa);
  CHECK(N.N.parity[1] == c.pb);
  CHECK(N.N.parity[2] == c.pa);
  CHECK(N.N.parity[3] == c.pb);
  REQUIRE(N.N_plain.shape() == N.N.shape());
  REQUIRE(fg_max_abs_entry(N.N.t) > 0.0);

  // T2-ii, clean half: the forbidden block was empty to begin with.
  INFO(c.label << " [T2-ii forbidden_ratio = " << N.forbidden_ratio << "]");
  CHECK(N.forbidden_ratio <= 1.0e-10);
  CHECK(fgf::parity_violation(N.N) <=
        1.0e-12 * std::max(1.0, fgf::max_abs(N.N)));

  // N_plain is N with the elementwise (-1)^{p(in_A) p(in_B)} mask.
  {
    tensor want = N.N.t;
    for (std::size_t n = 0; n < want.local_size(); ++n) {
      const mptensor::Index idx = want.global_index(n);
      if (c.pa[idx[0]] && c.pb[idx[1]]) {
        typename tensor::value_type v;
        want.get_value(idx, v);
        want.set_value(idx, -v);
      }
    }
    fg_check_allclose(N.N_plain, want, c.label + " [N_plain = masked N]");
  }

  // The magnitude of this window, taken from the closed path alone.
  double ref_scale = 0.0;
  {
    const tensor id_plain = fue_internal_identity<tensor>(c.pa, c.pb);
    const fg_ftensor<tensor> Oid =
        fgf::wrap_twosite_gate(id_plain, c.pa, c.pb);
    const auto halves =
        fgf::build_reduced_pair_halves(c.QAp, c.QBp, Oid, c.dir_e);
    const auto ref_id = fue_close_halves(c.env, halves);
    INFO(c.label << " [identity closure = " << ref_id << "]");
    REQUIRE(std::abs(ref_id) > 0.0);
    ref_scale = std::abs(ref_id);
    // The identity probe is itself a T2-i row.
    fue_check_rel(c.label + " [T2-i open vs closed, identity probe]",
                  fue_sum_product(N.N.t, Oid.t), ref_id, 1.0e-12);
  }

  for (const int op_seed : {13, 57, 91}) {
    fue_run_probe(c, N, op_seed, ref_scale, "probe");
  }
}

}  // namespace

TEST_CASE(
    "full update T2-i: sum_x N(x) O(x) equals the closed fold value for every "
    "parity-even wrap-form operator") {
  // Contract section 3.2 coverage: both directions, real and complex, d = 2
  // and d = 4, a mixed virtual ledger, and at least one nA != nB shape (the
  // last two rows: the two window sites carry different physical dimensions,
  // so their QR internal legs differ in size).
  fue_run_t2i_case<tenes::real_tensor>('h', "eo", 2, 2, 100, "real");
  fue_run_t2i_case<tenes::real_tensor>('v', "eo", 2, 2, 200, "real");
  fue_run_t2i_case<tenes::real_tensor>('h', "eo", 4, 4, 300, "real");
  fue_run_t2i_case<tenes::real_tensor>('v', "eo", 4, 4, 400, "real");
  fue_run_t2i_case<tenes::real_tensor>('h', "eeo", 2, 2, 500, "real");
  fue_run_t2i_case<tenes::complex_tensor>('h', "eo", 2, 2, 600, "complex");
  fue_run_t2i_case<tenes::complex_tensor>('v', "eo", 4, 4, 700, "complex");
  // nA != nB
  fue_run_t2i_case<tenes::real_tensor>('h', "eo", 2, 4, 800, "real");
  fue_run_t2i_case<tenes::complex_tensor>('v', "eo", 4, 2, 900, "complex");
}

namespace {

//! Elementwise reconstruction of N from the closed path: for every
//! parity-allowed index quadruple, close the fold with the delta operator
//! supported on that single element and compare with N's element. This is
//! T2-i without the averaging a dense probe does.
template <class tensor>
void fue_run_t2i_elementwise(char dir, const std::string& vps, int dA, int dB,
                             int seed, const char* type_name) {
  const fue_case<tensor> c = fue_make_case<tensor>(
      dir, fue_full_geom(dir), vps, dA, dB, 2, seed, false,
      std::string("T2-i-elementwise ") + type_name);
  INFO(c.label);
  const auto N = fue_build_N(c.env, c.QA, c.QB, c.dir_e);
  const std::size_t nA = c.pa.size();
  const std::size_t nB = c.pb.size();
  REQUIRE(nA * nB <= 32u);  // keep the number of closures bounded

  double worst = 0.0;
  double scale = 0.0;
  std::size_t probed = 0;
  std::size_t nonzero = 0;
  for (std::size_t i1 = 0; i1 < nA; ++i1) {
    for (std::size_t i2 = 0; i2 < nB; ++i2) {
      for (std::size_t o1 = 0; o1 < nA; ++o1) {
        for (std::size_t o2 = 0; o2 < nB; ++o2) {
          const bool odd = ((c.pa[i1] != c.pb[i2]) != (c.pa[o1] != c.pb[o2]));
          if (odd) {
            continue;
          }
          const mptensor::Index idx(i1, i2, o1, o2);
          tensor plain_op(mptensor::Shape(nA, nB, nA, nB));
          plain_op.set_value(idx, fue_one<tensor>());
          const fg_ftensor<tensor> O =
              fgf::wrap_twosite_gate(plain_op, c.pa, c.pb);
          const auto halves =
              fgf::build_reduced_pair_halves(c.QAp, c.QBp, O, c.dir_e);
          const auto ref = fue_close_halves(c.env, halves);
          const auto got = fue_sum_product(N.N.t, O.t);
          ++probed;
          if (std::abs(ref) > 0.0) {
            ++nonzero;
          }
          worst = std::max(worst, std::abs(got - ref));
          scale = std::max(scale, std::max(std::abs(got), std::abs(ref)));
        }
      }
    }
  }
  INFO(c.label << " [T2-i elementwise]: probed=" << probed
               << " nonzero=" << nonzero << " worst=" << worst
               << " scale=" << scale);
  REQUIRE(probed > 0);
  // Anti-vacuity: most allowed elements must actually carry signal.
  REQUIRE(nonzero * 2 >= probed);
  CHECK(worst <= 1.0e-12 * scale);
}

}  // namespace

TEST_CASE(
    "full update T2-i (elementwise): every allowed element of N equals the "
    "closed fold value of the corresponding delta operator") {
  // Two rows only: each one is nA^2 nB^2 / 2 = 128 separate closures, which
  // is the whole point (no averaging) but also the cost.
  fue_run_t2i_elementwise<tenes::real_tensor>('h', "eo", 2, 2, 1100, "real");
  fue_run_t2i_elementwise<tenes::complex_tensor>('v', "eo", 2, 2, 1200,
                                                 "complex");
}

// ---- T2-ii: the forbidden parity block guard ------------------------------

namespace {

//! Copy of the environment whose eT[slot] has an extra component on the
//! parity-ODD values of its fused (ket, bra) leg. A physical environment
//! never has one: the fused parity flowing through any cut of the double
//! layer is conserved. With it, the folded coefficient tensor grows a
//! forbidden block and build_full_update_environment() must refuse.
template <class tensor>
std::vector<tensor> fue_contaminate_env(const std::vector<tensor>& e, int slot,
                                        const fue::parity_vector& vp,
                                        int seed) {
  std::vector<tensor> bad = e;
  const fue::parity_vector pdd = fue::fuse(vp, vp);
  tensor& t = bad[slot];
  REQUIRE(t.shape().size() == 3);
  REQUIRE(t.shape()[2] == pdd.size());
  const double amp = fg_max_abs_entry(t);
  REQUIRE(amp > 0.0);
  std::size_t touched = 0;
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    if (!pdd[idx[2]]) {
      continue;
    }
    typename tensor::value_type v;
    t.get_value(idx, v);
    t.set_value(idx, v + amp * fue_one<tensor>() *
                             (0.5 + 0.25 * fue_draw(seed, idx, 3)));
    ++touched;
  }
  REQUIRE(touched > 0);
  return bad;
}

}  // namespace

TEST_CASE(
    "full update T2-ii: build_full_update_environment refuses an environment "
    "that carries a parity-odd component") {
  for (const char dir : {'h', 'v'}) {
    const fue_case<tenes::real_tensor> c = fue_make_case<tenes::real_tensor>(
        dir, fue_full_geom(dir), "eo", 2, 2, 2, 1400, false, "T2-ii real");
    INFO(c.label);
    // Control: the clean environment is accepted and its forbidden block is
    // empty. Without this, a function that always throws would pass.
    const auto clean = fue_build_N(c.env, c.QA, c.QB, c.dir_e);
    INFO(c.label << " [clean forbidden_ratio = " << clean.forbidden_ratio
                 << "]");
    CHECK(clean.forbidden_ratio <= 1.0e-10);

    // eT1 (slot 4) and eT3 (slot 6) sit on different sides of the window.
    const tenes::real_tensor id_blob = fgf::build_reduced_identity_pair(
        c.sites[c.a], c.sites[c.b], c.dir_e);
    const double clean_value = fue_close_blob(c.env, dir, id_blob);
    REQUIRE(std::abs(clean_value) > 0.0);
    for (const int slot : {4, 6}) {
      const std::vector<tenes::real_tensor> bad =
          fue_contaminate_env(c.env, slot, fg_pv("eo"), 77 + slot);
      INFO(c.label << " [contaminated slot " << slot << "]");
      // Anti-vacuity: the contamination must actually reach the network -
      // a perturbation that the closure annihilates would not create a
      // forbidden block either, and requiring a throw would be unfair.
      const double bad_value = fue_close_blob(bad, dir, id_blob);
      INFO(c.label << " [contaminated slot " << slot
                   << " closure]: clean=" << clean_value
                   << " bad=" << bad_value);
      REQUIRE(std::abs(bad_value - clean_value) >
              1.0e-6 * std::abs(clean_value));
      CHECK_THROWS_AS(fue_build_N(bad, c.QA, c.QB, c.dir_e),
                      std::runtime_error);
    }
  }
}

// ---- T2-iii: N_plain is Hermitian and positive semidefinite ---------------

namespace {

template <class tensor>
void fue_run_t2iii_case(char dir, const std::string& vps, int dA, int dB,
                        int seed, const char* type_name) {
  const fue_case<tensor> c =
      fue_make_case<tensor>(dir, fue_full_geom(dir), vps, dA, dB, 2, seed,
                            false, std::string("T2-iii ") + type_name);
  INFO(c.label);
  const auto N = fue_build_N(c.env, c.QA, c.QB, c.dir_e);
  const tensor& M = N.N_plain;

  // Hermiticity: M(a, b, a', b') == conj(M(a', b', a, b)).
  const tensor adj = mptensor::conj(mptensor::transpose(
      M, mptensor::Axes(2, 3, 0, 1)));
  fg_check_allclose(M, adj, c.label + " [T2-iii N_plain hermitian]");

  // Positive semidefiniteness of the (in_A, in_B) x (out_A, out_B) matrix.
  std::vector<double> w;
  const int info =
      mptensor::eigh(M, mptensor::Axes(0, 1), mptensor::Axes(2, 3), w);
  REQUIRE(info == 0);
  REQUIRE(!w.empty());
  double wmax = 0.0;
  double wmin = w[0];
  for (const double v : w) {
    wmax = std::max(wmax, v);
    wmin = std::min(wmin, v);
  }
  INFO(c.label << " [T2-iii spectrum]: min=" << wmin << " max=" << wmax);
  // Anti-vacuity: a zero metric would be trivially PSD.
  REQUIRE(wmax > 0.0);
  CHECK(wmin >= -1.0e-12 * wmax);
}

}  // namespace

TEST_CASE(
    "full update T2-iii: N_plain is hermitian and positive semidefinite in "
    "an exactly contracted environment") {
  fue_run_t2iii_case<tenes::real_tensor>('h', "eo", 2, 2, 1500, "real");
  fue_run_t2iii_case<tenes::real_tensor>('v', "eo", 2, 2, 1600, "real");
  fue_run_t2iii_case<tenes::real_tensor>('h', "eo", 4, 4, 1700, "real");
  fue_run_t2iii_case<tenes::complex_tensor>('h', "eo", 2, 2, 1800, "complex");
  fue_run_t2iii_case<tenes::complex_tensor>('v', "eo", 2, 4, 1900, "complex");
}

// ---- T2-iv: independent oracle (single-layer graded truth) ----------------

namespace {

/*! sum_x N(x) O(x) against the SINGLE-LAYER graded contraction of the same
 *  patch with the two window sites replaced by the pseudo-sites.
 *
 *  The pseudo-sites carry a dimension-1 even connecting bond, so contracting
 *  them in the patch is the same network the environment was built for; the
 *  operator O acts on their internal legs, which are exactly N's legs. This
 *  path uses fermion::tensordot / transpose / conj only - none of the folded
 *  machinery - so it is the one test here with power over the sign of N
 *  (contract section 4).
 */
template <class tensor>
void fue_run_t2iv_case(char dir, const fue_geom& g, const std::string& vps,
                       int dA, int dB, int seed, int op_seed,
                       const char* type_name) {
  const fue_case<tensor> c =
      fue_make_case<tensor>(dir, g, vps, dA, dB, 2, seed, false,
                            std::string("T2-iv ") + type_name);
  INFO(c.label);
  const auto N = fue_build_N(c.env, c.QA, c.QB, c.dir_e);
  CHECK(N.forbidden_ratio <= 1.0e-10);

  // The pseudo-patch: the same sites, except that the two window sites are
  // the QR factors carrying their internal leg in the physical slot.
  std::vector<fg_ftensor<tensor>> sites = c.sites;
  sites[c.a] = c.QAp;
  sites[c.b] = c.QBp;
  std::vector<int> labels;
  const fg_ftensor<tensor> psi = fg_patch_state(c.patch, sites, labels);
  const auto norm_psi = fg_graded_norm(psi);
  INFO(c.label << " [pseudo-patch graded norm = " << norm_psi << "]");
  REQUIRE(std::abs(norm_psi) > 0.0);

  const tensor plain_op = fue_random_even_plain<tensor>(c.pa, c.pb, op_seed);
  const fg_ftensor<tensor> O = fgf::wrap_twosite_gate(plain_op, c.pa, c.pb);
  const auto ref = fg_graded_twosite(psi, labels, c.a, c.b, O);
  const auto got = fue_sum_product(N.N.t, O.t);

  // Anti-vacuity: the signal floor comes from the truth side only (the
  // graded norm of the same pseudo-patch), never from N.
  const double floor =
      1.0e-6 * std::abs(norm_psi) * fg_max_abs_entry(O.t);
  INFO(c.label << " [signal]: ref=" << ref << " floor=" << floor);
  REQUIRE(std::abs(ref) > floor);

  fue_check_rel(c.label + " [T2-iv open N vs graded single-layer truth]", got,
                ref, 1.0e-12);

  // The probe must reach the odd sector of N's first leg.
  const auto flipped = fue_sum_product_flipped(N.N.t, O.t, c.pa);
  INFO(c.label << " [odd-sector reach]: plain=" << got
               << " odd-flipped=" << flipped);
  CHECK(std::abs(flipped - got) >
        1.0e-6 * (std::abs(got) + std::abs(flipped)));
}

}  // namespace

TEST_CASE(
    "full update T2-iv: N reproduces the single-layer graded truth on an "
    "exactly contracted patch") {
  // The truth here is the graded state of the WHOLE patch with the physical
  // legs open, so the patch is kept at 3x3: the two window placements per
  // direction between them leave no window leg trivial (the first has L_A
  // resp. T_A on the patch boundary, the second R_B resp. B_B), while a
  // single 4x3 patch would cost an order of magnitude more.
  const fue_geom h_left{3, 3, 0, 1};
  const fue_geom h_right{3, 3, 1, 1};
  const fue_geom v_top{3, 3, 1, 0};
  const fue_geom v_bottom{3, 3, 1, 1};
  // Mandated coverage: nontrivial virtual ledger, d = 4, both directions.
  fue_run_t2iv_case<tenes::real_tensor>('h', h_left, "eo", 2, 2, 2100, 23,
                                        "real");
  fue_run_t2iv_case<tenes::real_tensor>('h', h_right, "eo", 4, 4, 2150, 25,
                                        "real");
  fue_run_t2iv_case<tenes::real_tensor>('v', v_top, "eo", 2, 2, 2200, 29,
                                        "real");
  fue_run_t2iv_case<tenes::real_tensor>('v', v_bottom, "eo", 4, 4, 2250, 31,
                                        "real");
  fue_run_t2iv_case<tenes::real_tensor>('h', h_right, "eo", 4, 2, 2300, 37,
                                        "real");
  fue_run_t2iv_case<tenes::complex_tensor>('v', v_top, "eo", 2, 4, 2400, 43,
                                           "complex");
}

// ---- T2-v: bosonic degeneracy --------------------------------------------

namespace {

//! Split the fused (ket, bra) leg of a folded edge tensor into the two
//! legs the bosonic environment builder expects. mptensor flattens
//! column-major and doubled_pipeline_traced() fuses ket-first, so the plain
//! reshape already puts ket on axis 2 and bra on axis 3 - the order
//! Create_Environment_two_sites() contracts against Q1 and conj(Q1).
template <class tensor>
tensor fue_split_edge(const tensor& eT, std::size_t D) {
  const mptensor::Shape sh = eT.shape();
  REQUIRE(sh.size() == 3);
  REQUIRE(sh[2] == D * D);
  return mptensor::reshape(eT, mptensor::Shape(sh[0], sh[1], D, D));
}

template <class tensor>
void fue_run_t2v_case(char dir, const std::string& vps, int seed,
                      const char* type_name) {
  const fue_case<tensor> c =
      fue_make_case<tensor>(dir, fue_full_geom(dir), vps, 2, 2, 2, seed, true,
                            std::string("T2-v ") + type_name);
  INFO(c.label);
  // Premise: this really is the all-even control.
  for (int s = 0; s < c.patch.nsite(); ++s) {
    for (const auto& p : c.sites[s].parity) {
      for (std::size_t i = 0; i < p.size(); ++i) {
        REQUIRE(p[i] == false);
      }
    }
  }
  const auto N = fue_build_N(c.env, c.QA, c.QB, c.dir_e);
  CHECK(N.forbidden_ratio <= 1.0e-10);
  // The mask is the identity when no ledger has an odd value.
  fg_check_allclose(N.N_plain, N.N.t,
                    c.label + " [all-even: N_plain == N]");

  std::vector<tensor> s;
  for (int i = 4; i < 10; ++i) {
    s.push_back(fue_split_edge(c.env[i], c.D));
  }
  tensor got;
  if (dir == 'h') {
    // The window layout of the horizontal fold and of the bosonic builder
    // coincide slot for slot; only Q2 has to be rotated the way
    // Full_update_bond_horizontal() rotates it.
    got = tenes::itps::core::Create_Environment_two_sites(
        c.env[0], c.env[1], c.env[2], c.env[3], s[0], s[1], s[2], s[3], s[4],
        s[5], c.QA.t,
        mptensor::transpose(c.QB.t, mptensor::Axes(3, 0, 1, 2)));
  } else {
    // The bosonic builder only knows the horizontal window, so the vertical
    // one enters rotated by one slot - exactly the wiring
    // iTPS<tensor>::full_update() uses for source_leg == 3:
    //   bosonic (C1..C4)   = fold (C2, C3, C4, C1)
    //   bosonic (eT1..eT6) = fold (eT2, eT3, eT4, eT5, eT6, eT1)
    // and the two factors rotate with it: (l, t, r) -> (t, r, l).
    got = tenes::itps::core::Create_Environment_two_sites(
        c.env[1], c.env[2], c.env[3], c.env[0], s[1], s[2], s[3], s[4], s[5],
        s[0], mptensor::transpose(c.QA.t, mptensor::Axes(1, 2, 0, 3)),
        mptensor::transpose(c.QB.t, mptensor::Axes(3, 1, 2, 0)));
  }
  REQUIRE(fg_max_abs_entry(got) > 0.0);
  fg_check_allclose(N.N.t, got,
                    c.label + " [T2-v N vs Create_Environment_two_sites]");
}

}  // namespace

TEST_CASE(
    "full update T2-v: with all-even ledgers N equals the bosonic "
    "Create_Environment_two_sites") {
  fue_run_t2v_case<tenes::real_tensor>('h', "ee", 2600, "real");
  fue_run_t2v_case<tenes::real_tensor>('v', "ee", 2700, "real");
  fue_run_t2v_case<tenes::complex_tensor>('h', "ee", 2800, "complex");
}
