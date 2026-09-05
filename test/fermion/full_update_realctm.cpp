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

// ===== T2-vi / T2-vii / T3-viii: the REAL CTM environment ==================
//
// Contract: docs/superpowers/specs/2026-09-04-fermion-full-update-contract.md
// sections 3.2 (T2-vi, T2-vii), 3.3 (window selection) and 3.4 (T3-viii).
//
// Every other T2/T3 test closes the fold in an EXACTLY contracted patch. The
// three tests here use the environment the production driver uses:
//   build_reduced_density_tensors(Tn, finfo) -> Calc_CTM_Environment_density
// (the iTPS::update_CTM() path), on a parity-even random 2x2 unit cell with
// D = 2, d = 2 and mixed ledgers, with the ten window tensors picked exactly
// as the measurement path (contract section 3.3) picks them.
//
// Included into the test_fermion_layer TU after fermion/full_update_env.cpp
// and fermion/full_update_bond.cpp; it reuses their scalar judgment helpers
// (fue_check_rel), the pseudo-site packer (fue_pseudo_site), the internal
// identity (fue_internal_identity), the Schmidt counter (fub_bond_schmidt)
// and the exact-cutoff parameter set (fub_exact_params), plus the
// fold_geometry helpers (fg_max_abs_entry, fg_ledger_mixed, ...).
//
// Truth sources (never the machinery under test):
//   - T2-vi:   the measurement path
//              contract_reduced_pair_halves_density_CTM(
//                  build_reduced_pair_halves(Tn1(Y), Tn2(Y), O, dir)),
//              judged against the plain (bosonic-convention) quadratic form
//              of N_plain on the masked block. No reference value is baked
//              in anywhere: both sides are computed here, from the contract's
//              formulas, on the same environment.
//   - T2-vii:  build_pair_state / apply_pair_op, the graded primitives
//              fold_geometry already ties to the Fock oracle.
//   - T3-viii: build_pair_state of the input (identity gate) and the
//              measurement path (hopping gate).
//
// Premises asserted (the test is meaningless without them):
//   - the CTM converged (Calc_CTM_Environment_density returned fewer than
//     Max_CTM_Iteration sweeps) and the folded metric has no forbidden block
//     (forbidden_ratio <= 1e-10) - an unconverged CTM leaks parity, which
//     shows up as an N_plain error and would break the comparison for a
//     reason that is NOT a bug in N;
//   - both QR internal ledgers are mixed (an odd channel exists at all);
//   - the graded QR reconstructs the site (the axis tables were copied
//     right);
//   - the SVD assembly Tn1(X), Tn2(X) reproduces build_pair_state(Tn1, Tn2)
//     up to a scale (contract T2-vi step 1);
//   - every block Y has weight in its (a odd, beta odd) sector, so the
//     mask_{a beta} of step 3 actually does something;
//   - the hopping expectation on every block is non-negligible, so a sign
//     error confined to the odd channel of N_plain would move <G>.
//
// Element comparisons go through the GLOBAL index (contract section 0):
// mptensor evaluates transposes lazily, so local index n of two tensors of
// the same shape need not point at the same element.

#include <complex>

namespace {

// ---- scalar helpers -------------------------------------------------------

inline double fur_re(double v) { return v; }
inline double fur_re(const std::complex<double>& v) { return v.real(); }
inline double fur_im(double) { return 0.0; }
inline double fur_im(const std::complex<double>& v) { return v.imag(); }

//! Deterministic draw, distinct from fue_draw() / ib_random_value() so the
//! fixtures here do not accidentally reproduce another test's tensors.
inline double fur_draw(int seed, const mptensor::Index& idx, std::size_t rank) {
  const double w[6] = {1.9, 3.7, 5.9, 7.3, 9.1, 11.3};
  double x = 0.53 * (seed + 5) + 0.011 * seed * seed;
  for (std::size_t ax = 0; ax < rank; ++ax) {
    x += w[ax] * static_cast<double>(idx[ax]);
  }
  return 0.57 * std::sin(x) + 0.31 * std::cos(1.23 * x + 0.17 * seed);
}

inline void fur_set(tenes::real_tensor& t, const mptensor::Index& idx, int seed,
                    std::size_t rank) {
  t.set_value(idx, fur_draw(seed, idx, rank));
}

inline void fur_set(tenes::complex_tensor& t, const mptensor::Index& idx,
                    int seed, std::size_t rank) {
  t.set_value(idx, std::complex<double>(fur_draw(seed, idx, rank),
                                        fur_draw(seed + 900, idx, rank)));
}

template <class tensor>
tensor fur_from_real(const tenes::real_tensor& a);

template <>
inline tenes::real_tensor fur_from_real<tenes::real_tensor>(
    const tenes::real_tensor& a) {
  return a;
}

template <>
inline tenes::complex_tensor fur_from_real<tenes::complex_tensor>(
    const tenes::real_tensor& a) {
  return fg_to_complex(a);
}

//! Plain elements of exp(tau (c^dag_A c_B + c^dag_B c_A)) for d = 2, legs
//! (in1, in2, out1, out2), in the same element convention as ib_hop_plain(2).
template <class tensor>
tensor fur_hop_gate_plain(double tau) {
  tenes::real_tensor g(mptensor::Shape(2, 2, 2, 2));
  g.set_value(mptensor::Index(0, 0, 0, 0), 1.0);
  g.set_value(mptensor::Index(1, 1, 1, 1), 1.0);
  g.set_value(mptensor::Index(0, 1, 0, 1), std::cosh(tau));
  g.set_value(mptensor::Index(1, 0, 1, 0), std::cosh(tau));
  g.set_value(mptensor::Index(0, 1, 1, 0), std::sinh(tau));
  g.set_value(mptensor::Index(1, 0, 0, 1), std::sinh(tau));
  return fur_from_real<tensor>(g);
}

//! Parity-even random graded tensor over the given ledgers.
template <class tensor>
fg_ftensor<tensor> fur_random_even_block(const fgf::leg_parities& lp,
                                         int seed) {
  mptensor::Shape sh;
  for (const auto& p : lp) {
    sh.push(p.size());
  }
  tensor t(sh);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    if (fgf::count_odd(lp, idx) % 2 == 0) {
      fur_set(t, idx, seed, lp.size());
    }
  }
  return fg_ftensor<tensor>{t, lp};
}

//! Contract T2-vi step 3: (-1) on every element whose legs 0 and 1 are both
//! odd. Ledgers unchanged.
template <class tensor>
fg_ftensor<tensor> fur_mask_ab(const fg_ftensor<tensor>& Y) {
  REQUIRE(Y.rank() >= 2);
  fg_ftensor<tensor> r = Y;
  for (std::size_t n = 0; n < r.t.local_size(); ++n) {
    const mptensor::Index idx = r.t.global_index(n);
    if (Y.parity[0][idx[0]] && Y.parity[1][idx[1]]) {
      typename tensor::value_type v;
      r.t.get_value(idx, v);
      r.t.set_value(idx, -v);
    }
  }
  return r;
}

//! Largest |element| in the sector where legs 0 and 1 are both odd - the
//! sector fur_mask_ab() acts on.
template <class tensor>
double fur_max_abs_oddodd(const fg_ftensor<tensor>& Y) {
  double m = 0.0;
  for (std::size_t n = 0; n < Y.t.local_size(); ++n) {
    const mptensor::Index idx = Y.t.global_index(n);
    if (Y.parity[0][idx[0]] && Y.parity[1][idx[1]]) {
      typename tensor::value_type v;
      Y.t.get_value(idx, v);
      m = std::max(m, std::abs(v));
    }
  }
  return m;
}

/*! Least-squares proportionality: c = <ref, got> / <ref, ref>, and the
 *  deviation max |got - c ref| / max |ref|, both through the global index.
 *  This is the T3-viii judgment (and the T2-vi step-1 premise): the overall
 *  scale and phase are absorbed into c.
 */
template <class tensor>
double fur_lstsq_deviation(const tensor& got, const tensor& ref,
                           typename tensor::value_type& c) {
  using value_type = typename tensor::value_type;
  REQUIRE(got.shape() == ref.shape());
  value_type num = value_type(0.0);
  double den = 0.0;
  double ref_max = 0.0;
  for (std::size_t n = 0; n < ref.local_size(); ++n) {
    const mptensor::Index idx = ref.global_index(n);
    value_type r, g;
    ref.get_value(idx, r);
    got.get_value(idx, g);
    num += fgf::detail::scalar_conj(r) * g;
    den += std::abs(r) * std::abs(r);
    ref_max = std::max(ref_max, std::abs(r));
  }
  REQUIRE(den > 0.0);
  c = num / den;
  double dev = 0.0;
  mptensor::Index worst;
  for (std::size_t n = 0; n < ref.local_size(); ++n) {
    const mptensor::Index idx = ref.global_index(n);
    value_type r, g;
    ref.get_value(idx, r);
    got.get_value(idx, g);
    const double d = std::abs(g - c * r);
    if (d > dev) {
      dev = d;
      worst = idx;
    }
  }
  INFO("[lstsq] c=" << c << " max|got - c ref|=" << dev << " max|ref|="
                    << ref_max << " at " << ib_index_string(worst));
  return dev / ref_max;
}

//! Contract section 3.1 relative error of `got` against the reference:
//! max |got - ref| / max(1, max |ref|), through the global index.
template <class tensor>
double fur_rel_diff(const tensor& got, const tensor& ref) {
  REQUIRE(got.shape() == ref.shape());
  const double scale = std::max(1.0, fg_max_abs_entry(ref));
  double dev = 0.0;
  mptensor::Index worst;
  for (std::size_t n = 0; n < ref.local_size(); ++n) {
    const mptensor::Index idx = ref.global_index(n);
    typename tensor::value_type g, r;
    got.get_value(idx, g);
    ref.get_value(idx, r);
    const double d = std::abs(g - r);
    if (d > dev) {
      dev = d;
      worst = idx;
    }
  }
  INFO("[rel_diff] max|got - ref|=" << dev << " scale=" << scale << " at "
                                    << ib_index_string(worst));
  return dev / scale;
}

// ---- the real CTM environment ---------------------------------------------

template <class tensor>
struct fur_env {
  tenes::SquareLattice lattice{2, 2};
  fgf::FermionInfo finfo;
  tenes::itps::PEPS_Parameters params;
  std::vector<tensor> Tn;       // plain rank-5 site tensors
  std::vector<tensor> reduced;  // build_reduced_density_tensors(Tn, finfo)
  std::vector<tensor> C1, C2, C3, C4, eTt, eTr, eTb, eTl;
  int ctm_count = 0;
  std::string label;
};

/*! Parity-even random 2x2 unit cell (D = 2, d = 2, every ledger [even, odd])
 *  and its converged CTM environment, built through the production path.
 *
 *  Premise asserted here: the CTM converged before Max_CTM_Iteration. The
 *  initial corners of Calc_CTM_Environment_density() are made with all-ones
 *  boundary vectors and therefore carry a parity-odd component; only the
 *  converged environment is parity clean, and only a parity-clean
 *  environment gives the fold an empty forbidden block.
 */
template <class tensor>
fur_env<tensor> fur_make_env(int seed, int chi, const char* type_name) {
  fur_env<tensor> e;
  e.label = std::string("realctm ") + type_name +
            " seed=" + std::to_string(seed) + " chi=" + std::to_string(chi);
  INFO(e.label);
  const int N = e.lattice.N_UNIT;
  REQUIRE(N == 4);
  for (int s = 0; s < N; ++s) {
    e.lattice.physical_dims[s] = 2;
    e.lattice.virtual_dims[s] = {2, 2, 2, 2};
    e.lattice.initial_dirs[s] = {0.0};
    e.lattice.noises[s] = 0.0;
  }

  const fgf::parity_vector vp{false, true};
  const fgf::parity_vector phys = ib_phys_parity(2);
  e.finfo.enabled = true;
  e.finfo.phys.assign(N, phys);
  e.finfo.virt.assign(N, std::array<fgf::parity_vector, 4>{vp, vp, vp, vp});
  fgf::validate_neighbor_consistency(e.finfo, e.lattice);

  for (int s = 0; s < N; ++s) {
    const fgf::leg_parities lp = fgf::Tn_parity(e.finfo, s);
    tensor t(mptensor::Shape(2, 2, 2, 2, 2));
    for (std::size_t n = 0; n < t.local_size(); ++n) {
      const mptensor::Index idx = t.global_index(n);
      if (fgf::count_odd(lp, idx) % 2 == 0) {
        fur_set(t, idx, seed + 7 * s, 5);
      }
    }
    const fg_ftensor<tensor> wrapped{t, lp};
    REQUIRE(fgf::parity_violation(wrapped) == 0.0);
    REQUIRE(fgf::max_abs(wrapped) > 0.0);
    e.Tn.push_back(t);
  }

  e.params.print_level = tenes::PrintLevel::none;
  e.params.fermion = true;
  e.params.phys_parity.assign(N, phys);
  e.params.D = 2;
  e.params.CHI = chi;
  e.params.Use_RSVD = false;
  e.params.Max_CTM_Iteration = 2000;
  e.params.CTM_Convergence_Epsilon = 1.0e-12;

  e.reduced = fgf::build_reduced_density_tensors(e.Tn, e.finfo);
  e.C1.assign(N, tensor());
  e.C2.assign(N, tensor());
  e.C3.assign(N, tensor());
  e.C4.assign(N, tensor());
  e.eTt.assign(N, tensor());
  e.eTr.assign(N, tensor());
  e.eTb.assign(N, tensor());
  e.eTl.assign(N, tensor());
  e.ctm_count = tenes::itps::core::Calc_CTM_Environment_density(
      e.C1, e.C2, e.C3, e.C4, e.eTt, e.eTr, e.eTb, e.eTl, e.reduced, e.params,
      e.lattice);
  INFO(e.label << " [premise: CTM sweeps = " << e.ctm_count << " of at most "
               << e.params.Max_CTM_Iteration << "]");
  REQUIRE(e.ctm_count < e.params.Max_CTM_Iteration);
  return e;
}

// ---- one full-update window in that environment ---------------------------

template <class tensor>
struct fur_window {
  char dir = 'h';
  fgf::reduced_pair_direction dir_e = fgf::reduced_pair_direction::horizontal;
  int s1 = 0;
  int s2 = 0;
  std::vector<tensor> env;      // C1..C4, eT1..eT6 in argument order
  fg_ftensor<tensor> Tn1, Tn2;  // wrapped site tensors, legs (l, t, r, b, s)
  fg_ftensor<tensor> QA, RA, QB, RB;
  fgf::parity_vector pa, pb, ps1, ps2;
  fg_ftensor<tensor> X;  // the state's own block (a, beta, s1, s2)
  std::string label;
};

//! Contract section 3.3: the ten window tensors, s1 the raster-earlier site.
template <class tensor>
std::vector<tensor> fur_window_env(const fur_env<tensor>& e, int s1, int s2,
                                   char dir) {
  if (dir == 'h') {
    return {e.C1[s1],  e.C2[s2],  e.C3[s2],  e.C4[s1],  e.eTt[s1],
            e.eTt[s2], e.eTr[s2], e.eTb[s2], e.eTb[s1], e.eTl[s1]};
  }
  return {e.C1[s1],  e.C2[s1],  e.C3[s2],  e.C4[s2],  e.eTt[s1],
          e.eTr[s1], e.eTr[s2], e.eTb[s2], e.eTl[s2], e.eTl[s1]};
}

//! The measurement path: <psi(T1, T2)| O |psi(T1, T2)> closed in the window.
template <class tensor>
typename tensor::value_type fur_closed(const fur_window<tensor>& w,
                                       const fg_ftensor<tensor>& T1,
                                       const fg_ftensor<tensor>& T2,
                                       const fg_ftensor<tensor>& O) {
  const auto halves = fgf::build_reduced_pair_halves(T1, T2, O, w.dir_e);
  return fgf::contract_reduced_pair_halves_density_CTM(
      w.env[0], w.env[1], w.env[2], w.env[3], w.env[4], w.env[5], w.env[6],
      w.env[7], w.env[8], w.env[9], halves);
}

template <class tensor>
fue::full_update_environment<tensor> fur_build_N(const fur_window<tensor>& w) {
  return fue::build_full_update_environment(
      w.env[0], w.env[1], w.env[2], w.env[3], w.env[4], w.env[5], w.env[6],
      w.env[7], w.env[8], w.env[9], w.QA, w.QB, w.dir_e);
}

/*! Build the window whose raster-earlier site is s1, with the QR factors of
 *  contract T2-vi (bosonic axes) and the state's own block X.
 *
 *  Premises asserted: both internal ledgers are mixed; the graded QR
 *  reconstructs each site tensor (so the axis tables really were copied
 *  from the contract); the measurement-path norm is positive.
 */
template <class tensor>
fur_window<tensor> fur_make_window(const fur_env<tensor>& e, int s1, char dir,
                                   const std::string& tag) {
  fur_window<tensor> w;
  w.dir = dir;
  w.dir_e = dir == 'h' ? fgf::reduced_pair_direction::horizontal
                       : fgf::reduced_pair_direction::vertical;
  w.s1 = s1;
  w.s2 = dir == 'h' ? e.lattice.right(s1) : e.lattice.bottom(s1);
  REQUIRE(w.s2 != w.s1);
  w.label = tag + " " + e.label + " dir=" + dir +
            " s1=" + std::to_string(w.s1) + " s2=" + std::to_string(w.s2);
  INFO(w.label);
  w.env = fur_window_env(e, w.s1, w.s2, dir);
  w.Tn1 = fgf::wrap_Tn(e.Tn[w.s1], e.finfo, w.s1);
  w.Tn2 = fgf::wrap_Tn(e.Tn[w.s2], e.finfo, w.s2);
  REQUIRE(fgf::parity_violation(w.Tn1) == 0.0);
  REQUIRE(fgf::parity_violation(w.Tn2) == 0.0);

  // QR with the bosonic axes (contract T2-vi table).
  if (dir == 'h') {
    REQUIRE(fgf::qr(w.Tn1, mptensor::Axes(0, 1, 3), mptensor::Axes(2, 4), w.QA,
                    w.RA) == 0);
    REQUIRE(fgf::qr(w.Tn2, mptensor::Axes(1, 2, 3), mptensor::Axes(0, 4), w.QB,
                    w.RB) == 0);
  } else {
    REQUIRE(fgf::qr(w.Tn1, mptensor::Axes(0, 1, 2), mptensor::Axes(3, 4), w.QA,
                    w.RA) == 0);
    REQUIRE(fgf::qr(w.Tn2, mptensor::Axes(0, 2, 3), mptensor::Axes(1, 4), w.QB,
                    w.RB) == 0);
  }
  REQUIRE(w.QA.rank() == 4);
  REQUIRE(w.QB.rank() == 4);
  REQUIRE(w.RA.rank() == 3);
  REQUIRE(w.RB.rank() == 3);
  w.pa = w.QA.parity[3];
  w.pb = w.QB.parity[3];
  w.ps1 = w.Tn1.parity[4];
  w.ps2 = w.Tn2.parity[4];
  REQUIRE(fg_ledger_mixed(w.pa));
  REQUIRE(fg_ledger_mixed(w.pb));
  REQUIRE(fg_ledger_mixed(w.ps1));
  REQUIRE(fg_ledger_mixed(w.ps2));

  // Premise: Q R reproduces the site, in the contract's leg orders
  //   horizontal: QA(l,t,b,a) RA(a,r,s1)   QB(t,r,b,beta) RB(beta,l,s2)
  //   vertical:   QA(l,t,r,a) RA(a,b,s1)   QB(l,r,b,beta) RB(beta,t,s2)
  {
    fg_ftensor<tensor> back1, back2;
    if (dir == 'h') {
      back1 = fgf::transpose(
          fgf::tensordot(w.QA, w.RA, mptensor::Axes(3), mptensor::Axes(0)),
          mptensor::Axes(0, 1, 3, 2, 4));
      back2 = fgf::transpose(
          fgf::tensordot(w.QB, w.RB, mptensor::Axes(3), mptensor::Axes(0)),
          mptensor::Axes(3, 0, 1, 2, 4));
    } else {
      back1 = fgf::tensordot(w.QA, w.RA, mptensor::Axes(3), mptensor::Axes(0));
      back2 = fgf::transpose(
          fgf::tensordot(w.QB, w.RB, mptensor::Axes(3), mptensor::Axes(0)),
          mptensor::Axes(0, 3, 1, 2, 4));
    }
    REQUIRE(back1.parity == w.Tn1.parity);
    REQUIRE(back2.parity == w.Tn2.parity);
    const double d1 = fur_rel_diff(back1.t, w.Tn1.t);
    const double d2 = fur_rel_diff(back2.t, w.Tn2.t);
    INFO(w.label << " [premise: QR reconstruction] " << d1 << " / " << d2);
    REQUIRE(d1 <= 1.0e-12);
    REQUIRE(d2 <= 1.0e-12);
  }

  // X = transpose(tensordot(RA, RB, Axes(1), Axes(1)), Axes(0,2,1,3)):
  // legs (a, beta, s1, s2).
  w.X = fgf::transpose(
      fgf::tensordot(w.RA, w.RB, mptensor::Axes(1), mptensor::Axes(1)),
      mptensor::Axes(0, 2, 1, 3));
  REQUIRE(w.X.rank() == 4);
  const fgf::leg_parities block_ledgers{w.pa, w.pb, w.ps1, w.ps2};
  REQUIRE(w.X.parity == block_ledgers);
  REQUIRE(fgf::max_abs(w.X) > 0.0);

  // Premise: the measurement path gives a positive norm in this window.
  {
    const fg_ftensor<tensor> I = fgf::wrap_twosite_gate(
        fue_internal_identity<tensor>(w.ps1, w.ps2), w.ps1, w.ps2);
    const auto norm = fur_closed(w, w.Tn1, w.Tn2, I);
    INFO(w.label << " [premise: measurement-path norm = " << norm << "]");
    REQUIRE(std::isfinite(std::abs(norm)));
    REQUIRE(fur_re(norm) > 0.0);
    REQUIRE(std::abs(fur_im(norm)) <= 1.0e-10 * std::abs(norm));
  }
  return w;
}

/*! Contract T2-vi step 1: Tn1(Y), Tn2(Y) from the graded SVD of the block.
 *  Copied from the contract's table; whether the copy is right is asserted
 *  by fur_require_assembly() before any comparison uses it.
 */
template <class tensor>
void fur_assemble(const fur_window<tensor>& w, const fg_ftensor<tensor>& Y,
                  fg_ftensor<tensor>& T1, fg_ftensor<tensor>& T2) {
  REQUIRE(Y.rank() == 4);
  fg_ftensor<tensor> U, VT;
  std::vector<double> s;
  REQUIRE(fgf::svd(Y, mptensor::Axes(0, 2), mptensor::Axes(1, 3), U, s, VT) ==
          0);
  REQUIRE(!s.empty());
  fg_ftensor<tensor> R1 = U;  // (a, s1, m)
  R1.multiply_vector(s, 2);
  const fg_ftensor<tensor> R2 =
      fgf::transpose(VT, mptensor::Axes(1, 2, 0));  // (beta, s2, m)
  if (w.dir == 'h') {
    T1 = fgf::transpose(
        fgf::tensordot(w.QA, R1, mptensor::Axes(3), mptensor::Axes(0)),
        mptensor::Axes(0, 1, 4, 2, 3));
    T2 = fgf::transpose(
        fgf::tensordot(w.QB, R2, mptensor::Axes(3), mptensor::Axes(0)),
        mptensor::Axes(4, 0, 1, 2, 3));
  } else {
    T1 = fgf::transpose(
        fgf::tensordot(w.QA, R1, mptensor::Axes(3), mptensor::Axes(0)),
        mptensor::Axes(0, 1, 2, 4, 3));
    T2 = fgf::transpose(
        fgf::tensordot(w.QB, R2, mptensor::Axes(3), mptensor::Axes(0)),
        mptensor::Axes(0, 4, 1, 2, 3));
  }
  REQUIRE(T1.rank() == 5);
  REQUIRE(T2.rank() == 5);
  const int a_bond = w.dir == 'h' ? 2 : 3;
  const int b_bond = w.dir == 'h' ? 0 : 1;
  REQUIRE(T1.parity[a_bond] == T2.parity[b_bond]);
  for (int ax = 0; ax < 5; ++ax) {
    if (ax != a_bond) {
      REQUIRE(T1.parity[ax] == w.Tn1.parity[ax]);
    }
    if (ax != b_bond) {
      REQUIRE(T2.parity[ax] == w.Tn2.parity[ax]);
    }
  }
}

//! Contract T2-vi step 1, premise: with Y = X the assembly reproduces
//! build_pair_state(Tn1, Tn2) up to an overall scale.
template <class tensor>
void fur_require_assembly(const fur_window<tensor>& w) {
  fg_ftensor<tensor> T1, T2;
  fur_assemble(w, w.X, T1, T2);
  const fg_ftensor<tensor> got = fgf::build_pair_state(T1, T2, w.dir_e);
  const fg_ftensor<tensor> want = fgf::build_pair_state(w.Tn1, w.Tn2, w.dir_e);
  REQUIRE(got.shape() == want.shape());
  REQUIRE(got.parity == want.parity);
  REQUIRE(fgf::max_abs(want) > 0.0);
  typename tensor::value_type c;
  const double dev = fur_lstsq_deviation(got.t, want.t, c);
  INFO(w.label << " [premise: assembly Tn1(X), Tn2(X) vs build_pair_state] "
               << "dev=" << dev << " c=" << c);
  REQUIRE(std::abs(c) > 0.0);
  REQUIRE(dev <= 1.0e-12);
}

//! Bosonic-convention quadratic form <Yt| N_plain |Zt> of contract T2-vi
//! step 3, plain mptensor operations only.
template <class tensor>
typename tensor::value_type fur_quad(const tensor& N_plain, const tensor& Yt,
                                     const tensor& Zt) {
  const tensor NZ = mptensor::tensordot(N_plain, Zt, mptensor::Axes(0, 1),
                                        mptensor::Axes(0, 1));
  return mptensor::trace(NZ, mptensor::conj(Yt), mptensor::Axes(0, 1, 2, 3),
                         mptensor::Axes(0, 1, 2, 3));
}

// ---- T2-vi ---------------------------------------------------------------

/*! One block Y: the measurement path on Tn1(Y), Tn2(Y) against the plain
 *  quadratic form of N_plain on the masked block.
 */
template <class tensor>
void fur_run_t2vi_block(const fur_window<tensor>& w,
                        const fue::full_update_environment<tensor>& N,
                        const fg_ftensor<tensor>& I,
                        const fg_ftensor<tensor>& G,
                        const fg_ftensor<tensor>& hop,
                        const fg_ftensor<tensor>& Y, const std::string& what) {
  const std::string label = w.label + " Y=" + what;
  INFO(label);
  REQUIRE(Y.rank() == 4);
  const fgf::leg_parities block_ledgers{w.pa, w.pb, w.ps1, w.ps2};
  REQUIRE(Y.parity == block_ledgers);
  const double ymax = fgf::max_abs(Y);
  REQUIRE(ymax > 0.0);
  REQUIRE(fgf::parity_violation(Y) <= 1.0e-14 * ymax);
  // Premise: the mask of step 3 has something to act on.
  const double oddodd = fur_max_abs_oddodd(Y);
  INFO(label << " [premise: (a odd, beta odd) weight = " << oddodd << " of max "
             << ymax << "]");
  REQUIRE(oddodd > 1.0e-6 * ymax);

  // Step 1: the state psi(Y).
  fg_ftensor<tensor> T1, T2;
  fur_assemble(w, Y, T1, T2);

  // Step 2: the measurement path.
  const auto norm_meas = fur_closed(w, T1, T2, I);
  const auto g_meas = fur_closed(w, T1, T2, G);
  const auto hop_meas = fur_closed(w, T1, T2, hop);
  INFO(label << " [measurement path] norm=" << norm_meas
             << " <G>*norm=" << g_meas << " <hop>*norm=" << hop_meas);
  REQUIRE(std::isfinite(std::abs(norm_meas)));
  REQUIRE(std::abs(norm_meas) > 0.0);
  // Premise: the odd (hopping) channel carries weight on this block, so a
  // sign error confined to it would move <G> = cosh + sinh <hop>.
  REQUIRE(std::abs(hop_meas) > 1.0e-3 * std::abs(norm_meas));

  // Step 3: the plain quadratic form on the masked block.
  const fg_ftensor<tensor> Yt = fur_mask_ab(Y);
  const fg_ftensor<tensor> Th =
      fgf::tensordot(Y, G, mptensor::Axes(2, 3), mptensor::Axes(0, 1));
  REQUIRE(Th.parity == Y.parity);
  const fg_ftensor<tensor> Tht = fur_mask_ab(Th);
  const auto norm_quad = fur_quad(N.N_plain, Yt.t, Yt.t);
  const auto g_quad = fur_quad(N.N_plain, Yt.t, Tht.t);
  INFO(label << " [quadratic form] norm=" << norm_quad
             << " <G>*norm=" << g_quad);
  REQUIRE(std::abs(norm_quad) > 0.0);

  // Step 4.
  fue_check_rel(label + " [T2-vi norm]", norm_quad, norm_meas, 1.0e-10);
  const auto r_meas = g_meas / norm_meas;
  const auto r_quad = g_quad / norm_quad;
  fue_check_rel(label + " [T2-vi <G>]", r_quad, r_meas, 1.0e-10);
}

template <class tensor>
void fur_run_t2vi_window(const fur_env<tensor>& e, int s1, char dir, int seed) {
  const fur_window<tensor> w = fur_make_window(e, s1, dir, "T2-vi");
  INFO(w.label);
  fur_require_assembly(w);

  const auto N = fur_build_N(w);
  INFO(w.label << " [premise: forbidden_ratio = " << N.forbidden_ratio << "]");
  REQUIRE(N.forbidden_ratio <= 1.0e-10);
  REQUIRE(N.N_plain.shape() ==
          mptensor::Shape(w.pa.size(), w.pb.size(), w.pa.size(), w.pb.size()));

  const fg_ftensor<tensor> I = fgf::wrap_twosite_gate(
      fue_internal_identity<tensor>(w.ps1, w.ps2), w.ps1, w.ps2);
  // A strongly non-trivial hopping gate: the odd channel enters <G> with
  // weight sinh(tau), and tau = 0.3 keeps it well above the tolerance.
  const fg_ftensor<tensor> G =
      fgf::wrap_twosite_gate(fur_hop_gate_plain<tensor>(0.3), w.ps1, w.ps2);
  const fg_ftensor<tensor> hop = fgf::wrap_twosite_gate(
      fur_from_real<tensor>(ib_hop_plain(2)), w.ps1, w.ps2);

  // The blocks of contract T2-vi: X, Theta = X G, two random even blocks,
  // and X + 0.05 * random.
  const fgf::leg_parities lp{w.pa, w.pb, w.ps1, w.ps2};
  const fg_ftensor<tensor>& X = w.X;
  const fg_ftensor<tensor> Theta =
      fgf::tensordot(X, G, mptensor::Axes(2, 3), mptensor::Axes(0, 1));
  const fg_ftensor<tensor> R1 = fur_random_even_block<tensor>(lp, seed + 1);
  const fg_ftensor<tensor> R2 = fur_random_even_block<tensor>(lp, seed + 2);
  fg_ftensor<tensor> Xp = X;
  {
    fg_ftensor<tensor> R3 = fur_random_even_block<tensor>(lp, seed + 3);
    const double scale = 0.05 * fgf::max_abs(X) / fgf::max_abs(R3);
    for (std::size_t n = 0; n < Xp.t.local_size(); ++n) {
      const mptensor::Index idx = Xp.t.global_index(n);
      typename tensor::value_type x, r;
      Xp.t.get_value(idx, x);
      R3.t.get_value(idx, r);
      Xp.t.set_value(idx, x + scale * r);
    }
    // Anti-vacuity: the perturbation really moved X.
    REQUIRE(fur_rel_diff(Xp.t, X.t) > 1.0e-3);
  }

  fur_run_t2vi_block(w, N, I, G, hop, X, "X");
  fur_run_t2vi_block(w, N, I, G, hop, Theta, "Theta=X.G");
  fur_run_t2vi_block(w, N, I, G, hop, R1, "random-1");
  fur_run_t2vi_block(w, N, I, G, hop, R2, "random-2");
  fur_run_t2vi_block(w, N, I, G, hop, Xp, "X+0.05*random");
}

}  // namespace

TEST_CASE(
    "full update T2-vi: in the real CTM environment the plain quadratic form "
    "of N_plain reproduces the measurement path on general blocks") {
  {
    const fur_env<tenes::real_tensor> e =
        fur_make_env<tenes::real_tensor>(5100, 8, "real");
    fur_run_t2vi_window(e, 0, 'h', 5110);
    fur_run_t2vi_window(e, 3, 'v', 5120);
  }
}

// ---- T2-vii --------------------------------------------------------------

namespace {

template <class tensor>
void fur_run_t2vii_window(const fur_env<tensor>& e, int s1, char dir) {
  const fur_window<tensor> w = fur_make_window(e, s1, dir, "T2-vii");
  INFO(w.label);

  // The pseudo-sites of T2-i: the internal leg in the physical slot, an
  // even dimension-1 dummy in the connecting-bond slot.
  const fg_ftensor<tensor> QAp = fue_pseudo_site(w.QA, dir == 'h' ? 2 : 3);
  const fg_ftensor<tensor> QBp = fue_pseudo_site(w.QB, dir == 'h' ? 0 : 1);
  REQUIRE(fgf::parity_violation(QAp) <= 1.0e-13 * fgf::max_abs(QAp));
  REQUIRE(fgf::parity_violation(QBp) <= 1.0e-13 * fgf::max_abs(QBp));

  const fg_ftensor<tensor> want = fgf::build_pair_state(w.Tn1, w.Tn2, w.dir_e);
  const fg_ftensor<tensor> pseudo = fgf::build_pair_state(QAp, QBp, w.dir_e);
  REQUIRE(pseudo.rank() == 8);
  // The internal legs sit where the physical legs of a pair state sit.
  REQUIRE(pseudo.parity[3] == w.pa);
  REQUIRE(pseudo.parity[7] == w.pb);
  const fg_ftensor<tensor> got = fgf::apply_pair_op(pseudo, w.X);

  REQUIRE(got.shape() == want.shape());
  REQUIRE(got.parity == want.parity);
  REQUIRE(fgf::max_abs(want) > 0.0);
  // Premise: X reaches the (a odd, beta odd) sector, so the pairing of the
  // two internal legs through the pseudo-sites is exercised with a sign.
  const double oddodd = fur_max_abs_oddodd(w.X);
  INFO(w.label << " [premise: X (a odd, beta odd) weight = " << oddodd
               << " of max " << fgf::max_abs(w.X) << "]");
  REQUIRE(oddodd > 1.0e-6 * fgf::max_abs(w.X));

  const double dev = fur_rel_diff(got.t, want.t);
  INFO(w.label << " [T2-vii apply_pair_op(pair(QA', QB'), X) == "
                  "pair(Tn1, Tn2)] relative error = "
               << dev);
  CHECK(dev <= 1.0e-12);
}

}  // namespace

TEST_CASE(
    "full update T2-vii: the pseudo-site pair state with X applied equals "
    "the pair state of the sites") {
  {
    const fur_env<tenes::real_tensor> e =
        fur_make_env<tenes::real_tensor>(5300, 8, "real");
    fur_run_t2vii_window(e, 0, 'h');
    fur_run_t2vii_window(e, 3, 'v');
  }
}

// ---- T3-viii -------------------------------------------------------------

namespace {

//! Contract T3-iii shape: rank 5, clean parity, every non-bond leg keeps its
//! dimension and ledger, the bond ledger agrees at both ends.
template <class tensor>
void fur_check_output_ledgers(const fur_window<tensor>& w,
                              const fg_ftensor<tensor>& T1n,
                              const fg_ftensor<tensor>& T2n,
                              const std::string& label) {
  INFO(label << " [output ledgers]");
  REQUIRE(T1n.rank() == 5);
  REQUIRE(T2n.rank() == 5);
  const int a_bond = w.dir == 'h' ? 2 : 3;
  const int b_bond = w.dir == 'h' ? 0 : 1;
  for (int ax = 0; ax < 5; ++ax) {
    if (ax != a_bond) {
      CHECK(T1n.shape()[ax] == w.Tn1.shape()[ax]);
      CHECK(T1n.parity[ax] == w.Tn1.parity[ax]);
    }
    if (ax != b_bond) {
      CHECK(T2n.shape()[ax] == w.Tn2.shape()[ax]);
      CHECK(T2n.parity[ax] == w.Tn2.parity[ax]);
    }
  }
  REQUIRE(T1n.shape()[a_bond] == T2n.shape()[b_bond]);
  CHECK(T1n.parity[a_bond] == T2n.parity[b_bond]);
  REQUIRE(fgf::max_abs(T1n) > 0.0);
  REQUIRE(fgf::max_abs(T2n) > 0.0);
  REQUIRE(std::isfinite(fgf::max_abs(T1n)));
  REQUIRE(std::isfinite(fgf::max_abs(T2n)));
  INFO(label << " [parity violation] " << fgf::parity_violation(T1n) << " / "
             << fgf::parity_violation(T2n));
  CHECK(fgf::parity_violation(T1n) <= 1.0e-12 * fgf::max_abs(T1n));
  CHECK(fgf::parity_violation(T2n) <= 1.0e-12 * fgf::max_abs(T2n));
}

template <class tensor>
void fur_run_bond(const fur_window<tensor>& w, const fg_ftensor<tensor>& gate,
                  const tenes::itps::PEPS_Parameters& params,
                  fg_ftensor<tensor>& T1n, fg_ftensor<tensor>& T2n) {
  tenes::itps::core::Full_update_bond_fermion(
      w.env[0], w.env[1], w.env[2], w.env[3], w.env[4], w.env[5], w.env[6],
      w.env[7], w.env[8], w.env[9], w.Tn1, w.Tn2, gate, w.dir_e, params, T1n,
      T2n);
}

template <class tensor>
void fur_run_t3viii_window(const fur_env<tensor>& e, int s1, char dir) {
  const fur_window<tensor> w = fur_make_window(e, s1, dir, "T3-viii");
  INFO(w.label);

  // Premises on the environment: no forbidden block, and a positive
  // DEFINITE metric (otherwise "unchanged" is not an elementwise statement:
  // any state equal to the input in the N-norm would be an equally good
  // answer).
  {
    const auto N = fur_build_N(w);
    INFO(w.label << " [premise: forbidden_ratio = " << N.forbidden_ratio
                 << "]");
    REQUIRE(N.forbidden_ratio <= 1.0e-10);
    std::vector<double> wv;
    REQUIRE(mptensor::eigh(N.N_plain, mptensor::Axes(0, 1),
                           mptensor::Axes(2, 3), wv) == 0);
    REQUIRE(!wv.empty());
    double wmax = wv[0];
    double wmin = wv[0];
    for (const double v : wv) {
      wmax = std::max(wmax, v);
      wmin = std::min(wmin, v);
    }
    INFO(w.label << " [premise: metric spectrum] min=" << wmin
                 << " max=" << wmax << " size=" << wv.size());
    REQUIRE(wmax > 0.0);
    REQUIRE(wmin > 1.0e-10 * wmax);
  }

  const fg_ftensor<tensor> pair_old =
      fgf::build_pair_state(w.Tn1, w.Tn2, w.dir_e);
  // Premise: the state uses both parity sectors of the connecting bond, so
  // a sign error in either channel of the update would change it.
  {
    const fub_schmidt s0 = fub_bond_schmidt(pair_old);
    INFO(w.label << " [premise: bond Schmidt] total=" << s0.total
                 << " even=" << s0.even << " odd=" << s0.odd);
    REQUIRE(s0.even >= 1);
    REQUIRE(s0.odd >= 1);
  }

  const tenes::itps::PEPS_Parameters params = fub_exact_params(true);
  const fg_ftensor<tensor> I = fgf::wrap_twosite_gate(
      fue_internal_identity<tensor>(w.ps1, w.ps2), w.ps1, w.ps2);

  // ---- identity gate: the pair state is unchanged up to a scale --------
  {
    fg_ftensor<tensor> T1n, T2n;
    fur_run_bond(w, I, params, T1n, T2n);
    fur_check_output_ledgers(w, T1n, T2n, w.label + " [identity]");
    const fg_ftensor<tensor> pair_new =
        fgf::build_pair_state(T1n, T2n, w.dir_e);
    REQUIRE(pair_new.shape() == pair_old.shape());
    REQUIRE(pair_new.parity == pair_old.parity);
    typename tensor::value_type c;
    const double dev = fur_lstsq_deviation(pair_new.t, pair_old.t, c);
    INFO(w.label << " [T3-viii identity] |pair_new - c pair_old| / "
                    "max|pair_old| = "
                 << dev << " (c = " << c << ")");
    REQUIRE(std::abs(c) > 0.0);
    CHECK(dev <= 1.0e-10);
  }

  // ---- hopping gate: the bond's <G> goes up --------------------------
  {
    const fg_ftensor<tensor> G =
        fgf::wrap_twosite_gate(fur_hop_gate_plain<tensor>(0.01), w.ps1, w.ps2);
    fg_ftensor<tensor> T1n, T2n;
    fur_run_bond(w, G, params, T1n, T2n);
    fur_check_output_ledgers(w, T1n, T2n, w.label + " [hopping]");

    const auto norm_old = fur_closed(w, w.Tn1, w.Tn2, I);
    const auto g_old = fur_closed(w, w.Tn1, w.Tn2, G);
    const auto norm_new = fur_closed(w, T1n, T2n, I);
    const auto g_new = fur_closed(w, T1n, T2n, G);
    REQUIRE(fur_re(norm_old) > 0.0);
    REQUIRE(fur_re(norm_new) > 0.0);
    REQUIRE(std::abs(fur_im(norm_new)) <= 1.0e-10 * std::abs(norm_new));
    const auto r_old = g_old / norm_old;
    const auto r_new = g_new / norm_new;
    INFO(w.label << " [T3-viii hopping] <G>_old=" << r_old
                 << " <G>_new=" << r_new);
    REQUIRE(std::isfinite(fur_re(r_old)));
    REQUIRE(std::isfinite(fur_re(r_new)));
    REQUIRE(std::abs(fur_im(r_old)) <= 1.0e-10 * std::abs(r_old));
    REQUIRE(std::abs(fur_im(r_new)) <= 1.0e-10 * std::abs(r_new));
    // Anti-vacuity: the update really changed the state.
    {
      const fg_ftensor<tensor> pair_new =
          fgf::build_pair_state(T1n, T2n, w.dir_e);
      typename tensor::value_type c;
      const double dev = fur_lstsq_deviation(pair_new.t, pair_old.t, c);
      INFO(w.label << " [hopping: state moved] dev=" << dev);
      REQUIRE(dev > 1.0e-6);
    }
    // Not a bare ">": T3-iv measured the two ALS routes disagreeing on the
    // solution by up to 5e-8 (als_iterate stops on the cost, not on the
    // iterate), so a rise smaller than that would be round-off rather than
    // the update. The measured rise here is 1.6e-5 (h) and 2.5e-5 (v), so a
    // floor of 1e-7 sits two orders of magnitude above that noise and two
    // below the signal.
    CHECK(fur_re(r_new) > fur_re(r_old) + 1.0e-7);
  }
}

}  // namespace

TEST_CASE(
    "full update T3-viii: in the real CTM environment the identity gate "
    "leaves the pair state and the hopping gate raises the bond's <G>") {
  {
    const fur_env<tenes::real_tensor> e =
        fur_make_env<tenes::real_tensor>(5500, 8, "real");
    fur_run_t3viii_window(e, 0, 'h');
    fur_run_t3viii_window(e, 3, 'v');
  }
}
