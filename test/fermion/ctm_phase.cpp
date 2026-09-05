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

// ===== T8 / T2-vi complex: the overall phase of the fermionic CTM ==========
//
// Contract: docs/superpowers/specs/2026-09-04-fermion-full-update-contract.md
// section 2.5 (API), 3.6 (T8-i, T8-ii, T8-iv, T2-vi complex), 3.3 (window
// selection) and 3.2 (T2-vi).
//
// The fold-mode CTM (build_reduced_density_tensors ->
// Calc_CTM_Environment_density) converges to corners and edges that carry
// state-dependent phases (signs for real tensors); the phase differs from
// window to window and, for complex tensors, rotates from sweep to sweep.
// Measurements cancel it; the real-positivity guards of
// Check_Convergence_CTM_RDM and the Hermitian/positive metric of the full
// update do not. Section 2.5 (revision 2) adds a phase_invariant flag to the
// convergence check / CTM driver, and makes build_full_update_environment
// divide N by the phase of its own window norm sum_x N(x) I_wrap(x),
// reporting the removed phase in full_update_environment::phase. Nothing
// rotates the environment itself.
//
// Included into the test_fermion_layer TU after fermion/full_update_realctm.cpp
// and reuses its fixtures without re-declaring them: fur_env (the converged
// real-CTM environment record), fur_window / fur_make_window / fur_window_env
// (the ten window tensors of section 3.3, the QR factors and the state's own
// block X), fur_closed (the measurement path), fur_build_N, fur_mask_ab,
// fur_quad, fur_hop_gate_plain, fur_random_even_block, fur_lstsq_deviation,
// fur_rel_diff, fur_assemble, fur_require_assembly, and the earlier helpers
// fue_internal_identity, fue_sum_product, fue_pseudo_site, fue_close_halves,
// fue_random_even_plain, fue_check_rel, fue_make_case, fue_full_geom,
// fub_restrict_bond, fub_gate_plain, fub_exact_params, fub_check_output_parity,
// fub_require_definite_environment, fg_ledger_mixed, fg_max_abs_entry,
// ib_hop_plain, ib_index_string.
//
// Truth sources (never the machinery under test):
//   - T8-i:   the same call sequence on the same environments with and
//             without a phase; the requirement is agreement, plus the
//             documented first-call / invalid-environment behaviour.
//   - T8-ii:  the same window folded from three environments that differ
//             only by a phase on C1 or on one edge; the requirement is that
//             N, N_plain and forbidden_ratio agree and `phase` follows the
//             injected factor. The value of `phase` itself is pinned by the
//             closed measurement path on the pseudo-sites (contract 3.2,
//             T2-i note: closed = phase * sum N O).
//   - T8-iv:  the exactly contracted patch of fue_make_case (the T3-i input)
//             with C1 rotated by hand; the two outputs must be the same pair
//             state up to scale and phase.
//   - T2-vi complex: the measurement path on Tn1(Y), Tn2(Y) against the
//             plain quadratic form of N_plain, as in the real T2-vi, with
//             |norm| and the ratio <G> compared because the measurement path
//             carries the window's phase and N_plain does not.
//
// Premises asserted (the tests are meaningless without them):
//   - every CTM run converged (count < Max_CTM_Iteration);
//   - the folded metric has no forbidden block (forbidden_ratio <= 1e-10);
//   - every ledger that is supposed to be mixed is mixed;
//   - the reference series of T8-i really is a "slightly different"
//     environment: its second call reports a finite positive distance, and
//     the two epsilons used bracket it (so both return values occur);
//   - for T8-ii the window norm of the raw environment is finite and
//     nonzero, and the injected phases are not 1;
//   - for T8-iv the metric of the unmodified input is positive definite and
//     the gate really moves the state, so "(a) == (b)" is not "both equal
//     the input".
//
// Element comparisons go through the GLOBAL index (contract section 0):
// mptensor evaluates transposes lazily, so local index n of two tensors of
// the same shape need not point at the same element.

#include <complex>
#include <limits>

#include "../../src/iTPS/core/contract_density_ctm.hpp"
#include "../../src/iTPS/core/contract_itps_ctm.hpp"
#include "../../src/iTPS/core/ctm.hpp"
#include "../../src/iTPS/core/full_update_fermion.hpp"

namespace {

// ---- scalar helpers -------------------------------------------------------

inline std::complex<double> fcp_to_complex(double v) { return {v, 0.0}; }
inline std::complex<double> fcp_to_complex(const std::complex<double>& v) {
  return v;
}

//! A complex scalar as the tensor's value type. For real tensors the
//! imaginary part must vanish (a real environment only admits the phases
//! +1 and -1); that is asserted, not truncated.
template <class tensor>
typename tensor::value_type fcp_from_complex(const std::complex<double>& v);

template <>
inline double fcp_from_complex<tenes::real_tensor>(
    const std::complex<double>& v) {
  REQUIRE(v.imag() == 0.0);
  return v.real();
}

template <>
inline std::complex<double> fcp_from_complex<tenes::complex_tensor>(
    const std::complex<double>& v) {
  return v;
}

template <class tensor>
constexpr bool fcp_is_complex() {
  return std::is_same<typename tensor::value_type, std::complex<double>>::value;
}

//! Multiply every element of a tensor by a scalar, through the global index.
template <class tensor>
void fcp_scale_inplace(tensor& a, const std::complex<double>& factor) {
  const typename tensor::value_type f = fcp_from_complex<tensor>(factor);
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    const mptensor::Index idx = a.global_index(n);
    typename tensor::value_type v;
    a.get_value(idx, v);
    a.set_value(idx, f * v);
  }
}

template <class tensor>
void fcp_scale_all(std::vector<tensor>& v, const std::complex<double>& factor) {
  for (auto& t : v) {
    fcp_scale_inplace(t, factor);
  }
}

//! Exact elementwise identity (the "byte unchanged" judgment), through the
//! global index. NaN never compares equal, which is the right answer here.
template <class tensor>
bool fcp_identical(const tensor& a, const tensor& b) {
  if (!(a.shape() == b.shape())) {
    return false;
  }
  int same = 1;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    const mptensor::Index idx = a.global_index(n);
    typename tensor::value_type va, vb;
    a.get_value(idx, va);
    b.get_value(idx, vb);
    if (!(va == vb)) {
      same = 0;
    }
  }
  tenes::allreduce_sum(same, a.get_comm());
  return same == a.get_comm_size();
}

template <class tensor>
bool fcp_identical_all(const std::vector<tensor>& a,
                       const std::vector<tensor>& b) {
  if (a.size() != b.size()) {
    return false;
  }
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (!fcp_identical(a[i], b[i])) {
      return false;
    }
  }
  return true;
}

//! max |got - factor * ref| / max |ref|, through the global index. Pure
//! relative (no max(1, .) floor): the metric's scale is not 1 in general.
template <class tensor>
double fcp_scaled_rel_diff(const tensor& got, const tensor& ref,
                           const std::complex<double>& factor) {
  REQUIRE(got.shape() == ref.shape());
  const typename tensor::value_type f = fcp_from_complex<tensor>(factor);
  const double scale = fg_max_abs_entry(ref);
  REQUIRE(scale > 0.0);
  double dev = 0.0;
  mptensor::Index worst;
  for (std::size_t n = 0; n < ref.local_size(); ++n) {
    const mptensor::Index idx = ref.global_index(n);
    typename tensor::value_type g, r;
    got.get_value(idx, g);
    ref.get_value(idx, r);
    const double d = std::abs(g - f * r);
    if (d > dev) {
      dev = d;
      worst = idx;
    }
  }
  INFO("[scaled_rel_diff] factor=" << factor << " max|got - f ref|=" << dev
                                   << " max|ref|=" << scale << " at "
                                   << ib_index_string(worst));
  return dev / scale;
}

//! Hermiticity residual of a rank-4 metric N(a, b, a', b'):
//! max |N(a,b,a',b') - conj(N(a',b',a,b))| / max |N|.
template <class tensor>
double fcp_hermiticity_residual(const tensor& N) {
  REQUIRE(N.rank() == 4);
  const double scale = fg_max_abs_entry(N);
  REQUIRE(scale > 0.0);
  double dev = 0.0;
  for (std::size_t n = 0; n < N.local_size(); ++n) {
    const mptensor::Index idx = N.global_index(n);
    const mptensor::Index jdx(idx[2], idx[3], idx[0], idx[1]);
    typename tensor::value_type a, b;
    N.get_value(idx, a);
    N.get_value(jdx, b);
    dev = std::max(dev,
                   std::abs(fcp_to_complex(a) - std::conj(fcp_to_complex(b))));
  }
  return dev / scale;
}

//! Spectrum bounds of a rank-4 metric seen as a matrix (a b) x (a' b').
template <class tensor>
void fcp_spectrum(const tensor& N, double& wmin, double& wmax) {
  std::vector<double> w;
  REQUIRE(mptensor::eigh(N, mptensor::Axes(0, 1), mptensor::Axes(2, 3), w) ==
          0);
  REQUIRE(!w.empty());
  wmin = w[0];
  wmax = w[0];
  for (const double v : w) {
    wmin = std::min(wmin, v);
    wmax = std::max(wmax, v);
  }
}

//! Deterministic draw, distinct from the other fixtures' generators.
inline double fcp_draw(int seed, const mptensor::Index& idx, std::size_t rank) {
  const double w[6] = {2.3, 4.1, 6.7, 8.9, 10.3, 12.7};
  double x = 0.61 * (seed + 3) + 0.017 * seed * seed;
  for (std::size_t ax = 0; ax < rank; ++ax) {
    x += w[ax] * static_cast<double>(idx[ax]);
  }
  return 0.53 * std::sin(x) + 0.29 * std::cos(1.31 * x + 0.19 * seed);
}

inline void fcp_set(tenes::real_tensor& t, const mptensor::Index& idx, int seed,
                    std::size_t rank) {
  t.set_value(idx, fcp_draw(seed, idx, rank));
}

inline void fcp_set(tenes::complex_tensor& t, const mptensor::Index& idx,
                    int seed, std::size_t rank) {
  t.set_value(idx, std::complex<double>(fcp_draw(seed, idx, rank),
                                        fcp_draw(seed + 700, idx, rank)));
}

// ---- one-site norms --------------------------------------------------------

//! Trace of a rank-2 (distributed) RDM: the diagonal is gathered with
//! get_value and allreduce_sum, as contract section 2.5 prescribes.
template <class tensor>
std::complex<double> fcp_trace(const tensor& rdm) {
  using value_type = typename tensor::value_type;
  REQUIRE(rdm.rank() == 2);
  const std::size_t d0 = rdm.shape()[0];
  REQUIRE(rdm.shape()[1] == d0);
  std::vector<value_type> diag(d0, value_type(0.0));
  for (std::size_t n = 0; n < rdm.local_size(); ++n) {
    const mptensor::Index idx = rdm.global_index(n);
    if (idx[0] == idx[1]) {
      value_type v;
      rdm.get_value(idx, v);
      diag[idx[0]] = v;
    }
  }
  tenes::allreduce_sum(diag, rdm.get_comm());
  std::complex<double> t = 0.0;
  for (const value_type& v : diag) {
    t += fcp_to_complex(v);
  }
  return t;
}

//! One-site norm of every site: the trace of the one-site RDM in the fold
//! (is_density) or the bosonic convention.
template <class tensor>
std::vector<std::complex<double>> fcp_one_site_norms(const fur_env<tensor>& e,
                                                     bool is_density) {
  std::vector<std::complex<double>> norms;
  for (int i = 0; i < e.lattice.N_UNIT; ++i) {
    const tensor rdm =
        is_density ? tenes::itps::core::Contract_one_site_RDM_density_CTM(
                         e.C1[i], e.C2[i], e.C3[i], e.C4[i], e.eTt[i], e.eTr[i],
                         e.eTb[i], e.eTl[i], e.reduced[i])
                   : tenes::itps::core::Contract_one_site_RDM_iTPS_CTM(
                         e.C1[i], e.C2[i], e.C3[i], e.C4[i], e.eTt[i], e.eTr[i],
                         e.eTb[i], e.eTl[i], e.Tn[i]);
    norms.push_back(fcp_trace(rdm));
  }
  return norms;
}

//! Rotate every C1[i] by the conjugate phase of site 0's one-site norm and
//! assert that afterwards EVERY site's norm is real and positive. This is a
//! fixture normalisation for T8-i only: it hands the convergence check a
//! series whose reference phase is known. (The one-site norms of all sites
//! share one phase; two-site windows do not, which is why nothing in the
//! production code rotates the environment.)
template <class tensor>
void fcp_normalize_phase_by_hand(fur_env<tensor>& e, bool is_density,
                                 const std::string& label) {
  std::vector<std::complex<double>> norms = fcp_one_site_norms(e, is_density);
  REQUIRE(!norms.empty());
  REQUIRE(std::isfinite(std::abs(norms[0])));
  REQUIRE(std::abs(norms[0]) > 0.0);
  std::complex<double> phase = norms[0] / std::abs(norms[0]);
  if (!fcp_is_complex<tensor>()) {
    phase = std::complex<double>(phase.real() < 0.0 ? -1.0 : 1.0, 0.0);
  }
  INFO(label << " [fixture: one-site phase of site 0 = " << phase << "]");
  fcp_scale_all(e.C1, std::conj(phase));
  norms = fcp_one_site_norms(e, is_density);
  for (std::size_t i = 0; i < norms.size(); ++i) {
    INFO(label << " [fixture: normalised one-site norm " << i << " = "
               << norms[i] << "]");
    REQUIRE(std::isfinite(std::abs(norms[i])));
    REQUIRE(norms[i].real() > 0.0);
    REQUIRE(std::abs(norms[i].imag()) <= 1.0e-8 * std::abs(norms[i]));
  }
}

// ---- environments ----------------------------------------------------------

//! Parity-even random lx x ly unit cell (d = 2, or a checkerboard of d = 2
//! and d = 4 when `checkerboard_d` is set - that makes the two internal
//! QR legs of every nearest-neighbour window differ, nA != nB; every virtual
//! ledger even_first_parity(D), which is mixed for D >= 2) and its converged
//! fold CTM environment through the NEW driver
//! Calc_CTM_Environment_density(..., initialize = true, phase_invariant =
//! true). Premise asserted: the CTM converged before Max_CTM_Iteration.
template <class tensor>
fur_env<tensor> fcp_make_fold_env(int seed, int chi, int D, int lx, int ly,
                                  const char* type_name, int perturb_seed = 0,
                                  double perturb = 0.0,
                                  bool checkerboard_d = false) {
  fur_env<tensor> e;
  e.lattice = tenes::SquareLattice(lx, ly);
  e.label =
      std::string("ctm-phase fold ") + type_name +
      " seed=" + std::to_string(seed) + " chi=" + std::to_string(chi) +
      " D=" + std::to_string(D) + " cell=" + std::to_string(lx) + "x" +
      std::to_string(ly) +
      (perturb != 0.0 ? " perturbed(" + std::to_string(perturb) + ")" : "") +
      (checkerboard_d ? " d=2/4 checkerboard" : "");
  INFO(e.label);
  const int N = e.lattice.N_UNIT;
  REQUIRE(N == lx * ly);
  std::vector<int> phys_dim(N, 2);
  for (int s = 0; s < N; ++s) {
    if (checkerboard_d && (e.lattice.x(s) + e.lattice.y(s)) % 2 == 1) {
      phys_dim[s] = 4;
    }
    e.lattice.physical_dims[s] = phys_dim[s];
    e.lattice.virtual_dims[s] = {D, D, D, D};
    e.lattice.initial_dirs[s] = {0.0};
    e.lattice.noises[s] = 0.0;
  }

  const fgf::parity_vector vp = fgf::even_first_parity(D);
  REQUIRE(fg_ledger_mixed(vp));
  e.finfo.enabled = true;
  e.finfo.phys.clear();
  e.params.phys_parity.clear();
  for (int s = 0; s < N; ++s) {
    e.finfo.phys.push_back(ib_phys_parity(phys_dim[s]));
    e.params.phys_parity.push_back(ib_phys_parity(phys_dim[s]));
  }
  e.finfo.virt.assign(N, std::array<fgf::parity_vector, 4>{vp, vp, vp, vp});
  fgf::validate_neighbor_consistency(e.finfo, e.lattice);

  for (int s = 0; s < N; ++s) {
    const fgf::leg_parities lp = fgf::Tn_parity(e.finfo, s);
    tensor t(mptensor::Shape(D, D, D, D, phys_dim[s]));
    for (std::size_t n = 0; n < t.local_size(); ++n) {
      const mptensor::Index idx = t.global_index(n);
      if (fgf::count_odd(lp, idx) % 2 == 0) {
        fcp_set(t, idx, seed + 7 * s, 5);
        if (perturb != 0.0) {
          typename tensor::value_type v;
          t.get_value(idx, v);
          t.set_value(
              idx,
              v * typename tensor::value_type(
                      1.0 + perturb * fcp_draw(perturb_seed + 3 * s, idx, 5)));
        }
      }
    }
    const fg_ftensor<tensor> wrapped{t, lp};
    REQUIRE(fgf::parity_violation(wrapped) == 0.0);
    REQUIRE(fgf::max_abs(wrapped) > 0.0);
    e.Tn.push_back(t);
  }

  e.params.print_level = tenes::PrintLevel::none;
  e.params.fermion = true;
  e.params.D = D;
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
      e.lattice, true, true);
  INFO(e.label << " [premise: CTM sweeps = " << e.ctm_count << " of at most "
               << e.params.Max_CTM_Iteration << "]");
  REQUIRE(e.ctm_count < e.params.Max_CTM_Iteration);
  return e;
}

//! Plain random 2x2 bosonic unit cell (D = 2, d = 2) and its converged
//! bosonic CTM environment (Calc_CTM_Environment). Used by T8-i for the
//! is_density = false branch of Check_Convergence_CTM_RDM. e.reduced stays
//! empty: the bosonic one-site RDM is contracted from Tn itself.
template <class tensor>
fur_env<tensor> fcp_make_boson_env(int seed, int chi, const char* type_name,
                                   int perturb_seed = 0, double perturb = 0.0) {
  fur_env<tensor> e;
  e.lattice = tenes::SquareLattice(2, 2);
  e.label =
      std::string("ctm-phase boson ") + type_name +
      " seed=" + std::to_string(seed) + " chi=" + std::to_string(chi) +
      (perturb != 0.0 ? " perturbed(" + std::to_string(perturb) + ")" : "");
  INFO(e.label);
  const int N = e.lattice.N_UNIT;
  REQUIRE(N == 4);
  const int D = 2;
  for (int s = 0; s < N; ++s) {
    e.lattice.physical_dims[s] = 2;
    e.lattice.virtual_dims[s] = {D, D, D, D};
    e.lattice.initial_dirs[s] = {0.0};
    e.lattice.noises[s] = 0.0;
  }
  for (int s = 0; s < N; ++s) {
    tensor t(mptensor::Shape(D, D, D, D, 2));
    for (std::size_t n = 0; n < t.local_size(); ++n) {
      const mptensor::Index idx = t.global_index(n);
      fcp_set(t, idx, seed + 11 * s, 5);
      if (perturb != 0.0) {
        typename tensor::value_type v;
        t.get_value(idx, v);
        t.set_value(idx, v * typename tensor::value_type(
                                 1.0 + perturb * fcp_draw(perturb_seed + 5 * s,
                                                          idx, 5)));
      }
    }
    REQUIRE(fg_max_abs_entry(t) > 0.0);
    e.Tn.push_back(t);
  }
  e.params.print_level = tenes::PrintLevel::none;
  e.params.fermion = false;
  e.params.D = D;
  e.params.CHI = chi;
  e.params.Use_RSVD = false;
  e.params.Max_CTM_Iteration = 2000;
  e.params.CTM_Convergence_Epsilon = 1.0e-12;
  e.C1.assign(N, tensor());
  e.C2.assign(N, tensor());
  e.C3.assign(N, tensor());
  e.C4.assign(N, tensor());
  e.eTt.assign(N, tensor());
  e.eTr.assign(N, tensor());
  e.eTb.assign(N, tensor());
  e.eTl.assign(N, tensor());
  e.ctm_count = tenes::itps::core::Calc_CTM_Environment(
      e.C1, e.C2, e.C3, e.C4, e.eTt, e.eTr, e.eTb, e.eTl, e.Tn, e.params,
      e.lattice);
  INFO(e.label << " [premise: CTM sweeps = " << e.ctm_count << " of at most "
               << e.params.Max_CTM_Iteration << "]");
  REQUIRE(e.ctm_count < e.params.Max_CTM_Iteration);
  return e;
}

}  // namespace

// ---- T8-i: phase invariance of the convergence check ----------------------

namespace {

//! Two consecutive calls of Check_Convergence_CTM_RDM with a fresh history:
//! first on E0, then on E1. The history (rdm_old, has_rdm_old) is private to
//! the series, as the contract requires.
struct fcp_series {
  bool first = true;
  double first_dist = 0.0;
  bool second = true;
  double second_dist = 0.0;
};

template <class tensor>
fcp_series fcp_run_series(const fur_env<tensor>& E0, const fur_env<tensor>& E1,
                          bool is_density, double epsilon,
                          bool phase_invariant) {
  using value_type = typename tensor::value_type;
  std::vector<tenes::small_tensor<value_type>> rdm_old;
  bool has_rdm_old = false;
  fcp_series r;
  double dist = 0.0;
  r.first = tenes::itps::core::Check_Convergence_CTM_RDM(
      E0.C1, E0.C2, E0.C3, E0.C4, E0.eTt, E0.eTr, E0.eTb, E0.eTl,
      is_density ? E0.reduced : E0.Tn, E0.lattice, rdm_old, has_rdm_old,
      epsilon, is_density, dist, phase_invariant);
  r.first_dist = dist;
  dist = 0.0;
  r.second = tenes::itps::core::Check_Convergence_CTM_RDM(
      E1.C1, E1.C2, E1.C3, E1.C4, E1.eTt, E1.eTr, E1.eTb, E1.eTl,
      is_density ? E1.reduced : E1.Tn, E1.lattice, rdm_old, has_rdm_old,
      epsilon, is_density, dist, phase_invariant);
  r.second_dist = dist;
  return r;
}

//! A copy of the environment with every C1[i] multiplied by the phase.
template <class tensor>
fur_env<tensor> fcp_with_phase(const fur_env<tensor>& e,
                               const std::complex<double>& phase) {
  fur_env<tensor> r = e;
  fcp_scale_all(r.C1, phase);
  return r;
}

using fcp_phase_pair = std::pair<std::complex<double>, std::complex<double>>;

/*! The T8-i requirement on one pair (E0, E1) of positive-phase environments:
 *  the reference series E0 -> E1 and the transformed series
 *  p0 E0 -> p1 E1 (p0 != p1) agree call by call with phase_invariant = true.
 */
template <class tensor>
void fcp_check_t8i_pair(const fur_env<tensor>& E0, const fur_env<tensor>& E1,
                        bool is_density,
                        const std::vector<fcp_phase_pair>& phases,
                        const std::string& label) {
  INFO(label);
  // Bracket the reference distance: one epsilon above it (converged) and
  // one below it (not converged), so that BOTH return values are compared.
  // The reference distance itself is only used to place the brackets; the
  // requirement is agreement between the series.
  const fcp_series probe = fcp_run_series(E0, E1, is_density, 1.0, true);
  INFO(label << " [premise: reference series] first=" << probe.first
             << " first_dist=" << probe.first_dist << " second=" << probe.second
             << " second_dist=" << probe.second_dist);
  REQUIRE(probe.first == false);
  REQUIRE(std::isnan(probe.first_dist));
  REQUIRE(std::isfinite(probe.second_dist));
  REQUIRE(probe.second_dist > 0.0);
  const double eps_hi = 10.0 * probe.second_dist;
  const double eps_lo = 0.1 * probe.second_dist;

  for (const double eps : {eps_hi, eps_lo}) {
    const fcp_series ref = fcp_run_series(E0, E1, is_density, eps, true);
    INFO(label << " eps=" << eps << " [reference] second=" << ref.second
               << " dist=" << ref.second_dist);
    REQUIRE(ref.first == false);
    REQUIRE(std::isnan(ref.first_dist));
    REQUIRE(!std::isnan(ref.second_dist));
    // The brackets really bracket: the converged verdict at eps_hi and the
    // unconverged one at eps_lo. Without this, "the verdicts agree" could be
    // satisfied by an implementation that always says false.
    CHECK(ref.second == (eps == eps_hi));

    for (const auto& pp : phases) {
      const std::complex<double>& p0 = pp.first;
      const std::complex<double>& p1 = pp.second;
      REQUIRE(std::abs(std::abs(p0) - 1.0) <= 1.0e-15);
      REQUIRE(std::abs(std::abs(p1) - 1.0) <= 1.0e-15);
      REQUIRE(std::abs(p0 - p1) > 1.0e-3);
      const fur_env<tensor> T0 = fcp_with_phase(E0, p0);
      const fur_env<tensor> T1 = fcp_with_phase(E1, p1);
      const fcp_series tr = fcp_run_series(T0, T1, is_density, eps, true);
      INFO(label << " eps=" << eps << " phases=(" << p0 << ", " << p1
                 << ") [transformed] first=" << tr.first
                 << " first_dist=" << tr.first_dist << " second=" << tr.second
                 << " dist=" << tr.second_dist
                 << " reference dist=" << ref.second_dist);
      CHECK(tr.first == false);
      CHECK(std::isnan(tr.first_dist));
      CHECK(!std::isnan(tr.second_dist));
      CHECK(tr.second == ref.second);
      CHECK(std::abs(tr.second_dist - ref.second_dist) <=
            1.0e-12 *
                std::max(std::abs(tr.second_dist), std::abs(ref.second_dist)));
    }
  }
}

//! Contract T8-i, additional requirements: an all-zero C1, a NaN element and
//! an Inf element each leave the second call unconverged with a NaN
//! distance, whatever phase_invariant is.
template <class tensor>
void fcp_check_t8i_invalid(const fur_env<tensor>& E0, bool is_density,
                           const std::string& label) {
  INFO(label);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<std::pair<std::string, fur_env<tensor>>> bad;
  {
    fur_env<tensor> z = E0;
    fcp_scale_all(z.C1, 0.0);
    for (const auto& c : z.C1) {
      REQUIRE(fg_max_abs_entry(c) == 0.0);
    }
    bad.emplace_back("all C1 zero", z);
  }
  {
    fur_env<tensor> n = E0;
    const mptensor::Index origin(std::vector<std::size_t>(n.C1[0].rank(), 0));
    n.C1[0].set_value(origin, typename tensor::value_type(nan));
    bad.emplace_back("NaN in C1[0]", n);
  }
  {
    fur_env<tensor> n = E0;
    const mptensor::Index origin(std::vector<std::size_t>(n.eTr[1].rank(), 0));
    n.eTr[1].set_value(origin, typename tensor::value_type(inf));
    bad.emplace_back("Inf in eTr[1]", n);
  }
  for (const auto& item : bad) {
    for (const bool phase_invariant : {false, true}) {
      const fcp_series s = fcp_run_series(item.second, item.second, is_density,
                                          1.0e6, phase_invariant);
      INFO(label << " [" << item.first
                 << "] phase_invariant=" << phase_invariant
                 << " first=" << s.first << " first_dist=" << s.first_dist
                 << " second=" << s.second << " second_dist=" << s.second_dist);
      CHECK(s.first == false);
      CHECK(std::isnan(s.first_dist));
      CHECK(s.second == false);
      CHECK(std::isnan(s.second_dist));
    }
  }
}

//! Contract T8-i, last requirement: with phase_invariant = false (the
//! default) a real series with a -1 phase stays unconverged on the second
//! call too - the pre-existing behaviour is kept for the bosonic caller.
template <class tensor>
void fcp_check_t8i_default_real_sign(const fur_env<tensor>& E0,
                                     const fur_env<tensor>& E1, bool is_density,
                                     const std::string& label) {
  INFO(label);
  const fur_env<tensor> M0 = fcp_with_phase(E0, -1.0);
  const fur_env<tensor> M1 = fcp_with_phase(E1, -1.0);
  const fcp_series neg = fcp_run_series(M0, M1, is_density, 1.0e6, false);
  INFO(label << " [phase_invariant=false, both -1] first=" << neg.first
             << " second=" << neg.second << " dist=" << neg.second_dist);
  CHECK(neg.first == false);
  CHECK(std::isnan(neg.first_dist));
  CHECK(neg.second == false);
  CHECK(std::isnan(neg.second_dist));
  // ... and the positive series is accepted by the default path, so the
  // negative verdict above is about the sign and not about the environment.
  const fcp_series pos = fcp_run_series(E0, E1, is_density, 1.0e6, false);
  INFO(label << " [phase_invariant=false, both +1] first=" << pos.first
             << " second=" << pos.second << " dist=" << pos.second_dist);
  CHECK(pos.first == false);
  CHECK(!std::isnan(pos.second_dist));
  CHECK(pos.second == true);
}

template <class tensor>
void fcp_run_t8i(bool is_density, int seed, const char* type_name) {
  const std::string label =
      std::string("T8-i ") + type_name + (is_density ? " fold" : " boson");
  INFO(label);
  fur_env<tensor> E0 =
      is_density ? fcp_make_fold_env<tensor>(seed, 8, 2, 2, 2, type_name)
                 : fcp_make_boson_env<tensor>(seed, 8, type_name);
  // E1: the converged environment of a slightly different state (every Tn
  // element scaled by 1 + 0.01 * random), so that it is "a little different"
  // from E0 in a way a CTM iteration could produce - not an arbitrary
  // elementwise perturbation, which for complex tensors would rotate the
  // one-site trace by O(0.01) and be rejected for a reason unrelated to
  // the phase.
  fur_env<tensor> E1 =
      is_density
          ? fcp_make_fold_env<tensor>(seed, 8, 2, 2, 2, type_name, seed + 50,
                                      0.01)
          : fcp_make_boson_env<tensor>(seed, 8, type_name, seed + 50, 0.01);
  fcp_normalize_phase_by_hand(E0, is_density, label + " E0");
  fcp_normalize_phase_by_hand(E1, is_density, label + " E1");

  // Real tensors admit only the signs; complex ones two generic angles
  // (0.7 and 2.9 rad, both orders) and the signs.
  const std::complex<double> minus(-1.0, 0.0);
  const std::complex<double> plus(1.0, 0.0);
  std::vector<fcp_phase_pair> phases{{minus, plus}, {plus, minus}};
  if (fcp_is_complex<tensor>()) {
    phases.emplace_back(std::polar(1.0, 0.7), std::polar(1.0, 2.9));
    phases.emplace_back(std::polar(1.0, 2.9), std::polar(1.0, 0.7));
  }
  fcp_check_t8i_pair(E0, E1, is_density, phases, label);
  fcp_check_t8i_invalid(E0, is_density, label);
  if (!fcp_is_complex<tensor>()) {
    fcp_check_t8i_default_real_sign(E0, E1, is_density, label);
  }
}

}  // namespace

// One TEST_CASE per environment kind, so that a premise failure in one
// (e.g. the CTM not converging) does not hide the verdict on the others.
TEST_CASE(
    "ctm phase T8-i (fold real): Check_Convergence_CTM_RDM with "
    "phase_invariant is blind to a sign on C1 and rejects invalid input") {
  fcp_run_t8i<tenes::real_tensor>(true, 7100, "real");
}

TEST_CASE(
    "ctm phase T8-i (fold complex): Check_Convergence_CTM_RDM with "
    "phase_invariant is blind to a phase on C1 and rejects invalid input") {
  fcp_run_t8i<tenes::complex_tensor>(true, 7200, "complex");
}

TEST_CASE(
    "ctm phase T8-i (boson real): Check_Convergence_CTM_RDM with "
    "phase_invariant is blind to a sign on C1 for a bosonic Tn") {
  fcp_run_t8i<tenes::real_tensor>(false, 7300, "real");
}

TEST_CASE(
    "ctm phase T8-i (boson complex): Check_Convergence_CTM_RDM with "
    "phase_invariant is blind to a phase on C1 for a bosonic Tn") {
  fcp_run_t8i<tenes::complex_tensor>(false, 7400, "complex");
}

// ---- a window in a phased environment -------------------------------------

namespace {

/*! fur_make_window without its last premise (a real positive measurement-
 *  path norm). In a real CTM environment every window carries its own
 *  phase, so that premise cannot be met and the norm is only required to be
 *  finite and nonzero here. Everything else - the window tensors of section
 *  3.3, the QR axes of T2-vi, the reconstruction premise and the block X -
 *  is fur_make_window's, obtained by calling it on a copy of the environment
 *  whose C1[s1] has been rotated by hand so that THIS window's norm is real
 *  positive; the window's env slots are then replaced by the environment
 *  actually under test. (C1[s1] enters the window exactly once, so the
 *  rotation is a pure phase on the closed value.)
 */
template <class tensor>
fur_window<tensor> fcp_make_window(const fur_env<tensor>& e, int s1, char dir,
                                   const std::string& tag) {
  const int s2 = dir == 'h' ? e.lattice.right(s1) : e.lattice.bottom(s1);
  // The window's own norm in the environment under test.
  fur_env<tensor> helper = e;
  {
    const fg_ftensor<tensor> Tn1 = fgf::wrap_Tn(e.Tn[s1], e.finfo, s1);
    const fg_ftensor<tensor> Tn2 = fgf::wrap_Tn(e.Tn[s2], e.finfo, s2);
    const fg_ftensor<tensor> I = fgf::wrap_twosite_gate(
        fue_internal_identity<tensor>(Tn1.parity[4], Tn2.parity[4]),
        Tn1.parity[4], Tn2.parity[4]);
    const auto dir_e = dir == 'h' ? fgf::reduced_pair_direction::horizontal
                                  : fgf::reduced_pair_direction::vertical;
    const std::vector<tensor> env = fur_window_env(e, s1, s2, dir);
    const auto halves = fgf::build_reduced_pair_halves(Tn1, Tn2, I, dir_e);
    const auto norm = fgf::contract_reduced_pair_halves_density_CTM(
        env[0], env[1], env[2], env[3], env[4], env[5], env[6], env[7], env[8],
        env[9], halves);
    const std::complex<double> n = fcp_to_complex(norm);
    INFO(tag << " [window norm in the environment under test = " << n
             << ", arg " << std::arg(n) << "]");
    REQUIRE(std::isfinite(std::abs(n)));
    REQUIRE(std::abs(n) > 0.0);
    std::complex<double> phase = n / std::abs(n);
    if (!fcp_is_complex<tensor>()) {
      phase = std::complex<double>(phase.real() < 0.0 ? -1.0 : 1.0, 0.0);
    }
    // The helper copy: C1[s1] rotated so that fur_make_window's premise
    // holds. Only used to run the QR / reconstruction / X machinery.
    fcp_scale_inplace(helper.C1[s1], std::conj(phase));
  }
  fur_window<tensor> w = fur_make_window(helper, s1, dir, tag);
  REQUIRE(w.s2 == s2);
  w.env = fur_window_env(e, s1, s2, dir);
  return w;
}

//! The window norm of contract 2.5: sum_x N.t(x) I_wrap.t(x), the plain
//! product sum with the wrapped identity on the two internal legs.
template <class tensor>
std::complex<double> fcp_window_norm(
    const fur_window<tensor>& w,
    const fue::full_update_environment<tensor>& N) {
  const fg_ftensor<tensor> I = fgf::wrap_twosite_gate(
      fue_internal_identity<tensor>(w.pa, w.pb), w.pa, w.pb);
  REQUIRE(I.t.shape() == N.N.t.shape());
  return fcp_to_complex(fue_sum_product(N.N.t, I.t));
}

//! The closed measurement path on the pseudo-sites of T2-i with the
//! window's own ten environment tensors: what sum N O would be WITHOUT the
//! phase normalisation.
template <class tensor>
std::complex<double> fcp_closed_pseudo(const fur_window<tensor>& w,
                                       const fg_ftensor<tensor>& O) {
  const fg_ftensor<tensor> QAp = fue_pseudo_site(w.QA, w.dir == 'h' ? 2 : 3);
  const fg_ftensor<tensor> QBp = fue_pseudo_site(w.QB, w.dir == 'h' ? 0 : 1);
  const auto halves = fgf::build_reduced_pair_halves(QAp, QBp, O, w.dir_e);
  return fcp_to_complex(fue_close_halves(w.env, halves));
}

}  // namespace

// ---- T8-ii: phase invariance of N
// --------------------------------------------

namespace {

/*! Contract T8-ii on one window. Three environments that differ only by a
 *  phase are folded: (a) the CTM output, (b) every C1[i] times `factor`,
 *  (c) one edge of the window (eT1 = eTt[s1]) times `factor`. The window
 *  contracts each of C1[s1] and eTt[s1] exactly once, so (b) and (c) both
 *  multiply the unnormalised N by `factor` and the normalised N not at all.
 */
template <class tensor>
void fcp_run_t8ii_window(const fur_env<tensor>& e, int s1, char dir,
                         const std::complex<double>& factor,
                         bool require_unequal_internal,
                         const std::string& label) {
  const fur_window<tensor> wa = fcp_make_window(e, s1, dir, label);
  INFO(wa.label);
  REQUIRE(std::abs(std::abs(factor) - 1.0) <= 1.0e-15);
  REQUIRE(std::abs(factor - 1.0) > 0.1);
  INFO(wa.label << " [internal legs] nA=" << wa.pa.size()
                << " nB=" << wa.pb.size());
  if (require_unequal_internal) {
    // Contract T8-ii item 6: at least one window with nA != nB, so that a
    // leg-order error of I_wrap cannot hide behind a square shape.
    REQUIRE(wa.pa.size() != wa.pb.size());
  }

  // (b): every C1 of the environment. (c): eT1 = eTt[s1] of the window.
  fur_window<tensor> wb = wa;
  {
    fur_env<tensor> eb = fcp_with_phase(e, factor);
    wb.env = fur_window_env(eb, wa.s1, wa.s2, dir);
  }
  fur_window<tensor> wc = wa;
  fcp_scale_inplace(wc.env[4], factor);
  // Premise (anti-vacuity): (b) and (c) really are different inputs.
  REQUIRE(!fcp_identical(wb.env[0], wa.env[0]));
  REQUIRE(!fcp_identical(wc.env[4], wa.env[4]));

  const fue::full_update_environment<tensor> Na = fur_build_N(wa);
  const fue::full_update_environment<tensor> Nb = fur_build_N(wb);
  const fue::full_update_environment<tensor> Nc = fur_build_N(wc);
  INFO(wa.label << " [phases] a=" << Na.phase << " b=" << Nb.phase
                << " c=" << Nc.phase << " factor=" << factor);
  INFO(wa.label << " [premise: forbidden_ratio] a=" << Na.forbidden_ratio
                << " b=" << Nb.forbidden_ratio << " c=" << Nc.forbidden_ratio);
  REQUIRE(Na.forbidden_ratio <= 1.0e-10);
  REQUIRE(Nb.forbidden_ratio <= 1.0e-10);
  REQUIRE(Nc.forbidden_ratio <= 1.0e-10);
  REQUIRE(Na.N.rank() == 4);
  REQUIRE(Na.N.parity == fgf::leg_parities{wa.pa, wa.pb, wa.pa, wa.pb});
  REQUIRE(fg_max_abs_entry(Na.N.t) > 0.0);

  // (1) N.t and N_plain are the same for the three inputs.
  for (const auto* other : {&Nb, &Nc}) {
    const char which = other == &Nb ? 'b' : 'c';
    INFO(wa.label << " [T8-ii 1: N of (" << which << ") vs (a)]");
    REQUIRE(other->N.parity == Na.N.parity);
    CHECK(fcp_scaled_rel_diff(other->N.t, Na.N.t, 1.0) <= 1.0e-12);
    CHECK(fcp_scaled_rel_diff(other->N_plain, Na.N_plain, 1.0) <= 1.0e-12);
  }

  // (2) the reported phase follows the injected factor; |phase| = 1.
  CHECK(std::abs(std::abs(Na.phase) - 1.0) <= 1.0e-12);
  CHECK(std::abs(Nb.phase - factor * Na.phase) <= 1.0e-12);
  CHECK(std::abs(Nc.phase - factor * Na.phase) <= 1.0e-12);
  if (!fcp_is_complex<tensor>()) {
    // A real environment admits only the signs.
    CHECK(Na.phase.imag() == 0.0);
    CHECK(Nb.phase.imag() == 0.0);
    CHECK(Nc.phase.imag() == 0.0);
  }
  // Anti-vacuity: with factor != 1 at least one of the three phases is not 1.
  REQUIRE((std::abs(Na.phase - 1.0) > 0.1 || std::abs(Nb.phase - 1.0) > 0.1));

  // (3) the window norm of the normalised N is real positive, and
  // (4) N_plain is Hermitian and positive semidefinite,
  // (5) forbidden_ratio agrees - all three inputs.
  for (const auto* env : {&Na, &Nb, &Nc}) {
    const char which = env == &Na ? 'a' : (env == &Nb ? 'b' : 'c');
    const std::string vlabel = wa.label + " [variant " + which + "]";
    const std::complex<double> norm = fcp_window_norm(wa, *env);
    INFO(vlabel << " [T8-ii 3] window norm = " << norm);
    REQUIRE(std::isfinite(std::abs(norm)));
    CHECK(norm.real() > 0.0);
    CHECK(std::abs(norm.imag()) <= 1.0e-12 * std::abs(norm));

    const double herm = fcp_hermiticity_residual(env->N_plain);
    double wmin = 0.0;
    double wmax = 0.0;
    fcp_spectrum(env->N_plain, wmin, wmax);
    INFO(vlabel << " [T8-ii 4] hermiticity residual=" << herm
                << " wmin=" << wmin << " wmax=" << wmax);
    CHECK(herm <= 1.0e-12);
    CHECK(wmax > 0.0);
    CHECK(wmin >= -1.0e-12 * wmax);

    // The ratio is roundoff-level itself (contract T2-ii: <= 1e-10, and in
    // practice 1e-13..1e-16), so "agrees" is an absolute statement here.
    INFO(vlabel << " [T8-ii 5] forbidden_ratio=" << env->forbidden_ratio
                << " (a)=" << Na.forbidden_ratio);
    CHECK(std::abs(env->forbidden_ratio - Na.forbidden_ratio) <= 1.0e-12);
  }

  // (6) independent verification (and contract 3.2, T2-i note): the closed
  // measurement path on the pseudo-sites carries the environment's phase,
  // the normalised N does not, so closed(O) == phase * sum_x N(x) O(x). With
  // O = I_wrap this is the unnormalised window norm against
  // contract_reduced_pair_halves_density_CTM; it pins the VALUE of `phase`
  // to an independent path, not just its transformation law, and the leg
  // order / wrap mask / plain product sum of I_wrap against the closed fold.
  // Also checked with a dense random even operator, on all three inputs.
  {
    const fg_ftensor<tensor> I = fgf::wrap_twosite_gate(
        fue_internal_identity<tensor>(wa.pa, wa.pb), wa.pa, wa.pb);
    const fg_ftensor<tensor> R = fgf::wrap_twosite_gate(
        fue_random_even_plain<tensor>(wa.pa, wa.pb, 8100 + s1 + (dir == 'h')),
        wa.pa, wa.pb);
    const std::vector<std::pair<const fur_window<tensor>*,
                                const fue::full_update_environment<tensor>*>>
        variants{{&wa, &Na}, {&wb, &Nb}, {&wc, &Nc}};
    for (const auto& v : variants) {
      for (const auto* O : {&I, &R}) {
        const std::complex<double> closed = fcp_closed_pseudo(*v.first, *O);
        const std::complex<double> open =
            fcp_to_complex(fue_sum_product(v.second->N.t, O->t));
        const std::complex<double> want = v.second->phase * open;
        INFO(wa.label << " [T2-i note: closed = phase * sum N O] O="
                      << (O == &I ? "identity" : "random")
                      << " closed=" << closed << " phase*open=" << want);
        REQUIRE(std::abs(closed) > 0.0);
        CHECK(std::abs(closed - want) <= 1.0e-12 * std::abs(closed));
      }
    }
  }
}

template <class tensor>
void fcp_run_t8ii_state(int seed, int chi, int D, int lx, int ly, int v_site,
                        const char* type_name, bool checkerboard_d = false) {
  const fur_env<tensor> e = fcp_make_fold_env<tensor>(
      seed, chi, D, lx, ly, type_name, 0, 0.0, checkerboard_d);
  const std::string label = "T8-ii " + e.label;
  INFO(label);
  const std::complex<double> factor = fcp_is_complex<tensor>()
                                          ? std::polar(1.0, 1.3)
                                          : std::complex<double>(-1.0, 0.0);
  fcp_run_t8ii_window(e, 0, 'h', factor, checkerboard_d, label);
  fcp_run_t8ii_window(e, v_site, 'v', factor, checkerboard_d, label);

  // (7) the environments the fold must refuse: an all-zero C1 (no window
  // norm to normalise with), and a single NaN or Inf element in one of the
  // window's C or eT tensors. Note on the contract's "outside the support of
  // I_wrap": the fold is a dense contraction, so one non-finite input
  // element reaches EVERY element of N (NaN * 0 = NaN, Inf * 0 = NaN), the
  // window norm included; no injection position keeps the norm finite. The
  // non-finite check on N.t and the norm check are therefore both tripped
  // here, and this test cannot tell them apart (reported to the author).
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (const auto& sw :
         std::vector<std::pair<int, char>>{{0, 'h'}, {v_site, 'v'}}) {
      // fcp_make_window needs a nonzero norm to build the helper; take the
      // window geometry from the real environment and swap the slots.
      const fur_window<tensor> w0 =
          fcp_make_window(e, sw.first, sw.second, label);
      {
        fur_env<tensor> zero = e;
        fcp_scale_all(zero.C1, 0.0);
        for (const auto& c : zero.C1) {
          REQUIRE(fg_max_abs_entry(c) == 0.0);
        }
        fur_window<tensor> w = w0;
        w.env = fur_window_env(zero, w.s1, w.s2, sw.second);
        INFO(w.label << " [T8-ii 7: all C1 zero]");
        CHECK_THROWS_AS(fur_build_N(w), std::runtime_error);
      }
      {
        // NaN in one element of C2 (window slot 1).
        fur_window<tensor> w = w0;
        const mptensor::Index origin(
            std::vector<std::size_t>(w.env[1].rank(), 0));
        w.env[1].set_value(origin, typename tensor::value_type(nan));
        INFO(w.label << " [T8-ii 7: NaN in C2]");
        CHECK_THROWS_AS(fur_build_N(w), std::runtime_error);
      }
      {
        // Inf in one element of eT4 (window slot 7).
        fur_window<tensor> w = w0;
        const mptensor::Index origin(
            std::vector<std::size_t>(w.env[7].rank(), 0));
        w.env[7].set_value(origin, typename tensor::value_type(inf));
        INFO(w.label << " [T8-ii 7: Inf in eT4]");
        CHECK_THROWS_AS(fur_build_N(w), std::runtime_error);
      }
    }
  }
}

}  // namespace

TEST_CASE(
    "ctm phase T8-ii (a): on a real 2x2 D=3 fold environment N is invariant "
    "under a sign on C1 or on an edge and phase follows the sign") {
  fcp_run_t8ii_state<tenes::real_tensor>(6210, 12, 3, 2, 2, 3, "real");
}

TEST_CASE(
    "ctm phase T8-ii (b): on a complex 2x2 D=2 fold environment N is "
    "invariant under a phase on C1 or on an edge and is Hermitian") {
  fcp_run_t8ii_state<tenes::complex_tensor>(6300, 8, 2, 2, 2, 3, "complex");
}

TEST_CASE(
    "ctm phase T8-ii (c): on a complex 3x2 D=2 fold environment N is "
    "invariant under a phase on C1 or on an edge and is Hermitian") {
  fcp_run_t8ii_state<tenes::complex_tensor>(6410, 8, 2, 3, 2, 4, "complex");
}

TEST_CASE(
    "ctm phase T8-ii (d): with a d=2/4 checkerboard (nA != nB) the "
    "unnormalised window norm matches the closed measurement path") {
  fcp_run_t8ii_state<tenes::complex_tensor>(6500, 8, 2, 2, 2, 0, "complex",
                                            true);
}

// ---- T8-iv: the full-update kernel is blind to a phase on C1 --------------

namespace {

/*! The T3-i input (fue_make_case in the full geometry, bond-restricted
 *  sites, the a + b c^dag c gate, exact cutoffs), run twice: (a) as is and
 *  (b) with C1 multiplied by -1 (complex: e^{2i}). The window contracts C1
 *  exactly once, so the unnormalised metric picks up the factor; the kernel
 *  must normalise it away, throw nothing, and return the same pair state up
 *  to scale and phase.
 */
template <class tensor>
void fcp_run_t8iv_case(char dir, int d, int seed, const char* type_name) {
  const std::string tag = std::string("T8-iv ") + type_name;
  const fue_case<tensor> c = fue_make_case<tensor>(
      dir, fue_full_geom(dir), "eo", d, d, 2, seed, false, tag);
  INFO(c.label);
  // Premise: the unmodified metric is positive definite (otherwise "the same
  // pair state" is not an elementwise statement).
  fub_require_definite_environment(c, c.label);

  const int a_bond = dir == 'h' ? 2 : 3;
  const int b_bond = dir == 'h' ? 0 : 1;
  const fg_ftensor<tensor> Tn1 = fub_restrict_bond(c.sites[c.a], a_bond);
  const fg_ftensor<tensor> Tn2 = fub_restrict_bond(c.sites[c.b], b_bond);
  REQUIRE(fgf::max_abs(Tn1) > 0.0);
  REQUIRE(fgf::max_abs(Tn2) > 0.0);
  const tensor plain = fub_gate_plain<tensor>(d, 1.0, 0.8);
  const fg_ftensor<tensor> gate =
      fgf::wrap_twosite_gate(plain, Tn1.parity[4], Tn2.parity[4]);
  const tenes::itps::PEPS_Parameters params = fub_exact_params(true);
  const std::complex<double> factor = fcp_is_complex<tensor>()
                                          ? std::polar(1.0, 2.0)
                                          : std::complex<double>(-1.0, 0.0);

  auto run = [&](const std::vector<tensor>& env, fg_ftensor<tensor>& out1,
                 fg_ftensor<tensor>& out2) {
    tenes::itps::core::Full_update_bond_fermion(
        env[0], env[1], env[2], env[3], env[4], env[5], env[6], env[7], env[8],
        env[9], Tn1, Tn2, gate, c.dir_e, params, out1, out2);
  };

  // (a) as is.
  fg_ftensor<tensor> a1, a2;
  INFO(c.label << " [T8-iv a: unmodified]");
  REQUIRE_NOTHROW(run(c.env, a1, a2));
  fub_check_output_parity(c, a1, a2, Tn1, Tn2, c.label + " [T8-iv a]");
  const fg_ftensor<tensor> pair_a = fgf::build_pair_state(a1, a2, c.dir_e);
  REQUIRE(fgf::max_abs(pair_a) > 0.0);
  // Premise (anti-vacuity): the gate moved the state, so agreement between
  // (a) and (b) below is not "both still equal the input".
  {
    const fg_ftensor<tensor> pair_in = fgf::build_pair_state(Tn1, Tn2, c.dir_e);
    typename tensor::value_type cc;
    const double moved = fur_lstsq_deviation(pair_a.t, pair_in.t, cc);
    INFO(c.label << " [premise: the update moved the state] dev=" << moved);
    REQUIRE(moved > 1.0e-3);
  }

  // (b) C1 times the factor: no exception, the same pair state.
  std::vector<tensor> env = c.env;
  fcp_scale_inplace(env[0], factor);
  REQUIRE(!fcp_identical(env[0], c.env[0]));
  fg_ftensor<tensor> b1, b2;
  INFO(c.label << " [T8-iv b: C1 * " << factor << "]");
  REQUIRE_NOTHROW(run(env, b1, b2));
  fub_check_output_parity(c, b1, b2, Tn1, Tn2, c.label + " [T8-iv b]");
  const fg_ftensor<tensor> pair_b = fgf::build_pair_state(b1, b2, c.dir_e);
  REQUIRE(pair_b.shape() == pair_a.shape());
  REQUIRE(pair_b.parity == pair_a.parity);
  typename tensor::value_type cc;
  const double dev = fur_lstsq_deviation(pair_b.t, pair_a.t, cc);
  INFO(c.label << " [T8-iv] |pair_b - c pair_a| / max|pair_a| = " << dev
               << " (c = " << cc << ")");
  REQUIRE(std::abs(cc) > 0.0);
  CHECK(dev <= 1.0e-10);
}

}  // namespace

TEST_CASE(
    "ctm phase T8-iv: Full_update_bond_fermion returns the same pair state "
    "with or without a phase on C1 and throws nothing") {
  fcp_run_t8iv_case<tenes::real_tensor>('h', 2, 3100, "real d=2");
  fcp_run_t8iv_case<tenes::real_tensor>('v', 2, 3200, "real d=2");
  fcp_run_t8iv_case<tenes::complex_tensor>('h', 2, 3500, "complex d=2");
  fcp_run_t8iv_case<tenes::complex_tensor>('v', 4, 3600, "complex d=4");
}

// ---- T2-vi complex ---------------------------------------------------------

namespace {

/*! One block Y of contract T2-vi in a COMPLEX real-CTM environment: the
 *  measurement path on Tn1(Y), Tn2(Y) against the plain quadratic form of
 *  the phase-normalised N_plain on the masked block. The measurement path
 *  carries the window's phase, N_plain does not, so |norm| is compared with
 *  the (real positive) quadratic norm and <G> as the ratio.
 */
template <class tensor>
void fcp_run_t2vi_complex_block(const fur_window<tensor>& w,
                                const fue::full_update_environment<tensor>& N,
                                const fg_ftensor<tensor>& I,
                                const fg_ftensor<tensor>& G,
                                const fg_ftensor<tensor>& hop,
                                const fg_ftensor<tensor>& Y,
                                const std::string& what) {
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

  // Step 2: the measurement path (phased).
  const auto norm_meas = fur_closed(w, T1, T2, I);
  const auto g_meas = fur_closed(w, T1, T2, G);
  const auto hop_meas = fur_closed(w, T1, T2, hop);
  INFO(label << " [measurement path] norm=" << norm_meas << " (arg "
             << std::arg(fcp_to_complex(norm_meas)) << ")"
             << " <G>*norm=" << g_meas << " <hop>*norm=" << hop_meas);
  REQUIRE(std::isfinite(std::abs(norm_meas)));
  REQUIRE(std::abs(norm_meas) > 0.0);
  // Premise: the odd (hopping) channel carries weight on this block.
  REQUIRE(std::abs(hop_meas) > 1.0e-3 * std::abs(norm_meas));

  // Step 3: the plain quadratic form on the masked block (unphased).
  const fg_ftensor<tensor> Yt = fur_mask_ab(Y);
  const fg_ftensor<tensor> Th =
      fgf::tensordot(Y, G, mptensor::Axes(2, 3), mptensor::Axes(0, 1));
  REQUIRE(Th.parity == Y.parity);
  const fg_ftensor<tensor> Tht = fur_mask_ab(Th);
  const std::complex<double> norm_quad =
      fcp_to_complex(fur_quad(N.N_plain, Yt.t, Yt.t));
  const std::complex<double> g_quad =
      fcp_to_complex(fur_quad(N.N_plain, Yt.t, Tht.t));
  INFO(label << " [quadratic form] norm=" << norm_quad
             << " <G>*norm=" << g_quad);
  REQUIRE(std::abs(norm_quad) > 0.0);
  // N_plain is phase-normalised: its quadratic norm is real positive.
  CHECK(norm_quad.real() > 0.0);
  CHECK(std::abs(norm_quad.imag()) <= 1.0e-10 * std::abs(norm_quad));

  // Step 4: |norm| against |norm|, and the ratios.
  fue_check_rel(label + " [T2-vi complex |norm|]", std::abs(norm_quad),
                std::abs(norm_meas), 1.0e-10);
  const std::complex<double> r_meas =
      fcp_to_complex(g_meas) / fcp_to_complex(norm_meas);
  const std::complex<double> r_quad = g_quad / norm_quad;
  fue_check_rel(label + " [T2-vi complex <G>]", r_quad, r_meas, 1.0e-10);
}

template <class tensor>
void fcp_run_t2vi_complex_window(const fur_env<tensor>& e, int s1, char dir,
                                 int seed) {
  const fur_window<tensor> w = fcp_make_window(e, s1, dir, "T2-vi complex");
  INFO(w.label);
  fur_require_assembly(w);

  const auto N = fur_build_N(w);
  INFO(w.label << " [premise: forbidden_ratio = " << N.forbidden_ratio
               << ", phase = " << N.phase << "]");
  REQUIRE(N.forbidden_ratio <= 1.0e-10);
  REQUIRE(N.N_plain.shape() ==
          mptensor::Shape(w.pa.size(), w.pb.size(), w.pa.size(), w.pb.size()));
  REQUIRE(std::abs(std::abs(N.phase) - 1.0) <= 1.0e-12);

  const fg_ftensor<tensor> I = fgf::wrap_twosite_gate(
      fue_internal_identity<tensor>(w.ps1, w.ps2), w.ps1, w.ps2);
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
    REQUIRE(fur_rel_diff(Xp.t, X.t) > 1.0e-3);
  }

  fcp_run_t2vi_complex_block(w, N, I, G, hop, X, "X");
  fcp_run_t2vi_complex_block(w, N, I, G, hop, Theta, "Theta=X.G");
  fcp_run_t2vi_complex_block(w, N, I, G, hop, R1, "random-1");
  fcp_run_t2vi_complex_block(w, N, I, G, hop, R2, "random-2");
  fcp_run_t2vi_complex_block(w, N, I, G, hop, Xp, "X+0.05*random");
}

//! Contract T2-vi complex: the environment is the new CTM path only
//! (Calc_CTM_Environment_density with phase_invariant = true); no phase is
//! fixed anywhere.
template <class tensor>
void fcp_run_t2vi_complex(int seed, int chi) {
  const fur_env<tensor> e =
      fcp_make_fold_env<tensor>(seed, chi, 2, 2, 2, "complex");
  INFO("T2-vi complex " << e.label);
  fcp_run_t2vi_complex_window(e, 0, 'h', seed + 10);
  fcp_run_t2vi_complex_window(e, 3, 'v', seed + 20);
}

}  // namespace

TEST_CASE(
    "ctm phase T2-vi complex: in the complex CTM environment the plain "
    "quadratic form of the phase-normalised N_plain reproduces |norm| and <G> "
    "of the measurement path on general blocks") {
  fcp_run_t2vi_complex<tenes::complex_tensor>(7500, 8);
}

// ---- T8-v: phase normalisation of a gathered one-site RDM -----------------

namespace {

template <class V>
V fcp_rdm_value(double re, double im);
template <>
inline double fcp_rdm_value<double>(double re, double im) {
  REQUIRE(im == 0.0);
  return re;
}
template <>
inline std::complex<double> fcp_rdm_value<std::complex<double>>(double re,
                                                                double im) {
  return {re, im};
}

//! A random Hermitian positive definite d x d matrix, B B^dag + I, as a
//! non-distributed small_tensor. Its trace is real positive by
//! construction (the diagonal of B B^dag is a sum of |b|^2).
template <class V>
tenes::small_tensor<V> fcp_random_hpd(std::size_t d, int seed) {
  tenes::small_tensor<V> B(mptensor::Shape(d, d));
  for (std::size_t n = 0; n < B.local_size(); ++n) {
    const mptensor::Index idx = B.global_index(n);
    const double re = fcp_draw(seed, idx, 2);
    const double im =
        std::is_same<V, double>::value ? 0.0 : fcp_draw(seed + 300, idx, 2);
    B.set_value(idx, fcp_rdm_value<V>(re, im));
  }
  tenes::small_tensor<V> A = mptensor::tensordot(
      B, mptensor::conj(B), mptensor::Axes(1), mptensor::Axes(1));
  // The diagonal of B B^dag is real up to the roundoff of the BLAS kernel;
  // set it exactly real so that the trace has NO imaginary part and the
  // "nothing to do" branch of the contract (|phase - 1| <= 1e-14) is the one
  // exercised in (a).
  for (std::size_t i = 0; i < d; ++i) {
    V v;
    A.get_value(mptensor::Index(i, i), v);
    A.set_value(mptensor::Index(i, i),
                fcp_rdm_value<V>(fcp_to_complex(v).real() + 1.0, 0.0));
  }
  return A;
}

template <class V>
std::complex<double> fcp_small_trace(const tenes::small_tensor<V>& a) {
  std::complex<double> t = 0.0;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    const mptensor::Index idx = a.global_index(n);
    if (idx[0] == idx[1]) {
      V v;
      a.get_value(idx, v);
      t += fcp_to_complex(v);
    }
  }
  return t;
}

template <class V>
double fcp_small_frobenius(const tenes::small_tensor<V>& a) {
  double s = 0.0;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    V v;
    a.get_value(a.global_index(n), v);
    s += std::abs(v) * std::abs(v);
  }
  return std::sqrt(s);
}

template <class V>
bool fcp_small_identical(const tenes::small_tensor<V>& a,
                         const tenes::small_tensor<V>& b) {
  if (!(a.shape() == b.shape())) {
    return false;
  }
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    const mptensor::Index idx = a.global_index(n);
    V va, vb;
    a.get_value(idx, va);
    b.get_value(idx, vb);
    if (!(va == vb)) {
      return false;
    }
  }
  return true;
}

template <class V>
double fcp_small_rel_diff(const tenes::small_tensor<V>& got,
                          const tenes::small_tensor<V>& ref) {
  REQUIRE(got.shape() == ref.shape());
  double scale = 0.0;
  double dev = 0.0;
  for (std::size_t n = 0; n < ref.local_size(); ++n) {
    const mptensor::Index idx = ref.global_index(n);
    V g, r;
    got.get_value(idx, g);
    ref.get_value(idx, r);
    scale = std::max(scale, std::abs(r));
    dev = std::max(dev, std::abs(g - r));
  }
  REQUIRE(scale > 0.0);
  return dev / scale;
}

template <class V>
void fcp_run_t8v(const char* type_name, int seed) {
  const std::string label = std::string("T8-v ") + type_name;
  INFO(label);
  const std::size_t d = 4;
  const tenes::small_tensor<V> A = fcp_random_hpd<V>(d, seed);
  const std::complex<double> trA = fcp_small_trace(A);
  INFO(label << " [premise] trace(A) = " << trA);
  REQUIRE(trA.real() > 0.0);
  REQUIRE(trA.imag() == 0.0);
  const double frobA = fcp_small_frobenius(A);
  REQUIRE(frobA > 0.0);

  // (a) a real positive trace: no-op, exactly 1.
  {
    tenes::small_tensor<V> a = A;
    const std::complex<double> ph = tenes::itps::core::normalize_rdm_phase(a);
    INFO(label << " [T8-v a] returned " << ph);
    CHECK(ph == std::complex<double>(1.0, 0.0));
    CHECK(fcp_small_identical(a, A));
  }

  // (b) a known phase on the same matrix: back to A, the phase reported,
  // |trace| and the Frobenius norm kept.
  {
    const std::complex<double> factor = std::is_same<V, double>::value
                                            ? std::complex<double>(-1.0, 0.0)
                                            : std::polar(1.0, 2.3);
    tenes::small_tensor<V> b = A;
    for (std::size_t n = 0; n < b.local_size(); ++n) {
      const mptensor::Index idx = b.global_index(n);
      V v;
      b.get_value(idx, v);
      const std::complex<double> scaled = fcp_to_complex(v) * factor;
      b.set_value(idx, fcp_rdm_value<V>(scaled.real(), scaled.imag()));
    }
    const std::complex<double> tr_before = fcp_small_trace(b);
    const double frob_before = fcp_small_frobenius(b);
    // Premise: the injected phase really made the trace non-positive.
    REQUIRE(!(tr_before.real() > 0.0 &&
              std::abs(tr_before.imag()) <= 1.0e-14 * std::abs(tr_before)));
    const std::complex<double> ph = tenes::itps::core::normalize_rdm_phase(b);
    const std::complex<double> tr_after = fcp_small_trace(b);
    INFO(label << " [T8-v b] factor=" << factor << " returned " << ph
               << " trace before=" << tr_before << " after=" << tr_after);
    CHECK(std::abs(ph - factor) <= 1.0e-12);
    CHECK(fcp_small_rel_diff(b, A) <= 1.0e-12);
    CHECK(tr_after.real() > 0.0);
    CHECK(std::abs(tr_after.imag()) <= 1.0e-12 * std::abs(tr_after));
    CHECK(std::abs(std::abs(tr_after) - std::abs(tr_before)) <=
          1.0e-12 * std::abs(tr_before));
    CHECK(std::abs(fcp_small_frobenius(b) - frob_before) <=
          1.0e-12 * frob_before);
  }

  // (c) the matrices it must refuse.
  {
    tenes::small_tensor<V> zero(mptensor::Shape(d, d));
    for (std::size_t n = 0; n < zero.local_size(); ++n) {
      zero.set_value(zero.global_index(n), V(0.0));
    }
    INFO(label << " [T8-v c: zero matrix]");
    CHECK_THROWS_AS(tenes::itps::core::normalize_rdm_phase(zero),
                    std::runtime_error);
  }
  {
    // NaN off the diagonal only: the trace itself is finite and positive,
    // so a trace-only check would let it through.
    tenes::small_tensor<V> nan = A;
    nan.set_value(mptensor::Index(0, 1),
                  V(std::numeric_limits<double>::quiet_NaN()));
    const std::complex<double> tr = fcp_small_trace(nan);
    REQUIRE(std::isfinite(std::abs(tr)));
    REQUIRE(tr.real() > 0.0);
    INFO(label << " [T8-v c: NaN off-diagonal, trace = " << tr << "]");
    CHECK_THROWS_AS(tenes::itps::core::normalize_rdm_phase(nan),
                    std::runtime_error);
  }
  {
    tenes::small_tensor<V> inf = A;
    inf.set_value(mptensor::Index(1, 1),
                  V(std::numeric_limits<double>::infinity()));
    INFO(label << " [T8-v c: Inf on the diagonal]");
    CHECK_THROWS_AS(tenes::itps::core::normalize_rdm_phase(inf),
                    std::runtime_error);
  }
}

}  // namespace

TEST_CASE(
    "ctm phase T8-v: normalize_rdm_phase leaves a positive-trace RDM alone, "
    "removes a known phase and refuses zero or non-finite input") {
  fcp_run_t8v<double>("real", 8300);
  fcp_run_t8v<std::complex<double>>("complex", 8400);
}
