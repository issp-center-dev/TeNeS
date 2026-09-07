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

// ===== S-5: the warm start of iTPS::update_CTM() ===========================
//
// Contract:
// docs/superpowers/specs/2026-09-07-fermion-fast-full-update-contract.md
// sections 3 (R3) and 4 (S-5).
//
// R3: in a NON-fast fermionic run the CTM re-convergence after a bond update
// starts from the current environment, and falls back to the uniform-vector
// (cold) initialization when no environment has been built yet. The switch is
// not reachable from the input file, so - as S-5 says - this is where the
// cold/warm A/B has to be taken.
//
// The three cases:
//
//   1. same fixed point.  From one state and one converged environment,
//      perturbed, a cold and a warm re-convergence have to land on the same
//      environment, judged by the one-site reduced density matrices (which,
//      unlike the environment tensors themselves, do not depend on the gauge
//      the CTMRG happens to fix on its boundary indices).
//   2. really warm.  Case 1 alone would pass an implementation that took the
//      argument and ignored it. With a budget of zero CTM sweeps a warm start
//      keeps the converged environment while a cold start replaces it with
//      the uniform-vector one, and the two RDM sets are then far apart. This
//      is what shows the argument is honoured at all.
//   3. fallback.  A warm start on a state whose environment has never been
//      converged - straight after construction, and after a second
//      initialize_tensors(), which R3 names as an operation that invalidates
//      it - has to behave exactly like a cold start, not merely "not crash".
//
// Truth sources: no reference number is taken from update_CTM(). Case 1
// compares the two paths against each other; case 3 compares the warm path
// against the cold one; case 2 compares against the environment recorded
// before the call. The tolerances come from the measurements quoted at each
// case, taken with core::Calc_CTM_Environment_density(..., initialize=false)
// - the routine iTPS::update_CTM() already calls, with the flag R3 is about -
// on this very fixture (work/fermion/fastfu/testauthor/probe).
//
// ---------------------------------------------------------------------------
// ASSUMED SIGNATURE. The contract says only that update_CTM() "takes an
// argument that asks for a warm start"; it does not spell the argument. This
// file assumes
//
//     void iTPS<tensor>::update_CTM(bool warm_start);
//
// with false = the current cold behaviour. A defaulted parameter
// (bool warm_start = false) also compiles against this file; an enum or a
// differently named method does not. See report-testauthor.md.
// ---------------------------------------------------------------------------
//
// This is a separate executable rather than another file included into
// test_fermion_layer.cpp on purpose: until R3 lands the call below does not
// compile, and inside that translation unit a compile error would take the
// S-3 tests of fermion/fast_full_update.cpp - which do compile, and fail for
// their own reason - down with it.
//
// The bosonic half of S-5 ("a boson run is not changed by the warm start")
// is not tested here. The contract says the existing bosonic regression
// tests (AntiferroHeisenberg_real and the rest) carry it, and a unit-level
// version would have no power: a boson state that wrongly warm-started would
// still converge to the same fixed point, so the comparison would pass
// either way.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../doctest.h"
#include "../test_workdir.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include "../../src/SquareLattice.hpp"
#include "../../src/tensor.hpp"
#include "../../src/fermion/reduced_measure.hpp"
#include "../../src/iTPS/PEPS_Parameters.hpp"
#include "../../src/iTPS/iTPS.hpp"

namespace tenes::itps {
//! Local definition of the friend struct declared by iTPS (iTPS.hpp line 99).
//! test/test_fermion_layer.cpp and test/input.cpp have their own; this
//! translation unit is a separate executable, so the three do not collide.
struct iTPSTestAccessor {
  template <class tensor>
  static std::vector<tensor>& Tn(iTPS<tensor>& state) {
    return state.Tn;
  }
  template <class tensor>
  static tenes::fermion::FermionInfo& finfo(iTPS<tensor>& state) {
    return state.finfo;
  }
  //! Mutable, because case 2 has to shrink Max_CTM_Iteration after the
  //! environment has been converged; the parameter set is a copy the state
  //! owns, so there is no other way in.
  template <class tensor>
  static PEPS_Parameters& params(iTPS<tensor>& state) {
    return state.peps_parameters;
  }
  template <class tensor>
  static SquareLattice& lattice(iTPS<tensor>& state) {
    return state.lattice;
  }
};
}  // namespace tenes::itps

namespace {

using cws_tensor = tenes::real_tensor;
using cws_state = tenes::itps::iTPS<cws_tensor>;
using A = tenes::itps::iTPSTestAccessor;

// The fixture of fermion/fermion_guards.cpp: the schedule on which the
// fermionic fold CTM converges. eps is tightened from 1e-8 to 1e-10 so that
// "the two agree to within the CTM tolerance" is a statement about something
// much smaller than any real disagreement.
constexpr int cws_LX = 2;
constexpr int cws_LY = 2;
constexpr int cws_D = 2;
constexpr int cws_CHI = 8;
constexpr int cws_simple_steps = 200;
constexpr int cws_ctm_iteration_max = 400;
constexpr double cws_ctm_epsilon = 1.0e-10;
constexpr double cws_tau = 0.01;
constexpr const char* cws_ctm_warning = "CTM did not converge";

//! Redirect one stream into a string for the lifetime of the object (the
//! same device fermion/fermion_guards.cpp uses; a doctest cannot read the
//! process' own stdout).
class cws_capture {
 public:
  explicit cws_capture(std::ostream& stream) : stream_(stream) {
    saved_ = stream_.rdbuf(buffer_.rdbuf());
  }
  ~cws_capture() { stream_.rdbuf(saved_); }
  cws_capture(cws_capture const&) = delete;
  cws_capture& operator=(cws_capture const&) = delete;
  std::string str() const { return buffer_.str(); }

 private:
  std::ostream& stream_;
  std::ostringstream buffer_;
  std::streambuf* saved_ = nullptr;
};

cws_tensor cws_gate(double tau) {
  cws_tensor g(mptensor::Shape(2, 2, 2, 2));
  g.set_value(mptensor::Index(0, 0, 0, 0), 1.0);
  g.set_value(mptensor::Index(1, 1, 1, 1), 1.0);
  g.set_value(mptensor::Index(0, 1, 0, 1), std::cosh(tau));
  g.set_value(mptensor::Index(1, 0, 1, 0), std::cosh(tau));
  g.set_value(mptensor::Index(0, 1, 1, 0), std::sinh(tau));
  g.set_value(mptensor::Index(1, 0, 0, 1), std::sinh(tau));
  return g;
}

tenes::SquareLattice cws_make_lattice() {
  tenes::SquareLattice lattice(cws_LX, cws_LY);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    lattice.physical_dims[site] = 2;
    lattice.virtual_dims[site] = {cws_D, cws_D, cws_D, cws_D};
    lattice.initial_dirs[site] = {0.0};
    lattice.noises[site] = 1.0;
  }
  return lattice;
}

tenes::itps::PEPS_Parameters cws_make_params(
    const tenes::SquareLattice& lattice, const std::string& outdir) {
  tenes::itps::PEPS_Parameters params;
  params.fermion = true;
  params.phys_parity.assign(lattice.N_UNIT, std::vector<bool>{false, true});
  params.print_level = tenes::PrintLevel::warn;
  params.outdir = outdir;
  params.CHI = cws_CHI;
  params.Max_CTM_Iteration = cws_ctm_iteration_max;
  params.CTM_Convergence_Epsilon = cws_ctm_epsilon;
  params.Use_RSVD = false;
  params.seed = 11;
  return params;
}

std::vector<tenes::EvolutionOperator<cws_tensor>> cws_updates(
    const tenes::SquareLattice& lattice) {
  std::vector<tenes::EvolutionOperator<cws_tensor>> updates;
  const cws_tensor gate = cws_gate(cws_tau);
  for (int leg : {2, 1}) {
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      updates.push_back(tenes::make_twosite_EvolutionOperator<cws_tensor>(
          site, leg, 0, gate));
    }
  }
  return updates;
}

//! A fermionic state after `steps` simple-update sweeps. Deterministic: two
//! calls with the same arguments give bit-identical states, which is what
//! lets the cold and the warm run below be two separate objects with the
//! same history.
std::unique_ptr<cws_state> cws_build(const std::string& outdir, int steps) {
  const tenes::SquareLattice lattice = cws_make_lattice();
  const auto updates = cws_updates(lattice);
  auto state = std::make_unique<cws_state>(
      MPI_COMM_WORLD, cws_make_params(lattice, outdir), lattice, updates,
      tenes::EvolutionOperators<cws_tensor>{}, tenes::Operators<cws_tensor>{},
      tenes::Operators<cws_tensor>{}, tenes::Operators<cws_tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  for (int step = 0; step < steps; ++step) {
    for (const auto& up : updates) {
      state->simple_update(up);
    }
  }
  return state;
}

//! Multiply every element of every Tn by 1 + amp * f(index). Elementwise
//! scaling can only shrink the support, so the parity mask of the fermionic
//! state (every entry of odd total parity is zero) survives it untouched;
//! adding a random tensor would not.
void cws_perturb(cws_state& state, double amp) {
  auto& Tn = A::Tn(state);
  for (std::size_t site = 0; site < Tn.size(); ++site) {
    for (std::size_t n = 0; n < Tn[site].local_size(); ++n) {
      const mptensor::Index idx = Tn[site].global_index(n);
      double v;
      Tn[site].get_value(idx, v);
      double x = 0.0;
      for (std::size_t ax = 0; ax < idx.size(); ++ax) {
        x += (ax + 1.7) * static_cast<double>(idx[ax]);
      }
      Tn[site].set_value(idx, v * (1.0 + amp * std::sin(3.1 * x + 0.7 * site)));
    }
  }
}

//! The one-site reduced density matrices of the current environment,
//! normalized to unit trace. Gauge invariant, unlike the environment tensors
//! themselves: a CTMRG started from a different point converges to the same
//! fixed point only up to a gauge on the chi legs, so comparing C1 and eTl
//! element by element would fail on a correct implementation.
std::vector<std::vector<double>> cws_rdms(cws_state& state) {
  auto raw = state.measure_onesite_rdm();
  std::vector<std::vector<double>> out;
  for (auto& rdm : raw) {
    const auto shape = rdm.shape();
    double trace = 0.0;
    for (std::size_t i = 0; i < shape[0]; ++i) {
      double v;
      rdm.get_value(mptensor::Index(i, i), v);
      trace += v;
    }
    REQUIRE(std::abs(trace) > 0.0);
    std::vector<double> flat;
    for (std::size_t n = 0; n < rdm.local_size(); ++n) {
      double v;
      rdm.get_value(rdm.global_index(n), v);
      flat.push_back(v / trace);
    }
    out.push_back(flat);
  }
  return out;
}

double cws_rdm_distance(const std::vector<std::vector<double>>& a,
                        const std::vector<std::vector<double>>& b) {
  REQUIRE(a.size() == b.size());
  double worst = 0.0;
  for (std::size_t site = 0; site < a.size(); ++site) {
    REQUIRE(a[site].size() == b[site].size());
    for (std::size_t i = 0; i < a[site].size(); ++i) {
      worst = std::max(worst, std::abs(a[site][i] - b[site][i]));
    }
  }
  return worst;
}

void cws_cleanup(const std::string& outdir) {
  std::error_code ec;
  std::filesystem::remove_all(outdir, ec);
}

}  // namespace

TEST_CASE(
    "fermion CTM warm start S-5: a warm and a cold re-convergence reach the "
    "same fixed point") {
  // Measured on this fixture with core::Calc_CTM_Environment_density's
  // initialize flag (work/fermion/fastfu/testauthor/probe, 2026-09-07):
  // the worst difference between the two one-site RDM sets is 1.7e-11 for a
  // 1e-3 perturbation, 7.8e-11 for 1e-2 and 4.3e-11 for 1e-1, i.e. a few
  // hundred times CTM_Convergence_Epsilon, and it tracks that setting
  // linearly (at 1e-8 it was 1.7e-9 ... 5.0e-9). The bound below is 1000
  // times the epsilon: four decades above what was measured, and four
  // decades below the 1e-3 .. 1e-1 scale on which two CTMRGs that had gone
  // to DIFFERENT fixed points would differ.
  const double tol = 1.0e3 * cws_ctm_epsilon;
  const double amp = 1.0e-3;

  auto cold =
      cws_build("output_test_fermion_warm_start_same_cold", cws_simple_steps);
  auto warm =
      cws_build("output_test_fermion_warm_start_same_warm", cws_simple_steps);
  REQUIRE(A::finfo(*cold).enabled);

  std::string warnings;
  std::vector<std::vector<double>> before;
  {
    cws_capture out(std::cout);
    cws_capture err(std::cerr);
    cold->update_CTM(false);
    warm->update_CTM(false);
    before = cws_rdms(*warm);
    warnings = out.str() + err.str();
  }
  // Premise: the environment both runs start from is converged.
  INFO("while building the common environment: " << warnings);
  REQUIRE(warnings.find(cws_ctm_warning) == std::string::npos);

  cws_perturb(*cold, amp);
  cws_perturb(*warm, amp);

  std::string cold_warnings;
  std::string warm_warnings;
  {
    cws_capture out(std::cout);
    cws_capture err(std::cerr);
    cold->update_CTM(false);
    cold_warnings = out.str() + err.str();
  }
  {
    cws_capture out(std::cout);
    cws_capture err(std::cerr);
    warm->update_CTM(true);
    warm_warnings = out.str() + err.str();
  }
  // Premise: both re-convergences converged. Comparing two environments that
  // were stopped early would say nothing about the fixed point.
  INFO("cold re-convergence: " << cold_warnings);
  INFO("warm re-convergence: " << warm_warnings);
  REQUIRE(cold_warnings.find(cws_ctm_warning) == std::string::npos);
  REQUIRE(warm_warnings.find(cws_ctm_warning) == std::string::npos);

  const auto cold_rdm = cws_rdms(*cold);
  const auto warm_rdm = cws_rdms(*warm);

  // Premise: the perturbation actually moved the environment, otherwise both
  // runs would be reproducing an unchanged answer and the agreement below
  // would be empty. Measured: 6.6e-4 for amp = 1e-3.
  const double moved = cws_rdm_distance(before, cold_rdm);
  INFO("the perturbation moved the one-site RDMs by " << moved);
  REQUIRE(moved > 1.0e-5);

  const double distance = cws_rdm_distance(cold_rdm, warm_rdm);
  INFO("worst |rho_cold - rho_warm| = " << distance << " (tol " << tol
                                        << ", CTM_Convergence_Epsilon "
                                        << cws_ctm_epsilon << ")");
  CHECK(distance <= tol);

  cws_cleanup("output_test_fermion_warm_start_same_cold");
  cws_cleanup("output_test_fermion_warm_start_same_warm");
}

TEST_CASE(
    "fermion CTM warm start S-5: the warm start really reuses the "
    "environment") {
  // Without this case, an update_CTM(bool) that accepted the argument and
  // threw it away would pass every other check in this file: a cold restart
  // converges to the same fixed point, only more slowly.
  //
  // The discriminator is a budget of zero CTM sweeps. A warm start then has
  // nothing to do and keeps the converged environment; a cold start replaces
  // it with the uniform-vector initialization and does not iterate at all.
  // Measured on this fixture (probe, 2026-09-07): the two RDM sets are then
  // 6.20e-1 apart for a 1e-3 perturbation and 6.35e-1 for 1e-2, while the
  // warm result stays within 6.6e-4 (resp. 6.6e-3) of the environment it
  // started from. Both runs print the CTM warning (count = 0), which is why
  // this case does not ask for silence.
  const double amp = 1.0e-3;
  const double cold_warm_min = 1.0e-1;   // measured 6.2e-1
  const double warm_drift_max = 1.0e-2;  // measured 6.6e-4

  auto cold =
      cws_build("output_test_fermion_warm_start_reuse_cold", cws_simple_steps);
  auto warm =
      cws_build("output_test_fermion_warm_start_reuse_warm", cws_simple_steps);

  std::string warnings;
  std::vector<std::vector<double>> converged;
  {
    cws_capture out(std::cout);
    cws_capture err(std::cerr);
    cold->update_CTM(false);
    warm->update_CTM(false);
    converged = cws_rdms(*warm);
    warnings = out.str() + err.str();
  }
  INFO("while building the common environment: " << warnings);
  REQUIRE(warnings.find(cws_ctm_warning) == std::string::npos);

  cws_perturb(*cold, amp);
  cws_perturb(*warm, amp);
  A::params(*cold).Max_CTM_Iteration = 0;
  A::params(*warm).Max_CTM_Iteration = 0;

  {
    cws_capture out(std::cout);
    cws_capture err(std::cerr);
    cold->update_CTM(false);
    warm->update_CTM(true);
  }

  const auto cold_rdm = cws_rdms(*cold);
  const auto warm_rdm = cws_rdms(*warm);
  const double separation = cws_rdm_distance(cold_rdm, warm_rdm);
  const double drift = cws_rdm_distance(converged, warm_rdm);
  INFO("with a budget of zero sweeps: |rho_cold - rho_warm| = "
       << separation << " (must exceed " << cold_warm_min
       << "; equal values mean update_CTM ignored its argument), "
       << "|rho_warm - rho_converged| = " << drift << " (must stay below "
       << warm_drift_max << "; a larger drift means the warm start did not "
       << "keep the environment it was given)");
  CHECK(separation > cold_warm_min);
  CHECK(drift < warm_drift_max);

  cws_cleanup("output_test_fermion_warm_start_reuse_cold");
  cws_cleanup("output_test_fermion_warm_start_reuse_warm");
}

TEST_CASE(
    "fermion CTM warm start S-5: a warm start with no environment falls back "
    "to a cold one") {
  // R3's fallback. Both halves compare the warm result against the cold one
  // rather than merely asserting that nothing threw: "it happened to run" is
  // exactly the outcome a broken fallback would also produce.
  //
  // The two states of each half have identical histories, so a warm start
  // that correctly falls back reproduces the cold run bit for bit; the bound
  // is a machine-precision one, not a physical tolerance.
  const double tol = 1.0e-10;
  // Twenty simple-update steps are enough to leave the initial product-ish
  // state; this case is about the fallback, not about convergence, so it
  // does not need the long schedule.
  const int steps = 20;

  SUBCASE("straight after construction, before any environment was built") {
    auto cold = cws_build("output_test_fermion_warm_start_fresh_cold", steps);
    auto warm = cws_build("output_test_fermion_warm_start_fresh_warm", steps);
    {
      cws_capture out(std::cout);
      cws_capture err(std::cerr);
      REQUIRE_NOTHROW(cold->update_CTM(false));
      REQUIRE_NOTHROW(warm->update_CTM(true));
    }
    const double distance = cws_rdm_distance(cws_rdms(*cold), cws_rdms(*warm));
    INFO("worst |rho_cold - rho_warm| on a state with no environment = "
         << distance);
    CHECK(distance <= tol);
    cws_cleanup("output_test_fermion_warm_start_fresh_cold");
    cws_cleanup("output_test_fermion_warm_start_fresh_warm");
  }

  SUBCASE("after initialize_tensors() threw the environment away") {
    // R3: "the state of whether the environment is valid has to go back to
    // invalid on an operation that can change the shape of the environment
    // tensors (initializing the tensors, loading a checkpoint)".
    // initialize_tensors() clears C1..eTl and rebuilds them empty
    // (src/iTPS/tensors.cpp), so a warm start that still believed in the
    // environment would iterate from those.
    auto cold = cws_build("output_test_fermion_warm_start_reinit_cold", steps);
    auto warm = cws_build("output_test_fermion_warm_start_reinit_warm", steps);
    {
      cws_capture out(std::cout);
      cws_capture err(std::cerr);
      cold->update_CTM(false);
      warm->update_CTM(false);
      cold->initialize_tensors();
      warm->initialize_tensors();
      REQUIRE_NOTHROW(cold->update_CTM(false));
      REQUIRE_NOTHROW(warm->update_CTM(true));
    }
    // Premise: the two states really are the same state, so that any
    // difference below can only come from the warm-start argument.
    {
      const auto& a = A::Tn(*cold);
      const auto& b = A::Tn(*warm);
      REQUIRE(a.size() == b.size());
      double worst = 0.0;
      for (std::size_t site = 0; site < a.size(); ++site) {
        worst = std::max(worst, mptensor::max_abs(a[site] - b[site]));
      }
      INFO("the two re-initialized states differ by " << worst);
      REQUIRE(worst == 0.0);
    }
    const double distance = cws_rdm_distance(cws_rdms(*cold), cws_rdms(*warm));
    INFO("worst |rho_cold - rho_warm| after initialize_tensors() = "
         << distance);
    CHECK(distance <= tol);
    cws_cleanup("output_test_fermion_warm_start_reinit_cold");
    cws_cleanup("output_test_fermion_warm_start_reinit_warm");
  }
}
