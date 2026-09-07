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

// ===== S-3: which CTM moves the fermionic fast full update performs ========
//
// Contract:
// docs/superpowers/specs/2026-09-07-fermion-fast-full-update-contract.md
// sections 2, 3 (R2) and 4 (S-3).
//
// R2 fixes, for each source_leg, the two partial CTM moves that replace the
// full re-convergence after one full-update bond:
//
//   source_leg 0 (left)   Right_move_single (source column) then
//                         Left_move_single  (target column)
//   source_leg 1 (up)     Bottom_move_single(source row)    then
//                         Top_move_single   (target row)
//   source_leg 2 (right)  Left_move_single  (source column) then
//                         Right_move_single (target column)
//   source_leg 3 (down)   Top_move_single   (source row)    then
//                         Bottom_move_single(target row)
//
// and the rows/columns are those of the ORIGINAL source/target, not of the
// two sites after the driver has swapped them (it does, for source_leg 0 and
// 1; src/iTPS/full_update.cpp, the `if (source_leg == 0 || source_leg == 1)`
// block).
//
// S-3 requires this to be checked from the OUTSIDE, i.e. from the environment
// tensors, not by asking the implementation what it called. Two independent
// checks do that here.
//
//   (1) ffu_check_structure - which of the 8 x N_UNIT environment slots
//       changed at all. Contract section 2: each move rewrites a different
//       set of tensor KINDS ({C1,C4,eTl} / {C2,C3,eTr} / {C1,C2,eTt} /
//       {C3,C4,eTb}), and it rewrites the row/column NEXT TO the one it
//       absorbed (right(ix) / left(ix) / bottom(iy) / top(iy)). The expected
//       sets below are built from those two facts alone - nothing about how
//       the moves are implemented, and in particular nothing about which
//       tensor is fed to them. On the 3x3 unit cell used here left(ix),
//       ix and right(ix) are three different columns, so the pair
//       (kind set, column) identifies both the move and its argument.
//
//   (2) ffu_check_oracle - the numbers. The two moves of the R2 table are
//       replayed here, from the environment as it was BEFORE the bond, on
//       the state as it is AFTER it, and the result has to be the
//       environment the solver produced. This one also sees the ORDER of the
//       two moves (the second move reads what the first one wrote), which
//       (1) cannot see, and it sees a move that was called with the right
//       argument through a wrong route.
//
// Truth source: neither check takes a reference number from the routine
// under test. (1) is a set built from the contract's table, verified against
// src/iTPS/core/ctm_single.cpp; (2) calls core::*_move_single, which is the
// bosonic finite-temperature CTM move that this change does not touch, and
// which the contract (section 2) states the fermionic path already goes
// through - fermion-ness lives in the reduced tensor, not in the move.
//
// Included into the test_fermion_layer TU AFTER fermion/fermion_guards.cpp:
// it reuses that file's stream capture (fgd_stream_capture) and the TU's
// make_free_fermion_gate(). It uses the environment accessors added to
// iTPSTestAccessor at the top of test_fermion_layer.cpp.

#include <array>
#include <filesystem>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "../../src/fermion/reduced_measure.hpp"
#include "../../src/iTPS/core/ctm.hpp"

namespace {

using ffu_tensor = tenes::real_tensor;
using ffu_state = tenes::itps::iTPS<ffu_tensor>;

// ---- the fixture ----------------------------------------------------------
//
// A 3x3 unit cell: the smallest one in which left(ix), ix and right(ix) are
// three different columns, which is what lets check (1) tell
// "Left_move_single(source column)" from "Right_move_single(target column)".
// On the 2x2 cell used elsewhere in this suite the two would rewrite the very
// same column and mutation (b) of the contract would be invisible to (1).
//
// The schedule is the one fermion/fermion_guards.cpp established for a
// full-update run whose fold CTM is converged (200 simple-update steps,
// D = 2, chi = 8, iteration_max = 100, convergence_epsilon = 1e-8). An
// unconverged CTM leaks parity and makes build_full_update_environment throw
// on its forbidden-block check, which would fail these tests for a reason
// that has nothing to do with the move table; the premise is asserted below.
constexpr int ffu_LX = 3;
constexpr int ffu_LY = 3;
constexpr int ffu_D = 2;
constexpr int ffu_CHI = 8;
constexpr int ffu_simple_steps = 60;
constexpr int ffu_ctm_iteration_max = 100;
constexpr double ffu_ctm_epsilon = 1.0e-8;
constexpr double ffu_tau = 0.01;

//! The eight environment slots, in a fixed order, so that a changed slot can
//! be named as (kind, site).
enum ffu_kind {
  FFU_C1 = 0,
  FFU_C2,
  FFU_C3,
  FFU_C4,
  FFU_ETT,
  FFU_ETR,
  FFU_ETB,
  FFU_ETL,
  FFU_NKIND
};

//! Returned by value (not as a const char*) because doctest's INFO()
//! stringifies a raw pointer as an address, not as the text it points at.
std::string ffu_kind_name(int kind) {
  static const char* names[FFU_NKIND] = {"C1",  "C2",  "C3",  "C4",
                                         "eTt", "eTr", "eTb", "eTl"};
  return names[kind];
}

//! The four partial moves of contract section 2.
enum class ffu_move { left, right, top, bottom };

std::string ffu_move_name(ffu_move m) {
  switch (m) {
    case ffu_move::left:
      return "Left_move_single";
    case ffu_move::right:
      return "Right_move_single";
    case ffu_move::top:
      return "Top_move_single";
    case ffu_move::bottom:
      return "Bottom_move_single";
  }
  return "?";
}

//! A copy of the ten environment vectors.
struct ffu_env {
  std::array<std::vector<ffu_tensor>, FFU_NKIND> slot;
};

ffu_env ffu_snapshot(ffu_state& state) {
  using A = tenes::itps::iTPSTestAccessor;
  ffu_env env;
  env.slot[FFU_C1] = A::C1(state);
  env.slot[FFU_C2] = A::C2(state);
  env.slot[FFU_C3] = A::C3(state);
  env.slot[FFU_C4] = A::C4(state);
  env.slot[FFU_ETT] = A::eTt(state);
  env.slot[FFU_ETR] = A::eTr(state);
  env.slot[FFU_ETB] = A::eTb(state);
  env.slot[FFU_ETL] = A::eTl(state);
  return env;
}

//! max |a - b|, or -1 when the two do not even have the same shape (which
//! counts as "changed" and is reported as such).
double ffu_diff(const ffu_tensor& a, const ffu_tensor& b) {
  if (!(a.shape() == b.shape())) {
    return -1.0;
  }
  return mptensor::max_abs(a - b);
}

//! Apply one move of the contract's table to a copy of the environment.
void ffu_apply_move(ffu_move m, int coord, ffu_env& env,
                    const std::vector<ffu_tensor>& Tn_single,
                    const tenes::itps::PEPS_Parameters& params,
                    const tenes::SquareLattice& lattice) {
  namespace core = tenes::itps::core;
  auto& C1 = env.slot[FFU_C1];
  auto& C2 = env.slot[FFU_C2];
  auto& C3 = env.slot[FFU_C3];
  auto& C4 = env.slot[FFU_C4];
  auto& eTt = env.slot[FFU_ETT];
  auto& eTr = env.slot[FFU_ETR];
  auto& eTb = env.slot[FFU_ETB];
  auto& eTl = env.slot[FFU_ETL];
  switch (m) {
    case ffu_move::left:
      core::Left_move_single(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn_single,
                             coord, params, lattice);
      break;
    case ffu_move::right:
      core::Right_move_single(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn_single,
                              coord, params, lattice);
      break;
    case ffu_move::top:
      core::Top_move_single(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn_single,
                            coord, params, lattice);
      break;
    case ffu_move::bottom:
      core::Bottom_move_single(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn_single,
                               coord, params, lattice);
      break;
  }
}

//! The (kind, site) slots one move rewrites, straight from the right two
//! columns of the table in contract section 2:
//!
//!   Left_move_single(ix)   C1 C4 eTl   on column right(ix)
//!   Right_move_single(ix)  C2 C3 eTr   on column left(ix)
//!   Top_move_single(iy)    C1 C2 eTt   on row    bottom(iy)
//!   Bottom_move_single(iy) C3 C4 eTb   on row    top(iy)
//!
//! The row/column is expressed through lattice.right() / left() / bottom() /
//! top() rather than through arithmetic on the coordinate, so that the
//! expectation follows the geometry the solver itself uses.
std::set<int> ffu_move_kinds(ffu_move m) {
  switch (m) {
    case ffu_move::left:
      return {FFU_C1, FFU_C4, FFU_ETL};
    case ffu_move::right:
      return {FFU_C2, FFU_C3, FFU_ETR};
    case ffu_move::top:
      return {FFU_C1, FFU_C2, FFU_ETT};
    case ffu_move::bottom:
      return {FFU_C3, FFU_C4, FFU_ETB};
  }
  return {};
}

std::vector<int> ffu_move_sites(ffu_move m, int coord,
                                const tenes::SquareLattice& lattice) {
  std::vector<int> sites;
  switch (m) {
    case ffu_move::left:
      for (int iy = 0; iy < lattice.LY; ++iy) {
        sites.push_back(lattice.right(lattice.index(coord, iy)));
      }
      break;
    case ffu_move::right:
      for (int iy = 0; iy < lattice.LY; ++iy) {
        sites.push_back(lattice.left(lattice.index(coord, iy)));
      }
      break;
    case ffu_move::top:
      for (int ix = 0; ix < lattice.LX; ++ix) {
        sites.push_back(lattice.bottom(lattice.index(ix, coord)));
      }
      break;
    case ffu_move::bottom:
      for (int ix = 0; ix < lattice.LX; ++ix) {
        sites.push_back(lattice.top(lattice.index(ix, coord)));
      }
      break;
  }
  return sites;
}

std::set<std::pair<int, int>> ffu_expected_changed(
    ffu_move m, int coord, const tenes::SquareLattice& lattice) {
  std::set<std::pair<int, int>> out;
  for (int kind : ffu_move_kinds(m)) {
    for (int site : ffu_move_sites(m, coord, lattice)) {
      out.insert({kind, site});
    }
  }
  return out;
}

//! The two moves R2 prescribes for a bond named (source, source_leg).
std::array<std::pair<ffu_move, int>, 2> ffu_expected_moves(
    int source, int source_leg, const tenes::SquareLattice& lattice) {
  const int target = lattice.neighbor(source, source_leg);
  // Deliberately the ORIGINAL source and target: the driver swaps the two
  // sites for source_leg 0 and 1, and R2 says the moves must not follow that
  // swap.
  const int sx = lattice.x(source);
  const int sy = lattice.y(source);
  const int tx = lattice.x(target);
  const int ty = lattice.y(target);
  switch (source_leg) {
    case 0:
      return {std::make_pair(ffu_move::right, sx),
              std::make_pair(ffu_move::left, tx)};
    case 1:
      return {std::make_pair(ffu_move::bottom, sy),
              std::make_pair(ffu_move::top, ty)};
    case 2:
      return {std::make_pair(ffu_move::left, sx),
              std::make_pair(ffu_move::right, tx)};
    default:
      return {std::make_pair(ffu_move::top, sy),
              std::make_pair(ffu_move::bottom, ty)};
  }
}

std::string ffu_slot_list(const std::set<std::pair<int, int>>& slots) {
  std::ostringstream os;
  os << "{";
  bool first = true;
  for (const auto& s : slots) {
    if (!first) os << ", ";
    first = false;
    os << ffu_kind_name(s.first) << "[" << s.second << "]";
  }
  os << "}";
  return os.str();
}

// ---- building the state ---------------------------------------------------

tenes::SquareLattice ffu_make_lattice() {
  tenes::SquareLattice lattice(ffu_LX, ffu_LY);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    lattice.physical_dims[site] = 2;
    lattice.virtual_dims[site] = {ffu_D, ffu_D, ffu_D, ffu_D};
    lattice.initial_dirs[site] = {0.0};
    lattice.noises[site] = 1.0;
  }
  return lattice;
}

tenes::itps::PEPS_Parameters ffu_make_params(
    const tenes::SquareLattice& lattice, bool fast, const std::string& outdir) {
  tenes::itps::PEPS_Parameters params;
  params.fermion = true;
  params.phys_parity.assign(lattice.N_UNIT,
                            std::vector<bool>{false, true});  // spinless
  params.print_level = tenes::PrintLevel::warn;
  params.outdir = outdir;
  params.CHI = ffu_CHI;
  params.Max_CTM_Iteration = ffu_ctm_iteration_max;
  params.CTM_Convergence_Epsilon = ffu_ctm_epsilon;
  params.Use_RSVD = false;
  params.seed = 11;
  params.Full_Use_FastFullUpdate = fast;
  return params;
}

//! One gate per site per direction, always named from the raster-earlier end
//! (source_leg 2 and 1), which is what the simple update wants.
std::vector<tenes::EvolutionOperator<ffu_tensor>> ffu_make_simple_updates(
    const tenes::SquareLattice& lattice, const ffu_tensor& gate) {
  std::vector<tenes::EvolutionOperator<ffu_tensor>> updates;
  for (int leg : {2, 1}) {
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      updates.push_back(tenes::make_twosite_EvolutionOperator<ffu_tensor>(
          site, leg, 0, gate));
    }
  }
  return updates;
}

//! A fermionic state with a converged CTM, ready for one full-update bond.
//! Deterministic: the seed and the schedule are fixed, so two calls with the
//! same arguments produce bit-identical states.
std::unique_ptr<ffu_state> ffu_build_state(
    const tenes::SquareLattice& lattice,
    const tenes::itps::PEPS_Parameters& params,
    const std::vector<tenes::EvolutionOperator<ffu_tensor>>& simple_updates,
    const std::vector<tenes::EvolutionOperator<ffu_tensor>>& full_updates,
    std::string& warnings) {
  auto state = std::make_unique<ffu_state>(
      MPI_COMM_WORLD, params, lattice, simple_updates, full_updates,
      tenes::Operators<ffu_tensor>{}, tenes::Operators<ffu_tensor>{},
      tenes::Operators<ffu_tensor>{}, tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  {
    fgd_stream_capture out(std::cout);
    fgd_stream_capture err(std::cerr);
    for (int step = 0; step < ffu_simple_steps; ++step) {
      for (const auto& up : simple_updates) {
        state->simple_update(up);
      }
    }
    state->update_CTM();
    warnings = out.str() + err.str();
  }
  return state;
}

// ---- the two checks -------------------------------------------------------

//! Check (1): every environment slot that changed lies inside the two moves
//! of the contract's table, and each of the two moves left a mark on every
//! tensor kind it owns.
//!
//! Why the second half is "at least one site per kind" and not "every site of
//! the rewritten row/column": a move rewrites its whole row/column, but a
//! slot only changes OBSERVABLY when the two sites the bond touched reach it.
//! Calc_Next_CTM_single builds its two corners from disjoint inputs
//! (src/iTPS/core/ctm_single.cpp), so on a 3x3 cell with two moved sites some
//! corners of the rewritten row are recomputed from untouched tensors and come
//! out bit-identical. Measured with a conforming implementation: for
//! source_leg 1, C1[3] C2[5] C3[8] C4[6] are rewritten but unchanged, while
//! every eT of the row and the other two corners of each kind do change.
//! Demanding all of them would fail a correct implementation.
//!
//! The half that carries the detection power is the first one: a move on the
//! wrong row/column, an extra move, or a full CTM re-convergence all show up
//! as a slot outside `expected`. The per-kind half is what notices a move
//! that was not performed at all; the two moves of a pair always own disjoint
//! kinds (contract section 2), which the caller asserts.
void ffu_check_structure(const ffu_env& before, const ffu_env& after,
                         const std::set<std::pair<int, int>>& expected,
                         const std::array<std::set<int>, 2>& move_kinds,
                         const tenes::SquareLattice& lattice,
                         const std::string& label) {
  std::set<std::pair<int, int>> changed;
  for (int kind = 0; kind < FFU_NKIND; ++kind) {
    REQUIRE(before.slot[kind].size() ==
            static_cast<std::size_t>(lattice.N_UNIT));
    REQUIRE(after.slot[kind].size() == before.slot[kind].size());
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      const double d =
          ffu_diff(before.slot[kind][site], after.slot[kind][site]);
      if (d != 0.0) {
        changed.insert({kind, site});
      }
    }
  }

  std::set<std::pair<int, int>> extra;
  for (const auto& s : changed) {
    if (expected.count(s) == 0) extra.insert(s);
  }

  INFO(label << ": slots the table allows to change "
             << ffu_slot_list(expected));
  INFO(label << ": slots that actually changed      "
             << ffu_slot_list(changed));
  INFO(label << ": changed but not allowed " << ffu_slot_list(extra)
             << " (a move that the table does not ask for ran, or a move ran "
                "on the wrong row/column, or the whole CTM was re-converged)");
  CHECK(extra.empty());

  for (int which = 0; which < 2; ++which) {
    for (int kind : move_kinds[which]) {
      int hits = 0;
      for (const auto& s : changed) {
        if (s.first == kind) ++hits;
      }
      INFO(label << ": move " << (which + 1) << " of the pair owns "
                 << ffu_kind_name(kind) << " and changed " << hits << " of its "
                 << lattice.N_UNIT << " sites (0 means the move never ran)");
      CHECK(hits > 0);
    }
  }
}

//! Check (2): replay the two moves and compare the numbers.
void ffu_check_oracle(const ffu_env& before, const ffu_env& after,
                      const std::array<std::pair<ffu_move, int>, 2>& moves,
                      const std::vector<ffu_tensor>& Tn_single,
                      const tenes::itps::PEPS_Parameters& params,
                      const tenes::SquareLattice& lattice,
                      const std::string& label) {
  ffu_env replay = before;
  for (const auto& m : moves) {
    ffu_apply_move(m.first, m.second, replay, Tn_single, params, lattice);
  }

  double worst = 0.0;
  double scale = 0.0;
  std::string worst_where = "(none)";
  bool shape_mismatch = false;
  for (int kind = 0; kind < FFU_NKIND; ++kind) {
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      const ffu_tensor& want = replay.slot[kind][site];
      const ffu_tensor& got = after.slot[kind][site];
      scale = std::max(scale, mptensor::max_abs(want));
      if (!(want.shape() == got.shape())) {
        shape_mismatch = true;
        worst_where = ffu_kind_name(kind) + "[" + std::to_string(site) +
                      "] (shape mismatch)";
        continue;
      }
      const double d = mptensor::max_abs(got - want);
      if (d > worst) {
        worst = d;
        worst_where = ffu_kind_name(kind) + "[" + std::to_string(site) + "]";
      }
    }
  }

  // The replay runs the very same routine on the very same inputs, so an
  // implementation that does what R2 says reproduces it to the last bit; the
  // bound is only there so that a future reordering of a sum inside the move
  // does not turn the test red. It is nine orders below the size of the
  // tensors being compared.
  const double tol = 1.0e-11 * std::max(scale, 1.0);
  INFO(label << " [oracle]: replayed " << ffu_move_name(moves[0].first) << "("
             << moves[0].second << ") then " << ffu_move_name(moves[1].first)
             << "(" << moves[1].second << ")");
  INFO(label << " [oracle]: worst |got - replay| = " << worst << " at "
             << worst_where << " (tol " << tol << ", largest entry of the "
             << "replayed environment " << scale << ")");
  CHECK_FALSE(shape_mismatch);
  CHECK(worst <= tol);
  // Premise: the replay is not trivially equal to the input, otherwise the
  // comparison above would pass for any implementation that leaves the
  // environment alone.
  double moved = 0.0;
  for (int kind = 0; kind < FFU_NKIND; ++kind) {
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      const double d =
          ffu_diff(before.slot[kind][site], replay.slot[kind][site]);
      moved = std::max(moved, d == -1.0 ? scale : d);
    }
  }
  INFO(label << " [oracle premise]: the two moves move the environment by "
             << moved);
  REQUIRE(moved > 1.0e-12 * std::max(scale, 1.0));
}

}  // namespace

TEST_CASE(
    "fermion fast full update S-3: one bond runs exactly the two CTM moves "
    "of the contract's table") {
  using namespace tenes;
  using namespace tenes::itps;
  using A = iTPSTestAccessor;

  const SquareLattice lattice = ffu_make_lattice();
  const ffu_tensor gate = make_free_fermion_gate(ffu_tau);
  const auto simple_updates = ffu_make_simple_updates(lattice, gate);

  // The bond is always the one leaving the centre of the 3x3 cell, so that
  // its source column/row, the column/row to its left/above and the one to
  // its right/below are three different ones.
  const int source = lattice.index(1, 1);

  for (int source_leg = 0; source_leg < 4; ++source_leg) {
    const int target = lattice.neighbor(source, source_leg);
    const std::string label = "source_leg " + std::to_string(source_leg) +
                              " (site " + std::to_string(source) + " -> " +
                              std::to_string(target) + ")";
    INFO(label);

    std::vector<EvolutionOperator<ffu_tensor>> full_updates{
        make_twosite_EvolutionOperator<ffu_tensor>(source, source_leg, 0,
                                                   gate)};

    const PEPS_Parameters params =
        ffu_make_params(lattice, /*fast=*/true,
                        "output_test_fermion_fast_full_update_leg" +
                            std::to_string(source_leg));
    REQUIRE(params.Full_Use_FastFullUpdate == true);

    std::string warnings;
    auto state = ffu_build_state(lattice, params, simple_updates, full_updates,
                                 warnings);
    // Premise: the environment this test starts from is converged. An
    // unconverged fold CTM leaks parity and trips the forbidden-block guard
    // of build_full_update_environment, which would fail the checks below
    // for a reason that is not about the move table.
    INFO(label << ": warnings while preparing the state: " << warnings);
    REQUIRE(warnings.find("CTM did not converge") == std::string::npos);
    REQUIRE(A::finfo(*state).enabled);

    const ffu_env before = ffu_snapshot(*state);
    const std::vector<ffu_tensor> Tn_before = A::Tn(*state);

    REQUIRE_NOTHROW(state->full_update(full_updates[0]));

    const ffu_env after = ffu_snapshot(*state);

    // Premise: the bond update really changed the two sites, otherwise every
    // move below would be a no-op and both checks would be hollow.
    {
      const std::vector<ffu_tensor>& Tn_after = A::Tn(*state);
      double moved = 0.0;
      for (int site = 0; site < lattice.N_UNIT; ++site) {
        moved = std::max(
            moved, ffu_diff(Tn_before[site], Tn_after[site]) < 0.0
                       ? 1.0
                       : mptensor::max_abs(Tn_after[site] - Tn_before[site]));
      }
      INFO(label << ": the bond update moved Tn by " << moved);
      REQUIRE(moved > 1.0e-10);
    }

    const auto moves = ffu_expected_moves(source, source_leg, lattice);
    std::set<std::pair<int, int>> expected;
    for (const auto& m : moves) {
      const auto slots = ffu_expected_changed(m.first, m.second, lattice);
      expected.insert(slots.begin(), slots.end());
    }
    // Premise of check (1): the two moves rewrite disjoint kinds, so the two
    // of them together cannot look like one of them. (True for every leg by
    // the table of contract section 2; asserted so that a future edit of the
    // table above cannot silently weaken the check.)
    {
      const auto first =
          ffu_expected_changed(moves[0].first, moves[0].second, lattice);
      const auto second =
          ffu_expected_changed(moves[1].first, moves[1].second, lattice);
      REQUIRE(first.size() + second.size() == expected.size());
    }

    const std::array<std::set<int>, 2> move_kinds{
        ffu_move_kinds(moves[0].first), ffu_move_kinds(moves[1].first)};
    ffu_check_structure(before, after, expected, move_kinds, lattice, label);

    // The reduced single-layer tensors of the state AFTER the bond: contract
    // section 2 says this is what the fermionic CTM path feeds its moves
    // (build_reduced_density_tensors -> Make_single_tensor_density ->
    // *_move_single, exactly what iTPS::update_CTM() does in fermion mode).
    const std::vector<ffu_tensor> Tn_single =
        core::Make_single_tensor_density(fermion::build_reduced_density_tensors(
            A::Tn(*state), A::finfo(*state)));
    ffu_check_oracle(before, after, moves, Tn_single,
                     A::peps_parameters(*state), A::lattice(*state), label);

    state.reset();
    std::error_code ec;
    std::filesystem::remove_all(params.outdir, ec);
  }
}

TEST_CASE(
    "fermion fast full update S-3 control: with fastfullupdate = false the "
    "whole CTM is re-converged instead") {
  // The control for the test above: the same bond on the same state with
  // Full_Use_FastFullUpdate = false has to touch slots that no pair of
  // partial moves touches. Without it, a "fast" path that quietly fell back
  // to update_CTM() could not be told from one that ran the two moves - the
  // set of expected slots would just be a subset of everything.
  using namespace tenes;
  using namespace tenes::itps;

  const SquareLattice lattice = ffu_make_lattice();
  const ffu_tensor gate = make_free_fermion_gate(ffu_tau);
  const auto simple_updates = ffu_make_simple_updates(lattice, gate);
  const int source = lattice.index(1, 1);
  const int source_leg = 2;

  std::vector<EvolutionOperator<ffu_tensor>> full_updates{
      make_twosite_EvolutionOperator<ffu_tensor>(source, source_leg, 0, gate)};
  const PEPS_Parameters params = ffu_make_params(
      lattice, /*fast=*/false, "output_test_fermion_fast_full_update_control");
  REQUIRE(params.Full_Use_FastFullUpdate == false);

  std::string warnings;
  auto state =
      ffu_build_state(lattice, params, simple_updates, full_updates, warnings);
  INFO("warnings while preparing the state: " << warnings);
  REQUIRE(warnings.find("CTM did not converge") == std::string::npos);

  const ffu_env before = ffu_snapshot(*state);
  {
    fgd_stream_capture out(std::cout);
    fgd_stream_capture err(std::cerr);
    REQUIRE_NOTHROW(state->full_update(full_updates[0]));
  }
  const ffu_env after = ffu_snapshot(*state);

  const auto moves = ffu_expected_moves(source, source_leg, lattice);
  std::set<std::pair<int, int>> expected;
  for (const auto& m : moves) {
    const auto slots = ffu_expected_changed(m.first, m.second, lattice);
    expected.insert(slots.begin(), slots.end());
  }

  int outside = 0;
  for (int kind = 0; kind < FFU_NKIND; ++kind) {
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      if (expected.count({kind, site}) != 0) continue;
      if (ffu_diff(before.slot[kind][site], after.slot[kind][site]) != 0.0) {
        ++outside;
      }
    }
  }
  INFO("slots outside the two fast moves that the non-fast path changed: "
       << outside << " of "
       << (FFU_NKIND * lattice.N_UNIT - static_cast<int>(expected.size())));
  CHECK(outside > 0);

  state.reset();
  std::error_code ec;
  std::filesystem::remove_all(params.outdir, ec);
}
