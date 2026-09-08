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

// ===== Diagnosing a failed fermion full-update decomposition ================
//
// When a graded decomposition fails, tenes::fermion::svd() and qr() return
// the FIRST nonzero LAPACK info of the two block factorizations, so the
// caller cannot tell the even sector from the odd one - and the odd sector
// is the fermion-specific half, the one that goes empty or degenerate in
// the known failure modes. These tests pin the diagnostics channel that
// carries the per-sector facts out to the error message, the message
// builder that turns them into something a user can act on, and the bond /
// step context the sweep loop adds on top.
//
// Included into the test_fermion_layer TU after fermion/full_update_bond.cpp,
// whose fue_make_case / fub_gate_plain / fub_run_fermion helpers drive the
// one guard that a test can actually reach end to end.

#include <complex>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>

#include "../../src/exception.hpp"
#include "../../src/iTPS/full_update_diagnostics.hpp"

namespace {

using tenes::fermion::decomposition_diagnostics;

//! Substring test spelled out so a failure prints the whole haystack.
bool fdd_contains(const std::string& haystack, const std::string& needle) {
  return haystack.find(needle) != std::string::npos;
}

}  // namespace

// ---- the diagnostics channel of the graded decompositions ------------------

TEST_CASE("graded SVD reports the shape and the LAPACK info of each sector") {
  // Rows are (even, odd, even) and columns (even, odd, odd), so the even
  // block is 2x1 and the odd block is 1x2. Nothing else in the codebase
  // exposes that split.
  tenes::fermion::leg_parities p{{false, true, false}, {false, true, true}};
  ft a = make_even_ft(mptensor::Shape(3, 3), p, 61);
  ft u, vt;
  std::vector<double> s;
  decomposition_diagnostics diag;
  const int info = tenes::fermion::svd(a, mptensor::Axes(0), mptensor::Axes(1),
                                       u, s, vt, &diag);
  CHECK(info == 0);
  CHECK(diag.row_even == 2);
  CHECK(diag.col_even == 1);
  CHECK(diag.row_odd == 1);
  CHECK(diag.col_odd == 2);
  CHECK(diag.info_even == 0);
  CHECK(diag.info_odd == 0);
  CHECK(diag.from_svd == true);
  // A successful decomposition must not pay for the input scan: it is an
  // extra pass over the matrix (and, under MPI, an extra reduction) on a
  // path that runs once per bond per full-update step.
  CHECK(diag.scanned == false);
}

TEST_CASE("graded QR reports the shape and the LAPACK info of each sector") {
  tenes::fermion::leg_parities p{{false, true, false}, {false, true, true}};
  ft a = make_even_ft(mptensor::Shape(3, 3), p, 62);
  ft q, r;
  decomposition_diagnostics diag;
  const int info =
      tenes::fermion::qr(a, mptensor::Axes(0), mptensor::Axes(1), q, r, &diag);
  CHECK(info == 0);
  CHECK(diag.row_even == 2);
  CHECK(diag.col_even == 1);
  CHECK(diag.row_odd == 1);
  CHECK(diag.col_odd == 2);
  CHECK(diag.info_even == 0);
  CHECK(diag.info_odd == 0);
  CHECK(diag.from_svd == false);
  CHECK(diag.scanned == false);
}

TEST_CASE("truncated graded SVD forwards the diagnostics of the full one") {
  tenes::fermion::leg_parities p{{false, false, true}, {false, false, true}};
  ft a{tenes::real_tensor(mptensor::Shape(3, 3)), p};
  a.t.set_value(mptensor::Index(0, 0), 3.0);
  a.t.set_value(mptensor::Index(1, 1), 1.0);
  a.t.set_value(mptensor::Index(2, 2), 2.0);
  ft u, vt;
  std::vector<double> s;
  decomposition_diagnostics diag;
  const int info = tenes::fermion::svd_trunc(
      a, mptensor::Axes(0), mptensor::Axes(1), u, s, vt, 2, &diag);
  CHECK(info == 0);
  CHECK(diag.row_even == 2);
  CHECK(diag.col_even == 2);
  CHECK(diag.row_odd == 1);
  CHECK(diag.col_odd == 1);
  CHECK(diag.from_svd == true);
}

// ---- the input scan, and the condition that gates it ----------------------

TEST_CASE("the input scan runs only after a sector has failed") {
  tenes::real_tensor m(mptensor::Shape(2, 2));
  m.set_value(mptensor::Index(0, 0), 4.0);
  m.set_value(mptensor::Index(0, 1), -2.0);
  m.set_value(mptensor::Index(1, 1), std::numeric_limits<double>::quiet_NaN());

  decomposition_diagnostics clean;
  clean.note_failed_block(m, "even");
  CHECK(clean.scanned == false);
  CHECK(clean.has_nonfinite == false);
  CHECK(clean.max_abs == 0.0);

  decomposition_diagnostics failed;
  failed.info_odd = 1;
  failed.note_failed_block(m, "odd");
  CHECK(failed.scanned == true);
  CHECK(failed.has_nonfinite == true);
  // The NaN must not poison the magnitude that gets reported alongside it.
  CHECK(failed.max_abs == doctest::Approx(4.0));
}

// ---- describe(): the sentence that names the failing sector ---------------

TEST_CASE("describe() names each sector, its shape, and its LAPACK info") {
  decomposition_diagnostics d;
  d.row_even = 9;
  d.col_even = 9;
  d.row_odd = 7;
  d.col_odd = 7;
  d.info_odd = 3;
  const std::string text = d.describe();
  INFO(text);
  CHECK(fdd_contains(text, "even sector 9x9"));
  CHECK(fdd_contains(text, "odd sector 7x7"));
  CHECK(fdd_contains(text, "info=3"));
  CHECK(fdd_contains(text, "did not converge"));
}

TEST_CASE("describe() spells out a negative LAPACK info as a bad argument") {
  decomposition_diagnostics d;
  d.row_even = 4;
  d.col_even = 4;
  d.row_odd = 4;
  d.col_odd = 4;
  d.info_even = -6;
  const std::string text = d.describe();
  INFO(text);
  CHECK(fdd_contains(text, "info=-6"));
  CHECK(fdd_contains(text, "illegal argument 6"));
}

TEST_CASE("describe() calls an empty sector empty rather than 0x0") {
  // An empty odd sector is the shape the solver already warns about
  // ("kept an empty parity sector"); the message must connect to it.
  decomposition_diagnostics d;
  d.row_even = 4;
  d.col_even = 4;
  d.row_odd = 0;
  d.col_odd = 3;
  const std::string text = d.describe();
  INFO(text);
  CHECK(fdd_contains(text, "odd sector empty"));
}

TEST_CASE("describe() reports a non-finite input once the scan has run") {
  decomposition_diagnostics d;
  d.row_even = 2;
  d.col_even = 2;
  d.row_odd = 2;
  d.col_odd = 2;
  d.info_even = 1;
  d.scanned = true;
  d.has_nonfinite = true;
  d.max_abs = 1.25e3;
  const std::string text = d.describe();
  INFO(text);
  CHECK(fdd_contains(text, "non-finite"));
  CHECK(fdd_contains(text, "max|.|="));
}

TEST_CASE("describe() says the input was finite when the scan found nothing") {
  decomposition_diagnostics d;
  d.row_even = 2;
  d.col_even = 2;
  d.row_odd = 2;
  d.col_odd = 2;
  d.info_even = 1;
  d.scanned = true;
  d.has_nonfinite = false;
  d.max_abs = 3.5;
  const std::string text = d.describe();
  INFO(text);
  CHECK(fdd_contains(text, "all elements finite"));
}

// ---- the user-facing message ---------------------------------------------

TEST_CASE(
    "the full-update failure message names the decomposition, the sector, "
    "and the parameters to change") {
  decomposition_diagnostics d;
  d.row_even = 9;
  d.col_even = 9;
  d.row_odd = 7;
  d.col_odd = 7;
  d.info_odd = 3;
  const std::string msg = tenes::itps::fermion_full_update_failure_message(
      "balancing SVD", d.describe());
  INFO(msg);
  CHECK(fdd_contains(msg, "balancing SVD"));
  CHECK(fdd_contains(msg, "odd sector 7x7"));
  // The two upstream warnings that actually precede this failure in the
  // solver's own output, so the user knows where to look.
  CHECK(fdd_contains(msg, "CTM did not converge"));
  CHECK(fdd_contains(msg, "empty parity sector"));
  // The knobs, spelled exactly as input.toml spells them.
  CHECK(fdd_contains(msg, "parameter.ctm.iteration_max"));
  CHECK(fdd_contains(msg, "parameter.ctm.dimension"));
  CHECK(fdd_contains(msg, "parameter.simple_update.num_step"));
  CHECK(fdd_contains(msg, "parameter.full_update.tau"));
  // The recovery answer: nothing was saved, so arrange for a next time.
  CHECK(fdd_contains(msg, "parameter.general.tensor_save"));
  CHECK(fdd_contains(msg, "parameter.general.tensor_load"));
}

// ---- the bond and step context the sweep loop adds ------------------------

TEST_CASE("the bond context names both sites, the direction, and the step") {
  // Leg 2 is the right neighbour (SquareLattice::right), and the step is
  // reported one-based against the configured total.
  const std::string c = tenes::itps::full_update_bond_context(3, 4, 2, 6, 50);
  INFO(c);
  CHECK(fdd_contains(c, "site 3"));
  CHECK(fdd_contains(c, "site 4"));
  CHECK(fdd_contains(c, "right"));
  CHECK(fdd_contains(c, "step 7/50"));
}

TEST_CASE("the bond context names the other three legs correctly") {
  CHECK(fdd_contains(tenes::itps::full_update_bond_context(1, 0, 0, 0, 1),
                     "left"));
  CHECK(fdd_contains(tenes::itps::full_update_bond_context(1, 5, 1, 0, 1),
                     "top"));
  CHECK(fdd_contains(tenes::itps::full_update_bond_context(1, 5, 3, 0, 1),
                     "bottom"));
}

// ---- the guard as the solver actually reaches it --------------------------

namespace {

/*! Run the fermionic bond update with an all-zero gate. The two-site state
 *  it produces is exactly zero, every singular value with it, and the
 *  weight normalization has nothing to divide by - the one guard in
 *  Full_update_bond_fermion a test can reach without forcing LAPACK to
 *  fail. Returns the message; fails the test if nothing is thrown.
 */
template <class tensor>
std::string fdd_zero_gate_message(char dir, int d, int seed) {
  const fue_case<tensor> c =
      fue_make_case<tensor>(dir, fue_full_geom(dir), "eo", d, d, 2, seed, false,
                            std::string("zero gate ") + dir);
  const fg_ftensor<tensor> Tn1 = c.sites[c.a];
  const fg_ftensor<tensor> Tn2 = c.sites[c.b];
  // Premise: the inputs themselves are not degenerate, so the failure is
  // the gate's doing and not an accident of the fixture.
  REQUIRE(fgf::max_abs(Tn1) > 0.0);
  REQUIRE(fgf::max_abs(Tn2) > 0.0);

  const tensor plain = fub_gate_plain<tensor>(d, 0.0, 0.0);
  REQUIRE(mptensor::max_abs(plain) == 0.0);
  const fg_ftensor<tensor> gate =
      fgf::wrap_twosite_gate(plain, Tn1.parity[4], Tn2.parity[4]);

  const tenes::itps::PEPS_Parameters params = fub_exact_params(false);
  fg_ftensor<tensor> Tn1_new, Tn2_new;
  try {
    fub_run_fermion(c, Tn1, Tn2, gate, params, Tn1_new, Tn2_new);
  } catch (const tenes::runtime_error& e) {
    return std::string(e.what());
  }
  FAIL("Full_update_bond_fermion accepted an all-zero gate");
  return std::string();
}

}  // namespace

TEST_CASE(
    "a degenerate two-site state is reported as a tenes error with the "
    "recovery advice attached") {
  const std::string msg =
      fdd_zero_gate_message<tenes::real_tensor>('h', 2, 4700);
  INFO(msg);
  // Not "[UNEXPECTED ERROR]": tenes::runtime_error is what main.cpp prints
  // under "[ERROR]", i.e. as a condition of the run and not as a bug.
  CHECK(fdd_contains(msg, "fermion full update"));
  CHECK(fdd_contains(msg, "parameter.ctm.iteration_max"));
  CHECK(fdd_contains(msg, "parameter.general.tensor_save"));
}

// ---- the environment guards of the same full update -----------------------
//
// build_full_update_environment() refuses a window it cannot use, and its
// two numerical refusals are the ones the diagnostic logs of this branch
// show most often. They are conditions of the run like the decompositions
// above, so they must reach the user the same way: as tenes::runtime_error
// (printed under "[ERROR]", not "[UNEXPECTED ERROR]"), naming input.toml
// keys rather than PEPS_Parameters field names, and picking up the bond and
// step context that iTPS::full_update() appends.

TEST_CASE(
    "a parity-contaminated environment is refused as a tenes error naming "
    "the input.toml keys") {
  const fue_case<tenes::real_tensor> c = fue_make_case<tenes::real_tensor>(
      'h', fue_full_geom('h'), "eo", 2, 2, 2, 1400, false, "contaminated env");
  // Control: the clean window is accepted, so "throws" is a statement about
  // the contamination and not about the fixture.
  REQUIRE(fue_build_N(c.env, c.QA, c.QB, c.dir_e).forbidden_ratio <= 1.0e-10);

  const std::vector<tenes::real_tensor> bad =
      fue_contaminate_env(c.env, 4, fg_pv("eo"), 81);
  try {
    fue_build_N(bad, c.QA, c.QB, c.dir_e);
  } catch (const tenes::runtime_error& e) {
    const std::string msg = e.what();
    INFO(msg);
    CHECK(fdd_contains(msg, "parameter.ctm.iteration_max"));
    CHECK(fdd_contains(msg, "parameter.ctm.convergence_epsilon"));
    return;
  }
  FAIL("build_full_update_environment accepted a contaminated environment");
}

TEST_CASE(
    "a window with no norm left is refused as a tenes error naming the "
    "input.toml keys") {
  const fue_case<tenes::real_tensor> c = fue_make_case<tenes::real_tensor>(
      'h', fue_full_geom('h'), "eo", 2, 2, 2, 1401, false, "empty window");
  REQUIRE(fue_build_N(c.env, c.QA, c.QB, c.dir_e).forbidden_ratio <= 1.0e-10);

  // Zeroing one edge tensor leaves the parity ledgers intact - so the
  // forbidden-block guard has nothing to complain about - but collapses the
  // window norm, which is the shape a diverged state arrives in.
  std::vector<tenes::real_tensor> bad = c.env;
  bad[4] = tenes::real_tensor(bad[4].shape());
  REQUIRE(mptensor::max_abs(bad[4]) == 0.0);
  try {
    fue_build_N(bad, c.QA, c.QB, c.dir_e);
  } catch (const tenes::runtime_error& e) {
    const std::string msg = e.what();
    INFO(msg);
    CHECK(fdd_contains(msg, "parameter."));
    return;
  }
  FAIL("build_full_update_environment accepted a window with no norm");
}

// ===== What a real failure taught the message ==============================
//
// A run on oneAPI 2022.2.1 (icpx + MKL) failed the balancing SVD with
//
//   even sector 2x2 info=3 (LAPACK: did not converge), odd sector 1x1 info=0;
//   input max|.|=0.975684, all elements finite
//
// and every word of that was misleading. dgesvd's positive info counts the
// superdiagonals of an intermediate bidiagonal form that did not converge;
// a 2x2 block has one, so 3 cannot come out of a conforming LAPACK whatever
// the matrix held - the message pointed at parameter.ctm.* when the lead was
// the library. The three cases below pin the corrections: name an
// out-of-range info as such and stop recommending solver parameters for it;
// decide finiteness from the bit pattern, since the toolchain that produced
// this report defaults to -fp-model=fast, under which std::isfinite folds to
// true and "all elements finite" is worth nothing; and print a small failing
// block, because four numbers at full precision turn "it failed" into a
// standalone LAPACK reproducer.

// ---- 1. an info outside the LAPACK contract -------------------------------

TEST_CASE("describe() names an info that LAPACK cannot have returned") {
  decomposition_diagnostics d;
  d.row_even = 2;
  d.col_even = 2;
  d.row_odd = 1;
  d.col_odd = 1;
  // 9, not 3: under MPI, MIN(M,N)+1 = 3 is a legal ScaLAPACK code (see
  // below), so the out-of-range case needs a value neither convention can
  // justify.
  d.info_even = 9;
  const std::string text = d.describe();
  INFO(text);
  CHECK(fdd_contains(text, "out of range"));
  // The bound must be stated, so the reader can check the arithmetic.
  CHECK(fdd_contains(text, "at most 2"));
  CHECK(fdd_contains(text, "2x2"));
  // The old wording asserted a convergence failure that cannot have happened.
  CHECK_FALSE(fdd_contains(text, "did not converge"));
}

TEST_CASE(
    "describe() still reports an in-range info as a convergence failure") {
  decomposition_diagnostics d;
  d.row_even = 4;
  d.col_even = 4;
  d.row_odd = 2;
  d.col_odd = 2;
  d.info_even = 3;
  const std::string text = d.describe();
  INFO(text);
  CHECK(fdd_contains(text, "did not converge"));
  CHECK_FALSE(fdd_contains(text, "out of range"));
}

TEST_CASE("suspect_library() is what separates the two") {
  decomposition_diagnostics bad;
  bad.row_even = 2;
  bad.col_even = 2;
  bad.info_even = 9;
  CHECK(bad.suspect_library());

  decomposition_diagnostics ordinary;
  ordinary.row_even = 4;
  ordinary.col_even = 4;
  ordinary.info_even = 3;
  CHECK_FALSE(ordinary.suspect_library());

  decomposition_diagnostics clean;
  clean.row_even = 2;
  clean.col_even = 2;
  CHECK_FALSE(clean.suspect_library());
}

TEST_CASE(
    "the failure message stops recommending solver parameters once the "
    "library is the suspect") {
  const std::string msg = tenes::itps::fermion_full_update_failure_message(
      "balancing SVD", "even sector 2x2 info=9 (LAPACK: out of range)",
      tenes::itps::full_update_failure_lead::library);
  INFO(msg);
  CHECK(fdd_contains(msg, "LAPACK"));
  // The knobs are still listed - they are what to try if the library turns
  // out to be innocent - but the message must say so first.
  CHECK(fdd_contains(msg, "unlikely to help"));
  CHECK(fdd_contains(msg, "parameter.ctm.iteration_max"));
}

TEST_CASE(
    "MIN(M,N)+1 is ScaLAPACK's heterogeneity code, not a broken library") {
  // Reference ScaLAPACK, SRC/pdgesvd.f, on INFO:
  //   "> 0: if DBDSQR did not converge. If INFO = MIN(M,N) + 1, then
  //    PDGESVD has detected heterogeneity by finding that eigenvalues were
  //    not identical across the process grid."
  // A real run met exactly this: info=3 on a 2x2 under mpiexec -np 4, gone
  // at -np 1, with the block's two singular values agreeing to 2e-5. A
  // serial build never calls pdgesvd, so there the same value is out of
  // range for LAPACK's dgesvd and means what it meant before.
  decomposition_diagnostics d;
  d.row_even = 2;
  d.col_even = 2;
  d.row_odd = 1;
  d.col_odd = 1;
  d.info_even = 3;
  d.from_svd = true;
  const std::string text = d.describe();
  INFO(text);
#ifdef _NO_MPI
  CHECK(fdd_contains(text, "out of range"));
  CHECK(d.suspect_library());
  CHECK_FALSE(d.grid_heterogeneity());
#else
  CHECK(fdd_contains(text, "process grid"));
  CHECK(fdd_contains(text, "heterogeneity"));
  CHECK(d.grid_heterogeneity());
  CHECK_FALSE(d.suspect_library());
#endif
}

TEST_CASE("only an SVD can report ScaLAPACK's heterogeneity code") {
  // pdgeqrf documents INFO as 0 or -i only: it has no MIN(M,N)+1 code, so
  // a positive info from a QR is out of spec even under MPI. Reading the
  // SVD's convention into it would be the same mislabelling this change
  // exists to remove.
  decomposition_diagnostics from_qr;
  from_qr.row_even = 2;
  from_qr.col_even = 2;
  from_qr.info_even = 3;
  CHECK_FALSE(from_qr.grid_heterogeneity());
  CHECK(from_qr.suspect_library());

  decomposition_diagnostics from_svd;
  from_svd.row_even = 2;
  from_svd.col_even = 2;
  from_svd.info_even = 3;
  from_svd.from_svd = true;
#ifdef _NO_MPI
  CHECK_FALSE(from_svd.grid_heterogeneity());
#else
  CHECK(from_svd.grid_heterogeneity());
  CHECK_FALSE(from_svd.suspect_library());
#endif
}

TEST_CASE(
    "the failure message blames the process grid, not the state or the "
    "library, for a heterogeneity code") {
  const std::string msg = tenes::itps::fermion_full_update_failure_message(
      "balancing SVD", "even sector 2x2 info=3 (ScaLAPACK: heterogeneity)",
      tenes::itps::full_update_failure_lead::process_grid);
  INFO(msg);
  CHECK(fdd_contains(msg, "process grid"));
  // The one thing that is known to work, because it was measured.
  CHECK(fdd_contains(msg, "fewer MPI processes"));
  CHECK(fdd_contains(msg, "unlikely to help"));
  // Not the state, and not the library.
  CHECK_FALSE(fdd_contains(msg, "CTM did not converge"));
  CHECK_FALSE(fdd_contains(msg, "LAPACK and BLAS this binary"));
}

TEST_CASE("the failure message leads with the parameters when it should") {
  const std::string msg = tenes::itps::fermion_full_update_failure_message(
      "balancing SVD", "even sector 4x4 info=3 (LAPACK: did not converge)");
  INFO(msg);
  CHECK_FALSE(fdd_contains(msg, "unlikely to help"));
  CHECK(fdd_contains(msg, "parameter.ctm.iteration_max"));
}

// ---- 3. the failing block itself ------------------------------------------

TEST_CASE(
    "a small failing block is printed at a precision that reproduces it") {
  tenes::real_tensor block(mptensor::Shape(2, 2));
  block.set_value(mptensor::Index(0, 0), 0.1);
  block.set_value(mptensor::Index(0, 1), -0.25);
  block.set_value(mptensor::Index(1, 0), 0.5);
  block.set_value(mptensor::Index(1, 1), 2.0);

  decomposition_diagnostics d;
  d.row_even = 2;
  d.col_even = 2;
  d.row_odd = 1;
  d.col_odd = 1;
  d.info_even = 1;
  d.note_failed_block(block, "even");

  const std::string text = d.describe();
  INFO(text);
  CHECK(fdd_contains(text, "even block (2x2, row-major)"));
  // 17 significant digits, or the numbers cannot be fed back to LAPACK.
  CHECK(fdd_contains(text, "0.10000000000000001"));
  CHECK(fdd_contains(text, "-0.25"));
  CHECK(fdd_contains(text, "0.5"));
  CHECK(fdd_contains(text, "2"));
}

TEST_CASE("the odd block is named as such when the odd block failed") {
  tenes::real_tensor block(mptensor::Shape(1, 1));
  block.set_value(mptensor::Index(0, 0), 7.5);

  decomposition_diagnostics d;
  d.row_even = 2;
  d.col_even = 2;
  d.row_odd = 1;
  d.col_odd = 1;
  d.info_odd = 1;
  d.note_failed_block(block, "odd");

  const std::string text = d.describe();
  INFO(text);
  CHECK(fdd_contains(text, "odd block (1x1, row-major)"));
  CHECK(fdd_contains(text, "7.5"));
  CHECK_FALSE(fdd_contains(text, "even block"));
}

TEST_CASE("a block too large to read is summarized rather than printed") {
  // 5x5 is 25 elements; the dump is for blocks a human can retype into a
  // LAPACK call, not for pouring a matrix into a log.
  tenes::real_tensor block(mptensor::Shape(5, 5));
  block.set_value(mptensor::Index(0, 0), 1.0);

  decomposition_diagnostics d;
  d.row_even = 5;
  d.col_even = 5;
  d.info_even = 2;
  d.note_failed_block(block, "even");

  const std::string text = d.describe();
  INFO(text);
  CHECK(d.scanned == true);
  CHECK_FALSE(fdd_contains(text, "even block"));
  // The scan still reports what it always did.
  CHECK(fdd_contains(text, "max|.|="));
}

// ===== The scan must describe what LAPACK was handed =======================
//
// fermion::svd() and qr() do not decompose the parity-sorted matrix; they
// decompose a slice() copy of each diagonal block. The first version of
// these diagnostics scanned and dumped the SOURCE matrix instead, so a
// block that got corrupted on its way into LAPACK would be reported as
// healthy - which is exactly the question left open by a oneAPI 2022.2.1
// failure whose dumped 2x2 decomposed perfectly when fed to the same
// machine's dgesvd on its own. The recorder therefore takes the block, and
// the call sites sit next to the LAPACK call so the two cannot drift apart.

TEST_CASE("the recorder describes the tensor it is given, not its source") {
  // A 3x3 whose odd corner is huge and whose even block is small: if the
  // scan ever reads the whole matrix again, max_abs gives it away.
  tenes::real_tensor block(mptensor::Shape(2, 2));
  block.set_value(mptensor::Index(0, 0), 0.25);
  block.set_value(mptensor::Index(1, 1), -0.5);

  decomposition_diagnostics d;
  d.row_even = 2;
  d.col_even = 2;
  d.row_odd = 1;
  d.col_odd = 1;
  d.info_even = 3;
  d.note_failed_block(block, "even");

  CHECK(d.scanned == true);
  CHECK(d.max_abs == doctest::Approx(0.5));
  const std::string text = d.describe();
  INFO(text);
  CHECK(fdd_contains(text, "even block (2x2"));
  CHECK(fdd_contains(text, "0.25"));
}

TEST_CASE("the recorder finds a non-finite element in the block it is given") {
  tenes::real_tensor block(mptensor::Shape(2, 2));
  block.set_value(mptensor::Index(0, 0), 1.0);
  block.set_value(mptensor::Index(1, 1),
                  std::numeric_limits<double>::quiet_NaN());

  decomposition_diagnostics d;
  d.row_odd = 2;
  d.col_odd = 2;
  d.info_odd = 1;
  d.note_failed_block(block, "odd");

  CHECK(d.has_nonfinite == true);
  CHECK(d.max_abs == doctest::Approx(1.0));
  const std::string text = d.describe();
  INFO(text);
  CHECK(fdd_contains(text, "non-finite"));
  CHECK(fdd_contains(text, "odd block (2x2"));
}

TEST_CASE("the recorder does nothing for a decomposition that succeeded") {
  tenes::real_tensor block(mptensor::Shape(2, 2));
  block.set_value(mptensor::Index(0, 0), 3.0);

  decomposition_diagnostics d;
  d.row_even = 2;
  d.col_even = 2;
  d.note_failed_block(block, "even");

  CHECK(d.scanned == false);
  CHECK(d.max_abs == 0.0);
  CHECK(d.block_dump.empty());
}

TEST_CASE("a graded SVD that succeeds records no block at all") {
  tenes::fermion::leg_parities p{{false, true, false}, {false, true, true}};
  ft a = make_even_ft(mptensor::Shape(3, 3), p, 63);
  ft u, vt;
  std::vector<double> s;
  decomposition_diagnostics diag;
  REQUIRE(tenes::fermion::svd(a, mptensor::Axes(0), mptensor::Axes(1), u, s, vt,
                              &diag) == 0);
  CHECK(diag.block_dump.empty());
  CHECK(diag.scanned == false);
}
