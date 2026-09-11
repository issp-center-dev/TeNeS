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

// ===== Scalars reduced from a distributed tensor must be global ===========
//
// A fermionic Hubbard state measured through the same binary gave the same
// one-site observables at 1, 2 and 4 MPI ranks and a hopping that moved by
// 0.023 between 1 and 2.  Every two-site observable moved; no one-site one
// did.  The two-site CTM contraction ends in detail::trace_boundary_pairs(),
// which sums the elements it holds locally and returns -- a partial sum on
// every rank but one.
//
// ctest runs the MPI build with one rank, where a local sum is the global
// sum, so this file is registered a second time under two ranks.  These
// cases are written so that they pass trivially at one rank and only bite
// with more: the reference is assembled from get_value(), which is
// collective, so both sides of each CHECK are the same on every rank.

#include "../test_fermion_common.hpp"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/mpi.hpp"

namespace {

//! The entry at a global index: a closed form, so the reference below can
//! be assembled by arithmetic alone.  mptensor's get_value() leaves the
//! value untouched on ranks that do not own the element, so a reference
//! read through it would itself be a per-rank partial sum and the check
//! would pass whatever the code under test did.
double mr_entry(const mptensor::Index& idx, std::size_t n) {
  const double x = static_cast<double>(idx[0] + n * idx[1] + n * n * idx[2] +
                                       n * n * n * idx[3]);
  return std::cos(0.37 * x) + 0.5;
}

//! A rank-4 tensor filled from mr_entry, the same on every rank.
tenes::real_tensor mr_rank4(std::size_t n) {
  tenes::real_tensor a(mptensor::Shape(n, n, n, n));
  for (std::size_t k = 0; k < a.local_size(); ++k) {
    const mptensor::Index idx = a.global_index(k);
    a.set_value(idx, mr_entry(idx, n));
  }
  return a;
}

//! The double-delta trace sum_{i,k} a[i,i,k,k], from the closed form.
double mr_reference_trace(std::size_t n) {
  double value = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t k = 0; k < n; ++k) {
      value += mr_entry(mptensor::Index(i, i, k, k), n);
    }
  }
  return value;
}

}  // namespace

TEST_CASE(
    "trace_boundary_pairs is the global double-delta trace on every rank") {
  // n = 20 puts the matricized 400 x 400 tensor across more than one
  // block-cyclic block in each direction, so at two ranks no rank holds
  // every diagonal-pair element.
  const tenes::real_tensor a = mr_rank4(20);
  const double want = mr_reference_trace(20);
  const double got = tenes::fermion::detail::trace_boundary_pairs(a);
  INFO("ranks = " << a.get_comm_size() << ", rank = " << a.get_comm_rank());
  CHECK(got == doctest::Approx(want).epsilon(1.0e-12));
}

// ===== Parity guards decide on every rank alike ===========================
//
// The simple-update guard enforce_even_parity(), the layer checks of
// doubled_pipeline_traced() and validate_block_diagonal() compared the
// largest parity-odd magnitude of the PROCESS-LOCAL slice with a threshold.
// With two or more ranks only the rank that stored the offending element
// threw; the others went on into the next collective call and the run hung.
// The fix routes them through one collective predicate,
// require_even_parity(), which reduces both the odd-sector maximum and the
// scale over the ranks.  Contract:
// docs/superpowers/specs/2026-09-12-fermion-parity-guard-collective-contract.md
//
// Each guard is called on every rank and whatever it does is caught on
// every rank; only then are the outcomes counted with a reduction and
// checked.  A guard that throws on some ranks only therefore fails the
// checks instead of leaving the other ranks in a collective that never
// completes.  For the same reason every REQUIRE below is on a reduced
// quantity, so that when it stops a test case it stops it on every rank.
//
// The offending element is picked at run time from the local slices, and
// the premise that makes a rank-local guard disagree -- exactly one rank
// stores the violation, exactly one other rank stores the element that sets
// the scale -- is counted with a reduction, not assumed.  At one rank every
// decision is the old one; at two the cases bite.

namespace {

//! Leg dimension of the fixtures.  Six index values per leg make the
//! matricizations 36 rows tall, which spreads them over more than one
//! 16-row block of the block-cyclic distribution, so at two ranks every
//! rank stores part of both parity sectors.
constexpr std::size_t pg_dim = 6;

//! The one element above 0.5 in magnitude in every fixture, so that the
//! guards' threshold 1e-10 * max(1, max_abs) is 4e-10 rather than the bare
//! 1e-10 a rank without it would compute from its own slice.
constexpr double pg_scale = 4.0;
constexpr double pg_threshold = 1.0e-10 * pg_scale;

//! A violation far above the threshold.
constexpr double pg_gross = 0.25;

//! A violation between the bare 1e-10 and the threshold 4e-10: the
//! collective decision lets it pass, a decision scaled by a slice that does
//! not hold pg_scale rejects it.
constexpr double pg_faint = 2.0e-10;

//! Parity ledger of one virtual leg; the odd values are spread over the
//! leg so that both sectors appear in every block of the distribution.
tenes::fermion::parity_vector pg_leg() {
  return {false, true, true, false, true, false};
}

//! Column-major position of a global index (first leg fastest), as a
//! double so that it can travel through allreduce_max().
double pg_flat(const mptensor::Index& idx, const mptensor::Shape& shape) {
  double flat = 0.0;
  for (std::size_t ax = shape.size(); ax-- > 0;) {
    flat = flat * static_cast<double>(shape[ax]) + static_cast<double>(idx[ax]);
  }
  return flat;
}

//! Inverse of pg_flat().
mptensor::Index pg_unflat(double flat, const mptensor::Shape& shape) {
  std::size_t rest = static_cast<std::size_t>(flat);
  mptensor::Index idx;
  idx.resize(shape.size());
  for (std::size_t ax = 0; ax < shape.size(); ++ax) {
    idx[ax] = rest % shape[ax];
    rest /= shape[ax];
  }
  return idx;
}

//! Closed-form entry of size at most 0.5 for the parity-allowed positions.
double pg_entry(const mptensor::Index& idx, const mptensor::Shape& shape) {
  return 0.5 * std::cos(0.37 * pg_flat(idx, shape));
}

//! Closed-form round-off between 0.5e-12 and 1.5e-12 for the parity-odd
//! positions: what a working update leaves behind, far below any threshold.
double pg_residue(const mptensor::Index& idx, const mptensor::Shape& shape) {
  return 1.0e-12 * (1.0 + 0.5 * std::sin(0.53 * pg_flat(idx, shape)));
}

bool pg_is_odd(const tenes::fermion::leg_parities& legs,
               const mptensor::Index& idx) {
  return tenes::fermion::count_odd(legs, idx) % 2 == 1;
}

//! Number of ranks of t's communicator on which flag holds.
int pg_ranks_where(bool flag, const tenes::real_tensor& t) {
  int count = flag ? 1 : 0;
  tenes::allreduce_sum(count, t.get_comm());
  return count;
}

//! Whether this rank stores the element at idx.  A scan of the local
//! slice rather than get_value(), whose value is only written on the owner.
bool pg_stores(const tenes::real_tensor& t, const mptensor::Index& idx) {
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    if (t.global_index(n) == idx) {
      return true;
    }
  }
  return false;
}

//! Largest |element| of this rank's slice among the positions pred selects.
template <class Pred>
double pg_local_max(const tenes::real_tensor& t, Pred pred) {
  double v = 0.0;
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    if (pred(t.global_index(n))) {
      v = std::max(v, std::abs(t[n]));
    }
  }
  return v;
}

//! The first element, in local storage order, that pred selects on the
//! lowest (from_last == false) or the highest rank that stores one.  The
//! candidates travel through a reduction, so every rank returns the same
//! index, and the REQUIRE fails on every rank alike if no rank stores one.
template <class Pred>
mptensor::Index pg_pick(const tenes::real_tensor& t, bool from_last,
                        Pred pred) {
  std::vector<double> first(static_cast<std::size_t>(t.get_comm_size()),
                            -1.0);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    if (pred(idx)) {
      first[static_cast<std::size_t>(t.get_comm_rank())] =
          pg_flat(idx, t.shape());
      break;
    }
  }
  tenes::allreduce_max(first, t.get_comm());
  double chosen = -1.0;
  for (std::size_t k = 0; k < first.size() && chosen < 0.0; ++k) {
    chosen = first[from_last ? first.size() - 1 - k : k];
  }
  REQUIRE(chosen >= 0.0);
  return pg_unflat(chosen, t.shape());
}

//! Fill t from pg_entry() where allowed(idx), from odd_fill(idx) elsewhere;
//! then put pg_scale on an allowed element of the lowest rank that stores
//! one and `violation` on a forbidden element of the highest.  Every rank
//! runs the same code and set_value() writes on the owner only.
template <class Allowed, class OddFill>
void pg_fill(tenes::real_tensor& t, Allowed allowed, OddFill odd_fill,
             double violation, mptensor::Index& scale_at,
             mptensor::Index& odd_at) {
  const mptensor::Shape shape = t.shape();
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    t.set_value(idx, allowed(idx) ? pg_entry(idx, shape) : odd_fill(idx));
  }
  scale_at = pg_pick(t, false, allowed);
  odd_at = pg_pick(t, true, [&](const mptensor::Index& idx) {
    return !allowed(idx);
  });
  t.set_value(scale_at, pg_scale);
  t.set_value(odd_at, violation);

  // Each picked element lives on exactly one rank.  With two or more ranks
  // they live on different ones: then a rank that looks only at its own
  // slice either sees the violation without the scale or the scale without
  // the violation.  A layout that replicated the tensor on every rank, or
  // kept it on one, would make the rank-local guards agree with the
  // collective ones and the cases below hollow; this stops them there.
  const bool has_scale = pg_stores(t, scale_at);
  const bool has_odd = pg_stores(t, odd_at);
  INFO("ranks = " << t.get_comm_size() << ", rank = " << t.get_comm_rank());
  REQUIRE(pg_ranks_where(has_scale, t) == 1);
  REQUIRE(pg_ranks_where(has_odd, t) == 1);
  REQUIRE(pg_ranks_where(has_scale && has_odd, t) ==
          (t.get_comm_size() == 1 ? 1 : 0));
}

//! A graded fixture: see pg_fill().
struct pg_graded {
  ft a;
  mptensor::Index scale_at;  //!< Parity even; holds pg_scale.
  mptensor::Index odd_at;    //!< Parity odd; holds the violation.
};

template <class OddFill>
pg_graded pg_make_graded(const mptensor::Shape& shape,
                         const tenes::fermion::leg_parities& legs,
                         OddFill odd_fill, double violation) {
  pg_graded g{ft{tenes::real_tensor(shape), legs}, {}, {}};
  pg_fill(
      g.a.t, [&](const mptensor::Index& idx) { return !pg_is_odd(legs, idx); },
      odd_fill, violation, g.scale_at, g.odd_at);
  return g;
}

//! Rank-4 site of legs pg_leg(): odd sector at residue level except
//! `violation` on one element (0 and no residue: parity even).
pg_graded pg_make_site(double violation, bool residue) {
  const mptensor::Shape shape(pg_dim, pg_dim, pg_dim, pg_dim);
  const tenes::fermion::leg_parities legs{pg_leg(), pg_leg(), pg_leg(),
                                          pg_leg()};
  return pg_make_graded(
      shape, legs,
      [&](const mptensor::Index& idx) {
        return residue ? pg_residue(idx, shape) : 0.0;
      },
      violation);
}

//! What one guard call did, counted over the ranks.  The counts are the
//! same on every rank, so a CHECK on them agrees everywhere; the message is
//! this rank's own (empty where it did not throw std::runtime_error).
struct pg_tally {
  int ranks = 0;        //!< Size of the communicator.
  int threw = 0;        //!< Ranks on which it threw std::runtime_error.
  int threw_other = 0;  //!< Ranks on which it threw anything else.
  std::string what;
};

template <class Call>
pg_tally pg_run(const tenes::real_tensor& t, Call call) {
  bool threw = false;
  bool threw_other = false;
  pg_tally tally;
  try {
    call();
  } catch (const std::runtime_error& e) {
    threw = true;
    tally.what = e.what();
  } catch (...) {
    threw_other = true;
  }
  tally.ranks = t.get_comm_size();
  tally.threw = pg_ranks_where(threw, t);
  tally.threw_other = pg_ranks_where(threw_other, t);
  return tally;
}

bool pg_starts_with(const std::string& text, const std::string& prefix) {
  return text.compare(0, prefix.size(), prefix) == 0;
}

//! Whether some number written in text equals value to five significant
//! figures.  The contract fixes which numbers the message carries, not how
//! it formats them, so every numeric token is parsed.
bool pg_mentions(const std::string& text, double value) {
  const char* const s = text.c_str();
  std::size_t i = 0;
  while (i < text.size()) {
    const bool starts =
        std::isdigit(static_cast<unsigned char>(s[i])) ||
        (s[i] == '.' && i + 1 < text.size() &&
         std::isdigit(static_cast<unsigned char>(s[i + 1])));
    char* end = nullptr;
    const double v = starts ? std::strtod(s + i, &end) : 0.0;
    if (!starts || end == s + i) {
      ++i;
      continue;
    }
    if (std::abs(v - value) <= 1.0e-5 * std::abs(value)) {
      return true;
    }
    i = static_cast<std::size_t>(end - s);
  }
  return false;
}

}  // namespace

// ---- require_even_parity -------------------------------------------------

TEST_CASE(
    "require_even_parity: an odd element stored by one rank throws on every "
    "rank") {
  const pg_graded site = pg_make_site(pg_gross, true);
  const ft& a = site.a;
  const auto odd = [&](const mptensor::Index& idx) {
    return pg_is_odd(a.parity, idx);
  };
  // The premise a rank-local guard would trip over: one rank sees it.
  REQUIRE(pg_ranks_where(pg_local_max(a.t, odd) > pg_threshold, a.t) == 1);

  const char* const context = "mpi_reduction parity guard";
  const pg_tally tally =
      pg_run(a.t, [&] { tenes::fermion::require_even_parity(a, context); });
  INFO("ranks = " << tally.ranks << ", rank = " << a.t.get_comm_rank()
                  << ", message: " << tally.what);
  CHECK(tally.threw_other == 0);
  CHECK(tally.threw == tally.ranks);
  if (!tally.what.empty()) {
    CHECK(pg_starts_with(tally.what, context));
    CHECK(tally.what.find("not parity even") != std::string::npos);
    // Every rank reports the GLOBAL odd-sector maximum and the threshold
    // from the GLOBAL scale, though only one rank stores each of them.
    CHECK(pg_mentions(tally.what, pg_gross));
    CHECK(pg_mentions(tally.what, pg_threshold));
  }
}

TEST_CASE(
    "require_even_parity: below the threshold, or parity even, no rank "
    "throws") {
  const char* const context = "mpi_reduction parity guard";
  {
    const pg_graded site = pg_make_site(pg_faint, true);
    const ft& a = site.a;
    // The premise a guard scaled by its own slice would trip over: every
    // rank but one sees max_abs <= 1, hence a threshold of 1e-10, which the
    // reduced violation pg_faint exceeds.
    REQUIRE(pg_ranks_where(pg_local_max(a.t,
                                        [](const mptensor::Index&) {
                                          return true;
                                        }) > 1.0,
                           a.t) == 1);
    const pg_tally tally =
        pg_run(a.t, [&] { tenes::fermion::require_even_parity(a, context); });
    INFO("below the threshold; ranks = " << tally.ranks
                                         << ", rank = " << a.t.get_comm_rank()
                                         << ", message: " << tally.what);
    CHECK(tally.threw_other == 0);
    CHECK(tally.threw == 0);
  }
  {
    const pg_graded site = pg_make_site(0.0, false);
    const ft& a = site.a;
    const pg_tally tally =
        pg_run(a.t, [&] { tenes::fermion::require_even_parity(a, context); });
    INFO("parity even; ranks = " << tally.ranks
                                 << ", rank = " << a.t.get_comm_rank()
                                 << ", message: " << tally.what);
    CHECK(tally.threw_other == 0);
    CHECK(tally.threw == 0);
  }
}

// ---- enforce_even_parity (the simple-update guard) -----------------------

TEST_CASE(
    "enforce_even_parity: an odd element stored by one rank throws on every "
    "rank") {
  pg_graded site = pg_make_site(pg_gross, true);
  ft& a = site.a;
  const auto odd = [&](const mptensor::Index& idx) {
    return pg_is_odd(a.parity, idx);
  };
  REQUIRE(pg_ranks_where(pg_local_max(a.t, odd) > pg_threshold, a.t) == 1);

  const pg_tally tally =
      pg_run(a.t, [&] { tenes::fermion::enforce_even_parity(a); });
  INFO("ranks = " << tally.ranks << ", rank = " << a.t.get_comm_rank()
                  << ", message: " << tally.what);
  CHECK(tally.threw_other == 0);
  CHECK(tally.threw == tally.ranks);
  if (!tally.what.empty()) {
    CHECK(pg_mentions(tally.what, pg_gross));
    CHECK(pg_mentions(tally.what, pg_threshold));
  }
}

TEST_CASE(
    "enforce_even_parity: below the threshold no rank throws, and each rank "
    "clips its own odd elements") {
  for (const bool clean : {false, true}) {
    pg_graded site = clean ? pg_make_site(0.0, false)
                           : pg_make_site(pg_faint, true);
    ft& a = site.a;
    const mptensor::Shape shape = a.t.shape();
    const auto odd = [&](const mptensor::Index& idx) {
      return pg_is_odd(a.parity, idx);
    };
    const auto all = [](const mptensor::Index&) { return true; };
    INFO(std::string(clean ? "parity even" : "below the threshold")
         << "; ranks = " << a.t.get_comm_size()
         << ", rank = " << a.t.get_comm_rank());
    if (!clean) {
      REQUIRE(pg_ranks_where(pg_local_max(a.t, all) > 1.0, a.t) == 1);
      // Round-off on every rank that stores odd elements, at least two
      // ranks when there are two, so that "each rank clips its own" is
      // observed on more than one rank.
      REQUIRE(pg_ranks_where(pg_local_max(a.t, odd) > 0.0, a.t) >=
              std::min(a.t.get_comm_size(), 2));
    }

    const pg_tally tally =
        pg_run(a.t, [&] { tenes::fermion::enforce_even_parity(a); });
    INFO("message: " << tally.what);
    CHECK(tally.threw_other == 0);
    CHECK(tally.threw == 0);

    // The odd sector is zero on its owners; the even sector is what the
    // closed form (and pg_scale) put there, bit for bit.
    std::size_t odd_left = 0;
    std::size_t even_changed = 0;
    for (std::size_t n = 0; n < a.t.local_size(); ++n) {
      const mptensor::Index idx = a.t.global_index(n);
      if (odd(idx)) {
        odd_left += (a.t[n] != 0.0) ? 1 : 0;
      } else {
        const double want =
            (idx == site.scale_at) ? pg_scale : pg_entry(idx, shape);
        even_changed += (a.t[n] != want) ? 1 : 0;
      }
    }
    INFO("odd elements left on this rank = " << odd_left
                                              << ", even elements changed = "
                                              << even_changed);
    CHECK(pg_ranks_where(odd_left != 0, a.t) == 0);
    CHECK(pg_ranks_where(even_changed != 0, a.t) == 0);
  }
}

// ---- validate_block_diagonal (debug builds only) -------------------------

namespace {

//! Even-first sorted 36 x 36 matrix with 20 even rows and 17 even columns:
//! the block edges fall inside the 16-wide distribution blocks, so at two
//! ranks each rank stores part of both diagonal and off-diagonal blocks.
constexpr std::size_t pg_rows = pg_dim * pg_dim;
constexpr std::size_t pg_row_even = 20;
constexpr std::size_t pg_col_even = 17;

bool pg_diagonal_block(const mptensor::Index& idx) {
  return (idx[0] < pg_row_even) == (idx[1] < pg_col_even);
}

struct pg_blocked {
  tenes::real_tensor sorted;
  mptensor::Index scale_at;  //!< In a diagonal block; holds pg_scale.
  mptensor::Index off_at;    //!< In an off-diagonal block; the violation.
};

//! Off-diagonal blocks zero except `violation` on one element.
pg_blocked pg_make_blocked(double violation) {
  pg_blocked b{tenes::real_tensor(mptensor::Shape(pg_rows, pg_rows)), {}, {}};
  pg_fill(
      b.sorted, pg_diagonal_block,
      [](const mptensor::Index&) { return 0.0; }, violation, b.scale_at,
      b.off_at);
  return b;
}

}  // namespace

TEST_CASE(
    "validate_block_diagonal decides on every rank alike (debug builds)") {
  const char* const context = "mpi_reduction block check";
  const auto call = [&](const pg_blocked& b) {
    return pg_run(b.sorted, [&] {
      tenes::fermion::detail::validate_block_diagonal(b.sorted, pg_row_even,
                                                      pg_col_even, context);
    });
  };
  const auto off = [](const mptensor::Index& idx) {
    return !pg_diagonal_block(idx);
  };
  const auto all = [](const mptensor::Index&) { return true; };

  const pg_blocked gross = pg_make_blocked(pg_gross);
  REQUIRE(pg_ranks_where(pg_local_max(gross.sorted, off) > pg_threshold,
                         gross.sorted) == 1);
  const pg_blocked faint = pg_make_blocked(pg_faint);
  REQUIRE(pg_ranks_where(pg_local_max(faint.sorted, all) > 1.0,
                         faint.sorted) == 1);
  const pg_blocked clean = pg_make_blocked(0.0);

  const pg_tally t_gross = call(gross);
  const pg_tally t_faint = call(faint);
  const pg_tally t_clean = call(clean);
  INFO("ranks = " << t_gross.ranks
                  << ", rank = " << gross.sorted.get_comm_rank());
#ifndef NDEBUG
  {
    INFO("above the threshold; message: " << t_gross.what);
    CHECK(t_gross.threw_other == 0);
    CHECK(t_gross.threw == t_gross.ranks);
    if (!t_gross.what.empty()) {
      CHECK(pg_starts_with(t_gross.what, context));
    }
  }
  {
    INFO("below the threshold; message: " << t_faint.what);
    CHECK(t_faint.threw_other == 0);
    CHECK(t_faint.threw == 0);
  }
  {
    INFO("block diagonal; message: " << t_clean.what);
    CHECK(t_clean.threw_other == 0);
    CHECK(t_clean.threw == 0);
  }
#else
  // Release builds skip the scan: validate_block_diagonal() is a no-op
  // there, so no input throws and there is nothing to decide collectively.
  CHECK(t_gross.threw + t_gross.threw_other == 0);
  CHECK(t_faint.threw + t_faint.threw_other == 0);
  CHECK(t_clean.threw + t_clean.threw_other == 0);
#endif
}

// ---- the layer checks of doubled_pipeline_traced (debug builds only) -----
//
// Their decision is require_even_parity()'s, tested above with a violation
// that one rank stores.  doubled_pipeline_traced() itself is only ever
// handed a layer that is contaminated on EVERY rank: before the fix a rank
// that did not throw would go on into the contraction and wait there for
// the others, so a one-rank violation would hang this case instead of
// failing it.  What is checked here is that the layer checks go through
// require_even_parity(): the message on every rank starts with the layer's
// context and carries the GLOBAL odd-sector maximum pg_gross, which one
// rank stores while every other rank holds only pg_gross / 2, and the
// threshold from the GLOBAL scale.

TEST_CASE(
    "doubled_pipeline_traced: the layer checks report the global violation "
    "on every rank (debug builds)") {
  const mptensor::Shape shape(pg_dim, pg_dim, pg_dim, pg_dim, 2);
  const tenes::fermion::leg_parities legs{pg_leg(), pg_leg(), pg_leg(),
                                          pg_leg(), {false, true}};
  const pg_graded clean = pg_make_graded(
      shape, legs, [](const mptensor::Index&) { return 0.0; }, 0.0);
  const pg_graded dirty = pg_make_graded(
      shape, legs, [](const mptensor::Index&) { return 0.5 * pg_gross; },
      pg_gross);
  const auto odd = [&](const mptensor::Index& idx) {
    return pg_is_odd(legs, idx);
  };
  const tenes::real_tensor& t = dirty.a.t;
  // Safe to call on every build, fixed or not: every rank sees a violation.
  REQUIRE(pg_ranks_where(pg_local_max(t, odd) > pg_threshold, t) ==
          t.get_comm_size());
  REQUIRE(pg_ranks_where(pg_local_max(t, odd) == pg_gross, t) == 1);
  REQUIRE(pg_ranks_where(pg_local_max(clean.a.t, odd) > 0.0, t) == 0);

#ifndef NDEBUG
  const struct {
    const char* context;
    const ft& bra;
    const ft& ket;
  } layers[] = {{"doubled_pipeline_traced: bra layer", dirty.a, clean.a},
                {"doubled_pipeline_traced: ket layer", clean.a, dirty.a}};
  for (const auto& layer : layers) {
    const pg_tally tally = pg_run(t, [&] {
      static_cast<void>(tenes::fermion::detail::doubled_pipeline_traced(
          layer.bra, layer.ket));
    });
    INFO(std::string(layer.context)
         << "; ranks = " << tally.ranks << ", rank = " << t.get_comm_rank()
         << ", message: " << tally.what);
    CHECK(tally.threw_other == 0);
    CHECK(tally.threw == tally.ranks);
    if (!tally.what.empty()) {
      CHECK(pg_starts_with(tally.what, layer.context));
      CHECK(pg_mentions(tally.what, pg_gross));
      CHECK(pg_mentions(tally.what, pg_threshold));
    }
  }
#else
  // Release builds do not check the layers; nothing to decide.
#endif
}
