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

// ===== The graded decompositions stay on the tensor's own communicator ====
//
// tenes_itps_main() takes the communicator to run on. From the command line
// that is MPI_COMM_WORLD, but a library caller may hand it a subset of the
// processes -- two solvers side by side on halves of the machine, say. Every
// collective the solver makes must then be on the communicator its tensors
// live on, and never on MPI_COMM_WORLD.
//
// These cases pin that for the fermionic decompositions and for the
// auxiliary identities the full-update environment builds.
//
// How they detect the defect: MPI_COMM_WORLD is split into two groups, and
// the groups are given DIFFERENT work -- different matrix shapes and a
// different number of calls. A routine that is internally collective over
// all processes then either blocks (its partner is inside a different call,
// or none) or mixes the two groups' data, while a routine that only talks
// to its own group is unaffected by what the other group is doing. Giving
// both groups the SAME matrices would hide exactly the defect under test,
// because then the world-wide collectives would line up by accident.
//
// The oracles are self-contained: Q R reproduces the input, U diag(s) V^T
// reproduces the input, and the factors' communicator is the input's. The
// contract asks nothing about the values of the factors beyond that -- the
// serial correctness of the graded decompositions is pinned by
// fermion/local_decomposition.cpp and by the fold/full-update cases.
//
// Every check below is a CHECK rather than a REQUIRE inside the group-local
// part: the two groups deliberately run different code paths, so a REQUIRE
// that stopped one group's case early would leave the other group waiting in
// the closing barrier. Premises shared by all ranks are REQUIREd.
//
// A build with assertions on reports the defect earlier and more plainly
// than the group comparison does: mptensor's tensordot asserts
// a.get_comm() == b.get_comm(), and a communicator handle made by
// MPI_Comm_split is not MPI_COMM_WORLD even at one process. That is why
// there is a one-rank registration as well. Under NDEBUG the assertion is
// gone, and the two-rank registration -- different shapes, different call
// counts per group -- is the only thing left that can see it.

#define DOCTEST_CONFIG_IMPLEMENT
#include "../doctest.h"
#include "../test_workdir.hpp"

#include <cmath>
#include <string>
#include <vector>

#include "../../src/fermion/fops.hpp"
#include "../../src/fermion/ftensor.hpp"
#include "../../src/fermion/full_update_env.hpp"
#include "../../src/fermion/parity.hpp"
#include "../../src/fermion/reduced.hpp"
#include "../../src/initialize_mptensor.hpp"
#include "../../src/iTPS/iTPS.hpp"
#include "../../src/iTPS/load_toml.hpp"
#include "../../src/mpi.hpp"
#include "../../src/tensor.hpp"
#include "../../src/util/file.hpp"

namespace tenes::itps {
//! Reaches the site tensors, which iTPS keeps private. iTPS declares this
//! name a friend (iTPS.hpp); test_fermion_common.hpp declares a fuller one
//! for the executables that need more, and this file is not part of any of
//! them, so the two never meet in one program.
struct iTPSTestAccessor {
  template <class tensor>
  static std::vector<tensor> &Tn(iTPS<tensor> &state) {
    return state.Tn;
  }
  template <class tensor>
  static std::vector<tensor> &eTt(iTPS<tensor> &state) {
    return state.eTt;
  }
};
}  // namespace tenes::itps

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  tenes::initialize_mptensor();
  doctest::Context context(argc, argv);
  const int res = context.run();
  MPI_Finalize();
  return res;
}

namespace {

using ftr = tenes::fermion::ftensor<tenes::real_tensor>;

int sd_world_size() {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
}

int sd_world_rank() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

//! Whether two communicators hold the same processes in the same order.
//! Note that at one process every communicator is congruent to every other,
//! so this check only has power from two processes upwards -- which is
//! where the registration at two ranks comes in.
bool sd_same_comm(MPI_Comm a, MPI_Comm b) {
#ifdef _NO_MPI
  static_cast<void>(a);
  static_cast<void>(b);
  return true;
#else
  int result = MPI_UNEQUAL;
  MPI_Comm_compare(a, b, &result);
  return result == MPI_IDENT || result == MPI_CONGRUENT;
#endif
}

//! Parity ledger from a string of 'e' (even) and 'o' (odd).
tenes::fermion::parity_vector sd_pv(const std::string &s) {
  tenes::fermion::parity_vector p;
  for (std::size_t i = 0; i < s.size(); ++i) {
    p.push_back(s[i] == 'o');
  }
  return p;
}

//! Deterministic entry, fixed by the global index and a seed alone, so that
//! every process of a group fills the same matrix however it is distributed.
double sd_draw(int seed, const mptensor::Index &idx, std::size_t rank) {
  const double w[4] = {1.7, 3.1, 5.3, 7.9};
  double x = 0.37 * seed + 0.011 * seed * seed;
  for (std::size_t ax = 0; ax < rank; ++ax) {
    x += w[ax % 4] * static_cast<double>(idx[ax]);
  }
  return 0.73 * std::sin(x) + 0.31 * std::cos(1.7 * x + 0.19 * seed);
}

//! A parity-even graded tensor of the shape the ledgers imply, on @p comm.
ftr sd_even_ft(MPI_Comm comm, const tenes::fermion::leg_parities &legs,
               int seed) {
  mptensor::Shape shape;
  for (std::size_t ax = 0; ax < legs.size(); ++ax) {
    shape.push(legs[ax].size());
  }
  tenes::real_tensor t(comm, shape);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    const bool odd = tenes::fermion::count_odd(legs, idx) % 2 == 1;
    t.set_value(idx, odd ? 0.0 : sd_draw(seed, idx, legs.size()));
  }
  return ftr{t, legs};
}

//! Tensor of all ones on @p comm.
tenes::real_tensor sd_ones(MPI_Comm comm, const mptensor::Shape &sh) {
  tenes::real_tensor t(comm, sh);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    t.set_value(t.global_index(n), 1.0);
  }
  return t;
}

//! One job of a group: a rank-2 graded matrix and a seed.
struct sd_job {
  std::string rows;  //!< Ledger of the row leg.
  std::string cols;  //!< Ledger of the column leg.
  int seed;
};

/*! The jobs of group @p color.
 *
 *  The two lists differ in length AND in the shapes they name. That is the
 *  whole detection mechanism: a decomposition that is collective over all
 *  processes cannot survive one group asking for a 4x4 while the other asks
 *  for a 6x6, nor one group making three calls while the other makes one.
 */
std::vector<sd_job> sd_jobs(int color) {
  if (color == 0) {
    return {{"eoeo", "eeoo", 101}};
  }
  return {
      {"eooeoo", "eoeoeo", 211}, {"eo", "eo", 212}, {"eoeoeo", "eeoo", 213}};
}

//! Number of processes of @p comm on which @p flag holds.
int sd_count(bool flag, MPI_Comm comm) {
  int count = flag ? 1 : 0;
  tenes::allreduce_sum(count, comm);
  return count;
}

}  // namespace

// ---- controls: the jobs themselves are well posed ------------------------
//
// Run on MPI_COMM_WORLD, where the communicator question does not arise,
// these say that every matrix below decomposes and reassembles. Without
// them a failure of the cases that follow could be read as "the fixture is
// malformed" rather than "the decomposition left its communicator", so they
// come first: if the contract cases abort the process, this one has already
// reported.

TEST_CASE("control: every job of both groups decomposes on MPI_COMM_WORLD") {
  for (int color = 0; color < 2; ++color) {
    const std::vector<sd_job> jobs = sd_jobs(color);
    for (std::size_t j = 0; j < jobs.size(); ++j) {
      const tenes::fermion::leg_parities legs{sd_pv(jobs[j].rows),
                                              sd_pv(jobs[j].cols)};
      const ftr a = sd_even_ft(MPI_COMM_WORLD, legs, jobs[j].seed);
      INFO("group " << color << " job " << j << ": " << jobs[j].rows << " x "
                    << jobs[j].cols);

      ftr q, r;
      REQUIRE(tenes::fermion::qr(a, mptensor::Axes(0), mptensor::Axes(1), q,
                                 r) == 0);
      const ftr qr_recon =
          tenes::fermion::tensordot(q, r, mptensor::Axes(1), mptensor::Axes(0));
      CHECK(mptensor::max_abs(qr_recon.t - a.t) < 1.0e-12);

      ftr u, vt;
      std::vector<double> s;
      REQUIRE(tenes::fermion::svd(a, mptensor::Axes(0), mptensor::Axes(1), u, s,
                                  vt) == 0);
      ftr us = u;
      us.multiply_vector(s, 1);
      const ftr svd_recon = tenes::fermion::tensordot(us, vt, mptensor::Axes(1),
                                                      mptensor::Axes(0));
      CHECK(mptensor::max_abs(svd_recon.t - a.t) < 1.0e-12);

      ftr tu, tvt;
      std::vector<double> ts;
      const int dc = static_cast<int>(jobs[j].rows.size()) - 1;
      REQUIRE(tenes::fermion::svd_trunc(a, mptensor::Axes(0), mptensor::Axes(1),
                                        tu, ts, tvt, dc) == 0);
      CHECK(static_cast<int>(ts.size()) <= dc);
    }
  }
}

// ---- the graded QR and SVD ----------------------------------------------

TEST_CASE(
    "the graded QR and SVD of a tensor on a sub-communicator stay inside it") {
  const int world = sd_world_size();
  const int rank = sd_world_rank();
  INFO("world size = " << world << ", world rank = " << rank);
  if (world < 2) {
    MESSAGE(
        "one process only: every communicator is congruent to every other, so "
        "this case cannot tell a world-wide collective from a group-local "
        "one. It is registered a second time at two ranks, which is where it "
        "has power.");
  }

  const int color = rank % 2;
  MPI_Comm sub = MPI_COMM_WORLD;
#ifndef _NO_MPI
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &sub);
#endif

  const std::vector<sd_job> jobs = sd_jobs(color);
  int failures = 0;
  for (std::size_t j = 0; j < jobs.size(); ++j) {
    const tenes::fermion::leg_parities legs{sd_pv(jobs[j].rows),
                                            sd_pv(jobs[j].cols)};
    const ftr a = sd_even_ft(sub, legs, jobs[j].seed);
    INFO("group " << color << " job " << j << ": " << jobs[j].rows << " x "
                  << jobs[j].cols);

    {
      ftr q, r;
      const int info =
          tenes::fermion::qr(a, mptensor::Axes(0), mptensor::Axes(1), q, r);
      CHECK(info == 0);
      failures += (info != 0) ? 1 : 0;
      const bool q_here = sd_same_comm(q.t.get_comm(), sub);
      const bool r_here = sd_same_comm(r.t.get_comm(), sub);
      CHECK(q_here);
      CHECK(r_here);
      failures += (q_here ? 0 : 1) + (r_here ? 0 : 1);
      const ftr recon =
          tenes::fermion::tensordot(q, r, mptensor::Axes(1), mptensor::Axes(0));
      const double diff = mptensor::max_abs(recon.t - a.t);
      INFO("max |QR - A| = " << diff);
      CHECK(diff < 1.0e-12);
      failures += (diff < 1.0e-12) ? 0 : 1;
    }

    {
      ftr u, vt;
      std::vector<double> s;
      const int info = tenes::fermion::svd(a, mptensor::Axes(0),
                                           mptensor::Axes(1), u, s, vt);
      CHECK(info == 0);
      failures += (info != 0) ? 1 : 0;
      const bool u_here = sd_same_comm(u.t.get_comm(), sub);
      const bool vt_here = sd_same_comm(vt.t.get_comm(), sub);
      CHECK(u_here);
      CHECK(vt_here);
      failures += (u_here ? 0 : 1) + (vt_here ? 0 : 1);
      ftr us = u;
      us.multiply_vector(s, 1);
      const ftr recon = tenes::fermion::tensordot(us, vt, mptensor::Axes(1),
                                                  mptensor::Axes(0));
      const double diff = mptensor::max_abs(recon.t - a.t);
      INFO("max |U S Vt - A| = " << diff);
      CHECK(diff < 1.0e-12);
      failures += (diff < 1.0e-12) ? 0 : 1;
    }
  }

  // Counted over all processes, so the summary fails on every rank at once
  // rather than on whichever rank happened to hold the bad element.
  const int total = sd_count(failures != 0, MPI_COMM_WORLD);
  INFO("processes with at least one failed check = " << total);
  CHECK(total == 0);

#ifndef _NO_MPI
  MPI_Comm_free(&sub);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

TEST_CASE(
    "the truncated graded SVD of a tensor on a sub-communicator stays inside "
    "it") {
  const int world = sd_world_size();
  const int rank = sd_world_rank();
  INFO("world size = " << world << ", world rank = " << rank);

  const int color = rank % 2;
  MPI_Comm sub = MPI_COMM_WORLD;
#ifndef _NO_MPI
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &sub);
#endif

  // Different shapes and a different number of kept values per group, again
  // so that a world-wide collective cannot line up by accident.
  const std::vector<sd_job> jobs = sd_jobs(color);
  int failures = 0;
  for (std::size_t j = 0; j < jobs.size(); ++j) {
    const tenes::fermion::leg_parities legs{sd_pv(jobs[j].rows),
                                            sd_pv(jobs[j].cols)};
    const ftr a = sd_even_ft(sub, legs, jobs[j].seed + 1000);
    const int dc = static_cast<int>(jobs[j].rows.size()) - 1;
    INFO("group " << color << " job " << j << ": " << jobs[j].rows << " x "
                  << jobs[j].cols << ", dc = " << dc);

    ftr u, vt;
    std::vector<double> s;
    const int info = tenes::fermion::svd_trunc(a, mptensor::Axes(0),
                                               mptensor::Axes(1), u, s, vt, dc);
    CHECK(info == 0);
    failures += (info != 0) ? 1 : 0;
    const bool right_size = static_cast<int>(s.size()) <= dc;
    CHECK(right_size);
    failures += right_size ? 0 : 1;
    const bool u_here = sd_same_comm(u.t.get_comm(), sub);
    const bool vt_here = sd_same_comm(vt.t.get_comm(), sub);
    CHECK(u_here);
    CHECK(vt_here);
    failures += (u_here ? 0 : 1) + (vt_here ? 0 : 1);
    // The kept singular values are the leading ones of the full spectrum,
    // so the truncated product is within the discarded weight of A. The
    // bound is deliberately loose: what is under test is the communicator,
    // not the truncation, which local_decomposition.cpp already pins.
    ftr us = u;
    us.multiply_vector(s, 1);
    const ftr recon =
        tenes::fermion::tensordot(us, vt, mptensor::Axes(1), mptensor::Axes(0));
    const double diff = mptensor::max_abs(recon.t - a.t);
    const double scale = mptensor::max_abs(a.t);
    INFO("max |U S Vt - A| = " << diff << ", max |A| = " << scale);
    CHECK(diff <= scale);
    failures += (diff <= scale) ? 0 : 1;
  }

  const int total = sd_count(failures != 0, MPI_COMM_WORLD);
  INFO("processes with at least one failed check = " << total);
  CHECK(total == 0);

#ifndef _NO_MPI
  MPI_Comm_free(&sub);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

// ---- the auxiliary identities of the full-update environment -------------
//
// build_full_update_environment() builds three helper tensors that carry no
// data of their own: the two open-channel identity factors and the identity
// it closes the window with to take its norm. They are what the last
// sentence of the contract's section on sub-communicators refers to.
//
// The window here is the degenerate one: both sites sit alone, so all ten
// CTM slots are the all-ones tensor of dimension 1 and the only dimensions
// that matter are the two open channels, whose sizes differ per group. The
// environment is therefore not a physical one -- forbidden_tol is set so
// high that the parity guard cannot fire -- and no claim is made here about
// the VALUE of N. What is claimed is that the call completes on the group's
// own communicator and returns a tensor that lives there, with the shape
// the two channels imply.

namespace {

//! The ten trivial CTM slots of a window whose every environment leg has
//! dimension 1.
std::vector<tenes::real_tensor> sd_trivial_env(MPI_Comm comm) {
  std::vector<tenes::real_tensor> e;
  for (int i = 0; i < 4; ++i) {
    e.push_back(sd_ones(comm, mptensor::Shape(1, 1)));
  }
  for (int i = 0; i < 6; ++i) {
    e.push_back(sd_ones(comm, mptensor::Shape(1, 1, 1)));
  }
  return e;
}

//! An environment-side QR factor with legs (1, 1, 1, channel). Parity even,
//! because the layer checks of the fold refuse anything else in a debug
//! build. With three even legs of dimension 1 that leaves the even half of
//! the channel carrying the weight, which is enough for the window to have a
//! norm -- the only numerical property this case needs.
ftr sd_q_factor(MPI_Comm comm, const std::string &channel, int seed) {
  const tenes::fermion::parity_vector one = sd_pv("e");
  const tenes::fermion::leg_parities legs{one, one, one, sd_pv(channel)};
  return sd_even_ft(comm, legs, seed);
}

//! The open channel of group @p color: different sizes per group, so that a
//! world-wide collective inside the identity builders cannot line up.
std::string sd_channel_a(int color) { return color == 0 ? "eo" : "eooe"; }
std::string sd_channel_b(int color) { return color == 0 ? "eeo" : "eo"; }

}  // namespace

TEST_CASE(
    "control: the window of both groups builds an environment on "
    "MPI_COMM_WORLD") {
  for (int color = 0; color < 2; ++color) {
    const std::vector<tenes::real_tensor> e = sd_trivial_env(MPI_COMM_WORLD);
    const std::string ca = sd_channel_a(color);
    const std::string cb = sd_channel_b(color);
    const ftr QA = sd_q_factor(MPI_COMM_WORLD, ca, 501 + color);
    const ftr QB = sd_q_factor(MPI_COMM_WORLD, cb, 601 + color);
    INFO("group " << color << ": nA = " << ca.size() << ", nB = " << cb.size());
    const tenes::fermion::full_update_environment<tenes::real_tensor> env =
        tenes::fermion::build_full_update_environment(
            e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9], QA, QB,
            tenes::fermion::reduced_pair_direction::horizontal, 1.0e300);
    const mptensor::Shape got = env.N.t.shape();
    REQUIRE(got.size() == 4u);
    CHECK(got[0] == ca.size());
    CHECK(got[1] == cb.size());
    CHECK(got[2] == ca.size());
    CHECK(got[3] == cb.size());
    const double scale = mptensor::max_abs(env.N.t);
    INFO("max |N| = " << scale);
    CHECK(std::isfinite(scale));
    CHECK(scale > 0.0);
  }
}

TEST_CASE(
    "the full-update environment of a window on a sub-communicator stays "
    "inside it") {
  const int world = sd_world_size();
  const int rank = sd_world_rank();
  INFO("world size = " << world << ", world rank = " << rank);

  const int color = rank % 2;
  MPI_Comm sub = MPI_COMM_WORLD;
#ifndef _NO_MPI
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &sub);
#endif

  const std::vector<tenes::real_tensor> e = sd_trivial_env(sub);
  const std::string ca = sd_channel_a(color);
  const std::string cb = sd_channel_b(color);
  const ftr QA = sd_q_factor(sub, ca, 501 + color);
  const ftr QB = sd_q_factor(sub, cb, 601 + color);
  const std::size_t nA = ca.size();
  const std::size_t nB = cb.size();
  INFO("group " << color << ": nA = " << nA << ", nB = " << nB);

  // The window is deliberately unphysical, so the parity guard is disarmed;
  // this case is about the communicator, not about the environment's value.
  const double wide_open = 1.0e300;
  const tenes::fermion::full_update_environment<tenes::real_tensor> env =
      tenes::fermion::build_full_update_environment(
          e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9], QA, QB,
          tenes::fermion::reduced_pair_direction::horizontal, wide_open);

  int failures = 0;
  const bool here = sd_same_comm(env.N.t.get_comm(), sub);
  CHECK(here);
  failures += here ? 0 : 1;
  const mptensor::Shape got = env.N.t.shape();
  const bool shaped = got.size() == 4 && got[0] == nA && got[1] == nB &&
                      got[2] == nA && got[3] == nB;
  CHECK(shaped);
  failures += shaped ? 0 : 1;
  const double scale = mptensor::max_abs(env.N.t);
  INFO("max |N| = " << scale);
  const bool finite = std::isfinite(scale) && scale > 0.0;
  CHECK(finite);
  failures += finite ? 0 : 1;

  const int total = sd_count(failures != 0, MPI_COMM_WORLD);
  INFO("processes with at least one failed check = " << total);
  CHECK(total == 0);

#ifndef _NO_MPI
  MPI_Comm_free(&sub);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

// ---- reading a checkpoint back onto the group's own communicator --------
//
// Addendum 6. The decompositions above are not the only place a
// communicator can leak. mptensor's load() leaves a tensor's communicator
// alone, so whichever one the destination tensor was built with is the one
// the loaded data ends up on -- and a destination built without an explicit
// communicator is on MPI_COMM_WORLD. A solver handed a subset of the
// processes would then spread every loaded tensor over all of them.
//
// Each group saves its own checkpoint on its own communicator and reads it
// back there, so the two groups never touch the same directory and the
// number of processes that wrote a checkpoint is the number that reads it.

namespace {

//! A two-by-two unit cell at virtual dimension 2: the smallest thing that
//! goes through the whole save and load path.
tenes::SquareLattice sd_lattice() {
  return tenes::itps::gen_lattice(toml::parse_str(R"(
[tensor]
L_sub = [2, 2]
[[tensor.unitcell]]
index = []
physical_dim = 2
virtual_dim = 2
initial_state = [1.0, 0.0]
noise = 0.01
  )")
                                      .at("tensor"));
}

tenes::itps::iTPS<tenes::real_tensor> sd_make_itps(
    MPI_Comm comm, tenes::SquareLattice const &lattice,
    tenes::itps::PEPS_Parameters const &params) {
  using tenes::real_tensor;
  return tenes::itps::iTPS<real_tensor>(
      comm, params, lattice, tenes::EvolutionOperators<real_tensor>{},
      tenes::EvolutionOperators<real_tensor>{}, tenes::Operators<real_tensor>{},
      tenes::Operators<real_tensor>{}, tenes::Operators<real_tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
}

}  // namespace

TEST_CASE(
    "a checkpoint read on a sub-communicator gives tensors that live there") {
  const int world = sd_world_size();
  const int rank = sd_world_rank();
  INFO("world size = " << world << ", world rank = " << rank);

  const int color = rank % 2;
  MPI_Comm sub = MPI_COMM_WORLD;
#ifndef _NO_MPI
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &sub);
#endif

  // One directory per group: the checkpoint a group reads is the one it
  // wrote, on the communicator it wrote it with.
  const std::string tag = "output_subcomm_load_" + std::to_string(color);
  const tenes::SquareLattice lattice = sd_lattice();

  tenes::itps::PEPS_Parameters save_params;
  save_params.print_level = tenes::PrintLevel::none;
  save_params.outdir = tag;
  save_params.tensor_save_dir = tag + "/tensors";

  tenes::itps::PEPS_Parameters load_params = save_params;
  load_params.tensor_save_dir = "";
  load_params.tensor_load_dir = save_params.tensor_save_dir;

  int failures = 0;
  {
    auto saver = sd_make_itps(sub, lattice, save_params);
    const bool saved = saver.save_tensors();
    CHECK(saved);
    failures += saved ? 0 : 1;
  }

  // The load happens in the constructor.
  auto loader = sd_make_itps(sub, lattice, load_params);
  auto &Tn = tenes::itps::iTPSTestAccessor::Tn(loader);
  const bool any = !Tn.empty();
  CHECK(any);
  failures += any ? 0 : 1;
  if (any) {
    int wrong = 0;
    for (std::size_t i = 0; i < Tn.size(); ++i) {
      if (!sd_same_comm(Tn[i].get_comm(), sub)) {
        ++wrong;
      }
    }
    INFO("site tensors on the wrong communicator = " << wrong << " of "
                                                     << Tn.size());
    CHECK(wrong == 0);
    failures += (wrong == 0) ? 0 : 1;
  }

  const int total = sd_count(failures != 0, MPI_COMM_WORLD);
  INFO("processes with at least one failed check = " << total);
  CHECK(total == 0);

#ifndef _NO_MPI
  MPI_Comm_free(&sub);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

// ---- the same for a FERMIONIC checkpoint ---------------------------------
//
// Addendum 2, item 2. The case above saves a bosonic state, and a bosonic
// save writes the CTM environment straight out of the solver's own tensors,
// which already carry the right communicator. A fermionic save does not: the
// folded environment cannot be written as it stands, so save_tensors() makes
// a zero-filled placeholder for each of the ten CTM slots. Building those
// placeholders without the communicator is invisible to every other test --
// the reviewer's mutation check found exactly that hole.
//
// At more than one process the hole is visible in the files. A placeholder on
// MPI_COMM_WORLD asks mptensor to write the base file from WORLD rank 0 and
// the fragment for whatever rank it thinks it is, so the group whose member
// is not world rank 0 ends up with no base file and a fragment under the
// wrong number; the reload then cannot find the tensor. At one process the
// two communicators hold the same single process and nothing distinguishes
// them, which is why this case, like the others here, has its power at two.

namespace {

//! The parity ledger of the physical legs: one fermionic mode per site, so
//! basis state 0 is even and state 1 is odd.
std::vector<std::vector<bool>> sd_phys_parity(int n_unit) {
  return std::vector<std::vector<bool>>(static_cast<std::size_t>(n_unit),
                                        std::vector<bool>{false, true});
}

}  // namespace

TEST_CASE(
    "a fermionic checkpoint saved and read on a sub-communicator stays "
    "there") {
  const int world = sd_world_size();
  const int rank = sd_world_rank();
  INFO("world size = " << world << ", world rank = " << rank);

  const int color = rank % 2;
  MPI_Comm sub = MPI_COMM_WORLD;
#ifndef _NO_MPI
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &sub);
#endif

  const std::string tag =
      "output_subcomm_fermion_load_" + std::to_string(color);
  const std::string tensors = tag + "/tensors";
  const tenes::SquareLattice lattice = sd_lattice();

  tenes::itps::PEPS_Parameters save_params;
  save_params.print_level = tenes::PrintLevel::none;
  save_params.outdir = tag;
  save_params.tensor_save_dir = tensors;
  save_params.fermion = true;
  save_params.phys_parity = sd_phys_parity(lattice.N_UNIT);

  tenes::itps::PEPS_Parameters load_params = save_params;
  load_params.tensor_save_dir = "";
  load_params.tensor_load_dir = tensors;

  int failures = 0;
  {
    auto saver = sd_make_itps(sub, lattice, save_params);
    const bool saved = saver.save_tensors();
    CHECK(saved);
    failures += saved ? 0 : 1;
  }

  // Premise: this really was the fermionic path. The parity ledger is the
  // file only a fermionic save writes, and the CTM placeholders are what the
  // fermionic save puts where a bosonic one writes the environment.
  const bool ledger = tenes::util::path_exists(tensors + "/fermion.dat");
  CHECK(ledger);
  failures += ledger ? 0 : 1;

  // Written by this group, so numbered from this group's own rank 0. A
  // placeholder built on MPI_COMM_WORLD would put the base file and the
  // fragment somewhere else.
  const bool placeholder = tenes::util::path_exists(tensors + "/Et_0.dat");
  const bool fragment =
      tenes::util::path_exists(tensors + "/Et_0.dat.0000.bin");
  INFO("Et_0.dat present = " << placeholder
                             << ", Et_0.dat.0000.bin present = " << fragment);
  CHECK(placeholder);
  CHECK(fragment);
  failures += (placeholder ? 0 : 1) + (fragment ? 0 : 1);

  if (placeholder && fragment) {
    auto loader = sd_make_itps(sub, lattice, load_params);
    auto &Tn = tenes::itps::iTPSTestAccessor::Tn(loader);
    auto &eTt = tenes::itps::iTPSTestAccessor::eTt(loader);
    const bool any = !Tn.empty() && !eTt.empty();
    CHECK(any);
    failures += any ? 0 : 1;
    if (any) {
      int wrong = 0;
      for (std::size_t i = 0; i < Tn.size(); ++i) {
        if (!sd_same_comm(Tn[i].get_comm(), sub)) {
          ++wrong;
        }
      }
      for (std::size_t i = 0; i < eTt.size(); ++i) {
        if (!sd_same_comm(eTt[i].get_comm(), sub)) {
          ++wrong;
        }
      }
      INFO("tensors on the wrong communicator = " << wrong << " of "
                                                  << (Tn.size() + eTt.size()));
      CHECK(wrong == 0);
      failures += (wrong == 0) ? 0 : 1;
    }
  }

  const int total = sd_count(failures != 0, MPI_COMM_WORLD);
  INFO("processes with at least one failed check = " << total);
  CHECK(total == 0);

#ifndef _NO_MPI
  MPI_Comm_free(&sub);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}
