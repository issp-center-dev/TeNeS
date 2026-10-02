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

// ===== the physical-ledger check of the real-time simple update =============
//
// Contract: section 8 of
// docs/superpowers/specs/2026-10-01-fermion-real-time-evolution-contract.md.
// A fermion-mode real-time simple update (mode = "time") checks, once every
// gate of a step has been applied, that the physical parity ledger of every
// site is back at the input `parity`, and throws std::logic_error with the
// message of the ground-state simple update otherwise:
//
//   "fermion simple update invariant violated: physical ledger did not
//    return to its original value after a gate sweep"
//
//   A  a step whose gate chain changes a physical ledger on the way and
//      restores it (the two-hop chain of a next-nearest-neighbour hopping)
//      runs through iTPS::time_evolution() without an exception;
//   B  a state whose ledger does not match the input at the end of a step
//      makes time_evolution() throw that logic_error. Two ways to make it:
//      B1 the solver's copy of the input parity (peps_parameters.phys_parity)
//         is changed after construction, so the gates run normally and only
//         the comparison sees a difference;
//      B2 the ledger of a site no gate touches is changed after
//         construction (the labels of its two states exchanged);
//   C  the check sits at the end of the step, not after every gate: applying
//      the first gate of the chain alone leaves the ledger of the middle site
//      different from the input (so a per-gate check would refuse A), and the
//      second gate restores it.
//
// The gate chain is tenes_std's output (tool/tenes_std.py, mode = "time",
// tau = 0.01) for H = -(c^dag_0 c_3 + c^dag_3 c_0) on the (1, 1) bond of
// site 0 in a 2x2 cell of spinless fermions (parity [0, 1]). The std.toml it
// was generated from is test/data/fermion_realtime_ledger_std.toml; to
// regenerate, run
//
//     python3 tool/tenes_std.py test/data/fermion_realtime_ledger_std.toml \
//         -o input.toml
//
// and copy the two [[evolution.simple]] tables of input.toml (the
// [[evolution.full]] ones are the same gates) into rtl_input() below. The
// path is 0 -(right)-> 1 -(up)-> 3, and the first gate's out2 leg
// (dimension 8) carries the channel on site 1.
// The ledgers of the gates are inferred by infer_fermion_gate_ledgers(), as
// tenes does when it reads input.toml.
//
// Before the implementation (da46fa90) validate_fermion_constraints()
// refused mode = "time", so every case fails when it loads the input. With
// that guard lifted but no check in time_evolution() (a copy of
// src/iTPS/time_evolution.cpp without the call), B1 and B2 run through and
// fail while A and C pass; with the check moved after every gate, A fails.
//
// Its own executable: it drives the solver (time_evolution() measures with
// the mean-field environment and writes TE_*.dat into a scratch directory)
// and must not take the other fermion unit tests down if it aborts.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../test_fermion_common.hpp"

#include <filesystem>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../../src/exception.hpp"
#include "../../src/iTPS/load_toml.hpp"
#include "../../src/printlevel.hpp"

namespace {

using rtl_ptensor = tenes::complex_tensor;

const char* rtl_message =
    "fermion simple update invariant violated: physical ledger did not "
    "return to its original value after a gate sweep";

std::string rtl_input(std::string const& outdir) {
  return std::string(R"(
[parameter]
[parameter.general]
mode = "time"
fermion = true
is_real = false
output = ")") +
         outdir + R"("
measure_interval = 1
[parameter.simple_update]
tau = 0.01
num_step = 2
[parameter.full_update]
num_step = 0
[parameter.ctm]
meanfield_env = true
[parameter.random]
seed = 11

[tensor]
L_sub = [2, 2]
skew = 0
[[tensor.unitcell]]
index = []
physical_dim = 2
virtual_dim = 2
parity = [0, 1]
initial_state = [1.0, 0.0]
noise = 0.01

[observable]
[[observable.onesite]]
name = "n"
group = 0
sites = []
dim = 2
elements = """
1 1 1.0 0.0
"""

[evolution]
[[evolution.simple]]
group = 0
source_site = 0
source_leg = 2
dimensions = [2, 2, 2, 8]
elements = """
0 0 0 0  -1.4141782073286622 0.0
1 0 1 0  -1.4141782073286622 0.0
0 1 0 1  -1.4141782073286622 0.0
1 1 1 1  -1.4141782073286622 0.0
0 0 0 2  -3.535504443269654e-05 0.0
1 0 1 2  3.535504443269653e-05 0.0
0 1 0 3  -3.535504443269654e-05 0.0
1 1 1 3  3.535504443269653e-05 0.0
1 0 0 4  0.0 -0.009999833334166663
1 1 0 5  0.0 0.009999833334166663
0 0 1 6  0.0 0.009999833334166663
0 1 1 7  -0.0 -0.009999833334166663
"""
[[evolution.simple]]
group = 0
source_site = 1
source_leg = 1
dimensions = [8, 2, 2, 2]
elements = """
0 0 0 0  -0.7071067811865475 0.0
2 0 0 0  -0.7071067811865477 0.0
6 1 0 0  1.0 0.0
1 0 1 0  -0.7071067811865475 0.0
3 0 1 0  -0.7071067811865477 0.0
7 1 1 0  1.0 0.0
4 0 0 1  -1.0 0.0
0 1 0 1  -0.7071067811865477 0.0
2 1 0 1  0.7071067811865475 0.0
5 0 1 1  -1.0 0.0
1 1 1 1  -0.7071067811865477 0.0
3 1 1 1  0.7071067811865475 0.0
"""
)";
}

//! Everything tenes builds from input.toml before it constructs the solver,
//! in the order of src/iTPS/main.cpp (validate, then infer the ledgers).
struct rtl_input_data {
  tenes::itps::PEPS_Parameters params;
  tenes::SquareLattice lattice;
  tenes::EvolutionOperators<rtl_ptensor> simple_updates;
  tenes::EvolutionOperators<rtl_ptensor> full_updates;
  tenes::Operators<rtl_ptensor> onesite;
};

rtl_input_data rtl_load(std::string const& outdir) {
  using namespace tenes::itps;
  // mptensor initializes MPI on the first tensor it makes; make one before
  // anything that calls MPI directly.
  const tenes::real_tensor mpi_probe(MPI_COMM_WORLD, mptensor::Shape(1));
  static_cast<void>(mpi_probe);
  const auto toml_input = toml::parse_str(rtl_input(outdir));
  rtl_input_data in{gen_param(toml_input.at("parameter")),
                    gen_lattice(toml_input.at("tensor")),
                    {},
                    {},
                    {}};
  in.params.phys_parity = gen_phys_parity(toml_input.at("tensor"), in.lattice);
  in.params.print_level = tenes::PrintLevel::none;
  in.simple_updates =
      load_simple_updates<rtl_ptensor>(toml_input, MPI_COMM_WORLD);
  in.full_updates = load_full_updates<rtl_ptensor>(toml_input, MPI_COMM_WORLD);
  in.onesite =
      load_operators<rtl_ptensor>(toml_input, MPI_COMM_WORLD, in.lattice.N_UNIT,
                                  1, 0.0, "observable.onesite");
  validate_fermion_constraints(
      in.params, in.lattice, in.simple_updates, in.full_updates, in.onesite,
      tenes::Operators<rtl_ptensor>{}, tenes::Operators<rtl_ptensor>{},
      CorrelationParameter{});
  infer_fermion_gate_ledgers(in.params, in.lattice, in.simple_updates,
                             in.full_updates);
  return in;
}

tenes::itps::iTPS<rtl_ptensor> rtl_solver(rtl_input_data const& in) {
  return tenes::itps::iTPS<rtl_ptensor>(
      MPI_COMM_WORLD, in.params, in.lattice, in.simple_updates, in.full_updates,
      in.onesite, tenes::Operators<rtl_ptensor>{},
      tenes::Operators<rtl_ptensor>{}, tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
}

//! A fresh output directory for one case, removed by the destructor.
struct rtl_outdir {
  std::string path;
  explicit rtl_outdir(std::string p) : path(std::move(p)) {
    std::error_code ec;
    std::filesystem::remove_all(path, ec);
    std::filesystem::create_directories(path);
  }
  ~rtl_outdir() {
    std::error_code ec;
    std::filesystem::remove_all(path, ec);
  }
};

//! Premises shared by every case: the input really is a two-gate chain in
//! one group whose first gate leaves the middle site (1) on a ledger other
//! than its input parity, and whose second gate brings it back.
void rtl_require_chain(rtl_input_data const& in) {
  REQUIRE(in.params.fermion == true);
  REQUIRE(in.params.calcmode ==
          tenes::itps::PEPS_Parameters::CalculationMode::time_evolution);
  REQUIRE(in.params.num_simple_step.size() == 1);
  REQUIRE(in.params.num_simple_step[0] == 2);
  REQUIRE(in.simple_updates.size() == 2);
  for (auto const& up : in.simple_updates) {
    REQUIRE(up.group == 0);
    REQUIRE(up.fermion_legs.size() == 4);
  }
  REQUIRE(in.simple_updates[0].source_site == 0);
  REQUIRE(in.simple_updates[1].source_site == 1);
  // gate 0 writes site 1 through its out2 leg, gate 1 reads it through in1
  REQUIRE(in.simple_updates[0].fermion_legs[3] != in.params.phys_parity[1]);
  REQUIRE(in.simple_updates[1].fermion_legs[0] ==
          in.simple_updates[0].fermion_legs[3]);
  REQUIRE(in.simple_updates[1].fermion_legs[2] == in.params.phys_parity[1]);
}

//! Run time_evolution() and return the message of the std::logic_error it
//! throws ("" if it throws nothing). Any other exception fails the case.
std::string rtl_run(tenes::itps::iTPS<rtl_ptensor>& solver) {
  try {
    solver.time_evolution();
  } catch (std::logic_error const& e) {
    return std::string(e.what());
  } catch (std::exception const& e) {
    FAIL_CHECK("time_evolution() threw something other than a logic_error: "
               << std::string(e.what()));
    return std::string("(other exception)");
  }
  return std::string();
}

}  // namespace

TEST_CASE(
    "fermion real-time ledger A: a chain that changes the ledger within a "
    "step runs through") {
  rtl_outdir out("output_test_fermion_realtime_ledger_A");
  auto in = rtl_load(out.path);
  rtl_require_chain(in);
  auto solver = rtl_solver(in);
  const std::string thrown = rtl_run(solver);
  INFO("time_evolution() threw: " << thrown);
  CHECK(thrown.empty());
  auto const& finfo = tenes::itps::iTPSTestAccessor::finfo(solver);
  for (int site = 0; site < in.lattice.N_UNIT; ++site) {
    CHECK(finfo.phys[site] == in.params.phys_parity[site]);
  }
  // ... and it really took both steps: the time column of TE_onesite_obs.dat
  // has t = 0, 0.01 and 0.02 (measure_interval = 1).
  std::ifstream ifs(std::filesystem::path(out.path) / "TE_onesite_obs.dat");
  REQUIRE(ifs.good());
  std::set<std::string> times;
  std::string line;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream words(line);
    std::string time;
    words >> time;
    times.insert(time);
  }
  CHECK(times.size() == 3);
}

TEST_CASE(
    "fermion real-time ledger B1: an input parity that the ledger does not "
    "return to is caught at the end of the step") {
  rtl_outdir out("output_test_fermion_realtime_ledger_B1");
  auto in = rtl_load(out.path);
  rtl_require_chain(in);
  auto solver = rtl_solver(in);
  // The solver keeps its own copy of the parameters; the shared accessor
  // hands it out const, and the solver object itself is not const.
  auto& params = const_cast<tenes::itps::PEPS_Parameters&>(
      tenes::itps::iTPSTestAccessor::peps_parameters(solver));
  REQUIRE(params.phys_parity[2] == std::vector<bool>{false, true});
  params.phys_parity[2] = std::vector<bool>{true, false};
  const std::string thrown = rtl_run(solver);
  INFO("time_evolution() threw: " << thrown);
  CHECK(thrown == rtl_message);
}

TEST_CASE(
    "fermion real-time ledger B2: a ledger that differs from the input "
    "parity is caught at the end of the step") {
  rtl_outdir out("output_test_fermion_realtime_ledger_B2");
  auto in = rtl_load(out.path);
  rtl_require_chain(in);
  auto solver = rtl_solver(in);
  auto& finfo = tenes::itps::iTPSTestAccessor::finfo(solver);
  // Site 2 is on none of the two gates (they act on 0-1 and 1-3), so the
  // run itself never rewrites its ledger.
  REQUIRE(finfo.phys[2] == in.params.phys_parity[2]);
  finfo.phys[2] = tenes::fermion::parity_vector{true, false};
  const std::string thrown = rtl_run(solver);
  INFO("time_evolution() threw: " << thrown);
  CHECK(thrown == rtl_message);
}

TEST_CASE(
    "fermion real-time ledger C: the ledger is checked after the step, not "
    "after every gate") {
  rtl_outdir out("output_test_fermion_realtime_ledger_C");
  auto in = rtl_load(out.path);
  rtl_require_chain(in);
  auto solver = rtl_solver(in);
  auto const& finfo = tenes::itps::iTPSTestAccessor::finfo(solver);
  REQUIRE(finfo.phys[1] == in.params.phys_parity[1]);
  // After the first gate alone the middle site is on the gate's out2 ledger,
  // so a check after every gate would have refused case A ...
  solver.simple_update(in.simple_updates[0]);
  CHECK(finfo.phys[1] != in.params.phys_parity[1]);
  CHECK(finfo.phys[1] == in.simple_updates[0].fermion_legs[3]);
  // ... and the second gate brings it back, so the end-of-step check passes.
  solver.simple_update(in.simple_updates[1]);
  CHECK(finfo.phys[1] == in.params.phys_parity[1]);
}
