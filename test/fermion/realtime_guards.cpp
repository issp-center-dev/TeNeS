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

// ===== the input guards of fermion real-time evolution ======================
//
// Contract: section 6 of
// docs/superpowers/specs/2026-10-01-fermion-real-time-evolution-contract.md.
//
//   1. fermion mode with mode = "time" is accepted (no exception);
//   2. fermion mode with mode = "finite" is refused with tenes::input_error,
//      and the message names "finite-temperature";
//   3. fermion mode with mode = "time" and Simple_Gauge_Fix = true is still
//      refused, and the message names "Simple_Gauge_Fix".
//
// The cases call validate_fermion_constraints() directly, the way
// test/input.cpp does for the other fermion guards, with the parameters and
// the lattice built from TOML by gen_param() / gen_lattice() /
// gen_phys_parity() (src/iTPS/load_toml.hpp). The mode strings are the ones
// gen_param() parses ("time..." and "finite..." prefixes); each case asserts
// as a premise that the parsed calcmode is the one it is about, so a failure
// below is about the guard and not about the parsing.
//
// Before the implementation case 1 fails (the guard refuses every
// non-ground-state mode, "non-ground-state mode") and case 2 fails on the
// message (the same "non-ground-state mode" text names no
// "finite-temperature").
//
// Its own executable: a unit test that needs no tensors, no MPI and no
// solver, and nothing else in the fermion suite to fail along with it.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../doctest.h"
#include "../test_workdir.hpp"

#include <string>

#include "../../src/exception.hpp"
#include "../../src/iTPS/load_toml.hpp"
#include "../../src/iTPS/PEPS_Parameters.hpp"
#include "../../src/SquareLattice.hpp"
#include "../../src/tensor.hpp"

namespace {

// A 2x2 cell (no site is its own neighbour) with an odd physical state, so
// that fermion mode is fully on, and an even initial state; none of the other
// fermion guards can fire on it.
const char* rtg_cell_toml = R"(
[tensor]
L_sub = [2, 2]
skew = 0
[[tensor.unitcell]]
index = []
physical_dim = 2
virtual_dim = 2
parity = [0, 1]
initial_state = [1.0, 0.0]
noise = 0.0
)";

struct rtg_input {
  tenes::itps::PEPS_Parameters peps_parameters;
  tenes::SquareLattice lattice;
};

rtg_input rtg_make(std::string const& general_extra,
                   std::string const& simple_extra) {
  using namespace tenes::itps;
  auto tensor_toml = toml::parse_str(rtg_cell_toml);
  tenes::SquareLattice lattice = gen_lattice(tensor_toml.at("tensor"));
  const std::string param_text = std::string(R"(
[parameter]
[parameter.general]
fermion = true
)") + general_extra + R"(
[parameter.simple_update]
tau = 0.01
num_step = 10
)" + simple_extra;
  auto param_toml = toml::parse_str(param_text);
  PEPS_Parameters peps_parameters = gen_param(param_toml.at("parameter"));
  peps_parameters.phys_parity =
      gen_phys_parity(tensor_toml.at("tensor"), lattice);
  return rtg_input{peps_parameters, lattice};
}

void rtg_validate(rtg_input const& in) {
  using ptensor = tenes::complex_tensor;
  tenes::itps::validate_fermion_constraints(
      in.peps_parameters, in.lattice, tenes::EvolutionOperators<ptensor>{},
      tenes::EvolutionOperators<ptensor>{}, tenes::Operators<ptensor>{},
      tenes::Operators<ptensor>{}, tenes::Operators<ptensor>{},
      tenes::itps::CorrelationParameter{});
}

}  // namespace

TEST_CASE("fermion real-time guard 1: mode = \"time\" is accepted") {
  using tenes::itps::PEPS_Parameters;
  auto in = rtg_make("mode = \"time\"\n", "");
  REQUIRE(in.peps_parameters.fermion == true);
  REQUIRE(in.peps_parameters.calcmode ==
          PEPS_Parameters::CalculationMode::time_evolution);
  REQUIRE(in.peps_parameters.Simple_Gauge_Fix == false);
  try {
    rtg_validate(in);
  } catch (std::exception const& e) {
    FAIL_CHECK(
        "fermion mode refused mode = \"time\": " << std::string(e.what()));
  }
}

TEST_CASE(
    "fermion real-time guard 2: mode = \"finite\" is refused with "
    "input_error naming finite-temperature") {
  using tenes::itps::PEPS_Parameters;
  auto in = rtg_make("mode = \"finite\"\n", "");
  REQUIRE(in.peps_parameters.fermion == true);
  REQUIRE(in.peps_parameters.calcmode ==
          PEPS_Parameters::CalculationMode::finite_temperature);
  REQUIRE(in.peps_parameters.Simple_Gauge_Fix == false);
  try {
    rtg_validate(in);
    FAIL_CHECK("fermion mode accepted mode = \"finite\"");
  } catch (tenes::input_error const& e) {
    const std::string msg(e.what());
    INFO("message: " << msg);
    CHECK(msg.find("finite-temperature") != std::string::npos);
  } catch (std::exception const& e) {
    FAIL_CHECK(
        "mode = \"finite\" was refused with an exception that is not "
        "tenes::input_error: "
        << std::string(e.what()));
  }
}

TEST_CASE(
    "fermion real-time guard 3: mode = \"time\" with Simple_Gauge_Fix = true "
    "is still refused") {
  using tenes::itps::PEPS_Parameters;
  auto in = rtg_make("mode = \"time\"\n", "gauge_fix = true\n");
  REQUIRE(in.peps_parameters.fermion == true);
  REQUIRE(in.peps_parameters.calcmode ==
          PEPS_Parameters::CalculationMode::time_evolution);
  REQUIRE(in.peps_parameters.Simple_Gauge_Fix == true);
  try {
    rtg_validate(in);
    FAIL_CHECK(
        "fermion mode accepted Simple_Gauge_Fix = true with mode = "
        "\"time\"");
  } catch (tenes::input_error const& e) {
    const std::string msg(e.what());
    INFO("message: " << msg);
    CHECK(msg.find("Simple_Gauge_Fix") != std::string::npos);
  } catch (std::exception const& e) {
    FAIL_CHECK(
        "Simple_Gauge_Fix = true was refused with an exception that "
        "is not tenes::input_error: "
        << std::string(e.what()));
  }
}
