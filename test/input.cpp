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

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../src/tensor.hpp"
#include "../src/util/string.hpp"
#include "../src/load_toml.hpp"
#include "../src/mpi.hpp"

auto parse_str(std::string const &str) -> decltype(cpptoml::parse_file("")) {
  std::stringstream ss;
  ss << str;
  cpptoml::parser p{ss};
  return p.parse();
}

TEST_CASE("input") {
  using namespace tenes;
  using ptensor = complex_tensor;

  SUBCASE("parameter_default") {
    INFO("parameter_default");
    auto toml = parse_str(R"([parameter])");

    PEPS_Parameters peps_parameters = gen_param(toml->get_table("parameter"));

    CHECK(peps_parameters.CHI == 2);

    CHECK(peps_parameters.num_simple_step == 0);
    CHECK(peps_parameters.Inverse_lambda_cut == 1e-12);

    CHECK(peps_parameters.num_full_step == 0);
    CHECK(peps_parameters.Inverse_Env_cut == 1e-12);
    CHECK(peps_parameters.Full_Inverse_precision == 1e-12);
    CHECK(peps_parameters.Full_Convergence_Epsilon == 1e-6);
    CHECK(peps_parameters.Full_max_iteration == 100);
    CHECK(peps_parameters.Full_Gauge_Fix == true);
    CHECK(peps_parameters.Full_Use_FastFullUpdate == true);

    CHECK(peps_parameters.Inverse_projector_cut == 1e-12);
    CHECK(peps_parameters.CTM_Convergence_Epsilon == 1e-6);
    CHECK(peps_parameters.Max_CTM_Iteration == 100);
    CHECK(peps_parameters.CTM_Projector_corner == true);
    CHECK(peps_parameters.Use_RSVD == false);
    CHECK(peps_parameters.RSVD_Oversampling_factor == 2.0);

    CHECK(peps_parameters.seed == 11);
  }

  SUBCASE("parameter") {
    INFO("parameter");
    auto toml = parse_str(R"(
[parameter]
[parameter.tensor]
save_dir = "checkpoint"
load_dir = "checkpoint"

[parameter.simple_update]
num_step = 1000
lambda_cutoff = 1e-10

[parameter.full_update]
num_step = 1
inverse_precision = 1e-10
convergence_epsilon = 1e-10
env_cutoff = 1e-10
iteration_max = 100
gauge_fix = false
fastfullupdate = false

[parameter.ctm]
dimension = 16
projector_cutoff = 1e-10
convergence_epsilon = 1e-8
iteration_max = 10
projector_corner = false
use_rsvd = true
rsvd_oversampling_factor = 3.0

[parameter.random]
seed = 42)");

    PEPS_Parameters peps_parameters = gen_param(toml->get_table("parameter"));

    CHECK(peps_parameters.CHI == 16);

    CHECK(peps_parameters.num_simple_step == 1000);
    CHECK(peps_parameters.Inverse_lambda_cut == 1e-10);

    CHECK(peps_parameters.num_full_step == 1);
    CHECK(peps_parameters.Inverse_Env_cut == 1e-10);
    CHECK(peps_parameters.Full_Inverse_precision == 1e-10);
    CHECK(peps_parameters.Full_Convergence_Epsilon == 1e-10);
    CHECK(peps_parameters.Full_max_iteration == 100);
    CHECK(peps_parameters.Full_Gauge_Fix == false);
    CHECK(peps_parameters.Full_Use_FastFullUpdate == false);

    CHECK(peps_parameters.Inverse_projector_cut == 1e-10);
    CHECK(peps_parameters.CTM_Convergence_Epsilon == 1e-8);
    CHECK(peps_parameters.Max_CTM_Iteration == 10);
    CHECK(peps_parameters.CTM_Projector_corner == false);
    CHECK(peps_parameters.Use_RSVD == true);
    CHECK(peps_parameters.RSVD_Oversampling_factor == 3.0);

    CHECK(peps_parameters.seed == 42);
  }

  SUBCASE("tensor") {
    INFO("tensor");
    auto toml = parse_str(R"(
[tensor]
L_sub = [4, 1]
skew = 2
[[tensor.unitcell]]
index = [0, 2]
physical_dim = 2
virtual_dim = [4, 3, 4, 3]
initial_state = [1.0, 0.0]
noise = 0.01
[[tensor.unitcell]]
index = [1, 3]
physical_dim = 3
virtual_dim = [4, 1, 4, 1]
initial_state = [0.0, 1.0]
noise = 0.01
    )");
    SquareLattice lattice = gen_lattice(toml->get_table("tensor"));
    CHECK(lattice.LX == 4);
    CHECK(lattice.LY == 1);
    CHECK(lattice.skew == 2);
  }

  SUBCASE("evolution") {
    {
      INFO("simple_update");
      auto toml = parse_str(R"(
[evolution]
[[evolution.simple]]
source_site = 0
source_leg = 2
dimensions = [2,2,2,4]
elements = """
0 0 0 0 1.0 0.0
"""
      )");
      const auto simple_updates = tenes::load_simple_updates<ptensor>(toml);
      CHECK(simple_updates[0].source_site == 0);
      CHECK(simple_updates[0].source_leg == 2);
      auto &op = simple_updates[0].op;
      CHECK(op.shape() == mptensor::Shape{2, 2, 2, 4});
      std::complex<double> v = 0.0;
      op.get_value({0, 0, 0, 0}, v);
      CHECK(std::real(v) == 1.0);
      CHECK(std::imag(v) == 0.0);
    }
    {
      INFO("full_update");
      auto toml = parse_str(R"(
[evolution]
[[evolution.full]]
source_site = 0
source_leg = 2
dimensions = [2,2,2,4]
elements = """
0 0 0 0 0.0 1.0
"""
      )");
      const auto full_updates = tenes::load_full_updates<ptensor>(toml);
      CHECK(full_updates[0].source_site == 0);
      CHECK(full_updates[0].source_leg == 2);
      auto &op = full_updates[0].op;
      CHECK(op.shape() == mptensor::Shape{2, 2, 2, 4});
      std::complex<double> v = 0.0;
      op.get_value({0, 0, 0, 0}, v);
      CHECK(std::real(v) == 0.0);
      CHECK(std::imag(v) == 1.0);
    }
  }

  SUBCASE("observable") {
    {
      INFO("onesite");
      auto toml = parse_str(R"(
[observable]
[[observable.onesite]]
group = 0
sites = []
dim = 2
elements = """
0 0 1.0 0.0
"""
      )");
      const int nsites = 2;
      const int nbody = 1;
      auto onesites = load_operators<ptensor>(toml, nsites, nbody, 0.0,
                                              "observable.onesite");
      for (int i = 0; i < 2; ++i) {
        auto const &on = onesites[i];
        CHECK(on.group == 0);
        CHECK(on.source_site == i);
        CHECK(on.is_onesite());
        CHECK(on.op.shape() == mptensor::Shape{2, 2});
        std::complex<double> v = 0.0;
        on.op.get_value({0, 0}, v);
        CHECK(std::real(v) == 1.0);
        CHECK(std::imag(v) == 0.0);
      }
    }
    {
      INFO("twosite");
      auto toml = parse_str(R"(
[observable]
[[observable.twosite]]
group = 0
dim = [2,2]
bonds = """
0 1 0
1 2 1
"""
elements = """
0 0 0 0 0.0 1.0
"""
      )");
      const int nsites = 2;
      const int nbody = 2;
      auto twosites = load_operators<ptensor>(toml, nsites, nbody, 0.0,
                                              "observable.twosite");
      for (int i = 0; i < 2; ++i) {
        auto const &on = twosites[i];
        CHECK(on.group == 0);
        CHECK(on.source_site == i);
        CHECK(on.dx == std::vector<int>{i + 1});
        CHECK(on.dy == std::vector<int>{i});
        CHECK(on.op.shape() == mptensor::Shape{2, 2, 2, 2});
        std::complex<double> v = 0.0;
        on.op.get_value({0, 0, 0, 0}, v);
        CHECK(std::real(v) == 0.0);
        CHECK(std::imag(v) == 1.0);
      }
    }
  }

  SUBCASE("correlation") {}
}
