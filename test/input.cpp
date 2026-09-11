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
#include "test_workdir.hpp"

#include <filesystem>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "../src/tensor.hpp"
#include "../src/mpi.hpp"
#include "../src/fermion/fermion_info.hpp"
#include "../src/fermion/fops.hpp"
#include "../src/util/string.hpp"
#include "../src/arpack_solver.hpp"
#include "../src/iTPS/load_toml.hpp"
#include "../src/iTPS/iTPS.hpp"
#include "../src/iTPS/main.hpp"
#include "../src/iTPS/transfer_matrix.hpp"

toml::value parse_str(std::string const &str) { return toml::parse_str(str); }

namespace tenes::itps {
struct iTPSTestAccessor {
  template <class tensor>
  static std::vector<tensor> const &Tn(iTPS<tensor> const &state) {
    return state.Tn;
  }

  template <class tensor>
  static tenes::fermion::FermionInfo const &finfo(iTPS<tensor> const &state) {
    return state.finfo;
  }
};
}  // namespace tenes::itps

TEST_CASE("input") {
  using namespace tenes;
  using namespace tenes::itps;
  using ptensor = complex_tensor;

  SUBCASE("parameter_default") {
    INFO("parameter_default");
    auto toml = parse_str(R"([parameter])");

    PEPS_Parameters peps_parameters = gen_param(toml.at("parameter"));

    CHECK(peps_parameters.CHI == 2);

    CHECK(peps_parameters.num_simple_step.size() == 1);
    CHECK(peps_parameters.num_simple_step[0] == 0);
    CHECK(peps_parameters.Inverse_lambda_cut == 1e-12);

    CHECK(peps_parameters.num_full_step.size() == 1);
    CHECK(peps_parameters.num_full_step[0] == 0);
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

    PEPS_Parameters peps_parameters = gen_param(toml.at("parameter"));

    CHECK(peps_parameters.CHI == 16);

    CHECK(peps_parameters.num_simple_step.size() == 1);
    CHECK(peps_parameters.num_simple_step[0] == 1000);
    CHECK(peps_parameters.Inverse_lambda_cut == 1e-10);

    CHECK(peps_parameters.num_full_step.size() == 1);
    CHECK(peps_parameters.num_full_step[0] == 1);
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

  SUBCASE("tau is read from its own section") {
    INFO("tau");
    auto toml = parse_str(R"(
[parameter]
[parameter.simple_update]
tau = 0.1
[parameter.full_update]
tau = 0.01
)");

    PEPS_Parameters peps_parameters = gen_param(toml.at("parameter"));

    REQUIRE(peps_parameters.tau_simple_step.size() == 1);
    CHECK(peps_parameters.tau_simple_step[0] == 0.1);
    REQUIRE(peps_parameters.tau_full_step.size() == 1);
    CHECK(peps_parameters.tau_full_step[0] == 0.01);
  }

  SUBCASE("saved mode string") {
    INFO("saved mode string");

    auto count_occurrences = [](std::string const &filename,
                                std::string const &key) {
      std::ifstream ifs(filename);
      std::string line;
      int n = 0;
      while (std::getline(ifs, line)) {
        if (line.find(key) != std::string::npos) {
          ++n;
        }
      }
      return n;
    };

    struct {
      PEPS_Parameters::CalculationMode mode;
      const char *name;
    } cases[] = {
        {PEPS_Parameters::CalculationMode::ground_state, "ground state"},
        {PEPS_Parameters::CalculationMode::time_evolution, "time evolution"},
        {PEPS_Parameters::CalculationMode::finite_temperature,
         "finite temperature"},
    };

    for (auto const &c : cases) {
      PEPS_Parameters peps_parameters;
      peps_parameters.calcmode = c.mode;
      const std::string filename = "output_parameters_test.dat";
      peps_parameters.save(filename.c_str());

      CHECK(count_occurrences(filename, std::string("mode = ") + c.name) == 1);
      CHECK(count_occurrences(filename, "ground state") +
                count_occurrences(filename, "time evolution") +
                count_occurrences(filename, "finite temperature") ==
            1);
    }
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
    SquareLattice lattice = gen_lattice(toml.at("tensor"));
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
group = 0
source_site = 0
source_leg = 2
dimensions = [2,2,2,4]
elements = """
0 0 0 0 1.0 0.0
"""
      )");
      const auto simple_updates =
          tenes::itps::load_simple_updates<ptensor>(toml, MPI_COMM_WORLD);
      CHECK(simple_updates[0].source_site == 0);
      CHECK(simple_updates[0].source_leg == 2);
      CHECK(simple_updates[0].group == 0);
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
group = 0
source_site = 0
source_leg = 2
dimensions = [2,2,2,4]
elements = """
0 0 0 0 0.0 1.0
"""
      )");
      const auto full_updates =
          tenes::itps::load_full_updates<ptensor>(toml, MPI_COMM_WORLD);
      CHECK(full_updates[0].source_site == 0);
      CHECK(full_updates[0].source_leg == 2);
      CHECK(full_updates[0].group == 0);
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
      auto onesites = load_operators<ptensor>(toml, MPI_COMM_WORLD, nsites,
                                              nbody, 0.0, "observable.onesite");
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
      auto twosites = load_operators<ptensor>(toml, MPI_COMM_WORLD, nsites,
                                              nbody, 0.0, "observable.twosite");
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

  SUBCASE("iTPS without evolution operators") {
    INFO("iTPS without evolution operators");
    // measurement-only setup: no [[evolution.simple]] / [[evolution.full]]
    auto toml = parse_str(R"(
[tensor]
L_sub = [2, 2]
[[tensor.unitcell]]
index = []
physical_dim = 2
virtual_dim = 2
initial_state = [1.0, 0.0]
noise = 0.01
    )");
    PEPS_Parameters peps_parameters;
    peps_parameters.print_level = PrintLevel::none;
    peps_parameters.outdir = "output_itps_without_evolution";
    SquareLattice lattice = gen_lattice(toml.at("tensor"));

    CHECK_NOTHROW(iTPS<ptensor>(
        MPI_COMM_WORLD, peps_parameters, lattice, EvolutionOperators<ptensor>{},
        EvolutionOperators<ptensor>{}, Operators<ptensor>{},
        Operators<ptensor>{}, Operators<ptensor>{}, CorrelationParameter{},
        TransferMatrix_Parameters{}));
  }

  SUBCASE("fermion initialization masks odd-total Tn entries") {
    INFO("fermion initialization masks odd-total Tn entries");
    // L_sub = [2, 1] is left as-is: this subcase constructs iTPS<ptensor>
    // directly and never calls validate_fermion_constraints (only
    // itps_main does, in main.cpp), so the self-neighbour guard (a site of
    // this skew-0 one-row cell is its own vertical neighbour) is never on
    // this subcase's path.
    auto toml = parse_str(R"(
[tensor]
L_sub = [2, 1]
[[tensor.unitcell]]
index = []
physical_dim = 2
virtual_dim = 2
parity = [0, 1]
noise = 0.01
    )");
    PEPS_Parameters peps_parameters;
    peps_parameters.fermion = true;
    peps_parameters.print_level = PrintLevel::none;
    SquareLattice lattice = gen_lattice(toml.at("tensor"));
    peps_parameters.phys_parity = gen_phys_parity(toml.at("tensor"), lattice);

    iTPS<ptensor> state(MPI_COMM_WORLD, peps_parameters, lattice,
                        EvolutionOperators<ptensor>{},
                        EvolutionOperators<ptensor>{}, Operators<ptensor>{},
                        Operators<ptensor>{}, Operators<ptensor>{},
                        CorrelationParameter{}, TransferMatrix_Parameters{});
    const auto &fi = iTPSTestAccessor::finfo(state);
    REQUIRE(fi.enabled);
    const auto &tensors = iTPSTestAccessor::Tn(state);
    REQUIRE(tensors.size() == static_cast<std::size_t>(lattice.N_UNIT));
    for (int site = 0; site < lattice.N_UNIT; ++site) {
      auto ft = tenes::fermion::wrap_Tn(tensors[site], fi, site);
      CHECK(tenes::fermion::parity_violation(ft) == doctest::Approx(0.0));
    }
  }

  SUBCASE("correlation") {}

  SUBCASE("correlation_length eigensolver") {
    INFO("correlation_length eigensolver");
    auto toml_default = parse_str(R"([correlation_length])");
    auto p_default = gen_transfer_matrix_parameter(
        toml_default.at("correlation_length"), "correlation_length");
    CHECK(p_default.eigensolver == TransferMatrixEigensolver::automatic);

    auto toml_builtin = parse_str(R"(
[correlation_length]
eigensolver = "builtin"
)");
    auto p_builtin = gen_transfer_matrix_parameter(
        toml_builtin.at("correlation_length"), "correlation_length");
    CHECK(p_builtin.eigensolver == TransferMatrixEigensolver::builtin);

    auto toml_arpack = parse_str(R"(
[correlation_length]
eigensolver = "arpack"
)");
    if (tenes::arpack_available()) {
      auto p_arpack = gen_transfer_matrix_parameter(
          toml_arpack.at("correlation_length"), "correlation_length");
      CHECK(p_arpack.eigensolver == TransferMatrixEigensolver::arpack);
    } else {
      CHECK_THROWS_AS(
          gen_transfer_matrix_parameter(toml_arpack.at("correlation_length"),
                                        "correlation_length"),
          tenes::input_error);
    }

    auto toml_bad = parse_str(R"(
[correlation_length]
eigensolver = "lapack"
)");
    CHECK_THROWS_AS(
        gen_transfer_matrix_parameter(toml_bad.at("correlation_length"),
                                      "correlation_length"),
        tenes::input_error);
  }

  SUBCASE("correlation_length arnoldi defaults are automatic") {
    INFO("correlation_length arnoldi defaults are automatic");
    TransferMatrix_Parameters p;
    CHECK(p.arnoldi_maxdim == 0);      // 0 means automatic
    CHECK(p.arnoldi_restartdim == 0);  // 0 means automatic
    CHECK(p.arnoldi_maxiter == 0);     // 0 means automatic

    // ARPACK relies on restarts: max(2 * num_eigvals + 1, 25)
    CHECK(effective_arnoldi_maxdim(0, 4, true) == 25);
    CHECK(effective_arnoldi_maxdim(0, 12, true) == 25);
    CHECK(effective_arnoldi_maxdim(0, 15, true) == 31);
    // builtin solves in one large sweep: max(2 * num_eigvals + 1, 50)
    CHECK(effective_arnoldi_maxdim(0, 4, false) == 50);
    CHECK(effective_arnoldi_maxdim(0, 30, false) == 61);
    // explicit values are used as-is for both solvers
    CHECK(effective_arnoldi_maxdim(40, 4, true) == 40);
    CHECK(effective_arnoldi_maxdim(8, 15, false) == 8);

    // automatic maxiter: 10 restarts for ARPACK, none for builtin
    CHECK(effective_arnoldi_maxiter(0, true) == 10);
    CHECK(effective_arnoldi_maxiter(0, false) == 1);
    CHECK(effective_arnoldi_maxiter(3, true) == 3);
    CHECK(effective_arnoldi_maxiter(3, false) == 3);

    // automatic: max(num_eigvals + 1, maxdim / 2)
    CHECK(effective_arnoldi_restartdim(0, 4, 50) == 25);
    CHECK(effective_arnoldi_restartdim(0, 24, 31) == 25);
    // explicit values are used as-is
    CHECK(effective_arnoldi_restartdim(20, 4, 50) == 20);
  }

  SUBCASE("fermion parity input loads") {
    INFO("fermion parity input loads");
    // L_sub = [1, 1] is left as-is: this subcase only exercises
    // gen_param/gen_lattice/gen_phys_parity (the raw TOML parsing) and
    // never calls validate_fermion_constraints, so the self-neighbour guard
    // (which refuses a 1x1 cell) never runs here and the cell still loads
    // cleanly.
    auto param_toml = parse_str(R"(
[parameter]
[parameter.general]
fermion = true
)");
    auto tensor_toml = parse_str(R"(
[tensor]
L_sub = [1, 1]
[[tensor.unitcell]]
index = [0]
physical_dim = 2
virtual_dim = 2
parity = [0, 1]
)");
    PEPS_Parameters peps_parameters = gen_param(param_toml.at("parameter"));
    SquareLattice lattice = gen_lattice(tensor_toml.at("tensor"));
    peps_parameters.phys_parity =
        gen_phys_parity(tensor_toml.at("tensor"), lattice);
    CHECK(peps_parameters.fermion == true);
    REQUIRE(peps_parameters.phys_parity.size() == 1);
    CHECK(peps_parameters.phys_parity[0] == std::vector<bool>{false, true});
  }

  SUBCASE("fermion accepts mean-field environment") {
    INFO("fermion accepts mean-field environment");
    auto param_toml = parse_str(R"(
[parameter]
[parameter.general]
fermion = true
[parameter.ctm]
meanfield_env = true
)");
    // L_sub = [2, 2]: with a 1x1 cell the self-neighbour guard
    // (docs/superpowers/specs/2026-09-11-fermion-skew-guard-contract.md
    // section 2.1) would throw before MeanField_Env is even looked at, so
    // the subcase could not tell that the mean-field environment itself is
    // accepted. A 2x2 cell with a single site definition broadcast via
    // index = [] clears that guard; with the fermionic mean-field
    // measurement in place, MeanField_Env=true is a supported combination
    // and nothing else in this input is guarded.
    auto tensor_toml = parse_str(R"(
[tensor]
L_sub = [2, 2]
[[tensor.unitcell]]
index = []
physical_dim = 2
virtual_dim = 2
parity = [0, 1]
)");
    PEPS_Parameters peps_parameters = gen_param(param_toml.at("parameter"));
    SquareLattice lattice = gen_lattice(tensor_toml.at("tensor"));
    peps_parameters.phys_parity =
        gen_phys_parity(tensor_toml.at("tensor"), lattice);
    CHECK(peps_parameters.MeanField_Env == true);
    CHECK_NOTHROW(validate_fermion_constraints(
        peps_parameters, lattice, EvolutionOperators<ptensor>{},
        EvolutionOperators<ptensor>{}, Operators<ptensor>{},
        Operators<ptensor>{}, Operators<ptensor>{}, CorrelationParameter{}));
  }

  SUBCASE("fermion rejects odd one-site operator") {
    INFO("fermion rejects odd one-site operator");
    auto param_toml = parse_str(R"(
[parameter]
[parameter.general]
fermion = true
)");
    // L_sub = [2, 2]: same reasoning as "fermion accepts mean-field
    // environment" above -- a 1x1 cell would trip the self-neighbour guard
    // before the parity-odd one-site operator check this subcase is named
    // for.
    auto tensor_toml = parse_str(R"(
[tensor]
L_sub = [2, 2]
[[tensor.unitcell]]
index = []
physical_dim = 2
virtual_dim = 2
parity = [0, 1]
)");
    auto observable_toml = parse_str(R"(
[observable]
[[observable.onesite]]
group = 0
sites = [0]
dim = 2
elements = """
0 1 1.0 0.0
"""
)");
    PEPS_Parameters peps_parameters = gen_param(param_toml.at("parameter"));
    SquareLattice lattice = gen_lattice(tensor_toml.at("tensor"));
    peps_parameters.phys_parity =
        gen_phys_parity(tensor_toml.at("tensor"), lattice);
    auto onesite = load_operators<ptensor>(observable_toml, MPI_COMM_WORLD, 4,
                                           1, 0.0, "observable.onesite");
    CHECK_THROWS_AS(
        validate_fermion_constraints(
            peps_parameters, lattice, EvolutionOperators<ptensor>{},
            EvolutionOperators<ptensor>{}, onesite, Operators<ptensor>{},
            Operators<ptensor>{}, CorrelationParameter{}),
        tenes::input_error);
  }

  SUBCASE("fermion odd operator is rejected through itps_main load path") {
    INFO("fermion odd operator is rejected through itps_main load path");
    const std::string input_filename =
        "test_input_fermion_odd_operator_main_path.toml";
    const std::string outdir =
        "output_test_input_fermion_odd_operator_main_path";
    {
      std::ofstream ofs(input_filename);
      // L_sub = [2, 2]: a 1x1 cell trips the self-neighbour guard
      // (docs/superpowers/specs/2026-09-11-fermion-skew-guard-contract.md
      // section 2.1) before reaching the parity-odd one-site operator check
      // this subcase exercises; a single site definition broadcast via
      // index = [] keeps this a minimal fixture while clearing that guard.
      ofs << R"(
[parameter]
[parameter.general]
is_real = true
fermion = true
output = ")"
          << outdir << R"("

[tensor]
L_sub = [2, 2]
[[tensor.unitcell]]
index = []
physical_dim = 2
virtual_dim = 2
parity = [0, 1]

[observable]
[[observable.onesite]]
name = "odd"
group = 0
sites = [0]
dim = 2
elements = """
0 1 1.0 0.0
"""

[evolution]
)";
    }

    try {
      tenes::itps::itps_main(input_filename, MPI_COMM_WORLD, PrintLevel::none);
      FAIL("fermion odd operator was accepted through itps_main");
    } catch (const tenes::input_error &e) {
      CHECK(std::string(e.what()).find("parity-odd one-site operators") !=
            std::string::npos);
    }
    std::remove(input_filename.c_str());
  }
}

TEST_CASE("identity gates complete the bonds no Hamiltonian term gates") {
  using namespace tenes;
  using namespace tenes::itps;
  using ptensor = complex_tensor;
  MPI_Comm comm = MPI_COMM_WORLD;

  auto make_lattice = [&](std::string const &vdim) {
    auto tensor_toml = parse_str(R"(
[tensor]
L_sub = [2, 2]
[[tensor.unitcell]]
index = []
physical_dim = 2
virtual_dim = )" + vdim + R"(
)");
    return gen_lattice(tensor_toml.at("tensor"));
  };

  auto horizontal_gate = [&](int site, int group) {
    ptensor op(comm, mptensor::Shape(2, 2, 2, 2));
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        op.set_value(mptensor::Index(i, j, i, j), 0.5);
      }
    }
    return make_twosite_EvolutionOperator<ptensor>(site, 2, group, op);
  };

  auto count_leg = [](EvolutionOperators<ptensor> const &ops, int leg) {
    int n = 0;
    for (auto const &op : ops) {
      if (op.is_twosite() && op.source_leg == leg) ++n;
    }
    return n;
  };

  SUBCASE("ungated vertical bonds with D > 1 receive identity gates") {
    auto lattice = make_lattice("2");
    EvolutionOperators<ptensor> ops;
    for (int s = 0; s < 4; ++s) ops.push_back(horizontal_gate(s, 0));

    auto completed = complete_ungated_bonds<ptensor>(ops, lattice, comm);

    CHECK(completed.size() == 8);
    CHECK(count_leg(completed, 2) == 4);
    CHECK(count_leg(completed, 1) == 4);
    for (auto const &op : completed) {
      if (op.source_leg != 1) continue;
      CHECK(op.group == 0);
      for (int i1 = 0; i1 < 2; ++i1) {
        for (int i2 = 0; i2 < 2; ++i2) {
          for (int o1 = 0; o1 < 2; ++o1) {
            for (int o2 = 0; o2 < 2; ++o2) {
              typename ptensor::value_type v;
              op.op.get_value(mptensor::Index(i1, i2, o1, o2), v);
              const double expected = (i1 == o1 && i2 == o2) ? 1.0 : 0.0;
              CHECK(std::abs(v - expected) < 1e-15);
            }
          }
        }
      }
    }
  }

  SUBCASE("a D = 1 ungated leg needs nothing") {
    auto lattice = make_lattice("[2, 1, 2, 1]");
    EvolutionOperators<ptensor> ops;
    for (int s = 0; s < 4; ++s) ops.push_back(horizontal_gate(s, 0));
    auto completed = complete_ungated_bonds<ptensor>(ops, lattice, comm);
    CHECK(completed.size() == 4);
  }

  SUBCASE("a fully gated cell is returned unchanged") {
    auto lattice = make_lattice("2");
    EvolutionOperators<ptensor> ops;
    for (int s = 0; s < 4; ++s) ops.push_back(horizontal_gate(s, 0));
    for (int s = 0; s < 4; ++s) {
      ptensor op(comm, mptensor::Shape(2, 2, 2, 2));
      op.set_value(mptensor::Index(0, 0, 0, 0), 1.0);
      ops.push_back(make_twosite_EvolutionOperator<ptensor>(s, 1, 0, op));
    }
    auto completed = complete_ungated_bonds<ptensor>(ops, lattice, comm);
    CHECK(completed.size() == 8);
  }

  SUBCASE("a bond gated from the other end counts as gated") {
    // the vertical bond of site 0 (leg 1, up) is the same bond as the
    // bottom leg (3) of its upper neighbour
    auto lattice = make_lattice("2");
    EvolutionOperators<ptensor> ops;
    for (int s = 0; s < 4; ++s) ops.push_back(horizontal_gate(s, 0));
    for (int s = 0; s < 4; ++s) {
      ptensor op(comm, mptensor::Shape(2, 2, 2, 2));
      op.set_value(mptensor::Index(0, 0, 0, 0), 1.0);
      ops.push_back(make_twosite_EvolutionOperator<ptensor>(s, 3, 0, op));
    }
    auto completed = complete_ungated_bonds<ptensor>(ops, lattice, comm);
    CHECK(completed.size() == 8);
  }

  SUBCASE("a cell with only one-site gates is completed on the given comm") {
    // groups is empty here, so the generated gates fall back to group 0.
    // The communicator cannot be read off the (absent) two-site operators
    // either: it has to come from the caller.
    auto lattice = make_lattice("2");
    EvolutionOperators<ptensor> ops;
    for (int s = 0; s < 4; ++s) {
      ptensor op(comm, mptensor::Shape(2, 2));
      op.set_value(mptensor::Index(0, 0), 1.0);
      op.set_value(mptensor::Index(1, 1), 1.0);
      ops.push_back(make_onesite_EvolutionOperator<ptensor>(s, 0, op));
    }

    auto completed = complete_ungated_bonds<ptensor>(ops, lattice, comm);

    // every one of the 8 bonds of the 2x2 cell is ungated
    CHECK(completed.size() == 12);
    CHECK(count_leg(completed, 1) == 4);
    CHECK(count_leg(completed, 2) == 4);
    for (auto const &op : completed) {
      if (!op.is_twosite()) continue;
      CHECK(op.group == 0);
      CHECK(op.op.get_comm() == comm);
    }
  }

  SUBCASE("identity gates follow every group that is present") {
    auto lattice = make_lattice("2");
    EvolutionOperators<ptensor> ops;
    for (int s = 0; s < 4; ++s) ops.push_back(horizontal_gate(s, 0));
    for (int s = 0; s < 4; ++s) ops.push_back(horizontal_gate(s, 1));
    auto completed = complete_ungated_bonds<ptensor>(ops, lattice, comm);
    CHECK(completed.size() == 16);
    int g0 = 0, g1 = 0;
    for (auto const &op : completed) {
      if (op.source_leg == 1) (op.group == 0 ? g0 : g1)++;
    }
    CHECK(g0 == 4);
    CHECK(g1 == 4);
  }
}

namespace {

std::string fermion_cell_toml() {
  return R"(
[tensor]
L_sub = [2, 2]
skew = 0
[[tensor.unitcell]]
index = []
physical_dim = 2
virtual_dim = 2
parity = [0, 1]
noise = 0.01
)";
}

}  // namespace

TEST_CASE("fermion mode accepts the full update") {
  using namespace tenes;
  using namespace tenes::itps;
  using ptensor = complex_tensor;
  MPI_Comm comm = MPI_COMM_WORLD;

  auto tensor_toml = parse_str(fermion_cell_toml());
  SquareLattice lattice = gen_lattice(tensor_toml.at("tensor"));

  auto param_toml = parse_str(R"(
[parameter]
[parameter.general]
fermion = true
[parameter.full_update]
tau = 0.01
num_step = 1
)");
  PEPS_Parameters peps_parameters = gen_param(param_toml.at("parameter"));
  peps_parameters.phys_parity = gen_phys_parity(tensor_toml.at("tensor"), lattice);

  // Preconditions: the guard this subcase is about has to be reachable at
  // all. Without these the CHECK_NOTHROW below would also pass for an input
  // that simply has no full update in it.
  REQUIRE(peps_parameters.fermion == true);
  REQUIRE(peps_parameters.num_full_step.size() == 1);
  REQUIRE(peps_parameters.num_full_step[0] == 1);

  SUBCASE("a positive full-update step count is no longer refused") {
    INFO("a positive full-update step count is no longer refused");
    CHECK_NOTHROW(validate_fermion_constraints(
        peps_parameters, lattice, EvolutionOperators<ptensor>{},
        EvolutionOperators<ptensor>{}, Operators<ptensor>{},
        Operators<ptensor>{}, Operators<ptensor>{}, CorrelationParameter{}));
  }

  SUBCASE("a parity-odd full-update gate is still refused") {
    INFO("a parity-odd full-update gate is still refused");
    // Lifting the "full update" guard must not take the parity check on the
    // full-update gates with it. (0, 0, 0, 1) has one odd leg, so with the
    // physical ledger [0, 1] the gate is parity odd.
    ptensor op(comm, mptensor::Shape(2, 2, 2, 2));
    op.set_value(mptensor::Index(0, 0, 0, 1), 1.0);
    EvolutionOperators<ptensor> full_updates{
        make_twosite_EvolutionOperator<ptensor>(0, 2, 0, op)};
    CHECK_THROWS_AS(
        validate_fermion_constraints(
            peps_parameters, lattice, EvolutionOperators<ptensor>{},
            full_updates, Operators<ptensor>{}, Operators<ptensor>{},
            Operators<ptensor>{}, CorrelationParameter{}),
        tenes::input_error);
  }
}

TEST_CASE("fermion mode refuses the mean-field environment with a full update") {
  using namespace tenes;
  using namespace tenes::itps;
  using ptensor = complex_tensor;

  auto tensor_toml = parse_str(fermion_cell_toml());
  SquareLattice lattice = gen_lattice(tensor_toml.at("tensor"));

  auto make_parameters = [&](bool meanfield) {
    auto param_toml = parse_str(R"(
[parameter]
[parameter.general]
fermion = true
[parameter.full_update]
tau = 0.01
num_step = 1
[parameter.ctm]
dimension = 4
)");
    PEPS_Parameters p = gen_param(param_toml.at("parameter"));
    p.phys_parity = gen_phys_parity(tensor_toml.at("tensor"), lattice);
    p.print_level = PrintLevel::none;
    p.outdir = "output_test_input_fermion_full_meanfield";
    p.MeanField_Env = meanfield;
    return p;
  };

  auto build = [&](PEPS_Parameters const &p) {
    return iTPS<ptensor>(MPI_COMM_WORLD, p, lattice,
                         EvolutionOperators<ptensor>{},
                         EvolutionOperators<ptensor>{}, Operators<ptensor>{},
                         Operators<ptensor>{}, Operators<ptensor>{},
                         CorrelationParameter{}, TransferMatrix_Parameters{});
  };

  // Precondition / control: the very same configuration without the
  // mean-field environment must be accepted, otherwise the CHECK_THROWS
  // below would be passing because of the full update rather than because
  // of meanfield_env.
  auto plain = make_parameters(false);
  REQUIRE(plain.num_full_step[0] == 1);
  REQUIRE(plain.MeanField_Env == false);
  CHECK_NOTHROW(build(plain));

  auto meanfield = make_parameters(true);
  REQUIRE(meanfield.MeanField_Env == true);
  REQUIRE(meanfield.num_full_step[0] == 1);
  CHECK_THROWS_AS(build(meanfield), tenes::input_error);
}

// ===== Fermion mode and the shape of the unit cell =========================
//
// docs/superpowers/specs/2026-09-11-fermion-skew-guard-contract.md section
// 2.1 (evidence: docs/superpowers/notes/2026-09-11-fermion-skew-revisit.md).
// validate_fermion_constraints refuses a fermion-mode cell if and only if
// some site is its own nearest neighbour through the periodic + skew
// boundary. With T(x, y) = T(x + skew, y + LY) that is exactly LX == 1
// (horizontal), or LY == 1 with skew = 0 mod LX (vertical). Skewed cells as
// such are accepted: the "skew breaks the fermionic signs" measurement
// behind the old refusal was retracted (see the note), and
// test/fermion/skew_unfold.cpp now pins skewed cells to their unfolded
// skew-0 equivalents.
//
// The lattice keeps the sign of the input skew (skew % LX in C++ semantics),
// so [2,1] skew -1 holds skew = -1 and must be accepted like skew = 1.
//
// An input skew that is a non-zero multiple of LX (e.g. [2,1] skew 2) is the
// skew-0 lattice (contract section 8.1), so such cells reach the guard like
// any other and are among the cases below. Before that fix the SquareLattice
// constructor divided by zero (lcm(LX, 0) / 0) while building them: a SIGFPE
// on x86-64, a silent LY_noskew = 0 on arm64. The next two test cases pin
// the constructor itself.

namespace {

//! LY_noskew of an [lx, ly] cell with the given skew, from arithmetic: a
//! skew r = skew mod lx in [1, lx) repeats after lcm(lx, r) / r rows of
//! cells, a multiple of lx is no skew at all.
int expected_ly_noskew(int lx, int ly, int skew) {
  const int r = ((skew % lx) + lx) % lx;
  return r == 0 ? ly : ly * (std::lcm(lx, r) / r);
}

//! The site of an [lx, ly] cell with the given skew at global position
//! (x, y), from T(x, y) = T(x + skew, y + ly): (x, y + k ly) holds the site
//! at (x - k skew, y).
int skew_site_at(int lx, int ly, int skew, int x, int y) {
  const int k = (y >= 0 ? y : y - ly + 1) / ly;
  const int xs = (((x - skew * k) % lx) + lx) % lx;
  const int ys = ((y % ly) + ly) % ly;
  return xs + lx * ys;
}

}  // namespace

TEST_CASE("SquareLattice treats a skew that is a multiple of LX as no skew") {
  using tenes::SquareLattice;
  // Contract section 8.1. Every one of these used to divide by zero in the
  // constructor; on arm64 that left LY_noskew = N_UNIT_noskew = 0, which is
  // what the checks below see on the unfixed code there.
  const int cells[][3] = {{2, 2, 2},  {2, 1, 2}, {2, 1, -2}, {3, 1, 3},
                          {3, 2, 6},  {1, 2, 1}, {1, 1, 1},  {1, 2, 5},
                          {2, 2, -4}, {4, 1, 8}, {3, 3, -3}};
  for (const auto &c : cells) {
    const int lx = c[0];
    const int ly = c[1];
    const int skew = c[2];
    INFO("L_sub = [" << lx << ", " << ly << "], skew = " << skew);
    const SquareLattice lattice(lx, ly, skew);
    const SquareLattice plain(lx, ly, 0);
    CHECK(lattice.skew == 0);
    CHECK(lattice.LX == lx);
    CHECK(lattice.LY == ly);
    CHECK(lattice.N_UNIT == lx * ly);
    CHECK(lattice.LX_noskew == lx);
    CHECK(lattice.LY_noskew == ly);
    CHECK(lattice.N_UNIT_noskew == lx * ly);
    int neighbor_mismatch = 0;
    int other_mismatch = 0;
    int parity_mismatch = 0;
    for (int i = 0; i < lattice.N_UNIT; ++i) {
      for (int leg = 0; leg < 4; ++leg) {
        neighbor_mismatch += lattice.neighbor(i, leg) != plain.neighbor(i, leg);
      }
      for (int dx = -2; dx <= 2; ++dx) {
        for (int dy = -2; dy <= 2; ++dy) {
          other_mismatch += lattice.other(i, dx, dy) != plain.other(i, dx, dy);
        }
      }
      parity_mismatch += lattice.parity(i) != plain.parity(i);
    }
    int index_mismatch = 0;
    for (int x = -2 * lx; x < 3 * lx; ++x) {
      for (int y = -3 * ly; y < 3 * ly; ++y) {
        index_mismatch += lattice.index(x, y) != plain.index(x, y);
      }
    }
    CHECK(neighbor_mismatch == 0);
    CHECK(other_mismatch == 0);
    CHECK(index_mismatch == 0);
    CHECK(parity_mismatch == 0);
  }

  // The same through the input path.
  auto tensor_toml = parse_str(R"(
[tensor]
L_sub = [2, 2]
skew = 2
[[tensor.unitcell]]
index = []
physical_dim = 2
virtual_dim = 2
)");
  const SquareLattice from_input =
      tenes::itps::gen_lattice(tensor_toml.at("tensor"));
  CHECK(from_input.skew == 0);
  CHECK(from_input.LY_noskew == 2);
  CHECK(from_input.N_UNIT_noskew == 4);
}

TEST_CASE("SquareLattice keeps every other skew") {
  using tenes::SquareLattice;
  // Contract section 8.1: the fix must not touch these. The member keeps
  // skew % LX with the sign of the input (the guard test below relies on
  // it), and the neighbour and index maps follow T(x, y) = T(x + skew, y + LY).
  const int cells[][3] = {{2, 2, 1}, {2, 1, 1},  {2, 1, -1}, {3, 1, 1},
                          {3, 1, 2}, {3, 1, -1}, {3, 1, 4},  {3, 3, -4},
                          {4, 1, 2}, {4, 3, 6},  {2, 2, 7},  {4, 2, -2}};
  const int dx[4] = {-1, 0, 1, 0};
  const int dy[4] = {0, 1, 0, -1};
  for (const auto &c : cells) {
    const int lx = c[0];
    const int ly = c[1];
    const int skew = c[2];
    INFO("L_sub = [" << lx << ", " << ly << "], skew = " << skew);
    const SquareLattice lattice(lx, ly, skew);
    const int ly_noskew = expected_ly_noskew(lx, ly, skew);
    CHECK(lattice.skew == skew % lx);
    CHECK(lattice.LX_noskew == lx);
    CHECK(lattice.LY_noskew == ly_noskew);
    CHECK(lattice.N_UNIT_noskew == lx * ly_noskew);
    int neighbor_mismatch = 0;
    for (int i = 0; i < lattice.N_UNIT; ++i) {
      for (int leg = 0; leg < 4; ++leg) {
        neighbor_mismatch +=
            lattice.neighbor(i, leg) !=
            skew_site_at(lx, ly, skew, i % lx + dx[leg], i / lx + dy[leg]);
      }
    }
    int index_mismatch = 0;
    for (int x = -2 * lx; x < 3 * lx; ++x) {
      for (int y = -3 * ly_noskew; y < 3 * ly_noskew; ++y) {
        index_mismatch +=
            lattice.index(x, y) != skew_site_at(lx, ly, skew, x, y);
      }
    }
    CHECK(neighbor_mismatch == 0);
    CHECK(index_mismatch == 0);
  }
}

namespace {

std::string fermion_shaped_cell_toml(int lx, int ly, int skew) {
  std::ostringstream os;
  os << "[tensor]\n"
     << "L_sub = [" << lx << ", " << ly << "]\n"
     << "skew = " << skew << "\n"
     << R"([[tensor.unitcell]]
index = []
physical_dim = 2
virtual_dim = 2
parity = [0, 1]
noise = 0.01
)";
  return os.str();
}

//! True iff `text` contains a match of the ECMAScript regular expression
//! `pattern`, ignoring case.
bool contains_icase(std::string const &text, std::string const &pattern) {
  return std::regex_search(
      text, std::regex(pattern, std::regex::ECMAScript | std::regex::icase));
}

//! True iff `text` contains `word` followed, after characters that are
//! neither digits nor minus signs, by the given integers in order (each one
//! ending at a non-digit). "L_sub = [3, 1]", "L_sub=[3,1]" and
//! "tensor.L_sub (3 x 1)" all name L_sub 3 1; the check is case-insensitive.
bool names_numbers_after(std::string const &text, std::string const &word,
                         std::vector<int> const &numbers) {
  std::string pattern = word;
  for (std::size_t i = 0; i < numbers.size(); ++i) {
    pattern += (i == 0 ? "[^0-9-]*" : "[^0-9-]+");
    pattern += std::to_string(numbers[i]);
  }
  pattern += "(?![0-9])";
  return contains_icase(text, pattern);
}

}  // namespace

TEST_CASE("fermion mode refuses exactly the cells with a self-neighbour site") {
  using namespace tenes;
  using namespace tenes::itps;
  using ptensor = complex_tensor;

  struct cell {
    int lx;
    int ly;
    int skew;
  };

  auto validate = [](cell c) {
    auto tensor_toml = parse_str(fermion_shaped_cell_toml(c.lx, c.ly, c.skew));
    SquareLattice lattice = gen_lattice(tensor_toml.at("tensor"));
    auto param_toml = parse_str(R"(
[parameter]
[parameter.general]
fermion = true
)");
    PEPS_Parameters peps_parameters = gen_param(param_toml.at("parameter"));
    peps_parameters.phys_parity =
        gen_phys_parity(tensor_toml.at("tensor"), lattice);
    // Premises: the cell is what the input says, the lattice holds the skew
    // with the sign of the input, and the rest of the input is a valid
    // fermion input (parity metadata on every site), so that the only thing
    // left for the guard to judge is the shape.
    REQUIRE(lattice.LX == c.lx);
    REQUIRE(lattice.LY == c.ly);
    REQUIRE(lattice.skew == c.skew % c.lx);
    // The cell the guard judges is the whole cell (contract section 8.1; on
    // the unfixed code a skew that is a multiple of LX left LY_noskew = 0).
    CHECK(lattice.LY_noskew == expected_ly_noskew(c.lx, c.ly, c.skew));
    CHECK(lattice.N_UNIT_noskew ==
          c.lx * expected_ly_noskew(c.lx, c.ly, c.skew));
    REQUIRE(peps_parameters.fermion == true);
    REQUIRE(peps_parameters.phys_parity.size() ==
            static_cast<std::size_t>(c.lx * c.ly));
    validate_fermion_constraints(
        peps_parameters, lattice, EvolutionOperators<ptensor>{},
        EvolutionOperators<ptensor>{}, Operators<ptensor>{},
        Operators<ptensor>{}, Operators<ptensor>{}, CorrelationParameter{});
  };

  SUBCASE("cells in which no site is its own neighbour are accepted") {
    const cell accepted[] = {
        // Newly accepted: skewed cells with both sides >= 2,
        {2, 2, 1},
        {3, 2, 1},
        {3, 3, 2},
        {2, 3, 1},
        // one-row cells whose skew is not a multiple of LX (every bond of
        // such a row goes to a different site; [2,1] skew 1 is what
        // tenes_simple builds for a square lattice with W = 1),
        {2, 1, 1},
        {2, 1, -1},
        {3, 1, 1},
        {3, 1, 2},
        {3, 1, -1},
        {4, 1, 2},
        // and cells with |skew| >= LX.
        {2, 2, 7},
        {3, 1, 4},
        // Accepted before and still accepted.
        {2, 2, 0},
        {3, 3, 0},
        // A skew that is a multiple of LX on a cell two sites high: the
        // skew-0 cell (contract section 8.1).
        {2, 2, 2},
        {3, 2, 6},
    };
    for (const cell c : accepted) {
      INFO("L_sub = [" << c.lx << ", " << c.ly << "], skew = " << c.skew);
      CHECK_NOTHROW(validate(c));
    }
  }

  SUBCASE("cells with a self-neighbour site are refused, naming the cause") {
    const cell refused[] = {
        // LX == 1: horizontal self-neighbour.
        {1, 1, 0},
        {1, 2, 0},
        {1, 3, 0},
        // LY == 1 and skew = 0 mod LX: vertical self-neighbour.
        {2, 1, 0},
        {3, 1, 0},
        {4, 1, 0},
        // The same with a skew that is a non-zero multiple of LX, which
        // reaches the guard since contract section 8.1 (the lattice holds
        // skew 0, and the message names that).
        {2, 1, 2},
        {2, 1, -2},
        {3, 1, 3},
        {1, 1, 1},
        {1, 2, 5},
    };
    for (const cell c : refused) {
      INFO("L_sub = [" << c.lx << ", " << c.ly << "], skew = " << c.skew);
      std::string message;
      try {
        validate(c);
        FAIL_CHECK("a cell with a self-neighbour site was accepted");
        continue;
      } catch (const tenes::input_error &e) {
        message = e.what();
      }
      INFO("message: " << message);
      // Through the existing throw_fermion_guard wrapper.
      CHECK(message.find("fermion mode") != std::string::npos);
      // The cause: a site would be its own nearest neighbour.
      CHECK(contains_icase(message, "\\bown\\b"));
      CHECK(contains_icase(message, "neighbou?r"));
      // The shape: both L_sub values, and the skew the lattice holds.
      CHECK(names_numbers_after(message, "L_sub", {c.lx, c.ly}));
      CHECK(names_numbers_after(message, "skew", {c.skew % c.lx}));
    }
  }
}
