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
    // itps_main does, in main.cpp), so the new C2 tensor.L_sub-dimensions
    // guard is never on this subcase's path and does not need a 2x2 cell.
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
    // never calls validate_fermion_constraints, so the new C2
    // tensor.L_sub-dimensions guard never runs here and the 1x1 cell still
    // loads cleanly.
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
    // L_sub = [2, 2]: with a 1-wide cell the tensor.L_sub-dimensions guard
    // (C2, task-11-contract.md) would throw before MeanField_Env is even
    // looked at, so the subcase could not tell that the mean-field
    // environment itself is accepted. A 2x2 cell with a single site
    // definition broadcast via index = [] clears that guard; with the
    // fermionic mean-field measurement in place, MeanField_Env=true is a
    // supported combination and nothing else in this input is guarded.
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
    // environment" above -- a 1-wide cell would trip the new
    // tensor.L_sub-dimensions guard before the parity-odd one-site
    // operator check this subcase is named for.
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
      // L_sub = [2, 2]: a 1-wide cell now trips the new
      // tensor.L_sub-dimensions guard (C2, task-11-contract.md) before
      // reaching the parity-odd one-site operator check this subcase
      // exercises; a single site definition broadcast via index = []
      // keeps this a minimal fixture while clearing that guard.
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

//! Redirect std::cerr into a string for the lifetime of the object.
//! The fermionic fast-full-update fallback has to announce itself on the
//! standard error stream (contract 3.5, T6), and a doctest cannot read the
//! process' stderr, so the stream buffer is swapped instead.
class cerr_capture {
 public:
  cerr_capture() { saved_ = std::cerr.rdbuf(buffer_.rdbuf()); }
  ~cerr_capture() { std::cerr.rdbuf(saved_); }
  cerr_capture(cerr_capture const &) = delete;
  cerr_capture &operator=(cerr_capture const &) = delete;
  std::string str() const { return buffer_.str(); }

 private:
  std::ostringstream buffer_;
  std::streambuf *saved_ = nullptr;
};

//! (source_site, source_leg) of eight gates that cover the eight bonds of a
//! 2x2 unit cell exactly once, using every source_leg value.
//!
//! src/SquareLattice.cpp maps leg 0 to (x-1, y), leg 1 to (x, y+1), leg 2 to
//! (x+1, y) and leg 3 to (x, y-1), with site = x + 2 * y. The raster-earlier
//! site of a bond is its left (horizontal) or upper (vertical) end, so legs 2
//! and 3 name a bond from the earlier end and legs 0 and 1 from the later
//! one; the latter are the cases that force the driver to swap the two sites.
constexpr int kFullUpdateBonds[8][2] = {{0, 2}, {2, 2}, {0, 0}, {2, 0},
                                        {0, 1}, {2, 1}, {1, 3}, {3, 3}};

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

TEST_CASE("fermion mode falls back from the fast full update with a warning") {
  using namespace tenes;
  using namespace tenes::itps;

  // Precondition: fastfullupdate really is on by default, so the input
  // below exercises the fallback rather than the plain path.
  {
    auto defaults = gen_param(parse_str(R"([parameter])").at("parameter"));
    REQUIRE(defaults.Full_Use_FastFullUpdate == true);
  }

  const std::string input_filename =
      "test_input_fermion_fast_full_update.toml";
  const std::string outdir = "output_test_input_fermion_fast_full_update";

  // A parity-even hopping gate; the diagonal keeps the state from collapsing
  // and the off-diagonal block makes the bond genuinely correlated.
  const std::string gate_elements =
      "0 0 0 0  1.0 0.0\n"
      "0 1 0 1  1.0 0.0\n"
      "1 0 1 0  1.0 0.0\n"
      "1 1 1 1  1.0 0.0\n"
      "0 1 1 0  0.1 0.0\n"
      "1 0 0 1  0.1 0.0";

  {
    std::ostringstream ofs;
    ofs << R"(
[parameter]
[parameter.general]
is_real = true
fermion = true
output = ")"
        << outdir << R"("

[parameter.simple_update]
tau = 0.01
num_step = 20

[parameter.full_update]
tau = 0.01
num_step = 1

[parameter.ctm]
dimension = 8
convergence_epsilon = 1.0e-8
iteration_max = 20

[parameter.random]
seed = 11
)" << fermion_cell_toml()
        << R"(
[observable]
[[observable.onesite]]
name = "n"
group = 0
sites = []
dim = 2
elements = """
1 1  1.0 0.0
"""

[[observable.twosite]]
name = "bond_hamiltonian"
group = 0
bonds = """
0 1 0
1 1 0
2 1 0
3 1 0
0 0 1
1 0 1
2 0 1
3 0 1
"""
dim = [2, 2]
elements = """
0 1 1 0  -1.0 0.0
1 0 0 1  -1.0 0.0
"""

[evolution]
)";
    for (int leg : {2, 1}) {
      for (int site = 0; site < 4; ++site) {
        for (const char *kind : {"simple", "full"}) {
          ofs << "[[evolution." << kind << "]]\n"
              << "group = 0\n"
              << "source_site = " << site << "\n"
              << "source_leg = " << leg << "\n"
              << "dimensions = [2, 2, 2, 2]\n"
              << "elements = \"\"\"\n"
              << gate_elements << "\n\"\"\"\n";
        }
      }
    }
    std::ofstream out(input_filename.c_str());
    out << ofs.str();
  }

  std::string captured;
  {
    cerr_capture capture;
    CHECK_NOTHROW(
        tenes::itps::itps_main(input_filename, MPI_COMM_WORLD, PrintLevel::none));
    captured = capture.str();
  }
  INFO("captured stderr: " << captured);
  CHECK(captured.find("Full_Use_FastFullUpdate") != std::string::npos);

  std::remove(input_filename.c_str());
  std::error_code ec;
  std::filesystem::remove_all(outdir, ec);
}

TEST_CASE("a fermionic full update leaves the parity ledger consistent") {
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
[parameter.simple_update]
tau = 0.01
num_step = 20
[parameter.full_update]
tau = 0.01
num_step = 1
fastfullupdate = false
[parameter.ctm]
dimension = 8
convergence_epsilon = 1.0e-8
iteration_max = 20
[parameter.random]
seed = 11
)");
  PEPS_Parameters peps_parameters = gen_param(param_toml.at("parameter"));
  peps_parameters.phys_parity =
      gen_phys_parity(tensor_toml.at("tensor"), lattice);
  peps_parameters.print_level = PrintLevel::none;
  peps_parameters.outdir = "output_test_input_fermion_full_ledger";

  auto gate = [&]() {
    ptensor op(comm, mptensor::Shape(2, 2, 2, 2));
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        op.set_value(mptensor::Index(i, j, i, j), 1.0);
      }
    }
    op.set_value(mptensor::Index(0, 1, 1, 0), 0.1);
    op.set_value(mptensor::Index(1, 0, 0, 1), 0.1);
    return op;
  }();

  EvolutionOperators<ptensor> simple_updates;
  for (int leg : {2, 1}) {
    for (int site = 0; site < 4; ++site) {
      simple_updates.push_back(
          make_twosite_EvolutionOperator<ptensor>(site, leg, 0, gate));
    }
  }
  EvolutionOperators<ptensor> full_updates;
  for (const auto &bond : kFullUpdateBonds) {
    full_updates.push_back(
        make_twosite_EvolutionOperator<ptensor>(bond[0], bond[1], 0, gate));
  }

  // Preconditions. The bonds must be gated from the raster-later end at
  // least once (source_leg 0 and 1), because that is the path whose ledger
  // bookkeeping this test is about; and the eight gates must name eight
  // different bonds.
  auto bond_name = [&](int site, int leg) {
    const int x = site % lattice.LX;
    const int y = site / lattice.LX;
    if (leg == 2) return std::make_tuple(0, x, y);
    if (leg == 0) return std::make_tuple(0, (x - 1 + lattice.LX) % lattice.LX, y);
    if (leg == 1) return std::make_tuple(1, x, y);
    return std::make_tuple(1, x, (y - 1 + lattice.LY) % lattice.LY);
  };
  std::set<std::tuple<int, int, int>> named;
  bool has_leg0 = false;
  bool has_leg1 = false;
  for (const auto &up : full_updates) {
    named.insert(bond_name(up.source_site, up.source_leg));
    has_leg0 = has_leg0 || up.source_leg == 0;
    has_leg1 = has_leg1 || up.source_leg == 1;
  }
  REQUIRE(has_leg0);
  REQUIRE(has_leg1);
  REQUIRE(named.size() == 8);
  REQUIRE(peps_parameters.num_full_step[0] == 1);
  REQUIRE(peps_parameters.Full_Use_FastFullUpdate == false);

  iTPS<ptensor> state(comm, peps_parameters, lattice, simple_updates,
                      full_updates, Operators<ptensor>{}, Operators<ptensor>{},
                      Operators<ptensor>{}, CorrelationParameter{},
                      TransferMatrix_Parameters{});
  REQUIRE(iTPSTestAccessor::finfo(state).enabled);

  CHECK_NOTHROW(state.optimize());

  // The full update must actually have moved the state: a no-op would keep
  // every ledger trivially consistent and make the checks below hollow. The
  // baseline runs the identical simple update and stops there.
  {
    PEPS_Parameters su_only = peps_parameters;
    su_only.num_full_step.assign(su_only.num_full_step.size(), 0);
    iTPS<ptensor> baseline(comm, su_only, lattice, simple_updates,
                           full_updates, Operators<ptensor>{},
                           Operators<ptensor>{}, Operators<ptensor>{},
                           CorrelationParameter{},
                           TransferMatrix_Parameters{});
    baseline.optimize();
    const auto &after = iTPSTestAccessor::Tn(state);
    const auto &before = iTPSTestAccessor::Tn(baseline);
    REQUIRE(after.size() == before.size());
    double moved = 0.0;
    for (std::size_t site = 0; site < after.size(); ++site) {
      moved = std::max(moved, mptensor::max_abs(after[site] - before[site]));
    }
    INFO("largest change of Tn caused by the full update: " << moved);
    REQUIRE(moved > 1e-8);
  }

  const auto &fi = iTPSTestAccessor::finfo(state);
  CHECK_NOTHROW(tenes::fermion::validate_neighbor_consistency(fi, lattice));

  // The per-bond contract (T3-iii) asks for 1e-12 * max_abs on a single
  // bond update; this runs twenty simple-update steps and a full-update
  // sweep before looking, so the bound is relaxed by two decades. It is
  // still far tighter than anything a leaked odd component would satisfy.
  const auto &tensors = iTPSTestAccessor::Tn(state);
  REQUIRE(tensors.size() == static_cast<std::size_t>(lattice.N_UNIT));
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    auto ft = tenes::fermion::wrap_Tn(tensors[site], fi, site);
    const double scale = tenes::fermion::max_abs(ft);
    REQUIRE(scale > 0.0);
    CHECK(tenes::fermion::parity_violation(ft) <= 1e-10 * scale);
  }

  std::error_code ec;
  std::filesystem::remove_all(peps_parameters.outdir, ec);
}
