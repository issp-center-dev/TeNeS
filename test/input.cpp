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
#include "../src/mpi.hpp"
#include "../src/fermion/fops.hpp"
#include "../src/util/string.hpp"
#include "../src/arpack_solver.hpp"
#include "../src/iTPS/load_toml.hpp"
#include "../src/iTPS/iTPS.hpp"
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

  SUBCASE("fermion rejects mean-field environment") {
    INFO("fermion rejects mean-field environment");
    auto param_toml = parse_str(R"(
[parameter]
[parameter.general]
fermion = true
[parameter.ctm]
meanfield_env = true
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
    CHECK_THROWS_AS(
        validate_fermion_constraints(
            peps_parameters, lattice, EvolutionOperators<ptensor>{},
            EvolutionOperators<ptensor>{}, Operators<ptensor>{},
            Operators<ptensor>{}, Operators<ptensor>{}, CorrelationParameter{}),
        tenes::input_error);
  }

  SUBCASE("fermion rejects odd one-site operator") {
    INFO("fermion rejects odd one-site operator");
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
    auto onesite = load_operators<ptensor>(observable_toml, MPI_COMM_WORLD, 1,
                                           1, 0.0, "observable.onesite");
    CHECK_THROWS_AS(
        validate_fermion_constraints(
            peps_parameters, lattice, EvolutionOperators<ptensor>{},
            EvolutionOperators<ptensor>{}, onesite, Operators<ptensor>{},
            Operators<ptensor>{}, CorrelationParameter{}),
        tenes::input_error);
  }
}
