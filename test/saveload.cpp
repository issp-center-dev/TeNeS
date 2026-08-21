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

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"
#include "test_workdir.hpp"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#include "../src/tensor.hpp"
#include "../src/mpi.hpp"
#include "../src/exception.hpp"
#include "../src/iTPS/load_toml.hpp"
#include "../src/iTPS/iTPS.hpp"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  doctest::Context context(argc, argv);
  const int res = context.run();
  MPI_Finalize();
  return res;
}

namespace {

toml::value parse_str(std::string const &str) { return toml::parse_str(str); }

}  // namespace

TEST_CASE("save and load tensors") {
  using namespace tenes;
  using namespace tenes::itps;
  using ptensor = real_tensor;

  const std::string dir = "saveload_test_tensors";

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
  SquareLattice lattice = gen_lattice(toml.at("tensor"));

  PEPS_Parameters save_params;
  save_params.print_level = PrintLevel::none;
  save_params.outdir = "output_saveload_test";
  save_params.tensor_save_dir = dir;

  // the loading iTPS reads the tensors in its constructor
  PEPS_Parameters load_params = save_params;
  load_params.tensor_save_dir = "";
  load_params.tensor_load_dir = dir;

  auto make_itps = [&lattice](PEPS_Parameters const &params) {
    return iTPS<ptensor>(MPI_COMM_WORLD, params, lattice,
                         EvolutionOperators<ptensor>{},
                         EvolutionOperators<ptensor>{}, Operators<ptensor>{},
                         Operators<ptensor>{}, Operators<ptensor>{},
                         CorrelationParameter{}, TransferMatrix_Parameters{});
  };

  make_itps(save_params).save_tensors();

  SUBCASE("saved tensors can be loaded back") {
    CHECK_NOTHROW(make_itps(load_params));
  }

  SUBCASE("missing tensor file is an error") {
    REQUIRE(std::remove((dir + "/C1_0.dat").c_str()) == 0);
    CHECK_THROWS_AS(make_itps(load_params), tenes::load_error);
  }

  SUBCASE("missing lambda file is an error") {
    REQUIRE(std::remove((dir + "/lambda_2.dat").c_str()) == 0);
    CHECK_THROWS_AS(make_itps(load_params), tenes::load_error);
  }

  SUBCASE("truncated lambda file is an error") {
    {
      std::ofstream ofs(dir + "/lambda_1.dat");
      ofs << "1.0" << std::endl;
    }
    CHECK_THROWS_AS(make_itps(load_params), tenes::load_error);
  }
}
