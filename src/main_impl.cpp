#include <complex>

#include <cpptoml.h>

#include "Lattice.hpp"
#include "PEPS_Parameters.hpp"
#include "load_toml.cpp"
#include "operator.hpp"
#include "tenes.hpp"
#include "type.hpp"
#include "util/file.hpp"
#include "util/read_matrix.hpp"
#include "version.hpp"

namespace tenes {

int main_impl(int argc, char **argv) {
  // using ptensor = mptensor::Tensor<mptensor_matrix_type, std::complex<double>>;
  using ptensor = mptensor::Tensor<mptensor_matrix_type, double>;

  int mpisize = 0, mpirank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

  std::string usage =
      R"(TeNeS: TEnsor NEtwork Solver for 2D quantum lattice system
  
  Usage:
    tenes [--quiet] <input_toml>
    tenes --help
    tenes --version

  Options:
    -h --help       Show this help message.
    -v --version    Show the version.
    -q --quiet      Do not print any messages.
  )";

  if (argc == 1) {
    if (mpirank == 0)
      std::cout << usage << std::endl;
    return 0;
  }

  for (int i = 1; i < argc; ++i) {
    std::string opt = argv[i];
    if (opt == "-h" || opt == "--help") {
      if (mpirank == 0)
        std::cout << usage << std::endl;
      return 0;
    }
  }

  for (int i = 1; i < argc; ++i) {
    std::string opt = argv[i];
    if (opt == "-v" || opt == "--version") {
      if (mpirank == 0)
        std::cout << "TeNeS v" << TENES_VERSION << std::endl;
      return 0;
    }
  }

  using PrintLevel = PEPS_Parameters::PrintLevel;

  PrintLevel print_level = PrintLevel::info;
  std::string input_filename;
  for (int i = 1; i < argc; ++i) {
    std::string opt = argv[i];
    if (opt == "-q" || opt == "--quiet") {
      print_level = PrintLevel::none;
    } else {
      input_filename = opt;
    }
  }

  if (!file_exists(input_filename)) {
    if (mpirank == 0)
      std::cout << "ERROR: cannot find the input file: " << input_filename
                << std::endl;
    return 1;
  }

  auto input_toml = cpptoml::parse_file(input_filename);

  // Parameters
  auto toml_param = input_toml->get_table("parameter");
  PEPS_Parameters peps_parameters =
      (toml_param != nullptr ? gen_param(toml_param) : PEPS_Parameters());
  peps_parameters.print_level = print_level;
  peps_parameters.Bcast(MPI_COMM_WORLD);

  auto toml_lattice = input_toml->get_table("tensor");
  if (toml_lattice == nullptr) {
    // ERROR
    std::cout << "[tensor] not found" << std::endl;
    return 1;
  }
  Lattice lattice = gen_lattice(toml_lattice);
  lattice.Bcast(MPI_COMM_WORLD);

  // time evolution
  auto toml_evolution = input_toml->get_table("evolution");
  if (toml_evolution == nullptr) {
    // ERROR
    std::cout << "[evolution] not found" << std::endl;
    return 1;
  }

  const auto simple_updates = load_simple_updates<ptensor>(input_toml);
  const auto full_updates = load_full_updates<ptensor>(input_toml);

  // observable
  auto toml_observable = input_toml->get_table("observable");
  if (toml_observable == nullptr) {
    // ERROR
    std::cout << "[observable] not found" << std::endl;
    return 1;
  }

  // onesite observable
  const auto onesite_obs = load_operators<ptensor>(input_toml, lattice.N_UNIT, 1,
                                                  "observable.onesite");
  const auto twosite_obs = load_operators<ptensor>(input_toml, lattice.N_UNIT,
                                                   2, "observable.twosite");

  // correlation
  auto toml_correlation = input_toml->get_table("correlation");
  const auto corparam = (toml_correlation != nullptr
                             ? gen_corparam(toml_correlation, "correlation")
                             : CorrelationParameter());

  return tenes(MPI_COMM_WORLD, peps_parameters, lattice, simple_updates,
               full_updates, onesite_obs, twosite_obs,
               corparam);
}

} // end of namespace tenes
