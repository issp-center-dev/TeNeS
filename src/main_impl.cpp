#include <cpptoml.h>
#include "type.hpp"
#include "util/file.hpp"
#include "Lattice.hpp"
#include "PEPS_Parameters.hpp"
#include "edge.hpp"
#include "load_toml.cpp"
#include "tenes.hpp"
#include "util/read_matrix.hpp"
#include "version.hpp"

namespace tenes {

int main_impl(int argc, char **argv) {
  using ptensor = mptensor::Tensor<mptensor_matrix_type, double>;

  int mpisize=0, mpirank=0;
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

  if(argc==1){
    if(mpirank==0)
      std::cout << usage << std::endl;
    return 0;
  }

  for(int i=1; i<argc; ++i){
    std::string opt = argv[i];
    if(opt == "-h" || opt == "--help"){
      if(mpirank==0)
        std::cout << usage << std::endl;
      return 0;
    }
  }

  for(int i=1; i<argc; ++i){
    std::string opt = argv[i];
    if(opt == "-v" || opt == "--version"){
      if(mpirank==0)
        std::cout << "TeNeS v" << TENES_VERSION << std::endl;
      return 0;
    }
  }

  using PrintLevel = PEPS_Parameters::PrintLevel;

  PrintLevel print_level = PrintLevel::info;
  std::string input_filename;
  for(int i=1; i<argc; ++i){
    std::string opt = argv[i];
    if(opt == "-q" || opt == "--quiet"){
      print_level = PrintLevel::none;
    }else{
      input_filename = opt;
    }
  }

  if(!file_exists(input_filename)){
    if(mpirank==0)
      std::cout << "ERROR: cannot find the input file: "<< input_filename << std::endl;
    return 1;
  }

  auto input_toml = cpptoml::parse_file(input_filename);

  // Parameters
  auto toml_param = input_toml->get_table("parameter");
  PEPS_Parameters peps_parameters =
      (toml_param != nullptr ? gen_param(toml_param) : PEPS_Parameters());
  peps_parameters.print_level = print_level;
  peps_parameters.Bcast(MPI_COMM_WORLD);

  auto toml_lattice = input_toml->get_table("lattice");
  if (toml_lattice == nullptr) {
    // ERROR
    std::cout << "[lattice] not found" << std::endl;
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
  const auto simple_edges =
      gen_edges(toml_evolution, "simple_update", "evolution");
  const auto full_edges = gen_edges(toml_evolution, "full_update", "evolution");
  const auto evolutions =
      gen_matrices<ptensor>(toml_evolution, "matrix", "evolution");

  // observable
  auto toml_observable = input_toml->get_table("observable");
  if (toml_observable == nullptr) {
    // ERROR
    std::cout << "[observable] not found" << std::endl;
    return 1;
  }
  const auto lops =
      gen_matrices<ptensor>(toml_observable, "local_operator", "observable");
  const auto hams =
      gen_matrices<ptensor>(toml_observable, "hamiltonian", "observable");
  const auto ham_edges =
      gen_edges(toml_observable, "hamiltonian_bonds", "observable");

  // correlation
  auto toml_correlation = input_toml->get_table("correlation");
  const auto corparam = (toml_correlation != nullptr
                             ? gen_corparam(toml_correlation, "correlation")
                             : CorrelationParameter());

  return tenes(MPI_COMM_WORLD, peps_parameters, lattice, simple_edges,
               full_edges, ham_edges, evolutions, hams, lops, corparam);
}

} // end of namespace tenes
