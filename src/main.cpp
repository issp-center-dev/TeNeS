#define _USE_MATH_DEFINES
#include <random>
#include <sys/stat.h>

#include <cpptoml.h>

#include "Lattice.hpp"
#include "PEPS_Parameters.hpp"
#include "edge.hpp"
#include "util/read_matrix.hpp"
#include "tnsolve.hpp"


int main_impl(int argc, char **argv) {
  using ptensor = mptensor::Tensor<mptensor::scalapack::Matrix, double>;
  int mpisize, mpirank;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

  auto input_toml = cpptoml::parse_file("input.toml");

  // Parameters
  PEPS_Parameters peps_parameters;
  if (mpirank == 0) {
    peps_parameters.set(input_toml->get_table("parameter"));
  }
  peps_parameters.Bcast(MPI_COMM_WORLD);

  auto lattice_toml = input_toml->get_table("lattice");
  auto Lsub = lattice_toml->get_array_of<int64_t>("Lsub");
  if(!Lsub){
    std::cerr << "cannot find Lsub in the section [lattice]" << std::endl;
    return 1;
  }
  Lattice lattice((*Lsub)[0], (*Lsub)[1]);
  lattice.Bcast(MPI_COMM_WORLD);

  // time evolution
  auto evolution_toml = input_toml->get_table("evolution");
  auto str = evolution_toml->get_as<std::string>("simple_update");
  if(!str){
    std::cerr << "cannot find simple_update in the section [evolution]" << std::endl;
    return 1;
  }
  const auto simple_edges = make_edges(*str);

  str = evolution_toml->get_as<std::string>("full_update");
  if(!str){
    std::cerr << "cannot find full_update in the section [evolution]" << std::endl;
    return 1;
  }
  const auto full_edges = make_edges(*str);

  auto strs = evolution_toml->get_array_of<std::string>("matrix");
  if(!strs){
    std::cerr << "cannot find matrix in the section [evolution]" << std::endl;
    return 1;
  }
  const std::vector<ptensor> evolutions = util::read_matrix<ptensor>(*strs);

  // observable
  auto observable_toml = input_toml->get_table("observable");
  strs = observable_toml->get_array_of<std::string>("local_operator");
  if(!strs){
    std::cerr << "cannot find local_operator in the section [observable]" << std::endl;
    return 1;
  }
  const std::vector<ptensor> lops = util::read_matrix<ptensor>(*strs);

  strs = observable_toml->get_array_of<std::string>("hamiltonian");
  if(!strs){
    std::cerr << "cannot find hamiltonian in the section [observable]" << std::endl;
    return 1;
  }
  const std::vector<ptensor> hams = util::read_matrix<ptensor>(*strs);

  str = observable_toml->get_as<std::string>("hamiltonian_bonds");
  if(!str){
    std::cerr << "cannot find hamiltonian_bonds in the section [obsevable]" << std::endl;
    return 1;
  }
  const auto ham_edges = make_edges(*str);

  return tnsolve(MPI_COMM_WORLD, peps_parameters, lattice, simple_edges, full_edges, ham_edges, evolutions, hams, lops);
}


int main(int argc, char **argv){
  MPI_Init(&argc, &argv);
  int status = main_impl(argc, argv);
  MPI_Finalize();
  return status;
}
