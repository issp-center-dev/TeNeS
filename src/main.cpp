#define _USE_MATH_DEFINES
#include <random>
#include <sys/stat.h>

#include <toml11/toml.hpp>

#include "Lattice.hpp"
#include "PEPS_Parameters.hpp"
#include "edge.hpp"
#include "hamiltonian.hpp"
#include "observable.hpp"
#include "tnsolve.hpp"


int main(int argc, char **argv) {
  using ptensor = mptensor::Tensor<mptensor::scalapack::Matrix, double>;
  int mpisize, mpirank;

  /* MPI initialization */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

  toml::table input_toml = toml::parse("input.toml");

  // Parameters
  PEPS_Parameters peps_parameters;

  Lattice lattice(2,2);

  if (mpirank == 0) {
    peps_parameters.set(input_toml);
  }

  peps_parameters.Bcast(MPI_COMM_WORLD);

  lattice.Bcast(MPI_COMM_WORLD);

  std::vector<ptensor> hams = load_hamiltonians(input_toml);
  std::vector<ptensor> lops = load_local_operators(input_toml);

  auto bonds_str = toml::find<std::string>(toml::find(input_toml, "bond"), "simple_update");
  auto fullbonds_str = toml::find<std::string>(toml::find(input_toml, "bond"), "full_update");
  Edges simple_edges = make_edges(bonds_str);
  Edges full_edges = make_edges(fullbonds_str);

  tnsolve(MPI_COMM_WORLD, peps_parameters, lattice, simple_edges, full_edges, hams, lops);

  MPI_Finalize();
  return 0;
}
