#define _USE_MATH_DEFINES
#include <random>
#include <sys/stat.h>

#include <mptensor/complex.hpp>
#include <mptensor/rsvd.hpp>
#include <mptensor/tensor.hpp>
#include <toml11/toml.hpp>

#include "Lattice.hpp"
#include "PEPS_Basics.hpp"
#include "PEPS_Parameters.hpp"
#include "Parameters.hpp"
#include "edge.hpp"
#include "hamiltonian.hpp"
#include "tnsolve.hpp"

using namespace mptensor;
using ptensor = Tensor<scalapack::Matrix, double>;

int main(int argc, char **argv) {
  int mpisize, mpirank;

  /* MPI initialization */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

  toml::table input_toml = toml::parse("input.toml");

  // Parameters
  PEPS_Parameters peps_parameters;
  Parameters local_parameters;

  Lattice lattice(2,2);

  if (mpirank == 0) {
    local_parameters.set(input_toml);
    peps_parameters.set(input_toml);
  }

  local_parameters.Bcast_parameters(MPI_COMM_WORLD);
  peps_parameters.Bcast_parameters(MPI_COMM_WORLD);

  lattice.Bcast_parameters(MPI_COMM_WORLD);

  const double hx = 1.4;

  // std::vector<ptensor> hams = load_hamiltonians(input_toml);
  std::vector<ptensor> hams = {Set_Hamiltonian(hx)};

  auto bonds_str = toml::find<std::string>(toml::find(input_toml, "bond"), "simple_update");
  auto fullbonds_str = toml::find<std::string>(toml::find(input_toml, "bond"), "full_update");
  Edges simple_edges = make_edges(bonds_str);
  Edges full_edges = make_edges(fullbonds_str);

  tnsolve(MPI_COMM_WORLD, peps_parameters, local_parameters, lattice, simple_edges, full_edges, hams);


  MPI_Finalize();
  return 0;
}
