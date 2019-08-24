#define _USE_MATH_DEFINES
#include <random>
#include <sys/stat.h>

#include <toml11/toml.hpp>

#include "Lattice.hpp"
#include "PEPS_Parameters.hpp"
#include "edge.hpp"
#include "util/read_matrix.hpp"
#include "tnsolve.hpp"


int main(int argc, char **argv) {
  using toml::find;
  using ptensor = mptensor::Tensor<mptensor::scalapack::Matrix, double>;
  int mpisize, mpirank;

  /* MPI initialization */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

  toml::value input_toml = toml::parse("input.toml");

  // Parameters
  PEPS_Parameters peps_parameters;

  std::vector<int> Lsub = toml::find<std::vector<int>>(toml::find(input_toml, "lattice"), "Lsub");

  Lattice lattice(Lsub[0], Lsub[1]);

  if (mpirank == 0) {
    peps_parameters.set(input_toml);
  }

  peps_parameters.Bcast(MPI_COMM_WORLD);

  lattice.Bcast(MPI_COMM_WORLD);

  // time evolution
  auto str = find<std::string>(find(input_toml, "evolution"), "simple_update");
  const auto simple_edges = make_edges(str);

  str = find<std::string>(find(input_toml, "evolution"), "full_update");
  const auto full_edges = make_edges(str);

  auto strs = find<std::vector<std::string>>(find(input_toml, "evolution"), "matrix");
  const std::vector<ptensor> evolutions = util::read_matrix<ptensor>(strs);

  // observables
  strs = find<std::vector<std::string>>(find(input_toml, "observable"), "local_operator");
  const std::vector<ptensor> lops = util::read_matrix<ptensor>(strs);

  strs = find<std::vector<std::string>>(find(input_toml, "observable"), "hamiltonian");
  const std::vector<ptensor> hams = util::read_matrix<ptensor>(strs);

  str = find<std::string>(find(input_toml, "observable"), "hamiltonian_bonds");
  const auto ham_edges = make_edges(str);


  tnsolve(MPI_COMM_WORLD, peps_parameters, lattice, simple_edges, full_edges, ham_edges, evolutions, hams, lops);

  MPI_Finalize();
  return 0;
}
