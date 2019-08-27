#ifndef TENES_HPP
#define TENES_HPP

#include <mpi.h>
#include <vector>

class Lattice;
class PEPS_Parameters;
struct Edge;
using Edges = std::vector<Edge>;

template <class tensor>
int tenes(MPI_Comm comm, PEPS_Parameters peps_parameters, Lattice lattice,
          Edges simple_edges, Edges full_edges, Edges ham_edges,
          std::vector<tensor> timeevolutions, std::vector<tensor> hamiltonians,
          std::vector<tensor> local_operators);

#endif // TENES_HPP
