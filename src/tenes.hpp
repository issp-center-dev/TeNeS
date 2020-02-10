#ifndef TENES_HPP
#define TENES_HPP

#include <vector>

#include "mpi.hpp"

#include "operator.hpp"

namespace tenes {

class Lattice;
class PEPS_Parameters;
struct Edge;
using Edges = std::vector<Edge>;
struct CorrelationParameter;

template <class tensor>
int tenes(MPI_Comm comm, PEPS_Parameters peps_parameters, Lattice lattice,
          Operators<tensor> simple_updates, Operators<tensor> full_updates,
          Edges ham_edges, std::vector<tensor> hamiltonians,
          std::vector<tensor> local_operators, CorrelationParameter corparam);

} // end of namespace tenes

#endif // TENES_HPP
