#ifndef TNSOLVE_HPP
#define TNSOLVE_HPP

#include <vector>
#include <mpi.h>

class Lattice;
class PEPS_Parameters;
struct Edge;
using Edges = std::vector<Edge>;

template<class tensor>
int tnsolve(MPI_Comm comm,
            PEPS_Parameters peps_parameters,
            Lattice lattice,
            Edges simple_edges,
            Edges full_edges,
            Edges ham_edges,
            std::vector<tensor> timeevolutions,
            std::vector<tensor> hamiltonians,
            std::vector<tensor> local_operators
            );

#endif // TNSOLVE_HPP
