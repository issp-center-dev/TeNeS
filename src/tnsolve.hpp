#ifndef TNSOLVE_HPP
#define TNSOLVE_HPP

#include <vector>
#include <mpi.h>

class Lattice;
class PEPS_Parameters;
class Edge;
using Edges = std::vector<Edge>;

template<class tensor>
int tnsolve(MPI_Comm comm,
            PEPS_Parameters peps_parameters,
            Lattice lattice,
            Edges simple_edges,
            Edges full_edges,
            std::vector<tensor> hams,
            std::vector<tensor> lops
            );

#endif // TNSOLVE_HPP
