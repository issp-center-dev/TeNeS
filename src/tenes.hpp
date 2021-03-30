/* TeNeS - Massively parallel tensor network solver /
/ Copyright (C) 2019- The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify /
/ it under the terms of the GNU General Public License as published by /
/ the Free Software Foundation, either version 3 of the License, or /
/ (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, /
/ but WITHOUT ANY WARRANTY; without even the implied warranty of /
/ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the /
/ GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License /
/ along with this program. If not, see http://www.gnu.org/licenses/. */

#ifndef TENES_SRC_TENES_HPP_
#define TENES_SRC_TENES_HPP_

#include <vector>

#include "mpi.hpp"

#include "operator.hpp"

namespace tenes {

class SquareLattice;
class PEPS_Parameters;
struct Edge;
using Edges = std::vector<Edge>;
struct CorrelationParameter;
struct CorrelationLengthCalculator_Parameters;

template <class tensor>
int tenes(MPI_Comm comm, PEPS_Parameters peps_parameters, SquareLattice lattice,
          NNOperators<tensor> simple_updates, NNOperators<tensor> full_updates,
          Operators<tensor> onesite_operators,
          Operators<tensor> twosite_operators, CorrelationParameter corparam,
          CorrelationLengthCalculator_Parameters clength_param);

}  // end of namespace tenes

#endif  // TENES_SRC_TENES_HPP_
