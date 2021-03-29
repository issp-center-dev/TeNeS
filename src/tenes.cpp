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

#include "tenes.hpp"
#include "iTPS/iTPS.hpp"

namespace tenes {

template <class tensor>
int tenes(MPI_Comm comm, PEPS_Parameters peps_parameters, Lattice lattice,
          NNOperators<tensor> simple_updates, NNOperators<tensor> full_updates,
          Operators<tensor> onesite_operators,
          Operators<tensor> twosite_operators, CorrelationParameter corparam,
          CorrelationLengthCalculator_Parameters clength_param) {
  iTPS<tensor> tns(comm, peps_parameters, lattice, simple_updates, full_updates,
                   onesite_operators, twosite_operators, corparam,
                   clength_param);
  tns.optimize();
  tns.save_tensors();
  if (peps_parameters.to_measure) {
    tns.measure();
  }
  tns.summary();
  return 0;
}

// template specialization
template int tenes<real_tensor>(
    MPI_Comm comm, PEPS_Parameters peps_parameters, Lattice lattice,
    NNOperators<real_tensor> simple_updates,
    NNOperators<real_tensor> full_updates,
    Operators<real_tensor> onesite_operators,
    Operators<real_tensor> twosite_operators, CorrelationParameter corparam,
    CorrelationLengthCalculator_Parameters clength_param);

template int tenes<complex_tensor>(
    MPI_Comm comm, PEPS_Parameters peps_parameters, Lattice lattice,
    NNOperators<complex_tensor> simple_updates,
    NNOperators<complex_tensor> full_updates,
    Operators<complex_tensor> onesite_operators,
    Operators<complex_tensor> twosite_operators, CorrelationParameter corparam,
    CorrelationLengthCalculator_Parameters clength_param);

}  // end of namespace tenes
