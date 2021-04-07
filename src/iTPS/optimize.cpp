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

#include "iTPS.hpp"
#include "../printlevel.hpp"
#include "../tensor.hpp"

namespace tenes {

template <class ptensor>
void iTPS<ptensor>::optimize() {
  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "Start simple update" << std::endl;
  }
  simple_update();

  if (peps_parameters.num_full_step > 0) {
    if (peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "Start full update" << std::endl;
    }
    full_update();
  }
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // end of namespace tenes
