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

/*
 *
 Basic routines independent on unit cell structures.
 Using mptensor libraries
 (Test version)
 2015 Dec.  Tsuyoshi Okubo
*/

#ifndef TENES_SRC_ITPS_CORE_SIMPLE_UPDATE_HPP_
#define TENES_SRC_ITPS_CORE_SIMPLE_UPDATE_HPP_

#include <vector>

namespace tenes {
namespace itps {
class PEPS_Parameters;

namespace core {

template <class tensor>
void Simple_update_bond(const tensor &Tn1, const tensor &Tn2,
                        const std::vector<std::vector<double>> &lambda1,
                        const std::vector<std::vector<double>> &lambda2,
                        const tensor &op12, const int connect1,
                        const PEPS_Parameters peps_parameters, tensor &Tn1_new,
                        tensor &Tn2_new, std::vector<double> &lambda_c);
}  // end of namespace core
}  // namespace itps
}  // namespace tenes

#endif  // TENES_SRC_ITPS_CORE_SIMPLE_UPDATE_HPP_
