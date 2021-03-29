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

#ifndef TENES_SRC_ITPS_CORE_FULL_UPDATE_HPP_
#define TENES_SRC_ITPS_CORE_FULL_UPDATE_HPP_

namespace tenes {

class PEPS_Parameters;

template <class tensor>
void Full_update_bond_horizontal(
    const tensor &C1, const tensor &C2,
    const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2,
    const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6,
    const tensor &Tn1, const tensor &Tn2,
    const tensor &op12, const PEPS_Parameters peps_parameters,
    tensor &Tn1_new, tensor &Tn2_new);

template <class tensor>
void Full_update_bond(
    const tensor &C1, const tensor &C2,
    const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2,
    const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6,
    const tensor &Tn1, const tensor &Tn2,
    const tensor &op12, const int connect1,
    const PEPS_Parameters peps_parameters, tensor &Tn1_new,
    tensor &Tn2_new);

}  // end of namespace tenes

#endif  // TENES_SRC_ITPS_CORE_FULL_UPDATE_HPP_
