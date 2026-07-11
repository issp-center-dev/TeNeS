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

#ifndef TENES_SRC_ITPS_CORRELATION_LENGTH_HPP_
#define TENES_SRC_ITPS_CORRELATION_LENGTH_HPP_

namespace tenes::itps {

/*! @brief Computes the correlation length from transfer matrix eigenvalues.
 *
 *  @param[in] e0_abs absolute value of the largest eigenvalue
 *  @param[in] e1_abs absolute value of the second largest eigenvalue
 *  @param[in] L      length of the unit cell along the transfer direction
 *
 *  @return xi = -L / log(e1_abs / e0_abs).
 *          Returns +infinity when the two eigenvalues are degenerate
 *          (e1_abs >= e0_abs) and 0.0 when e0_abs is not positive.
 */
double calc_correlation_length(double e0_abs, double e1_abs, int L);

}  // namespace tenes::itps

#endif  // TENES_SRC_ITPS_CORRELATION_LENGTH_HPP_
