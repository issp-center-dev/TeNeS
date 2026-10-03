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

#ifndef TENES_SRC_ITPS_MEASURE_WINDOW_HPP_
#define TENES_SRC_ITPS_MEASURE_WINDOW_HPP_

#include <algorithm>
#include <cstdlib>
#include <sstream>

#include "../exception.hpp"
#include "../operator.hpp"

namespace tenes::itps {

//! Number of sites along an edge of the window the measurement contracts.
constexpr int measure_window_size = 4;

/*!
 * @brief Reject observables that do not fit into the measurement window.
 *
 * measure_twosite() and measure_multisite() contract a window of at most
 * 4 x 4 sites. An observable that does not fit cannot be measured; it used
 * to be left out of the result with a warning. Throws tenes::input_error
 * for the first one found.
 *
 * Called when the input is read, so that a calculation does not start, and
 * again by the measurement, for a state that was not made from an input file.
 */
template <class tensor>
void validate_measure_window(const Operators<tensor> &twosite_operators,
                             const Operators<tensor> &multisite_operators) {
  constexpr int nmax = measure_window_size;
  for (auto const &op : twosite_operators) {
    if (op.dx.empty() ||
        (std::abs(op.dx[0]) < nmax && std::abs(op.dy[0]) < nmax)) {
      continue;
    }
    std::stringstream ss;
    ss << "ERROR: twosite observable \"" << op.name << "\" (group "
       << op.group << ") has the bond source_site = " << op.source_site
       << ", dx = " << op.dx[0] << ", dy = " << op.dy[0]
       << ", which is outside the " << nmax << "x" << nmax
       << " measurement window and cannot be measured.\n"
       << "       |dx| and |dy| must not exceed " << nmax - 1 << ".";
    throw tenes::input_error(ss.str());
  }
  for (auto const &op : multisite_operators) {
    int mindx = 0, maxdx = 0, mindy = 0, maxdy = 0;
    for (auto dx : op.dx) {
      mindx = std::min(mindx, dx);
      maxdx = std::max(maxdx, dx);
    }
    for (auto dy : op.dy) {
      mindy = std::min(mindy, dy);
      maxdy = std::max(maxdy, dy);
    }
    if (maxdx - mindx < nmax && maxdy - mindy < nmax) {
      continue;
    }
    std::stringstream ss;
    ss << "ERROR: multisite observable \"" << op.name << "\" (group "
       << op.group << ") has sites around source_site = " << op.source_site
       << " that do not fit into the " << nmax << "x" << nmax
       << " measurement window, and cannot be measured.";
    throw tenes::input_error(ss.str());
  }
}

}  // namespace tenes::itps

#endif  // TENES_SRC_ITPS_MEASURE_WINDOW_HPP_
