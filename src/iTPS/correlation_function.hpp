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

/*! @file
 *  @brief Result and parameter types of the correlation-function
 *         measurement @f$\langle A(0) B(r) \rangle@f$ along the axes.
 */

#ifndef TENES_SRC_CORRELATION_HPP_
#define TENES_SRC_CORRELATION_HPP_

#include <tuple>
#include <vector>

namespace tenes::itps {

//! One measured correlation value
//! @f$\langle A(x_0, y_0)\, B(x_0 + dx, y_0 + dy) \rangle@f$.
struct Correlation {
  int left_index;  //!< site of the left operator
  int right_dx;    //!< x displacement of the right operator
  int right_dy;    //!< y displacement of the right operator
  int left_op;     //!< one-site operator index of A
  int right_op;    //!< one-site operator index of B
  double real;     //!< real part of the measured value
  double imag;     //!< imaginary part of the measured value
};

//! Settings of the correlation-function measurement, read from the
//! [correlation] table of input.toml.
struct CorrelationParameter {
  //! Maximum displacement along an axis (r_max); 0 disables the
  //! measurement.
  int r_max;
  //! Pairs of one-site operator indices to correlate (operators).
  std::vector<std::tuple<int, int>> operators;
  //! Disabled measurement (r_max = 0).
  CorrelationParameter() : r_max(0) {}
  //! Measure the given operator pairs up to distance r_max.
  CorrelationParameter(int r_max, std::vector<std::tuple<int, int>> const& ops)
      : r_max(r_max), operators(ops) {}
};

}  // namespace tenes::itps

#endif  // TENES_SRC_CORRELATION_HPP_
