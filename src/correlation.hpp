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

#ifndef CORRELATION_HPP
#define CORRELATION_HPP

#include <tuple>
#include <vector>

namespace tenes {

struct Correlation {
  int left_index;
  int right_dx, right_dy;
  int left_op, right_op;
  double real, imag;
};

struct CorrelationParameter {
  int r_max;
  std::vector<std::tuple<int, int>> operators;
  CorrelationParameter() : r_max(0) {}
  CorrelationParameter(int r_max, std::vector<std::tuple<int, int>> const& ops)
      : r_max(r_max), operators(ops) {}
};

}  // end of namespace tenes

#endif  // CORRELATION_HPP
