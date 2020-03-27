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

#ifndef TENES_TIMER_HPP
#define TENES_TIMER_HPP

#include <chrono>

namespace tenes {

template <class CT = std::chrono::high_resolution_clock>
class Timer {
 public:
  using clock_type = CT;

  Timer() : start(clock_type::now()) {}
  void reset() { start = clock_type::now(); }
  double elapsed() const {
    const auto elapsed_time =
        std::chrono::duration_cast<std::chrono::nanoseconds>(clock_type::now() - start);
    return elapsed_time.count() * 1.0e-9;
  }

 private:
  typename clock_type::time_point start;
};

}  // end of namespace tenes

#endif  // TENES_TIMER_HPP
