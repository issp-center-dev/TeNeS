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

#ifndef TENES_SRC_UTIL_TYPE_TRAITS_HPP_
#define TENES_SRC_UTIL_TYPE_TRAITS_HPP_

#include <complex>
#include <type_traits>

namespace tenes {

template <class to, class from>
to convert_complex(from const &v) {
  if constexpr (std::is_arithmetic_v<to>) {
    // convert to real
    return std::real(v);
  } else {
    // convert to complex
    return to(std::real(v), std::imag(v));
  }
}

}  // end of namespace tenes

#endif  // TENES_SRC_UTIL_TYPE_TRAITS_HPP_
