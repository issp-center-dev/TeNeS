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

#include <mptensor/tensor.hpp>

namespace tenes {

namespace traits {

struct get_value_type {
  template <class T>
  static typename T::value_type check(typename T::value_type *);
  template <class T>
  static T check(...);
};

template <class to, class from>
to convert_complex(
    from const &v,
    typename std::enable_if<std::is_arithmetic<to>::value>::type *) {
  // convert to real
  return std::real(v);
}

template <class to, class from>
to convert_complex(
    from const &v,
    typename std::enable_if<!std::is_floating_point<to>::value>::type *) {
  // convert to complex
  return to(std::real(v), std::imag(v));
}

}  // end of namespace traits

template <class T>
struct get_value_type {
  using value_type = decltype(traits::get_value_type::check<T>(nullptr));
};

template <class to, class from>
to convert_complex(from const &v) {
  return traits::convert_complex<to>(v, nullptr);
}

}  // end of namespace tenes

#endif  // TENES_SRC_UTIL_TYPE_TRAITS_HPP_
