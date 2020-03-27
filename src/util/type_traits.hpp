#ifndef UTIL_TYPE_TRAITS_HPP
#define UTIL_TYPE_TRAITS_HPP

#include <complex>
#include <type_traits>

#include <mptensor/tensor.hpp>

namespace tenes {

namespace traits {

struct get_value_type {
  template <class T>
  static typename T::value_type check(typename T::value_type *);
  template <class T> static T check(...);
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

} // end of namespace traits

template <class T> struct get_value_type {
  using value_type = decltype(traits::get_value_type::check<T>(nullptr));
};

template <class to, class from> to convert_complex(from const &v) {
  return traits::convert_complex<to>(v, nullptr);
}

} // end of namespace tenes

#endif // UTIL_TYPE_TRAITS_HPP
