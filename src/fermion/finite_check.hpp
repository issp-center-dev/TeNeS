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
 *  @brief Deciding whether a tensor element is finite.
 *
 *  Its own header so that test_fermion_finite_math can compile just this
 *  predicate with -ffast-math and check that it still works.
 */

#ifndef TENES_SRC_FERMION_FINITE_CHECK_HPP_
#define TENES_SRC_FERMION_FINITE_CHECK_HPP_

#include <complex>
#include <cstdint>
#include <cstring>

namespace tenes {
namespace fermion {
namespace detail {

/*! @brief True iff v is neither infinite nor NaN.
 *
 * Decided from the bit pattern rather than with std::isfinite, because
 * under finite-math-only (-ffast-math, or icpx's default -fp-model=fast)
 * the compiler may assume no operand is Inf or NaN and fold every
 * floating-point predicate to true. This one is the last line of defence
 * for a diverged state, so it must not be something the optimizer is
 * allowed to reason away. IEEE-754 binary64: an all-ones exponent field
 * is Inf (zero mantissa) or NaN (nonzero mantissa), and nothing else is.
 *
 * Reading the bits is not by itself enough. clang recognises the
 * mask-and-compare as a floating-point predicate and rewrites it back into
 * one -- the IR for the plain version is `fabs(v) != Inf` -- and under
 * finite-math-only the argument arrives with nofpclass(nan inf), which
 * folds that to a constant. Apple clang 21 at -O1 and above emits
 * `ret i1 true` for the whole function; GCC 16 does not. The volatile
 * stops the rewrite by making the value an opaque integer, and the store
 * and load it costs are on a path that only runs once a decomposition has
 * already failed (fops.hpp, note_failed_block).
 */
inline bool is_finite_value(double v) {
  static_assert(sizeof(double) == sizeof(std::uint64_t),
                "is_finite_value assumes IEEE-754 binary64");
  std::uint64_t bits = 0;
  std::memcpy(&bits, &v, sizeof(bits));
  volatile std::uint64_t opaque = bits;
  return (opaque & 0x7ff0000000000000ULL) != 0x7ff0000000000000ULL;
}

//! is_finite_value() for a complex element: both parts must be finite.
template <class T>
inline bool is_finite_value(const std::complex<T>& v) {
  return is_finite_value(static_cast<double>(v.real())) &&
         is_finite_value(static_cast<double>(v.imag()));
}

}  // namespace detail
}  // namespace fermion
}  // namespace tenes

#endif  // TENES_SRC_FERMION_FINITE_CHECK_HPP_
