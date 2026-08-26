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
 *  @brief The graded tensor type ::tenes::fermion::ftensor.
 *
 *  An ftensor pairs a plain mptensor tensor with a per-leg parity ledger
 *  (::tenes::fermion::leg_parities). All fermionic sign bookkeeping in TeNeS
 *  is expressed through this pairing: the tensor elements are stored
 *  unsigned, and the graded operations in fops.hpp / sign_sweep.hpp consult
 *  the ledger to apply the Koszul signs where legs are reordered or
 *  contracted.
 */

#ifndef TENES_SRC_FERMION_FTENSOR_HPP_
#define TENES_SRC_FERMION_FTENSOR_HPP_

#include "parity.hpp"

namespace tenes {
namespace fermion {
namespace detail {

/*!
 * @brief Koszul sign of one tensor element under a leg permutation.
 *
 * The sign of a permutation acting on graded legs is @f$(-1)^k@f$ where
 * @f$k@f$ counts the inversion pairs of the permutation whose two legs are
 * both odd for the given index: each such pair is one crossing of two odd
 * legs, contributing one factor of @f$-1@f$.
 *
 * @param[in] parity Parity ledgers of all legs (pre-transpose order).
 * @param[in] idx Element index (pre-transpose order).
 * @param[in] axes Transpose permutation, as passed to mptensor transpose.
 * @return +1 or -1.
 */
inline int transpose_sign(const leg_parities& parity,
                          const mptensor::Index& idx,
                          const mptensor::Axes& axes) {
  int k = 0;
  for (std::size_t x = 0; x < axes.size(); ++x) {
    for (std::size_t y = x + 1; y < axes.size(); ++y) {
      if (axes[x] > axes[y] && parity[axes[y]][idx[axes[y]]] &&
          parity[axes[x]][idx[axes[x]]]) {
        ++k;
      }
    }
  }
  return (k % 2 == 0) ? 1 : -1;
}

}  // namespace detail

/*! @brief A graded (fermionic) tensor: a plain tensor plus per-leg parity.
 *
 * The elements are stored exactly as in the underlying mptensor tensor
 * @c t; no signs are baked into them. The parity ledger @c parity tells,
 * for each leg and each index value, whether that basis state is fermion
 * odd, and is what the graded operations (tensordot, transpose, conj, ...)
 * consult to generate the Koszul signs.
 *
 * Most members simply forward to the wrapped tensor. The graded exceptions
 * are transpose(), which applies the swap-sign sweep before permuting the
 * legs, and the free functions of fops.hpp that take ftensor arguments.
 *
 * @tparam tensor mptensor tensor type (::tenes::real_tensor or
 *         ::tenes::complex_tensor).
 */
template <class tensor>
struct ftensor {
  //! Element type of the wrapped tensor.
  using value_type = typename tensor::value_type;

  //! Tensor elements, stored without any fermionic signs baked in.
  tensor t;
  //! One parity ledger per leg, in leg order; sizes match t.shape().
  leg_parities parity;

  ftensor() = default;
  //! Wrap an existing tensor with the given per-leg parity ledger.
  ftensor(const tensor& t_, const leg_parities& parity_)
      : t(t_), parity(parity_) {}
  //! Construct an all-even (purely bosonic) tensor of the given shape.
  ftensor(typename tensor::comm_type comm, const mptensor::Shape& shape)
      : t(comm, shape), parity(shape.size()) {
    for (std::size_t ax = 0; ax < shape.size(); ++ax) {
      parity[ax].assign(shape[ax], false);
    }
  }
  //! Shape of the wrapped tensor.
  mptensor::Shape shape() const { return t.shape(); }
  //! Number of legs.
  int rank() const { return static_cast<int>(parity.size()); }
  //! Communicator of the wrapped tensor.
  typename tensor::comm_type get_comm() const { return t.get_comm(); }
  //! Number of elements stored on this process.
  std::size_t local_size() const { return t.local_size(); }
  //! Global index of the n-th locally stored element.
  mptensor::Index global_index(std::size_t n) const {
    return t.global_index(n);
  }
  //! Read one element (collective, as in mptensor).
  template <class D>
  void get_value(const mptensor::Index& idx, D& v) const {
    t.get_value(idx, v);
  }
  //! Write one element (as in mptensor; only the owning process stores it).
  template <class D>
  void set_value(const mptensor::Index& idx, D v) {
    t.set_value(idx, v);
  }
  //! Divide all elements by a scalar. Parities are unchanged.
  ftensor& operator/=(double v) {
    t /= v;
    return *this;
  }
  //! Multiply a diagonal vector onto one leg. Parities are unchanged.
  template <typename D>
  ftensor& multiply_vector(const std::vector<D>& vec, std::size_t n_axes) {
    t.multiply_vector(vec, n_axes);
    return *this;
  }
  //! Multiply diagonal vectors onto two legs. Parities are unchanged.
  template <typename D0, typename D1>
  ftensor& multiply_vector(const std::vector<D0>& vec0, std::size_t n_axes0,
                           const std::vector<D1>& vec1, std::size_t n_axes1) {
    t.multiply_vector(vec0, n_axes0, vec1, n_axes1);
    return *this;
  }
  //! Multiply diagonal vectors onto three legs. Parities are unchanged.
  template <typename D0, typename D1, typename D2>
  ftensor& multiply_vector(const std::vector<D0>& vec0, std::size_t n_axes0,
                           const std::vector<D1>& vec1, std::size_t n_axes1,
                           const std::vector<D2>& vec2, std::size_t n_axes2) {
    t.multiply_vector(vec0, n_axes0, vec1, n_axes1, vec2, n_axes2);
    return *this;
  }
  /*!
   * @brief Graded in-place transpose.
   *
   * Applies the Koszul sign of the permutation (see
   * detail::transpose_sign()) to the elements, then permutes both the
   * tensor legs and the parity ledgers. Defined in sign_sweep.hpp.
   *
   * @param[in] axes Permutation, as in mptensor transpose.
   */
  ftensor& transpose(const mptensor::Axes& axes);
};

}  // namespace fermion
}  // namespace tenes

#endif  // TENES_SRC_FERMION_FTENSOR_HPP_
