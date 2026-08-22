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

#ifndef TENES_SRC_FERMION_FTENSOR_HPP_
#define TENES_SRC_FERMION_FTENSOR_HPP_

#include "parity.hpp"

namespace tenes {
namespace fermion {
namespace detail {

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

template <class tensor>
struct ftensor {
  using value_type = typename tensor::value_type;

  tensor t;
  leg_parities parity;

  ftensor() = default;
  ftensor(const tensor& t_, const leg_parities& parity_)
      : t(t_), parity(parity_) {}
  ftensor(typename tensor::comm_type comm, const mptensor::Shape& shape)
      : t(comm, shape), parity(shape.size()) {
    for (std::size_t ax = 0; ax < shape.size(); ++ax) {
      parity[ax].assign(shape[ax], false);
    }
  }
  mptensor::Shape shape() const { return t.shape(); }
  int rank() const { return static_cast<int>(parity.size()); }
  typename tensor::comm_type get_comm() const { return t.get_comm(); }
  std::size_t local_size() const { return t.local_size(); }
  mptensor::Index global_index(std::size_t n) const {
    return t.global_index(n);
  }
  template <class D>
  void get_value(const mptensor::Index& idx, D& v) const {
    t.get_value(idx, v);
  }
  template <class D>
  void set_value(const mptensor::Index& idx, D v) {
    t.set_value(idx, v);
  }
  ftensor& operator/=(double v) {
    t /= v;
    return *this;
  }
  template <typename D>
  ftensor& multiply_vector(const std::vector<D>& vec, std::size_t n_axes) {
    t.multiply_vector(vec, n_axes);
    return *this;
  }
  template <typename D0, typename D1>
  ftensor& multiply_vector(const std::vector<D0>& vec0, std::size_t n_axes0,
                           const std::vector<D1>& vec1, std::size_t n_axes1) {
    t.multiply_vector(vec0, n_axes0, vec1, n_axes1);
    return *this;
  }
  template <typename D0, typename D1, typename D2>
  ftensor& multiply_vector(const std::vector<D0>& vec0, std::size_t n_axes0,
                           const std::vector<D1>& vec1, std::size_t n_axes1,
                           const std::vector<D2>& vec2, std::size_t n_axes2) {
    t.multiply_vector(vec0, n_axes0, vec1, n_axes1, vec2, n_axes2);
    return *this;
  }
  ftensor& transpose(const mptensor::Axes& axes);
};

}  // namespace fermion
}  // namespace tenes

#endif  // TENES_SRC_FERMION_FTENSOR_HPP_
