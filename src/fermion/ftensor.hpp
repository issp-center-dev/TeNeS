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
  tensor t;
  leg_parities parity;

  mptensor::Shape shape() const { return t.shape(); }
  int rank() const { return static_cast<int>(parity.size()); }
  ftensor& transpose(const mptensor::Axes& axes) {
    for (std::size_t n = 0; n < t.local_size(); ++n) {
      auto idx = t.global_index(n);
      if (detail::transpose_sign(parity, idx, axes) < 0) {
        typename tensor::value_type v;
        t.get_value(idx, v);
        t.set_value(idx, -v);
      }
    }
    t.transpose(axes);
    leg_parities next;
    next.reserve(axes.size());
    for (std::size_t i = 0; i < axes.size(); ++i) {
      next.push_back(parity[axes[i]]);
    }
    parity = next;
    return *this;
  }
};

}  // namespace fermion
}  // namespace tenes

#endif  // TENES_SRC_FERMION_FTENSOR_HPP_
