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

#ifndef TENES_SRC_FERMION_PARITY_HPP_
#define TENES_SRC_FERMION_PARITY_HPP_

#include <cstddef>
#include <vector>

#include <mptensor/tensor.hpp>

namespace tenes {
namespace fermion {

using parity_vector = std::vector<bool>;
using leg_parities = std::vector<parity_vector>;

// Fuse two leg parities into the parity of the combined index, using the
// COLUMN-MAJOR flattening convention of mptensor (ScaLAPACK heritage):
// mptensor::reshape combines adjacent legs as index = i + a.size() * j, i.e.
// the FIRST leg is the fastest-varying one. A row-major fuse here silently
// mislabels every fused index whose leg parity patterns are asymmetric.
inline parity_vector fuse(const parity_vector& a, const parity_vector& b) {
  parity_vector r(a.size() * b.size());
  for (std::size_t j = 0; j < b.size(); ++j) {
    for (std::size_t i = 0; i < a.size(); ++i) {
      r[i + a.size() * j] = a[i] != b[j];
    }
  }
  return r;
}

inline int count_odd(const leg_parities& p, const mptensor::Index& idx) {
  int m = 0;
  for (std::size_t ax = 0; ax < p.size(); ++ax) {
    m += p[ax][idx[ax]] ? 1 : 0;
  }
  return m;
}

inline std::vector<std::size_t> parity_sort_perm(const parity_vector& p) {
  std::vector<std::size_t> perm;
  perm.reserve(p.size());
  for (std::size_t i = 0; i < p.size(); ++i) {
    if (!p[i]) {
      perm.push_back(i);
    }
  }
  for (std::size_t i = 0; i < p.size(); ++i) {
    if (p[i]) {
      perm.push_back(i);
    }
  }
  return perm;
}

}  // namespace fermion
}  // namespace tenes

#endif  // TENES_SRC_FERMION_PARITY_HPP_
