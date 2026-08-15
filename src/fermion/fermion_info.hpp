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

#ifndef TENES_SRC_FERMION_FERMION_INFO_HPP_
#define TENES_SRC_FERMION_FERMION_INFO_HPP_

#include <array>
#include <stdexcept>
#include <vector>

#include "ftensor.hpp"
#include "parity.hpp"

namespace tenes {
namespace fermion {

struct FermionInfo {
  bool enabled = false;
  std::vector<parity_vector> phys;
  std::vector<std::array<parity_vector, 4>> virt;
  std::vector<std::array<parity_vector, 2>> C_par[4];
  std::vector<std::array<parity_vector, 4>> eT_par[4];
};

inline parity_vector even_first_parity(std::size_t dim) {
  parity_vector p(dim, false);
  const std::size_t neven = (dim + 1) / 2;
  for (std::size_t i = neven; i < dim; ++i) {
    p[i] = true;
  }
  return p;
}

inline leg_parities Tn_parity(const FermionInfo& fi, int site) {
  if (site < 0 || static_cast<std::size_t>(site) >= fi.phys.size() ||
      static_cast<std::size_t>(site) >= fi.virt.size()) {
    throw std::runtime_error("FermionInfo: site is out of range");
  }
  return {fi.virt[site][0], fi.virt[site][1], fi.virt[site][2],
          fi.virt[site][3], fi.phys[site]};
}

template <class tensor>
ftensor<tensor> wrap_Tn(const tensor& Tn, const FermionInfo& fi, int site) {
  return ftensor<tensor>{Tn, Tn_parity(fi, site)};
}

template <class tensor>
void unwrap_Tn(const ftensor<tensor>& ft, tensor& Tn, FermionInfo&, int) {
  if (ft.parity.size() != 5) {
    throw std::runtime_error("unwrap_Tn expects a five-leg Tn ftensor");
  }
  Tn = ft.t;
}

}  // namespace fermion
}  // namespace tenes

#endif  // TENES_SRC_FERMION_FERMION_INFO_HPP_
