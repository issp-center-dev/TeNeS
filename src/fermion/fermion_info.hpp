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
#include <sstream>
#include <vector>

#include "../SquareLattice.hpp"
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

inline void validate_neighbor_consistency(const FermionInfo& fi,
                                          const SquareLattice& lattice) {
  if (!fi.enabled) {
    return;
  }
  if (fi.virt.size() != static_cast<std::size_t>(lattice.N_UNIT)) {
    throw std::runtime_error("FermionInfo: virtual ledger size mismatch");
  }
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    for (int leg = 0; leg < 4; ++leg) {
      const int neighbor = lattice.neighbor(site, leg);
      const int neighbor_leg = (leg + 2) % 4;
      if (fi.virt[site][leg] != fi.virt[neighbor][neighbor_leg]) {
        std::stringstream ss;
        ss << "FermionInfo: inconsistent virtual parity ledger at site " << site
           << " leg " << leg << " versus neighbor " << neighbor << " leg "
           << neighbor_leg;
        throw std::runtime_error(ss.str());
      }
    }
  }
}

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

inline leg_parities C_parity(const FermionInfo& fi, int corner, int site) {
  if (corner < 0 || corner >= 4 || site < 0 ||
      static_cast<std::size_t>(site) >= fi.C_par[corner].size()) {
    throw std::runtime_error("FermionInfo: C index is out of range");
  }
  return {fi.C_par[corner][site][0], fi.C_par[corner][site][1]};
}

template <class tensor>
ftensor<tensor> wrap_C(const tensor& C, const FermionInfo& fi, int corner,
                       int site) {
  return ftensor<tensor>{C, C_parity(fi, corner, site)};
}

template <class tensor>
void unwrap_C(const ftensor<tensor>& ft, tensor& C, FermionInfo& fi, int corner,
              int site) {
  if (ft.parity.size() != 2) {
    throw std::runtime_error("unwrap_C expects a two-leg C ftensor");
  }
  C = ft.t;
  fi.C_par[corner][site] = {ft.parity[0], ft.parity[1]};
}

inline leg_parities eT_parity(const FermionInfo& fi, int edge, int site) {
  if (edge < 0 || edge >= 4 || site < 0 ||
      static_cast<std::size_t>(site) >= fi.eT_par[edge].size()) {
    throw std::runtime_error("FermionInfo: eT index is out of range");
  }
  return {fi.eT_par[edge][site][0], fi.eT_par[edge][site][1],
          fi.eT_par[edge][site][2], fi.eT_par[edge][site][3]};
}

template <class tensor>
ftensor<tensor> wrap_eT(const tensor& eT, const FermionInfo& fi, int edge,
                        int site) {
  return ftensor<tensor>{eT, eT_parity(fi, edge, site)};
}

template <class tensor>
void unwrap_eT(const ftensor<tensor>& ft, tensor& eT, FermionInfo& fi, int edge,
               int site) {
  if (ft.parity.size() != 4) {
    throw std::runtime_error("unwrap_eT expects a four-leg eT ftensor");
  }
  eT = ft.t;
  fi.eT_par[edge][site] = {ft.parity[0], ft.parity[1], ft.parity[2],
                           ft.parity[3]};
}

}  // namespace fermion
}  // namespace tenes

#endif  // TENES_SRC_FERMION_FERMION_INFO_HPP_
