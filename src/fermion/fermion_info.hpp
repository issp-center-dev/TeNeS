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
 *  @brief Per-site parity ledgers of the iTPS and the Tn wrap/unwrap layer.
 *
 *  ::tenes::fermion::FermionInfo carries the parity ledgers of every site
 *  tensor in the unit cell. wrap_Tn() combines a plain center tensor Tn
 *  with these ledgers into a graded ::tenes::fermion::ftensor; unwrap_Tn()
 *  extracts the plain tensor back. The wrapped Tn is rank 5 with the leg
 *  order shared by the whole solver:
 *
 *  @verbatim
          t(1)
           |
    l(0) - Tn - r(2)
           |  \
          b(3) s(4 = physical)
    @endverbatim
 */

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

/*! @brief Parity ledgers of every site tensor in the unit cell.
 *
 * When @c enabled is false the run is bosonic and the ledgers are unused.
 * The virtual ledgers are shared between neighbors: leg @c l of one site and
 * leg @c (l+2)%4 of its neighbor describe the same bond and must be equal
 * (checked by validate_neighbor_consistency()).
 */
struct FermionInfo {
  //! True iff the run treats the model fermionically.
  bool enabled = false;
  //! Physical-leg parity ledger of each site, indexed by site.
  std::vector<parity_vector> phys;
  //! Virtual-leg ledgers of each site, indexed by site then by leg
  //! (0=left, 1=top, 2=right, 3=bottom).
  std::vector<std::array<parity_vector, 4>> virt;
};

/*!
 * @brief Check that shared bonds carry identical parity ledgers.
 *
 * @param[in] fi Ledger set to check; a disabled ledger passes trivially.
 * @param[in] lattice Unit-cell geometry supplying the neighbor map.
 * @throw std::runtime_error If the ledger count does not match the unit
 *        cell, or if any bond's two ends disagree.
 */
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

/*!
 * @brief Even-first parity ledger for a leg of the given dimension.
 *
 * The first ceil(dim/2) index values are even, the rest odd. This is the
 * layout the graded decompositions produce for virtual legs, and the
 * default assumed when a ledger has to be invented (e.g. on initial state
 * construction).
 *
 * @param[in] dim Leg dimension.
 */
inline parity_vector even_first_parity(std::size_t dim) {
  parity_vector p(dim, false);
  const std::size_t neven = (dim + 1) / 2;
  for (std::size_t i = neven; i < dim; ++i) {
    p[i] = true;
  }
  return p;
}

/*!
 * @brief Ledgers of one site's Tn in the wrapped leg order (l, t, r, b, s).
 *
 * @param[in] fi Ledger set; must be populated for @p site.
 * @param[in] site Site index within the unit cell.
 * @throw std::runtime_error If @p site is out of range.
 */
inline leg_parities Tn_parity(const FermionInfo& fi, int site) {
  if (site < 0 || static_cast<std::size_t>(site) >= fi.phys.size() ||
      static_cast<std::size_t>(site) >= fi.virt.size()) {
    throw std::runtime_error("FermionInfo: site is out of range");
  }
  return {fi.virt[site][0], fi.virt[site][1], fi.virt[site][2],
          fi.virt[site][3], fi.phys[site]};
}

/*!
 * @brief Wrap a plain center tensor into a graded ftensor.
 *
 * @param[in] Tn Rank-5 center tensor with legs (l, t, r, b, s).
 * @param[in] fi Ledger set of the unit cell.
 * @param[in] site Site index of @p Tn within the unit cell.
 * @return @p Tn paired with the ledgers Tn_parity(fi, site).
 */
template <class tensor>
ftensor<tensor> wrap_Tn(const tensor& Tn, const FermionInfo& fi, int site) {
  return ftensor<tensor>{Tn, Tn_parity(fi, site)};
}

/*!
 * @brief Extract the plain tensor from a wrapped Tn.
 *
 * The inverse of wrap_Tn() for the element data: since ftensor stores its
 * elements unsigned, this is a plain copy-out. The ledger is not written
 * back — the FermionInfo argument only mirrors the wrap_Tn() signature.
 *
 * @param[in] ft Wrapped rank-5 center tensor.
 * @param[out] Tn Receives the plain tensor.
 * @throw std::runtime_error If @p ft is not rank 5.
 */
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
