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

#ifndef TENES_SRC_FERMION_REDUCED_HPP_
#define TENES_SRC_FERMION_REDUCED_HPP_

#include <stdexcept>

#include "fops.hpp"
#include "ftensor.hpp"

namespace tenes::fermion {

enum class reduced_pair_direction { horizontal, vertical };

struct reduced_variant {
  unsigned bra_mask = 0;
  unsigned ket_mask = 0;
  bool ket_first = true;
  bool phys_twist = false;
};

namespace detail {

inline int leg_bit(int leg) { return 1 << leg; }

template <class tensor>
void apply_spatial_parities(ftensor<tensor>& a,
                            const std::vector<int>& bra_axes,
                            const std::vector<int>& ket_axes,
                            const std::vector<int>& leg_ids,
                            const reduced_variant& variant) {
  for (std::size_t i = 0; i < leg_ids.size(); ++i) {
    const int bit = leg_bit(leg_ids[i]);
    if ((variant.bra_mask & bit) != 0) {
      apply_parity(a, bra_axes[i]);
    }
    if ((variant.ket_mask & bit) != 0) {
      apply_parity(a, ket_axes[i]);
    }
  }
}

template <class tensor>
ftensor<tensor> apply_reduced_two_site_op(const ftensor<tensor>& psi,
                                          const ftensor<tensor>& op) {
  ftensor<tensor> applied =
      tensordot(psi, op, mptensor::Axes(3, 7), mptensor::Axes(0, 1));
  return transpose(applied, mptensor::Axes(0, 1, 2, 6, 3, 4, 5, 7));
}

template <class tensor>
tensor fuse_doubled_external_legs(const ftensor<tensor>& doubled,
                                  const std::vector<int>& leg_ids,
                                  const reduced_variant& variant) {
  const std::size_t nlegs = leg_ids.size();
  ftensor<tensor> prepared = doubled;
  std::vector<int> bra_axes;
  std::vector<int> ket_axes;
  for (std::size_t ax = 0; ax < nlegs; ++ax) {
    bra_axes.push_back(static_cast<int>(ax));
    ket_axes.push_back(static_cast<int>(nlegs + ax));
  }
  apply_spatial_parities(prepared, bra_axes, ket_axes, leg_ids, variant);

  mptensor::Axes interleaved;
  for (std::size_t ax = 0; ax < nlegs; ++ax) {
    if (variant.ket_first) {
      interleaved.push(nlegs + ax);
      interleaved.push(ax);
    } else {
      interleaved.push(ax);
      interleaved.push(nlegs + ax);
    }
  }
  ftensor<tensor> ordered = transpose(prepared, interleaved);
  mptensor::Shape sh;
  for (std::size_t ax = 0; ax < nlegs; ++ax) {
    sh.push(ordered.shape()[2 * ax] * ordered.shape()[2 * ax + 1]);
  }
  return mptensor::reshape(ordered.t, sh);
}

}  // namespace detail

template <class tensor>
tensor build_reduced_op(const ftensor<tensor>& Tn,
                        const reduced_variant& variant) {
  if (Tn.rank() != 5) {
    throw std::runtime_error("build_reduced_op expects a five-leg Tn ftensor");
  }
  ftensor<tensor> doubled =
      tensordot(conj(Tn), Tn, mptensor::Axes(), mptensor::Axes());
  detail::apply_spatial_parities(doubled, {0, 1, 2, 3}, {5, 6, 7, 8},
                                 {0, 1, 2, 3}, variant);
  if (variant.phys_twist) {
    apply_parity(doubled, 4);
  }
  mptensor::Axes interleaved;
  for (int ax = 0; ax < 4; ++ax) {
    if (variant.ket_first) {
      interleaved.push(5 + ax);
      interleaved.push(ax);
    } else {
      interleaved.push(ax);
      interleaved.push(5 + ax);
    }
  }
  interleaved.push(9);
  interleaved.push(4);
  ftensor<tensor> ordered = transpose(doubled, interleaved);
  mptensor::Shape sh;
  for (std::size_t ax = 0; ax < 4; ++ax) {
    sh.push(ordered.shape()[2 * ax] * ordered.shape()[2 * ax + 1]);
  }
  sh.push(ordered.shape()[8]);
  sh.push(ordered.shape()[9]);
  return mptensor::reshape(ordered.t, sh);
}

template <class tensor>
tensor build_reduced_op(const ftensor<tensor>& Tn) {
  return build_reduced_op(Tn, reduced_variant{});
}

template <class tensor>
tensor build_reduced(const ftensor<tensor>& Tn,
                     const reduced_variant& variant) {
  return mptensor::contract(build_reduced_op(Tn, variant), mptensor::Axes(4),
                            mptensor::Axes(5));
}

template <class tensor>
tensor build_reduced_pair(const ftensor<tensor>& TnA,
                          const ftensor<tensor>& TnB,
                          const ftensor<tensor>& op12,
                          reduced_pair_direction direction,
                          const reduced_variant& variant) {
  ftensor<tensor> ket;
  std::vector<int> leg_ids;
  switch (direction) {
    case reduced_pair_direction::horizontal:
      ket = tensordot(TnA, TnB, mptensor::Axes(2), mptensor::Axes(0));
      leg_ids = {0, 1, 3, 1, 2, 3};
      break;
    case reduced_pair_direction::vertical:
      ket = tensordot(TnA, TnB, mptensor::Axes(3), mptensor::Axes(1));
      leg_ids = {0, 1, 2, 0, 2, 3};
      break;
    default:
      throw std::runtime_error("build_reduced_pair: invalid direction");
  }
  ftensor<tensor> op_ket = detail::apply_reduced_two_site_op(ket, op12);
  ftensor<tensor> bra = conj(ket);
  if (variant.phys_twist) {
    apply_parity(bra, 3);
    apply_parity(bra, 7);
  }
  ftensor<tensor> doubled =
      tensordot(bra, op_ket, mptensor::Axes(3, 7), mptensor::Axes(3, 7));
  return detail::fuse_doubled_external_legs(doubled, leg_ids, variant);
}

template <class tensor>
tensor build_reduced(const ftensor<tensor>& Tn) {
  return build_reduced(Tn, reduced_variant{});
}

template <class tensor>
tensor build_reduced_pair(const ftensor<tensor>& TnA,
                          const ftensor<tensor>& TnB,
                          const ftensor<tensor>& op12,
                          reduced_pair_direction direction) {
  return build_reduced_pair(TnA, TnB, op12, direction, reduced_variant{});
}

}  // namespace tenes::fermion

#endif  // TENES_SRC_FERMION_REDUCED_HPP_
