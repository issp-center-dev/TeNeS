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
#include <vector>

#include "fops.hpp"
#include "ftensor.hpp"

namespace tenes::fermion {

enum class reduced_pair_direction { horizontal, vertical };

namespace detail {

constexpr int joint_bit(int x, int y) { return x * 3 + (y < x ? y : y - 1); }

// This joint-swap pattern is gauge-equivalent to YASTN's fuse_layers()
// convention, empirically pinned down using the oracle as arbiter.
constexpr unsigned kDoubledJointMask =
    (1u << joint_bit(0, 3)) | (1u << joint_bit(1, 0)) |
    (1u << joint_bit(2, 3)) | (1u << joint_bit(3, 0));

template <class tensor>
void apply_joint_swaps(ftensor<tensor>& a, const std::vector<int>& bra_axes,
                       const std::vector<int>& ket_axes,
                       const std::vector<int>& leg_ids) {
  for (int x = 0; x < 4; ++x) {
    for (int y = 0; y < 4; ++y) {
      if (x == y) {
        continue;
      }
      if ((kDoubledJointMask & (1u << joint_bit(x, y))) == 0) {
        continue;
      }
      for (std::size_t ix = 0; ix < leg_ids.size(); ++ix) {
        if (leg_ids[ix] != x) {
          continue;
        }
        for (std::size_t iy = 0; iy < leg_ids.size(); ++iy) {
          if (leg_ids[iy] != y) {
            continue;
          }
          apply_swap(a, ket_axes[ix], bra_axes[iy]);
          apply_swap(a, bra_axes[ix], bra_axes[iy]);
        }
      }
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
                                  const std::vector<int>& leg_ids) {
  const std::size_t nlegs = leg_ids.size();
  ftensor<tensor> prepared = doubled;
  std::vector<int> bra_axes;
  std::vector<int> ket_axes;
  for (std::size_t ax = 0; ax < nlegs; ++ax) {
    bra_axes.push_back(static_cast<int>(ax));
    ket_axes.push_back(static_cast<int>(nlegs + ax));
  }
  apply_joint_swaps(prepared, bra_axes, ket_axes, leg_ids);

  mptensor::Axes interleaved;
  for (std::size_t ax = 0; ax < nlegs; ++ax) {
    interleaved.push(nlegs + ax);
    interleaved.push(ax);
  }
  ftensor<tensor> ordered = transpose(prepared, interleaved);
  mptensor::Shape sh;
  for (std::size_t ax = 0; ax < nlegs; ++ax) {
    sh.push(ordered.shape()[2 * ax] * ordered.shape()[2 * ax + 1]);
  }
  return mptensor::reshape(ordered.t, sh);
}

template <class tensor>
void scale_in_place(tensor& a, double scale) {
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    const auto idx = a.global_index(n);
    typename tensor::value_type v;
    a.get_value(idx, v);
    a.set_value(idx, scale * v);
  }
}

template <class tensor>
void apply_fused_leg_gauge(tensor& a, const parity_vector& leg_parity,
                           std::size_t ax, bool ket_odd_bra_even) {
  std::vector<double> sign(leg_parity.size() * leg_parity.size(), 1.0);
  for (std::size_t bra = 0; bra < leg_parity.size(); ++bra) {
    for (std::size_t ket = 0; ket < leg_parity.size(); ++ket) {
      const bool flip = ket_odd_bra_even
                            ? (leg_parity[ket] && !leg_parity[bra])
                            : (!leg_parity[ket] && leg_parity[bra]);
      sign[ket + leg_parity.size() * bra] = flip ? -1.0 : 1.0;
    }
  }
  a.multiply_vector(sign, ax);
}

template <class tensor>
tensor doubled_pipeline(const ftensor<tensor>& braTn,
                        const ftensor<tensor>& ketT) {
  // Shared double-layer pipeline: outer product of bra and ket layers,
  // frozen joint-swap dressing, (ket, bra) interleave, column-major fusion.
  // Output legs: ([l lb], [t tb], [r rb], [b bb], s_ket, s_bra).
  ftensor<tensor> doubled =
      tensordot(conj(braTn), ketT, mptensor::Axes(), mptensor::Axes());
  apply_joint_swaps(doubled, {0, 1, 2, 3}, {5, 6, 7, 8}, {0, 1, 2, 3});
  mptensor::Axes interleaved;
  for (int ax = 0; ax < 4; ++ax) {
    interleaved.push(5 + ax);
    interleaved.push(ax);
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
tensor fuse_doubled_cluster(const ftensor<tensor>& doubled,
                            const std::vector<int>& leg_ids) {
  ftensor<tensor> prepared = doubled;
  constexpr std::size_t kExternalLegs = 6;
  const std::vector<int> cluster_axes = {0, 1, 2, 4, 5, 6};
  std::vector<int> bra_axes;
  std::vector<int> ket_axes;
  for (const int ax : cluster_axes) {
    bra_axes.push_back(ax);
    ket_axes.push_back(ax + 8);
  }
  apply_joint_swaps(prepared, bra_axes, ket_axes, leg_ids);

  mptensor::Axes interleaved;
  for (std::size_t i = 0; i < kExternalLegs; ++i) {
    interleaved.push(ket_axes[i]);
    interleaved.push(bra_axes[i]);
  }
  interleaved.push(11);
  interleaved.push(15);
  interleaved.push(3);
  interleaved.push(7);
  ftensor<tensor> ordered = transpose(prepared, interleaved);

  mptensor::Shape sh;
  for (std::size_t ax = 0; ax < kExternalLegs; ++ax) {
    sh.push(ordered.shape()[2 * ax] * ordered.shape()[2 * ax + 1]);
  }
  sh.push(ordered.shape()[12]);
  sh.push(ordered.shape()[13]);
  sh.push(ordered.shape()[14]);
  sh.push(ordered.shape()[15]);

  tensor fused = mptensor::reshape(ordered.t, sh);
  return mptensor::contract(fused, mptensor::Axes(6, 7), mptensor::Axes(8, 9));
}

}  // namespace detail

template <class tensor>
tensor build_reduced_op(const ftensor<tensor>& Tn) {
  if (Tn.rank() != 5) {
    throw std::runtime_error("build_reduced_op expects a five-leg Tn ftensor");
  }
  return detail::doubled_pipeline(Tn, Tn);
}

// Doubled impurity with a (possibly parity-odd) single-site operator inserted
// into the KET layer BEFORE doubling. Odd operators anticommute past virtual
// legs during the interleave, so post-doubling bosonic insertion is wrong for
// them; the pre-doubling insertion lets the graded pipeline generate those
// signs. For parity-even operators both insertion points agree.
template <class tensor>
tensor build_reduced_channel(const ftensor<tensor>& Tn,
                             const ftensor<tensor>& op_channel) {
  ftensor<tensor> ket_op =
      tensordot(Tn, op_channel, mptensor::Axes(4), mptensor::Axes(0));
  return mptensor::contract(detail::doubled_pipeline(Tn, ket_op),
                            mptensor::Axes(4), mptensor::Axes(5));
}

template <class tensor>
tensor build_reduced(const ftensor<tensor>& Tn) {
  return mptensor::contract(build_reduced_op(Tn), mptensor::Axes(4),
                            mptensor::Axes(5));
}

// variant bit 0: doubled-cluster global sign +1 (0) or -1 (1).
inline int g_reduced_pair_variant = 0;

template <class tensor>
tensor doubled_cluster(const ftensor<tensor>& TnA, const ftensor<tensor>& TnB,
                       const ftensor<tensor>& op12,
                       reduced_pair_direction direction, int variant = 0);

template <class tensor>
tensor build_reduced_pair(const ftensor<tensor>& TnA,
                          const ftensor<tensor>& TnB,
                          const ftensor<tensor>& op12,
                          reduced_pair_direction direction) {
  return doubled_cluster(TnA, TnB, op12, direction, g_reduced_pair_variant);
}

template <class tensor>
tensor doubled_cluster(const ftensor<tensor>& TnA, const ftensor<tensor>& TnB,
                       const ftensor<tensor>& op12,
                       reduced_pair_direction direction, int variant) {
  ftensor<tensor> ket_ab;
  std::vector<int> leg_ids;
  switch (direction) {
    case reduced_pair_direction::horizontal:
      ket_ab = tensordot(TnA, TnB, mptensor::Axes(2), mptensor::Axes(0));
      leg_ids = {0, 1, 3, 1, 2, 3};
      break;
    case reduced_pair_direction::vertical:
      ket_ab = tensordot(TnA, TnB, mptensor::Axes(3), mptensor::Axes(1));
      leg_ids = {0, 1, 2, 0, 2, 3};
      break;
    default:
      throw std::runtime_error("doubled_cluster: invalid direction");
  }

  ftensor<tensor> ket_op = detail::apply_reduced_two_site_op(ket_ab, op12);
  ftensor<tensor> doubled =
      tensordot(conj(ket_ab), ket_op, mptensor::Axes(), mptensor::Axes());
  tensor ret = detail::fuse_doubled_cluster(doubled, leg_ids);
  // Gauge alignment with the single-site doubling convention; measured
  // directly via comparison against build_reduced_op/build_reduced-based
  // direct composition, not derived analytically.
  if (direction == reduced_pair_direction::horizontal) {
    detail::apply_fused_leg_gauge(ret, TnA.parity[3], 2, true);
    detail::apply_fused_leg_gauge(ret, TnB.parity[3], 5, false);
  } else {
    detail::apply_fused_leg_gauge(ret, TnA.parity[0], 0, true);
    detail::apply_fused_leg_gauge(ret, TnB.parity[0], 3, false);
  }
  if ((variant & 1) != 0) {
    detail::scale_in_place(ret, -1.0);
  }
  return ret;
}

}  // namespace tenes::fermion

#endif  // TENES_SRC_FERMION_REDUCED_HPP_
