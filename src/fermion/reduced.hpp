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

struct JointSwapForms {
  SwapForm cross;
  SwapForm bra;
};

inline JointSwapForms joint_swap_forms(const std::vector<int>& bra_axes,
                                       const std::vector<int>& ket_axes,
                                       const std::vector<int>& leg_ids) {
  if (bra_axes.size() != leg_ids.size() || ket_axes.size() != leg_ids.size()) {
    throw std::runtime_error("joint_swap_forms: axis and leg id size mismatch");
  }
  for (std::size_t i = 0; i < leg_ids.size(); ++i) {
    if (bra_axes[i] < 0 || ket_axes[i] < 0) {
      throw std::runtime_error("joint_swap_forms: negative axis");
    }
  }

  JointSwapForms forms;
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
          forms.cross.toggle(ket_axes[ix], bra_axes[iy]);
          forms.bra.toggle(bra_axes[ix], bra_axes[iy]);
        }
      }
    }
  }
  return forms;
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
tensor doubled_pipeline(const ftensor<tensor>& bra_Tn,
                        const ftensor<tensor>& ket_Tn) {
  // Shared double-layer pipeline: outer product of bra and ket layers,
  // frozen joint-swap dressing, (ket, bra) interleave, column-major fusion.
  // Output legs: ([l lb], [t tb], [r rb], [b bb], s_ket, s_bra).
  // The two arguments differ only where an operator has been inserted into
  // the ket layer beforehand; for the plain reduced tensor they are equal.
  const auto forms = joint_swap_forms({0, 1, 2, 3}, {5, 6, 7, 8}, {0, 1, 2, 3});
  ftensor<tensor> bra = conj(bra_Tn);
  apply_swap_form(bra, forms.bra);
  ftensor<tensor> doubled =
      tensordot(bra, ket_Tn, mptensor::Axes(), mptensor::Axes());
  mptensor::Axes interleaved;
  for (int ax = 0; ax < 4; ++ax) {
    interleaved.push(5 + ax);
    interleaved.push(ax);
  }
  interleaved.push(9);
  interleaved.push(4);
  transpose_with_swap_form(doubled, forms.cross, interleaved);
  mptensor::Shape sh;
  for (std::size_t ax = 0; ax < 4; ++ax) {
    sh.push(doubled.shape()[2 * ax] * doubled.shape()[2 * ax + 1]);
  }
  sh.push(doubled.shape()[8]);
  sh.push(doubled.shape()[9]);
  return mptensor::reshape(doubled.t, sh);
}

template <class tensor>
tensor fuse_doubled_cluster(const ftensor<tensor>& bra_pair,
                            const ftensor<tensor>& ket_pair,
                            const std::vector<int>& leg_ids) {
  constexpr std::size_t kExternalLegs = 6;
  const std::vector<int> cluster_axes = {0, 1, 2, 4, 5, 6};
  std::vector<int> bra_axes;
  std::vector<int> ket_axes;
  for (const int ax : cluster_axes) {
    bra_axes.push_back(ax);
    ket_axes.push_back(ax + 8);
  }
  const auto forms = joint_swap_forms(bra_axes, ket_axes, leg_ids);
  ftensor<tensor> bra = bra_pair;
  apply_swap_form(bra, forms.bra);
  ftensor<tensor> doubled =
      tensordot(bra, ket_pair, mptensor::Axes(), mptensor::Axes());

  mptensor::Axes interleaved;
  for (std::size_t i = 0; i < kExternalLegs; ++i) {
    interleaved.push(ket_axes[i]);
    interleaved.push(bra_axes[i]);
  }
  interleaved.push(11);
  interleaved.push(15);
  interleaved.push(3);
  interleaved.push(7);
  transpose_with_swap_form(doubled, forms.cross, interleaved);

  mptensor::Shape sh;
  for (std::size_t ax = 0; ax < kExternalLegs; ++ax) {
    sh.push(doubled.shape()[2 * ax] * doubled.shape()[2 * ax + 1]);
  }
  sh.push(doubled.shape()[12]);
  sh.push(doubled.shape()[13]);
  sh.push(doubled.shape()[14]);
  sh.push(doubled.shape()[15]);

  tensor fused = mptensor::reshape(doubled.t, sh);
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

template <class tensor>
tensor build_reduced(const ftensor<tensor>& Tn) {
  return mptensor::contract(build_reduced_op(Tn), mptensor::Axes(4),
                            mptensor::Axes(5));
}

template <class tensor>
ftensor<tensor> build_pair_state(const ftensor<tensor>& TnA,
                                 const ftensor<tensor>& TnB,
                                 reduced_pair_direction direction) {
  if (TnA.rank() != 5 || TnB.rank() != 5) {
    throw std::runtime_error("build_pair_state expects five-leg Tn ftensors");
  }
  switch (direction) {
    case reduced_pair_direction::horizontal:
      // (l_A, t_A, b_A, s_A, t_B, r_B, b_B, s_B)
      return tensordot(TnA, TnB, mptensor::Axes(2), mptensor::Axes(0));
    case reduced_pair_direction::vertical:
      // (l_A, t_A, r_A, s_A, l_B, r_B, b_B, s_B)
      return tensordot(TnA, TnB, mptensor::Axes(3), mptensor::Axes(1));
  }
  throw std::runtime_error("build_pair_state: invalid direction");
}

// Apply a two-site operator op12 (in_A, in_B, out_A, out_B) to the physical
// legs (3, 7) of a pair state and restore the original leg order.
template <class tensor>
ftensor<tensor> apply_pair_op(const ftensor<tensor>& pair,
                              const ftensor<tensor>& op12) {
  if (pair.rank() != 8) {
    throw std::runtime_error("apply_pair_op expects an eight-leg pair state");
  }
  if (op12.rank() != 4) {
    throw std::runtime_error("apply_pair_op expects a four-leg operator");
  }
  ftensor<tensor> applied =
      tensordot(pair, op12, mptensor::Axes(3, 7), mptensor::Axes(0, 1));
  return transpose(applied, mptensor::Axes(0, 1, 2, 6, 3, 4, 5, 7));
}

template <class tensor>
tensor build_reduced_pair(const ftensor<tensor>& TnA,
                          const ftensor<tensor>& TnB,
                          const ftensor<tensor>& op12,
                          reduced_pair_direction direction) {
  const ftensor<tensor> ket_ab = build_pair_state(TnA, TnB, direction);
  std::vector<int> leg_ids;
  switch (direction) {
    case reduced_pair_direction::horizontal:
      leg_ids = {0, 1, 3, 1, 2, 3};
      break;
    case reduced_pair_direction::vertical:
      leg_ids = {0, 1, 2, 0, 2, 3};
      break;
    default:
      throw std::runtime_error("doubled_cluster: invalid direction");
  }

  ftensor<tensor> ket_op = apply_pair_op(ket_ab, op12);
  tensor ret = detail::fuse_doubled_cluster(conj(ket_ab), ket_op, leg_ids);
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
  return ret;
}

}  // namespace tenes::fermion

#endif  // TENES_SRC_FERMION_REDUCED_HPP_
