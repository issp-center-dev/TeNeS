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
 *  @brief Relay of a two-site operator channel along the virtual bonds of a
 *         measurement window.
 *
 *  A fermionic two-site operator is split by a graded SVD into channels of
 *  definite parity. With the channel fixed, its k leg is a dimension-1
 *  auxiliary leg kappa. For an odd channel kappa is a string that is relayed
 *  from the source to the target along a nearest-neighbour path of the
 *  window: each site on the path fuses kappa into its entry and exit bonds,
 *  so that every tensor of the network is parity even again and the folded
 *  window can be handed to the bosonic density kernels. When the string flows
 *  against the canonical bond orientation (right-to-left or bottom-to-top),
 *  the upstream site carries the corresponding kappa supertrace sign. See
 *  docs/superpowers/specs/2026-09-25-fermion-longrange-measure-design.md.
 */

#ifndef TENES_SRC_FERMION_RELAY_HPP_
#define TENES_SRC_FERMION_RELAY_HPP_

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "fops.hpp"
#include "ftensor.hpp"
#include "reduced.hpp"

namespace tenes::fermion {

namespace detail {

//! True iff leg is one of the four virtual legs (l, t, r, b).
inline bool valid_relay_leg(int leg) { return 0 <= leg && leg < 4; }

//! Dimension-1 identity on the incoming and outgoing kappa string legs.
template <class tensor>
ftensor<tensor> relay_identity_string(typename tensor::comm_type comm,
                                      bool kappa_parity) {
  parity_vector kappa{kappa_parity};
  ftensor<tensor> S{tensor(comm, mptensor::Shape(1, 1)), {kappa, kappa}};
  S.set_value(mptensor::Index(0, 0), typename tensor::value_type(1));
  return S;
}

//! Apply the (-1)^(p_bond * p_kappa) sign for moving kappa past a bond leg.
template <class tensor>
void apply_relay_crossing_mask(ftensor<tensor>& a, int axis,
                               const parity_vector& bond_parity,
                               bool kappa_parity) {
  std::vector<double> mask(bond_parity.size(), 1.0);
  if (kappa_parity) {
    for (std::size_t i = 0; i < bond_parity.size(); ++i) {
      if (bond_parity[i]) {
        mask[i] = -1.0;
      }
    }
  }
  a.multiply_vector(mask, axis);
}

//! Multiply every local element of an ftensor by a scalar.
template <class tensor>
void multiply_ftensor_scalar(ftensor<tensor>& a,
                             typename tensor::value_type factor) {
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    a.t[n] *= factor;
  }
}

//! Append a virtual leg and its optional adjacent kappa legs to an axis list.
inline void relay_push_virtual_with_kappa(mptensor::Axes& axes, int leg,
                                          int entry_leg, int entry_axis,
                                          int exit_leg, int exit_axis) {
  axes.push(leg);
  if (leg == entry_leg) {
    axes.push(entry_axis);
  }
  if (leg == exit_leg) {
    axes.push(exit_axis);
  }
}

//! Move kappa next to its virtual bond legs and fuse those pairs explicitly.
//!
//! Do not replace this with the shared tenes::fermion::reshape() from
//! fops.hpp: kappa has dimension 1, so dimensions alone cannot distinguish
//! whether it should be fused into the preceding bond or left as the next leg.
//! A future reshape() substitution would therefore risk fusing the dimension-1
//! kappa leg into the wrong neighbour.
template <class tensor>
ftensor<tensor> fuse_relay_kappa(const ftensor<tensor>& a,
                                 const mptensor::Axes& axes,
                                 const std::vector<bool>& has_kappa) {
  const ftensor<tensor> moved = transpose(a, axes);
  ftensor<tensor> ret;
  mptensor::Shape sh;
  std::size_t pos = 0;
  for (int leg = 0; leg < 4; ++leg) {
    std::size_t dim = moved.shape()[pos++];
    if (has_kappa[leg]) {
      dim *= moved.shape()[pos++];
    }
    sh.push(dim);
  }
  sh.push(moved.shape()[pos]);
  ret.t = mptensor::reshape(moved.t, sh);
  pos = 0;
  for (int leg = 0; leg < 4; ++leg) {
    parity_vector p = moved.parity[pos++];
    if (has_kappa[leg]) {
      p = fuse(p, moved.parity[pos++]);
    }
    ret.parity.push_back(p);
  }
  ret.parity.push_back(moved.parity[pos]);
  return ret;
}

}  // namespace detail

//! Cell of a measurement window in Contract_* orientation:
//! row 0 is the top row, col 0 the left column.
struct window_cell {
  int row;
  int col;
};

//! Order in which relay_path() walks the two axes.
enum class relay_order { x_first, y_first };

//! Role of a site on the relay path.
enum class relay_role { source, middle, target };

/*!
 * @brief Cells from source to target inclusive, one nearest-neighbour step
 *        at a time.
 *
 * x_first walks along the source row to the target column, then along that
 * column; y_first the other way round.
 *
 * @throw std::invalid_argument if source == target.
 */
inline std::vector<window_cell> relay_path(window_cell source,
                                           window_cell target,
                                           relay_order order) {
  if (source.row == target.row && source.col == target.col) {
    throw std::invalid_argument("relay_path: source equals target");
  }
  std::vector<window_cell> path;
  path.push_back(source);
  window_cell cur = source;
  const auto walk_col = [&]() {
    while (cur.col != target.col) {
      cur.col += cur.col < target.col ? 1 : -1;
      path.push_back(cur);
    }
  };
  const auto walk_row = [&]() {
    while (cur.row != target.row) {
      cur.row += cur.row < target.row ? 1 : -1;
      path.push_back(cur);
    }
  };
  if (order == relay_order::x_first) {
    walk_col();
    walk_row();
  } else {
    walk_row();
    walk_col();
  }
  return path;
}

/*!
 * @brief Virtual leg (0 = l, 1 = t, 2 = r, 3 = b) of @p from pointing at the
 *        adjacent cell @p to (row - 1 is up = t).
 *
 * @throw std::invalid_argument if the cells are not adjacent.
 */
inline int relay_leg(window_cell from, window_cell to) {
  const int dr = to.row - from.row;
  const int dc = to.col - from.col;
  if (dr == 0 && dc == 1) {
    return 2;
  }
  if (dr == 0 && dc == -1) {
    return 0;
  }
  if (dr == -1 && dc == 0) {
    return 1;
  }
  if (dr == 1 && dc == 0) {
    return 3;
  }
  throw std::invalid_argument("relay_leg: cells are not adjacent");
}

/*!
 * @brief One channel of a two-site operator.
 *
 * u has legs (in, out, kappa) and vt (kappa, in, out); kappa has dimension 1
 * and a definite parity; the singular value is folded into u.
 */
template <class tensor>
struct relay_channel {
  ftensor<tensor> u;
  ftensor<tensor> vt;
};

/*!
 * @brief Graded SVD of a wrap_twosite_gate()-loaded operator (legs in1, in2,
 *        out1, out2), split into dimension-1 channels.
 *
 * Channels with a zero singular value may be omitted.
 */
template <class tensor>
std::vector<relay_channel<tensor>> relay_channels(const ftensor<tensor>& op12) {
  if (op12.rank() != 4) {
    throw std::invalid_argument("relay_channels expects a rank-4 operator");
  }
  ftensor<tensor> u, vt;
  std::vector<double> s;
  const int info =
      svd(op12, mptensor::Axes(0, 2), mptensor::Axes(1, 3), u, s, vt);
  if (info != 0) {
    throw std::invalid_argument("relay_channels: operator SVD failed");
  }
  u.multiply_vector(s, 2);

  std::vector<relay_channel<tensor>> channels;
  for (std::size_t k = 0; k < s.size(); ++k) {
    if (s[k] == 0.0) {
      continue;
    }
    channels.push_back({slice(u, 2, k, k + 1), slice(vt, 0, k, k + 1)});
  }
  return channels;
}

/*!
 * @brief One folded site of the relay network.
 *
 * @return Plain rank-6 tensor ([l lb], [t tb], [r rb], [b bb], s_ket,
 *         s_bra).
 * @param[in] Tn Wrapped rank-5 center tensor.
 * @param[in] role source: entry_leg == -1; target: exit_leg == -1; middle:
 *            both in 0..3 and distinct.
 * @param[in] entry_leg Leg through which the string enters.
 * @param[in] exit_leg Leg through which the string leaves; the crossing sign
 *            is applied on it (source, middle). If this exit is l or t, an
 *            odd kappa also contributes the bond-orientation supertrace sign.
 * @param[in] channel The channel whose kappa is relayed.
 * @throw std::invalid_argument on an inconsistent role / leg combination.
 */
template <class tensor>
tensor build_relay_site(const ftensor<tensor>& Tn, relay_role role,
                        int entry_leg, int exit_leg,
                        const relay_channel<tensor>& channel) {
  if (Tn.rank() != 5 || channel.u.rank() != 3 || channel.vt.rank() != 3) {
    throw std::invalid_argument(
        "build_relay_site expects a rank-5 site and rank-3 channel factors");
  }
  if (channel.u.shape()[2] != 1 || channel.vt.shape()[0] != 1) {
    throw std::invalid_argument(
        "build_relay_site expects dimension-1 kappa channel factors");
  }
  if (channel.u.parity[2] != channel.vt.parity[0]) {
    throw std::invalid_argument(
        "build_relay_site expects matching kappa parities");
  }
  const bool kappa_parity = channel.u.parity[2][0];

  ftensor<tensor> ket;
  mptensor::Axes axes;
  std::vector<bool> has_kappa(4, false);
  int crossing_axis = -1;
  parity_vector crossing_parity;
  switch (role) {
    case relay_role::source:
      if (entry_leg != -1 || !detail::valid_relay_leg(exit_leg)) {
        throw std::invalid_argument("build_relay_site: invalid source legs");
      }
      ket = tensordot(Tn, channel.u, mptensor::Axes(4), mptensor::Axes(0));
      has_kappa[exit_leg] = true;
      for (int leg = 0; leg < 4; ++leg) {
        detail::relay_push_virtual_with_kappa(axes, leg, -1, -1, exit_leg, 5);
      }
      axes.push(4);
      crossing_axis = exit_leg;
      crossing_parity = Tn.parity[exit_leg];
      break;
    case relay_role::target:
      if (!detail::valid_relay_leg(entry_leg) || exit_leg != -1) {
        throw std::invalid_argument("build_relay_site: invalid target legs");
      }
      ket = tensordot(Tn, channel.vt, mptensor::Axes(4), mptensor::Axes(1));
      has_kappa[entry_leg] = true;
      for (int leg = 0; leg < 4; ++leg) {
        detail::relay_push_virtual_with_kappa(axes, leg, entry_leg, 4, -1, -1);
      }
      axes.push(5);
      break;
    case relay_role::middle:
      if (!detail::valid_relay_leg(entry_leg) ||
          !detail::valid_relay_leg(exit_leg) || entry_leg == exit_leg) {
        throw std::invalid_argument("build_relay_site: invalid middle legs");
      }
      ket = tensordot(
          Tn,
          detail::relay_identity_string<tensor>(Tn.get_comm(), kappa_parity),
          mptensor::Axes(), mptensor::Axes());
      has_kappa[entry_leg] = true;
      has_kappa[exit_leg] = true;
      for (int leg = 0; leg < 4; ++leg) {
        detail::relay_push_virtual_with_kappa(axes, leg, entry_leg, 5, exit_leg,
                                              6);
      }
      axes.push(4);
      crossing_axis = exit_leg;
      crossing_parity = Tn.parity[exit_leg];
      break;
    default:
      throw std::invalid_argument("build_relay_site: invalid role");
  }

  ftensor<tensor> fused = detail::fuse_relay_kappa(ket, axes, has_kappa);
  if (crossing_axis >= 0 && exit_leg <= 1 && kappa_parity) {
    // Design §3.3: bond-orientation supertrace sign for l/t exits.
    detail::multiply_ftensor_scalar(fused, typename tensor::value_type(-1));
  }
  if (crossing_axis >= 0) {
    // Design §3.3: crossing sign for kappa moving past the exit bond.
    detail::apply_relay_crossing_mask(fused, crossing_axis, crossing_parity,
                                      kappa_parity);
  }
  return detail::doubled_pipeline(Tn, fused);
}

/*!
 * @brief Folded window for one channel.
 *
 * @return [row][col] rank-6 tensors. Cells on the path come from
 *         build_relay_site(); the others are detail::doubled_pipeline(Tn,
 *         Tn).
 */
template <class tensor>
std::vector<std::vector<tensor>> build_relay_window(
    const std::vector<std::vector<ftensor<tensor>>>& Tn, window_cell source,
    window_cell target, const relay_channel<tensor>& channel,
    relay_order order = relay_order::x_first) {
  if (Tn.empty() || Tn[0].empty()) {
    throw std::invalid_argument("build_relay_window: empty window");
  }
  const std::size_t ncol = Tn[0].size();
  for (const auto& row : Tn) {
    if (row.size() != ncol) {
      throw std::invalid_argument("build_relay_window: ragged window");
    }
  }
  const auto in_window = [&](window_cell c) {
    return 0 <= c.row && static_cast<std::size_t>(c.row) < Tn.size() &&
           0 <= c.col && static_cast<std::size_t>(c.col) < ncol;
  };
  if (!in_window(source) || !in_window(target)) {
    throw std::invalid_argument("build_relay_window: endpoint out of range");
  }

  std::vector<std::vector<tensor>> window(Tn.size());
  for (std::size_t row = 0; row < Tn.size(); ++row) {
    window[row].resize(ncol);
    for (std::size_t col = 0; col < ncol; ++col) {
      window[row][col] = detail::doubled_pipeline(Tn[row][col], Tn[row][col]);
    }
  }

  const std::vector<window_cell> path = relay_path(source, target, order);
  for (std::size_t i = 0; i < path.size(); ++i) {
    const window_cell c = path[i];
    if (!in_window(c)) {
      throw std::invalid_argument("build_relay_window: path out of range");
    }
    relay_role role = relay_role::middle;
    int entry = -1;
    int exit = -1;
    if (i == 0) {
      role = relay_role::source;
      exit = relay_leg(path[i], path[i + 1]);
    } else if (i + 1 == path.size()) {
      role = relay_role::target;
      entry = relay_leg(path[i], path[i - 1]);
    } else {
      entry = relay_leg(path[i], path[i - 1]);
      exit = relay_leg(path[i], path[i + 1]);
    }
    window[c.row][c.col] =
        build_relay_site(Tn[c.row][c.col], role, entry, exit, channel);
  }
  return window;
}

}  // namespace tenes::fermion

#endif  // TENES_SRC_FERMION_RELAY_HPP_
