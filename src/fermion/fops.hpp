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

#ifndef TENES_SRC_FERMION_FOPS_HPP_
#define TENES_SRC_FERMION_FOPS_HPP_

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "ftensor.hpp"

namespace tenes {
namespace fermion {
namespace detail {

inline bool contains_axis(const mptensor::Axes& axes, std::size_t ax) {
  for (std::size_t i = 0; i < axes.size(); ++i) {
    if (axes[i] == ax) {
      return true;
    }
  }
  return false;
}

inline mptensor::Axes tensordot_left_perm(std::size_t rank,
                                          const mptensor::Axes& axes) {
  mptensor::Axes perm;
  for (std::size_t ax = 0; ax < rank; ++ax) {
    if (!contains_axis(axes, ax)) {
      perm.push(ax);
    }
  }
  for (std::size_t i = 0; i < axes.size(); ++i) {
    perm.push(axes[i]);
  }
  return perm;
}

inline mptensor::Axes tensordot_right_perm(std::size_t rank,
                                           const mptensor::Axes& axes) {
  mptensor::Axes perm;
  for (std::size_t i = axes.size(); i > 0; --i) {
    perm.push(axes[i - 1]);
  }
  for (std::size_t ax = 0; ax < rank; ++ax) {
    if (!contains_axis(axes, ax)) {
      perm.push(ax);
    }
  }
  return perm;
}

template <class tensor>
ftensor<tensor> apply_transpose_sign_mask(const ftensor<tensor>& a,
                                          const mptensor::Axes& axes) {
  ftensor<tensor> ret = a;
  for (std::size_t n = 0; n < ret.t.local_size(); ++n) {
    auto idx = ret.t.global_index(n);
    if (transpose_sign(ret.parity, idx, axes) < 0) {
      typename tensor::value_type v;
      ret.t.get_value(idx, v);
      ret.t.set_value(idx, -v);
    }
  }
  return ret;
}

inline void validate_contracted_parity(const leg_parities& parity_a,
                                       const leg_parities& parity_b,
                                       const mptensor::Axes& axes_a,
                                       const mptensor::Axes& axes_b) {
  if (axes_a.size() != axes_b.size()) {
    throw std::runtime_error(
        "fermion tensordot parity check: axis count mismatch");
  }
  for (std::size_t i = 0; i < axes_a.size(); ++i) {
    if (parity_a[axes_a[i]] != parity_b[axes_b[i]]) {
      throw std::runtime_error(
          "fermion tensordot parity check: contracted parity mismatch");
    }
  }
}

inline leg_parities free_leg_parities(const leg_parities& parity,
                                      const mptensor::Axes& contracted) {
  leg_parities ret;
  for (std::size_t ax = 0; ax < parity.size(); ++ax) {
    if (!contains_axis(contracted, ax)) {
      ret.push_back(parity[ax]);
    }
  }
  return ret;
}

}  // namespace detail

template <class tensor>
void apply_swap(ftensor<tensor>& a, int ax1, int ax2) {
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    auto idx = a.t.global_index(n);
    if (a.parity[ax1][idx[ax1]] && a.parity[ax2][idx[ax2]]) {
      typename tensor::value_type v;
      a.t.get_value(idx, v);
      a.t.set_value(idx, -v);
    }
  }
}

template <class tensor>
void apply_parity(ftensor<tensor>& a, int ax) {
  std::vector<double> sign(a.parity[ax].size());
  for (std::size_t i = 0; i < sign.size(); ++i) {
    sign[i] = a.parity[ax][i] ? -1.0 : 1.0;
  }
  a.t.multiply_vector(sign, ax);
}

template <class tensor>
double parity_violation(const ftensor<tensor>& a) {
  double v = 0.0;
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    auto idx = a.t.global_index(n);
    if (count_odd(a.parity, idx) % 2 == 1) {
      typename tensor::value_type x;
      a.t.get_value(idx, x);
      v = std::max(v, std::abs(x));
    }
  }
  return v;
}

template <class tensor>
ftensor<tensor> transpose(const ftensor<tensor>& a,
                          const mptensor::Axes& axes) {
  ftensor<tensor> ret = a;
  ret.transpose(axes);
  return ret;
}

template <class tensor>
ftensor<tensor> tensordot(const ftensor<tensor>& a, const ftensor<tensor>& b,
                          const mptensor::Axes& axes_a,
                          const mptensor::Axes& axes_b) {
  detail::validate_contracted_parity(a.parity, b.parity, axes_a, axes_b);
  ftensor<tensor> a_masked = detail::apply_transpose_sign_mask(
      a, detail::tensordot_left_perm(a.parity.size(), axes_a));
  ftensor<tensor> b_masked = detail::apply_transpose_sign_mask(
      b, detail::tensordot_right_perm(b.parity.size(), axes_b));
  ftensor<tensor> ret;
  ret.t = mptensor::tensordot(a_masked.t, b_masked.t, axes_a, axes_b);
  ret.parity = detail::free_leg_parities(a.parity, axes_a);
  leg_parities b_free = detail::free_leg_parities(b.parity, axes_b);
  ret.parity.insert(ret.parity.end(), b_free.begin(), b_free.end());
  return ret;
}

template <class tensor>
typename tensor::value_type trace(const ftensor<tensor>& a,
                                  const ftensor<tensor>& b,
                                  const mptensor::Axes& axes_a,
                                  const mptensor::Axes& axes_b) {
  detail::validate_contracted_parity(a.parity, b.parity, axes_a, axes_b);
  ftensor<tensor> a_masked = detail::apply_transpose_sign_mask(
      a, detail::tensordot_left_perm(a.parity.size(), axes_a));
  ftensor<tensor> b_masked = detail::apply_transpose_sign_mask(
      b, detail::tensordot_right_perm(b.parity.size(), axes_b));
  return mptensor::trace(a_masked.t, b_masked.t, axes_a, axes_b);
}

}  // namespace fermion
}  // namespace tenes

#endif  // TENES_SRC_FERMION_FOPS_HPP_
