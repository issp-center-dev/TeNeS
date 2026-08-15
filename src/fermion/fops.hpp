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
#include <complex>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
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

inline int fermion_svd_log_limit() {
  const char* raw = std::getenv("TENES_FERMION_SVD_LOG_LIMIT");
  if (raw == nullptr) {
    return 0;
  }
  return std::max(0, std::atoi(raw));
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

inline parity_vector fuse_axes(const leg_parities& parity,
                               const mptensor::Axes& axes) {
  if (axes.size() == 0) {
    return parity_vector{false};
  }
  parity_vector ret = parity[axes[0]];
  for (std::size_t i = 1; i < axes.size(); ++i) {
    ret = fuse(ret, parity[axes[i]]);
  }
  return ret;
}

inline mptensor::Shape shape_from_axes(const mptensor::Shape& shape,
                                       const mptensor::Axes& axes) {
  mptensor::Shape ret;
  for (std::size_t i = 0; i < axes.size(); ++i) {
    ret.push(shape[axes[i]]);
  }
  return ret;
}

inline std::size_t product_shape(const mptensor::Shape& shape) {
  std::size_t ret = 1;
  for (std::size_t i = 0; i < shape.size(); ++i) {
    ret *= shape[i];
  }
  return ret;
}

inline std::size_t count_even(const parity_vector& parity) {
  std::size_t ret = 0;
  for (bool p : parity) {
    if (!p) {
      ++ret;
    }
  }
  return ret;
}

inline double scalar_conj(double v) { return v; }

inline std::complex<double> scalar_conj(std::complex<double> v) {
  return std::conj(v);
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

template <class tensor>
ftensor<tensor> conj(const ftensor<tensor>& a) {
  ftensor<tensor> ret = a;
  for (std::size_t n = 0; n < ret.t.local_size(); ++n) {
    auto idx = ret.t.global_index(n);
    const int m = count_odd(ret.parity, idx);
    const int sign = ((m * (m - 1) / 2) % 2 == 0) ? 1 : -1;
    typename tensor::value_type v;
    ret.t.get_value(idx, v);
    ret.t.set_value(idx, static_cast<double>(sign) * detail::scalar_conj(v));
  }
  return ret;
}

template <class tensor>
ftensor<tensor> slice(const ftensor<tensor>& a, std::size_t ax, std::size_t b,
                      std::size_t e) {
  ftensor<tensor> ret;
  ret.t = mptensor::slice(a.t, ax, b, e);
  ret.parity = a.parity;
  ret.parity[ax] =
      parity_vector(a.parity[ax].begin() + b, a.parity[ax].begin() + e);
  return ret;
}

template <class tensor>
ftensor<tensor> extend(const ftensor<tensor>& a, const mptensor::Shape& sh) {
  ftensor<tensor> ret;
  ret.t = mptensor::extend(a.t, sh);
  ret.parity = a.parity;
  for (std::size_t ax = 0; ax < ret.parity.size(); ++ax) {
    if (ret.parity[ax].size() > sh[ax]) {
      throw std::runtime_error("fermion extend: new shape is smaller");
    }
    ret.parity[ax].resize(sh[ax], false);
  }
  return ret;
}

template <class tensor>
ftensor<tensor> reshape(const ftensor<tensor>& a, const mptensor::Shape& sh) {
  ftensor<tensor> ret;
  ret.t = mptensor::reshape(a.t, sh);
  ret.parity.clear();
  std::size_t old_ax = 0;
  for (std::size_t new_ax = 0; new_ax < sh.size(); ++new_ax) {
    if (old_ax >= a.parity.size()) {
      throw std::runtime_error("fermion reshape: too many output axes");
    }
    std::size_t dim = a.parity[old_ax].size();
    parity_vector fused = a.parity[old_ax];
    ++old_ax;
    while (dim < sh[new_ax] && old_ax < a.parity.size()) {
      fused = fuse(fused, a.parity[old_ax]);
      dim *= a.parity[old_ax].size();
      ++old_ax;
    }
    if (dim != sh[new_ax]) {
      throw std::runtime_error(
          "fermion reshape: only adjacent leg fusion is supported");
    }
    ret.parity.push_back(fused);
  }
  if (old_ax != a.parity.size()) {
    throw std::runtime_error("fermion reshape: unused input axes");
  }
  return ret;
}

template <class tensor>
double max_abs(const ftensor<tensor>& a) {
  return mptensor::max_abs(a.t);
}

template <class tensor>
int svd(const ftensor<tensor>& a, std::vector<double>& s) {
  return mptensor::svd(a.t, s);
}

inline std::vector<double>& parity_cleanup_observations() {
  static std::vector<double> values;
  return values;
}

template <class tensor>
double verify_and_clean_even_parity(ftensor<tensor>& a, const char* context) {
  const double v = parity_violation(a);
  const double scale = std::max(1.0, max_abs(a));
  const double threshold = 1.0e-10 * scale;
  parity_cleanup_observations().push_back(v);
  if (v > threshold) {
    std::stringstream ss;
    ss << context << " produced odd-parity elements: max_abs=" << v
       << " threshold=" << threshold;
    throw std::runtime_error(ss.str());
  }
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    const auto index = a.t.global_index(n);
    if (count_odd(a.parity, index) % 2 == 1) {
      a.t.set_value(index, typename tensor::value_type{});
    }
  }
  return v;
}

template <class tensor>
tensor make_perm_matrix(const std::vector<std::size_t>& perm) {
  tensor ret(mptensor::Shape(perm.size(), perm.size()));
  for (std::size_t n = 0; n < ret.local_size(); ++n) {
    auto idx = ret.global_index(n);
    ret.set_value(idx, perm[idx[0]] == idx[1] ? 1.0 : 0.0);
  }
  return ret;
}

template <class tensor>
int qr(const ftensor<tensor>& a, const mptensor::Axes& rows,
       const mptensor::Axes& cols, ftensor<tensor>& q, ftensor<tensor>& r) {
  ftensor<tensor> a_ordered = transpose(a, rows + cols);
  mptensor::Shape row_shape = detail::shape_from_axes(a.shape(), rows);
  mptensor::Shape col_shape = detail::shape_from_axes(a.shape(), cols);
  const std::size_t drow = detail::product_shape(row_shape);
  const std::size_t dcol = detail::product_shape(col_shape);
  typename tensor::comm_type comm = a.t.get_comm();
  tensor mat = mptensor::reshape(a_ordered.t, mptensor::Shape(drow, dcol));

  parity_vector row_parity = detail::fuse_axes(a.parity, rows);
  parity_vector col_parity = detail::fuse_axes(a.parity, cols);
  std::vector<std::size_t> row_perm = parity_sort_perm(row_parity);
  std::vector<std::size_t> col_perm = parity_sort_perm(col_parity);
  tensor prow = make_perm_matrix<tensor>(row_perm);
  tensor pcol = make_perm_matrix<tensor>(col_perm);
  tensor sorted =
      mptensor::tensordot(prow, mat, mptensor::Axes(1), mptensor::Axes(0));
  sorted =
      mptensor::tensordot(sorted, pcol, mptensor::Axes(1), mptensor::Axes(1));

  const std::size_t row_even = detail::count_even(row_parity);
  const std::size_t col_even = detail::count_even(col_parity);
  const std::size_t row_odd = drow - row_even;
  const std::size_t col_odd = dcol - col_even;
  const std::size_t size_even = std::min(row_even, col_even);
  const std::size_t size_odd = std::min(row_odd, col_odd);
  const std::size_t size = size_even + size_odd;

  tensor q_sorted(comm, mptensor::Shape(drow, size));
  tensor r_sorted(comm, mptensor::Shape(size, dcol));
  int info = 0;
  if (size_even > 0) {
    tensor even_block = mptensor::slice(sorted, mptensor::Index(0, 0),
                                        mptensor::Index(row_even, col_even));
    tensor qe, re;
    info =
        mptensor::qr(even_block, mptensor::Axes(0), mptensor::Axes(1), qe, re);
    q_sorted.set_slice(qe, mptensor::Index(0, 0),
                       mptensor::Index(row_even, size_even));
    r_sorted.set_slice(re, mptensor::Index(0, 0),
                       mptensor::Index(size_even, col_even));
  }
  if (size_odd > 0) {
    tensor odd_block =
        mptensor::slice(sorted, mptensor::Index(row_even, col_even),
                        mptensor::Index(drow, dcol));
    tensor qo, ro;
    int info_odd =
        mptensor::qr(odd_block, mptensor::Axes(0), mptensor::Axes(1), qo, ro);
    if (info == 0) {
      info = info_odd;
    }
    q_sorted.set_slice(qo, mptensor::Index(row_even, size_even),
                       mptensor::Index(drow, size));
    r_sorted.set_slice(ro, mptensor::Index(size_even, col_even),
                       mptensor::Index(size, dcol));
  }

  tensor q_mat =
      mptensor::tensordot(prow, q_sorted, mptensor::Axes(0), mptensor::Axes(0));
  tensor r_mat =
      mptensor::tensordot(r_sorted, pcol, mptensor::Axes(1), mptensor::Axes(0));

  mptensor::Shape q_shape = row_shape;
  q_shape.push(size);
  q.t = mptensor::reshape(q_mat, q_shape);
  mptensor::Shape r_shape;
  r_shape.push(size);
  for (std::size_t i = 0; i < col_shape.size(); ++i) {
    r_shape.push(col_shape[i]);
  }
  r.t = mptensor::reshape(r_mat, r_shape);

  parity_vector internal(size, false);
  for (std::size_t i = size_even; i < size; ++i) {
    internal[i] = true;
  }
  q.parity.clear();
  for (std::size_t i = 0; i < rows.size(); ++i) {
    q.parity.push_back(a.parity[rows[i]]);
  }
  q.parity.push_back(internal);
  r.parity.clear();
  r.parity.push_back(internal);
  for (std::size_t i = 0; i < cols.size(); ++i) {
    r.parity.push_back(a.parity[cols[i]]);
  }
  return info;
}

template <class tensor>
int svd(const ftensor<tensor>& a, const mptensor::Axes& rows,
        const mptensor::Axes& cols, ftensor<tensor>& u, std::vector<double>& s,
        ftensor<tensor>& vt) {
  ftensor<tensor> a_ordered = transpose(a, rows + cols);
  mptensor::Shape row_shape = detail::shape_from_axes(a.shape(), rows);
  mptensor::Shape col_shape = detail::shape_from_axes(a.shape(), cols);
  const std::size_t drow = detail::product_shape(row_shape);
  const std::size_t dcol = detail::product_shape(col_shape);
  typename tensor::comm_type comm = a.t.get_comm();
  tensor mat = mptensor::reshape(a_ordered.t, mptensor::Shape(drow, dcol));

  parity_vector row_parity = detail::fuse_axes(a.parity, rows);
  parity_vector col_parity = detail::fuse_axes(a.parity, cols);
  std::vector<std::size_t> row_perm = parity_sort_perm(row_parity);
  std::vector<std::size_t> col_perm = parity_sort_perm(col_parity);
  tensor prow = make_perm_matrix<tensor>(row_perm);
  tensor pcol = make_perm_matrix<tensor>(col_perm);
  tensor sorted =
      mptensor::tensordot(prow, mat, mptensor::Axes(1), mptensor::Axes(0));
  sorted =
      mptensor::tensordot(sorted, pcol, mptensor::Axes(1), mptensor::Axes(1));

  const std::size_t row_even = detail::count_even(row_parity);
  const std::size_t col_even = detail::count_even(col_parity);
  const std::size_t row_odd = drow - row_even;
  const std::size_t col_odd = dcol - col_even;
  const std::size_t size_even = std::min(row_even, col_even);
  const std::size_t size_odd = std::min(row_odd, col_odd);
  const std::size_t size = size_even + size_odd;

  tensor u_sorted(comm, mptensor::Shape(drow, size));
  tensor vt_sorted(comm, mptensor::Shape(size, dcol));
  std::vector<double> s_even, s_odd;
  int info = 0;
  if (size_even > 0) {
    tensor even_block = mptensor::slice(sorted, mptensor::Index(0, 0),
                                        mptensor::Index(row_even, col_even));
    tensor ue, vte;
    info = mptensor::svd(even_block, mptensor::Axes(0), mptensor::Axes(1), ue,
                         s_even, vte);
    u_sorted.set_slice(ue, mptensor::Index(0, 0),
                       mptensor::Index(row_even, size_even));
    vt_sorted.set_slice(vte, mptensor::Index(0, 0),
                        mptensor::Index(size_even, col_even));
  }
  if (size_odd > 0) {
    tensor odd_block =
        mptensor::slice(sorted, mptensor::Index(row_even, col_even),
                        mptensor::Index(drow, dcol));
    tensor uo, vto;
    int info_odd = mptensor::svd(odd_block, mptensor::Axes(0),
                                 mptensor::Axes(1), uo, s_odd, vto);
    if (info == 0) {
      info = info_odd;
    }
    u_sorted.set_slice(uo, mptensor::Index(row_even, size_even),
                       mptensor::Index(drow, size));
    vt_sorted.set_slice(vto, mptensor::Index(size_even, col_even),
                        mptensor::Index(size, dcol));
  }

  tensor u_mat =
      mptensor::tensordot(prow, u_sorted, mptensor::Axes(0), mptensor::Axes(0));
  tensor vt_mat = mptensor::tensordot(vt_sorted, pcol, mptensor::Axes(1),
                                      mptensor::Axes(0));

  mptensor::Shape u_shape = row_shape;
  u_shape.push(size);
  u.t = mptensor::reshape(u_mat, u_shape);
  mptensor::Shape vt_shape;
  vt_shape.push(size);
  for (std::size_t i = 0; i < col_shape.size(); ++i) {
    vt_shape.push(col_shape[i]);
  }
  vt.t = mptensor::reshape(vt_mat, vt_shape);

  s = s_even;
  s.insert(s.end(), s_odd.begin(), s_odd.end());
  parity_vector internal(size, false);
  for (std::size_t i = size_even; i < size; ++i) {
    internal[i] = true;
  }
  u.parity.clear();
  for (std::size_t i = 0; i < rows.size(); ++i) {
    u.parity.push_back(a.parity[rows[i]]);
  }
  u.parity.push_back(internal);
  vt.parity.clear();
  vt.parity.push_back(internal);
  for (std::size_t i = 0; i < cols.size(); ++i) {
    vt.parity.push_back(a.parity[cols[i]]);
  }
  return info;
}

template <class tensor>
int svd_trunc(const ftensor<tensor>& a, const mptensor::Axes& rows,
              const mptensor::Axes& cols, ftensor<tensor>& u,
              std::vector<double>& s, ftensor<tensor>& vt, int dc) {
  ftensor<tensor> full_u, full_vt;
  std::vector<double> full_s;
  int info = svd(a, rows, cols, full_u, full_s, full_vt);
  const int nkeep = std::min<int>(dc, full_s.size());
  std::vector<std::size_t> order(full_s.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    order[i] = i;
  }
  const parity_vector& full_internal = full_u.parity.back();
  static int svd_log_count = 0;
  const int svd_log_limit = detail::fermion_svd_log_limit();
  if (svd_log_count < svd_log_limit) {
    std::cerr << "TENES_FERMION_SVD call=" << svd_log_count << " dc=" << dc
              << " full_size=" << full_s.size() << " even=[";
    bool first = true;
    for (std::size_t i = 0; i < full_s.size(); ++i) {
      if (!full_internal[i]) {
        std::cerr << (first ? "" : ",") << full_s[i];
        first = false;
      }
    }
    std::cerr << "] odd=[";
    first = true;
    for (std::size_t i = 0; i < full_s.size(); ++i) {
      if (full_internal[i]) {
        std::cerr << (first ? "" : ",") << full_s[i];
        first = false;
      }
    }
    std::cerr << "]\n";
  }
  ++svd_log_count;
  std::stable_sort(order.begin(), order.end(),
                   [&](std::size_t lhs, std::size_t rhs) {
                     if (full_s[lhs] != full_s[rhs]) {
                       return full_s[lhs] > full_s[rhs];
                     }
                     if (full_internal[lhs] != full_internal[rhs]) {
                       return !full_internal[lhs];
                     }
                     return lhs < rhs;
                   });
  order.resize(nkeep);
  std::stable_sort(order.begin(), order.end(),
                   [&](std::size_t lhs, std::size_t rhs) {
                     if (full_internal[lhs] != full_internal[rhs]) {
                       return !full_internal[lhs];
                     }
                     return full_s[lhs] > full_s[rhs];
                   });

  mptensor::Shape row_shape = detail::shape_from_axes(a.shape(), rows);
  mptensor::Shape col_shape = detail::shape_from_axes(a.shape(), cols);
  const std::size_t drow = detail::product_shape(row_shape);
  const std::size_t dcol = detail::product_shape(col_shape);
  tensor full_u_mat =
      mptensor::reshape(full_u.t, mptensor::Shape(drow, full_s.size()));
  tensor full_vt_mat =
      mptensor::reshape(full_vt.t, mptensor::Shape(full_s.size(), dcol));
  tensor selector(a.t.get_comm(), mptensor::Shape(full_s.size(), nkeep));
  for (std::size_t n = 0; n < selector.local_size(); ++n) {
    auto idx = selector.global_index(n);
    selector.set_value(idx, order[idx[1]] == idx[0] ? 1.0 : 0.0);
  }
  tensor u_mat = mptensor::tensordot(full_u_mat, selector, mptensor::Axes(1),
                                     mptensor::Axes(0));
  tensor vt_mat = mptensor::tensordot(selector, full_vt_mat, mptensor::Axes(0),
                                      mptensor::Axes(0));

  mptensor::Shape u_shape = row_shape;
  u_shape.push(nkeep);
  u.t = mptensor::reshape(u_mat, u_shape);
  mptensor::Shape vt_shape;
  vt_shape.push(nkeep);
  for (std::size_t i = 0; i < col_shape.size(); ++i) {
    vt_shape.push(col_shape[i]);
  }
  vt.t = mptensor::reshape(vt_mat, vt_shape);

  s.resize(nkeep);
  parity_vector internal(nkeep, false);
  for (int i = 0; i < nkeep; ++i) {
    s[i] = full_s[order[i]];
    internal[i] = full_internal[order[i]];
  }
  u.parity.clear();
  for (std::size_t i = 0; i < rows.size(); ++i) {
    u.parity.push_back(a.parity[rows[i]]);
  }
  u.parity.push_back(internal);
  vt.parity.clear();
  vt.parity.push_back(internal);
  for (std::size_t i = 0; i < cols.size(); ++i) {
    vt.parity.push_back(a.parity[cols[i]]);
  }
  return info;
}

}  // namespace fermion
}  // namespace tenes

#endif  // TENES_SRC_FERMION_FOPS_HPP_
