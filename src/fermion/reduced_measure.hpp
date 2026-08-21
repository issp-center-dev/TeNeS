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

#ifndef TENES_SRC_FERMION_REDUCED_MEASURE_HPP_
#define TENES_SRC_FERMION_REDUCED_MEASURE_HPP_

#include <stdexcept>
#include <vector>

#include "fermion_info.hpp"
#include "reduced.hpp"

namespace tenes::fermion {

template <class tensor>
tensor lambda_dressed_tensor(const tensor& Tn,
                             const std::vector<std::vector<double>>& lambda) {
  tensor dressed = Tn;
  for (int leg = 0; leg < 4; ++leg) {
    dressed.multiply_vector(lambda[leg], leg);
  }
  return dressed;
}

template <class tensor>
std::vector<tensor> lambda_dressed_tensors(
    const std::vector<tensor>& Tn,
    const std::vector<std::vector<std::vector<double>>>& lambda) {
  std::vector<tensor> dressed;
  dressed.reserve(Tn.size());
  for (int site = 0; site < static_cast<int>(Tn.size()); ++site) {
    dressed.push_back(lambda_dressed_tensor(Tn[site], lambda[site]));
  }
  return dressed;
}

template <class tensor>
std::vector<tensor> build_reduced_density_tensors(const std::vector<tensor>& Tn,
                                                  const FermionInfo& finfo) {
  std::vector<tensor> reduced;
  reduced.reserve(Tn.size());
  for (int site = 0; site < static_cast<int>(Tn.size()); ++site) {
    reduced.push_back(build_reduced_op(wrap_Tn(Tn[site], finfo, site)));
  }
  return reduced;
}

template <class tensor>
tensor build_reduced_identity_pair(const ftensor<tensor>& TnA,
                                   const ftensor<tensor>& TnB,
                                   reduced_pair_direction direction) {
  if (direction == reduced_pair_direction::horizontal) {
    return mptensor::tensordot(build_reduced(TnA), build_reduced(TnB),
                               mptensor::Axes(2), mptensor::Axes(0));
  }
  return mptensor::tensordot(build_reduced(TnA), build_reduced(TnB),
                             mptensor::Axes(3), mptensor::Axes(1));
}

namespace detail {

inline mptensor::Axes all_axes(int rank) {
  mptensor::Axes axes;
  for (int ax = 0; ax < rank; ++ax) {
    axes.push(ax);
  }
  return axes;
}

}  // namespace detail

// Mean-field norm of a two-site pair state: <pair|pair> as a graded full
// contraction. No environment tensors: the lambda weights on the open legs
// are expected to be multiplied into TnA / TnB beforehand (the same dressing
// the bosonic mean-field path applies in measure_twosite).
template <class tensor>
typename tensor::value_type contract_pair_MF(const ftensor<tensor>& pair) {
  const mptensor::Axes axes = detail::all_axes(pair.rank());
  return trace(conj(pair), pair, axes, axes);
}

// Mean-field expectation value (unnormalized) of op12 on a pair state.
// op12 must be loaded with wrap_twosite_gate (input-leg swap only): this is
// the single-layer convention pinned by the Fock-verified direct path
// (r2_expect_two in the R5 test and the mf_* oracle cases), not the blob
// convention (wrap_reduced_pair_op).
template <class tensor>
typename tensor::value_type contract_pair_MF(const ftensor<tensor>& pair,
                                             const ftensor<tensor>& op12) {
  const mptensor::Axes axes = detail::all_axes(pair.rank());
  return trace(conj(pair), apply_pair_op(pair, op12), axes, axes);
}

namespace detail {

template <class tensor>
typename tensor::value_type trace_boundary_pairs(const tensor& a) {
  if (a.rank() != 4) {
    throw std::runtime_error("expected four boundary legs");
  }
  typename tensor::value_type value = 0.0;
  mptensor::Index idx;
  idx.resize(4);
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    a.global_index_fast(n, idx);
    if (idx[0] == idx[1] && idx[2] == idx[3]) {
      value += a[n];
    }
  }
  return value;
}

}  // namespace detail

template <class tensor>
typename tensor::value_type contract_reduced_pair_horizontal_density_CTM(
    const tensor& C1, const tensor& C2, const tensor& C3, const tensor& C4,
    const tensor& eT1, const tensor& eT2, const tensor& eT3, const tensor& eT4,
    const tensor& eT5, const tensor& eT6, const tensor& blob) {
  using mptensor::Axes;
  const tensor left_lower = tensordot(
      eT5,
      tensordot(C1, tensordot(C4, eT6, Axes(1), Axes(0)), Axes(0), Axes(1)),
      Axes(1), Axes(1));
  tensor work = tensordot(blob, left_lower, Axes(0, 2), Axes(3, 1));
  work = tensordot(eT1, work, Axes(0, 2), Axes(5, 0));

  const tensor right_lower = tensordot(
      eT4,
      tensordot(C2, tensordot(C3, eT3, Axes(0), Axes(1)), Axes(1), Axes(1)),
      Axes(0), Axes(1));
  work = tensordot(work, right_lower, Axes(2, 3), Axes(3, 1));
  work = tensordot(eT2, work, Axes(1, 2), Axes(4, 1));
  return detail::trace_boundary_pairs(work);
}

template <class tensor>
typename tensor::value_type contract_reduced_pair_vertical_density_CTM(
    const tensor& C1, const tensor& C2, const tensor& C3, const tensor& C4,
    const tensor& eT1, const tensor& eT2, const tensor& eT3, const tensor& eT4,
    const tensor& eT5, const tensor& eT6, const tensor& blob) {
  using mptensor::Axes;
  const tensor top_left = tensordot(
      eT6,
      tensordot(C1, tensordot(C2, eT1, Axes(0), Axes(1)), Axes(1), Axes(1)),
      Axes(1), Axes(0));
  tensor work = tensordot(blob, top_left, Axes(0, 1), Axes(1, 3));
  work = tensordot(eT2, work, Axes(0, 2), Axes(5, 0));

  const tensor bottom_right = tensordot(
      eT5,
      tensordot(C3, tensordot(C4, eT4, Axes(0), Axes(1)), Axes(1), Axes(1)),
      Axes(0), Axes(1));
  work = tensordot(work, bottom_right, Axes(1, 3), Axes(1, 3));
  work = tensordot(eT3, work, Axes(1, 2), Axes(4, 1));
  return detail::trace_boundary_pairs(work);
}

template <class tensor>
typename tensor::value_type contract_reduced_pair_density_CTM(
    const tensor& C1, const tensor& C2, const tensor& C3, const tensor& C4,
    const tensor& eT1, const tensor& eT2, const tensor& eT3, const tensor& eT4,
    const tensor& eT5, const tensor& eT6, const tensor& blob,
    reduced_pair_direction direction) {
  if (direction == reduced_pair_direction::horizontal) {
    return contract_reduced_pair_horizontal_density_CTM(
        C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6, blob);
  }
  return contract_reduced_pair_vertical_density_CTM(C1, C2, C3, C4, eT1, eT2,
                                                    eT3, eT4, eT5, eT6, blob);
}

}  // namespace tenes::fermion

#endif  // TENES_SRC_FERMION_REDUCED_MEASURE_HPP_
