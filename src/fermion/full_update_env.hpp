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
 *  @brief Open-channel folded two-site environment for fermionic full update.
 */

#ifndef TENES_SRC_FERMION_FULL_UPDATE_ENV_HPP_
#define TENES_SRC_FERMION_FULL_UPDATE_ENV_HPP_

#include <algorithm>
#include <cmath>
#include <complex>
#include <sstream>
#include <stdexcept>
#include <type_traits>

#include "reduced_measure.hpp"

namespace tenes::fermion {

template <class tensor>
struct full_update_environment {
  ftensor<tensor> N;
  tensor N_plain;
  double forbidden_ratio = 0.0;
  std::complex<double> phase = {1.0, 0.0};
};

namespace detail {

template <class tensor>
ftensor<tensor> insert_even_dummy_leg(const ftensor<tensor>& Q, int slot) {
  if (Q.rank() != 4) {
    throw std::runtime_error("build_full_update_environment expects rank-4 Q");
  }
  const mptensor::Shape q_shape = Q.t.shape();
  mptensor::Shape shape;
  leg_parities parity;
  for (int ax = 0; ax < 5; ++ax) {
    if (ax == slot) {
      shape.push(1);
      parity.push_back(parity_vector{false});
    } else {
      const int src = ax < slot ? ax : ax - 1;
      shape.push(q_shape[src]);
      parity.push_back(Q.parity[src]);
    }
  }
  return ftensor<tensor>{mptensor::reshape(Q.t, shape), parity};
}

template <class tensor>
ftensor<tensor> make_left_identity_factor(const parity_vector& p) {
  const std::size_t n = p.size();
  ftensor<tensor> I{tensor(mptensor::Shape(n, n, n, n)), {p, p, p, p}};
  for (std::size_t in = 0; in < n; ++in) {
    for (std::size_t out = 0; out < n; ++out) {
      I.t.set_value(mptensor::Index(in, out, in, out),
                    typename tensor::value_type(1.0));
    }
  }
  return reshape(I, mptensor::Shape(n, n, n * n));
}

template <class tensor>
ftensor<tensor> make_right_identity_factor(const parity_vector& p) {
  const std::size_t n = p.size();
  ftensor<tensor> I{tensor(mptensor::Shape(n, n, n, n)), {p, p, p, p}};
  for (std::size_t in = 0; in < n; ++in) {
    for (std::size_t out = 0; out < n; ++out) {
      I.t.set_value(mptensor::Index(in, out, in, out),
                    typename tensor::value_type(1.0));
    }
  }
  return reshape(I, mptensor::Shape(n * n, n, n));
}

template <class tensor>
void project_even(ftensor<tensor>& a) {
  mptensor::Index idx;
  idx.resize(a.rank());
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    a.t.global_index_fast(n, idx);
    if (count_odd(a.parity, idx) % 2 == 1) {
      a.t[n] = typename tensor::value_type(0.0);
    }
  }
}

template <class tensor>
void apply_input_pair_mask(tensor& a, const parity_vector& p0,
                           const parity_vector& p1) {
  mptensor::Index idx;
  idx.resize(a.shape().size());
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    a.global_index_fast(n, idx);
    if (p0[idx[0]] && p1[idx[1]]) {
      a[n] = -a[n];
    }
  }
}

}  // namespace detail

template <class tensor>
full_update_environment<tensor> build_full_update_environment(
    const tensor& C1, const tensor& C2, const tensor& C3, const tensor& C4,
    const tensor& eT1, const tensor& eT2, const tensor& eT3, const tensor& eT4,
    const tensor& eT5, const tensor& eT6, const ftensor<tensor>& QA,
    const ftensor<tensor>& QB, reduced_pair_direction direction,
    double forbidden_tol = 1.0e-8) {
  if (direction != reduced_pair_direction::horizontal &&
      direction != reduced_pair_direction::vertical) {
    throw std::runtime_error("build_full_update_environment: invalid direction");
  }
  if (QA.rank() != 4 || QB.rank() != 4) {
    throw std::runtime_error(
        "build_full_update_environment expects rank-4 QR factors");
  }

  ftensor<tensor> QAp;
  ftensor<tensor> QBp;
  if (direction == reduced_pair_direction::horizontal) {
    QAp = detail::insert_even_dummy_leg(QA, 2);
    QBp = detail::insert_even_dummy_leg(QB, 0);
  } else {
    QAp = detail::insert_even_dummy_leg(QA, 3);
    QBp = detail::insert_even_dummy_leg(QB, 1);
  }

  const parity_vector& pA = QA.parity[3];
  const parity_vector& pB = QB.parity[3];
  const std::size_t nA = pA.size();
  const std::size_t nB = pB.size();
  const ftensor<tensor> u = detail::make_left_identity_factor<tensor>(pA);
  const ftensor<tensor> vt = detail::make_right_identity_factor<tensor>(pB);

  const reduced_pair_halves<tensor> halves =
      build_reduced_pair_halves_from_factors(QAp, QBp, u, vt, direction);
  tensor left, right;
  absorb_reduced_pair_halves(C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6,
                             halves, left, right);

  tensor M = mptensor::tensordot(left, right, mptensor::Axes(0, 2),
                                 mptensor::Axes(0, 2));
  M = mptensor::reshape(M, mptensor::Shape(nA, nA, nB, nB));
  ftensor<tensor> Ntilde{M, {pA, pA, pB, pB}};
  const double forbidden_abs = parity_violation(Ntilde);
  const double scale = max_abs(Ntilde);
  const double forbidden_ratio = scale > 0.0 ? forbidden_abs / scale : 0.0;
  if (forbidden_ratio > forbidden_tol) {
    std::stringstream ss;
    ss << "build_full_update_environment: forbidden parity block ratio "
       << forbidden_ratio << " exceeds " << forbidden_tol
       << " (forbidden max_abs=" << forbidden_abs << ", max_abs=" << scale
       << "); the CTM may not have converged (check iteration_max and "
          "CTM_Convergence_Epsilon)";
    throw std::runtime_error(ss.str());
  }
  detail::project_even(Ntilde);

  full_update_environment<tensor> result;
  result.N = transpose(Ntilde, mptensor::Axes(0, 2, 1, 3));

  tensor identity(mptensor::Shape(nA, nB, nA, nB));
  for (std::size_t a = 0; a < nA; ++a) {
    for (std::size_t b = 0; b < nB; ++b) {
      identity.set_value(mptensor::Index(a, b, a, b),
                         typename tensor::value_type(1.0));
    }
  }
  const ftensor<tensor> identity_wrap = wrap_twosite_gate(identity, pA, pB);

  double local_finite = 1.0;
  for (std::size_t local = 0; local < result.N.t.local_size(); ++local) {
    const auto value = result.N.t[local];
    if (!std::isfinite(std::real(value)) || !std::isfinite(std::imag(value))) {
      local_finite = 0.0;
    }
  }
  std::vector<double> finite_collective{local_finite};
  tenes::allreduce_min(finite_collective, result.N.t.get_comm());
  const double N_max_abs = tenes::fermion::max_abs(result.N);
  const mptensor::Axes all4(0, 1, 2, 3);
  const std::complex<double> norm =
      mptensor::trace(result.N.t, identity_wrap.t, all4, all4);
  const double norm_abs = std::abs(norm);
  if (finite_collective[0] == 0.0 || !std::isfinite(N_max_abs) ||
      !std::isfinite(norm.real()) || !std::isfinite(norm.imag()) ||
      norm_abs <= 1.0e-12 * N_max_abs * nA * nB) {
    std::stringstream ss;
    ss << "build_full_update_environment: invalid window norm " << norm
       << " (max_abs=" << N_max_abs << ", nA=" << nA << ", nB=" << nB
       << ", elements_finite=" << finite_collective[0] << ")";
    throw std::runtime_error(ss.str());
  }
  const std::complex<double> phase = norm / norm_abs;
  if (std::abs(phase - std::complex<double>(1.0, 0.0)) > 1.0e-14) {
    for (std::size_t local = 0; local < result.N.t.local_size(); ++local) {
      if constexpr (std::is_same_v<typename tensor::value_type,
                                   std::complex<double>>) {
        result.N.t[local] *= std::conj(phase);
      } else {
        result.N.t[local] *= phase.real();
      }
    }
    result.phase = phase;
  }
  result.N_plain = result.N.t;
  detail::apply_input_pair_mask(result.N_plain, result.N.parity[0],
                                result.N.parity[1]);
  result.forbidden_ratio = forbidden_ratio;
  return result;
}

}  // namespace tenes::fermion

#endif  // TENES_SRC_FERMION_FULL_UPDATE_ENV_HPP_
