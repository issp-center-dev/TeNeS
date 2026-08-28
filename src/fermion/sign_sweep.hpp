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
 *  @brief Fermionic sign machinery: swap forms and elementwise sign sweeps.
 *
 *  A single exchange of two graded legs multiplies an element by
 *  @f$(-1)^{p_a p_b}@f$, which is bilinear in the two parities. The sign of
 *  a whole leg permutation is therefore a quadratic form over
 *  @f$\mathbb{F}_2@f$: @f$(-1)^{\sum_{(x,y)\in\mathrm{Inv}(\pi)} p_x p_y}@f$
 *  where the sum runs over the inversion pairs of the permutation. The pair
 *  set does not depend on how the permutation is decomposed into adjacent
 *  transpositions, so it is a canonical form.
 *
 *  ::tenes::fermion::SwapForm holds such a pair set; toggling the same pair
 *  twice cancels it (@f$(-1)^x(-1)^x = 1@f$), which is what lets several
 *  sign-generating steps be merged and applied in one elementwise pass over
 *  the tensor (apply_sign_sweep()). The pass itself has two evaluation
 *  strategies, selected by ::tenes::fermion::SignEval: a precomputed
 *  per-parity-mask sign table, or direct evaluation of the form per element.
 */

#ifndef TENES_SRC_FERMION_SIGN_SWEEP_HPP_
#define TENES_SRC_FERMION_SIGN_SWEEP_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include "ftensor.hpp"

namespace tenes {
namespace fermion {

/*! @brief Evaluation strategy for the sign sweep.
 *
 * @c table precomputes a @f$2^{\mathrm{rank}}@f$ sign table indexed by the
 * per-leg parity mask of an element; @c direct evaluates the quadratic form
 * per element. @c automatic picks the table when the rank allows it and the
 * table is not larger than the local tensor slice
 * (detail::use_sign_table()).
 */
enum class SignEval { automatic, table, direct };

/*! @brief A set of leg pairs representing a product of swap-gate signs.
 *
 * Applying the form to an element multiplies it by @f$(-1)^{p_x p_y}@f$ for
 * every pair @f$(x, y)@f$ in the set — the canonical pair-set form of a
 * sign-generating leg reordering (see the file description). The set is the
 * exponent of a quadratic form over @f$\mathbb{F}_2@f$, so membership is
 * modulo 2: toggle() implements the symmetric difference, and merging two
 * forms by toggling all pairs of one into the other composes their signs.
 */
class SwapForm {
 public:
  /*!
   * @brief Flip membership of the pair (ax1, ax2) in the set.
   *
   * Toggling a pair already present removes it — the two sign factors
   * cancel. The pair is stored unordered (normalized to ax1 < ax2), and
   * toggling a pair with ax1 == ax2 is a no-op rather than a parity gauge:
   * use apply_parity() for @f$(-1)^{p}@f$ on a single leg.
   */
  void toggle(int ax1, int ax2) {
    if (ax1 == ax2) {
      return;
    }
    if (ax2 < ax1) {
      std::swap(ax1, ax2);
    }
    const std::pair<int, int> term(ax1, ax2);
    std::vector<std::pair<int, int>>::iterator it =
        std::lower_bound(terms_.begin(), terms_.end(), term);
    if (it != terms_.end() && *it == term) {
      terms_.erase(it);
    } else {
      terms_.insert(it, term);
    }
  }

  //! True iff the form is the identity (no sign anywhere).
  bool empty() const { return terms_.empty(); }

  //! The pair set, sorted, with each pair normalized to first < second.
  const std::vector<std::pair<int, int>>& terms() const { return terms_; }

 private:
  //! Sorted pair set; each pair normalized to first < second.
  std::vector<std::pair<int, int>> terms_;
};

/*! @brief A diagonal factor multiplied onto one leg during a sign sweep.
 *
 * Unlike a SwapForm term this is linear, not bilinear: element n picks up
 * factor[idx[axis]] regardless of the other legs. This is a general sign-sweep
 * utility retained for callers that need a per-leg diagonal factor.
 */
struct LegGauge {
  //! Leg the factor acts on.
  int axis;
  //! Per-index-value factor; size must equal the leg dimension.
  std::vector<double> factor;
};

namespace detail {
//! Largest tensor rank for which a sign table may be built (table size
//! 2^rank; beyond this the shifts would be undefined and the memory
//! impractical, so direct evaluation is used instead).
constexpr std::size_t kMaxTableRank = 24;
//! Local sizes below this are swept serially with global_index_fast.
//! Measured by work/twosite-perf/task2_threshold_bench.cpp on this task:
//! with 8 OpenMP threads, global_index_fast serial evaluation is faster up
//! to local_size 4096, while the l2g/OpenMP sweep wins at 8192 and above.
constexpr std::size_t kSerialSweepThreshold = 8192;

//! Whether SignEval::automatic picks the table: rank must allow it and the
//! 2^rank table must not exceed the local slice being swept.
inline bool use_sign_table(std::size_t rank, std::size_t local_size) {
  if (rank > kMaxTableRank) {
    return false;
  }
  const std::size_t table_size = std::size_t(1) << rank;
  return table_size <= local_size;
}

/*!
 * @brief Evaluate a swap form on one element index.
 *
 * @param[in] terms Pair set of the form (SwapForm::terms()).
 * @param[in] idx Element index, one entry per leg.
 * @param[in] parity_bits All leg ledgers flattened into one array.
 * @param[in] parity_offsets Start of each leg's ledger in @p parity_bits.
 * @return +1 or -1: @f$(-1)^{\#\{(x,y)\in terms : p_x = p_y = 1\}}@f$.
 */
inline int direct_sign(const std::vector<std::pair<int, int>>& terms,
                       const std::size_t* idx,
                       const std::vector<std::uint8_t>& parity_bits,
                       const std::vector<std::size_t>& parity_offsets) {
  int exponent = 0;
  for (std::size_t i = 0; i < terms.size(); ++i) {
    const int ax1 = terms[i].first;
    const int ax2 = terms[i].second;
    if (parity_bits[parity_offsets[ax1] + idx[ax1]] != 0 &&
        parity_bits[parity_offsets[ax2] + idx[ax2]] != 0) {
      exponent ^= 1;
    }
  }
  return exponent == 0 ? 1 : -1;
}

//! Position of the (single) set bit, i.e. the axis a one-bit mask names.
inline int lowest_bit_axis(std::uint32_t bit) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_ctz(bit);
#else
  int axis = 0;
  while ((std::uint32_t(1) << axis) != bit) {
    ++axis;
  }
  return axis;
#endif
}

//! True iff an odd number of bits is set.
inline bool popcount_odd(std::uint32_t bits) {
#if defined(__GNUC__) || defined(__clang__)
  return (__builtin_popcount(bits) & 1) != 0;
#else
  bool odd = false;
  while (bits != 0) {
    odd = !odd;
    bits &= bits - 1;
  }
  return odd;
#endif
}

//! True iff the permutation is the identity.
inline bool is_identity_axes(const mptensor::Axes& axes) {
  for (std::size_t i = 0; i < axes.size(); ++i) {
    if (axes[i] != i) {
      return false;
    }
  }
  return true;
}

//! Throw std::runtime_error(context) unless axes is a permutation of
//! 0..rank-1.
inline void validate_axes(const mptensor::Axes& axes, std::size_t rank,
                          const char* context) {
  if (axes.size() != rank) {
    throw std::runtime_error(context);
  }
  std::vector<std::uint8_t> seen(rank, 0);
  for (std::size_t i = 0; i < axes.size(); ++i) {
    if (axes[i] >= rank || seen[axes[i]] != 0) {
      throw std::runtime_error(context);
    }
    seen[axes[i]] = 1;
  }
}

/*!
 * @brief Canonical pair set of a leg permutation: its inversion pairs.
 *
 * Applying the returned form realizes the Koszul sign of transposing by
 * @p axes (the same sign detail::transpose_sign() evaluates per element).
 * Pairs are expressed in pre-transpose axis labels.
 */
inline SwapForm transpose_sign_form(const mptensor::Axes& axes) {
  SwapForm form;
  for (std::size_t x = 0; x < axes.size(); ++x) {
    for (std::size_t y = x + 1; y < axes.size(); ++y) {
      if (axes[x] > axes[y]) {
        form.toggle(static_cast<int>(axes[y]), static_cast<int>(axes[x]));
      }
    }
  }
  return form;
}
}  // namespace detail

/*!
 * @brief One elementwise pass applying a swap form and leg gauges in place.
 *
 * Every element is multiplied by the form's sign (see SwapForm) and by the
 * gauge factors of its index values. Bundling all pending sign work into a
 * single sweep is the point of the canonical form: the pass over the
 * (distributed) tensor dominates the cost, so callers merge their forms
 * with SwapForm::toggle() and sweep once.
 *
 * The evaluation strategy follows @p eval; note that an explicit
 * SignEval::table request falls back to direct evaluation when the rank
 * exceeds detail::kMaxTableRank. Small local slices (below
 * detail::kSerialSweepThreshold) are swept serially, larger ones with
 * OpenMP.
 *
 * @param[in,out] a Tensor to sweep; parities are unchanged.
 * @param[in] form Pair set to apply.
 * @param[in] gauges Per-leg diagonal factors to apply in the same pass.
 * @param[in] eval Evaluation strategy (default: automatic).
 * @throw std::runtime_error If the ledger does not match the tensor shape,
 *        or a form/gauge axis is out of range.
 */
template <class tensor>
void apply_sign_sweep(ftensor<tensor>& a, const SwapForm& form,
                      const std::vector<LegGauge>& gauges,
                      SignEval eval = SignEval::automatic) {
  const std::size_t rank = a.t.shape().size();
  if (a.parity.size() != rank) {
    throw std::runtime_error("fermion sign sweep: parity rank mismatch");
  }
  for (std::size_t ax = 0; ax < rank; ++ax) {
    if (a.parity[ax].size() != a.t.shape()[ax]) {
      throw std::runtime_error("fermion sign sweep: parity dimension mismatch");
    }
  }
  for (std::size_t i = 0; i < gauges.size(); ++i) {
    if (gauges[i].axis < 0 ||
        static_cast<std::size_t>(gauges[i].axis) >= rank) {
      throw std::runtime_error("fermion sign sweep: gauge axis out of range");
    }
  }
  const std::size_t n_local = a.t.local_size();
  const std::vector<std::pair<int, int>>& terms = form.terms();
  const bool has_form = !terms.empty();
  for (std::size_t i = 0; i < terms.size(); ++i) {
    if (terms[i].first < 0 || terms[i].second < 0 ||
        static_cast<std::size_t>(terms[i].first) >= rank ||
        static_cast<std::size_t>(terms[i].second) >= rank) {
      throw std::runtime_error("fermion sign sweep: form axis out of range");
    }
  }
  if (!has_form && gauges.empty()) {
    return;
  }
  const bool use_table =
      (eval == SignEval::table && rank <= detail::kMaxTableRank) ||
      (eval == SignEval::automatic && detail::use_sign_table(rank, n_local));
  const bool build_table = has_form && use_table;

  std::vector<std::size_t> parity_offsets;
  std::vector<std::uint8_t> parity_bits;
  std::vector<std::uint32_t> shifted_parity_bits;
  std::vector<int> form_axes;
  if (has_form) {
    parity_offsets.assign(rank + 1, 0);
    for (std::size_t ax = 0; ax < rank; ++ax) {
      parity_offsets[ax + 1] = parity_offsets[ax] + a.parity[ax].size();
    }
    parity_bits.assign(parity_offsets[rank], 0);
    if (build_table) {
      shifted_parity_bits.assign(parity_offsets[rank], 0);
    }
    for (std::size_t ax = 0; ax < rank; ++ax) {
      for (std::size_t i = 0; i < a.parity[ax].size(); ++i) {
        if (a.parity[ax][i]) {
          parity_bits[parity_offsets[ax] + i] = 1;
          if (build_table) {
            shifted_parity_bits[parity_offsets[ax] + i] = std::uint32_t(1)
                                                          << ax;
          }
        }
      }
    }

    std::vector<std::uint8_t> seen(rank, 0);
    for (std::size_t i = 0; i < terms.size(); ++i) {
      if (seen[terms[i].first] == 0) {
        form_axes.push_back(terms[i].first);
        seen[terms[i].first] = 1;
      }
      if (seen[terms[i].second] == 0) {
        form_axes.push_back(terms[i].second);
        seen[terms[i].second] = 1;
      }
    }
  }

  std::vector<std::int8_t> signs;
  if (build_table) {
    std::vector<std::uint32_t> adjacent(rank, 0);
    for (std::size_t i = 0; i < terms.size(); ++i) {
      adjacent[terms[i].first] |= std::uint32_t(1) << terms[i].second;
      adjacent[terms[i].second] |= std::uint32_t(1) << terms[i].first;
    }
    signs.assign(std::size_t(1) << rank, 1);
    for (std::uint32_t mask = 1; mask < signs.size(); ++mask) {
      const std::uint32_t bit = mask & (~mask + 1);
      const std::uint32_t rest = mask ^ bit;
      const int ax = detail::lowest_bit_axis(bit);
      int sign = signs[rest];
      if (detail::popcount_odd(adjacent[ax] & rest)) {
        sign = -sign;
      }
      signs[mask] = static_cast<std::int8_t>(sign);
    }
  }

  if (n_local < detail::kSerialSweepThreshold) {
    mptensor::Index idx;
    idx.resize(rank);
    if (!has_form) {
      for (std::size_t n = 0; n < n_local; ++n) {
        a.t.global_index_fast(n, idx);
        for (std::size_t i = 0; i < gauges.size(); ++i) {
          const LegGauge& g = gauges[i];
          a.t[n] *= g.factor[idx[g.axis]];
        }
      }
    } else if (build_table) {
      for (std::size_t n = 0; n < n_local; ++n) {
        a.t.global_index_fast(n, idx);
        std::uint32_t mask = 0;
        for (std::size_t i = 0; i < form_axes.size(); ++i) {
          const int ax = form_axes[i];
          mask |= shifted_parity_bits[parity_offsets[ax] + idx[ax]];
        }
        if (signs[mask] < 0) {
          a.t[n] = -a.t[n];
        }
        for (std::size_t i = 0; i < gauges.size(); ++i) {
          const LegGauge& g = gauges[i];
          a.t[n] *= g.factor[idx[g.axis]];
        }
      }
    } else {
      for (std::size_t n = 0; n < n_local; ++n) {
        a.t.global_index_fast(n, idx);
        if (detail::direct_sign(terms, &idx[0], parity_bits, parity_offsets) <
            0) {
          a.t[n] = -a.t[n];
        }
        for (std::size_t i = 0; i < gauges.size(); ++i) {
          const LegGauge& g = gauges[i];
          a.t[n] *= g.factor[idx[g.axis]];
        }
      }
    }
    return;
  }

  a.t.prep_local_to_global();
  a.t.make_l2g_map();
  if (!has_form) {
#ifndef _NO_OMP
#pragma omp parallel default(shared)
#endif
    {
      std::vector<std::size_t> idx(rank);
#ifndef _NO_OMP
#pragma omp for
#endif
      for (std::size_t n = 0; n < n_local; ++n) {
        a.t.global_index_l2g_map(n, idx.data());
        for (std::size_t i = 0; i < gauges.size(); ++i) {
          const LegGauge& g = gauges[i];
          a.t[n] *= g.factor[idx[g.axis]];
        }
      }
    }
  } else if (build_table) {
#ifndef _NO_OMP
#pragma omp parallel default(shared)
#endif
    {
      std::vector<std::size_t> idx(rank);
#ifndef _NO_OMP
#pragma omp for
#endif
      for (std::size_t n = 0; n < n_local; ++n) {
        a.t.global_index_l2g_map(n, idx.data());
        std::uint32_t mask = 0;
        for (std::size_t i = 0; i < form_axes.size(); ++i) {
          const int ax = form_axes[i];
          mask |= shifted_parity_bits[parity_offsets[ax] + idx[ax]];
        }
        if (signs[mask] < 0) {
          a.t[n] = -a.t[n];
        }
        for (std::size_t i = 0; i < gauges.size(); ++i) {
          const LegGauge& g = gauges[i];
          a.t[n] *= g.factor[idx[g.axis]];
        }
      }
    }
  } else {
    // Explicit table requests for rank > kMaxTableRank also take this path:
    // constructing the table would require undefined shifts and impractical
    // memory, so direct evaluation is the defined fallback.
#ifndef _NO_OMP
#pragma omp parallel default(shared)
#endif
    {
      std::vector<std::size_t> idx(rank);
#ifndef _NO_OMP
#pragma omp for
#endif
      for (std::size_t n = 0; n < n_local; ++n) {
        a.t.global_index_l2g_map(n, idx.data());
        if (detail::direct_sign(terms, idx.data(), parity_bits,
                                parity_offsets) < 0) {
          a.t[n] = -a.t[n];
        }
        for (std::size_t i = 0; i < gauges.size(); ++i) {
          const LegGauge& g = gauges[i];
          a.t[n] *= g.factor[idx[g.axis]];
        }
      }
    }
  }
}

/*!
 * @brief Apply a swap form in place: apply_sign_sweep() without gauges.
 */
template <class tensor>
void apply_swap_form(ftensor<tensor>& a, const SwapForm& form,
                     SignEval eval = SignEval::automatic) {
  apply_sign_sweep(a, form, std::vector<LegGauge>(), eval);
}

/*!
 * @brief Apply leg gauges to a PLAIN tensor in place.
 *
 * Same sweep as apply_sign_sweep() restricted to gauges, for callers whose
 * tensor is not (or no longer) wrapped in an ftensor — no ledger is needed
 * since gauges act on index values directly.
 *
 * @param[in,out] a Plain tensor.
 * @param[in] gauges Per-leg diagonal factors.
 * @throw std::runtime_error If a gauge axis is out of range.
 */
template <class tensor>
void apply_leg_gauges(tensor& a, const std::vector<LegGauge>& gauges) {
  const std::size_t rank = a.shape().size();
  const std::size_t n_local = a.local_size();
  for (std::size_t i = 0; i < gauges.size(); ++i) {
    if (gauges[i].axis < 0 ||
        static_cast<std::size_t>(gauges[i].axis) >= rank) {
      throw std::runtime_error("fermion sign sweep: gauge axis out of range");
    }
  }
  if (gauges.empty()) {
    return;
  }
  if (n_local < detail::kSerialSweepThreshold) {
    mptensor::Index idx;
    idx.resize(rank);
    for (std::size_t n = 0; n < n_local; ++n) {
      a.global_index_fast(n, idx);
      for (std::size_t i = 0; i < gauges.size(); ++i) {
        const LegGauge& g = gauges[i];
        a[n] *= g.factor[idx[g.axis]];
      }
    }
    return;
  }
  a.prep_local_to_global();
  a.make_l2g_map();
#ifndef _NO_OMP
#pragma omp parallel default(shared)
#endif
  {
    std::vector<std::size_t> idx(rank);
#ifndef _NO_OMP
#pragma omp for
#endif
    for (std::size_t n = 0; n < n_local; ++n) {
      a.global_index_l2g_map(n, idx.data());
      for (std::size_t i = 0; i < gauges.size(); ++i) {
        const LegGauge& g = gauges[i];
        a[n] *= g.factor[idx[g.axis]];
      }
    }
  }
}

/*!
 * @brief Graded transpose with an extra swap form merged into its sweep.
 *
 * Merges @p form with the Koszul form of @p axes
 * (detail::transpose_sign_form()), sweeps once, then permutes the tensor
 * legs and the parity ledgers. Equivalent to apply_swap_form() followed by
 * a graded transpose, at the cost of a single sweep. @p form is expressed
 * in pre-transpose axis labels.
 *
 * @param[in,out] a Tensor to transpose.
 * @param[in] form Extra pair set to apply in the same sweep.
 * @param[in] axes Permutation, as in mptensor transpose.
 * @param[in] eval Evaluation strategy (default: automatic).
 * @throw std::runtime_error If @p axes is not a permutation of the legs.
 */
template <class tensor>
void transpose_with_swap_form(ftensor<tensor>& a, const SwapForm& form,
                              const mptensor::Axes& axes,
                              SignEval eval = SignEval::automatic) {
  const std::size_t rank = a.t.shape().size();
  detail::validate_axes(axes, rank,
                        "fermion sign sweep: transpose axes out of range");
  if (detail::is_identity_axes(axes)) {
    apply_swap_form(a, form, eval);
    return;
  }

  SwapForm combined = form;
  const SwapForm transpose_form = detail::transpose_sign_form(axes);
  const std::vector<std::pair<int, int>>& transpose_terms =
      transpose_form.terms();
  for (std::size_t i = 0; i < transpose_terms.size(); ++i) {
    combined.toggle(transpose_terms[i].first, transpose_terms[i].second);
  }
  apply_swap_form(a, combined, eval);
  a.t.transpose(axes);
  leg_parities next;
  next.reserve(axes.size());
  for (std::size_t i = 0; i < axes.size(); ++i) {
    next.push_back(a.parity[axes[i]]);
  }
  a.parity = next;
}

template <class tensor>
ftensor<tensor>& ftensor<tensor>::transpose(const mptensor::Axes& axes) {
  transpose_with_swap_form(*this, SwapForm{}, axes);
  return *this;
}

}  // namespace fermion
}  // namespace tenes

#endif  // TENES_SRC_FERMION_SIGN_SWEEP_HPP_
