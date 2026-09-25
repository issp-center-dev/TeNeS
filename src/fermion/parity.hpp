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
 *  @brief Z2 parity bookkeeping for graded (fermionic) tensor legs.
 *
 *  A graded tensor leg is described by a ::tenes::fermion::parity_vector
 *  that records, for every index value of the leg, whether that basis state
 *  carries odd fermion parity. This header provides the primitive operations
 *  on such parity ledgers: fusing two legs into one (fuse()), counting the
 *  odd legs selected by a tensor index (count_odd()), and building the
 *  even-first sorting permutation used by the block decompositions
 *  (parity_sort_perm()).
 */

#ifndef TENES_SRC_FERMION_PARITY_HPP_
#define TENES_SRC_FERMION_PARITY_HPP_

#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <mptensor/tensor.hpp>

#include "../mpi.hpp"

namespace tenes {
namespace fermion {

//! Per-leg parity ledger: parity_vector[i] is true iff index value i is odd.
using parity_vector = std::vector<bool>;
//! Parity ledgers for all legs of a tensor, in leg order.
using leg_parities = std::vector<parity_vector>;

/*!
 * @brief Fuse two leg parities into the parity of the combined index.
 *
 * Uses the COLUMN-MAJOR flattening convention of mptensor (ScaLAPACK
 * heritage): mptensor::reshape combines adjacent legs as
 * index = i + a.size() * j, i.e. the FIRST leg is the fastest-varying one.
 * A row-major fuse here silently mislabels every fused index whose leg
 * parity patterns are asymmetric.
 *
 * @param[in] a Parity ledger of the first (fastest-varying) leg.
 * @param[in] b Parity ledger of the second leg.
 * @return Parity ledger of the fused leg of dimension a.size() * b.size();
 *         a fused index is odd iff exactly one of its two factors is odd.
 */
inline parity_vector fuse(const parity_vector& a, const parity_vector& b) {
  parity_vector r(a.size() * b.size());
  for (std::size_t j = 0; j < b.size(); ++j) {
    for (std::size_t i = 0; i < a.size(); ++i) {
      r[i + a.size() * j] = a[i] != b[j];
    }
  }
  return r;
}

/*!
 * @brief Count how many legs selected by a tensor index carry odd parity.
 *
 * @param[in] p Parity ledgers of all legs.
 * @param[in] idx Tensor index (one index value per leg).
 * @return Number of legs whose selected index value is odd.
 */
inline int count_odd(const leg_parities& p, const mptensor::Index& idx) {
  int m = 0;
  for (std::size_t ax = 0; ax < p.size(); ++ax) {
    m += p[ax][idx[ax]] ? 1 : 0;
  }
  return m;
}

/*!
 * @brief Build the even-first sorting permutation of a leg.
 *
 * The block decompositions (fermion qr() / svd()) sort a fused leg so that
 * all even index values precede all odd ones; this returns the permutation
 * that realizes the sort. Within each parity sector the original order is
 * kept (the sort is stable).
 *
 * @param[in] p Parity ledger of the leg.
 * @return perm such that perm[k] is the original position of the k-th
 *         index value after the even-first sort.
 */
inline std::vector<std::size_t> parity_sort_perm(const parity_vector& p) {
  std::vector<std::size_t> perm;
  perm.reserve(p.size());
  for (std::size_t i = 0; i < p.size(); ++i) {
    if (!p[i]) {
      perm.push_back(i);
    }
  }
  for (std::size_t i = 0; i < p.size(); ++i) {
    if (p[i]) {
      perm.push_back(i);
    }
  }
  return perm;
}

//! Parity class of a one-site operator.
enum class op_parity { even, odd, mixed };

/*!
 * @brief Parity class of a one-site operator op[in, out].
 *
 * even: every nonzero element keeps the parity of the physical index; odd:
 * every nonzero element flips it; mixed: both kinds occur. The zero
 * operator is even. Collective over the tensor's communicator: every rank
 * returns the same value.
 *
 * @param[in] op Plain rank-2 one-site operator.
 * @param[in] phys Physical-leg parity ledger.
 */
template <class tensor>
op_parity operator_parity(const tensor& op, const parity_vector& phys) {
  if (op.rank() != 2) {
    throw std::invalid_argument("operator_parity expects a rank-2 operator");
  }
  if (op.shape()[0] != phys.size() || op.shape()[1] != phys.size()) {
    throw std::invalid_argument(
        "operator_parity: physical parity size does not match operator");
  }

  double local_even = 0.0;
  double local_odd = 0.0;
  mptensor::Index idx;
  idx.resize(2);
  for (std::size_t n = 0; n < op.local_size(); ++n) {
    op.global_index_fast(n, idx);
    const double a = std::abs(op[n]);
    if (phys[idx[0]] == phys[idx[1]]) {
      local_even = std::max(local_even, a);
    } else {
      local_odd = std::max(local_odd, a);
    }
  }

  std::vector<double> reduced{local_even, local_odd};
  tenes::allreduce_max(reduced, op.get_comm());
  const bool has_even = reduced[0] > 0.0;
  const bool has_odd = reduced[1] > 0.0;
  if (has_even && has_odd) {
    return op_parity::mixed;
  }
  return has_odd ? op_parity::odd : op_parity::even;
}

}  // namespace fermion
}  // namespace tenes

#endif  // TENES_SRC_FERMION_PARITY_HPP_
