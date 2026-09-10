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

// ===== Scalars reduced from a distributed tensor must be global ===========
//
// A fermionic Hubbard state measured through the same binary gave the same
// one-site observables at 1, 2 and 4 MPI ranks and a hopping that moved by
// 0.023 between 1 and 2.  Every two-site observable moved; no one-site one
// did.  The two-site CTM contraction ends in detail::trace_boundary_pairs(),
// which sums the elements it holds locally and returns -- a partial sum on
// every rank but one.
//
// ctest runs the MPI build with one rank, where a local sum is the global
// sum, so this file is registered a second time under two ranks.  These
// cases are written so that they pass trivially at one rank and only bite
// with more: the reference is assembled from get_value(), which is
// collective, so both sides of each CHECK are the same on every rank.

#include "../test_fermion_common.hpp"

#include <cmath>

namespace {

//! The entry at a global index: a closed form, so the reference below can
//! be assembled by arithmetic alone.  mptensor's get_value() leaves the
//! value untouched on ranks that do not own the element, so a reference
//! read through it would itself be a per-rank partial sum and the check
//! would pass whatever the code under test did.
double mr_entry(const mptensor::Index& idx, std::size_t n) {
  const double x = static_cast<double>(idx[0] + n * idx[1] + n * n * idx[2] +
                                       n * n * n * idx[3]);
  return std::cos(0.37 * x) + 0.5;
}

//! A rank-4 tensor filled from mr_entry, the same on every rank.
tenes::real_tensor mr_rank4(std::size_t n) {
  tenes::real_tensor a(mptensor::Shape(n, n, n, n));
  for (std::size_t k = 0; k < a.local_size(); ++k) {
    const mptensor::Index idx = a.global_index(k);
    a.set_value(idx, mr_entry(idx, n));
  }
  return a;
}

//! The double-delta trace sum_{i,k} a[i,i,k,k], from the closed form.
double mr_reference_trace(std::size_t n) {
  double value = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t k = 0; k < n; ++k) {
      value += mr_entry(mptensor::Index(i, i, k, k), n);
    }
  }
  return value;
}

}  // namespace

TEST_CASE(
    "trace_boundary_pairs is the global double-delta trace on every rank") {
  // n = 20 puts the matricized 400 x 400 tensor across more than one
  // block-cyclic block in each direction, so at two ranks no rank holds
  // every diagonal-pair element.
  const tenes::real_tensor a = mr_rank4(20);
  const double want = mr_reference_trace(20);
  const double got = tenes::fermion::detail::trace_boundary_pairs(a);
  INFO("ranks = " << a.get_comm_size() << ", rank = " << a.get_comm_rank());
  CHECK(got == doctest::Approx(want).epsilon(1.0e-12));
}
