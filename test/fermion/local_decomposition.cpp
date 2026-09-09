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

// ===== Small blocks are decomposed on one rank, not over the grid ==========
//
// The graded decompositions split every matrix by parity and hand LAPACK the
// two diagonal blocks, so where the bosonic path decomposes one D x D matrix
// this one decomposes two of about half that side. A fermionic full update at
// D = 3, d = 4 ends up asking ScaLAPACK to factorize a 2x2 and a 1x1.
//
// That is not merely wasteful. pdgesvd cross-checks the singular values it
// computed on each process of the grid and returns INFO = MIN(M,N)+1 when
// they do not agree; on nearly degenerate blocks - which small blocks often
// are - the ranks stop agreeing and the run ends. A 2x2 matrix lives inside a
// single 16x16 block-cyclic block anyway, so spreading it over a grid buys
// nothing and costs exactly this failure mode.
//
// So blocks up to a fixed element budget are gathered onto every rank and
// factorized there. gather() is an MPI_Allreduce, so every rank starts from
// the same matrix, and identical inputs to the same routine cannot disagree.
//
// These cases pin the budget, its runtime override, and - the part that
// matters - that a block within the budget really went through the
// non-distributed routine rather than through ScaLAPACK.

#include "../test_fermion_common.hpp"

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

constexpr const char* kBudgetEnv = "TENES_FERMION_LOCAL_DECOMP_MAX";

//! Sets TENES_FERMION_LOCAL_DECOMP_MAX for the life of the object and puts
//! the environment back afterwards, so one case cannot leak its budget into
//! the next. A null value means "unset".
class fld_budget {
 public:
  explicit fld_budget(const char* value) {
    const char* previous = std::getenv(kBudgetEnv);
    had_previous_ = previous != nullptr;
    if (had_previous_) {
      previous_ = previous;
    }
    if (value == nullptr) {
      unsetenv(kBudgetEnv);
    } else {
      setenv(kBudgetEnv, value, 1);
    }
  }
  ~fld_budget() {
    if (had_previous_) {
      setenv(kBudgetEnv, previous_.c_str(), 1);
    } else {
      unsetenv(kBudgetEnv);
    }
  }
  fld_budget(const fld_budget&) = delete;
  fld_budget& operator=(const fld_budget&) = delete;

 private:
  bool had_previous_;
  std::string previous_;
};

//! The entry of the test matrices, fixed by the global index alone so that
//! every rank holds the same matrix however a tensor is distributed, and so
//! that the same matrix can be built in either tensor type.
double fld_entry(const mptensor::Index& idx, std::size_t cols) {
  const double x = static_cast<double>(idx[0] * cols + idx[1] + 1);
  return std::sin(x) + 0.25 * std::cos(3.0 * x);
}

//! fld_entry as the solver's tensor type: distributed under MPI.
tenes::real_tensor fld_matrix(std::size_t rows, std::size_t cols) {
  tenes::real_tensor m(mptensor::Shape(rows, cols));
  for (std::size_t n = 0; n < m.local_size(); ++n) {
    const mptensor::Index idx = m.global_index(n);
    m.set_value(idx, fld_entry(idx, cols));
  }
  return m;
}

//! fld_entry as a non-distributed tensor. Built from the closed form rather
//! than from fld_matrix, so a reference decomposition taken from it owes
//! nothing to the gathering the code under test does.
tenes::small_tensor<double> fld_local_matrix(std::size_t rows,
                                             std::size_t cols) {
  tenes::small_tensor<double> m(mptensor::Shape(rows, cols));
  for (std::size_t n = 0; n < m.local_size(); ++n) {
    const mptensor::Index idx = m.global_index(n);
    m.set_value(idx, fld_entry(idx, cols));
  }
  return m;
}

//! max |a - b| over two rank-2 tensors of the same shape.
double fld_max_diff(const tenes::real_tensor& a, const tenes::real_tensor& b) {
  double worst = 0.0;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    const mptensor::Index idx = a.global_index(n);
    double va = 0.0, vb = 0.0;
    a.get_value(idx, va);
    b.get_value(idx, vb);
    worst = std::max(worst, std::abs(va - vb));
  }
  return worst;
}

}  // namespace

// ---- the budget ------------------------------------------------------------

TEST_CASE("the local-decomposition budget defaults to 4096 elements") {
  fld_budget budget(nullptr);
  CHECK(tenes::fermion::detail::local_decomposition_max_elements() == 4096u);
}

TEST_CASE("a block is decomposed locally exactly up to the budget") {
  fld_budget budget(nullptr);
  using tenes::fermion::detail::fits_local_decomposition;
  // 4096 is where the largest blocks a D = 8, d = 4 fermionic full update
  // produces land: the QR of a site tensor is D^3 x D*d, so its blocks are
  // about 256 x 16, and the truncated SVD of Theta is D*d^2 square, so its
  // blocks are about 64 x 64.
  CHECK(fits_local_decomposition(256, 16));
  CHECK(fits_local_decomposition(64, 64));
  CHECK_FALSE(fits_local_decomposition(64, 65));
  // D = 12 would put 864 x 24 into the QR; that one stays on the grid.
  CHECK_FALSE(fits_local_decomposition(864, 24));
}

TEST_CASE("TENES_FERMION_LOCAL_DECOMP_MAX moves the budget") {
  fld_budget budget("4");
  CHECK(tenes::fermion::detail::local_decomposition_max_elements() == 4u);
  CHECK(tenes::fermion::detail::fits_local_decomposition(2, 2));
  CHECK_FALSE(tenes::fermion::detail::fits_local_decomposition(2, 3));
}

TEST_CASE("TENES_FERMION_LOCAL_DECOMP_MAX=0 keeps every block on the grid") {
  fld_budget budget("0");
  CHECK(tenes::fermion::detail::local_decomposition_max_elements() == 0u);
  CHECK_FALSE(tenes::fermion::detail::fits_local_decomposition(1, 1));
}

// ---- which routine a block actually reaches --------------------------------

TEST_CASE(
    "a block within the budget is factorized by the non-distributed SVD") {
  fld_budget budget(nullptr);
  tenes::real_tensor m = fld_matrix(5, 3);
  tenes::real_tensor u, vt;
  std::vector<double> s;
  const int info = tenes::fermion::detail::block_svd(m, u, s, vt);
  REQUIRE(info == 0);

  // The reference is the routine the local path must have used: gather the
  // block onto every rank, hand it to LAPACK. Equality here is exact, not
  // approximate, and that is the whole point - an approximate check would
  // also pass if the block had gone to pdgesvd, which is what this path
  // exists to avoid.
  tenes::small_tensor<double> local = fld_local_matrix(5, 3);
  tenes::small_tensor<double> local_u, local_vt;
  std::vector<double> s_reference;
  REQUIRE(mptensor::svd(local, mptensor::Axes(0), mptensor::Axes(1), local_u,
                        s_reference, local_vt) == 0);
  REQUIRE(s.size() == s_reference.size());
  for (std::size_t i = 0; i < s.size(); ++i) {
    CHECK(s[i] == s_reference[i]);
  }
}

TEST_CASE("a block within the budget is factorized by the non-distributed QR") {
  fld_budget budget(nullptr);
  tenes::real_tensor m = fld_matrix(5, 3);
  tenes::real_tensor q, r;
  const int info = tenes::fermion::detail::block_qr(m, q, r);
  REQUIRE(info == 0);

  tenes::small_tensor<double> local = fld_local_matrix(5, 3);
  tenes::small_tensor<double> local_q, local_r;
  REQUIRE(mptensor::qr(local, mptensor::Axes(0), mptensor::Axes(1), local_q,
                       local_r) == 0);
  REQUIRE(r.shape() == local_r.shape());
  for (std::size_t n = 0; n < local_r.local_size(); ++n) {
    const mptensor::Index idx = local_r.global_index(n);
    double got = 0.0, want = 0.0;
    r.get_value(idx, got);
    local_r.get_value(idx, want);
    CHECK(got == want);
  }
}

// ---- both sides of the budget still decompose correctly --------------------

TEST_CASE("block_svd reconstructs its block whichever path it takes") {
  const char* setting = nullptr;
  SUBCASE("within the budget") { setting = nullptr; }
  SUBCASE("budget disabled") { setting = "0"; }
  fld_budget budget(setting);

  tenes::real_tensor m = fld_matrix(4, 6);
  tenes::real_tensor u, vt;
  std::vector<double> s;
  REQUIRE(tenes::fermion::detail::block_svd(m, u, s, vt) == 0);
  REQUIRE(s.size() == 4u);
  tenes::real_tensor us = u;
  us.multiply_vector(s, 1);
  const tenes::real_tensor recon =
      mptensor::tensordot(us, vt, mptensor::Axes(1), mptensor::Axes(0));
  CHECK(fld_max_diff(recon, fld_matrix(4, 6)) < 1.0e-12);
}

TEST_CASE("block_qr reconstructs its block whichever path it takes") {
  const char* setting = nullptr;
  SUBCASE("within the budget") { setting = nullptr; }
  SUBCASE("budget disabled") { setting = "0"; }
  fld_budget budget(setting);

  tenes::real_tensor m = fld_matrix(6, 4);
  tenes::real_tensor q, r;
  REQUIRE(tenes::fermion::detail::block_qr(m, q, r) == 0);
  const tenes::real_tensor recon =
      mptensor::tensordot(q, r, mptensor::Axes(1), mptensor::Axes(0));
  CHECK(fld_max_diff(recon, fld_matrix(6, 4)) < 1.0e-12);
}

// ---- the graded layer on top -----------------------------------------------

TEST_CASE("the graded SVD gives the same answer on both sides of the budget") {
  tenes::fermion::leg_parities p{{false, true, false}, {false, true, true}};
  ft a = make_even_ft(mptensor::Shape(3, 3), p, 71);

  ft u_local, vt_local, u_grid, vt_grid;
  std::vector<double> s_local, s_grid;
  {
    fld_budget budget(nullptr);
    REQUIRE(tenes::fermion::svd(a, mptensor::Axes(0), mptensor::Axes(1),
                                u_local, s_local, vt_local) == 0);
  }
  {
    fld_budget budget("0");
    REQUIRE(tenes::fermion::svd(a, mptensor::Axes(0), mptensor::Axes(1), u_grid,
                                s_grid, vt_grid) == 0);
  }
  REQUIRE(s_local.size() == s_grid.size());
  for (std::size_t i = 0; i < s_local.size(); ++i) {
    CHECK(s_local[i] == doctest::Approx(s_grid[i]).epsilon(1.0e-12));
  }
  CHECK(tenes::fermion::parity_violation(u_local) == doctest::Approx(0.0));
  CHECK(tenes::fermion::parity_violation(vt_local) == doctest::Approx(0.0));
  CHECK(u_local.parity[1] == u_grid.parity[1]);

  ft us = u_local;
  us.multiply_vector(s_local, 1);
  const ft recon = tenes::fermion::tensordot(us, vt_local, mptensor::Axes(1),
                                             mptensor::Axes(0));
  CHECK(fld_max_diff(recon.t, a.t) < 1.0e-12);
}

TEST_CASE("the graded QR gives the same answer on both sides of the budget") {
  tenes::fermion::leg_parities p{{false, true, false}, {false, true, true}};
  ft a = make_even_ft(mptensor::Shape(3, 3), p, 61);

  ft q_local, r_local, q_grid, r_grid;
  {
    fld_budget budget(nullptr);
    REQUIRE(tenes::fermion::qr(a, mptensor::Axes(0), mptensor::Axes(1), q_local,
                               r_local) == 0);
  }
  {
    fld_budget budget("0");
    REQUIRE(tenes::fermion::qr(a, mptensor::Axes(0), mptensor::Axes(1), q_grid,
                               r_grid) == 0);
  }
  CHECK(tenes::fermion::parity_violation(q_local) == doctest::Approx(0.0));
  CHECK(tenes::fermion::parity_violation(r_local) == doctest::Approx(0.0));
  CHECK(q_local.parity[1] == q_grid.parity[1]);

  const ft recon = tenes::fermion::tensordot(
      q_local, r_local, mptensor::Axes(1), mptensor::Axes(0));
  CHECK(fld_max_diff(recon.t, a.t) < 1.0e-12);
}
