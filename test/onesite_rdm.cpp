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

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <random>
#include <vector>

#include "../src/tensor.hpp"
#include "../src/iTPS/core/contract_itps_ctm.hpp"
#include "../src/iTPS/core/contract_density_ctm.hpp"
#include "../src/iTPS/core/contract_itps_mf.hpp"

/*
 * Consistency of the one-site RDM contraction functions against the
 * existing scalar contraction functions.
 *
 * The three functions under test return the (unnormalized) d x d
 * one-site reduced density matrix and must satisfy, for any op:
 *
 *  - iTPS CTM:
 *      trace(RDM, op, Axes(0,1), Axes(1,0))
 *        == Contract_one_site_iTPS_CTM(C1..C4, eT1..eT4, Tn1, op)
 *  - density CTM:
 *      trace(op, RDM, Axes(0,1), Axes(0,1))
 *        == Contract_one_site_density_CTM(C1..C4, eT1..eT4, Tn1, op)
 *  - MF (Tn1 has the mean-field environment already absorbed,
 *    the same precondition as the existing scalar function):
 *      trace(op, RDM, Axes(0,1), Axes(0,1))
 *        == Contract_one_site_iTPS_MF(Tn1, op)
 */

namespace {

using mptensor::Axes;
using mptensor::Index;
using mptensor::Shape;

template <class value_type>
value_type random_value(std::mt19937 &gen,
                        std::uniform_real_distribution<double> &dist) {
  if constexpr (std::is_same_v<value_type, double>) {
    return dist(gen);
  } else {
    return value_type(dist(gen), dist(gen));
  }
}

//! Fill a tensor with deterministic pseudo-random values.
template <class tensor>
void fill_random(tensor &A, std::mt19937 &gen) {
  using value_type = typename tensor::value_type;
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  const Shape shape = A.shape();
  const std::size_t rank = shape.size();
  std::vector<std::size_t> idx(rank, 0);
  while (true) {
    Index index;
    for (std::size_t i = 0; i < rank; ++i) {
      index.push(idx[i]);
    }
    A.set_value(index, random_value<value_type>(gen, dist));
    std::size_t k = 0;
    while (k < rank) {
      ++idx[k];
      if (idx[k] < shape[k]) {
        break;
      }
      idx[k] = 0;
      ++k;
    }
    if (k == rank) {
      break;
    }
  }
}

//! Build the three test operators: identity, diagonal, non-symmetric.
template <class tensor>
std::vector<tensor> make_operators(int d) {
  using value_type = typename tensor::value_type;
  std::vector<tensor> ops;

  tensor identity(Shape(d, d));
  tensor diagonal(Shape(d, d));
  tensor nonsym(Shape(d, d));
  std::mt19937 gen(54321);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (int i = 0; i < d; ++i) {
    for (int j = 0; j < d; ++j) {
      identity.set_value(Index(i, j),
                         (i == j) ? value_type(1.0) : value_type(0.0));
      diagonal.set_value(
          Index(i, j), (i == j) ? value_type(0.3 - 0.9 * i) : value_type(0.0));
      // dense and non-symmetric (complex entries in the complex case)
      nonsym.set_value(Index(i, j), random_value<value_type>(gen, dist));
    }
  }
  ops.push_back(identity);
  ops.push_back(diagonal);
  ops.push_back(nonsym);
  return ops;
}

template <class value_type>
double relative_difference(value_type a, value_type b) {
  const double scale =
      std::max({std::abs(a), std::abs(b), std::numeric_limits<double>::min()});
  return std::abs(a - b) / scale;
}

//! Sum of the diagonal elements of a d x d tensor.
template <class tensor>
typename tensor::value_type diagonal_sum(const tensor &A, int d) {
  typename tensor::value_type sum = 0.0;
  for (int i = 0; i < d; ++i) {
    typename tensor::value_type v;
    const bool has = A.get_value(Index(i, i), v);
    REQUIRE(has);
    sum += v;
  }
  return sum;
}

const int chi = 3;
const int D = 2;
const int d = 2;
const double tol = 1.0e-12;

template <class tensor>
void check_rdm_itps_ctm() {
  std::mt19937 gen(12345);

  tensor C1(Shape(chi, chi)), C2(Shape(chi, chi)), C3(Shape(chi, chi)),
      C4(Shape(chi, chi));
  tensor eT1(Shape(chi, chi, D, D)), eT2(Shape(chi, chi, D, D)),
      eT3(Shape(chi, chi, D, D)), eT4(Shape(chi, chi, D, D));
  tensor Tn1(Shape(D, D, D, D, d));
  fill_random(C1, gen);
  fill_random(C2, gen);
  fill_random(C3, gen);
  fill_random(C4, gen);
  fill_random(eT1, gen);
  fill_random(eT2, gen);
  fill_random(eT3, gen);
  fill_random(eT4, gen);
  fill_random(Tn1, gen);

  const tensor rdm = tenes::itps::core::Contract_one_site_RDM_iTPS_CTM(
      C1, C2, C3, C4, eT1, eT2, eT3, eT4, Tn1);
  REQUIRE(rdm.rank() == 2);
  REQUIRE(rdm.shape() == Shape(d, d));

  const auto ops = make_operators<tensor>(d);
  for (std::size_t iop = 0; iop < ops.size(); ++iop) {
    CAPTURE(iop);
    const auto expected = tenes::itps::core::Contract_one_site_iTPS_CTM(
        C1, C2, C3, C4, eT1, eT2, eT3, eT4, Tn1, ops[iop]);
    const auto actual = trace(rdm, ops[iop], Axes(0, 1), Axes(1, 0));
    CHECK(relative_difference(expected, actual) <= tol);
  }

  // trace(RDM) equals the scalar contraction with the identity operator
  const auto norm = tenes::itps::core::Contract_one_site_iTPS_CTM(
      C1, C2, C3, C4, eT1, eT2, eT3, eT4, Tn1, ops[0]);
  CHECK(relative_difference(diagonal_sum(rdm, d), norm) <= tol);
}

template <class tensor>
void check_rdm_density_ctm() {
  std::mt19937 gen(23456);

  tensor C1(Shape(chi, chi)), C2(Shape(chi, chi)), C3(Shape(chi, chi)),
      C4(Shape(chi, chi));
  tensor eT1(Shape(chi, chi, D)), eT2(Shape(chi, chi, D)),
      eT3(Shape(chi, chi, D)), eT4(Shape(chi, chi, D));
  tensor Tn1(Shape(D, D, D, D, d, d));
  fill_random(C1, gen);
  fill_random(C2, gen);
  fill_random(C3, gen);
  fill_random(C4, gen);
  fill_random(eT1, gen);
  fill_random(eT2, gen);
  fill_random(eT3, gen);
  fill_random(eT4, gen);
  fill_random(Tn1, gen);

  const tensor rdm = tenes::itps::core::Contract_one_site_RDM_density_CTM(
      C1, C2, C3, C4, eT1, eT2, eT3, eT4, Tn1);
  REQUIRE(rdm.rank() == 2);
  REQUIRE(rdm.shape() == Shape(d, d));

  const auto ops = make_operators<tensor>(d);
  for (std::size_t iop = 0; iop < ops.size(); ++iop) {
    CAPTURE(iop);
    const auto expected = tenes::itps::core::Contract_one_site_density_CTM(
        C1, C2, C3, C4, eT1, eT2, eT3, eT4, Tn1, ops[iop]);
    const auto actual = trace(ops[iop], rdm, Axes(0, 1), Axes(0, 1));
    CHECK(relative_difference(expected, actual) <= tol);
  }

  const auto norm = tenes::itps::core::Contract_one_site_density_CTM(
      C1, C2, C3, C4, eT1, eT2, eT3, eT4, Tn1, ops[0]);
  CHECK(relative_difference(diagonal_sum(rdm, d), norm) <= tol);
}

template <class tensor>
void check_rdm_itps_mf() {
  std::mt19937 gen(34567);

  // The MF one-site contraction requires Tn1 with the mean-field
  // environment absorbed; the algebraic identity checked here holds
  // for any tensor with that layout.
  tensor Tn1(Shape(D, D, D, D, d));
  fill_random(Tn1, gen);

  const tensor rdm = tenes::itps::core::Contract_one_site_RDM_iTPS_MF(Tn1);
  REQUIRE(rdm.rank() == 2);
  REQUIRE(rdm.shape() == Shape(d, d));

  const auto ops = make_operators<tensor>(d);
  for (std::size_t iop = 0; iop < ops.size(); ++iop) {
    CAPTURE(iop);
    const auto expected =
        tenes::itps::core::Contract_one_site_iTPS_MF(Tn1, ops[iop]);
    const auto actual = trace(ops[iop], rdm, Axes(0, 1), Axes(0, 1));
    CHECK(relative_difference(expected, actual) <= tol);
  }

  const auto norm = tenes::itps::core::Contract_one_site_iTPS_MF(Tn1, ops[0]);
  CHECK(relative_difference(diagonal_sum(rdm, d), norm) <= tol);
}

}  // namespace

TEST_CASE("one-site RDM contraction (real)") {
  using tensor = tenes::real_tensor;
  SUBCASE("iTPS CTM") { check_rdm_itps_ctm<tensor>(); }
  SUBCASE("density CTM") { check_rdm_density_ctm<tensor>(); }
  SUBCASE("iTPS MF") { check_rdm_itps_mf<tensor>(); }
}

TEST_CASE("one-site RDM contraction (complex)") {
  using tensor = tenes::complex_tensor;
  SUBCASE("iTPS CTM") { check_rdm_itps_ctm<tensor>(); }
  SUBCASE("density CTM") { check_rdm_density_ctm<tensor>(); }
  SUBCASE("iTPS MF") { check_rdm_itps_mf<tensor>(); }
}
