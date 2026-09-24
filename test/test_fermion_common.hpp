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

//! @file
//! Definitions shared by the translation units of the test_fermion_layer
//! executable.  The fermion test cases used to be one 18,000-line translation
//! unit (test_fermion_layer.cpp including every test/fermion/*.cpp at its
//! end), which took about 100 s to compile on a single core while the other
//! three cores of a CI runner had nothing left to do.  They are separate
//! translation units now, so what they share has to live in a header.

#ifndef TENES_TEST_FERMION_COMMON_HPP
#define TENES_TEST_FERMION_COMMON_HPP

#include "doctest.h"
#include "test_workdir.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>

#include "../src/fermion/fermion_info.hpp"
#include "../src/fermion/fops.hpp"
#include "../src/fermion/ftensor.hpp"
#include "../src/fermion/parity.hpp"
#include "../src/fermion/reduced.hpp"
#include "../src/fermion/reduced_measure.hpp"
#include "../src/SquareLattice.hpp"
#include "../src/tensor.hpp"
#include "../src/iTPS/PEPS_Parameters.hpp"
#include "../src/iTPS/iTPS.hpp"
#include "../src/iTPS/core/ctm.hpp"
#include "../src/iTPS/core/contract.hpp"
#include "../src/iTPS/core/contract_itps_ctm.hpp"
#include "../src/iTPS/core/simple_update.hpp"

using namespace tenes::fermion;

using ft = tenes::fermion::ftensor<tenes::real_tensor>;

namespace tenes::itps {
struct iTPSTestAccessor {
  template <class tensor>
  static std::vector<tensor>& Tn(iTPS<tensor>& state) {
    return state.Tn;
  }

  template <class tensor>
  static std::vector<std::vector<std::vector<double>>>& lambda_tensor(
      iTPS<tensor>& state) {
    return state.lambda_tensor;
  }

  template <class tensor>
  static tenes::fermion::FermionInfo& finfo(iTPS<tensor>& state) {
    return state.finfo;
  }

  template <class tensor>
  static std::vector<std::string>& twosite_operator_names(iTPS<tensor>& state) {
    return state.twosite_operator_names;
  }

  template <class tensor>
  static void update_reduced_density_environment(iTPS<tensor>& state) {
    // Bare Tn, matching measure.cpp: the CTM provides the environment.
    const std::vector<tensor> reduced =
        tenes::fermion::build_reduced_density_tensors(state.Tn, state.finfo);
    core::Calc_CTM_Environment_density(
        state.C1, state.C2, state.C3, state.C4, state.eTt, state.eTr, state.eTb,
        state.eTl, reduced, state.peps_parameters, state.lattice);
  }

  // ---- the CTM environment -------------------------------------------------
  //
  // The ten environment slots are private members of iTPS (iTPS.hpp, after the
  // `private:` on line 318). fermion/fast_full_update.cpp needs to read them
  // before and after a single full-update bond to see which CTM move ran, so
  // they are exposed here rather than by widening iTPS' own interface.
  template <class tensor>
  static std::vector<tensor>& C1(iTPS<tensor>& state) {
    return state.C1;
  }
  template <class tensor>
  static std::vector<tensor>& C2(iTPS<tensor>& state) {
    return state.C2;
  }
  template <class tensor>
  static std::vector<tensor>& C3(iTPS<tensor>& state) {
    return state.C3;
  }
  template <class tensor>
  static std::vector<tensor>& C4(iTPS<tensor>& state) {
    return state.C4;
  }
  template <class tensor>
  static std::vector<tensor>& eTt(iTPS<tensor>& state) {
    return state.eTt;
  }
  template <class tensor>
  static std::vector<tensor>& eTr(iTPS<tensor>& state) {
    return state.eTr;
  }
  template <class tensor>
  static std::vector<tensor>& eTb(iTPS<tensor>& state) {
    return state.eTb;
  }
  template <class tensor>
  static std::vector<tensor>& eTl(iTPS<tensor>& state) {
    return state.eTl;
  }

  //! The parameter set and the geometry the state was built with. The
  //! oracle of fermion/fast_full_update.cpp replays the two CTM moves with
  //! exactly the objects the solver used, so that a discrepancy cannot come
  //! from the test having rebuilt them differently.
  template <class tensor>
  static PEPS_Parameters const& peps_parameters(iTPS<tensor>& state) {
    return state.peps_parameters;
  }
  template <class tensor>
  static SquareLattice const& lattice(iTPS<tensor>& state) {
    return state.lattice;
  }
};
}  // namespace tenes::itps

//! A random real ftensor with the given shape and leg parities.
inline ft make_random_ft(const mptensor::Shape& sh,
                         const tenes::fermion::leg_parities& p, unsigned seed) {
  tenes::real_tensor t(sh);
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    t.set_value(t.global_index(n), dist(gen));
  }
  return ft{t, p};
}

inline ft make_even_ft(const mptensor::Shape& sh,
                       const tenes::fermion::leg_parities& p, unsigned seed) {
  ft a = make_random_ft(sh, p, seed);
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    auto idx = a.t.global_index(n);
    if (tenes::fermion::count_odd(p, idx) % 2 == 1) {
      a.t.set_value(idx, 0.0);
    }
  }
  return a;
}

inline tenes::real_tensor make_free_fermion_gate(double tau) {
  tenes::real_tensor gate(mptensor::Shape(2, 2, 2, 2));
  gate.set_value(mptensor::Index(0, 0, 0, 0), 1.0);
  gate.set_value(mptensor::Index(1, 1, 1, 1), 1.0);
  gate.set_value(mptensor::Index(0, 1, 0, 1), std::cosh(tau));
  gate.set_value(mptensor::Index(1, 0, 1, 0), std::cosh(tau));
  gate.set_value(mptensor::Index(0, 1, 1, 0), std::sinh(tau));
  gate.set_value(mptensor::Index(1, 0, 0, 1), std::sinh(tau));
  return gate;
}

inline std::vector<double> sorted_desc(std::vector<double> values) {
  std::sort(values.begin(), values.end(), std::greater<double>());
  return values;
}

inline double lambda_relative_diff(const std::vector<double>& a,
                                   const std::vector<double>& b) {
  double diff = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    const double scale = std::max({1.0e-300, std::abs(a[i]), std::abs(b[i])});
    diff = std::max(diff, std::abs(a[i] - b[i]) / scale);
  }
  return diff;
}

inline std::string vector_to_string(const std::vector<double>& values) {
  std::ostringstream os;
  os << "[";
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i != 0) {
      os << ",";
    }
    os << std::setprecision(17) << values[i];
  }
  os << "]";
  return os.str();
}

// Apply c_mode (dagger=false) or c^dag_mode to basis state g.
// Returns {g', sign}; sign = 0 means the result vanishes.
inline std::pair<int, double> electron_mode_op(int g, int mode, bool dagger) {
  const int bit = 3 - mode;  // mode 0 (1up) is the highest bit
  const int occ = (g >> bit) & 1;
  if (dagger == (occ == 1)) {
    return {0, 0.0};
  }
  int string_count = 0;
  for (int m = 0; m < mode; ++m) {
    string_count += (g >> (3 - m)) & 1;
  }
  const double sign = (string_count % 2 == 0) ? 1.0 : -1.0;
  return {g ^ (1 << bit), sign};
}

// 16x16 matrix of the electron bond Hamiltonian
//   h = -t sum_sigma (c^dag_{1 sigma} c_{2 sigma} + h.c.)
//       + (U/4) (n_{1up} n_{1dn} + n_{2up} n_{2dn})   [U split over 4 bonds]
//       - (mu/4) (n_1 + n_2)
// in the ordered Fock basis, indexed by (i1, i2) as i1 * 4 + i2.
inline std::array<std::array<double, 16>, 16> electron_bond_hamiltonian(
    double t, double u, double mu) {
  auto local_to_bits = [](int i) {  // i = n_up + 2 n_dn -> (n_up, n_dn)
    return std::pair<int, int>{i & 1, (i >> 1) & 1};
  };
  auto pair_to_g = [&](int i1, int i2) {
    const auto [a, b] = local_to_bits(i1);
    const auto [c, d] = local_to_bits(i2);
    return a << 3 | b << 2 | c << 1 | d;
  };
  std::array<int, 16> g_of_index{};
  std::array<int, 16> index_of_g{};
  for (int i1 = 0; i1 < 4; ++i1) {
    for (int i2 = 0; i2 < 4; ++i2) {
      const int idx = i1 * 4 + i2;
      g_of_index[idx] = pair_to_g(i1, i2);
      index_of_g[g_of_index[idx]] = idx;
    }
  }

  std::array<std::array<double, 16>, 16> h{};
  for (int idx = 0; idx < 16; ++idx) {
    const int g = g_of_index[idx];
    // diagonal: U and mu terms
    const int a = (g >> 3) & 1, b = (g >> 2) & 1, c = (g >> 1) & 1, d = g & 1;
    h[idx][idx] += 0.25 * u * (a * b + c * d) - 0.25 * mu * (a + b + c + d);
    // hopping: -t (c^dag_{1s} c_{2s} + c^dag_{2s} c_{1s});
    // modes: 1up=0, 1dn=1, 2up=2, 2dn=3
    const int hop_pairs[4][2] = {{0, 2}, {2, 0}, {1, 3}, {3, 1}};
    for (const auto& mp : hop_pairs) {
      const auto [g1, s1] = electron_mode_op(g, mp[1], false);
      if (s1 == 0.0) {
        continue;
      }
      const auto [g2, s2] = electron_mode_op(g1, mp[0], true);
      if (s2 == 0.0) {
        continue;
      }
      h[index_of_g[g2]][idx] += -t * s1 * s2;
    }
  }
  return h;
}

// Gate tensor exp(-tau h) as op[in1][in2][out1][out2] via a Taylor series
// (norm(tau h) << 1 for the parameters used here).
inline tenes::real_tensor electron_gate(double t, double u, double mu,
                                        double tau) {
  const auto h = electron_bond_hamiltonian(t, u, mu);
  std::array<std::array<double, 16>, 16> gate{};
  std::array<std::array<double, 16>, 16> term{};
  for (int i = 0; i < 16; ++i) {
    gate[i][i] = 1.0;
    term[i][i] = 1.0;
  }
  for (int order = 1; order <= 20; ++order) {
    std::array<std::array<double, 16>, 16> next{};
    for (int i = 0; i < 16; ++i) {
      for (int k = 0; k < 16; ++k) {
        if (term[i][k] == 0.0) {
          continue;
        }
        const double w = term[i][k] * (-tau) / order;
        for (int j = 0; j < 16; ++j) {
          next[i][j] += w * h[k][j];
        }
      }
    }
    term = next;
    for (int i = 0; i < 16; ++i) {
      for (int j = 0; j < 16; ++j) {
        gate[i][j] += term[i][j];
      }
    }
  }
  tenes::real_tensor op(mptensor::Shape(4, 4, 4, 4));
  for (int i1 = 0; i1 < 4; ++i1) {
    for (int i2 = 0; i2 < 4; ++i2) {
      for (int o1 = 0; o1 < 4; ++o1) {
        for (int o2 = 0; o2 < 4; ++o2) {
          const double v = gate[o1 * 4 + o2][i1 * 4 + i2];
          if (std::abs(v) > 1.0e-16) {
            op.set_value(mptensor::Index(i1, i2, o1, o2), v);
          }
        }
      }
    }
  }
  return op;
}

#endif  // TENES_TEST_FERMION_COMMON_HPP
