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
//! Fermion mode on skewed unit cells: covariance under unfolding.

// ===== A skewed cell is the same lattice as its unfolded skew-0 cell =======
//
// Contract: docs/superpowers/specs/2026-09-11-fermion-skew-guard-contract.md
// section 4; evidence in
// docs/superpowers/notes/2026-09-11-fermion-skew-revisit.md (experiments 2
// to 4). This is the test that keeps skewed fermion cells correct: it builds
// the iTPS directly and never goes through validate_fermion_constraints.
//
// A cell [LX, LY] with skew s (s not a multiple of LX) tiles the plane
// exactly like the cell [LX, LY * m] with skew 0, m = lcm(LX, s') / s' with
// s' = s mod LX in [1, LX). With TeNeS' convention T(x, y) = T(x + s, y + LY)
// (src/SquareLattice.hpp) the site at (x, y + k LY) holds the tensor of
// (x - k s, y), so unfolded site (x, y) holds skew-cell site
//
//     ((x - s' floor(y / LY)) mod LX,  y mod LY).
//
// That map is computed here by arithmetic alone (sku_skew_site_at), never
// through SquareLattice, so every check compares the code under test against
// geometry it did not produce.
//
// Every skew-cell bond (site, leg) has m images in the unfolded cell. Since no
// site of these cells is its own nearest neighbour, the images touch pairwise
// disjoint site tensors and bond weights, so applying a gate to all images one
// after the other leaves the unfolded state equal to the unfolding of the
// skew state after that gate. The checks:
//
//   geometry  neighbor(), other() and index() of the skewed SquareLattice
//             agree with the image map (exact, integers).
//   check 1   simple update: bit-identical Tn, bond weights and virtual
//             parity ledgers after sku_sweeps sweeps.
//   check 3   mean-field measurement: identical to the last bit, norms
//             included (it reads nothing but Tn and the bond weights).
//   check 2   CTM measurement: equal up to the finite-chi CTMRG residual.
//   check 4   one full-update bond, horizontal and vertical: equal up to the
//             same residual, and the other images stay untouched.
//
// Why checks 1 and 3 are exact comparisons (==, no ulp allowance): the two
// cells run the same arithmetic on the same numbers in the same order; only
// the storage slots differ. Observed bit-identical on HEAD with one and with
// four OpenMP threads. A BLAS whose kernels depend on the alignment of their
// operands (MKL without CNR mode) could break bit equality; that would show
// up here as differences of a few ulp, never as the O(1) differences a
// geometry error produces (with the skew branch of the neighbour map
// removed: hundreds of differing tensor elements and mean-field values off
// by up to 0.09, or a graded decomposition refusing a tensor outright).
//
// The state of checks 1 and 3 (sku_evolved_pair). Its bond ledgers must not
// all stay equal, or the ledger comparisons of check 1 cannot fail: a first
// version used D = 2, where every ledger stays the even-first [0, 1], and a
// simple update that wrote the ledger of a boundary bond to a site found by
// in-cell coordinates (ignoring the skew) passed it. Now: D = 3, a gate with
// hopping, a weak repulsion and a chemical potential off the particle-hole
// symmetric point, 10 sweeps at tau = 0.5, and a seed for which the update
// moves ledgers to [0, 1, 1] in every cell and leaves the bonds that cross
// the skewed boundary with at least two different ledgers. Both are
// asserted as premises (see sku_check_simple_update); with them, the same
// mutation fails every check-1 and check-3 case.
//
// Why check 2 is not exact, and the state it uses. At finite chi the two
// cells converge to slightly different CTM fixed points: a left move of the
// skew cell absorbs a skewed column that covers every site, so each iteration
// updates the environment of each site LX times, while the unfolded cell
// updates each column once. The difference is a truncation effect: it does
// not change with the convergence threshold and shrinks with chi (the note
// cited above, experiment 2). How large it is depends on how entangled the
// state is. Measured on HEAD at chi = 8, D = 2:
//
//   state                                  CTM residual     min |hopping|
//   random Tn (entries in [-1, 1])         7e-6 .. 6e-5     0.011
//   the same after 4 simple-update sweeps  4e-5 .. 3e-4     0.011  ([2,1]
//   (tau 0.3, V 0.8, mu 0.3)                                skew 1: no CTM
//                                                           convergence)
//   random Tn, odd virtual legs * 0.3      8e-14 .. 3e-11  4.4e-4
//
// A sign error or a swapped pair changes a value by the size of the
// observable, so what matters is observable / residual: 2e2 for the plain
// random state, 1e7 for the last one. Checks 2 and 4 therefore use the last
// one (sku_short_ranged_pair): the same random, parity-even, site-distinct
// Tn, with every element multiplied by sku_odd_scale per odd virtual index.
// Its CTM converges in 6 to 7 sweeps, which is also what keeps this file
// fast. The premises asserted on it: the CTM converged in both cells, on
// every bond the two densities differ and the hopping is non-zero, each by
// at least 1e3 times the tolerance. On the two bonds check 4 updates, the
// bond basis is reversed at both ends ([0, 1] becomes the odd-first [1, 0]),
// a change of basis, not of state; the full update returns even-first
// ledgers, so it must rewrite that bond's ledger at both ends, and a ledger
// written to the wrong slot leaves a stale one (asserted as a premise).
//
// What these checks can and cannot see. Both cells run the SAME code, so a
// mistake that is not skew-specific (say a sign on every vertical bond)
// corrupts both sides identically and passes here; the end-to-end and
// convention tests own that. What is covered is the geometry only a skewed
// cell exercises: the skewed neighbour map, index() resolving a y outside
// the cell through the skew, and the CTM moves walking a skewed column that
// is longer than the cell (LY_noskew > LY). The LX = 3 cells matter: with
// LX = 2 a skew applied with the wrong sign (x + s instead of x - s) is
// invisible, because +1 and -1 coincide mod 2, and no other test in the
// suite has a skewed cell with LX > 2.
//
// Observables. The one-site density n, and two two-site observables on every
// site and in all four directions (so both ends of every bond, horizontal and
// vertical, and every bond that wraps through the skewed boundary):
//   n (1 - n)  not symmetric under exchanging its sites, so a pair measured
//              with its ends swapped changes value (by n_i - n_j);
//   hopping    -(c^dag_1 c_2 + h.c.), non-zero (odd, odd) elements, so a
//              lost fermionic sign changes the value (by twice the value).

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../test_fermion_common.hpp"

#include <array>
#include <chrono>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace skf = tenes::fermion;
using sku_tensor = tenes::real_tensor;
using sku_state = tenes::itps::iTPS<sku_tensor>;
using sku_access = tenes::itps::iTPSTestAccessor;

constexpr int sku_d = 2;  //!< physical dimension (spinless fermion)

//! Virtual bond dimension of the state of checks 1 and 3 (the evolved state)
//! and of checks 2 and 4 (the short-ranged state), see the file comment.
constexpr int sku_D_evolved = 3;
constexpr int sku_D_short = 2;

//! Simple-update schedule of check 1 (and of the state check 3 measures):
//! sweeps over all bonds and the imaginary-time step (also the step of the
//! full-update gate of check 4). Every site tensor moves by O(1) from its
//! random start, and bond ledgers move away from even-first (both asserted).
constexpr int sku_sweeps = 10;
constexpr double sku_tau = 0.5;

//! Bond Hamiltonian of the gate: hopping, nearest-neighbour repulsion and a
//! chemical potential,
//!   h = -t (c^dag_1 c_2 + h.c.) + V n_1 n_2 - mu (n_1 + n_2).
//! mu != V / 2 keeps the gate off the particle-hole symmetric point, where
//! exactly degenerate Schmidt values could decide a ledger by rounding.
constexpr double sku_t = 1.0;
constexpr double sku_V = 0.2;
constexpr double sku_mu = 0.2;

//! Seeds of the random site tensors (site s draws from seed + 97 s). The
//! evolved state's seed is one for which, with the schedule above, the
//! ledger premises of check 1 hold in every cell, and keep holding for 9 to
//! 11 sweeps and tau from 0.49 to 0.51 (so no ledger is decided by a near
//! tie); the premises are asserted, so a seed that stops working shows up.
constexpr unsigned sku_seed_evolved = 3001;
constexpr unsigned sku_seed_short = 1009;

//! Suppression of the odd virtual components in the state of checks 2 and 4
//! (see the file comment).
constexpr double sku_odd_scale = 0.3;

//! CTM settings. D = 2 makes the reduced (doubled) legs D^2 = 4.
constexpr int sku_chi = 8;
constexpr int sku_ctm_iteration_max = 100;
constexpr double sku_ctm_epsilon = 1.0e-12;

//! Full update: a fixed number of ALS sweeps. The ALS loop stops as soon as
//! its cost changes by less than Full_Convergence_Epsilon; two cells whose
//! environments differ by the CTM residual can stop one sweep apart, which
//! would make the comparison jump from the residual to the size of one ALS
//! step on some platform and not on another. With the threshold at zero both
//! cells run exactly sku_als_sweeps sweeps.
constexpr int sku_als_sweeps = 10;

//! Tolerances, anchored on HEAD 03cd7cd6 plus this change's edits to src/
//! (none on these code paths; Debug, g++-16, macOS arm64, one OpenMP
//! thread). Largest differences observed over the four cells:
//!   CTM measurement   3.0e-11  ([3,1] skew 1, two-site)
//!   full update       1.3e-10  ([3,1] skew 1, horizontal bond, source site)
//! Both tolerances leave a factor of 30 over that. A geometry error shows up
//! at the size of the observables; the smallest one this state offers is the
//! density difference across a bond of [2,1] skew 1, 1.7e-4, i.e. 1.7e5 times
//! the tolerance. (Measured with the mutations of contract section 5: 2.5e-5
//! to 0.9 in check 2, and 3e-4 to 1.6e-3 in check 4 where the full update
//! does not refuse the broken environment outright.)
constexpr double sku_tol_ctm = 1.0e-9;
constexpr double sku_tol_full_update = 4.0e-9;

struct sku_cell {
  int lx;
  int ly;
  int skew;  //!< as written in the input; may be negative
};

std::string sku_cell_name(sku_cell c) {
  return "[" + std::to_string(c.lx) + "," + std::to_string(c.ly) + "] skew " +
         std::to_string(c.skew);
}

int sku_mod(int a, int b) { return ((a % b) + b) % b; }

int sku_floor_div(int a, int b) {
  const int q = a / b;
  return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q;
}

//! The skew-cell site that holds the tensor at global position (x, y):
//! T(x, y + k LY) = T(x - k s, y).
int sku_skew_site_at(sku_cell c, int x, int y) {
  const int k = sku_floor_div(y, c.ly);
  const int xs = sku_mod(x - c.skew * k, c.lx);
  const int ys = sku_mod(y, c.ly);
  return xs + c.lx * ys;
}

//! The skew-cell site at the other end of bond (site, leg), from the same
//! arithmetic; the fixtures use it rather than SquareLattice::neighbor(),
//! which is under test.
int sku_other_end(sku_cell c, int site, int leg) {
  const int dx[4] = {-1, 0, 1, 0};
  const int dy[4] = {0, 1, 0, -1};
  return sku_skew_site_at(c, site % c.lx + dx[leg], site / c.lx + dy[leg]);
}

//! The unfolding of a skewed cell, from arithmetic alone.
struct sku_unfolding {
  sku_cell cell{};
  int m = 0;            //!< images per skew-cell site
  int ly_unfolded = 0;  //!< LY * m
  //! preimage[u]: the skew-cell site whose tensor unfolded site u holds.
  std::vector<int> preimage;
  //! images[s]: the unfolded sites holding skew-cell site s, ascending.
  std::vector<std::vector<int>> images;
  std::string label;
};

sku_unfolding sku_unfold(sku_cell c) {
  sku_unfolding u;
  u.cell = c;
  const int s = sku_mod(c.skew, c.lx);
  REQUIRE(s != 0);
  const int lcm = c.lx / std::gcd(c.lx, s) * s;
  u.m = lcm / s;
  u.ly_unfolded = c.ly * u.m;
  u.images.assign(c.lx * c.ly, std::vector<int>{});
  for (int y = 0; y < u.ly_unfolded; ++y) {
    for (int x = 0; x < c.lx; ++x) {
      const int site = sku_skew_site_at(c, x, y);
      u.preimage.push_back(site);
      u.images[site].push_back(x + c.lx * y);
    }
  }
  for (const auto& im : u.images) {
    REQUIRE(static_cast<int>(im.size()) == u.m);
  }
  u.label = sku_cell_name(c) + " -> [" + std::to_string(c.lx) + "," +
            std::to_string(u.ly_unfolded) + "] skew 0";
  return u;
}

//! The cells of this file. [2,1] skew 1 is what tenes_simple builds for a
//! square lattice with W = 1; [3,1] skew 1 and [3,1] skew -1 are LX = 3 cells
//! (where the sign of the skew matters; the negative one also exercises the
//! C++ member keeping the sign of the input, here -1 % 3 = -1); [2,2] skew 1
//! has both sides >= 2.
const std::vector<sku_cell>& sku_cells() {
  static const std::vector<sku_cell> cells{
      {2, 1, 1}, {3, 1, 1}, {3, 1, -1}, {2, 2, 1}};
  return cells;
}

tenes::SquareLattice sku_lattice(int lx, int ly, int skew, int D) {
  tenes::SquareLattice lattice(lx, ly, skew);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    lattice.physical_dims[site] = sku_d;
    lattice.virtual_dims[site] = {D, D, D, D};
    lattice.initial_dirs[site] = {0.0};
    lattice.noises[site] = 1.0;
  }
  return lattice;
}

tenes::itps::PEPS_Parameters sku_params(int n_unit, bool meanfield) {
  tenes::itps::PEPS_Parameters p;
  p.fermion = true;
  p.is_real = true;
  p.phys_parity.assign(n_unit, skf::parity_vector{false, true});
  p.print_level = tenes::PrintLevel::none;
  p.outdir = "output_test_fermion_skew_unfold";
  p.CHI = sku_chi;
  p.Max_CTM_Iteration = sku_ctm_iteration_max;
  p.CTM_Convergence_Epsilon = sku_ctm_epsilon;
  p.Use_RSVD = false;
  p.MeanField_Env = meanfield;
  p.Full_max_iteration = sku_als_sweeps;
  p.Full_Convergence_Epsilon = 0.0;
  return p;
}

// ---- operators --------------------------------------------------------------
// Two-site tensors are indexed (in_1, in_2, out_1, out_2), matching
// electron_gate() in test_fermion_common.hpp. In the basis |n_1 n_2> with
// site 1 first in the Jordan-Wigner order, <10| c^dag_1 c_2 |01> = +1.

//! exp(-tau h) in closed form: 1 on |00>, e^{tau mu} [[cosh, sinh],
//! [sinh, cosh]](tau t) on {|01>, |10>}, e^{-tau (V - 2 mu)} on |11>.
sku_tensor sku_gate(double tau) {
  sku_tensor g(mptensor::Shape(2, 2, 2, 2));
  const double one = std::exp(tau * sku_mu);
  g.set_value(mptensor::Index(0, 0, 0, 0), 1.0);
  g.set_value(mptensor::Index(0, 1, 0, 1), one * std::cosh(tau * sku_t));
  g.set_value(mptensor::Index(1, 0, 1, 0), one * std::cosh(tau * sku_t));
  g.set_value(mptensor::Index(0, 1, 1, 0), one * std::sinh(tau * sku_t));
  g.set_value(mptensor::Index(1, 0, 0, 1), one * std::sinh(tau * sku_t));
  g.set_value(mptensor::Index(1, 1, 1, 1),
              std::exp(-tau * (sku_V - 2.0 * sku_mu)));
  return g;
}

sku_tensor sku_number() {
  sku_tensor n(mptensor::Shape(2, 2));
  n.set_value(mptensor::Index(1, 1), 1.0);
  return n;
}

//! n (1 - n): 1 on |10> only.
sku_tensor sku_particle_hole() {
  sku_tensor o(mptensor::Shape(2, 2, 2, 2));
  o.set_value(mptensor::Index(1, 0, 1, 0), 1.0);
  return o;
}

//! -(c^dag_1 c_2 + c^dag_2 c_1).
sku_tensor sku_hopping() {
  sku_tensor o(mptensor::Shape(2, 2, 2, 2));
  o.set_value(mptensor::Index(0, 1, 1, 0), -1.0);
  o.set_value(mptensor::Index(1, 0, 0, 1), -1.0);
  return o;
}

constexpr int sku_group_n_1mn = 0;
constexpr int sku_group_hop = 1;
constexpr int sku_directions[4][2] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};

struct sku_observables {
  tenes::Operators<sku_tensor> onesite;
  tenes::Operators<sku_tensor> twosite;
};

sku_observables sku_make_observables(int n_unit) {
  sku_observables obs;
  const sku_tensor n = sku_number();
  const sku_tensor ph = sku_particle_hole();
  const sku_tensor hop = sku_hopping();
  for (int site = 0; site < n_unit; ++site) {
    obs.onesite.emplace_back("n", 0, site, n);
    for (const auto& d : sku_directions) {
      obs.twosite.emplace_back("n1mn", sku_group_n_1mn, site, d[0], d[1], ph);
      obs.twosite.emplace_back("hop", sku_group_hop, site, d[0], d[1], hop);
    }
  }
  return obs;
}

std::unique_ptr<sku_state> sku_make_state(const tenes::SquareLattice& lattice,
                                          bool meanfield) {
  const sku_observables obs = sku_make_observables(lattice.N_UNIT);
  return std::make_unique<sku_state>(
      MPI_COMM_WORLD, sku_params(lattice.N_UNIT, meanfield), lattice,
      tenes::EvolutionOperators<sku_tensor>{},
      tenes::EvolutionOperators<sku_tensor>{}, obs.onesite, obs.twosite,
      tenes::Operators<sku_tensor>{}, tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
}

// ---- tensor comparisons -----------------------------------------------------

struct sku_diff {
  double max_abs = 0.0;
  std::size_t n_differ = 0;
};

sku_diff sku_compare(const sku_tensor& a, const sku_tensor& b) {
  REQUIRE(a.rank() == b.rank());
  for (std::size_t ax = 0; ax < a.rank(); ++ax) {
    REQUIRE(a.shape()[ax] == b.shape()[ax]);
  }
  REQUIRE(a.local_size() == b.local_size());
  sku_diff d;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    if (a[n] != b[n]) {
      ++d.n_differ;
      d.max_abs = std::max(d.max_abs, std::abs(a[n] - b[n]));
    }
  }
  return d;
}

double sku_max_abs(const sku_tensor& a) {
  double m = 0.0;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    m = std::max(m, std::abs(a[n]));
  }
  return m;
}

double sku_seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
      .count();
}

// ---- the pair of states -----------------------------------------------------

struct sku_pair {
  sku_unfolding map;
  tenes::SquareLattice skew_lattice;
  tenes::SquareLattice unfolded_lattice;
  std::unique_ptr<sku_state> skew;
  std::unique_ptr<sku_state> unfolded;
  std::vector<sku_tensor> skew_initial;  //!< Tn of the skew cell as seeded

  sku_pair(sku_cell c, bool meanfield, int D)
      : map(sku_unfold(c)),
        skew_lattice(sku_lattice(c.lx, c.ly, c.skew, D)),
        unfolded_lattice(sku_lattice(c.lx, map.ly_unfolded, 0, D)),
        skew(sku_make_state(skew_lattice, meanfield)),
        unfolded(sku_make_state(unfolded_lattice, meanfield)) {}
};

//! Random, parity-even, site-distinct Tn and site-distinct bond weights on
//! the skew cell, copied onto every image in the unfolded cell together with
//! the virtual ledgers. Every element is multiplied by odd_scale per odd
//! virtual index it carries.
void sku_seed(sku_pair& p, double odd_scale, unsigned seed) {
  auto& Ts = sku_access::Tn(*p.skew);
  auto& Tu = sku_access::Tn(*p.unfolded);
  auto& Ls = sku_access::lambda_tensor(*p.skew);
  auto& Lu = sku_access::lambda_tensor(*p.unfolded);
  const auto& fs = sku_access::finfo(*p.skew);
  auto& fu = sku_access::finfo(*p.unfolded);
  const int n_skew = p.skew_lattice.N_UNIT;
  const int D = p.skew_lattice.virtual_dims[0][0];
  for (int s = 0; s < n_skew; ++s) {
    const skf::leg_parities parity = skf::Tn_parity(fs, s);
    sku_tensor t(mptensor::Shape(D, D, D, D, sku_d));
    std::mt19937 gen(seed + 97u * static_cast<unsigned>(s));
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (std::size_t n = 0; n < t.local_size(); ++n) {
      const mptensor::Index idx = t.global_index(n);
      double v = dist(gen);
      if (skf::count_odd(parity, idx) % 2 != 0) {
        v = 0.0;
      }
      for (int leg = 0; leg < 4; ++leg) {
        if (parity[leg][idx[leg]]) {
          v *= odd_scale;
        }
      }
      t.set_value(idx, v);
    }
    REQUIRE(skf::parity_violation(skf::ftensor<sku_tensor>{t, parity}) == 0.0);
    Ts[s] = t;
  }
  // One weight vector per skew-cell bond, written to both of its ends.
  int bond = 0;
  for (int s = 0; s < n_skew; ++s) {
    for (int leg : {2, 1}) {
      const double w = 0.35 + 0.07 * bond++;
      std::vector<double> weights(D, 1.0);
      for (int k = 1; k < D; ++k) {
        weights[k] = weights[k - 1] * w;
      }
      Ls[s][leg] = weights;
      Ls[sku_other_end(p.map.cell, s, leg)][(leg + 2) % 4] = weights;
    }
  }
  for (int u = 0; u < p.unfolded_lattice.N_UNIT; ++u) {
    Tu[u] = Ts[p.map.preimage[u]];
    Lu[u] = Ls[p.map.preimage[u]];
    fu.virt[u] = fs.virt[p.map.preimage[u]];
  }
  // The unfolded cell has skew 0, so its neighbour map is not what this
  // file tests: a premise. On the skewed cell the same consistency is a
  // property of its neighbour map, so a CHECK, and the checks run on.
  REQUIRE_NOTHROW(skf::validate_neighbor_consistency(fu, p.unfolded_lattice));
  CHECK_NOTHROW(skf::validate_neighbor_consistency(fs, p.skew_lattice));
  p.skew_initial = Ts;
}

//! The gate list of the skew cell, one gate per bond, in sweep order. Leg 1
//! names a vertical bond from its lower (raster-later) end, which makes the
//! solver swap the two sites; leg 2 names a horizontal bond from its left
//! end, which it does not. Both branches run.
std::vector<std::pair<int, int>> sku_bond_order(int n_skew) {
  std::vector<std::pair<int, int>> order;
  for (int leg : {2, 1}) {
    for (int s = 0; s < n_skew; ++s) {
      order.emplace_back(s, leg);
    }
  }
  return order;
}

//! Check 1's state (and check 3's): full-amplitude random Tn evolved by
//! sku_sweeps simple-update sweeps; the unfolded cell applies every gate to
//! all its images in a row.
std::unique_ptr<sku_pair> sku_evolved_pair(sku_cell c, bool meanfield) {
  auto p = std::make_unique<sku_pair>(c, meanfield, sku_D_evolved);
  sku_seed(*p, 1.0, sku_seed_evolved);
  const sku_tensor gate = sku_gate(sku_tau);
  for (int sweep = 0; sweep < sku_sweeps; ++sweep) {
    for (const auto& bond : sku_bond_order(p->skew_lattice.N_UNIT)) {
      const int s = bond.first;
      const int leg = bond.second;
      p->skew->simple_update(
          tenes::make_twosite_EvolutionOperator(s, leg, 0, gate));
      for (int u : p->map.images[s]) {
        p->unfolded->simple_update(
            tenes::make_twosite_EvolutionOperator(u, leg, 0, gate));
      }
    }
  }
  return p;
}

//! The two bonds of check 4, as (skew-cell site, leg). Horizontal: the bond
//! that wraps periodically in x. Vertical: a bond of the top row, which
//! wraps through the skewed boundary (and is an inner bond of the unfolded
//! cell).
std::array<std::pair<int, int>, 2> sku_full_update_bonds(sku_cell c) {
  return {{{c.lx - 1, 2}, {c.lx * (c.ly - 1), 1}}};
}

//! Reverses the basis of one virtual leg of a site tensor (index k becomes
//! D - 1 - k).
sku_tensor sku_reverse_leg(const sku_tensor& t, int leg) {
  sku_tensor out(t.shape());
  const std::size_t dim = t.shape()[leg];
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    mptensor::Index idx = t.global_index(n);
    const double v = t[n];
    idx[leg] = dim - 1 - idx[leg];
    out.set_value(idx, v);
  }
  return out;
}

//! Checks 2 and 4's state: random Tn with the odd virtual components
//! suppressed, so that its CTM converges fast and the finite-chi residual is
//! small next to the observables.
//!
//! On the two bonds check 4 updates, the bond basis is then reversed at both
//! ends, together with its ledger and weights: [0, 1] becomes the odd-first
//! [1, 0]. That is a change of basis on the bond, not of the state. A full
//! update returns its bond ledger in even-first order, so it necessarily
//! rewrites the ledger at both ends, and a ledger written to the wrong slot
//! leaves a stale one behind (asserted as a premise in check 4); with every
//! ledger the even-first [0, 1] such a write would change nothing.
std::unique_ptr<sku_pair> sku_short_ranged_pair(sku_cell c) {
  auto p = std::make_unique<sku_pair>(c, false, sku_D_short);
  sku_seed(*p, sku_odd_scale, sku_seed_short);
  auto& Ts = sku_access::Tn(*p->skew);
  auto& Ls = sku_access::lambda_tensor(*p->skew);
  auto& fs = sku_access::finfo(*p->skew);
  for (const auto& bond : sku_full_update_bonds(c)) {
    const int ends[2][2] = {
        {bond.first, bond.second},
        {sku_other_end(c, bond.first, bond.second), (bond.second + 2) % 4}};
    for (const auto& end : ends) {
      const int s = end[0];
      const int leg = end[1];
      Ts[s] = sku_reverse_leg(Ts[s], leg);
      std::reverse(Ls[s][leg].begin(), Ls[s][leg].end());
      std::reverse(fs.virt[s][leg].begin(), fs.virt[s][leg].end());
    }
  }
  auto& Tu = sku_access::Tn(*p->unfolded);
  auto& Lu = sku_access::lambda_tensor(*p->unfolded);
  auto& fu = sku_access::finfo(*p->unfolded);
  for (int u = 0; u < p->unfolded_lattice.N_UNIT; ++u) {
    Tu[u] = Ts[p->map.preimage[u]];
    Lu[u] = Ls[p->map.preimage[u]];
    fu.virt[u] = fs.virt[p->map.preimage[u]];
  }
  for (int s = 0; s < p->skew_lattice.N_UNIT; ++s) {
    REQUIRE(skf::parity_violation(
                skf::ftensor<sku_tensor>{Ts[s], skf::Tn_parity(fs, s)}) == 0.0);
  }
  REQUIRE_NOTHROW(skf::validate_neighbor_consistency(fu, p->unfolded_lattice));
  CHECK_NOTHROW(skf::validate_neighbor_consistency(fs, p->skew_lattice));
  p->skew_initial = Ts;
  return p;
}

//! The fermionic branch of iTPS::update_CTM() (a cold start), returning the
//! number of CTM sweeps so that convergence can be asserted.
int sku_converge_ctm(sku_state& state) {
  const std::vector<sku_tensor> reduced = skf::build_reduced_density_tensors(
      sku_access::Tn(state), sku_access::finfo(state));
  return tenes::itps::core::Calc_CTM_Environment_density(
      sku_access::C1(state), sku_access::C2(state), sku_access::C3(state),
      sku_access::C4(state), sku_access::eTt(state), sku_access::eTr(state),
      sku_access::eTb(state), sku_access::eTl(state), reduced,
      sku_access::peps_parameters(state), sku_access::lattice(state), true,
      true);
}

//! Largest difference between every unfolded measurement and the one of the
//! skew-cell site/bond it images. `include_norms` also compares the norm
//! rows, which only the mean-field path fixes independently of the cell (a
//! CTM normalises its environment per cell, and the norm is not an
//! observable).
struct sku_measure_diff {
  double onesite = 0.0;
  double twosite = 0.0;
  std::size_t n_compared = 0;
};

sku_measure_diff sku_compare_measurements(sku_pair& p, bool include_norms) {
  const auto one_s = p.skew->measure_onesite();
  const auto one_u = p.unfolded->measure_onesite();
  const auto two_s = p.skew->measure_twosite();
  const auto two_u = p.unfolded->measure_twosite();
  REQUIRE(one_s.size() == 2);  // n, norms
  REQUIRE(one_u.size() == 2);
  REQUIRE(two_s.size() == 3);  // n(1-n), hopping, norms
  REQUIRE(two_u.size() == 3);
  sku_measure_diff d;
  const std::size_t n_one = include_norms ? 2 : 1;
  const std::size_t n_two = include_norms ? 3 : 2;
  for (std::size_t g = 0; g < n_one; ++g) {
    REQUIRE(one_u[g].size() ==
            static_cast<std::size_t>(p.unfolded_lattice.N_UNIT));
    for (int u = 0; u < p.unfolded_lattice.N_UNIT; ++u) {
      const double a = one_u[g][u];
      const double b = one_s[g][p.map.preimage[u]];
      INFO(p.map.label << " one-site row " << g << " unfolded site " << u);
      REQUIRE(std::isfinite(a));
      REQUIRE(std::isfinite(b));
      d.onesite = std::max(d.onesite, std::abs(a - b));
      ++d.n_compared;
    }
  }
  for (std::size_t g = 0; g < n_two; ++g) {
    REQUIRE(two_u[g].size() == two_s[g].size() * p.map.m);
    for (const auto& entry : two_u[g]) {
      const tenes::itps::Bond bond = entry.first;
      const double a = entry.second;
      const tenes::itps::Bond image{p.map.preimage[bond.source_site], bond.dx,
                                    bond.dy};
      INFO(p.map.label << " two-site row " << g << " bond (" << bond.source_site
                       << ", " << bond.dx << ", " << bond.dy << ")");
      REQUIRE(two_s[g].count(image) == 1);
      const double b = two_s[g].at(image);
      REQUIRE(std::isfinite(a));
      REQUIRE(std::isfinite(b));
      d.twosite = std::max(d.twosite, std::abs(a - b));
      ++d.n_compared;
    }
  }
  return d;
}

}  // namespace

// ---- geometry ---------------------------------------------------------------

TEST_CASE("fermion skew unfolding: the skewed lattice maps are the image map") {
  for (const sku_cell c : sku_cells()) {
    const sku_unfolding map = sku_unfold(c);
    INFO(map.label);
    const tenes::SquareLattice skew = sku_lattice(c.lx, c.ly, c.skew, 2);
    const tenes::SquareLattice unfolded =
        sku_lattice(c.lx, map.ly_unfolded, 0, 2);
    CHECK(skew.LY_noskew == map.ly_unfolded);
    // neighbor(): the neighbour of an image is an image of the neighbour.
    for (int u = 0; u < unfolded.N_UNIT; ++u) {
      for (int leg = 0; leg < 4; ++leg) {
        INFO("unfolded site " << u << " leg " << leg);
        CHECK(map.preimage[unfolded.neighbor(u, leg)] ==
              skew.neighbor(map.preimage[u], leg));
      }
    }
    // other(): displacements up to two steps from every skew-cell site.
    for (int s = 0; s < skew.N_UNIT; ++s) {
      for (int dx = -2; dx <= 2; ++dx) {
        for (int dy = -2; dy <= 2; ++dy) {
          INFO("skew site " << s << " displacement (" << dx << ", " << dy
                            << ")");
          CHECK(skew.other(s, dx, dy) ==
                sku_skew_site_at(c, skew.x(s) + dx, skew.y(s) + dy));
        }
      }
    }
    // index(): coordinates outside the cell, as the CTM moves use them
    // (y up to LY_noskew - 1) and beyond.
    for (int y = -2 * map.ly_unfolded; y < 2 * map.ly_unfolded; ++y) {
      for (int x = -c.lx; x < 2 * c.lx; ++x) {
        INFO("index(" << x << ", " << y << ")");
        CHECK(skew.index(x, y) == sku_skew_site_at(c, x, y));
      }
    }
  }
}

// ---- check 1: simple update -------------------------------------------------

namespace {

//! Renders a parity ledger as its 0/1 string, e.g. "011".
std::string sku_ledger_string(const skf::parity_vector& ledger) {
  std::string out;
  for (const bool odd : ledger) {
    out += odd ? '1' : '0';
  }
  return out;
}

void sku_check_simple_update(sku_cell c) {
  const auto t0 = std::chrono::steady_clock::now();
  auto pp = sku_evolved_pair(c, false);
  sku_pair& p = *pp;
  INFO(p.map.label);
  auto& Ts = sku_access::Tn(*p.skew);
  auto& Tu = sku_access::Tn(*p.unfolded);
  auto& Ls = sku_access::lambda_tensor(*p.skew);
  auto& Lu = sku_access::lambda_tensor(*p.unfolded);
  auto& fs = sku_access::finfo(*p.skew);
  auto& fu = sku_access::finfo(*p.unfolded);

  // Premises: the update did something - every skew-cell tensor moved by a
  // sizeable fraction of its scale - and the sites are distinct, so an image
  // holding the wrong site cannot go unnoticed.
  double min_change = std::numeric_limits<double>::infinity();
  for (int s = 0; s < p.skew_lattice.N_UNIT; ++s) {
    const double scale =
        std::max(sku_max_abs(Ts[s]), sku_max_abs(p.skew_initial[s]));
    REQUIRE(scale > 0.0);
    min_change = std::min(
        min_change, sku_compare(Ts[s], p.skew_initial[s]).max_abs / scale);
  }
  double min_distinct = std::numeric_limits<double>::infinity();
  for (int s = 0; s < p.skew_lattice.N_UNIT; ++s) {
    for (int r = s + 1; r < p.skew_lattice.N_UNIT; ++r) {
      min_distinct = std::min(min_distinct, sku_compare(Ts[s], Ts[r]).max_abs);
    }
  }
  // Premises of the ledger comparison. Every ledger starts even-first; if
  // they all stayed so, a ledger written to the wrong slot would go
  // unnoticed. (a) Some bond ledger has moved away from even-first. (b) The
  // bonds that cross the skewed boundary (vertical bonds of the top row) do
  // not all carry the same ledger: a site found by in-cell coordinates,
  // ignoring the skew, is the bottom end of another boundary bond, so that
  // is the mix-up (b) makes visible.
  const skf::parity_vector even_first = skf::even_first_parity(sku_D_evolved);
  int n_bonds_moved = 0;
  std::set<std::string> boundary_ledgers;
  std::string ledgers;
  for (int s = 0; s < p.skew_lattice.N_UNIT; ++s) {
    for (const int leg : {2, 1}) {
      const std::string ledger = sku_ledger_string(fs.virt[s][leg]);
      n_bonds_moved += fs.virt[s][leg] != even_first ? 1 : 0;
      if (leg == 1 && p.skew_lattice.y(s) == c.ly - 1) {
        boundary_ledgers.insert(ledger);
      }
      ledgers +=
          " (" + std::to_string(s) + "," + std::to_string(leg) + ")=" + ledger;
    }
  }
  std::cout << std::setprecision(3) << "skew unfolding SU " << p.map.label
            << ": min relative change of Tn " << min_change
            << ", min distance between sites " << min_distinct
            << ", bond ledgers (site,leg)=ledger" << ledgers << " ("
            << sku_seconds_since(t0) << " s)" << std::endl;
  CHECK(min_change > 0.1);
  CHECK(min_distinct > 0.1);
  CHECK(n_bonds_moved >= 1);
  CHECK(boundary_ledgers.size() >= 2);

  CHECK_NOTHROW(skf::validate_neighbor_consistency(fs, p.skew_lattice));
  CHECK_NOTHROW(skf::validate_neighbor_consistency(fu, p.unfolded_lattice));
  std::size_t n_tensor_mismatch = 0;
  std::size_t n_weight_mismatch = 0;
  std::size_t n_ledger_mismatch = 0;
  for (int u = 0; u < p.unfolded_lattice.N_UNIT; ++u) {
    const int s = p.map.preimage[u];
    INFO("unfolded site " << u << " images skew site " << s);
    const sku_diff d = sku_compare(Tu[u], Ts[s]);
    CHECK(d.n_differ == 0);
    CHECK(d.max_abs == 0.0);
    n_tensor_mismatch += d.n_differ;
    for (int leg = 0; leg < 4; ++leg) {
      INFO("leg " << leg);
      CHECK(Lu[u][leg] == Ls[s][leg]);
      CHECK(fu.virt[u][leg] == fs.virt[s][leg]);
      n_weight_mismatch += Lu[u][leg] == Ls[s][leg] ? 0 : 1;
      n_ledger_mismatch += fu.virt[u][leg] == fs.virt[s][leg] ? 0 : 1;
    }
  }
  std::cout << "skew unfolding SU " << p.map.label
            << ": mismatching Tn elements " << n_tensor_mismatch
            << ", bond weights " << n_weight_mismatch << ", ledgers "
            << n_ledger_mismatch << std::endl;
}

// ---- check 3: mean-field measurement ----------------------------------------

void sku_check_mean_field(sku_cell c) {
  const auto t0 = std::chrono::steady_clock::now();
  auto pp = sku_evolved_pair(c, true);
  sku_pair& p = *pp;
  INFO(p.map.label);
  const sku_measure_diff d = sku_compare_measurements(p, true);
  std::cout << std::setprecision(3) << "skew unfolding MF " << p.map.label
            << ": max |diff| one-site " << d.onesite << ", two-site "
            << d.twosite << " over " << d.n_compared << " values ("
            << sku_seconds_since(t0) << " s)" << std::endl;
  CHECK(d.onesite == 0.0);
  CHECK(d.twosite == 0.0);
}

}  // namespace

// One test case per cell, so that an exception in one cell (a graded
// decomposition refusing a tensor whose ledger went stale, say) does not
// hide what the other cells do.

TEST_CASE("fermion skew unfolding: the simple update is exact, [2,1] skew 1") {
  sku_check_simple_update({2, 1, 1});
}

TEST_CASE("fermion skew unfolding: the simple update is exact, [3,1] skew 1") {
  sku_check_simple_update({3, 1, 1});
}

TEST_CASE("fermion skew unfolding: the simple update is exact, [3,1] skew -1") {
  sku_check_simple_update({3, 1, -1});
}

TEST_CASE("fermion skew unfolding: the simple update is exact, [2,2] skew 1") {
  sku_check_simple_update({2, 2, 1});
}

TEST_CASE("fermion skew unfolding: mean-field values are exact, [2,1] skew 1") {
  sku_check_mean_field({2, 1, 1});
}

TEST_CASE("fermion skew unfolding: mean-field values are exact, [3,1] skew 1") {
  sku_check_mean_field({3, 1, 1});
}

TEST_CASE(
    "fermion skew unfolding: mean-field values are exact, [3,1] skew -1") {
  sku_check_mean_field({3, 1, -1});
}

TEST_CASE("fermion skew unfolding: mean-field values are exact, [2,2] skew 1") {
  sku_check_mean_field({2, 2, 1});
}

// ---- checks 2 and 4: CTM measurement and one full-update bond ---------------

namespace {

void sku_check_ctm_and_full_update(sku_cell c) {
  const auto t0 = std::chrono::steady_clock::now();
  auto pp = sku_short_ranged_pair(c);
  sku_pair& p = *pp;
  INFO(p.map.label);

  // Premise: both environments converged. A CHECK rather than a REQUIRE, so
  // that a skew cell whose CTM no longer converges still has its values
  // compared below (they are the more informative failure).
  const int count_skew = sku_converge_ctm(*p.skew);
  const int count_unfolded = sku_converge_ctm(*p.unfolded);
  {
    INFO("premise: CTM sweeps " << count_skew << " (skew) and "
                                << count_unfolded << " (unfolded) of at most "
                                << sku_ctm_iteration_max);
    CHECK(count_skew < sku_ctm_iteration_max);
    CHECK(count_unfolded < sku_ctm_iteration_max);
  }
  const double t_ctm = sku_seconds_since(t0);

  // Premises of the two-site observables, on the unfolded cell (which does
  // not depend on the skew geometry): on every bond the two densities differ
  // (so n (1 - n) sees a swapped pair) and the hopping is not negligible (so
  // it sees a lost sign), both far above the tolerance.
  {
    const auto one = p.unfolded->measure_onesite();
    const auto two = p.unfolded->measure_twosite();
    double min_density_gap = std::numeric_limits<double>::infinity();
    double min_hopping = std::numeric_limits<double>::infinity();
    double min_density = std::numeric_limits<double>::infinity();
    double max_density = 0.0;
    for (const double n : one[0]) {
      min_density = std::min(min_density, n);
      max_density = std::max(max_density, n);
    }
    for (const auto& entry : two[sku_group_hop]) {
      const tenes::itps::Bond bond = entry.first;
      const double v = entry.second;
      const int other =
          p.unfolded_lattice.other(bond.source_site, bond.dx, bond.dy);
      min_density_gap = std::min(
          min_density_gap, std::abs(one[0][bond.source_site] - one[0][other]));
      min_hopping = std::min(min_hopping, std::abs(v));
    }
    std::cout << std::setprecision(3) << "skew unfolding CTM " << p.map.label
              << ": premise densities in [" << min_density << ", "
              << max_density << "], min density gap across a bond "
              << min_density_gap << ", min |hopping| " << min_hopping
              << std::endl;
    CHECK(min_density_gap > 1.0e3 * sku_tol_ctm);
    CHECK(min_hopping > 1.0e3 * sku_tol_ctm);
  }

  // Check 2.
  const sku_measure_diff d = sku_compare_measurements(p, false);
  std::cout << std::setprecision(3) << "skew unfolding CTM " << p.map.label
            << ": CTM sweeps " << count_skew << " / " << count_unfolded
            << ", max |diff| one-site " << d.onesite << ", two-site "
            << d.twosite << " over " << d.n_compared << " values (" << t_ctm
            << " s CTM, " << sku_seconds_since(t0) << " s so far)" << std::endl;
  CHECK(d.onesite <= sku_tol_ctm);
  CHECK(d.twosite <= sku_tol_ctm);

  // Check 4: one full-update bond, on the skew cell and on the FIRST image
  // in the unfolded cell, each from its converged environment above.
  const sku_tensor gate = sku_gate(sku_tau);
  for (const auto& fu_bond : sku_full_update_bonds(c)) {
    const int s = fu_bond.first;
    const int leg = fu_bond.second;
    INFO("full update on skew bond (" << s << ", leg " << leg << ")");
    sku_state skew_fu = *p.skew;
    sku_state unfolded_fu = *p.unfolded;
    const int s_nb = p.skew_lattice.neighbor(s, leg);
    const int u = p.map.images[s][0];
    const int u_nb = p.unfolded_lattice.neighbor(u, leg);
    CHECK(p.map.preimage[u_nb] == s_nb);

    const std::vector<sku_tensor> before_u = sku_access::Tn(unfolded_fu);
    const skf::parity_vector ledger_before =
        sku_access::finfo(skew_fu).virt[s][leg];
    skew_fu.full_update(tenes::make_twosite_EvolutionOperator(s, leg, 0, gate));
    unfolded_fu.full_update(
        tenes::make_twosite_EvolutionOperator(u, leg, 0, gate));
    const auto& Ts = sku_access::Tn(skew_fu);
    const auto& Tu = sku_access::Tn(unfolded_fu);

    // Premise: the bond really changed. Measured on the two site tensors
    // contracted over the updated bond, which does not depend on the basis
    // or gauge the update leaves on that bond.
    const auto pair_product = [&](const std::vector<sku_tensor>& T) {
      return mptensor::tensordot(T[s], T[s_nb], mptensor::Axes(leg),
                                 mptensor::Axes((leg + 2) % 4));
    };
    const sku_tensor product_after = pair_product(Ts);
    const double change =
        sku_compare(product_after, pair_product(sku_access::Tn(*p.skew)))
            .max_abs /
        sku_max_abs(product_after);
    const double diff_source = sku_compare(Tu[u], Ts[s]).max_abs;
    const double diff_target = sku_compare(Tu[u_nb], Ts[s_nb]).max_abs;
    const skf::parity_vector ledger_after =
        sku_access::finfo(skew_fu).virt[s][leg];
    std::cout << std::setprecision(3) << "skew unfolding FU " << p.map.label
              << " bond (" << s << ", leg " << leg
              << "): relative change of the pair " << change << ", ledger "
              << sku_ledger_string(ledger_before) << " -> "
              << sku_ledger_string(ledger_after) << ", max |diff| source "
              << diff_source << ", target " << diff_target << std::endl;
    CHECK(change > 1.0e-2);
    // Premise of the ledger comparison below: the update rewrote the bond's
    // ledger, so a ledger written to the wrong slot leaves a stale one.
    CHECK(ledger_after != ledger_before);
    CHECK(diff_source <= sku_tol_full_update);
    CHECK(diff_target <= sku_tol_full_update);
    for (int l = 0; l < 4; ++l) {
      INFO("ledger leg " << l);
      CHECK(sku_access::finfo(unfolded_fu).virt[u][l] ==
            sku_access::finfo(skew_fu).virt[s][l]);
      CHECK(sku_access::finfo(unfolded_fu).virt[u_nb][l] ==
            sku_access::finfo(skew_fu).virt[s_nb][l]);
    }
    // The other images of the bond are untouched.
    for (std::size_t k = 1; k < p.map.images[s].size(); ++k) {
      const int v = p.map.images[s][k];
      const int v_nb = p.unfolded_lattice.neighbor(v, leg);
      INFO("other image (" << v << ", " << v_nb << ")");
      CHECK(sku_compare(Tu[v], before_u[v]).n_differ == 0);
      CHECK(sku_compare(Tu[v_nb], before_u[v_nb]).n_differ == 0);
    }
  }
  std::cout << "skew unfolding CTM+FU " << p.map.label << ": "
            << sku_seconds_since(t0) << " s" << std::endl;
}

}  // namespace

TEST_CASE("fermion skew unfolding: CTM and full update, [2,1] skew 1") {
  sku_check_ctm_and_full_update({2, 1, 1});
}

TEST_CASE("fermion skew unfolding: CTM and full update, [3,1] skew 1") {
  sku_check_ctm_and_full_update({3, 1, 1});
}

TEST_CASE("fermion skew unfolding: CTM and full update, [3,1] skew -1") {
  sku_check_ctm_and_full_update({3, 1, -1});
}

TEST_CASE("fermion skew unfolding: CTM and full update, [2,2] skew 1") {
  sku_check_ctm_and_full_update({2, 2, 1});
}
