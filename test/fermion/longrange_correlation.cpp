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
//! Task T3 of docs/superpowers/plans/2026-09-25-fermion-longrange-measure.md:
//! the correlation function (r_max > 0) in fermion mode with the CTM
//! environment (design
//! docs/superpowers/specs/2026-09-25-fermion-longrange-measure-design.md,
//! sections 3.3, 3.4 and 4.2). Behaviour contract:
//! work/fermion-longrange/t3/contract.md, items 1 to 8.
//!
//! Conventions fixed for every case in this file:
//!   - One-site operators use TeNeS' layout op[in, out] = <out|A|in>; two-site
//!     operators op[in_s, in_t, out_s, out_t], source (left / lower) first.
//!   - d = 2 basis {|0>, |1>}, ledger [e, o]. d = 4 basis {|0>, |up>, |dn>,
//!     |updn> = c+_up c+_dn |0>}, ledger [e, o, o, e].
//!   - A correlation row {left, (r, 0)} holds <A_left B_right^r(left)> and a
//!     row {left, (0, r)} holds <A_left B_top^r(left)>, as in the bosonic
//!     measure_correlation().
//!
//! Truth sources (never measure_correlation() nor the r_max guards):
//!   - the window measurement measure_twosite() of the two-site operator
//!     product_twosite_op(A, B) at (r, 0) and (0, r), r <= 3 (task T2,
//!     verified before this task; contract item 4 names it);
//!   - the existing nearest-neighbour path of measure_twosite() (bundled-k)
//!     for r = 1, with the product table written out by the test itself
//!     (contract item 5);
//!   - the "direct chain": T1's relay site tensors (build_relay_site(), with
//!     the roles and legs build_relay_window() gives a straight path) on the
//!     row (column) from the left (lower) site to the right (upper) one,
//!     closed with the solver's environment through the density correlation
//!     kernels StartCorrelation / Transfer / FinishCorrelation_density_CTM,
//!     as design section 3.4 describes, assembled here. It is pinned to the
//!     window measurement at r <= 3 by a [truth] case and is the reference
//!     at r = 4, 5, where no window exists (contract item 6);
//!   - for a pair of different parity, the exact 0.0 of the contract; a
//!     premise shows that the plain contraction of such a pair is far from 0
//!     on the states used, so the check is not satisfied by accident.
//!
//! Every state is rank-consistent: the tensor entries are a hash of the
//! global index (no random stream over local elements), and the environment
//! noise is scaled by a global maximum, so every case means the same at any
//! rank count. The case T3-8 is registered a second time at two MPI ranks
//! (contract item 8, test/CMakeLists.txt).
//!
//! The window measurement runs on a solver of its own without [correlation]
//! (r_max = 0), with the same Tn and an element-wise copy of the same
//! environment as the solver that measures the correlation function: the
//! reference never passes through the r_max guards or the T3 code.
//!
//! Cases whose name carries [truth] check the reference side only, and cases
//! whose name carries [kept] pin a rejection that the stub already performs;
//! both are expected to pass against the stub. Every other case must fail
//! until T3 is implemented.

#include "../test_fermion_common.hpp"

#include <algorithm>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "../../src/fermion/relay.hpp"
#include "../../src/iTPS/load_toml.hpp"
#include "../../src/mpi.hpp"

namespace {

namespace rf = tenes::fermion;
using lc_acc = tenes::itps::iTPSTestAccessor;
using tenes::complex_tensor;
using tenes::real_tensor;
using tenes::itps::Bond;
using tenes::itps::Correlation;
using tenes::itps::CorrelationParameter;
using lc_cplx = std::complex<double>;

template <class tensor>
using lc_state = tenes::itps::iTPS<tensor>;

// ---- small helpers ---------------------------------------------------------

template <class tensor>
const char* lc_type_name() {
  return std::is_same<tensor, complex_tensor>::value ? "complex" : "real";
}

//! A communicator of this process alone: the "one rank" solver of contract
//! item 8. Without MPI every communicator is that.
MPI_Comm lc_self_comm() {
#ifdef _NO_MPI
  return MPI_COMM_WORLD;
#else
  return MPI_COMM_SELF;
#endif
}

//! Rank and size of MPI_COMM_WORLD. A tensor is made first because mptensor
//! initializes MPI on first use and this file has no MPI_Init of its own.
std::pair<int, int> lc_world_rank_size() {
  const real_tensor probe(MPI_COMM_WORLD, mptensor::Shape(1));
  return {probe.get_comm_rank(), probe.get_comm_size()};
}

// |got - want| <= rtol * max(|want|, scale). scale is a magnitude that does
// not cancel (lc_ref::scale), so that a small result is not judged by a
// tolerance proportional to itself.
bool lc_check_close(const std::string& label, lc_cplx got, lc_cplx want,
                    double scale, double rtol) {
  const double tol = rtol * std::max(std::abs(want), scale);
  const double diff = std::abs(got - want);
  INFO(label << ": got=" << got << " want=" << want << " |diff|=" << diff
             << " tol=" << tol);
  CHECK(diff <= tol);
  return diff <= tol;
}

inline rf::parity_vector lc_phys(int d) {
  if (d == 2) {
    return {false, true};
  }
  if (d == 4) {
    return {false, true, true, false};
  }
  throw std::runtime_error("lc_phys: unsupported physical dimension");
}

// ---- one-site operators as element tables ---------------------------------

// m[in * d + out] = <out|A|in>.
struct lc_onesite {
  int d = 2;
  std::vector<lc_cplx> m;
  bool odd = false;
};

inline lc_onesite lc_onesite_op(const std::string& name, int d) {
  lc_onesite o;
  o.d = d;
  o.m.assign(d * d, 0.0);
  auto set = [&](int in, int out, double v) { o.m[in * d + out] = v; };
  if (name == "I") {
    for (int i = 0; i < d; ++i) {
      set(i, i, 1.0);
    }
    return o;
  }
  if (d == 2) {
    if (name == "c") {
      set(1, 0, 1.0);
      o.odd = true;
    } else if (name == "cdag") {
      set(0, 1, 1.0);
      o.odd = true;
    } else if (name == "n") {
      set(1, 1, 1.0);
    } else {
      throw std::runtime_error("lc_onesite_op: unknown d=2 operator " + name);
    }
    return o;
  }
  if (d != 4) {
    throw std::runtime_error("lc_onesite_op: unsupported dimension");
  }
  if (name == "c_up") {
    set(1, 0, 1.0);  // c_up |up> = |0>
    set(3, 2, 1.0);  // c_up |updn> = |dn>
    o.odd = true;
  } else if (name == "c_dn") {
    set(2, 0, 1.0);   // c_dn |dn> = |0>
    set(3, 1, -1.0);  // c_dn |updn> = -|up>
    o.odd = true;
  } else if (name == "cdag_up") {
    set(0, 1, 1.0);
    set(2, 3, 1.0);
    o.odd = true;
  } else if (name == "cdag_dn") {
    set(0, 2, 1.0);
    set(1, 3, -1.0);
    o.odd = true;
  } else if (name == "n") {
    set(1, 1, 1.0);
    set(2, 2, 1.0);
    set(3, 3, 2.0);
  } else if (name == "Sz") {
    set(1, 1, 0.5);
    set(2, 2, -0.5);
  } else {
    throw std::runtime_error("lc_onesite_op: unknown d=4 operator " + name);
  }
  return o;
}

inline void lc_set_scalar(real_tensor& t, const mptensor::Index& idx,
                          lc_cplx v) {
  REQUIRE(v.imag() == 0.0);
  t.set_value(idx, v.real());
}

inline void lc_set_scalar(complex_tensor& t, const mptensor::Index& idx,
                          lc_cplx v) {
  t.set_value(idx, v);
}

template <class tensor>
tensor lc_onesite_tensor(const lc_onesite& o, MPI_Comm comm) {
  tensor t(comm, mptensor::Shape(o.d, o.d));
  for (int in = 0; in < o.d; ++in) {
    for (int out = 0; out < o.d; ++out) {
      const lc_cplx v = o.m[in * o.d + out];
      if (v != 0.0) {
        lc_set_scalar(t, mptensor::Index(in, out), v);
      }
    }
  }
  return t;
}

//! The contract's product written out by the test (design section 4.2):
//! op4[i_s, i_t, o_s, o_t] = (-1)^{p_B p(i_s)} A[i_s, o_s] B[i_t, o_t].
//! Used by contract item 5 only, so that the nearest-neighbour anchor does
//! not go through product_twosite_op().
template <class tensor>
tensor lc_product_table(const lc_onesite& A, const lc_onesite& B,
                        MPI_Comm comm) {
  const int ds = A.d;
  const int dt = B.d;
  const rf::parity_vector ps = lc_phys(ds);
  tensor op(comm, mptensor::Shape(ds, dt, ds, dt));
  for (int is = 0; is < ds; ++is) {
    for (int it = 0; it < dt; ++it) {
      for (int os = 0; os < ds; ++os) {
        for (int ot = 0; ot < dt; ++ot) {
          const double sign = (B.odd && ps[is]) ? -1.0 : 1.0;
          const lc_cplx v = sign * A.m[is * ds + os] * B.m[it * dt + ot];
          if (v != 0.0) {
            lc_set_scalar(op, mptensor::Index(is, it, os, ot), v);
          }
        }
      }
    }
  }
  return op;
}

// ---- rank-consistent deterministic entries ---------------------------------

inline std::uint64_t lc_mix(std::uint64_t x) {  // splitmix64 finalizer
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

//! A number in [-1, 1) that depends on (seed, tag, linear global index,
//! part) only, never on the rank that owns the element.
inline double lc_unit(unsigned seed, int tag, std::size_t lin, int part) {
  std::uint64_t k = lc_mix(static_cast<std::uint64_t>(seed) * 1000003ULL +
                           static_cast<std::uint64_t>(tag));
  k = lc_mix(k ^ (static_cast<std::uint64_t>(lin) * 2ULL +
                  static_cast<std::uint64_t>(part)));
  return 2.0 * (static_cast<double>(k >> 11) * 0x1.0p-53) - 1.0;
}

//! Row-major position of a global index.
inline std::size_t lc_linear(const mptensor::Index& idx,
                             const mptensor::Shape& shape) {
  std::size_t lin = 0;
  for (std::size_t ax = 0; ax < shape.size(); ++ax) {
    lin = lin * shape[ax] + idx[ax];
  }
  return lin;
}

inline std::size_t lc_total(const mptensor::Shape& shape) {
  std::size_t n = 1;
  for (std::size_t ax = 0; ax < shape.size(); ++ax) {
    n *= shape[ax];
  }
  return n;
}

inline double lc_det_value(real_tensor*, unsigned seed, int tag,
                           std::size_t lin, double scale) {
  return scale * lc_unit(seed, tag, lin, 0);
}

inline lc_cplx lc_det_value(complex_tensor*, unsigned seed, int tag,
                            std::size_t lin, double scale) {
  return scale *
         lc_cplx(lc_unit(seed, tag, lin, 0), lc_unit(seed, tag, lin, 1));
}

template <class tensor>
double lc_global_max_abs(const tensor& t) {
  std::vector<double> m{0.0};
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    m[0] = std::max(m[0], std::abs(t[n]));
  }
  tenes::allreduce_max(m, t.get_comm());
  return m[0];
}

//! Parity-even, site-distinct Tn whose entries are a hash of (seed, site,
//! global index), every element multiplied by odd_scale per odd virtual
//! index (the T2 solver fixture, made rank-consistent).
template <class tensor>
void lc_seed_Tn(lc_state<tensor>& state, double odd_scale, unsigned seed) {
  auto& Tn = lc_acc::Tn(state);
  const auto& fi = lc_acc::finfo(state);
  for (std::size_t s = 0; s < Tn.size(); ++s) {
    const rf::leg_parities parity = rf::Tn_parity(fi, static_cast<int>(s));
    mptensor::Shape sh;
    for (const auto& leg : parity) {
      sh.push(leg.size());
    }
    tensor t(Tn[s].get_comm(), sh);
    for (std::size_t n = 0; n < t.local_size(); ++n) {
      const mptensor::Index idx = t.global_index(n);
      if (rf::count_odd(parity, idx) % 2 != 0) {
        continue;
      }
      double scale = 1.0;
      for (int leg = 0; leg < 4; ++leg) {
        if (parity[leg][idx[leg]]) {
          scale *= odd_scale;
        }
      }
      t.set_value(idx,
                  lc_det_value(static_cast<tensor*>(nullptr), seed,
                               static_cast<int>(s), lc_linear(idx, sh), scale));
    }
    REQUIRE(rf::parity_violation(rf::ftensor<tensor>{t, parity}) == 0.0);
    Tn[s] = t;
  }
}

template <class tensor>
std::vector<std::vector<tensor>*> lc_env_slots(lc_state<tensor>& state) {
  return {&lc_acc::C1(state),  &lc_acc::C2(state),  &lc_acc::C3(state),
          &lc_acc::C4(state),  &lc_acc::eTt(state), &lc_acc::eTr(state),
          &lc_acc::eTb(state), &lc_acc::eTl(state)};
}

//! Additive noise on every environment tensor, of relative size eps:
//! t += eps * max|t| * u. Additive, so that the elements the converged CTM
//! left exactly zero (a parity block, say) become non-zero too.
template <class tensor>
void lc_perturb_env(lc_state<tensor>& state, double eps, unsigned seed) {
  const auto slots = lc_env_slots(state);
  for (std::size_t k = 0; k < slots.size(); ++k) {
    for (std::size_t i = 0; i < slots[k]->size(); ++i) {
      tensor& t = (*slots[k])[i];
      const double scale = eps * lc_global_max_abs(t);
      const mptensor::Shape sh = t.shape();
      const int tag = static_cast<int>(1000 + 16 * k + i);
      for (std::size_t n = 0; n < t.local_size(); ++n) {
        t[n] += lc_det_value(static_cast<tensor*>(nullptr), seed, tag,
                             lc_linear(t.global_index(n), sh), scale);
      }
    }
  }
}

//! Copy every environment tensor of `from` into `to`, whatever the two
//! communicators: the full element list is assembled by a sum over the ranks
//! of `from` (every element has exactly one owner), then each rank of `to`
//! takes its own elements.
template <class tensor>
void lc_copy_env(lc_state<tensor>& from, lc_state<tensor>& to) {
  const auto src = lc_env_slots(from);
  const auto dst = lc_env_slots(to);
  const auto to_comm = lc_acc::Tn(to)[0].get_comm();
  for (std::size_t k = 0; k < src.size(); ++k) {
    REQUIRE(src[k]->size() == dst[k]->size());
    for (std::size_t i = 0; i < src[k]->size(); ++i) {
      const tensor& a = (*src[k])[i];
      const mptensor::Shape sh = a.shape();
      std::vector<typename tensor::value_type> full(lc_total(sh), 0.0);
      for (std::size_t n = 0; n < a.local_size(); ++n) {
        full[lc_linear(a.global_index(n), sh)] = a[n];
      }
      tenes::allreduce_sum(full, a.get_comm());
      tensor b(to_comm, sh);
      for (std::size_t n = 0; n < b.local_size(); ++n) {
        b[n] = full[lc_linear(b.global_index(n), sh)];
      }
      (*dst[k])[i] = b;
    }
  }
}

// ---- the solver fixture ----------------------------------------------------

constexpr int lc_D = 2;
constexpr int lc_chi = 4;
constexpr int lc_ctm_iteration_max = 30;
constexpr double lc_ctm_epsilon = 1.0e-10;
//! Odd virtual components are scaled by this per odd index (as in T2).
constexpr double lc_odd_scale = 0.6;
constexpr unsigned lc_seed = 3031;
constexpr double lc_env_noise = 0.3;

//! Contract items 4, 5 and 6.
constexpr double lc_rtol = 1.0e-10;
//! Contract item 8 (one rank against two).
constexpr double lc_rtol_mpi = 1.0e-12;
//! Premise: a compared value carries signal far above the tolerance.
constexpr double lc_min_signal = 1.0e-6;

tenes::SquareLattice lc_lattice(int lx, int ly, int d) {
  tenes::SquareLattice lattice(lx, ly, 0);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    lattice.physical_dims[site] = d;
    lattice.virtual_dims[site] = {lc_D, lc_D, lc_D, lc_D};
    lattice.initial_dirs[site] = {0.0};
    lattice.noises[site] = 1.0;
  }
  return lattice;
}

template <class tensor>
tenes::itps::PEPS_Parameters lc_params(int n_unit, int d, bool meanfield,
                                       const std::string& outdir) {
  tenes::itps::PEPS_Parameters p;
  p.fermion = true;
  p.is_real = std::is_same<tensor, real_tensor>::value;
  p.phys_parity.assign(n_unit, lc_phys(d));
  p.print_level = tenes::PrintLevel::none;
  p.outdir = outdir;
  p.CHI = lc_chi;
  p.Max_CTM_Iteration = lc_ctm_iteration_max;
  p.CTM_Convergence_Epsilon = lc_ctm_epsilon;
  p.Use_RSVD = false;
  p.MeanField_Env = meanfield;
  return p;
}

//! One correlation case: the one-site groups (every group on every site),
//! the [correlation] pairs (indices into `names`) and r_max.
struct lc_case {
  std::string label;
  int d = 2;
  std::vector<std::string> names;
  std::vector<std::pair<int, int>> pairs;
  int r_max = 3;
  double env_noise = 0.0;
  unsigned seed = lc_seed;
};

inline bool lc_same_parity(const lc_case& c, const std::pair<int, int>& p) {
  return lc_onesite_op(c.names[p.first], c.d).odd ==
         lc_onesite_op(c.names[p.second], c.d).odd;
}

inline CorrelationParameter lc_corparam(const lc_case& c) {
  std::vector<std::tuple<int, int>> ops;
  for (const auto& p : c.pairs) {
    ops.emplace_back(p.first, p.second);
  }
  return CorrelationParameter(c.r_max, ops);
}

//! A solver of a case on `comm`: the one-site groups on every site, the
//! given two-site observables and, with_correlation, the [correlation] pairs
//! (otherwise r_max = 0); Tn seeded, no CTM run yet.
template <class tensor>
std::unique_ptr<lc_state<tensor>> lc_make(
    const lc_case& c, MPI_Comm comm, const tenes::Operators<tensor>& twosite,
    bool with_correlation,
    const std::string& outdir = "output_test_fermion_longrange_correlation",
    tenes::itps::TransferMatrix_Parameters tmatrix =
        tenes::itps::TransferMatrix_Parameters{}) {
  const tenes::SquareLattice lattice = lc_lattice(2, 2, c.d);
  tenes::Operators<tensor> onesite;
  for (std::size_t g = 0; g < c.names.size(); ++g) {
    const tensor op =
        lc_onesite_tensor<tensor>(lc_onesite_op(c.names[g], c.d), comm);
    for (int s = 0; s < lattice.N_UNIT; ++s) {
      onesite.emplace_back(c.names[g], static_cast<int>(g), s, op);
    }
  }
  auto state = std::make_unique<lc_state<tensor>>(
      comm, lc_params<tensor>(lattice.N_UNIT, c.d, false, outdir), lattice,
      tenes::EvolutionOperators<tensor>{}, tenes::EvolutionOperators<tensor>{},
      onesite, twosite, tenes::Operators<tensor>{},
      with_correlation ? lc_corparam(c) : CorrelationParameter{}, tmatrix);
  lc_seed_Tn(*state, lc_odd_scale, c.seed);
  return state;
}

//! The window observables of a case: product_twosite_op(A, B) of every
//! same-parity pair at (r, 0) and (0, r), every site, r = 1 .. r_window.
//! window_group[p] is the group of pair p (-1 for a pair of different
//! parity).
template <class tensor>
tenes::Operators<tensor> lc_window_ops(const lc_case& c, MPI_Comm comm,
                                       int n_unit, int r_window,
                                       std::vector<int>& window_group) {
  const rf::parity_vector phys = lc_phys(c.d);
  tenes::Operators<tensor> twosite;
  window_group.assign(c.pairs.size(), -1);
  int group = 0;
  for (std::size_t p = 0; p < c.pairs.size(); ++p) {
    if (!lc_same_parity(c, c.pairs[p])) {
      continue;
    }
    const lc_onesite A = lc_onesite_op(c.names[c.pairs[p].first], c.d);
    const lc_onesite B = lc_onesite_op(c.names[c.pairs[p].second], c.d);
    const tensor op4 = rf::product_twosite_op(
        lc_onesite_tensor<tensor>(A, comm), lc_onesite_tensor<tensor>(B, comm),
        phys, phys, B.odd);
    window_group[p] = group;
    for (int s = 0; s < n_unit; ++s) {
      for (int r = 1; r <= r_window; ++r) {
        twosite.emplace_back("pair" + std::to_string(p), group, s, r, 0, op4);
        twosite.emplace_back("pair" + std::to_string(p), group, s, 0, r, op4);
      }
    }
    ++group;
  }
  return twosite;
}

//! Two solvers of one state. `ref` holds the window observables and no
//! [correlation] (r_max = 0), so that its measure_twosite() is the T2 path
//! alone and runs whatever the r_max guards do; `corr` holds the
//! [correlation] pairs and nothing else. Both carry the same Tn (a hash of
//! the global index) and, after lc_build, the same environment: the CTM of
//! `ref`, perturbed when env_noise > 0, copied element by element into
//! `corr` after `corr` ran its own CTM (so that nothing update_CTM() may
//! prepare is missing there).
template <class tensor>
struct lc_fixture {
  std::unique_ptr<lc_state<tensor>> ref;
  std::unique_ptr<lc_state<tensor>> corr;
  std::vector<int> window_group;
  int r_window = 0;
};

template <class tensor>
lc_fixture<tensor> lc_build(const lc_case& c, MPI_Comm ref_comm,
                            MPI_Comm corr_comm) {
  lc_fixture<tensor> fx;
  fx.r_window = std::min(c.r_max, 3);
  const int n_unit = lc_lattice(2, 2, c.d).N_UNIT;
  fx.ref = lc_make<tensor>(
      c, ref_comm,
      lc_window_ops<tensor>(c, ref_comm, n_unit, fx.r_window, fx.window_group),
      false);
  fx.ref->update_CTM();
  if (c.env_noise > 0.0) {
    lc_perturb_env(*fx.ref, c.env_noise, c.seed + 1);
  }
  fx.corr = lc_make<tensor>(c, corr_comm, tenes::Operators<tensor>{}, true);
  fx.corr->update_CTM();
  lc_copy_env(*fx.ref, *fx.corr);
  return fx;
}

// ---- the direct chain ------------------------------------------------------

//! A reference value. lc_chain_all() sets scale = sum_k |c_k| / |norm|;
//! lc_chain_refs() raises it to at least |A| |B| (largest elements of the
//! two one-site operators), the size <A_s B_t> has for a normalized state.
//! Neither cancels: a pair such as (n, n) has one channel, so the first is
//! |value| itself, and a value of 1e-4 then carries rounding errors of 1e-17
//! in absolute terms (1e-13 relative, measured at two ranks) that a
//! tolerance relative to the value alone would count against it.
template <class tensor>
struct lc_ref {
  lc_cplx value = 0.0;  //!< sum_k c_k / norm
  double scale = 0.0;   //!< a magnitude of the value that does not cancel
  std::size_t nchannel = 0;
};

inline double lc_max_abs(const lc_onesite& o) {
  double m = 0.0;
  for (const lc_cplx& x : o.m) {
    m = std::max(m, std::abs(x));
  }
  return m;
}

//! Sites of the straight chain from `left` over r steps to the right
//! (horizontal) or to the top (vertical).
inline std::vector<int> lc_chain_sites(const tenes::SquareLattice& lat,
                                       int left, int r, bool vertical) {
  std::vector<int> sites(r + 1);
  sites[0] = left;
  for (int i = 1; i <= r; ++i) {
    sites[i] = vertical ? lat.top(sites[i - 1]) : lat.right(sites[i - 1]);
  }
  REQUIRE(sites[r] == lat.other(left, vertical ? 0 : r, vertical ? r : 0));
  return sites;
}

//! A chain of folded rank-6 tensors (chain order, left / lower end first,
//! already rotated for a vertical chain) closed with the solver's
//! environment through the density correlation kernels. The environment
//! assignment is the one of the bosonic measure_correlation_ctm().
template <class tensor>
typename tensor::value_type lc_close_chain(lc_state<tensor>& st,
                                           const std::vector<int>& sites,
                                           bool vertical,
                                           const std::vector<tensor>& chain,
                                           const tensor& op_first,
                                           const tensor& op_last) {
  namespace core = tenes::itps::core;
  const int r = static_cast<int>(sites.size()) - 1;
  const int s0 = sites[0];
  const int sr = sites[r];
  tensor A;
  if (!vertical) {
    core::StartCorrelation_density_CTM(
        A, lc_acc::C1(st)[s0], lc_acc::C4(st)[s0], lc_acc::eTt(st)[s0],
        lc_acc::eTb(st)[s0], lc_acc::eTl(st)[s0], chain[0], op_first);
    for (int i = 1; i < r; ++i) {
      core::Transfer_density_CTM(A, lc_acc::eTt(st)[sites[i]],
                                 lc_acc::eTb(st)[sites[i]], chain[i]);
    }
    return core::FinishCorrelation_density_CTM(
        A, lc_acc::C2(st)[sr], lc_acc::C3(st)[sr], lc_acc::eTt(st)[sr],
        lc_acc::eTr(st)[sr], lc_acc::eTb(st)[sr], chain[r], op_last);
  }
  core::StartCorrelation_density_CTM(A, lc_acc::C4(st)[s0], lc_acc::C3(st)[s0],
                                     lc_acc::eTl(st)[s0], lc_acc::eTr(st)[s0],
                                     lc_acc::eTb(st)[s0], chain[0], op_first);
  for (int i = 1; i < r; ++i) {
    core::Transfer_density_CTM(A, lc_acc::eTl(st)[sites[i]],
                               lc_acc::eTr(st)[sites[i]], chain[i]);
  }
  return core::FinishCorrelation_density_CTM(
      A, lc_acc::C1(st)[sr], lc_acc::C2(st)[sr], lc_acc::eTl(st)[sr],
      lc_acc::eTt(st)[sr], lc_acc::eTr(st)[sr], chain[r], op_last);
}

template <class tensor>
tensor lc_orient(const tensor& t, bool vertical) {
  return vertical ? transpose(t, mptensor::Axes(3, 0, 1, 2, 4, 5)) : t;
}

//! The norm chain (plain reduced tensors, identity operators).
template <class tensor>
typename tensor::value_type lc_chain_norm(lc_state<tensor>& st,
                                          const std::vector<int>& sites,
                                          bool vertical) {
  const auto& fi = lc_acc::finfo(st);
  const auto& Tn = lc_acc::Tn(st);
  const auto comm = Tn[0].get_comm();
  std::vector<tensor> chain;
  for (const int s : sites) {
    chain.push_back(
        lc_orient(rf::build_reduced_op(rf::wrap_Tn(Tn[s], fi, s)), vertical));
  }
  const int d0 = static_cast<int>(fi.phys[sites.front()].size());
  const int dr = static_cast<int>(fi.phys[sites.back()].size());
  return lc_close_chain(
      st, sites, vertical, chain,
      lc_onesite_tensor<tensor>(lc_onesite_op("I", d0), comm),
      lc_onesite_tensor<tensor>(lc_onesite_op("I", dr), comm));
}

//! Design section 3.4, assembled by the test, for r = 1 .. r_max at once.
//! Per graded-SVD channel of op4, T1's relay site tensors with the roles and
//! legs build_relay_window() gives a straight path (horizontal: source exit
//! r, middle entry l and exit r, target entry l; vertical: source exit t,
//! middle entry b and exit t, target entry b), rotated for a vertical chain,
//! go through StartCorrelation, Transfer and FinishCorrelation_density_CTM;
//! the channel values are summed and divided by the same chain of plain
//! reduced tensors. The [truth] case T3-0 pins this to the window
//! measurement at r <= 3.
template <class tensor>
std::vector<lc_ref<tensor>> lc_chain_all(lc_state<tensor>& st, int left,
                                         int r_max, bool vertical,
                                         const tensor& op4) {
  namespace core = tenes::itps::core;
  using value_type = typename tensor::value_type;
  const tenes::SquareLattice& lat = lc_acc::lattice(st);
  const auto& fi = lc_acc::finfo(st);
  const auto& Tn = lc_acc::Tn(st);
  const auto comm = Tn[0].get_comm();
  const std::vector<int> sites = lc_chain_sites(lat, left, r_max, vertical);
  const int s0 = sites[0];
  for (const int s : sites) {
    // One channel set serves every r.
    REQUIRE(fi.phys[s] == fi.phys[s0]);
  }
  const int exit_leg = vertical ? 1 : 2;
  const int entry_leg = vertical ? 3 : 0;
  const tensor id = lc_onesite_tensor<tensor>(
      lc_onesite_op("I", static_cast<int>(fi.phys[s0].size())), comm);

  const auto start = [&](tensor& A, const tensor& t) {
    if (!vertical) {
      core::StartCorrelation_density_CTM(
          A, lc_acc::C1(st)[s0], lc_acc::C4(st)[s0], lc_acc::eTt(st)[s0],
          lc_acc::eTb(st)[s0], lc_acc::eTl(st)[s0], t, id);
    } else {
      core::StartCorrelation_density_CTM(
          A, lc_acc::C4(st)[s0], lc_acc::C3(st)[s0], lc_acc::eTl(st)[s0],
          lc_acc::eTr(st)[s0], lc_acc::eTb(st)[s0], t, id);
    }
  };
  const auto transfer = [&](tensor& A, int s, const tensor& t) {
    if (!vertical) {
      core::Transfer_density_CTM(A, lc_acc::eTt(st)[s], lc_acc::eTb(st)[s], t);
    } else {
      core::Transfer_density_CTM(A, lc_acc::eTl(st)[s], lc_acc::eTr(st)[s], t);
    }
  };
  const auto finish = [&](const tensor& A, int s,
                          const tensor& t) -> value_type {
    if (!vertical) {
      return core::FinishCorrelation_density_CTM(
          A, lc_acc::C2(st)[s], lc_acc::C3(st)[s], lc_acc::eTt(st)[s],
          lc_acc::eTr(st)[s], lc_acc::eTb(st)[s], t, id);
    }
    return core::FinishCorrelation_density_CTM(
        A, lc_acc::C1(st)[s], lc_acc::C2(st)[s], lc_acc::eTl(st)[s],
        lc_acc::eTt(st)[s], lc_acc::eTr(st)[s], t, id);
  };
  const auto wrap = [&](int s) { return rf::wrap_Tn(Tn[s], fi, s); };

  // The norm chain.
  std::map<int, tensor> reduced;
  const auto red = [&](int s) -> const tensor& {
    if (reduced.count(s) == 0) {
      reduced[s] = lc_orient(rf::build_reduced_op(wrap(s)), vertical);
    }
    return reduced.at(s);
  };
  std::vector<value_type> norm(r_max + 1, 0.0);
  {
    tensor A;
    start(A, red(s0));
    for (int r = 1; r <= r_max; ++r) {
      norm[r] = finish(A, sites[r], red(sites[r]));
      if (r < r_max) {
        transfer(A, sites[r], red(sites[r]));
      }
    }
  }

  const auto channels =
      rf::relay_channels(rf::wrap_twosite_gate(op4, fi.phys[s0], fi.phys[s0]));
  REQUIRE(!channels.empty());
  std::vector<value_type> sum(r_max + 1, 0.0);
  std::vector<double> abs_sum(r_max + 1, 0.0);
  for (const auto& ch : channels) {
    std::map<int, tensor> target;
    std::map<int, tensor> middle;
    tensor A;
    start(A, lc_orient(rf::build_relay_site(wrap(s0), rf::relay_role::source,
                                            -1, exit_leg, ch),
                       vertical));
    for (int r = 1; r <= r_max; ++r) {
      const int s = sites[r];
      if (target.count(s) == 0) {
        target[s] =
            lc_orient(rf::build_relay_site(wrap(s), rf::relay_role::target,
                                           entry_leg, -1, ch),
                      vertical);
      }
      const value_type ck = finish(A, s, target.at(s));
      sum[r] += ck;
      abs_sum[r] += std::abs(ck);
      if (r < r_max) {
        if (middle.count(s) == 0) {
          middle[s] =
              lc_orient(rf::build_relay_site(wrap(s), rf::relay_role::middle,
                                             entry_leg, exit_leg, ch),
                        vertical);
        }
        transfer(A, s, middle.at(s));
      }
    }
  }
  std::vector<lc_ref<tensor>> refs(r_max);
  for (int r = 1; r <= r_max; ++r) {
    REQUIRE(std::abs(norm[r]) > 0.0);
    refs[r - 1].value = lc_cplx(sum[r] / norm[r]);
    refs[r - 1].scale = abs_sum[r] / std::abs(norm[r]);
    refs[r - 1].nchannel = channels.size();
  }
  return refs;
}

template <class tensor>
lc_ref<tensor> lc_chain(lc_state<tensor>& st, int left, int r, bool vertical,
                        const tensor& op4) {
  return lc_chain_all(st, left, r, vertical, op4).back();
}

//! The plain contraction of <A_left B_right>: the reduced tensors with A and
//! B on the physical legs of the ends, no string. For a pair of different
//! parity this is what an implementation that forgets the parity rule would
//! report.
template <class tensor>
lc_cplx lc_plain_chain(lc_state<tensor>& st, int left, int r, bool vertical,
                       const lc_onesite& A, const lc_onesite& B) {
  const auto& fi = lc_acc::finfo(st);
  const auto& Tn = lc_acc::Tn(st);
  const auto comm = Tn[0].get_comm();
  const std::vector<int> sites =
      lc_chain_sites(lc_acc::lattice(st), left, r, vertical);
  std::vector<tensor> chain;
  for (const int s : sites) {
    chain.push_back(
        lc_orient(rf::build_reduced_op(rf::wrap_Tn(Tn[s], fi, s)), vertical));
  }
  const auto raw = lc_close_chain(st, sites, vertical, chain,
                                  lc_onesite_tensor<tensor>(A, comm),
                                  lc_onesite_tensor<tensor>(B, comm));
  return lc_cplx(raw / lc_chain_norm(st, sites, vertical));
}

// ---- correlation rows ------------------------------------------------------

//! (left_index, right_dx, right_dy, left_op, right_op)
using lc_key = std::tuple<int, int, int, int, int>;

inline std::string lc_key_name(const lc_key& k) {
  return "left " + std::to_string(std::get<0>(k)) + " (dx,dy)=(" +
         std::to_string(std::get<1>(k)) + "," + std::to_string(std::get<2>(k)) +
         ") ops [" + std::to_string(std::get<3>(k)) + ", " +
         std::to_string(std::get<4>(k)) + "]";
}

inline std::map<lc_key, std::vector<lc_cplx>> lc_rows(
    const std::vector<Correlation>& corr) {
  std::map<lc_key, std::vector<lc_cplx>> rows;
  for (const Correlation& c : corr) {
    rows[{c.left_index, c.right_dx, c.right_dy, c.left_op, c.right_op}]
        .push_back(lc_cplx(c.real, c.imag));
  }
  return rows;
}

//! The row of a key, required to be there exactly once.
inline lc_cplx lc_row(const std::map<lc_key, std::vector<lc_cplx>>& rows,
                      const lc_key& k) {
  INFO("correlation row " << lc_key_name(k));
  REQUIRE(rows.count(k) == 1);
  REQUIRE(rows.at(k).size() == 1);
  return rows.at(k)[0];
}

template <class tensor>
std::vector<Correlation> lc_measure(lc_state<tensor>& st) {
  std::vector<Correlation> corr;
  REQUIRE_NOTHROW(corr = st.measure_correlation());
  return corr;
}

//! Every (left, pair, r, direction) of the case, horizontal then vertical.
struct lc_item {
  int left;
  std::size_t pair;
  int r;
  bool vertical;
};

inline std::vector<lc_item> lc_items(const lc_case& c, int n_unit) {
  std::vector<lc_item> items;
  for (int s = 0; s < n_unit; ++s) {
    for (std::size_t p = 0; p < c.pairs.size(); ++p) {
      for (int r = 1; r <= c.r_max; ++r) {
        for (const bool v : {false, true}) {
          items.push_back({s, p, r, v});
        }
      }
    }
  }
  return items;
}

inline lc_key lc_item_key(const lc_case& c, const lc_item& it) {
  return {it.left, it.vertical ? 0 : it.r, it.vertical ? it.r : 0,
          c.pairs[it.pair].first, c.pairs[it.pair].second};
}

template <class tensor>
tensor lc_pair_op4(const lc_case& c, std::size_t p, MPI_Comm comm) {
  const lc_onesite A = lc_onesite_op(c.names[c.pairs[p].first], c.d);
  const lc_onesite B = lc_onesite_op(c.names[c.pairs[p].second], c.d);
  const rf::parity_vector phys = lc_phys(c.d);
  return rf::product_twosite_op(lc_onesite_tensor<tensor>(A, comm),
                                lc_onesite_tensor<tensor>(B, comm), phys, phys,
                                B.odd);
}

//! Direct-chain references of every same-parity item of a case.
template <class tensor>
std::map<lc_key, lc_ref<tensor>> lc_chain_refs(const lc_case& c,
                                               lc_state<tensor>& st) {
  const auto comm = lc_acc::Tn(st)[0].get_comm();
  std::map<lc_key, lc_ref<tensor>> refs;
  std::vector<tensor> op4(c.pairs.size());
  for (std::size_t p = 0; p < c.pairs.size(); ++p) {
    if (lc_same_parity(c, c.pairs[p])) {
      op4[p] = lc_pair_op4<tensor>(c, p, comm);
    }
  }
  for (int left = 0; left < lc_acc::lattice(st).N_UNIT; ++left) {
    for (std::size_t p = 0; p < c.pairs.size(); ++p) {
      if (!lc_same_parity(c, c.pairs[p])) {
        continue;
      }
      for (const bool v : {false, true}) {
        const double op_scale =
            lc_max_abs(lc_onesite_op(c.names[c.pairs[p].first], c.d)) *
            lc_max_abs(lc_onesite_op(c.names[c.pairs[p].second], c.d));
        const auto all = lc_chain_all(st, left, c.r_max, v, op4[p]);
        for (int r = 1; r <= c.r_max; ++r) {
          lc_ref<tensor> ref = all[r - 1];
          ref.scale = std::max(ref.scale, op_scale);
          refs[lc_item_key(c, lc_item{left, p, r, v})] = ref;
        }
      }
    }
  }
  return refs;
}

//! Contract items 2 and 3 on one case: one row per (left, pair, r,
//! direction) and no other; a pair of different parity is exactly 0.0.
void lc_check_shape_and_zeros(const lc_case& c, int n_unit,
                              const std::vector<Correlation>& corr) {
  const auto rows = lc_rows(corr);
  const auto items = lc_items(c, n_unit);
  CHECK(corr.size() == items.size());
  for (const lc_item& it : items) {
    const lc_key k = lc_item_key(c, it);
    INFO(c.label << " " << lc_key_name(k));
    CHECK(rows.count(k) == 1);
    if (rows.count(k) == 1) {
      CHECK(rows.at(k).size() == 1);
      if (!lc_same_parity(c, c.pairs[it.pair])) {
        // Exactly zero, real and imaginary parts (contract item 3).
        CHECK(rows.at(k)[0].real() == 0.0);
        CHECK(rows.at(k)[0].imag() == 0.0);
      }
    }
  }
}

//! Contract item 4 (and 5, 6): every same-parity row against the window
//! measurement (r <= r_window) and the direct chain (every r).
template <class tensor>
void lc_run_window_agreement(const lc_case& c) {
  INFO(c.label);
  auto fx = lc_build<tensor>(c, MPI_COMM_WORLD, MPI_COMM_WORLD);
  const int n_unit = lc_acc::lattice(*fx.ref).N_UNIT;

  // References first, so that a throw of the code under test comes after
  // them.
  const auto refs = lc_chain_refs(c, *fx.ref);
  double min_signal = std::numeric_limits<double>::infinity();
  for (const auto& [k, ref] : refs) {
    min_signal = std::min(min_signal, std::abs(ref.value));
  }
  {
    INFO("premise: every compared value carries signal; smallest |value| "
         << min_signal);
    REQUIRE(min_signal > lc_min_signal);
  }
  std::vector<std::map<Bond, typename tensor::value_type>> window;
  REQUIRE_NOTHROW(window = fx.ref->measure_twosite());

  const auto corr = lc_measure(*fx.corr);
  lc_check_shape_and_zeros(c, n_unit, corr);
  const auto rows = lc_rows(corr);
  double max_rel_window = 0.0;
  double max_rel_chain = 0.0;
  for (const lc_item& it : lc_items(c, n_unit)) {
    if (!lc_same_parity(c, c.pairs[it.pair])) {
      continue;
    }
    const lc_key k = lc_item_key(c, it);
    const std::string what = c.label + " " + lc_key_name(k) + " (" +
                             c.names[c.pairs[it.pair].first] + ", " +
                             c.names[c.pairs[it.pair].second] + ")";
    INFO(what);
    const lc_cplx got = lc_row(rows, k);
    if (std::is_same<tensor, real_tensor>::value) {
      CHECK(got.imag() == 0.0);
    }
    const lc_ref<tensor>& ref = refs.at(k);
    if (it.r <= fx.r_window) {
      const int g = fx.window_group[it.pair];
      REQUIRE(g >= 0);
      REQUIRE(static_cast<int>(window.size()) > g);
      const Bond b{it.left, std::get<1>(k), std::get<2>(k)};
      REQUIRE(window[g].count(b) == 1);
      const lc_cplx want = window[g].at(b);
      lc_check_close(what + " [correlation vs window measure_twosite]", got,
                     want, ref.scale, lc_rtol);
      max_rel_window =
          std::max(max_rel_window,
                   std::abs(got - want) / std::max(std::abs(want), ref.scale));
    }
    lc_check_close(what + " [correlation vs direct chain]", got, ref.value,
                   ref.scale, lc_rtol);
    max_rel_chain =
        std::max(max_rel_chain, std::abs(got - ref.value) /
                                    std::max(std::abs(ref.value), ref.scale));
  }
  std::cout << std::setprecision(3) << "longrange-corr " << c.label
            << ": max relative |correlation - window| " << max_rel_window
            << ", |correlation - chain| " << max_rel_chain << std::endl;
}

const std::vector<std::string> lc_names2 = {"n", "cdag", "c"};
const std::vector<std::string> lc_names4 = {"n",       "cdag_up", "c_up",
                                            "cdag_dn", "c_dn",    "Sz"};

//! d = 2 pairs of contract item 4: (c+, c), (c, c+), (n, n).
const std::vector<std::pair<int, int>> lc_pairs2 = {{1, 2}, {2, 1}, {0, 0}};
//! d = 4 pairs of contract item 4: (c+_up, c_up), (c+_dn, c_dn), (Sz, Sz).
const std::vector<std::pair<int, int>> lc_pairs4 = {{1, 2}, {3, 4}, {5, 5}};
//! More d = 4 pairs, on the perturbed environment: (c_up, c+_up), the
//! asymmetric even pair (n, Sz) and (c+_dn, c_dn) once more.
const std::vector<std::pair<int, int>> lc_pairs4_more = {
    {2, 1}, {0, 5}, {3, 4}};

}  // namespace

// ============================================================================
// The reference side
// ============================================================================

// The direct chain of design section 3.4 (T1's relay site tensors on the
// straight chain, the density correlation kernels) against the window
// measurement of T2 at r = 1, 2, 3, horizontal and vertical, on converged and
// perturbed environments. Neither side is T3 code: this pins the environment
// assignment, orientation and norm of the reference used by items 4 to 6
// and 8, and makes it the reference at r = 4, 5 (item 6), where no window
// exists.
TEST_CASE(
    "longrange-corr [truth] T3-0: the direct chain equals the verified window "
    "measurement at r = 1, 2, 3") {
  const auto run = [](auto tag, const lc_case& c) {
    using tensor = decltype(tag);
    INFO(c.label);
    auto fx = lc_build<tensor>(c, MPI_COMM_WORLD, MPI_COMM_WORLD);
    auto& st = *fx.ref;
    const auto refs = lc_chain_refs(c, st);
    const auto window = st.measure_twosite();
    double max_rel = 0.0;
    for (const lc_item& it : lc_items(c, lc_acc::lattice(st).N_UNIT)) {
      const lc_key k = lc_item_key(c, it);
      const int g = fx.window_group[it.pair];
      REQUIRE(g >= 0);
      const Bond b{it.left, std::get<1>(k), std::get<2>(k)};
      REQUIRE(window[g].count(b) == 1);
      const lc_ref<tensor>& ref = refs.at(k);
      const lc_cplx want = window[g].at(b);
      lc_check_close(c.label + " " + lc_key_name(k) + " [chain vs window]",
                     ref.value, want, ref.scale, 1.0e-12);
      max_rel = std::max(max_rel, std::abs(ref.value - want) /
                                      std::max(std::abs(want), ref.scale));
    }
    std::cout << std::setprecision(3) << "longrange-corr [truth] " << c.label
              << ": max relative |chain - window| " << max_rel << std::endl;
  };
  run(real_tensor{}, lc_case{"T3-0 d=2 real perturbed", 2, lc_names2, lc_pairs2,
                             3, lc_env_noise, lc_seed + 100});
  run(real_tensor{}, lc_case{"T3-0 d=4 real converged",
                             4,
                             lc_names4,
                             {{1, 2}, {0, 5}},
                             3,
                             0.0,
                             lc_seed + 101});
  run(complex_tensor{}, lc_case{"T3-0 d=2 complex perturbed",
                                2,
                                lc_names2,
                                {{1, 2}, {2, 1}},
                                3,
                                lc_env_noise,
                                lc_seed + 102});
}

// ============================================================================
// Contract item 1: input acceptance and rejection
// ============================================================================

namespace {

using lc_gtensor = real_tensor;

struct lc_guard_input {
  tenes::itps::PEPS_Parameters params;
  tenes::SquareLattice lattice = lc_lattice(2, 2, 2);
  tenes::Operators<lc_gtensor> onesite;
  CorrelationParameter corparam;

  //! One-site groups `names` on every site; the [correlation] pairs and
  //! r_max as given.
  lc_guard_input(int d, bool meanfield, const std::vector<std::string>& names,
                 int r_max, const std::vector<std::tuple<int, int>>& pairs)
      : corparam(r_max, pairs) {
    lattice = lc_lattice(2, 2, d);
    params = lc_params<lc_gtensor>(lattice.N_UNIT, d, meanfield,
                                   "output_test_fermion_longrange_correlation");
    for (std::size_t g = 0; g < names.size(); ++g) {
      for (int s = 0; s < lattice.N_UNIT; ++s) {
        onesite.emplace_back(names[g], static_cast<int>(g), s,
                             lc_onesite_tensor<lc_gtensor>(
                                 lc_onesite_op(names[g], d), MPI_COMM_WORLD));
      }
    }
  }

  void validate() const {
    tenes::itps::validate_fermion_constraints(
        params, lattice, tenes::EvolutionOperators<lc_gtensor>{},
        tenes::EvolutionOperators<lc_gtensor>{}, onesite,
        tenes::Operators<lc_gtensor>{}, tenes::Operators<lc_gtensor>{},
        corparam);
  }
};

void lc_check_accepts(const lc_guard_input& in, const std::string& what) {
  INFO(what);
  try {
    in.validate();
  } catch (const std::exception& e) {
    FAIL_CHECK(what << " was rejected: " << std::string(e.what()));
  }
}

void lc_check_rejects(const lc_guard_input& in, const std::string& what) {
  INFO(what);
  try {
    in.validate();
    FAIL_CHECK(what << " was accepted");
  } catch (const tenes::input_error&) {
  } catch (const std::exception& e) {
    FAIL_CHECK(what << " threw something other than tenes::input_error: "
                    << std::string(e.what()));
  }
}

}  // namespace

TEST_CASE(
    "longrange-corr T3-1a: r_max > 0 is accepted with the CTM environment, "
    "for odd, even and mixed-parity pairs") {
  lc_check_accepts(
      lc_guard_input(2, false, lc_names2, 3,
                     {{0, 0}, {1, 2}, {2, 1}, {1, 1}, {0, 2}, {1, 0}}),
      "d = 2, r_max = 3, pairs [0,0] [1,2] [2,1] [1,1] [0,2] [1,0]");
  lc_check_accepts(lc_guard_input(2, false, lc_names2, 5, {{1, 2}}),
                   "d = 2, r_max = 5, pair [1,2]");
  lc_check_accepts(lc_guard_input(2, false, {"n"}, 1, {{0, 0}}),
                   "d = 2, r_max = 1, even operators only");
  lc_check_accepts(lc_guard_input(4, false, lc_names4, 3,
                                  {{1, 2}, {3, 4}, {5, 5}, {0, 5}, {5, 2}}),
                   "d = 4, r_max = 3, pairs [1,2] [3,4] [5,5] [0,5] [5,2]");
}

// Task T5 of the same plan (contract item 1) lifts the mean-field
// restriction: this case used to pin that the mean-field environment rejects
// r_max > 0 at load and at measurement until T5 ("[kept] T3-1b ... stays
// rejected") and now requires that both accept it and that the measurement
// run. The values are checked in test/fermion/longrange_mf.cpp. Red until T5
// is implemented.
TEST_CASE(
    "longrange-corr T3-1b: with meanfield_env, r_max > 0 is accepted at load "
    "and measured (task T5)") {
  lc_check_accepts(lc_guard_input(2, true, lc_names2, 3, {{1, 2}}),
                   "mean field, r_max = 3, pair [1,2]");
  lc_check_accepts(lc_guard_input(2, true, {"n"}, 2, {{0, 0}}),
                   "mean field, r_max = 2, pair [0,0] (even only)");

  const lc_case c{"T3-1b mean field", 2, lc_names2, {{0, 0}, {1, 2}}, 2, 0.0,
                  lc_seed + 10};
  const tenes::SquareLattice lattice = lc_lattice(2, 2, 2);
  tenes::Operators<real_tensor> onesite;
  for (std::size_t g = 0; g < c.names.size(); ++g) {
    for (int s = 0; s < lattice.N_UNIT; ++s) {
      onesite.emplace_back(c.names[g], static_cast<int>(g), s,
                           lc_onesite_tensor<real_tensor>(
                               lc_onesite_op(c.names[g], 2), MPI_COMM_WORLD));
    }
  }
  lc_state<real_tensor> state(
      MPI_COMM_WORLD,
      lc_params<real_tensor>(lattice.N_UNIT, 2, true,
                             "output_test_fermion_longrange_correlation"),
      lattice, tenes::EvolutionOperators<real_tensor>{},
      tenes::EvolutionOperators<real_tensor>{}, onesite,
      tenes::Operators<real_tensor>{}, tenes::Operators<real_tensor>{},
      lc_corparam(c), tenes::itps::TransferMatrix_Parameters{});
  lc_seed_Tn(state, lc_odd_scale, c.seed);
  std::vector<Correlation> rows;
  try {
    rows = state.measure_correlation();
  } catch (const std::exception& e) {
    FAIL_CHECK("measure_correlation with meanfield_env and r_max > 0 threw: "
               << std::string(e.what()));
  }
  // 4 left sites x 2 pairs x r = 1, 2 x 2 directions.
  CHECK(rows.size() == 32);
  for (const Correlation& row : rows) {
    CHECK(std::isfinite(row.real));
    CHECK(std::isfinite(row.imag));
  }
}

// Green against the stub by construction (the stub rejects every r_max > 0);
// pins that accepting r_max > 0 does not extend to pairs that name no
// one-site operator, or a one-site operator of mixed parity (contract item
// 1: "each pair names defined one-site operators of definite parity").
TEST_CASE(
    "longrange-corr [kept] T3-1c: pairs naming an undefined or a mixed-parity "
    "one-site operator are rejected") {
  lc_check_rejects(lc_guard_input(2, false, lc_names2, 3, {{0, 3}}),
                   "pair [0,3] with one-site groups 0..2");
  lc_check_rejects(lc_guard_input(2, false, lc_names2, 3, {{7, 1}}),
                   "pair [7,1] with one-site groups 0..2");
  lc_check_rejects(lc_guard_input(2, false, lc_names2, 3, {{-1, 0}}),
                   "pair [-1,0]");

  lc_guard_input mixed(2, false, {"n"}, 3, {{0, 1}});
  lc_onesite o = lc_onesite_op("c", 2);
  const lc_onesite n = lc_onesite_op("n", 2);
  for (int k = 0; k < 4; ++k) {
    o.m[k] += 0.5 * n.m[k];
  }
  for (int s = 0; s < mixed.lattice.N_UNIT; ++s) {
    mixed.onesite.emplace_back(
        "mixed", 1, s, lc_onesite_tensor<lc_gtensor>(o, MPI_COMM_WORLD));
  }
  lc_check_rejects(mixed, "pair [0,1] where group 1 has mixed parity");
}

// ============================================================================
// Contract items 2 and 3: output shape, pairs of different parity
// ============================================================================

TEST_CASE(
    "longrange-corr T3-2: one row per left site, pair, r = 1..r_max and "
    "direction, as for bosons") {
  const lc_case c{"T3-2 d=2 real",
                  2,
                  lc_names2,
                  {{0, 0}, {1, 2}, {2, 1}, {1, 1}, {0, 2}, {1, 0}},
                  3,
                  lc_env_noise,
                  lc_seed + 20};
  auto fx = lc_build<real_tensor>(c, MPI_COMM_WORLD, MPI_COMM_WORLD);
  const int n_unit = lc_acc::lattice(*fx.corr).N_UNIT;
  const auto corr = lc_measure(*fx.corr);
  lc_check_shape_and_zeros(c, n_unit, corr);
  for (const Correlation& row : corr) {
    INFO("row left " << row.left_index << " (" << row.right_dx << ", "
                     << row.right_dy << ") ops [" << row.left_op << ", "
                     << row.right_op << "]");
    CHECK(((row.right_dx >= 1 && row.right_dy == 0) ||
           (row.right_dx == 0 && row.right_dy >= 1)));
    CHECK(std::isfinite(row.real));
    CHECK(row.imag == 0.0);  // real tensors
  }
}

namespace {

//! Contract item 3 on one case: the plain contraction of every
//! different-parity item is far from 0 (premise), and the solver reports
//! exactly 0.0 for it, in a row of its own.
template <class tensor>
void lc_run_mixed_zero(const lc_case& c) {
  INFO(c.label);
  auto fx = lc_build<tensor>(c, MPI_COMM_WORLD, MPI_COMM_WORLD);
  auto& st = *fx.ref;
  const int n_unit = lc_acc::lattice(st).N_UNIT;
  double min_plain = std::numeric_limits<double>::infinity();
  std::size_t n_mixed = 0;
  for (const lc_item& it : lc_items(c, n_unit)) {
    if (lc_same_parity(c, c.pairs[it.pair])) {
      continue;
    }
    ++n_mixed;
    const lc_cplx plain =
        lc_plain_chain(st, it.left, it.r, it.vertical,
                       lc_onesite_op(c.names[c.pairs[it.pair].first], c.d),
                       lc_onesite_op(c.names[c.pairs[it.pair].second], c.d));
    min_plain = std::min(min_plain, std::abs(plain));
  }
  REQUIRE(n_mixed > 0);
  {
    INFO(
        "premise: the plain contraction of a different-parity pair is not "
        "zero on this state; smallest |value| "
        << min_plain);
    REQUIRE(min_plain > lc_min_signal);
  }
  const auto corr = lc_measure(*fx.corr);
  lc_check_shape_and_zeros(c, n_unit, corr);
}

}  // namespace

TEST_CASE(
    "longrange-corr T3-3: a pair of different parity is exactly 0.0, real "
    "and imaginary parts, and its rows are written") {
  lc_run_mixed_zero<real_tensor>(lc_case{"T3-3 d=2 real perturbed",
                                         2,
                                         lc_names2,
                                         {{0, 2}, {1, 0}, {1, 2}},
                                         3,
                                         lc_env_noise,
                                         lc_seed + 30});
  lc_run_mixed_zero<real_tensor>(lc_case{"T3-3 d=4 real perturbed",
                                         4,
                                         lc_names4,
                                         {{5, 2}, {3, 0}, {1, 5}},
                                         2,
                                         lc_env_noise,
                                         lc_seed + 31});
  lc_run_mixed_zero<complex_tensor>(lc_case{"T3-3 d=2 complex perturbed",
                                            2,
                                            lc_names2,
                                            {{0, 1}, {2, 0}},
                                            2,
                                            lc_env_noise,
                                            lc_seed + 32});
}

// ============================================================================
// Contract item 4: the correlation function equals the window measurement
// ============================================================================

TEST_CASE(
    "longrange-corr T3-4a: r = 1, 2, 3 equal measure_twosite of "
    "product_twosite_op at (r, 0) and (0, r), d = 2") {
  lc_run_window_agreement<real_tensor>(lc_case{"T3-4a d=2 real converged", 2,
                                               lc_names2, lc_pairs2, 3, 0.0,
                                               lc_seed + 40});
  lc_run_window_agreement<real_tensor>(lc_case{"T3-4a d=2 real perturbed", 2,
                                               lc_names2, lc_pairs2, 3,
                                               lc_env_noise, lc_seed + 41});
  // complex_tensor (Review Focus 3 of T2; here <c+_s c_t> and <c_s c+_t>
  // are complex and not conjugate-symmetric in s <-> t).
  lc_run_window_agreement<complex_tensor>(lc_case{"T3-4a d=2 complex perturbed",
                                                  2, lc_names2, lc_pairs2, 3,
                                                  lc_env_noise, lc_seed + 42});
}

TEST_CASE(
    "longrange-corr T3-4b: r = 1, 2, 3 equal measure_twosite of "
    "product_twosite_op at (r, 0) and (0, r), d = 4") {
  lc_run_window_agreement<real_tensor>(lc_case{"T3-4b d=4 real converged", 4,
                                               lc_names4, lc_pairs4, 3, 0.0,
                                               lc_seed + 43});
  lc_run_window_agreement<real_tensor>(lc_case{"T3-4b d=4 real perturbed", 4,
                                               lc_names4, lc_pairs4_more, 3,
                                               lc_env_noise, lc_seed + 44});
}

// ============================================================================
// Contract item 5: r = 1 equals the existing nearest-neighbour path
// ============================================================================

TEST_CASE(
    "longrange-corr T3-5: r = 1 equals the existing nearest-neighbour "
    "measure_twosite of the product") {
  const auto run = [](const lc_case& c) {
    INFO(c.label);
    using tensor = real_tensor;
    // lc_build's pair of solvers, but the reference one measures the
    // nearest-neighbour observables built from the test's own product table.
    const int n_unit = lc_lattice(2, 2, c.d).N_UNIT;
    std::vector<tensor> op4;
    tenes::Operators<tensor> twosite;
    for (std::size_t p = 0; p < c.pairs.size(); ++p) {
      op4.push_back(lc_product_table<tensor>(
          lc_onesite_op(c.names[c.pairs[p].first], c.d),
          lc_onesite_op(c.names[c.pairs[p].second], c.d), MPI_COMM_WORLD));
      for (int s = 0; s < n_unit; ++s) {
        twosite.emplace_back("nn", static_cast<int>(p), s, 1, 0, op4[p]);
        twosite.emplace_back("nn", static_cast<int>(p), s, 0, 1, op4[p]);
      }
    }
    auto ref = lc_make<tensor>(c, MPI_COMM_WORLD, twosite, false);
    ref->update_CTM();
    lc_perturb_env(*ref, c.env_noise, c.seed + 1);
    auto corr =
        lc_make<tensor>(c, MPI_COMM_WORLD, tenes::Operators<tensor>{}, true);
    corr->update_CTM();
    lc_copy_env(*ref, *corr);

    const auto nn = ref->measure_twosite();
    std::map<lc_key, lc_cplx> want;
    std::map<lc_key, double> scale;
    for (int s = 0; s < n_unit; ++s) {
      for (std::size_t p = 0; p < c.pairs.size(); ++p) {
        for (const bool v : {false, true}) {
          const lc_key k{s, v ? 0 : 1, v ? 1 : 0, c.pairs[p].first,
                         c.pairs[p].second};
          const Bond b{s, std::get<1>(k), std::get<2>(k)};
          INFO(lc_key_name(k));
          REQUIRE(nn[p].count(b) == 1);
          want[k] = nn[p].at(b);
          scale[k] = std::max(
              lc_chain(*ref, s, 1, v, op4[p]).scale,
              lc_max_abs(lc_onesite_op(c.names[c.pairs[p].first], c.d)) *
                  lc_max_abs(lc_onesite_op(c.names[c.pairs[p].second], c.d)));
          INFO("premise: the value carries signal: " << want[k]);
          REQUIRE(std::abs(want[k]) > lc_min_signal);
        }
      }
    }
    const auto rows = lc_rows(lc_measure(*corr));
    for (const auto& kv : want) {
      const std::string what = c.label + " " + lc_key_name(kv.first);
      INFO(what);
      lc_check_close(what + " [r = 1 vs nearest-neighbour path]",
                     lc_row(rows, kv.first), kv.second, scale.at(kv.first),
                     lc_rtol);
    }
  };
  run(lc_case{"T3-5 d=2", 2, lc_names2, lc_pairs2, 1, lc_env_noise,
              lc_seed + 50});
  run(lc_case{"T3-5 d=4",
              4,
              lc_names4,
              {{1, 2}, {4, 3}, {5, 5}},
              1,
              lc_env_noise,
              lc_seed + 51});
}

// ============================================================================
// Contract item 6: r_max beyond two unit cells (Review Focus 1)
// ============================================================================

TEST_CASE(
    "longrange-corr T3-6: 2x2 cell with r_max = 5, r = 2, 3 equal the window "
    "measurement and r = 4, 5 the direct chain (Review Focus 1)") {
  {
    const tenes::SquareLattice lattice = lc_lattice(2, 2, 2);
    for (int s = 0; s < lattice.N_UNIT; ++s) {
      INFO("premise: r = 4 returns to the left site's sublattice, site " << s);
      REQUIRE(lattice.other(s, 4, 0) == s);
      REQUIRE(lattice.other(s, 0, 4) == s);
    }
  }
  lc_run_window_agreement<real_tensor>(
      lc_case{"T3-6 d=2 real perturbed r_max=5", 2, lc_names2, lc_pairs2, 5,
              lc_env_noise, lc_seed + 60});
  lc_run_window_agreement<real_tensor>(
      lc_case{"T3-6 d=2 real converged r_max=5",
              2,
              lc_names2,
              {{1, 2}},
              5,
              0.0,
              lc_seed + 61});
  lc_run_window_agreement<real_tensor>(
      lc_case{"T3-6 d=4 real perturbed r_max=5",
              4,
              lc_names4,
              {{1, 2}, {5, 5}},
              5,
              lc_env_noise,
              lc_seed + 62});
}

// ============================================================================
// Contract item 7: the correlation length stays disabled, with a warning
// ============================================================================

TEST_CASE(
    "longrange-corr T3-7: measure() with r_max > 0 writes the correlation "
    "function, disables the correlation length and warns") {
  const std::string outdir = "output_test_fermion_longrange_correlation_t3_7";
  const int rank = lc_world_rank_size().first;
  if (rank == 0) {
    std::error_code ec;
    std::filesystem::remove_all(outdir, ec);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  const lc_case c{"T3-7", 2, lc_names2, {{0, 0}, {1, 2}}, 2, 0.0, lc_seed + 70};
  tenes::itps::TransferMatrix_Parameters tmatrix;
  tmatrix.to_calculate = true;
  auto state =
      lc_make<real_tensor>(c, MPI_COMM_WORLD, tenes::Operators<real_tensor>{},
                           true, outdir, tmatrix);
  std::ostringstream err;
  std::streambuf* old = std::cerr.rdbuf(err.rdbuf());
  std::string failure;
  try {
    state->measure();
  } catch (const std::exception& e) {
    failure = e.what();
  }
  std::cerr.rdbuf(old);
  INFO("measure() threw: " << failure);
  REQUIRE(failure.empty());
  INFO("captured stderr: " << err.str());
  CHECK_FALSE(lc_acc::tmatrix_param(*state).to_calculate);
  if (rank == 0) {
    // One line that warns about the correlation length.
    std::istringstream lines(err.str());
    std::string line;
    bool warned = false;
    while (std::getline(lines, line)) {
      warned = warned || (line.find("WARNING") != std::string::npos &&
                          line.find("correlation_length") != std::string::npos);
    }
    CHECK(warned);
    CHECK_FALSE(std::filesystem::exists(outdir + "/correlation_length.dat"));
    REQUIRE(std::filesystem::exists(outdir + "/correlation.dat"));
    std::ifstream ifs(outdir + "/correlation.dat");
    int data_lines = 0;
    while (std::getline(ifs, line)) {
      if (!line.empty() && line[0] != '#') {
        ++data_lines;
      }
    }
    // 4 left sites x 2 pairs x r = 1, 2 x 2 directions.
    CHECK(data_lines == 4 * 2 * 2 * 2);
  }
}

// ============================================================================
// Contract item 8: one rank against two
// ============================================================================
//
// Registered twice: with the rest of test_fermion_longrange (ctest at one
// rank) and alone at two ranks (test_fermion_longrange_mpi2). The same state
// is built in one run on MPI_COMM_WORLD (distributed at two ranks) and on a
// communicator of this process alone (the one-rank solvers), all with the
// environment of the one-rank reference solver copied element by element.
// The two-rank correlation function must equal the one-rank one (1e-12), the
// one-rank window measurement and the one-rank direct chain (1e-10), and be
// the same on every rank; a pair of different parity must be exactly 0.0 on
// every rank.
// At one rank the two solvers coincide in layout, so the case then repeats
// item 4; the rank-count comparison bites at two.

namespace {

template <class tensor>
void lc_run_mpi(const lc_case& c) {
  INFO(c.label);
  const std::pair<int, int> rank_size = lc_world_rank_size();
  INFO("MPI_COMM_WORLD rank " << rank_size.first << " of " << rank_size.second);
  // One-rank pair of solvers (reference and correlation) and the
  // distributed correlation solver, all with the environment of the
  // one-rank reference solver.
  auto self = lc_build<tensor>(c, lc_self_comm(), lc_self_comm());
  auto world =
      lc_make<tensor>(c, MPI_COMM_WORLD, tenes::Operators<tensor>{}, true);
  world->update_CTM();
  lc_copy_env(*self.ref, *world);
  const int n_unit = lc_acc::lattice(*world).N_UNIT;

  // One-rank references.
  const auto refs = lc_chain_refs(c, *self.ref);
  std::vector<std::map<Bond, typename tensor::value_type>> window;
  REQUIRE_NOTHROW(window = self.ref->measure_twosite());
  double min_signal = std::numeric_limits<double>::infinity();
  for (const auto& [k, ref] : refs) {
    min_signal = std::min(min_signal, std::abs(ref.value));
  }
  {
    INFO("premise: smallest |value| " << min_signal);
    REQUIRE(min_signal > lc_min_signal);
  }

  const auto corr_world = lc_measure(*world);
  const auto corr_self = lc_measure(*self.corr);
  lc_check_shape_and_zeros(c, n_unit, corr_world);
  const auto rows_world = lc_rows(corr_world);
  const auto rows_self = lc_rows(corr_self);
  for (const lc_item& it : lc_items(c, n_unit)) {
    if (!lc_same_parity(c, c.pairs[it.pair])) {
      continue;
    }
    const lc_key k = lc_item_key(c, it);
    const std::string what = c.label + " " + lc_key_name(k);
    INFO(what);
    const lc_ref<tensor>& ref = refs.at(k);
    const lc_cplx got = lc_row(rows_world, k);
    lc_check_close(what + " [world vs one-rank correlation]", got,
                   lc_row(rows_self, k), ref.scale, lc_rtol_mpi);
    lc_check_close(what + " [world vs one-rank direct chain]", got, ref.value,
                   ref.scale, lc_rtol);
    if (it.r <= self.r_window) {
      const int g = self.window_group[it.pair];
      REQUIRE(g >= 0);
      const Bond b{it.left, std::get<1>(k), std::get<2>(k)};
      REQUIRE(window[g].count(b) == 1);
      lc_check_close(what + " [world vs one-rank window]", got, window[g].at(b),
                     ref.scale, lc_rtol);
    }
  }

  // The same rows, in the same key order, on every rank. The row count is
  // compared first, so that the element-wise reduction below is entered with
  // the same length everywhere.
  std::vector<double> count{static_cast<double>(corr_world.size())};
  std::vector<double> count_max = count;
  std::vector<double> count_min = count;
  tenes::allreduce_max(count_max, MPI_COMM_WORLD);
  tenes::allreduce_min(count_min, MPI_COMM_WORLD);
  REQUIRE(count_max[0] == count_min[0]);
  std::vector<double> values;
  std::vector<double> keys;
  for (const auto& [k, v] : rows_world) {
    keys.push_back(std::get<0>(k) + 10.0 * std::get<1>(k) +
                   100.0 * std::get<2>(k) + 1000.0 * std::get<3>(k) +
                   10000.0 * std::get<4>(k));
    for (const lc_cplx& x : v) {
      values.push_back(x.real());
      values.push_back(x.imag());
    }
  }
  std::vector<double> len{static_cast<double>(values.size()),
                          static_cast<double>(keys.size())};
  std::vector<double> len_max = len;
  std::vector<double> len_min = len;
  tenes::allreduce_max(len_max, MPI_COMM_WORLD);
  tenes::allreduce_min(len_min, MPI_COMM_WORLD);
  REQUIRE(len_max == len_min);
  std::vector<double> keys_max = keys;
  std::vector<double> keys_min = keys;
  tenes::allreduce_max(keys_max, MPI_COMM_WORLD);
  tenes::allreduce_min(keys_min, MPI_COMM_WORLD);
  CHECK(keys_max == keys_min);
  std::vector<double> values_max = values;
  std::vector<double> values_min = values;
  tenes::allreduce_max(values_max, MPI_COMM_WORLD);
  tenes::allreduce_min(values_min, MPI_COMM_WORLD);
  double max_spread = 0.0;
  for (std::size_t i = 0; i < values.size(); ++i) {
    max_spread = std::max(max_spread, values_max[i] - values_min[i]);
  }
  double max_scale = 0.0;
  for (const auto& [k, ref] : refs) {
    max_scale = std::max(max_scale, std::max(std::abs(ref.value), ref.scale));
  }
  INFO("largest spread of a correlation value over the ranks " << max_spread);
  CHECK(max_spread <= lc_rtol_mpi * max_scale);
}

}  // namespace

TEST_CASE(
    "longrange-corr T3-8: the correlation function at two MPI ranks equals "
    "the one-rank values and is the same on every rank") {
  lc_run_mpi<real_tensor>(lc_case{"T3-8 d=2 real perturbed",
                                  2,
                                  lc_names2,
                                  {{1, 2}, {2, 1}, {0, 0}, {0, 2}},
                                  3,
                                  lc_env_noise,
                                  lc_seed + 80});
  lc_run_mpi<complex_tensor>(lc_case{"T3-8 d=2 complex perturbed",
                                     2,
                                     lc_names2,
                                     {{1, 2}, {0, 0}, {1, 0}},
                                     3,
                                     lc_env_noise,
                                     lc_seed + 81});
  lc_run_mpi<real_tensor>(lc_case{"T3-8 d=4 real perturbed",
                                  4,
                                  lc_names4,
                                  {{1, 2}, {5, 5}},
                                  2,
                                  lc_env_noise,
                                  lc_seed + 82});
}

// ============================================================================
// Contract item 9: the factors of a product operator
// ============================================================================
//
// relay_product_source(A) and relay_product_target(B) are the u and vt of the
// single channel of the wrapped product A_s B_t (design section 3.4, "the
// factorization of the correlation function"). Two checks, both against T2
// and T1 code only: the channel contracted over kappa (graded) rebuilds
// wrap_twosite_gate(product_twosite_op(A, B)), as relay T1-2 rebuilds
// op12 from relay_channels(); and on open patches, the relay value of that
// one channel equals the sum over relay_channels() of the same operator.

namespace {

//! A parity-definite one-site operator with every allowed element filled
//! from the hash (complex entries for complex_tensor): odd moves between
//! the parity sectors, even stays inside them.
template <class tensor>
tensor lc_generic_op(int d, bool odd, unsigned seed) {
  const rf::parity_vector ph = lc_phys(d);
  const mptensor::Shape sh(d, d);
  tensor t(MPI_COMM_WORLD, sh);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    if ((ph[idx[0]] != ph[idx[1]]) == odd) {
      t.set_value(idx, lc_det_value(static_cast<tensor*>(nullptr), seed, 7,
                                    lc_linear(idx, sh), 1.0));
    }
  }
  return t;
}

//! Every element of a (possibly distributed) tensor in row-major order, the
//! same on every rank.
template <class tensor>
std::vector<typename tensor::value_type> lc_full(const tensor& a) {
  const mptensor::Shape sh = a.shape();
  std::vector<typename tensor::value_type> full(lc_total(sh), 0.0);
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    full[lc_linear(a.global_index(n), sh)] = a[n];
  }
  tenes::allreduce_sum(full, a.get_comm());
  return full;
}

//! One (A, B) pair of an item-9 case.
template <class tensor>
struct lc_factor_pair {
  std::string label;
  tensor A;
  tensor B;
  int ds;
  int dt;
  bool odd;
};

template <class tensor>
std::vector<lc_factor_pair<tensor>> lc_factor_pairs() {
  const MPI_Comm w = MPI_COMM_WORLD;
  const auto named = [&](const std::string& name, int d) {
    return lc_onesite_tensor<tensor>(lc_onesite_op(name, d), w);
  };
  std::vector<lc_factor_pair<tensor>> pairs;
  // Named operators, A != B.
  pairs.push_back({"d=2 (c+, c)", named("cdag", 2), named("c", 2), 2, 2, true});
  pairs.push_back({"d=2 (c, c+)", named("c", 2), named("cdag", 2), 2, 2, true});
  pairs.push_back({"d=2 (n, n)", named("n", 2), named("n", 2), 2, 2, false});
  pairs.push_back(
      {"d=4 (c+_up, c_dn)", named("cdag_up", 4), named("c_dn", 4), 4, 4, true});
  pairs.push_back(
      {"d=4 (c_dn, c+_up)", named("c_dn", 4), named("cdag_up", 4), 4, 4, true});
  pairs.push_back({"d=4 (n, Sz)", named("n", 4), named("Sz", 4), 4, 4, false});
  // Generic parity-definite operators (every allowed element non-zero).
  pairs.push_back({"d=2 generic odd", lc_generic_op<tensor>(2, true, 901),
                   lc_generic_op<tensor>(2, true, 902), 2, 2, true});
  pairs.push_back({"d=2 generic even", lc_generic_op<tensor>(2, false, 903),
                   lc_generic_op<tensor>(2, false, 904), 2, 2, false});
  pairs.push_back({"d=4 generic odd", lc_generic_op<tensor>(4, true, 905),
                   lc_generic_op<tensor>(4, true, 906), 4, 4, true});
  pairs.push_back({"d=4 generic even", lc_generic_op<tensor>(4, false, 907),
                   lc_generic_op<tensor>(4, false, 908), 4, 4, false});
  // Different physical dimensions on the two ends.
  pairs.push_back({"d=2 -> d=4 generic odd",
                   lc_generic_op<tensor>(2, true, 909),
                   lc_generic_op<tensor>(4, true, 910), 2, 4, true});
  pairs.push_back({"d=4 -> d=2 generic even",
                   lc_generic_op<tensor>(4, false, 911),
                   lc_generic_op<tensor>(2, false, 912), 4, 2, false});
  return pairs;
}

template <class tensor>
double lc_max_abs_full(const tensor& t) {
  double m = 0.0;
  for (const auto& x : lc_full(t)) {
    m = std::max(m, std::abs(x));
  }
  return m;
}

//! Contract item 9, first half: the factors' legs and ledgers, and the
//! graded contraction over kappa against the wrapped product.
template <class tensor>
void lc_run_factor_rebuild() {
  for (const auto& fp : lc_factor_pairs<tensor>()) {
    const std::string label =
        fp.label + " " + std::string(lc_type_name<tensor>());
    INFO(label);
    const rf::parity_vector ps = lc_phys(fp.ds);
    const rf::parity_vector pt = lc_phys(fp.dt);
    const rf::ftensor<tensor> op12 = rf::wrap_twosite_gate(
        rf::product_twosite_op(fp.A, fp.B, ps, pt, fp.odd), ps, pt);
    const double op_max = std::max(1.0, lc_max_abs_full(op12.t));
    {
      INFO("premise: the product is not zero");
      REQUIRE(lc_max_abs_full(op12.t) > 0.1);
    }

    rf::ftensor<tensor> u;
    rf::ftensor<tensor> vt;
    REQUIRE_NOTHROW(u = rf::relay_product_source(fp.A, ps, fp.odd));
    REQUIRE_NOTHROW(vt = rf::relay_product_target(fp.B, pt, fp.odd));
    REQUIRE(u.rank() == 3);
    REQUIRE(vt.rank() == 3);
    CHECK(u.shape() == mptensor::Shape(fp.ds, fp.ds, 1));
    CHECK(vt.shape() == mptensor::Shape(1, fp.dt, fp.dt));
    REQUIRE(u.parity.size() == 3);
    REQUIRE(vt.parity.size() == 3);
    CHECK(u.parity[0] == ps);
    CHECK(u.parity[1] == ps);
    CHECK(vt.parity[1] == pt);
    CHECK(vt.parity[2] == pt);
    REQUIRE(u.parity[2].size() == 1);
    REQUIRE(vt.parity[0].size() == 1);
    CHECK(u.parity[2][0] == fp.odd);
    CHECK(vt.parity[0][0] == fp.odd);
    // Each factor has a definite total parity (even, counting kappa).
    std::vector<double> viol{rf::parity_violation(u), rf::parity_violation(vt)};
    tenes::allreduce_max(viol, MPI_COMM_WORLD);
    CHECK(viol[0] <= 1.0e-14 * op_max);
    CHECK(viol[1] <= 1.0e-14 * op_max);

    // (in1, out1, kappa) x (kappa, in2, out2) -> (in1, out1, in2, out2),
    // then graded to (in1, in2, out1, out2), as relay T1-2 does.
    const rf::ftensor<tensor> uv =
        rf::tensordot(u, vt, mptensor::Axes(2), mptensor::Axes(0));
    const rf::ftensor<tensor> back =
        rf::transpose(uv, mptensor::Axes(0, 2, 1, 3));
    REQUIRE(back.t.shape() == op12.t.shape());
    const auto got = lc_full(back.t);
    const auto want = lc_full(op12.t);
    double max_dev = 0.0;
    for (std::size_t i = 0; i < want.size(); ++i) {
      max_dev = std::max(max_dev, std::abs(got[i] - want[i]));
    }
    INFO("max |u vt - op12| = " << max_dev << " (tol " << 1.0e-12 * op_max
                                << ")");
    CHECK(max_dev <= 1.0e-12 * op_max);
  }
}

//! An open nrow x ncol patch: every perimeter leg of dimension 1 (even),
//! inner bonds D = 2 [e, o], physical dimension d on every site, parity-even
//! entries from the hash (the same at any rank count).
template <class tensor>
std::vector<std::vector<rf::ftensor<tensor>>> lc_open_patch(int nrow, int ncol,
                                                            int d,
                                                            unsigned seed) {
  const rf::parity_vector edge{false};
  const rf::parity_vector eo{false, true};
  std::vector<std::vector<rf::ftensor<tensor>>> grid(nrow);
  for (int r = 0; r < nrow; ++r) {
    for (int c = 0; c < ncol; ++c) {
      rf::leg_parities lp{c > 0 ? eo : edge, r > 0 ? eo : edge,
                          c + 1 < ncol ? eo : edge, r + 1 < nrow ? eo : edge,
                          lc_phys(d)};
      mptensor::Shape sh;
      for (const auto& leg : lp) {
        sh.push(leg.size());
      }
      tensor t(MPI_COMM_WORLD, sh);
      for (std::size_t n = 0; n < t.local_size(); ++n) {
        const mptensor::Index idx = t.global_index(n);
        if (rf::count_odd(lp, idx) % 2 == 0) {
          t.set_value(idx, lc_det_value(static_cast<tensor*>(nullptr), seed,
                                        r * ncol + c, lc_linear(idx, sh), 1.0));
        }
      }
      grid[r].push_back(rf::ftensor<tensor>{t, lp});
    }
  }
  return grid;
}

//! The exact value of a folded open window: the density window kernel with
//! a trivial environment (every corner and edge a single 1), identity
//! operators on every cell.
template <class tensor>
typename tensor::value_type lc_contract_open(
    const std::vector<std::vector<tensor>>& w) {
  const int nrow = static_cast<int>(w.size());
  const int ncol = static_cast<int>(w[0].size());
  tensor C(MPI_COMM_WORLD, mptensor::Shape(1, 1));
  C.set_value(mptensor::Index(0, 0), 1.0);
  tensor eT(MPI_COMM_WORLD, mptensor::Shape(1, 1, 1));
  eT.set_value(mptensor::Index(0, 0, 0), 1.0);
  std::vector<std::vector<tensor>> ident(nrow);
  std::vector<std::vector<const tensor*>> cells(nrow), ops(nrow);
  for (int r = 0; r < nrow; ++r) {
    for (int c = 0; c < ncol; ++c) {
      ident[r].push_back(lc_onesite_tensor<tensor>(
          lc_onesite_op("I", static_cast<int>(w[r][c].shape()[4])),
          MPI_COMM_WORLD));
    }
  }
  for (int r = 0; r < nrow; ++r) {
    for (int c = 0; c < ncol; ++c) {
      cells[r].push_back(&w[r][c]);
      ops[r].push_back(&ident[r][c]);
    }
  }
  const std::vector<const tensor*> corners(4, &C);
  const std::vector<const tensor*> rows(nrow, &eT);
  const std::vector<const tensor*> cols(ncol, &eT);
  return tenes::itps::core::Contract_density_CTM(corners, cols, rows, cols,
                                                 rows, cells, ops);
}

//! Contract item 9, second half: on open patches, the relay value of the
//! product channel equals the sum over relay_channels() of the same
//! operator, for straight, reversed and bent paths.
template <class tensor>
void lc_run_factor_patch(int d) {
  struct patch_case {
    int nrow;
    int ncol;
    rf::window_cell src;
    rf::window_cell tgt;
  };
  const patch_case cases[] = {{1, 3, {0, 0}, {0, 2}},   // (2, 0)
                              {1, 3, {0, 2}, {0, 0}},   // (-2, 0)
                              {3, 1, {2, 0}, {0, 0}},   // (0, 2)
                              {2, 3, {1, 0}, {0, 2}},   // (2, 1), bent
                              {2, 3, {0, 2}, {1, 0}},   // (-2, -1), bent
                              {2, 2, {0, 0}, {1, 1}}};  // (1, -1)
  const rf::parity_vector ph = lc_phys(d);
  std::vector<lc_factor_pair<tensor>> pairs;
  for (auto& fp : lc_factor_pairs<tensor>()) {
    if (fp.ds == d && fp.dt == d) {
      pairs.push_back(fp);
    }
  }
  REQUIRE(!pairs.empty());
  for (const patch_case& pc : cases) {
    const auto grid =
        lc_open_patch<tensor>(pc.nrow, pc.ncol, d, lc_seed + 90 + pc.nrow);
    for (const auto& fp : pairs) {
      const std::string label =
          fp.label + " " + std::string(lc_type_name<tensor>()) + " patch " +
          std::to_string(pc.nrow) + "x" + std::to_string(pc.ncol) + " (" +
          std::to_string(pc.src.row) + "," + std::to_string(pc.src.col) +
          ") -> (" + std::to_string(pc.tgt.row) + "," +
          std::to_string(pc.tgt.col) + ")";
      INFO(label);
      // Reference: every channel of relay_channels (T1).
      const auto channels = rf::relay_channels(rf::wrap_twosite_gate(
          rf::product_twosite_op(fp.A, fp.B, ph, ph, fp.odd), ph, ph));
      REQUIRE(!channels.empty());
      typename tensor::value_type want = 0.0;
      double scale = 0.0;
      for (const auto& ch : channels) {
        const auto v =
            lc_contract_open(rf::build_relay_window(grid, pc.src, pc.tgt, ch));
        want += v;
        scale += std::abs(v);
      }
      {
        INFO("premise: the relay value carries signal: " << want);
        REQUIRE(std::abs(want) > 1.0e-6 * scale);
        REQUIRE(scale > 0.0);
      }
      rf::relay_channel<tensor> product;
      REQUIRE_NOTHROW(product.u = rf::relay_product_source(fp.A, ph, fp.odd));
      REQUIRE_NOTHROW(product.vt = rf::relay_product_target(fp.B, ph, fp.odd));
      std::vector<std::vector<tensor>> w;
      REQUIRE_NOTHROW(
          w = rf::build_relay_window(grid, pc.src, pc.tgt, product));
      const typename tensor::value_type got = lc_contract_open(w);
      lc_check_close(label + " [product channel vs relay_channels]",
                     lc_cplx(got), lc_cplx(want), scale, 1.0e-12);
    }
  }
}

}  // namespace

TEST_CASE(
    "longrange-corr T3-9a: relay_product_source (x) relay_product_target "
    "rebuilds the wrapped product A_s B_t") {
  lc_run_factor_rebuild<real_tensor>();
  lc_run_factor_rebuild<complex_tensor>();
}

TEST_CASE(
    "longrange-corr T3-9b: on open patches the product channel gives the "
    "relay value of relay_channels") {
  lc_run_factor_patch<real_tensor>(2);
  lc_run_factor_patch<real_tensor>(4);
  lc_run_factor_patch<complex_tensor>(2);
  lc_run_factor_patch<complex_tensor>(4);
}

// ============================================================================
// Contract item 11: groups defined on some sites only
// ============================================================================
//
// As for bosons (measure_correlation_ctm), a [correlation] pair whose
// one-site group is missing on a site gives no row with that site as the
// left (A) or the right (B) end; the input is accepted. Only an index below
// zero or not below the number of one-site groups is rejected.

namespace {

//! Sites on which each one-site group of the item-11 case is defined
//! (2x2 cell; group 0 = n everywhere, 1 = c+ on 0 and 3, 2 = c on 0, 1, 2).
const std::vector<std::vector<int>> lc_partial_sites = {
    {0, 1, 2, 3}, {0, 3}, {0, 1, 2}};

inline bool lc_partial_defined(int group, int site) {
  const auto& s = lc_partial_sites[group];
  return std::find(s.begin(), s.end(), site) != s.end();
}

template <class tensor>
tenes::Operators<tensor> lc_partial_onesite(int d) {
  tenes::Operators<tensor> onesite;
  for (std::size_t g = 0; g < lc_partial_sites.size(); ++g) {
    for (const int s : lc_partial_sites[g]) {
      onesite.emplace_back(lc_names2[g], static_cast<int>(g), s,
                           lc_onesite_tensor<tensor>(
                               lc_onesite_op(lc_names2[g], d), MPI_COMM_WORLD));
    }
  }
  return onesite;
}

}  // namespace

TEST_CASE(
    "longrange-corr T3-11a: pairs of groups defined on some sites only are "
    "accepted; negative and too large indices are rejected") {
  const std::vector<std::tuple<int, int>> pairs = {
      {1, 2}, {2, 1}, {0, 1}, {0, 0}};
  lc_guard_input partial(2, false, {}, 3, pairs);
  partial.onesite = lc_partial_onesite<lc_gtensor>(2);
  lc_check_accepts(partial,
                   "groups 1 (c+ on 0, 3) and 2 (c on 0, 1, 2) in pairs");

  lc_guard_input only_nowhere(2, false, {}, 2, {{0, 1}});
  only_nowhere.onesite = lc_partial_onesite<lc_gtensor>(2);
  // Group 1 is missing on sites 1 and 2: rows from there are skipped.
  lc_check_accepts(only_nowhere, "pair [0, 1] with group 1 on sites 0, 3");

  for (const auto& bad :
       std::vector<std::tuple<int, int>>{{0, 3}, {3, 0}, {-1, 0}, {0, -1}}) {
    lc_guard_input in(2, false, {}, 3, {bad});
    in.onesite = lc_partial_onesite<lc_gtensor>(2);
    lc_check_rejects(in, "pair [" + std::to_string(std::get<0>(bad)) + ", " +
                             std::to_string(std::get<1>(bad)) +
                             "] with three one-site groups");
  }
}

TEST_CASE(
    "longrange-corr T3-11b: rows whose left or right site lacks the group are "
    "left out, the others are measured as usual") {
  using tensor = real_tensor;
  const lc_case c{"T3-11b d=2 real perturbed",
                  2,
                  lc_names2,
                  {{1, 2}, {2, 1}, {0, 1}, {0, 0}, {2, 2}},
                  3,
                  lc_env_noise,
                  lc_seed + 110};
  // Reference solver: every group on every site, no [correlation].
  auto ref =
      lc_make<tensor>(c, MPI_COMM_WORLD, tenes::Operators<tensor>{}, false);
  ref->update_CTM();
  lc_perturb_env(*ref, c.env_noise, c.seed + 1);
  // Solver under test: the partial groups and the [correlation] pairs.
  const tenes::SquareLattice lattice = lc_lattice(2, 2, c.d);
  lc_state<tensor> corr(
      MPI_COMM_WORLD,
      lc_params<tensor>(lattice.N_UNIT, c.d, false,
                        "output_test_fermion_longrange_correlation"),
      lattice, tenes::EvolutionOperators<tensor>{},
      tenes::EvolutionOperators<tensor>{}, lc_partial_onesite<tensor>(c.d),
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{}, lc_corparam(c),
      tenes::itps::TransferMatrix_Parameters{});
  lc_seed_Tn(corr, lc_odd_scale, c.seed);
  corr.update_CTM();
  lc_copy_env(*ref, corr);

  const auto refs = lc_chain_refs(c, *ref);
  std::map<lc_key, bool> expected;  // key -> same parity
  int skipped_left = 0;
  int skipped_right = 0;
  for (const lc_item& it : lc_items(c, lattice.N_UNIT)) {
    const lc_key k = lc_item_key(c, it);
    const int target = lattice.other(it.left, std::get<1>(k), std::get<2>(k));
    const bool has_left = lc_partial_defined(std::get<3>(k), it.left);
    const bool has_right = lc_partial_defined(std::get<4>(k), target);
    skipped_left += has_left ? 0 : 1;
    skipped_right += (has_left && !has_right) ? 1 : 0;
    if (has_left && has_right) {
      expected[k] = lc_same_parity(c, c.pairs[it.pair]);
    }
  }
  {
    INFO("premise: rows are skipped for a missing left group ("
         << skipped_left << ") and for a missing right group (" << skipped_right
         << "), and rows remain (" << expected.size() << ")");
    REQUIRE(skipped_left > 0);
    REQUIRE(skipped_right > 0);
    REQUIRE(!expected.empty());
  }

  const auto rows = lc_rows(lc_measure(corr));
  CHECK(rows.size() == expected.size());
  for (const auto& kv : rows) {
    INFO("row " << lc_key_name(kv.first));
    CHECK(expected.count(kv.first) == 1);
  }
  for (const auto& [k, same] : expected) {
    const std::string what = c.label + " " + lc_key_name(k);
    INFO(what);
    const lc_cplx got = lc_row(rows, k);
    if (!same) {
      CHECK(got.real() == 0.0);
      CHECK(got.imag() == 0.0);
      continue;
    }
    const lc_ref<tensor>& r = refs.at(k);
    lc_check_close(what + " [correlation vs direct chain]", got, r.value,
                   r.scale, lc_rtol);
  }
}

// ============================================================================
// Design section 5.1: the measurement-side guard of correlation.operators
// ============================================================================
//
// validate_fermion_ctm_measurement() is the second guard of design section
// 5.1: an iTPS built directly (as a library caller does, without
// validate_fermion_constraints) must not let measure_correlation() index its
// operator tables with a negative index or one not below the number of
// one-site groups. The guard is called on its own first: if it lets the pair
// through, the case fails there, before measure_correlation() could index
// past the end of a vector (undefined behaviour, an abort in a Debug build)
// and take the other cases of this binary with it.

TEST_CASE(
    "longrange-corr T3-12: measure_correlation rejects a negative or too "
    "large one-site group index without validate_fermion_constraints") {
  using tensor = real_tensor;
  const lc_case base{"T3-12", 2, lc_names2, {}, 2, 0.0, lc_seed + 120};
  const auto make = [&](const std::vector<std::pair<int, int>>& pairs) {
    lc_case c = base;
    c.pairs = pairs;
    return lc_make<tensor>(c, MPI_COMM_WORLD, tenes::Operators<tensor>{}, true);
  };
  {
    // Control: in-range pairs pass the guard, a group on some sites only
    // included (contract item 11).
    auto state = make({{0, 0}, {1, 2}, {2, 0}});
    CHECK_NOTHROW(lc_acc::validate_fermion_ctm_measurement(*state));
    lc_state<tensor> partial(
        MPI_COMM_WORLD,
        lc_params<tensor>(4, 2, false,
                          "output_test_fermion_longrange_correlation"),
        lc_lattice(2, 2, 2), tenes::EvolutionOperators<tensor>{},
        tenes::EvolutionOperators<tensor>{}, lc_partial_onesite<tensor>(2),
        tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
        CorrelationParameter(2, {{1, 2}, {0, 1}}),
        tenes::itps::TransferMatrix_Parameters{});
    CHECK_NOTHROW(lc_acc::validate_fermion_ctm_measurement(partial));
  }
  for (const auto& bad : std::vector<std::pair<int, int>>{
           {0, 3}, {3, 0}, {7, 1}, {-1, 0}, {0, -1}}) {
    const std::string what = "pair [" + std::to_string(bad.first) + ", " +
                             std::to_string(bad.second) +
                             "] with three one-site groups";
    INFO(what);
    auto state = make({{0, 0}, bad});
    // The guard alone first (see above): fail here rather than let
    // measure_correlation() run into the out-of-range index.
    REQUIRE_THROWS_AS(lc_acc::validate_fermion_ctm_measurement(*state),
                      tenes::input_error);
    CHECK_THROWS_AS(state->measure_correlation(), tenes::input_error);
  }
}
