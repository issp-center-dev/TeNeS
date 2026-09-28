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
//! Task T1 of docs/superpowers/plans/2026-09-25-fermion-longrange-measure.md:
//! the string relay of src/fermion/relay.hpp (design
//! docs/superpowers/specs/2026-09-25-fermion-longrange-measure-design.md,
//! sections 3, 3.3 and 3.5). Behaviour contract:
//! work/fermion-longrange/t1/contract.md, items 1-9.
//!
//! Conventions fixed for every case in this file:
//!   - A patch is "NROW x NCOL" = rows x columns, in the Contract_* window
//!     orientation (row 0 on top, column 0 on the left). Site index
//!     s = col + ncol * row, which is fock_oracle.py's x + lx * y with
//!     x = col, y = row, lx = ncol.
//!   - A displacement (dx, dy) follows twosite_obs.cpp: dx = target.col -
//!     source.col, dy = source.row - target.row (dy > 0 is up).
//!   - The patch is open: every leg on the outer perimeter has dimension 1
//!     and an even ledger; the internal bonds carry the case's ledger.
//!   - A two-site operator is the 4-leg plain matrix
//!     op4[i_s, i_t, o_s, o_t] = sum (-1)^{p_B p(i_s)} A[i_s, o_s] B[i_t, o_t]
//!     (plan, Global Constraints; design section 4.2), source first, loaded
//!     with wrap_twosite_gate(). Leg 1 of op12 acts on the source.
//!
//! Truth sources (never the relay under test):
//!   - the single-layer graded contraction <psi|O|psi> / <psi|psi> of the
//!     same site tensors (tenes::fermion::tensordot / transpose / conj /
//!     trace; nothing folded), and
//!   - literals from test/fermion/fock_oracle.py (contract item 4).
//! The exact window contraction of the relay value is a plain mptensor
//! contraction written here (physical legs traced, fused bonds contracted);
//! it does not call any builder.
//!
//! Cases whose name carries [truth] check the reference side only and are
//! expected to pass against the stub; every other case must fail until
//! relay.hpp is implemented.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../test_fermion_common.hpp"

#include <complex>
#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "../../src/fermion/relay.hpp"

namespace {

namespace rf = tenes::fermion;

template <class tensor>
using rw_ft = rf::ftensor<tensor>;

using rf::relay_order;
using rf::relay_role;
using rf::window_cell;

// ---- small helpers ---------------------------------------------------------

inline std::string rw_cell_string(window_cell c) {
  return "(" + std::to_string(c.row) + "," + std::to_string(c.col) + ")";
}

inline bool rw_same_cell(window_cell a, window_cell b) {
  return a.row == b.row && a.col == b.col;
}

inline const char* rw_order_name(relay_order order) {
  return order == relay_order::x_first ? "x_first" : "y_first";
}

inline mptensor::Axes rw_all_axes(int rank) {
  mptensor::Axes axes;
  for (int ax = 0; ax < rank; ++ax) {
    axes.push(ax);
  }
  return axes;
}

template <class tensor>
double rw_max_abs_entry(const tensor& a) {
  double m = 0.0;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    typename tensor::value_type v;
    a.get_value(a.global_index(n), v);
    m = std::max(m, std::abs(v));
  }
  return m;
}

// Scalar judgment: |got - want| <= 1e-12 * max(|want|, scale). The scale is
// the operator's largest matrix element, a bound on |<O>| that does not
// cancel (a tolerance proportional to a cancelling result is meaningless).
template <class V>
void rw_check_close(const std::string& label, V got, V want, double scale) {
  const double tol = 1.0e-12 * std::max(std::abs(want), scale);
  INFO(label << ": got=" << got << " want=" << want
             << " |diff|=" << std::abs(got - want) << " tol=" << tol);
  CHECK(std::abs(got - want) <= tol);
}

// Elementwise: |diff| <= 1e-12 * max(1, max |element| of either tensor).
template <class tensor>
void rw_check_allclose(const tensor& got, const tensor& want,
                       const std::string& label) {
  {
    INFO(label << ": shape mismatch");
    REQUIRE(got.shape() == want.shape());
  }
  const double scale =
      std::max(1.0, std::max(rw_max_abs_entry(got), rw_max_abs_entry(want)));
  const double tol = 1.0e-12 * scale;
  double max_dev = 0.0;
  for (std::size_t n = 0; n < want.local_size(); ++n) {
    const mptensor::Index idx = want.global_index(n);
    typename tensor::value_type w, g;
    want.get_value(idx, w);
    got.get_value(idx, g);
    max_dev = std::max(max_dev, std::abs(g - w));
  }
  INFO(label << ": max |got-want| = " << max_dev << " (tol " << tol << ")");
  CHECK(max_dev <= tol);
}

// ---- deterministic parity-even site tensors -------------------------------

// Same formula as fock_oracle.deterministic_tensor(site, parities, seed) and
// fold_geometry.cpp's fg_det_entry:
// x = (site + 2) * (1 + seed + sum_ax (ax + 3 + seed % 5) * idx[ax]),
// value 0.19 sin(x) + 0.13 cos(0.37 x); only parity-even elements are set.
inline double rw_det_entry(int site, int seed, const mptensor::Index& idx) {
  double x = 1.0 + seed;
  for (int ax = 0; ax < 5; ++ax) {
    x += static_cast<double>((ax + 3 + seed % 5) * idx[ax]);
  }
  x *= site + 2;
  return 0.19 * std::sin(x) + 0.13 * std::cos(0.37 * x);
}

inline void rw_set_entry(tenes::real_tensor& t, const mptensor::Index& idx,
                         int site, int seed) {
  t.set_value(idx, rw_det_entry(site, seed, idx));
}

// Complex: an independent deterministic imaginary part (seed + 1000), as in
// fold_geometry.cpp.
inline void rw_set_entry(tenes::complex_tensor& t, const mptensor::Index& idx,
                         int site, int seed) {
  t.set_value(idx, std::complex<double>(rw_det_entry(site, seed, idx),
                                        rw_det_entry(site, seed + 1000, idx)));
}

inline rf::parity_vector rw_pv(const std::string& s) {
  rf::parity_vector v;
  for (char c : s) {
    v.push_back(c == 'o');
  }
  return v;
}

inline rf::parity_vector rw_phys(int d) {
  if (d == 2) {
    return {false, true};
  }
  if (d == 4) {
    return {false, true, true, false};
  }
  throw std::runtime_error("rw_phys: unsupported physical dimension");
}

// ---- open patch geometry --------------------------------------------------

struct rw_patch {
  int nrow;
  int ncol;
  // When nonempty, only these internal bonds (site pairs, smaller site
  // first) carry the case ledger; every other internal bond is dimension 1
  // and even. Used by the reduced 3x3 Fock-oracle patch only.
  std::vector<std::pair<int, int>> keep;

  int nsite() const { return nrow * ncol; }
  int site(int row, int col) const { return col + ncol * row; }
  int site(window_cell c) const { return site(c.row, c.col); }
  int row_of(int s) const { return s / ncol; }
  int col_of(int s) const { return s % ncol; }

  // Neighbour across virtual leg 0 = l, 1 = t, 2 = r, 3 = b; -1 if outside.
  int neighbour(int s, int leg) const {
    const int r = row_of(s);
    const int c = col_of(s);
    switch (leg) {
      case 0:
        return c > 0 ? s - 1 : -1;
      case 1:
        return r > 0 ? s - ncol : -1;
      case 2:
        return c + 1 < ncol ? s + 1 : -1;
      case 3:
        return r + 1 < nrow ? s + ncol : -1;
      default:
        return -1;
    }
  }

  bool nontrivial(int s, int leg) const {
    const int n = neighbour(s, leg);
    if (n < 0) {
      return false;
    }
    if (keep.empty()) {
      return true;
    }
    const std::pair<int, int> b{std::min(s, n), std::max(s, n)};
    return std::find(keep.begin(), keep.end(), b) != keep.end();
  }

  // Leg labels: site * 5 + leg (legs 0..3 virtual, 4 physical). Partner
  // label across an internal bond (trivial or not), -1 if none.
  int partner(int label) const {
    const int s = label / 5;
    const int leg = label % 5;
    if (leg == 4) {
      return -1;
    }
    const int n = neighbour(s, leg);
    return n < 0 ? -1 : n * 5 + (leg + 2) % 4;
  }
};

template <class tensor>
rw_ft<tensor> rw_make_site(const rw_patch& p, int s,
                           const rf::parity_vector& vp,
                           const rf::parity_vector& phys, int seed) {
  const rf::parity_vector edge{false};
  rf::leg_parities lp;
  for (int leg = 0; leg < 4; ++leg) {
    lp.push_back(p.nontrivial(s, leg) ? vp : edge);
  }
  lp.push_back(phys);
  tensor t(mptensor::Shape(lp[0].size(), lp[1].size(), lp[2].size(),
                           lp[3].size(), lp[4].size()));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    if (rf::count_odd(lp, idx) % 2 == 0) {
      rw_set_entry(t, idx, s, seed);
    }
  }
  return rw_ft<tensor>{t, lp};
}

// ---- generic network assembly ---------------------------------------------

inline tenes::real_tensor rw_dot(const tenes::real_tensor& a,
                                 const tenes::real_tensor& b,
                                 const mptensor::Axes& axes_a,
                                 const mptensor::Axes& axes_b) {
  return mptensor::tensordot(a, b, axes_a, axes_b);
}

inline tenes::complex_tensor rw_dot(const tenes::complex_tensor& a,
                                    const tenes::complex_tensor& b,
                                    const mptensor::Axes& axes_a,
                                    const mptensor::Axes& axes_b) {
  return mptensor::tensordot(a, b, axes_a, axes_b);
}

template <class tensor>
rw_ft<tensor> rw_dot(const rw_ft<tensor>& a, const rw_ft<tensor>& b,
                     const mptensor::Axes& axes_a,
                     const mptensor::Axes& axes_b) {
  return rf::tensordot(a, b, axes_a, axes_b);
}

// Contract (tensor, leg labels) nodes over every patch bond whose two ends
// are present; nodes are merged in the given order (graded for ftensor,
// plain for tensor). Remaining labels come back in result axis order.
template <class TT>
TT rw_assemble(const rw_patch& patch,
               const std::vector<std::pair<TT, std::vector<int>>>& nodes,
               std::vector<int>& labels_out) {
  TT acc = nodes[0].first;
  std::vector<int> labels = nodes[0].second;
  for (std::size_t i = 1; i < nodes.size(); ++i) {
    const auto& node = nodes[i];
    mptensor::Axes axes_a;
    mptensor::Axes axes_b;
    std::vector<bool> used_a(labels.size(), false);
    std::vector<bool> used_b(node.second.size(), false);
    for (std::size_t j = 0; j < node.second.size(); ++j) {
      const int p = patch.partner(node.second[j]);
      if (p < 0) {
        continue;
      }
      const auto it = std::find(labels.begin(), labels.end(), p);
      if (it == labels.end()) {
        continue;
      }
      const std::size_t k = static_cast<std::size_t>(it - labels.begin());
      axes_a.push(static_cast<int>(k));
      axes_b.push(static_cast<int>(j));
      used_a[k] = true;
      used_b[j] = true;
    }
    std::vector<int> new_labels;
    for (std::size_t k = 0; k < labels.size(); ++k) {
      if (!used_a[k]) {
        new_labels.push_back(labels[k]);
      }
    }
    for (std::size_t j = 0; j < node.second.size(); ++j) {
      if (!used_b[j]) {
        new_labels.push_back(node.second[j]);
      }
    }
    acc = rw_dot(acc, node.first, axes_a, axes_b);
    labels = new_labels;
  }
  labels_out = labels;
  return acc;
}

// Value of a fully contracted plain network: every remaining leg must have
// dimension 1.
template <class tensor>
typename tensor::value_type rw_scalar(const tensor& a) {
  const mptensor::Shape sh = a.shape();
  mptensor::Index idx;
  idx.resize(sh.size());
  for (std::size_t ax = 0; ax < sh.size(); ++ax) {
    REQUIRE(sh[ax] == 1);
    idx[ax] = 0;
  }
  typename tensor::value_type v;
  a.get_value(idx, v);
  return v;
}

// ---- operators -------------------------------------------------------------

// One-site operator in TeNeS' op[in, out] layout, m[in * d + out]
// = <out|A|in>. d = 4 basis: 0 = |0>, 1 = |up> = c+_up|0>,
// 2 = |dn> = c+_dn|0>, 3 = |updn> = c+_up c+_dn|0>.
struct rw_onesite {
  std::vector<double> m;
  bool odd;
};

inline rw_onesite rw_onesite_op(const std::string& name, int d) {
  rw_onesite o{std::vector<double>(d * d, 0.0), false};
  auto set = [&](int in, int out, double v) { o.m[in * d + out] = v; };
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
      throw std::runtime_error("rw_onesite_op: unknown d=2 operator " + name);
    }
    return o;
  }
  if (d != 4) {
    throw std::runtime_error("rw_onesite_op: unsupported dimension");
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
  } else {
    throw std::runtime_error("rw_onesite_op: unknown d=4 operator " + name);
  }
  return o;
}

// coef * A_s B_t.
struct rw_term {
  std::complex<double> coef;
  std::string a;
  std::string b;
};

// Complex coefficients of the non-Hermitian "cplx" operator: a conj() that
// goes missing or lands on the wrong layer changes its expectation value,
// which a real (or Hermitian) operator could hide.
const std::complex<double> rw_alpha(0.8, 0.6);
const std::complex<double> rw_beta(-0.3, 0.5);
const std::complex<double> rw_gamma(0.4, -0.7);

// Operator kinds:
//   cdagc  c+_s c_t (d = 2)
//   hop    sum_sigma (c+_s c_t + c+_t c_s), with c+_t c_s = -c_s c+_t
//   hopnn  hop + n_s n_t
//   cplx   sum_sigma (alpha c+_s c_t + beta c+_t c_s) + gamma n_s n_t
inline std::vector<rw_term> rw_op_terms(const std::string& kind, int d) {
  std::vector<std::pair<std::string, std::string>> flavours;
  if (d == 2) {
    flavours.push_back({"cdag", "c"});
  } else {
    flavours.push_back({"cdag_up", "c_up"});
    flavours.push_back({"cdag_dn", "c_dn"});
  }
  std::vector<rw_term> terms;
  if (kind == "cdagc") {
    if (d != 2) {
      throw std::runtime_error("rw_op_terms: cdagc is d = 2 only");
    }
    terms.push_back({1.0, "cdag", "c"});
    return terms;
  }
  if (kind == "hop" || kind == "hopnn") {
    for (const auto& f : flavours) {
      terms.push_back({1.0, f.first, f.second});
      terms.push_back({-1.0, f.second, f.first});
    }
    if (kind == "hopnn") {
      terms.push_back({1.0, "n", "n"});
    }
    return terms;
  }
  if (kind == "cplx") {
    for (const auto& f : flavours) {
      terms.push_back({rw_alpha, f.first, f.second});
      terms.push_back({-rw_beta, f.second, f.first});
    }
    terms.push_back({rw_gamma, "n", "n"});
    return terms;
  }
  throw std::runtime_error("rw_op_terms: unknown kind " + kind);
}

inline bool rw_phys_odd(int d, int i) { return rw_phys(d)[i]; }

// op4[i_s, i_t, o_s, o_t] = sum coef (-1)^{p_B p(i_s)} A[i_s,o_s] B[i_t,o_t],
// flattened as ((i_s * d + i_t) * d + o_s) * d + o_t.
inline std::vector<std::complex<double>> rw_op_elements(const std::string& kind,
                                                        int d) {
  std::vector<std::complex<double>> el(d * d * d * d, 0.0);
  for (const rw_term& term : rw_op_terms(kind, d)) {
    const rw_onesite A = rw_onesite_op(term.a, d);
    const rw_onesite B = rw_onesite_op(term.b, d);
    for (int is = 0; is < d; ++is) {
      for (int it = 0; it < d; ++it) {
        for (int os = 0; os < d; ++os) {
          for (int ot = 0; ot < d; ++ot) {
            const double sign = (B.odd && rw_phys_odd(d, is)) ? -1.0 : 1.0;
            el[((is * d + it) * d + os) * d + ot] +=
                term.coef * sign * A.m[is * d + os] * B.m[it * d + ot];
          }
        }
      }
    }
  }
  return el;
}

inline void rw_set_op(tenes::real_tensor& t, const mptensor::Index& idx,
                      std::complex<double> v) {
  REQUIRE(v.imag() == 0.0);
  t.set_value(idx, v.real());
}

inline void rw_set_op(tenes::complex_tensor& t, const mptensor::Index& idx,
                      std::complex<double> v) {
  t.set_value(idx, v);
}

template <class tensor>
tensor rw_op_plain(const std::string& kind, int d) {
  const auto el = rw_op_elements(kind, d);
  tensor op(mptensor::Shape(d, d, d, d));
  for (int is = 0; is < d; ++is) {
    for (int it = 0; it < d; ++it) {
      for (int os = 0; os < d; ++os) {
        for (int ot = 0; ot < d; ++ot) {
          const std::complex<double> v = el[((is * d + it) * d + os) * d + ot];
          if (v != 0.0) {
            rw_set_op(op, mptensor::Index(is, it, os, ot), v);
          }
        }
      }
    }
  }
  return op;
}

template <class tensor>
rw_ft<tensor> rw_op12(const std::string& kind, int d) {
  const rf::parity_vector phys = rw_phys(d);
  return rf::wrap_twosite_gate(rw_op_plain<tensor>(kind, d), phys, phys);
}

// ---- the state of one patch -------------------------------------------------

template <class tensor>
struct rw_state {
  using value_type = typename tensor::value_type;
  std::string label;
  rw_patch patch;
  int d = 2;
  rf::parity_vector phys;
  std::vector<rw_ft<tensor>> sites;              // by site index
  std::vector<std::vector<rw_ft<tensor>>> grid;  // [row][col]
  bool with_truth = false;
  rw_ft<tensor> psi;            // graded patch state, physical legs only
  std::vector<int> psi_labels;  // site * 5 + 4 per axis of psi
  value_type norm_graded = 0.0;
  value_type norm_window = 0.0;  // exact contraction of the identity window
};

// Graded patch state with the dimension-1 even boundary legs dropped: they
// carry no sign, so a plain reshape onto the physical legs is exact.
template <class tensor>
void rw_build_psi(rw_state<tensor>& st) {
  std::vector<std::pair<rw_ft<tensor>, std::vector<int>>> nodes;
  for (int s = 0; s < st.patch.nsite(); ++s) {
    nodes.push_back(
        {st.sites[s], {s * 5 + 0, s * 5 + 1, s * 5 + 2, s * 5 + 3, s * 5 + 4}});
  }
  std::vector<int> labels;
  const rw_ft<tensor> full = rw_assemble(st.patch, nodes, labels);
  mptensor::Shape sh;
  rf::leg_parities lp;
  std::vector<int> kept;
  for (std::size_t ax = 0; ax < labels.size(); ++ax) {
    if (labels[ax] % 5 == 4) {
      sh.push(full.shape()[ax]);
      lp.push_back(full.parity[ax]);
      kept.push_back(labels[ax]);
    } else {
      REQUIRE(full.shape()[ax] == 1);
      REQUIRE(!full.parity[ax][0]);
    }
  }
  st.psi = rw_ft<tensor>{mptensor::reshape(full.t, sh), lp};
  st.psi_labels = kept;
  const mptensor::Axes axes = rw_all_axes(st.psi.rank());
  st.norm_graded = rf::trace(rf::conj(st.psi), st.psi, axes, axes);
}

// Exact contraction of a folded window: each rank-6 cell
// ([l lb], [t tb], [r rb], [b bb], s_ket, s_bra) has its physical pair traced
// (the identity one-site operator the kernels are given), then the fused
// bonds are contracted plainly. Shapes are REQUIREd first so that a wrong
// shape fails this case instead of tripping an mptensor assertion.
template <class tensor>
typename tensor::value_type rw_contract_window(
    const rw_patch& p, const std::vector<std::vector<tensor>>& w,
    const std::string& label) {
  INFO(label);
  REQUIRE(w.size() == static_cast<std::size_t>(p.nrow));
  for (int r = 0; r < p.nrow; ++r) {
    REQUIRE(w[r].size() == static_cast<std::size_t>(p.ncol));
    for (int c = 0; c < p.ncol; ++c) {
      INFO("cell " << rw_cell_string({r, c}));
      REQUIRE(w[r][c].shape().size() == 6);
      REQUIRE(w[r][c].shape()[4] == w[r][c].shape()[5]);
    }
  }
  for (int s = 0; s < p.nsite(); ++s) {
    const auto& t = w[p.row_of(s)][p.col_of(s)];
    for (int leg = 0; leg < 4; ++leg) {
      const int n = p.neighbour(s, leg);
      INFO("cell " << rw_cell_string({p.row_of(s), p.col_of(s)}) << " leg "
                   << leg);
      if (n < 0) {
        REQUIRE(t.shape()[leg] == 1);
      } else {
        const auto& u = w[p.row_of(n)][p.col_of(n)];
        REQUIRE(t.shape()[leg] == u.shape()[(leg + 2) % 4]);
      }
    }
  }
  std::vector<std::pair<tensor, std::vector<int>>> nodes;
  for (int s = 0; s < p.nsite(); ++s) {
    const tensor traced = mptensor::contract(
        w[p.row_of(s)][p.col_of(s)], mptensor::Axes(4), mptensor::Axes(5));
    nodes.push_back({traced, {s * 5 + 0, s * 5 + 1, s * 5 + 2, s * 5 + 3}});
  }
  std::vector<int> rest;
  return rw_scalar(rw_assemble(p, nodes, rest));
}

template <class tensor>
std::vector<std::vector<tensor>> rw_identity_window(
    const rw_state<tensor>& st) {
  std::vector<std::vector<tensor>> w(st.patch.nrow);
  for (int r = 0; r < st.patch.nrow; ++r) {
    for (int c = 0; c < st.patch.ncol; ++c) {
      w[r].push_back(
          rf::detail::doubled_pipeline(st.grid[r][c], st.grid[r][c]));
    }
  }
  return w;
}

template <class tensor>
rw_state<tensor> rw_make_state(const rw_patch& p, int d, const std::string& vps,
                               int seed, bool with_truth) {
  rw_state<tensor> st;
  st.label = "patch " + std::to_string(p.nrow) + "x" + std::to_string(p.ncol) +
             " (rows x cols) d=" + std::to_string(d) + " vp=" + vps +
             " seed=" + std::to_string(seed) +
             (p.keep.empty() ? "" : " reduced-bond") +
             (std::is_same<tensor, tenes::complex_tensor>::value ? " complex"
                                                                 : " real");
  st.patch = p;
  st.d = d;
  st.phys = rw_phys(d);
  const rf::parity_vector vp = rw_pv(vps);
  for (int s = 0; s < p.nsite(); ++s) {
    st.sites.push_back(rw_make_site<tensor>(p, s, vp, st.phys, seed));
  }
  st.grid.resize(p.nrow);
  for (int r = 0; r < p.nrow; ++r) {
    for (int c = 0; c < p.ncol; ++c) {
      st.grid[r].push_back(st.sites[p.site(r, c)]);
    }
  }
  st.norm_window = rw_contract_window(p, rw_identity_window(st),
                                      st.label + " [identity window]");
  REQUIRE(std::abs(st.norm_window) > 0.0);
  st.with_truth = with_truth;
  if (with_truth) {
    rw_build_psi(st);
  }
  return st;
}

// ---- single-layer graded truth --------------------------------------------

// Unnormalized <psi| gate(a, b) |psi>: gate legs (in_a, in_b, out_a, out_b).
template <class tensor>
typename tensor::value_type rw_graded_apply(const rw_state<tensor>& st, int a,
                                            int b, const rw_ft<tensor>& gate) {
  const rw_ft<tensor>& psi = st.psi;
  const int rank = psi.rank();
  int ax_a = -1;
  int ax_b = -1;
  for (std::size_t k = 0; k < st.psi_labels.size(); ++k) {
    if (st.psi_labels[k] == a * 5 + 4) {
      ax_a = static_cast<int>(k);
    }
    if (st.psi_labels[k] == b * 5 + 4) {
      ax_b = static_cast<int>(k);
    }
  }
  REQUIRE(ax_a >= 0);
  REQUIRE(ax_b >= 0);
  rw_ft<tensor> applied = rf::tensordot(psi, gate, mptensor::Axes(ax_a, ax_b),
                                        mptensor::Axes(0, 1));
  // applied: psi's free legs in order, then (out_a, out_b).
  mptensor::Axes perm;
  int run = 0;
  for (int ax = 0; ax < rank; ++ax) {
    if (ax == ax_a) {
      perm.push(rank - 2);
    } else if (ax == ax_b) {
      perm.push(rank - 1);
    } else {
      perm.push(run++);
    }
  }
  applied = rf::transpose(applied, perm);
  const mptensor::Axes axes = rw_all_axes(rank);
  return rf::trace(rf::conj(psi), applied, axes, axes);
}

// Normalized <O> with op12's leg 1 on the source. Site tensors stay in
// geometric (= Jordan-Wigner) order; a source that comes later is expressed
// on the gate by the graded transpose (1,0,3,2), the convention of the
// existing nearest-neighbour path (fold_geometry.cpp, source=second).
template <class tensor>
typename tensor::value_type rw_graded_value(const rw_state<tensor>& st,
                                            window_cell src, window_cell tgt,
                                            const rw_ft<tensor>& op12) {
  REQUIRE(st.with_truth);
  const int s = st.patch.site(src);
  const int t = st.patch.site(tgt);
  if (s < t) {
    return rw_graded_apply(st, s, t, op12) / st.norm_graded;
  }
  return rw_graded_apply(st, t, s,
                         rf::transpose(op12, mptensor::Axes(1, 0, 3, 2))) /
         st.norm_graded;
}

// ---- relay value under test ----------------------------------------------

// sum_k [exact contraction of window k] / [exact contraction of the identity
// window] (contract item 3).
template <class tensor>
typename tensor::value_type rw_relay_value(
    const rw_state<tensor>& st, window_cell src, window_cell tgt,
    const std::vector<rf::relay_channel<tensor>>& channels, relay_order order,
    const std::string& label) {
  REQUIRE(!channels.empty());
  typename tensor::value_type sum = 0.0;
  for (std::size_t k = 0; k < channels.size(); ++k) {
    const auto w =
        rf::build_relay_window(st.grid, src, tgt, channels[k], order);
    sum += rw_contract_window(
        st.patch, w,
        label + " channel " + std::to_string(k) + " " + rw_order_name(order));
  }
  return sum / st.norm_window;
}

// ---- nearest-neighbour bundled-k reference ----------------------------------

// Existing bundled-k construction (build_reduced_pair_halves, A = left or
// top site, source-second through the graded transpose (1,0,3,2) of the
// gate), closed exactly with the rest of the patch's reduced tensors.
template <class tensor>
typename tensor::value_type rw_bundled_value(const rw_state<tensor>& st,
                                             window_cell src, window_cell tgt,
                                             const rw_ft<tensor>& op12) {
  const rw_patch& p = st.patch;
  REQUIRE(std::abs(src.row - tgt.row) + std::abs(src.col - tgt.col) == 1);
  const bool horizontal = src.row == tgt.row;
  const bool source_first = horizontal ? src.col < tgt.col : src.row < tgt.row;
  const window_cell A = source_first ? src : tgt;
  const window_cell B = source_first ? tgt : src;
  const rw_ft<tensor> gate =
      source_first ? op12 : rf::transpose(op12, mptensor::Axes(1, 0, 3, 2));
  const auto dir = horizontal ? rf::reduced_pair_direction::horizontal
                              : rf::reduced_pair_direction::vertical;
  const auto halves = rf::build_reduced_pair_halves(
      st.grid[A.row][A.col], st.grid[B.row][B.col], gate, dir);
  const int a = p.site(A);
  const int b = p.site(B);
  std::vector<std::pair<tensor, std::vector<int>>> nodes;
  for (int s = 0; s < p.nsite(); ++s) {
    const std::vector<int> legs{s * 5 + 0, s * 5 + 1, s * 5 + 2, s * 5 + 3};
    if (s == a) {
      nodes.push_back({halves.PA, legs});
    } else if (s == b) {
      nodes.push_back({halves.PB, legs});
    } else {
      nodes.push_back({rf::build_reduced(st.sites[s]), legs});
    }
  }
  std::vector<int> rest;
  return rw_scalar(rw_assemble(p, nodes, rest)) / st.norm_window;
}

// ---- the case table of contract item 3 -------------------------------------

struct rw_case {
  const char* name;
  int nrow;
  int ncol;
  window_cell src;
  window_cell tgt;
  int dx() const { return tgt.col - src.col; }
  int dy() const { return src.row - tgt.row; }
  bool bend() const { return dx() != 0 && dy() != 0; }
  std::string label() const {
    return std::string(name) + " [" + std::to_string(nrow) + "x" +
           std::to_string(ncol) + " rows x cols, src " + rw_cell_string(src) +
           " -> tgt " + rw_cell_string(tgt) + ", (dx,dy)=(" +
           std::to_string(dx()) + "," + std::to_string(dy()) + ")]";
  }
};

// Contract item 3, in full, followed by extra 3x3 cases (marked extra).
// The contract's "3x2 (2,1)" needs three columns: it is the 2 rows x 3
// columns patch here.
const rw_case rw_cases[] = {
    {"straight (2,0)", 1, 3, {0, 0}, {0, 2}},
    {"straight (3,0)", 1, 4, {0, 0}, {0, 3}},
    {"straight (0,2)", 3, 1, {2, 0}, {0, 0}},
    {"straight (-2,0)", 1, 3, {0, 2}, {0, 0}},
    {"corner (1,1)", 2, 2, {1, 0}, {0, 1}},
    {"corner (1,-1)", 2, 2, {0, 0}, {1, 1}},
    {"corner (-1,1)", 2, 2, {1, 1}, {0, 0}},
    {"corner (-1,-1)", 2, 2, {0, 1}, {1, 0}},
    {"corner (2,1)", 2, 3, {1, 0}, {0, 2}},
    {"corner (-2,1)", 2, 3, {1, 2}, {0, 0}},
    {"3x3 straight through the centre: top-mid -> bottom-mid",
     3,
     3,
     {0, 1},
     {2, 1}},
    {"3x3 bend at the centre: left-mid -> top-mid", 3, 3, {1, 0}, {0, 1}},
    {"extra 3x3 straight through the centre: right-mid -> left-mid",
     3,
     3,
     {1, 2},
     {1, 0}},
    {"extra 3x3 corner to corner: bottom-left -> top-right",
     3,
     3,
     {2, 0},
     {0, 2}},
    {"extra 3x3 source at the centre: centre -> bottom-right",
     3,
     3,
     {1, 1},
     {2, 2}},
    {"extra 3x3 target at the centre: top-right -> centre",
     3,
     3,
     {0, 2},
     {1, 1}},
    {"extra 3x3 bend at the centre (y_first): bottom-mid -> right-mid",
     3,
     3,
     {2, 1},
     {1, 2}},
    {"extra 3x3 bend at the centre (x_first): left-mid -> bottom-mid",
     3,
     3,
     {1, 0},
     {2, 1}},
};

// Cases of the thinner sweeps (D = 3 ledger, d = 4 complex): the (1,1)
// corner, the (2,1) corner with a three-leg middle site, the two mandatory
// 3x3 paths and the corner-to-corner path.
const int rw_subset[] = {4, 8, 10, 11, 13};

template <class tensor>
using rw_state_cache = std::map<std::pair<int, int>, rw_state<tensor>>;

template <class tensor>
rw_state<tensor>& rw_cached_state(rw_state_cache<tensor>& cache, int nrow,
                                  int ncol, int d, const std::string& vps,
                                  int seed, bool with_truth) {
  const std::pair<int, int> key{nrow, ncol};
  auto it = cache.find(key);
  if (it == cache.end()) {
    it = cache
             .emplace(key, rw_make_state<tensor>(rw_patch{nrow, ncol, {}}, d,
                                                 vps, seed, with_truth))
             .first;
  }
  return it->second;
}

// Contract item 3: relay value (default order = x_first) vs graded truth.
template <class tensor>
void rw_run_exact(const rw_case& c, int d, const std::string& vps, int seed,
                  const std::vector<std::string>& kinds,
                  rw_state_cache<tensor>& cache) {
  rw_state<tensor>& st =
      rw_cached_state(cache, c.nrow, c.ncol, d, vps, seed, true);
  for (const std::string& kind : kinds) {
    const std::string label = c.label() + " op=" + kind + " " + st.label;
    INFO(label);
    const rw_ft<tensor> op12 = rw_op12<tensor>(kind, d);
    const double op_max = rw_max_abs_entry(op12.t);
    const auto truth = rw_graded_value(st, c.src, c.tgt, op12);
    const auto channels = rf::relay_channels(op12);
    REQUIRE(!channels.empty());
    typename tensor::value_type sum = 0.0;
    for (std::size_t k = 0; k < channels.size(); ++k) {
      // Default argument: the production order is x_first.
      const auto w = rf::build_relay_window(st.grid, c.src, c.tgt, channels[k]);
      sum += rw_contract_window(st.patch, w,
                                label + " channel " + std::to_string(k));
    }
    rw_check_close(label + " [relay (default order) vs graded truth]",
                   sum / st.norm_window, truth, op_max);
  }
}

// Contract item 5: x_first and y_first agree.
template <class tensor>
void rw_run_order(const rw_case& c, int d, const std::string& vps, int seed,
                  const std::vector<std::string>& kinds,
                  rw_state_cache<tensor>& cache) {
  rw_state<tensor>& st =
      rw_cached_state(cache, c.nrow, c.ncol, d, vps, seed, false);
  for (const std::string& kind : kinds) {
    const std::string label = c.label() + " op=" + kind + " " + st.label;
    INFO(label);
    const rw_ft<tensor> op12 = rw_op12<tensor>(kind, d);
    const double op_max = rw_max_abs_entry(op12.t);
    const auto channels = rf::relay_channels(op12);
    const auto vx =
        rw_relay_value(st, c.src, c.tgt, channels, relay_order::x_first, label);
    const auto vy =
        rw_relay_value(st, c.src, c.tgt, channels, relay_order::y_first, label);
    rw_check_close(label + " [y_first vs x_first]", vy, vx, op_max);
  }
}

// ---- nearest-neighbour cases (contract item 6) ---------------------------

struct rw_nn_case {
  const char* name;
  int nrow;
  int ncol;
  window_cell src;
  window_cell tgt;
};

const rw_nn_case rw_nn_cases[] = {
    {"2x2 source left", 2, 2, {0, 0}, {0, 1}},
    {"2x2 source right", 2, 2, {1, 1}, {1, 0}},
    {"2x2 source top", 2, 2, {0, 1}, {1, 1}},
    {"2x2 source bottom", 2, 2, {1, 0}, {0, 0}},
    {"3x3 centre, source left: centre -> right-mid", 3, 3, {1, 1}, {1, 2}},
    {"3x3 centre, source right: centre -> left-mid", 3, 3, {1, 1}, {1, 0}},
    {"3x3 centre, source top: centre -> bottom-mid", 3, 3, {1, 1}, {2, 1}},
    {"3x3 centre, source bottom: centre -> top-mid", 3, 3, {1, 1}, {0, 1}},
};

inline std::string rw_nn_label(const rw_nn_case& c) {
  return std::string(c.name) + " [src " + rw_cell_string(c.src) + " -> tgt " +
         rw_cell_string(c.tgt) + "]";
}

template <class tensor>
void rw_run_nn(const rw_nn_case& c, int d,
               const std::vector<std::string>& kinds,
               rw_state_cache<tensor>& cache, bool truth_only) {
  rw_state<tensor>& st =
      rw_cached_state(cache, c.nrow, c.ncol, d, "eo", 0, truth_only);
  for (const std::string& kind : kinds) {
    const std::string label = rw_nn_label(c) + " op=" + kind + " " + st.label;
    INFO(label);
    const rw_ft<tensor> op12 = rw_op12<tensor>(kind, d);
    const double op_max = rw_max_abs_entry(op12.t);
    const auto bundled = rw_bundled_value(st, c.src, c.tgt, op12);
    if (truth_only) {
      rw_check_close(label + " [bundled-k vs graded truth]", bundled,
                     rw_graded_value(st, c.src, c.tgt, op12), op_max);
      continue;
    }
    const auto channels = rf::relay_channels(op12);
    for (const relay_order order :
         {relay_order::x_first, relay_order::y_first}) {
      rw_check_close(
          label + " [relay " + rw_order_name(order) + " vs bundled-k]",
          rw_relay_value(st, c.src, c.tgt, channels, order, label), bundled,
          op_max);
    }
  }
}

// ---- Fock-oracle anchors (contract item 4) ----------------------------------
//
// Generated from test/fermion/fock_oracle.py (unchanged), run inside
// test/fermion/ with `python3 - <<'EOF'` and the script below (seed 0,
// internal bonds [e, o] = [False, True], physical [e, o]):
//
//   from fock_oracle import Oracle, Patch, deterministic_tensor, make_case
//
//   class Sub(Patch):
//       # Only the bonds in KEEP carry [e, o]; every other leg is dim 1, even.
//       def internal_bonds(self):
//           return [b for b in Patch.internal_bonds(self) if b in KEEP]
//
//   KEEP = [(0, "r", 1, "l"), (0, "b", 3, "t"), (1, "b", 4, "t"),
//           (3, "r", 4, "l"), (4, "r", 5, "l"), (4, "b", 7, "t")]
//
//   def report(name, oracle, pairs):
//       n = oracle.norm()
//       print("%s norm = %.17e" % (name, n))
//       for s, t in pairs:
//           print("%s one_body(%d, %d) / norm = %.17e"
//                 % (name, s, t, oracle.one_body(s, t) / n))
//
//   for lx, ly, pairs in [(3, 1, [(0, 2)]), (2, 2, [(2, 1), (0, 3)]),
//                         (3, 2, [(3, 2), (5, 0)])]:
//       patch, tensors, lp = make_case(lx, ly, [False, True], 0)
//       report("make_case(%d, %d)" % (lx, ly), Oracle(patch, tensors, lp),
//              pairs)
//
//   eo, e = [False, True], [False]
//   lp = [[e, e, e, e, eo] for _ in range(9)]
//   for a, al, b, bl in KEEP:
//       lp[a]["ltrb".index(al)] = eo
//       lp[b]["ltrb".index(bl)] = eo
//   tensors = [deterministic_tensor(s, lp[s], 0) for s in range(9)]
//   report("Sub(3, 3)", Oracle(Sub(3, 3), tensors, lp), [(3, 1)])
//
// one_body(i, j) = <psi| c+_i c_j |psi> with site index x + lx * y, so
// one_body(s, t) / norm is <c+_s c_t> for source s and target t.
//
// The 3x3 anchor cannot use the full [e, o] patch: the oracle keeps one
// Fock mode per site and two per internal bond, 9 + 2 * 12 = 33 modes, a
// dense vector of 2^33 doubles. Only the six KEEP bonds carry [e, o] (21
// modes, about 10 s); the rest are dimension 1 and even, which the oracle
// needs no mode for. The centre still has four nontrivial legs, and the
// top-left plaquette is a loop: x_first bends at the centre (entry l, exit
// t), y_first bends at the corner (0,0) (entry b, exit r).

struct rw_anchor {
  const char* name;
  int nrow;
  int ncol;
  bool reduced;
  window_cell src;
  window_cell tgt;
  double norm;
  double value;
};

const rw_anchor rw_anchors[] = {
    {"(2,0) 1x3, make_case(3, 1) one_body(0, 2)",
     1,
     3,
     false,
     {0, 0},
     {0, 2},
     2.67265746667580693e-05,
     -1.10806860734770890e-01},
    {"(1,1) 2x2, make_case(2, 2) one_body(2, 1)",
     2,
     2,
     false,
     {1, 0},
     {0, 1},
     9.97335335163249894e-06,
     1.42176220745016541e-01},
    {"(1,-1) 2x2, make_case(2, 2) one_body(0, 3)",
     2,
     2,
     false,
     {0, 0},
     {1, 1},
     9.97335335163249894e-06,
     3.28573406583586441e-02},
    {"(2,1) 2x3, make_case(3, 2) one_body(3, 2)",
     2,
     3,
     false,
     {1, 0},
     {0, 2},
     3.51258625440347983e-08,
     5.31061262706638906e-02},
    {"(-2,1) 2x3, make_case(3, 2) one_body(5, 0)",
     2,
     3,
     false,
     {1, 2},
     {0, 0},
     3.51258625440347983e-08,
     2.94524622779692054e-02},
    {"3x3 bend (1,1) left-mid -> top-mid, Sub(3, 3) one_body(3, 1)",
     3,
     3,
     true,
     {1, 0},
     {0, 1},
     2.00317454667758503e-15,
     -8.83440955729557498e-02},
};

inline rw_patch rw_anchor_patch(const rw_anchor& a) {
  rw_patch p{a.nrow, a.ncol, {}};
  if (a.reduced) {
    p.keep = {{0, 1}, {0, 3}, {1, 4}, {3, 4}, {4, 5}, {4, 7}};
  }
  return p;
}

// ---- expected path, written out independently of relay_path ----------------

inline std::vector<window_cell> rw_expected_path(window_cell src,
                                                 window_cell tgt,
                                                 relay_order order) {
  std::vector<window_cell> path{src};
  window_cell cur = src;
  auto walk_col = [&]() {
    while (cur.col != tgt.col) {
      cur.col += tgt.col > cur.col ? 1 : -1;
      path.push_back(cur);
    }
  };
  auto walk_row = [&]() {
    while (cur.row != tgt.row) {
      cur.row += tgt.row > cur.row ? 1 : -1;
      path.push_back(cur);
    }
  };
  if (order == relay_order::x_first) {
    walk_col();
    walk_row();
  } else {
    walk_row();
    walk_col();
  }
  return path;
}

inline int rw_expected_leg(window_cell from, window_cell to) {
  if (to.row == from.row && to.col == from.col - 1) {
    return 0;
  }
  if (to.row == from.row - 1 && to.col == from.col) {
    return 1;
  }
  if (to.row == from.row && to.col == from.col + 1) {
    return 2;
  }
  if (to.row == from.row + 1 && to.col == from.col) {
    return 3;
  }
  return -1;
}

// A channel built here, from the graded SVD directly, so that the input
// checks of build_relay_site do not depend on relay_channels.
template <class tensor>
rf::relay_channel<tensor> rw_manual_channel(const rw_ft<tensor>& op12,
                                            bool want_odd) {
  rw_ft<tensor> u, vt;
  std::vector<double> s;
  REQUIRE(rf::svd(op12, mptensor::Axes(0, 2), mptensor::Axes(1, 3), u, s, vt) ==
          0);
  u.multiply_vector(s, 2);
  for (std::size_t k = 0; k < s.size(); ++k) {
    if (s[k] > 1.0e-12 && u.parity[2][k] == want_odd) {
      return {rf::slice(u, 2, k, k + 1), rf::slice(vt, 0, k, k + 1)};
    }
  }
  FAIL("rw_manual_channel: no channel of the requested parity");
  return {};
}

}  // namespace

// ============================================================================
// Contract item 1: relay_path and relay_leg
// ============================================================================

TEST_CASE("relay T1-1a: relay_path walks x_first / y_first as specified") {
  // Literal paths first (the bend of contract item 3's 3x3 case).
  {
    const auto px = rf::relay_path({1, 0}, {0, 1}, relay_order::x_first);
    REQUIRE(px.size() == 3);
    CHECK(rw_same_cell(px[0], {1, 0}));
    CHECK(rw_same_cell(px[1], {1, 1}));
    CHECK(rw_same_cell(px[2], {0, 1}));
    const auto py = rf::relay_path({1, 0}, {0, 1}, relay_order::y_first);
    REQUIRE(py.size() == 3);
    CHECK(rw_same_cell(py[0], {1, 0}));
    CHECK(rw_same_cell(py[1], {0, 0}));
    CHECK(rw_same_cell(py[2], {0, 1}));
  }
  // Every ordered pair of distinct cells of a 4x4 window, both orders.
  for (int sr = 0; sr < 4; ++sr) {
    for (int sc = 0; sc < 4; ++sc) {
      for (int tr = 0; tr < 4; ++tr) {
        for (int tc = 0; tc < 4; ++tc) {
          if (sr == tr && sc == tc) {
            continue;
          }
          const window_cell src{sr, sc};
          const window_cell tgt{tr, tc};
          for (const relay_order order :
               {relay_order::x_first, relay_order::y_first}) {
            INFO("src " << rw_cell_string(src) << " tgt " << rw_cell_string(tgt)
                        << " " << rw_order_name(order));
            const auto path = rf::relay_path(src, tgt, order);
            const std::size_t len = static_cast<std::size_t>(
                std::abs(tr - sr) + std::abs(tc - sc) + 1);
            REQUIRE(path.size() == len);
            CHECK(rw_same_cell(path.front(), src));
            CHECK(rw_same_cell(path.back(), tgt));
            for (std::size_t i = 1; i < path.size(); ++i) {
              CHECK(std::abs(path[i].row - path[i - 1].row) +
                        std::abs(path[i].col - path[i - 1].col) ==
                    1);
            }
            // x_first: along the source row to the target column first;
            // y_first: along the source column to the target row first.
            const std::size_t turn = static_cast<std::size_t>(
                order == relay_order::x_first ? std::abs(tc - sc)
                                              : std::abs(tr - sr));
            for (std::size_t i = 0; i < path.size(); ++i) {
              if (order == relay_order::x_first) {
                CHECK(path[i].row == (i <= turn ? sr : path[i].row));
                CHECK(path[i].col == (i >= turn ? tc : path[i].col));
              } else {
                CHECK(path[i].col == (i <= turn ? sc : path[i].col));
                CHECK(path[i].row == (i >= turn ? tr : path[i].row));
              }
            }
            const auto want = rw_expected_path(src, tgt, order);
            for (std::size_t i = 0; i < path.size(); ++i) {
              CHECK(rw_same_cell(path[i], want[i]));
            }
          }
        }
      }
    }
  }
}

TEST_CASE("relay T1-1b: relay_path rejects source == target") {
  for (const window_cell c :
       {window_cell{0, 0}, window_cell{1, 2}, window_cell{3, 3}}) {
    for (const relay_order order :
         {relay_order::x_first, relay_order::y_first}) {
      INFO("cell " << rw_cell_string(c) << " " << rw_order_name(order));
      CHECK_THROWS_AS(rf::relay_path(c, c, order), std::invalid_argument);
    }
  }
}

TEST_CASE("relay T1-1c: relay_leg returns l/t/r/b = 0/1/2/3") {
  const window_cell from{1, 1};
  CHECK(rf::relay_leg(from, {1, 2}) == 2);
  CHECK(rf::relay_leg(from, {1, 0}) == 0);
  CHECK(rf::relay_leg(from, {0, 1}) == 1);
  CHECK(rf::relay_leg(from, {2, 1}) == 3);
  // Opposite ends of a bond see each other through opposite legs.
  for (const window_cell to : {window_cell{1, 2}, window_cell{1, 0},
                               window_cell{0, 1}, window_cell{2, 1}}) {
    INFO("to " << rw_cell_string(to));
    CHECK(rf::relay_leg(to, from) == (rf::relay_leg(from, to) + 2) % 4);
  }
}

TEST_CASE("relay T1-1d: relay_leg rejects non-adjacent cells") {
  const window_cell from{1, 1};
  for (const window_cell to :
       {window_cell{1, 1}, window_cell{0, 0}, window_cell{2, 2},
        window_cell{0, 2}, window_cell{1, 3}, window_cell{3, 1}}) {
    INFO("to " << rw_cell_string(to));
    CHECK_THROWS_AS(rf::relay_leg(from, to), std::invalid_argument);
  }
}

// ============================================================================
// Contract item 2: relay_channels
// ============================================================================

namespace {

template <class tensor>
void rw_run_channels(const std::string& kind, int d, bool want_even,
                     bool want_odd) {
  const std::string label =
      "op=" + kind + " d=" + std::to_string(d) +
      (std::is_same<tensor, tenes::complex_tensor>::value ? " complex"
                                                          : " real");
  INFO(label);
  const rf::parity_vector phys = rw_phys(d);
  const rw_ft<tensor> op12 = rw_op12<tensor>(kind, d);
  const double op_max = rw_max_abs_entry(op12.t);
  const auto channels = rf::relay_channels(op12);
  REQUIRE(!channels.empty());
  int n_even = 0;
  int n_odd = 0;
  rw_ft<tensor> sum;
  for (std::size_t k = 0; k < channels.size(); ++k) {
    INFO("channel " << k);
    const auto& ch = channels[k];
    REQUIRE(ch.u.rank() == 3);
    REQUIRE(ch.vt.rank() == 3);
    CHECK(ch.u.shape() == mptensor::Shape(d, d, 1));
    CHECK(ch.vt.shape() == mptensor::Shape(1, d, d));
    REQUIRE(ch.u.parity[2].size() == 1);
    REQUIRE(ch.vt.parity[0].size() == 1);
    CHECK(ch.u.parity[2][0] == ch.vt.parity[0][0]);
    CHECK(ch.u.parity[0] == phys);
    CHECK(ch.u.parity[1] == phys);
    CHECK(ch.vt.parity[1] == phys);
    CHECK(ch.vt.parity[2] == phys);
    // Each factor has a definite total parity (even, counting kappa).
    CHECK(rf::parity_violation(ch.u) <= 1.0e-14 * op_max);
    CHECK(rf::parity_violation(ch.vt) <= 1.0e-14 * op_max);
    const rw_ft<tensor> uv =
        rf::tensordot(ch.u, ch.vt, mptensor::Axes(2), mptensor::Axes(0));
    if (rw_max_abs_entry(uv.t) > 1.0e-12 * op_max) {
      if (ch.u.parity[2][0]) {
        ++n_odd;
      } else {
        ++n_even;
      }
    }
    if (k == 0) {
      sum = uv;
    } else {
      REQUIRE(uv.shape() == sum.shape());
      sum.t += uv.t;
    }
  }
  // (in1, out1, in2, out2) -> (in1, in2, out1, out2), graded.
  const rw_ft<tensor> back = rf::transpose(sum, mptensor::Axes(0, 2, 1, 3));
  rw_check_allclose(back.t, op12.t, label + " [sum_k u_k vt_k vs op12]");
  if (want_even) {
    CHECK(n_even >= 1);
  }
  if (want_odd) {
    CHECK(n_odd >= 1);
  }
}

}  // namespace

TEST_CASE(
    "relay T1-2: relay_channels splits op12 into dimension-1 channels that "
    "sum back to it") {
  // d = 2 hopping is odd x odd only (both terms move one fermion across);
  // hopping + nn adds the even channel.
  rw_run_channels<tenes::real_tensor>("hop", 2, false, true);
  rw_run_channels<tenes::real_tensor>("hopnn", 2, true, true);
  rw_run_channels<tenes::real_tensor>("hopnn", 4, true, true);
  rw_run_channels<tenes::complex_tensor>("hop", 2, false, true);
  rw_run_channels<tenes::complex_tensor>("cplx", 2, true, true);
  rw_run_channels<tenes::complex_tensor>("cplx", 4, true, true);
}

// ============================================================================
// Contract item 3: exact window contraction vs graded truth
// ============================================================================

TEST_CASE(
    "relay T1-3a: relay equals graded truth, d=2 real D=2, all cases "
    "(c+_s c_t, hopping)") {
  rw_state_cache<tenes::real_tensor> cache;
  for (const rw_case& c : rw_cases) {
    rw_run_exact<tenes::real_tensor>(c, 2, "eo", 0, {"cdagc", "hop"}, cache);
  }
}

TEST_CASE(
    "relay T1-3b: relay equals graded truth, d=2 real D=3 [e,o,o] seed 173, "
    "subset (hopping)") {
  rw_state_cache<tenes::real_tensor> cache;
  for (const int i : rw_subset) {
    rw_run_exact<tenes::real_tensor>(rw_cases[i], 2, "eoo", 173, {"hop"},
                                     cache);
  }
}

TEST_CASE(
    "relay T1-3c: relay equals graded truth, d=2 complex D=2, all cases "
    "(c+_s c_t, complex non-Hermitian)") {
  rw_state_cache<tenes::complex_tensor> cache;
  for (const rw_case& c : rw_cases) {
    rw_run_exact<tenes::complex_tensor>(c, 2, "eo", 0, {"cdagc", "cplx"},
                                        cache);
  }
}

TEST_CASE(
    "relay T1-3d: relay equals graded truth, d=4 real D=2, all cases "
    "(hopping + nn)") {
  rw_state_cache<tenes::real_tensor> cache;
  for (const rw_case& c : rw_cases) {
    rw_run_exact<tenes::real_tensor>(c, 4, "eo", 0, {"hopnn"}, cache);
  }
}

TEST_CASE(
    "relay T1-3e: relay equals graded truth, d=4 complex D=2, subset "
    "(complex non-Hermitian hopping + nn)") {
  rw_state_cache<tenes::complex_tensor> cache;
  for (const int i : rw_subset) {
    rw_run_exact<tenes::complex_tensor>(rw_cases[i], 4, "eo", 0, {"cplx"},
                                        cache);
  }
}

namespace {

// build_relay_window's cells are build_relay_site() on the path (roles and
// legs from the documented path) and doubled_pipeline(Tn, Tn) elsewhere.
template <class tensor>
void rw_run_composition(const rw_case& c, int d, const std::string& kind,
                        relay_order order, bool use_default) {
  const rw_patch p{c.nrow, c.ncol, {}};
  const rw_state<tensor> st = rw_make_state<tensor>(p, d, "eo", 0, false);
  const std::string label =
      c.label() + " op=" + kind + " " + rw_order_name(order) +
      (use_default ? " (default argument)" : "") + " " + st.label;
  INFO(label);
  const auto channels = rf::relay_channels(rw_op12<tensor>(kind, d));
  REQUIRE(!channels.empty());
  const auto path = rw_expected_path(c.src, c.tgt, order);
  for (std::size_t k = 0; k < channels.size(); ++k) {
    INFO("channel " << k);
    const auto w =
        use_default
            ? rf::build_relay_window(st.grid, c.src, c.tgt, channels[k])
            : rf::build_relay_window(st.grid, c.src, c.tgt, channels[k], order);
    REQUIRE(w.size() == static_cast<std::size_t>(c.nrow));
    for (int r = 0; r < c.nrow; ++r) {
      REQUIRE(w[r].size() == static_cast<std::size_t>(c.ncol));
      for (int col = 0; col < c.ncol; ++col) {
        const window_cell cell{r, col};
        std::size_t on = path.size();
        for (std::size_t i = 0; i < path.size(); ++i) {
          if (rw_same_cell(path[i], cell)) {
            on = i;
          }
        }
        tensor want;
        std::string what;
        if (on == path.size()) {
          want = rf::detail::doubled_pipeline(st.grid[r][col], st.grid[r][col]);
          what = "off path: doubled_pipeline(Tn, Tn)";
        } else {
          const std::size_t last = path.size() - 1;
          const relay_role role = on == 0      ? relay_role::source
                                  : on == last ? relay_role::target
                                               : relay_role::middle;
          const int entry = on == 0 ? -1 : rw_expected_leg(cell, path[on - 1]);
          const int exit =
              on == last ? -1 : rw_expected_leg(cell, path[on + 1]);
          want = rf::build_relay_site(st.grid[r][col], role, entry, exit,
                                      channels[k]);
          what = "on path: build_relay_site(entry " + std::to_string(entry) +
                 ", exit " + std::to_string(exit) + ")";
        }
        rw_check_allclose(w[r][col], want,
                          label + " cell " + rw_cell_string(cell) + " " + what);
      }
    }
  }
}

}  // namespace

TEST_CASE(
    "relay T1-3f: build_relay_window is build_relay_site on the documented "
    "path and doubled_pipeline elsewhere") {
  // 3x3 bend: x_first bends at the centre, y_first at the corner (0,0).
  const rw_case& bend = rw_cases[11];
  rw_run_composition<tenes::real_tensor>(bend, 2, "hopnn", relay_order::x_first,
                                         true);
  rw_run_composition<tenes::real_tensor>(bend, 2, "hopnn", relay_order::x_first,
                                         false);
  rw_run_composition<tenes::real_tensor>(bend, 2, "hopnn", relay_order::y_first,
                                         false);
  // (-2,1) on 2 rows x 3 columns: source right of and below the target,
  // complex channels.
  rw_run_composition<tenes::complex_tensor>(rw_cases[9], 2, "cplx",
                                            relay_order::x_first, false);
  rw_run_composition<tenes::complex_tensor>(rw_cases[9], 2, "cplx",
                                            relay_order::y_first, false);
}

// ============================================================================
// Contract item 4: Fock-oracle anchors
// ============================================================================

TEST_CASE("relay T1-4: relay <c+_s c_t> equals the Fock oracle anchors") {
  for (const rw_anchor& a : rw_anchors) {
    const rw_state<tenes::real_tensor> st = rw_make_state<tenes::real_tensor>(
        rw_anchor_patch(a), 2, "eo", 0, false);
    const std::string label = std::string("anchor ") + a.name + " " + st.label;
    INFO(label);
    const auto op12 = rw_op12<tenes::real_tensor>("cdagc", 2);
    const auto channels = rf::relay_channels(op12);
    for (const relay_order order :
         {relay_order::x_first, relay_order::y_first}) {
      rw_check_close(label + " [relay " + rw_order_name(order) + " vs oracle]",
                     rw_relay_value(st, a.src, a.tgt, channels, order, label),
                     a.value, 1.0);
    }
  }
}

// ============================================================================
// Contract item 5: path independence on every bend
// ============================================================================

TEST_CASE("relay T1-5a: x_first equals y_first on every bend, d=2 real") {
  rw_state_cache<tenes::real_tensor> cache;
  for (const rw_case& c : rw_cases) {
    if (c.bend()) {
      rw_run_order<tenes::real_tensor>(c, 2, "eo", 0, {"cdagc", "hop"}, cache);
    }
  }
}

TEST_CASE("relay T1-5b: x_first equals y_first on every bend, d=2 complex") {
  rw_state_cache<tenes::complex_tensor> cache;
  for (const rw_case& c : rw_cases) {
    if (c.bend()) {
      rw_run_order<tenes::complex_tensor>(c, 2, "eo", 0, {"cdagc", "cplx"},
                                          cache);
    }
  }
}

TEST_CASE("relay T1-5c: x_first equals y_first on every bend, d=4 real") {
  rw_state_cache<tenes::real_tensor> cache;
  for (const rw_case& c : rw_cases) {
    if (c.bend()) {
      rw_run_order<tenes::real_tensor>(c, 4, "eo", 0, {"hopnn"}, cache);
    }
  }
}

TEST_CASE("relay T1-5d: x_first equals y_first on every bend, d=2 real D=3") {
  rw_state_cache<tenes::real_tensor> cache;
  for (const int i : rw_subset) {
    if (rw_cases[i].bend()) {
      rw_run_order<tenes::real_tensor>(rw_cases[i], 2, "eoo", 173, {"hop"},
                                       cache);
    }
  }
}

// ============================================================================
// Contract item 6: nearest neighbours vs the bundled-k construction
// ============================================================================

TEST_CASE("relay T1-6a: nearest neighbours equal bundled-k, d=2 real") {
  rw_state_cache<tenes::real_tensor> cache;
  for (const rw_nn_case& c : rw_nn_cases) {
    rw_run_nn<tenes::real_tensor>(c, 2, {"cdagc", "hop", "hopnn"}, cache,
                                  false);
  }
}

TEST_CASE("relay T1-6b: nearest neighbours equal bundled-k, d=2 complex") {
  rw_state_cache<tenes::complex_tensor> cache;
  for (const rw_nn_case& c : rw_nn_cases) {
    rw_run_nn<tenes::complex_tensor>(c, 2, {"cdagc", "cplx"}, cache, false);
  }
}

TEST_CASE("relay T1-6c: nearest neighbours equal bundled-k, d=4 real") {
  rw_state_cache<tenes::real_tensor> cache;
  for (const rw_nn_case& c : rw_nn_cases) {
    rw_run_nn<tenes::real_tensor>(c, 4, {"hopnn"}, cache, false);
  }
}

TEST_CASE("relay T1-6d: nearest neighbours equal bundled-k, d=4 complex 2x2") {
  rw_state_cache<tenes::complex_tensor> cache;
  for (int i = 0; i < 4; ++i) {
    rw_run_nn<tenes::complex_tensor>(rw_nn_cases[i], 4, {"cplx"}, cache, false);
  }
}

// ============================================================================
// Contract item 7: input checks of build_relay_site
// ============================================================================

namespace {

template <class tensor>
void rw_run_input_checks() {
  const rf::parity_vector phys = rw_phys(2);
  // A freestanding site: every virtual leg [e, o], so every leg is valid.
  const rw_patch p{3, 3, {}};
  const rw_ft<tensor> Tn = rw_make_site<tensor>(p, 4, rw_pv("eo"), phys, 0);
  const auto op12 = rw_op12<tensor>("hop", 2);
  const rf::relay_channel<tensor> ch = rw_manual_channel(op12, true);

  // Consistent combinations are accepted and give the rank-6 fold.
  for (const int exit : {0, 1, 2, 3}) {
    INFO("source exit " << exit);
    CHECK_NOTHROW(rf::build_relay_site(Tn, relay_role::source, -1, exit, ch));
  }
  for (const int entry : {0, 1, 2, 3}) {
    INFO("target entry " << entry);
    CHECK_NOTHROW(rf::build_relay_site(Tn, relay_role::target, entry, -1, ch));
  }
  for (const auto& io : std::vector<std::pair<int, int>>{
           {0, 2}, {2, 0}, {1, 3}, {3, 1}, {0, 1}, {3, 2}}) {
    INFO("middle entry " << io.first << " exit " << io.second);
    CHECK_NOTHROW(
        rf::build_relay_site(Tn, relay_role::middle, io.first, io.second, ch));
  }

  // source with entry_leg != -1.
  for (const int entry : {0, 1, 2, 3}) {
    INFO("source entry " << entry);
    CHECK_THROWS_AS(rf::build_relay_site(Tn, relay_role::source, entry,
                                         (entry + 2) % 4, ch),
                    std::invalid_argument);
  }
  // target with exit_leg != -1.
  for (const int exit : {0, 1, 2, 3}) {
    INFO("target exit " << exit);
    CHECK_THROWS_AS(
        rf::build_relay_site(Tn, relay_role::target, (exit + 2) % 4, exit, ch),
        std::invalid_argument);
  }
  // middle with entry == exit.
  for (const int leg : {0, 1, 2, 3}) {
    INFO("middle entry == exit == " << leg);
    CHECK_THROWS_AS(rf::build_relay_site(Tn, relay_role::middle, leg, leg, ch),
                    std::invalid_argument);
  }
  // middle with a leg out of 0..3 (including -1).
  for (const auto& io : std::vector<std::pair<int, int>>{
           {-1, 2}, {2, -1}, {4, 2}, {2, 4}, {-2, 0}, {0, 7}}) {
    INFO("middle entry " << io.first << " exit " << io.second);
    CHECK_THROWS_AS(
        rf::build_relay_site(Tn, relay_role::middle, io.first, io.second, ch),
        std::invalid_argument);
  }
}

}  // namespace

TEST_CASE(
    "relay T1-7: build_relay_site throws std::invalid_argument on an "
    "inconsistent role / leg combination") {
  rw_run_input_checks<tenes::real_tensor>();
  rw_run_input_checks<tenes::complex_tensor>();
}

// ============================================================================
// [truth] reference-side self-checks (expected green against the stub)
// ============================================================================

TEST_CASE(
    "relay [truth] identity window contraction equals the graded norm on "
    "every patch") {
  const std::pair<int, int> shapes[] = {{1, 3}, {1, 4}, {3, 1},
                                        {2, 2}, {2, 3}, {3, 3}};
  for (const auto& sh : shapes) {
    const rw_patch p{sh.first, sh.second, {}};
    for (const int d : {2, 4}) {
      const auto st = rw_make_state<tenes::real_tensor>(p, d, "eo", 0, true);
      rw_check_close(st.label + " [identity window vs graded norm]",
                     st.norm_window, st.norm_graded, std::abs(st.norm_graded));
      const auto stc =
          rw_make_state<tenes::complex_tensor>(p, d, "eo", 0, true);
      rw_check_close(stc.label + " [identity window vs graded norm]",
                     stc.norm_window, stc.norm_graded,
                     std::abs(stc.norm_graded));
    }
    const auto st3 = rw_make_state<tenes::real_tensor>(p, 2, "eoo", 173, true);
    rw_check_close(st3.label + " [identity window vs graded norm]",
                   st3.norm_window, st3.norm_graded, std::abs(st3.norm_graded));
  }
}

TEST_CASE("relay [truth] graded truth equals the Fock oracle anchors") {
  for (const rw_anchor& a : rw_anchors) {
    const auto st =
        rw_make_state<tenes::real_tensor>(rw_anchor_patch(a), 2, "eo", 0, true);
    const std::string label = std::string("anchor ") + a.name + " " + st.label;
    rw_check_close(label + " [graded norm vs oracle norm]", st.norm_graded,
                   a.norm, std::abs(a.norm));
    rw_check_close(label + " [identity window vs oracle norm]", st.norm_window,
                   a.norm, std::abs(a.norm));
    const auto op12 = rw_op12<tenes::real_tensor>("cdagc", 2);
    rw_check_close(label + " [graded <c+_s c_t> vs oracle]",
                   rw_graded_value(st, a.src, a.tgt, op12), a.value, 1.0);
  }
}

TEST_CASE(
    "relay [truth] the case table covers negative displacements, bends and "
    "every signal is large enough to see a sign error") {
  bool src_right = false;
  bool src_below = false;
  for (const rw_case& c : rw_cases) {
    src_right = src_right || c.src.col > c.tgt.col;
    src_below = src_below || c.src.row > c.tgt.row;
  }
  CHECK(src_right);
  CHECK(src_below);
  // Signal floor: |<O>| >= 1e-3 * max|op element| for every (case, op) of
  // contract item 3, so that a flipped or dropped odd channel (a change of
  // the order of |<O>|) is far outside the 1e-12 tolerance. The odd part
  // alone (hopping) is checked for the mixed operators.
  auto floor_check = [](const std::string& label, double v, double op_max) {
    INFO(label << ": |<O>|=" << v << " floor=" << 1.0e-3 * op_max);
    CHECK(v >= 1.0e-3 * op_max);
  };
  {
    rw_state_cache<tenes::real_tensor> c2, c4, c3;
    rw_state_cache<tenes::complex_tensor> z2;
    for (const rw_case& c : rw_cases) {
      auto& s2 = rw_cached_state(c2, c.nrow, c.ncol, 2, "eo", 0, true);
      for (const char* kind : {"cdagc", "hop"}) {
        const auto op = rw_op12<tenes::real_tensor>(kind, 2);
        floor_check(c.label() + " d=2 real " + kind,
                    std::abs(rw_graded_value(s2, c.src, c.tgt, op)),
                    rw_max_abs_entry(op.t));
      }
      auto& z = rw_cached_state(z2, c.nrow, c.ncol, 2, "eo", 0, true);
      for (const char* kind : {"cdagc", "cplx"}) {
        const auto op = rw_op12<tenes::complex_tensor>(kind, 2);
        floor_check(c.label() + " d=2 complex " + kind,
                    std::abs(rw_graded_value(z, c.src, c.tgt, op)),
                    rw_max_abs_entry(op.t));
      }
      auto& s4 = rw_cached_state(c4, c.nrow, c.ncol, 4, "eo", 0, true);
      for (const char* kind : {"hopnn", "hop"}) {
        const auto op = rw_op12<tenes::real_tensor>(kind, 4);
        floor_check(c.label() + " d=4 real " + kind,
                    std::abs(rw_graded_value(s4, c.src, c.tgt, op)),
                    rw_max_abs_entry(op.t));
      }
    }
    for (const int i : rw_subset) {
      const rw_case& c = rw_cases[i];
      auto& s3 = rw_cached_state(c3, c.nrow, c.ncol, 2, "eoo", 173, true);
      const auto op = rw_op12<tenes::real_tensor>("hop", 2);
      floor_check(c.label() + " d=2 real D=3 hop",
                  std::abs(rw_graded_value(s3, c.src, c.tgt, op)),
                  rw_max_abs_entry(op.t));
    }
  }
}

TEST_CASE(
    "relay [truth] graded truth: applying op12 at (source, target) directly "
    "equals the source-second transpose convention") {
  // Both forms are graded-covariant rewrites of one another; the anchors
  // decide which is right, this pins that the truth does not depend on it.
  rw_state_cache<tenes::complex_tensor> cache;
  for (const rw_case& c : rw_cases) {
    auto& st = rw_cached_state(cache, c.nrow, c.ncol, 2, "eo", 0, true);
    const auto op12 = rw_op12<tenes::complex_tensor>("cplx", 2);
    const int s = st.patch.site(c.src);
    const int t = st.patch.site(c.tgt);
    const auto direct = rw_graded_apply(st, s, t, op12) / st.norm_graded;
    rw_check_close(c.label() + " [direct vs source-second transpose]", direct,
                   rw_graded_value(st, c.src, c.tgt, op12),
                   rw_max_abs_entry(op12.t));
  }
}

TEST_CASE(
    "relay [truth] nearest-neighbour bundled-k reference equals the graded "
    "truth") {
  rw_state_cache<tenes::real_tensor> r2, r4;
  rw_state_cache<tenes::complex_tensor> z2, z4;
  for (const rw_nn_case& c : rw_nn_cases) {
    rw_run_nn<tenes::real_tensor>(c, 2, {"cdagc", "hop", "hopnn"}, r2, true);
    rw_run_nn<tenes::complex_tensor>(c, 2, {"cdagc", "cplx"}, z2, true);
    rw_run_nn<tenes::real_tensor>(c, 4, {"hopnn"}, r4, true);
  }
  for (int i = 0; i < 4; ++i) {
    rw_run_nn<tenes::complex_tensor>(rw_nn_cases[i], 4, {"cplx"}, z4, true);
  }
}

TEST_CASE(
    "relay [truth] the source-first product operators reproduce the "
    "hopping tables of fold_geometry.cpp") {
  // d = 2: c+_1 c_2 + c+_2 c_1 has exactly (0,1,1,0) = (1,0,0,1) = 1.
  {
    const auto hop = rw_op_plain<tenes::real_tensor>("hop", 2);
    tenes::real_tensor want(mptensor::Shape(2, 2, 2, 2));
    want.set_value(mptensor::Index(0, 1, 1, 0), 1.0);
    want.set_value(mptensor::Index(1, 0, 0, 1), 1.0);
    rw_check_allclose(hop, want, "d=2 hop");
  }
  // d = 4: the Jordan-Wigner table of fold_geometry.cpp's ib_hop_plain(4)
  // (c_A = c x 1, c_B = P x c), written out again here.
  {
    tenes::real_tensor want(mptensor::Shape(4, 4, 4, 4));
    double c[2][4][4] = {};  // c[s][out][in]
    c[0][0][1] = 1.0;
    c[0][2][3] = 1.0;
    c[1][0][2] = 1.0;
    c[1][1][3] = -1.0;
    const double ps[4] = {1.0, -1.0, -1.0, 1.0};
    for (int i1 = 0; i1 < 4; ++i1) {
      for (int i2 = 0; i2 < 4; ++i2) {
        for (int o1 = 0; o1 < 4; ++o1) {
          for (int o2 = 0; o2 < 4; ++o2) {
            double v = 0.0;
            for (int s = 0; s < 2; ++s) {
              v += c[s][i1][o1] * ps[i1] * c[s][o2][i2];
              v += ps[o1] * c[s][o1][i1] * c[s][i2][o2];
            }
            if (v != 0.0) {
              want.set_value(mptensor::Index(i1, i2, o1, o2), v);
            }
          }
        }
      }
    }
    rw_check_allclose(rw_op_plain<tenes::real_tensor>("hop", 4), want,
                      "d=4 hop");
  }
}
