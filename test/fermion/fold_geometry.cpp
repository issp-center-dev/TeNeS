// Tests for the fermionic folding geometry: the folded (double-layer
// reduced) contraction must agree with the single-layer graded contraction
// on arbitrary open patches, including the former blind-spot shapes where a
// site carries three or more nontrivial legs (2x3, 3x2, 3x3).
//
// Contract: work/fermion/ctm-fold-fix/CONTRACT.md (rev. 2.1).
// This file is included into the test_fermion_layer TU. The ib_* plain
// operator builders below were moved here VERBATIM from the retired
// impurity_blob.cpp (2026-08-28; its direct==naive equivalence suite is
// superseded by T13, which is grounded in the graded truth). The contract
// pins the d = 4 hopping to that ib_hop_plain(4), so the moved bodies are
// bit-identical to the originals.
//
// Shape naming convention (fixed for every test in this file):
// "LXxLY" means lx columns times ly rows, site = x + lx * y (the Fock
// oracle's convention, test/fermion/fock_oracle.py make_case(lx, ly, ...)).
// The patch is open: every leg on the outer perimeter has dimension 1 with
// an even ledger; interior bonds carry the case's virtual ledger.
//
// Truth sources (never the folded machinery under test):
//   - the single-layer graded contraction of the same site tensors
//     (tenes::fermion::tensordot / trace / conj / transpose), or
//   - literals burned in from test/fermion/fock_oracle.py (T5).

namespace {

namespace fgf = tenes::fermion;

template <class tensor>
using fg_ftensor = fgf::ftensor<tensor>;

// ---- plain operator builders (moved verbatim from impurity_blob.cpp) ------

tenes::fermion::parity_vector ib_phys_parity(int d) {
  if (d == 2) {
    return {false, true};
  }
  if (d == 4) {
    return {false, true, true, false};
  }
  throw std::runtime_error("ib_phys_parity: unsupported physical dimension");
}

double ib_random_value(int seed, const mptensor::Index& idx, std::size_t rank) {
  const double w[5] = {3.1, 5.3, 7.9, 11.7, 13.1};
  double x = 0.7 * (seed + 1) + 0.013 * seed * seed;
  for (std::size_t ax = 0; ax < rank; ++ax) {
    x += w[ax] * static_cast<double>(idx[ax]);
  }
  return 0.5 * std::sin(x) + 0.25 * std::cos(1.7 * x + 0.3 * seed);
}

tenes::real_tensor ib_identity_plain(int d) {
  tenes::real_tensor op(mptensor::Shape(d, d, d, d));
  for (int a = 0; a < d; ++a) {
    for (int b = 0; b < d; ++b) {
      op.set_value(mptensor::Index(a, b, a, b), 1.0);
    }
  }
  return op;
}

tenes::real_tensor ib_hop_plain(int d) {
  tenes::real_tensor op(mptensor::Shape(d, d, d, d));
  if (d == 2) {
    // c^dag_A c_B + c^dag_B c_A, element table from the contract.
    op.set_value(mptensor::Index(0, 1, 1, 0), 1.0);
    op.set_value(mptensor::Index(1, 0, 0, 1), 1.0);
    return op;
  }
  // d = 4 electron, basis 0=|0>, 1=|up>, 2=|dn>, 3=|updn>:
  // sum_s (c^dag_As c_Bs + c^dag_Bs c_As) with the Jordan-Wigner matrix
  // representation c_A = c x 1, c_B = P x c that reproduces the d = 2
  // table above in the spinless case. c[s][out][in].
  double c[2][4][4] = {};
  c[0][0][1] = 1.0;   // c_up |up> -> |0>
  c[0][2][3] = 1.0;   // c_up |updn> -> |dn>
  c[1][0][2] = 1.0;   // c_dn |dn> -> |0>
  c[1][1][3] = -1.0;  // c_dn |updn> -> -|up>
  const double ps[4] = {1.0, -1.0, -1.0, 1.0};
  for (int i1 = 0; i1 < 4; ++i1) {
    for (int i2 = 0; i2 < 4; ++i2) {
      for (int o1 = 0; o1 < 4; ++o1) {
        for (int o2 = 0; o2 < 4; ++o2) {
          double v = 0.0;
          for (int s = 0; s < 2; ++s) {
            // (c^dag_s P) x c_s : c^dag_A c_B
            v += c[s][i1][o1] * ps[i1] * c[s][o2][i2];
            // (P c_s) x c^dag_s : c^dag_B c_A
            v += ps[o1] * c[s][o1][i1] * c[s][i2][o2];
          }
          if (v != 0.0) {
            op.set_value(mptensor::Index(i1, i2, o1, o2), v);
          }
        }
      }
    }
  }
  return op;
}

tenes::real_tensor ib_nn_plain(int d) {
  tenes::real_tensor op(mptensor::Shape(d, d, d, d));
  if (d == 2) {
    op.set_value(mptensor::Index(1, 1, 1, 1), 1.0);
    return op;
  }
  const double nv[4] = {0.0, 1.0, 1.0, 2.0};
  for (int i1 = 0; i1 < 4; ++i1) {
    for (int i2 = 0; i2 < 4; ++i2) {
      if (nv[i1] * nv[i2] != 0.0) {
        op.set_value(mptensor::Index(i1, i2, i1, i2), nv[i1] * nv[i2]);
      }
    }
  }
  return op;
}

tenes::real_tensor ib_random_even_plain(int d, int seed) {
  const tenes::fermion::parity_vector phys = ib_phys_parity(d);
  tenes::real_tensor op(mptensor::Shape(d, d, d, d));
  for (std::size_t n = 0; n < op.local_size(); ++n) {
    const mptensor::Index idx = op.global_index(n);
    const bool odd =
        ((phys[idx[0]] != phys[idx[1]]) != (phys[idx[2]] != phys[idx[3]]));
    if (!odd) {
      op.set_value(idx, ib_random_value(seed, idx, 4));
    }
  }
  return op;
}

std::string ib_index_string(const mptensor::Index& idx) {
  std::string s = "(";
  for (std::size_t i = 0; i < idx.size(); ++i) {
    if (i > 0) {
      s += ",";
    }
    s += std::to_string(idx[i]);
  }
  s += ")";
  return s;
}

// ---- deterministic parity-even site tensors -------------------------------

// Same formula as fock_oracle.deterministic_tensor(site, parities, seed):
// x = (site + 2) * (1 + seed + sum_ax (ax + 3 + seed % 5) * idx[ax]),
// value 0.19 sin(x) + 0.13 cos(0.37 x); only parity-even elements are set.
inline double fg_det_entry(int site, int seed, const mptensor::Index& idx) {
  double x = 1.0 + seed;
  for (int ax = 0; ax < 5; ++ax) {
    x += static_cast<double>((ax + 3 + seed % 5) * idx[ax]);
  }
  x *= site + 2;
  return 0.19 * std::sin(x) + 0.13 * std::cos(0.37 * x);
}

inline void fg_set_entry(tenes::real_tensor& t, const mptensor::Index& idx,
                         int site, int seed) {
  t.set_value(idx, fg_det_entry(site, seed, idx));
}

// Complex variant: an independent deterministic imaginary part (seed
// offset 1000), still parity-even because only even elements are filled.
inline void fg_set_entry(tenes::complex_tensor& t, const mptensor::Index& idx,
                         int site, int seed) {
  t.set_value(idx, std::complex<double>(fg_det_entry(site, seed, idx),
                                        fg_det_entry(site, seed + 1000, idx)));
}

// ---- open patch geometry --------------------------------------------------

struct fg_patch {
  int lx;
  int ly;
  int nsite() const { return lx * ly; }
  int site(int x, int y) const { return x + lx * y; }
  int x_of(int s) const { return s % lx; }
  int y_of(int s) const { return s / lx; }

  // Leg labels: site * 5 + leg, legs 0=l, 1=t, 2=r, 3=b, 4=phys.
  int partner(int label) const {
    const int s = label / 5;
    const int leg = label % 5;
    const int x = x_of(s);
    const int y = y_of(s);
    switch (leg) {
      case 0:
        return x > 0 ? (s - 1) * 5 + 2 : -1;
      case 1:
        return y > 0 ? (s - lx) * 5 + 3 : -1;
      case 2:
        return x + 1 < lx ? (s + 1) * 5 + 0 : -1;
      case 3:
        return y + 1 < ly ? (s + lx) * 5 + 1 : -1;
      default:
        return -1;
    }
  }
};

fgf::parity_vector fg_pv(const std::string& s) {
  fgf::parity_vector v;
  for (char c : s) {
    v.push_back(c == 'o');
  }
  return v;
}

fgf::leg_parities fg_site_parity(const fg_patch& p, int s,
                                 const fgf::parity_vector& vp,
                                 const fgf::parity_vector& phys) {
  const fgf::parity_vector edge{false};
  const int x = p.x_of(s);
  const int y = p.y_of(s);
  return {x > 0 ? vp : edge, y > 0 ? vp : edge, x + 1 < p.lx ? vp : edge,
          y + 1 < p.ly ? vp : edge, phys};
}

template <class tensor>
fg_ftensor<tensor> fg_make_site(const fg_patch& p, int s,
                                const fgf::parity_vector& vp,
                                const fgf::parity_vector& phys, int seed) {
  const fgf::leg_parities lp = fg_site_parity(p, s, vp, phys);
  tensor t(mptensor::Shape(lp[0].size(), lp[1].size(), lp[2].size(),
                           lp[3].size(), lp[4].size()));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    if (fgf::count_odd(lp, idx) % 2 == 0) {
      fg_set_entry(t, idx, s, seed);
    }
  }
  return fg_ftensor<tensor>{t, lp};
}

// Freestanding single-bond site (T6/T7/T8/T13): every virtual leg carries
// the same nontrivial ledger vp (no dimension-1 boundary).
template <class tensor>
fg_ftensor<tensor> fg_make_bond_site(int site_id, const fgf::parity_vector& vp,
                                     const fgf::parity_vector& phys, int seed) {
  const fgf::leg_parities lp{vp, vp, vp, vp, phys};
  tensor t(
      mptensor::Shape(vp.size(), vp.size(), vp.size(), vp.size(), phys.size()));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    if (fgf::count_odd(lp, idx) % 2 == 0) {
      fg_set_entry(t, idx, site_id, seed);
    }
  }
  return fg_ftensor<tensor>{t, lp};
}

// ---- generic network assembly ---------------------------------------------

inline tenes::real_tensor fg_dot(const tenes::real_tensor& a,
                                 const tenes::real_tensor& b,
                                 const mptensor::Axes& axes_a,
                                 const mptensor::Axes& axes_b) {
  return mptensor::tensordot(a, b, axes_a, axes_b);
}

inline tenes::complex_tensor fg_dot(const tenes::complex_tensor& a,
                                    const tenes::complex_tensor& b,
                                    const mptensor::Axes& axes_a,
                                    const mptensor::Axes& axes_b) {
  return mptensor::tensordot(a, b, axes_a, axes_b);
}

template <class tensor>
fg_ftensor<tensor> fg_dot(const fg_ftensor<tensor>& a,
                          const fg_ftensor<tensor>& b,
                          const mptensor::Axes& axes_a,
                          const mptensor::Axes& axes_b) {
  return fgf::tensordot(a, b, axes_a, axes_b);
}

// Contract a list of (tensor, leg labels) nodes over all patch bonds whose
// two ends are present; nodes are merged in the given order. Remaining
// labels (boundary and physical legs) are returned through labels_out in
// the axis order of the result.
template <class TT>
TT fg_assemble(const fg_patch& patch,
               const std::vector<std::pair<TT, std::vector<int>>>& nodes,
               std::vector<int>& labels_out) {
  TT acc = nodes[0].first;
  std::vector<int> labels = nodes[0].second;
  for (std::size_t i = 1; i < nodes.size(); ++i) {
    const auto& node = nodes[i];
    mptensor::Axes axes_a;
    mptensor::Axes axes_b;
    std::vector<bool> contracted(node.second.size(), false);
    for (std::size_t j = 0; j < node.second.size(); ++j) {
      const int p = patch.partner(node.second[j]);
      if (p < 0) {
        continue;
      }
      const auto it = std::find(labels.begin(), labels.end(), p);
      if (it == labels.end()) {
        continue;
      }
      axes_a.push(static_cast<int>(it - labels.begin()));
      axes_b.push(static_cast<int>(j));
      contracted[j] = true;
    }
    std::vector<int> new_labels;
    for (std::size_t k = 0; k < labels.size(); ++k) {
      bool used = false;
      for (std::size_t m = 0; m < axes_a.size(); ++m) {
        if (axes_a[m] == k) {
          used = true;
        }
      }
      if (!used) {
        new_labels.push_back(labels[k]);
      }
    }
    for (std::size_t j = 0; j < node.second.size(); ++j) {
      if (!contracted[j]) {
        new_labels.push_back(node.second[j]);
      }
    }
    acc = fg_dot(acc, node.first, axes_a, axes_b);
    labels = new_labels;
  }
  labels_out = labels;
  return acc;
}

// Value of a fully contracted plain network: every remaining leg must have
// dimension 1.
template <class tensor>
typename tensor::value_type fg_scalar(const tensor& a) {
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

inline mptensor::Axes fg_all_axes(int rank) {
  mptensor::Axes axes;
  for (int ax = 0; ax < rank; ++ax) {
    axes.push(ax);
  }
  return axes;
}

// ---- single-layer graded truths -------------------------------------------

template <class tensor>
fg_ftensor<tensor> fg_patch_state(const fg_patch& p,
                                  const std::vector<fg_ftensor<tensor>>& sites,
                                  std::vector<int>& labels_out) {
  std::vector<std::pair<fg_ftensor<tensor>, std::vector<int>>> nodes;
  for (int s = 0; s < p.nsite(); ++s) {
    nodes.push_back(
        {sites[s], {s * 5 + 0, s * 5 + 1, s * 5 + 2, s * 5 + 3, s * 5 + 4}});
  }
  return fg_assemble(p, nodes, labels_out);
}

template <class tensor>
typename tensor::value_type fg_graded_norm(const fg_ftensor<tensor>& psi) {
  const mptensor::Axes axes = fg_all_axes(psi.rank());
  return fgf::trace(fgf::conj(psi), psi, axes, axes);
}

// Unnormalized <psi| gate(a, b) |psi>; the gate carries legs
// (in_a, in_b, out_a, out_b) in geometric (= JW site) order and is already
// wrapped (wrap_twosite_gate, plus the source-second transpose if any).
template <class tensor>
typename tensor::value_type fg_graded_twosite(const fg_ftensor<tensor>& psi,
                                              const std::vector<int>& labels,
                                              int site_a, int site_b,
                                              const fg_ftensor<tensor>& gate) {
  const int rank = psi.rank();
  int ax_a = -1;
  int ax_b = -1;
  for (std::size_t k = 0; k < labels.size(); ++k) {
    if (labels[k] == site_a * 5 + 4) {
      ax_a = static_cast<int>(k);
    }
    if (labels[k] == site_b * 5 + 4) {
      ax_b = static_cast<int>(k);
    }
  }
  REQUIRE(ax_a >= 0);
  REQUIRE(ax_b >= 0);
  fg_ftensor<tensor> applied = fgf::tensordot(
      psi, gate, mptensor::Axes(ax_a, ax_b), mptensor::Axes(0, 1));
  // applied legs: psi's free legs in original relative order, then
  // (out_a, out_b); permute out_a/out_b back to the physical positions.
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
  applied = fgf::transpose(applied, perm);
  const mptensor::Axes axes = fg_all_axes(rank);
  return fgf::trace(fgf::conj(psi), applied, axes, axes);
}

// ---- folded (plain) networks ----------------------------------------------

template <class tensor>
typename tensor::value_type fg_folded_norm(const fg_patch& p,
                                           const std::vector<tensor>& reduced) {
  std::vector<std::pair<tensor, std::vector<int>>> nodes;
  for (int s = 0; s < p.nsite(); ++s) {
    nodes.push_back({reduced[s], {s * 5 + 0, s * 5 + 1, s * 5 + 2, s * 5 + 3}});
  }
  std::vector<int> rest;
  return fg_scalar(fg_assemble(p, nodes, rest));
}

// Close a rank-6 pair blob on bond (a, b) with the exactly contracted
// environment of the remaining sites' build_reduced tensors.
template <class tensor>
typename tensor::value_type fg_blob_closure(const fg_patch& p,
                                            const std::vector<tensor>& reduced,
                                            int a, int b, char dir,
                                            const tensor& blob) {
  std::vector<std::pair<tensor, std::vector<int>>> nodes;
  for (int s = 0; s < p.nsite(); ++s) {
    if (s == a || s == b) {
      continue;
    }
    nodes.push_back({reduced[s], {s * 5 + 0, s * 5 + 1, s * 5 + 2, s * 5 + 3}});
  }
  const std::vector<int> blob_labels =
      dir == 'h' ? std::vector<int>{a * 5 + 0, a * 5 + 1, a * 5 + 3,
                                    b * 5 + 1, b * 5 + 2, b * 5 + 3}
                 : std::vector<int>{a * 5 + 0, a * 5 + 1, a * 5 + 2,
                                    b * 5 + 0, b * 5 + 2, b * 5 + 3};
  nodes.push_back({blob, blob_labels});
  std::vector<int> rest;
  return fg_scalar(fg_assemble(p, nodes, rest));
}

// ---- operators -------------------------------------------------------------

// Pairing c_A c_B (A = JW-earlier site) in the JW representation
// c_A = c x 1, c_B = P x c (P = diag(1, -1), the hopping convention):
// the single nonzero plain element is (in_A, in_B, out_A, out_B) =
// (1, 1, 0, 0) = -1, because c_B |11> = P|1> x c|1> = -|10> and then
// c_A |10> = |00>. The oracle's pairing(a, b) applies c_b then c_a with a
// the JW-earlier site; pairing(b, a) is the negative (anticommutation).
inline tenes::real_tensor fg_pairing_ca_cb_plain() {
  tenes::real_tensor op(mptensor::Shape(2, 2, 2, 2));
  op.set_value(mptensor::Index(1, 1, 0, 0), -1.0);
  return op;
}

inline tenes::complex_tensor fg_to_complex(const tenes::real_tensor& a) {
  tenes::complex_tensor t(a.shape());
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    const mptensor::Index idx = a.global_index(n);
    double v;
    a.get_value(idx, v);
    t.set_value(idx, std::complex<double>(v, 0.0));
  }
  return t;
}

inline tenes::real_tensor fg_op_plain_real(const std::string& kind, int d,
                                           int op_seed) {
  if (kind == "identity") {
    return ib_identity_plain(d);
  }
  if (kind == "hopping") {
    return ib_hop_plain(d);
  }
  if (kind == "nn") {
    return ib_nn_plain(d);
  }
  if (kind == "pairing") {
    return fg_pairing_ca_cb_plain();
  }
  if (kind == "random") {
    return ib_random_even_plain(d, op_seed);
  }
  throw std::runtime_error("fg_op_plain_real: unknown operator kind");
}

template <class tensor>
tensor fg_op_plain(const std::string& kind, int d, int op_seed);

template <>
tenes::real_tensor fg_op_plain<tenes::real_tensor>(const std::string& kind,
                                                   int d, int op_seed) {
  return fg_op_plain_real(kind, d, op_seed);
}

template <>
tenes::complex_tensor fg_op_plain<tenes::complex_tensor>(
    const std::string& kind, int d, int op_seed) {
  if (kind == "random") {
    // Genuinely complex parity-even elements from two real draws.
    const tenes::real_tensor re = ib_random_even_plain(d, op_seed);
    const tenes::real_tensor im = ib_random_even_plain(d, op_seed + 50);
    tenes::complex_tensor op(re.shape());
    for (std::size_t n = 0; n < re.local_size(); ++n) {
      const mptensor::Index idx = re.global_index(n);
      double vr, vi;
      re.get_value(idx, vr);
      im.get_value(idx, vi);
      op.set_value(idx, std::complex<double>(vr, vi));
    }
    return op;
  }
  return fg_to_complex(fg_op_plain_real(kind, d, op_seed));
}

template <class tensor>
double fg_max_abs_entry(const tensor& a) {
  double m = 0.0;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    typename tensor::value_type v;
    a.get_value(a.global_index(n), v);
    m = std::max(m, std::abs(v));
  }
  return m;
}

// ---- judgment helpers (contract formulas) ----------------------------------

// Scalar: |got - ref| <= 1e-10 * max(|ref_norm| * op_max, |ref|).
template <class V>
void fg_check_scalar(const std::string& label, V got, V ref, V ref_norm,
                     double op_max) {
  const double tol =
      1.0e-10 * std::max(std::abs(ref_norm) * op_max, std::abs(ref));
  INFO(label << ": got=" << got << " ref=" << ref
             << " |diff|=" << std::abs(got - ref) << " tol=" << tol);
  CHECK(std::abs(got - ref) <= tol);
}

// Signal floor: |ref| >= 1e-3 * |ref_norm| * op_max (seed quality gate).
template <class V>
void fg_require_signal(const std::string& label, V ref, V ref_norm,
                       double op_max) {
  INFO(label << ": ref=" << ref << " ref_norm=" << ref_norm
             << " op_max=" << op_max);
  REQUIRE(std::abs(ref) >= 1.0e-3 * std::abs(ref_norm) * op_max);
}

// Elementwise: |diff| <= 1e-12 * max(1, max |element| of either tensor).
template <class tensor>
void fg_check_allclose(const tensor& got, const tensor& want,
                       const std::string& label) {
  REQUIRE(got.shape() == want.shape());
  const double scale =
      std::max(1.0, std::max(fg_max_abs_entry(got), fg_max_abs_entry(want)));
  const double tol = 1.0e-12 * scale;
  double max_dev = 0.0;
  mptensor::Index worst;
  for (std::size_t n = 0; n < want.local_size(); ++n) {
    const mptensor::Index idx = want.global_index(n);
    typename tensor::value_type w, g;
    want.get_value(idx, w);
    got.get_value(idx, g);
    const double dev = std::abs(g - w);
    if (dev > max_dev) {
      max_dev = dev;
      worst = idx;
    }
  }
  INFO(label << ": max |got-want| = " << max_dev << " (tol " << tol
             << ") at index " << ib_index_string(worst));
  CHECK(max_dev <= tol);
}

// ---- T1-T4: full-patch folded norm vs graded norm --------------------------

template <class tensor>
void fg_run_norm_case(int lx, int ly, int d, const std::string& vps, int seed,
                      bool also_boson) {
  const std::string label = "shape=" + std::to_string(lx) + "x" +
                            std::to_string(ly) + " d=" + std::to_string(d) +
                            " vp=" + vps + " seed=" + std::to_string(seed);
  INFO(label);
  const fg_patch p{lx, ly};
  const fgf::parity_vector vp = fg_pv(vps);
  const fgf::parity_vector phys =
      also_boson ? fgf::parity_vector{false, false} : ib_phys_parity(d);
  std::vector<fg_ftensor<tensor>> sites;
  for (int s = 0; s < p.nsite(); ++s) {
    sites.push_back(fg_make_site<tensor>(p, s, vp, phys, seed));
  }
  std::vector<int> labels;
  const fg_ftensor<tensor> psi = fg_patch_state(p, sites, labels);
  const auto norm_graded = fg_graded_norm(psi);
  REQUIRE(std::abs(norm_graded) > 0.0);
  std::vector<tensor> reduced;
  for (int s = 0; s < p.nsite(); ++s) {
    reduced.push_back(fgf::build_reduced(sites[s]));
  }
  const auto norm_folded = fg_folded_norm(p, reduced);
  fg_check_scalar(label + " [folded norm vs graded norm]", norm_folded,
                  norm_graded, norm_graded, 1.0);
  if (also_boson) {
    // All-even control (T4): the plain (sign-free boson) contraction of the
    // same network must agree with both.
    std::vector<std::pair<tensor, std::vector<int>>> nodes;
    for (int s = 0; s < p.nsite(); ++s) {
      nodes.push_back(
          {sites[s].t,
           {s * 5 + 0, s * 5 + 1, s * 5 + 2, s * 5 + 3, s * 5 + 4}});
    }
    std::vector<int> rest;
    const tensor psi_plain = fg_assemble(p, nodes, rest);
    const mptensor::Axes axes = fg_all_axes(psi_plain.shape().size());
    const auto norm_boson = mptensor::trace(psi_plain, psi_plain, axes, axes);
    fg_check_scalar(label + " [boson norm vs graded norm]", norm_boson,
                    norm_graded, norm_graded, 1.0);
    fg_check_scalar(label + " [boson norm vs folded norm]", norm_boson,
                    norm_folded, norm_graded, 1.0);
  }
}

}  // namespace

TEST_CASE(
    "fold geometry T1: folded norm equals graded norm on open patches "
    "(d=2, D=2, shapes LXxLY = columns x rows)") {
  const std::pair<int, int> shapes[] = {{1, 3}, {2, 2}, {2, 3}, {3, 2}, {3, 3}};
  for (const auto& sh : shapes) {
    for (const int seed : {0, 173}) {
      fg_run_norm_case<tenes::real_tensor>(sh.first, sh.second, 2, "eo", seed,
                                           false);
    }
  }
}

TEST_CASE("fold geometry T2: folded norm equals graded norm at D=3 ledgers") {
  const std::pair<int, int> shapes[] = {{2, 3}, {3, 3}};
  for (const auto& sh : shapes) {
    for (const char* vps : {"eeo", "eoo"}) {
      for (const int seed : {0, 173}) {
        fg_run_norm_case<tenes::real_tensor>(sh.first, sh.second, 2, vps, seed,
                                             false);
      }
    }
  }
}

TEST_CASE("fold geometry T3: folded norm equals graded norm at d=4 (2x3)") {
  // d=4 [e,e,o,o] full-patch is excluded from CI by the contract's T3 row;
  // [e,e,o,o] appears only in the single-bond T7/T8 rows.
  for (const char* vps : {"eo", "eeo"}) {
    fg_run_norm_case<tenes::real_tensor>(2, 3, 4, vps, 0, false);
  }
}

TEST_CASE(
    "fold geometry T4: all-even ledger makes graded, folded and plain boson "
    "norms coincide (2x3, d=2)") {
  fg_run_norm_case<tenes::real_tensor>(2, 3, 2, "ee", 0, true);
}

// ---- T5: Fock-oracle anchors ------------------------------------------------

namespace {

struct fg_t5_bond {
  int a;
  int b;
  char dir;
  double hop;  // (one_body(a,b) + one_body(b,a)) / norm
};

struct fg_t5_anchor {
  int seed;
  double norm;
  fg_t5_bond bonds[7];
  // raw (unnormalized) <psi| c_a c_b |psi> values, a = JW-earlier site
  int pair_a[2];
  int pair_b[2];
  char pair_dir[2];
  double pair_raw[2];
  int npair;
};

// Generated from test/fermion/fock_oracle.py (run inside test/fermion/):
//   patch, tensors, lp = make_case(3, 2, [False, True], seed)
//   oracle = Oracle(patch, tensors, lp)
//   norm; (one_body(a,b)+one_body(b,a))/norm per internal bond;
//   oracle.pairing(a, b) raw.
// The two seed-0 pairing values match the contract's verified examples.
const fg_t5_anchor fg_t5_anchors[] = {
    {0,
     3.51258625440347983e-08,
     {{0, 1, 'h', 2.62440901358146905e-02},
      {0, 3, 'v', -6.76011456137395189e-02},
      {1, 2, 'h', -8.82886250302460790e-02},
      {1, 4, 'v', -1.18660107881122398e-01},
      {2, 5, 'v', 3.65591319316701135e-01},
      {3, 4, 'h', 4.31414486641046035e-01},
      {4, 5, 'h', 4.89815128479346995e-01}},
     {4, 1},
     {5, 4},
     {'h', 'v'},
     {2.97909261178536610e-09, 3.44039889340022481e-09},
     2},
    {173,
     5.30237345364552974e-09,
     {{0, 1, 'h', -1.52604887557228036e-01},
      {0, 3, 'v', 2.37418698623122273e-02},
      {1, 2, 'h', -7.44498160631679728e-02},
      {1, 4, 'v', 4.91466219855048680e-01},
      {2, 5, 'v', 1.80237359513616607e-01},
      {3, 4, 'h', -6.15215587899567359e-02},
      {4, 5, 'h', 1.51645226978584981e-02}},
     {4, 0},
     {5, 0},
     {'h', 0},
     {2.02334911785438392e-09, 0.0},
     1},
};

}  // namespace

TEST_CASE(
    "fold geometry T5: Fock oracle anchors on the 3x2 patch "
    "(3 columns x 2 rows, site = x + 3y)") {
  const fg_patch p{3, 2};
  const fgf::parity_vector vp = fg_pv("eo");
  const fgf::parity_vector phys = fg_pv("eo");
  for (const auto& anchor : fg_t5_anchors) {
    const std::string base = "T5 seed=" + std::to_string(anchor.seed);
    std::vector<fg_ftensor<tenes::real_tensor>> sites;
    for (int s = 0; s < p.nsite(); ++s) {
      sites.push_back(
          fg_make_site<tenes::real_tensor>(p, s, vp, phys, anchor.seed));
    }
    std::vector<int> labels;
    const auto psi = fg_patch_state(p, sites, labels);
    const double norm_graded = fg_graded_norm(psi);
    fg_check_scalar(base + " [graded norm vs oracle]", norm_graded, anchor.norm,
                    anchor.norm, 1.0);

    std::vector<tenes::real_tensor> reduced;
    for (int s = 0; s < p.nsite(); ++s) {
      reduced.push_back(fgf::build_reduced(sites[s]));
    }
    const double norm_folded = fg_folded_norm(p, reduced);
    fg_check_scalar(base + " [folded norm vs oracle]", norm_folded, anchor.norm,
                    anchor.norm, 1.0);

    const auto hop_gate = fgf::wrap_twosite_gate(ib_hop_plain(2), phys, phys);
    for (const auto& bond : anchor.bonds) {
      const double raw =
          fg_graded_twosite(psi, labels, bond.a, bond.b, hop_gate);
      fg_check_scalar(base + " [graded hop(" + std::to_string(bond.a) + "," +
                          std::to_string(bond.b) + ") vs oracle]",
                      raw / norm_graded, bond.hop, 1.0, 1.0);
    }

    const auto pair_gate =
        fgf::wrap_twosite_gate(fg_pairing_ca_cb_plain(), phys, phys);
    for (int k = 0; k < anchor.npair; ++k) {
      const double raw = fg_graded_twosite(psi, labels, anchor.pair_a[k],
                                           anchor.pair_b[k], pair_gate);
      fg_check_scalar(base + " [graded pairing(" +
                          std::to_string(anchor.pair_a[k]) + "," +
                          std::to_string(anchor.pair_b[k]) + ") raw vs oracle]",
                      raw, anchor.pair_raw[k], anchor.norm, 1.0);
    }
  }
}

// ---- T6: graded SVD reconstruction of gates --------------------------------

namespace {

template <class tensor>
void fg_run_t6_case(const std::string& label, int d, const tensor& plain) {
  INFO(label);
  const fgf::parity_vector phys = ib_phys_parity(d);
  const fg_ftensor<tensor> gate = fgf::wrap_twosite_gate(plain, phys, phys);
  fg_ftensor<tensor> u, vt;
  std::vector<double> s;
  REQUIRE(fgf::svd(gate, mptensor::Axes(0, 2), mptensor::Axes(1, 3), u, s,
                   vt) == 0);
  fg_ftensor<tensor> us = u;
  us.multiply_vector(s, 2);
  const fg_ftensor<tensor> recon =
      fgf::tensordot(us, vt, mptensor::Axes(2), mptensor::Axes(0));
  const fg_ftensor<tensor> want =
      fgf::transpose(gate, mptensor::Axes(0, 2, 1, 3));
  fg_check_allclose(recon.t, want.t, label);
}

}  // namespace

TEST_CASE(
    "fold geometry T6: U diag(s) VT reconstructs the graded transpose "
    "(0,2,1,3) of the gate") {
  fg_run_t6_case<tenes::real_tensor>("T6 hop d=2 real", 2, ib_hop_plain(2));
  fg_run_t6_case<tenes::real_tensor>("T6 hop d=4 real", 4, ib_hop_plain(4));
  fg_run_t6_case<tenes::real_tensor>("T6 random-even d=2 real op_seed=31", 2,
                                     ib_random_even_plain(2, 31));
  fg_run_t6_case<tenes::complex_tensor>("T6 hop d=2 complex", 2,
                                        fg_to_complex(ib_hop_plain(2)));
  fg_run_t6_case<tenes::complex_tensor>("T6 hop d=4 complex", 4,
                                        fg_to_complex(ib_hop_plain(4)));
  fg_run_t6_case<tenes::complex_tensor>(
      "T6 random-even d=2 complex op_seed=31", 2,
      fg_op_plain<tenes::complex_tensor>("random", 2, 31));
}

// ---- T7/T8: single-bond svd-split anchors ----------------------------------

namespace {

// One row of the T7/T8 matrix: freestanding two-site bond, every virtual
// leg carrying the ledger vps; site tensors are parity-even deterministic
// (seeds 11 and 12). The reference is apply_pair_op(build_pair_state(...))
// - existing primitives only, no blob API.
template <class tensor>
void fg_run_t78_case(const std::string& label, int d, const std::string& vps,
                     char dir, const std::string& op_kind, int op_seed) {
  INFO(label << " site_seeds=(11,12)");
  const fgf::parity_vector vp = fg_pv(vps);
  const fgf::parity_vector phys = ib_phys_parity(d);
  const fg_ftensor<tensor> TnA = fg_make_bond_site<tensor>(0, vp, phys, 11);
  const fg_ftensor<tensor> TnB = fg_make_bond_site<tensor>(1, vp, phys, 12);
  const tensor plain = fg_op_plain<tensor>(op_kind, d, op_seed);
  const fg_ftensor<tensor> gate = fgf::wrap_twosite_gate(plain, phys, phys);
  const auto dir_e = dir == 'h' ? fgf::reduced_pair_direction::horizontal
                                : fgf::reduced_pair_direction::vertical;
  const fg_ftensor<tensor> ref =
      fgf::apply_pair_op(fgf::build_pair_state(TnA, TnB, dir_e), gate);

  fg_ftensor<tensor> u, vt;
  std::vector<double> s;
  REQUIRE(fgf::svd(gate, mptensor::Axes(0, 2), mptensor::Axes(1, 3), u, s,
                   vt) == 0);
  fg_ftensor<tensor> us = u;
  us.multiply_vector(s, 2);
  // A' legs (l, t, r, b, out_A, k); B' legs (l, t, r, b, k, out_B).
  const fg_ftensor<tensor> Ap =
      fgf::tensordot(TnA, us, mptensor::Axes(4), mptensor::Axes(0));
  const fg_ftensor<tensor> Bp =
      fgf::tensordot(TnB, vt, mptensor::Axes(4), mptensor::Axes(1));

  // T7: graded contraction over (bond, k) as two separate axes.
  const fg_ftensor<tensor> got7 =
      dir == 'h'
          ? fgf::tensordot(Ap, Bp, mptensor::Axes(2, 5), mptensor::Axes(0, 4))
          : fgf::tensordot(Ap, Bp, mptensor::Axes(3, 5), mptensor::Axes(1, 4));
  fg_check_allclose(got7.t, ref.t, label + " [T7 split vs apply_pair_op]");

  // T8: bundled-k fuse (bond first, fused = bond + D_bond * k), crossing
  // mask (-1)^{p(bond) p(k)} on the A-side fused leg only, then a single
  // fused-axis graded contraction.
  const std::size_t D = vp.size();
  const std::size_t K = s.size();
  const fgf::parity_vector pk = u.parity[2];
  std::vector<double> mask(D * K, 1.0);
  for (std::size_t i = 0; i < mask.size(); ++i) {
    if (vp[i % D] && pk[i / D]) {
      mask[i] = -1.0;
    }
  }
  fg_ftensor<tensor> got8;
  if (dir == 'h') {
    fg_ftensor<tensor> Af =
        fgf::reshape(fgf::transpose(Ap, mptensor::Axes(0, 1, 2, 5, 3, 4)),
                     mptensor::Shape(D, D, D * K, D, phys.size()));
    Af.multiply_vector(mask, 2);
    const fg_ftensor<tensor> Bf =
        fgf::reshape(fgf::transpose(Bp, mptensor::Axes(0, 4, 1, 2, 3, 5)),
                     mptensor::Shape(D * K, D, D, D, phys.size()));
    got8 = fgf::tensordot(Af, Bf, mptensor::Axes(2), mptensor::Axes(0));
  } else {
    fg_ftensor<tensor> Af =
        fgf::reshape(fgf::transpose(Ap, mptensor::Axes(0, 1, 2, 3, 5, 4)),
                     mptensor::Shape(D, D, D, D * K, phys.size()));
    Af.multiply_vector(mask, 3);
    const fg_ftensor<tensor> Bf =
        fgf::reshape(fgf::transpose(Bp, mptensor::Axes(0, 1, 4, 2, 3, 5)),
                     mptensor::Shape(D, D * K, D, D, phys.size()));
    got8 = fgf::tensordot(Af, Bf, mptensor::Axes(3), mptensor::Axes(1));
  }
  fg_check_allclose(got8.t, ref.t, label + " [T8 fused vs apply_pair_op]");
}

struct fg_t78_row {
  int d;
  const char* vps;
  char dir;
  const char* op;
  bool complex_case;
};

// The ten mandatory rows (T7 and T8 share them; T8 must not be thinned).
const fg_t78_row fg_t78_rows[] = {
    {2, "eo", 'h', "hopping", false},   {2, "eo", 'v', "hopping", false},
    {2, "eeo", 'h', "random", false},   {2, "eo", 'h', "random", true},
    {4, "eo", 'h', "hopping", false},   {4, "eo", 'v', "hopping", false},
    {4, "eeo", 'h', "random", false},   {4, "eeoo", 'h', "hopping", false},
    {4, "eeoo", 'v', "hopping", false}, {4, "eo", 'v', "nn", false},
};

}  // namespace

TEST_CASE(
    "fold geometry T7/T8: single-bond svd-split and bundled-k fuse match "
    "apply_pair_op elementwise") {
  int row_id = 0;
  for (const auto& row : fg_t78_rows) {
    ++row_id;
    const std::string label =
        "T7/T8 row " + std::to_string(row_id) + ": d=" + std::to_string(row.d) +
        " vp=" + row.vps + " dir=" + row.dir + " op=" + row.op +
        (row.complex_case ? " complex" : " real") + " op_seed=31";
    if (row.complex_case) {
      fg_run_t78_case<tenes::complex_tensor>(label, row.d, row.vps, row.dir,
                                             row.op, 31);
    } else {
      fg_run_t78_case<tenes::real_tensor>(label, row.d, row.vps, row.dir,
                                          row.op, 31);
    }
  }
}

// ---- T9-T12: blob closure vs graded truth ----------------------------------

namespace {

// Deterministic bond choice per shape: kept cheap (a dimension-1 patch
// boundary leg on the pair keeps the rank-16 reference fuse small) while
// still placing at least one pair site with >= 3 nontrivial legs inside the
// blob for every blind-spot shape (3x3 even contains the 4-leg center).
std::pair<int, int> fg_pick_bond(const fg_patch& p, char dir) {
  if (dir == 'h') {
    const int a = p.site(0, p.ly / 2);
    return {a, a + 1};
  }
  const int a = p.site(p.lx / 2, 0);
  return {a, a + p.lx};
}

// One closure case: blob (direct and naive, loaded with wrap_twosite_gate)
// closed in the exactly contracted build_reduced environment, judged
// against the single-layer graded truth. For op == identity the T12 norm
// checks run as well.
template <class tensor>
void fg_run_closure_case(int lx, int ly, int d, const std::string& vps,
                         char dir, const std::string& op_kind, int seed,
                         int op_seed, bool source_second) {
  const std::string label =
      "shape=" + std::to_string(lx) + "x" + std::to_string(ly) +
      " d=" + std::to_string(d) + " vp=" + vps + " dir=" + dir +
      " op=" + op_kind + " seed=" + std::to_string(seed) +
      " op_seed=" + std::to_string(op_seed) +
      (source_second ? " source=second" : " source=first");
  INFO(label);
  const fg_patch p{lx, ly};
  const fgf::parity_vector vp = fg_pv(vps);
  const fgf::parity_vector phys = ib_phys_parity(d);
  std::vector<fg_ftensor<tensor>> sites;
  for (int s = 0; s < p.nsite(); ++s) {
    sites.push_back(fg_make_site<tensor>(p, s, vp, phys, seed));
  }
  std::vector<int> labels;
  const fg_ftensor<tensor> psi = fg_patch_state(p, sites, labels);
  const auto norm_raw = fg_graded_norm(psi);
  REQUIRE(std::abs(norm_raw) > 0.0);

  const std::pair<int, int> bond = fg_pick_bond(p, dir);
  const int a = bond.first;
  const int b = bond.second;
  const tensor plain = fg_op_plain<tensor>(op_kind, d, op_seed);
  const double op_max = fg_max_abs_entry(plain);
  // Production convention: site tensors always go in geometric order; a
  // source on the geometrically-second site is expressed on the gate only,
  // via the graded transpose (1,0,3,2) after wrap_twosite_gate.
  fg_ftensor<tensor> gate = fgf::wrap_twosite_gate(plain, phys, phys);
  if (source_second) {
    gate = fgf::transpose(gate, mptensor::Axes(1, 0, 3, 2));
  }
  const auto value_ref = fg_graded_twosite(psi, labels, a, b, gate);
  if (op_kind != "identity") {
    fg_require_signal(label, value_ref, norm_raw, op_max);
  }

  std::vector<tensor> reduced;
  for (int s = 0; s < p.nsite(); ++s) {
    reduced.push_back(fgf::build_reduced(sites[s]));
  }
  const auto dir_e = dir == 'h' ? fgf::reduced_pair_direction::horizontal
                                : fgf::reduced_pair_direction::vertical;
  const tensor blob_direct =
      fgf::build_reduced_pair_direct(sites[a], sites[b], gate, dir_e);
  const tensor blob_naive =
      fgf::build_reduced_pair_naive(sites[a], sites[b], gate, dir_e);
  fg_check_scalar(label + " [direct closure vs graded truth]",
                  fg_blob_closure(p, reduced, a, b, dir, blob_direct),
                  value_ref, norm_raw, op_max);
  fg_check_scalar(label + " [naive closure vs graded truth]",
                  fg_blob_closure(p, reduced, a, b, dir, blob_naive), value_ref,
                  norm_raw, op_max);

  if (op_kind == "identity") {
    // T12: norm path == value path. The identity-pair blob and the folded
    // norm are reference-side (expected GREEN); the id-op blobs above are
    // the value path under test.
    const auto folded_norm = fg_folded_norm(p, reduced);
    fg_check_scalar(label + " [T12 folded norm vs graded norm]", folded_norm,
                    norm_raw, norm_raw, 1.0);
    fg_check_scalar(label + " [T12 graded id value vs graded norm]", value_ref,
                    norm_raw, norm_raw, 1.0);
    const tensor id_pair =
        fgf::build_reduced_identity_pair(sites[a], sites[b], dir_e);
    fg_check_scalar(label + " [T12 identity-pair closure vs graded norm]",
                    fg_blob_closure(p, reduced, a, b, dir, id_pair), norm_raw,
                    norm_raw, 1.0);
  }
}

}  // namespace

TEST_CASE(
    "fold geometry T9: blob closure equals graded truth on open patches "
    "(d=2)") {
  const std::pair<int, int> shapes[] = {{2, 2}, {2, 3}, {3, 2}, {3, 3}};
  for (const auto& sh : shapes) {
    for (const char* vps : {"eo", "eoo"}) {
      for (const char dir : {'h', 'v'}) {
        for (const char* op : {"identity", "hopping", "random"}) {
          for (const int seed : {0, 173}) {
            fg_run_closure_case<tenes::real_tensor>(sh.first, sh.second, 2, vps,
                                                    dir, op, seed, 31, false);
          }
        }
      }
    }
  }
  // Source-second registration, one horizontal and one vertical case.
  fg_run_closure_case<tenes::real_tensor>(2, 3, 2, "eo", 'h', "random", 0, 32,
                                          true);
  fg_run_closure_case<tenes::real_tensor>(3, 2, 2, "eo", 'v', "random", 0, 32,
                                          true);
}

TEST_CASE("fold geometry T10: blob closure equals graded truth at d=4 (2x3)") {
  for (const char* vps : {"eo", "eeo"}) {
    for (const char dir : {'h', 'v'}) {
      for (const char* op : {"identity", "hopping"}) {
        fg_run_closure_case<tenes::real_tensor>(2, 3, 4, vps, dir, op, 0, 31,
                                                false);
      }
    }
  }
}

TEST_CASE(
    "fold geometry T11: complex blob closure equals graded truth "
    "(2x3, d=2, horizontal, random-even)") {
  fg_run_closure_case<tenes::complex_tensor>(2, 3, 2, "eo", 'h', "random", 0,
                                             31, false);
}

// ---- T13: direct == naive ---------------------------------------------------

namespace {

template <class tensor>
void fg_run_t13_case(int d, char dir, const std::string& op_kind, int op_seed,
                     const char* type_name) {
  const std::string label = std::string("T13 d=") + std::to_string(d) +
                            " dir=" + dir + " op=" + op_kind + " " + type_name +
                            " op_seed=" + std::to_string(op_seed) +
                            " site_seeds=(21,22) vp=eo(mixed)";
  INFO(label);
  const fgf::parity_vector vp = fg_pv("eo");
  const fgf::parity_vector phys = ib_phys_parity(d);
  const fg_ftensor<tensor> TnA = fg_make_bond_site<tensor>(0, vp, phys, 21);
  const fg_ftensor<tensor> TnB = fg_make_bond_site<tensor>(1, vp, phys, 22);
  const tensor plain = fg_op_plain<tensor>(op_kind, d, op_seed);
  // New loading convention: wrap_twosite_gate (input swap only).
  const fg_ftensor<tensor> gate = fgf::wrap_twosite_gate(plain, phys, phys);
  const auto dir_e = dir == 'h' ? fgf::reduced_pair_direction::horizontal
                                : fgf::reduced_pair_direction::vertical;
  const tensor blob_direct =
      fgf::build_reduced_pair_direct(TnA, TnB, gate, dir_e);
  const tensor blob_naive =
      fgf::build_reduced_pair_naive(TnA, TnB, gate, dir_e);
  fg_check_allclose(blob_direct, blob_naive, label + " [direct vs naive]");
}

}  // namespace

TEST_CASE(
    "fold geometry T13: build_reduced_pair_direct equals "
    "build_reduced_pair_naive elementwise (wrap_twosite_gate loading)") {
  for (const int d : {2, 4}) {
    const std::vector<const char*> ops =
        d == 2
            ? std::vector<const char*>{"identity", "hopping", "nn", "pairing",
                                       "random"}
            : std::vector<const char*>{"identity", "hopping", "nn", "random"};
    for (const char dir : {'h', 'v'}) {
      for (const char* op : ops) {
        fg_run_t13_case<tenes::real_tensor>(d, dir, op, 31, "real");
        fg_run_t13_case<tenes::complex_tensor>(d, dir, op, 33, "complex");
      }
    }
  }
}

TEST_CASE(
    "fold geometry T13b: representative 2x3 closure - direct, naive and "
    "graded truth coincide") {
  // One environment-closed representative (contract T13, last clause):
  // 2x3, d=2, [e,o], horizontal, hopping, seed 0.
  const std::string label = "T13b shape=2x3 d=2 vp=eo dir=h op=hopping seed=0";
  INFO(label);
  const fg_patch p{2, 3};
  const fgf::parity_vector vp = fg_pv("eo");
  const fgf::parity_vector phys = fg_pv("eo");
  std::vector<fg_ftensor<tenes::real_tensor>> sites;
  for (int s = 0; s < p.nsite(); ++s) {
    sites.push_back(fg_make_site<tenes::real_tensor>(p, s, vp, phys, 0));
  }
  std::vector<int> labels;
  const auto psi = fg_patch_state(p, sites, labels);
  const double norm_raw = fg_graded_norm(psi);
  const std::pair<int, int> bond = fg_pick_bond(p, 'h');
  const auto gate = fgf::wrap_twosite_gate(ib_hop_plain(2), phys, phys);
  const double value_ref =
      fg_graded_twosite(psi, labels, bond.first, bond.second, gate);
  std::vector<tenes::real_tensor> reduced;
  for (int s = 0; s < p.nsite(); ++s) {
    reduced.push_back(fgf::build_reduced(sites[s]));
  }
  const auto dir_e = fgf::reduced_pair_direction::horizontal;
  const double closure_direct =
      fg_blob_closure(p, reduced, bond.first, bond.second, 'h',
                      fgf::build_reduced_pair_direct(
                          sites[bond.first], sites[bond.second], gate, dir_e));
  const double closure_naive =
      fg_blob_closure(p, reduced, bond.first, bond.second, 'h',
                      fgf::build_reduced_pair_naive(
                          sites[bond.first], sites[bond.second], gate, dir_e));
  fg_check_scalar(label + " [direct closure vs naive closure]", closure_direct,
                  closure_naive, norm_raw, 1.0);
  fg_check_scalar(label + " [direct closure vs graded truth]", closure_direct,
                  value_ref, norm_raw, 1.0);
  fg_check_scalar(label + " [naive closure vs graded truth]", closure_naive,
                  value_ref, norm_raw, 1.0);
}

// ---- T14: unsupported-path guards ------------------------------------------

namespace {

constexpr const char* fg_guard_message =
    "fermion CTM measurement supports nearest-neighbor two-site observables "
    "only";

tenes::SquareLattice fg_guard_lattice() {
  tenes::SquareLattice lattice(2, 2);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    lattice.physical_dims[site] = 2;
    lattice.virtual_dims[site] = {2, 2, 2, 2};
    lattice.initial_dirs[site] = {0.0};
    lattice.noises[site] = 1.0;
  }
  return lattice;
}

tenes::itps::PEPS_Parameters fg_guard_params(bool mean_field,
                                             const std::string& outdir) {
  tenes::itps::PEPS_Parameters params;
  params.fermion = true;
  params.phys_parity.assign(4, tenes::fermion::parity_vector{false, true});
  params.MeanField_Env = mean_field;
  params.print_level = tenes::PrintLevel::none;
  params.outdir = outdir;
  params.CHI = 4;
  params.Max_CTM_Iteration = 2;
  params.CTM_Convergence_Epsilon = 1.0e-8;
  params.Use_RSVD = false;
  params.seed = 3;
  return params;
}

tenes::real_tensor fg_guard_hopping() {
  tenes::real_tensor hopping(mptensor::Shape(2, 2, 2, 2));
  hopping.set_value(mptensor::Index(0, 1, 1, 0), 1.0);
  hopping.set_value(mptensor::Index(1, 0, 0, 1), 1.0);
  return hopping;
}

// Inject a translation-invariant parity-even state and (for CTM) refresh
// the reduced-tensor environment; mirrors the wrapped-bond translation test.
void fg_guard_inject_state(tenes::itps::iTPS<tenes::real_tensor>& state,
                           bool mean_field) {
  auto& finfo = tenes::itps::iTPSTestAccessor::finfo(state);
  auto& Tn = tenes::itps::iTPSTestAccessor::Tn(state);
  tenes::real_tensor uniform(mptensor::Shape(2, 2, 2, 2, 2));
  const auto parity = tenes::fermion::Tn_parity(finfo, 0);
  for (std::size_t n = 0; n < uniform.local_size(); ++n) {
    const auto idx = uniform.global_index(n);
    if (tenes::fermion::count_odd(parity, idx) % 2 == 0) {
      uniform.set_value(idx, 0.125 * static_cast<double>(1 + n));
    }
  }
  for (auto& site_tensor : Tn) {
    site_tensor = uniform;
  }
  auto& lambda = tenes::itps::iTPSTestAccessor::lambda_tensor(state);
  for (auto& site_lambda : lambda) {
    for (auto& leg_lambda : site_lambda) {
      leg_lambda = {1.0, 0.75};
    }
  }
  if (!mean_field) {
    tenes::itps::iTPSTestAccessor::update_reduced_density_environment(state);
  }
}

}  // namespace

// T14a-T14c ran skipped (doctest::skip) while the guard was unimplemented:
// back then the measurement fell through to the bosonic contraction, whose
// tensors mismatch the fermionic reduced environment, and an mptensor
// assert SIGABRTed the whole test binary. The guard was implemented on
// 2026-08-28, so the cases now run unconditionally.
TEST_CASE(
    "fold geometry T14a: fermion CTM rejects a distance-2 two-site "
    "observable") {
  using tensor = tenes::real_tensor;
  const tenes::SquareLattice lattice = fg_guard_lattice();
  const auto params = fg_guard_params(false, "output_test_fold_geometry_t14a");
  tenes::Operators<tensor> twosite_ops;
  twosite_ops.emplace_back("hop_d2", 0, 0, 2, 0, fg_guard_hopping());
  tenes::itps::iTPS<tensor> state(
      MPI_COMM_WORLD, params, lattice, tenes::EvolutionOperators<tensor>{},
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      twosite_ops, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  fg_guard_inject_state(state, false);
  CHECK_THROWS_WITH(state.measure_twosite(),
                    doctest::Contains(fg_guard_message));
}

TEST_CASE("fold geometry T14b: fermion CTM rejects a multisite observable") {
  using tensor = tenes::real_tensor;
  const tenes::SquareLattice lattice = fg_guard_lattice();
  const auto params = fg_guard_params(false, "output_test_fold_geometry_t14b");
  // Three-site chain operator n0 n1 n2 as an explicit rank-6 tensor.
  tenes::real_tensor op3(mptensor::Shape(2, 2, 2, 2, 2, 2));
  for (int n0 = 0; n0 < 2; ++n0) {
    for (int n1 = 0; n1 < 2; ++n1) {
      for (int n2 = 0; n2 < 2; ++n2) {
        op3.set_value(mptensor::Index(n0, n1, n2, n0, n1, n2),
                      static_cast<double>(n0 * n1 * n2));
      }
    }
  }
  tenes::Operators<tensor> multisite_ops;
  multisite_ops.emplace_back("nnn", 0, 0, std::vector<int>{1, 2},
                             std::vector<int>{0, 0}, op3);
  tenes::itps::iTPS<tensor> state(
      MPI_COMM_WORLD, params, lattice, tenes::EvolutionOperators<tensor>{},
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, multisite_ops,
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  fg_guard_inject_state(state, false);
  CHECK_THROWS_WITH(state.measure_multisite(),
                    doctest::Contains(fg_guard_message));
}

TEST_CASE(
    "fold geometry T14c: fermion CTM rejects correlation measurement with "
    "r_max > 0") {
  using tensor = tenes::real_tensor;
  const tenes::SquareLattice lattice = fg_guard_lattice();
  const auto params = fg_guard_params(false, "output_test_fold_geometry_t14c");
  tenes::real_tensor number(mptensor::Shape(2, 2));
  number.set_value(mptensor::Index(1, 1), 1.0);
  tenes::Operators<tensor> onesite_ops;
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    onesite_ops.emplace_back("n", 0, site, number);
  }
  const tenes::itps::CorrelationParameter corparam(
      2, std::vector<std::tuple<int, int>>{{0, 0}});
  tenes::itps::iTPS<tensor> state(
      MPI_COMM_WORLD, params, lattice, tenes::EvolutionOperators<tensor>{},
      tenes::EvolutionOperators<tensor>{}, onesite_ops,
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{}, corparam,
      tenes::itps::TransferMatrix_Parameters{});
  fg_guard_inject_state(state, false);
  CHECK_THROWS_WITH(state.measure_correlation(),
                    doctest::Contains(fg_guard_message));
}

TEST_CASE(
    "fold geometry T14d: mean-field environment still measures "
    "nearest-neighbor two-site observables") {
  using tensor = tenes::real_tensor;
  const tenes::SquareLattice lattice = fg_guard_lattice();
  const auto params = fg_guard_params(true, "output_test_fold_geometry_t14d");
  tenes::Operators<tensor> twosite_ops;
  twosite_ops.emplace_back("hop_d2", 0, 0, 1, 0, fg_guard_hopping());
  twosite_ops.emplace_back("hop_d2", 0, 0, 0, 1, fg_guard_hopping());
  tenes::itps::iTPS<tensor> state(
      MPI_COMM_WORLD, params, lattice, tenes::EvolutionOperators<tensor>{},
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      twosite_ops, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  fg_guard_inject_state(state, true);
  std::vector<std::map<tenes::itps::Bond, double>> measured;
  CHECK_NOTHROW(measured = state.measure_twosite());
  REQUIRE(measured.size() >= 1);
  for (const auto& entry : measured[0]) {
    CHECK(std::isfinite(entry.second));
  }
}

TEST_CASE(
    "fold geometry T14e: fermion CTM rejects a same-site (0,0) two-site "
    "observable") {
  // Pins the dx = dy = 0 corner of the is_nearest_neighbor test (a two-site
  // operator on a single site is NOT a nearest-neighbor bond); the guard
  // implementation rejects it, but no case fixed that until the full-branch
  // review flagged the gap.
  using tensor = tenes::real_tensor;
  const tenes::SquareLattice lattice = fg_guard_lattice();
  const auto params = fg_guard_params(false, "output_test_fold_geometry_t14e");
  tenes::Operators<tensor> twosite_ops;
  twosite_ops.emplace_back("hop_d2", 0, 0, 0, 0, fg_guard_hopping());
  tenes::itps::iTPS<tensor> state(
      MPI_COMM_WORLD, params, lattice, tenes::EvolutionOperators<tensor>{},
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      twosite_ops, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});
  fg_guard_inject_state(state, false);
  CHECK_THROWS_WITH(state.measure_twosite(),
                    doctest::Contains(fg_guard_message));
}

// ---- addendum: doubled_pipeline with asymmetric bra/ket layers -------------
//
// The upcoming blob rewrite calls detail::doubled_pipeline(bra, ket) with
// bra != ket AND with a bond leg whose dimension and parity ledger differ
// between the two layers (the ket side carries the gate's bundled-k leg).
// This test pins that usage independently of the implementation: the
// reference below rebuilds the pipeline from existing primitives only
// (f::conj, f::tensordot, f::apply_swap, f::transpose, mptensor::reshape),
// with the six-pair joint mask written out as a local constant.

namespace {

template <class tensor>
fg_ftensor<tensor> fgdp_make_site(const fgf::leg_parities& lp, int site_id,
                                  int seed) {
  mptensor::Shape sh;
  for (const auto& leg : lp) {
    sh.push(leg.size());
  }
  tensor t(sh);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    if (fgf::count_odd(lp, idx) % 2 == 0) {
      fg_set_entry(t, idx, site_id, seed);
    }
  }
  return fg_ftensor<tensor>{t, lp};
}

template <class tensor>
void fgdp_run_asymmetric_case(const char* type_name) {
  const std::string label =
      std::string("doubled_pipeline asymmetric bra/ket ") + type_name +
      " site_seeds=(41,42)";
  INFO(label);
  // bra: every leg D=2 [e,o], physical d=2 [e,o].  ket: identical except the
  // RIGHT virtual leg, widened to D=4 [e,o,o,e] (dimension AND ledger differ
  // between the layers on that one leg).  Both parity-even, deterministic.
  const fgf::parity_vector eo = fg_pv("eo");
  const fgf::parity_vector eooe = fg_pv("eooe");
  const fgf::leg_parities bra_lp{eo, eo, eo, eo, eo};
  const fgf::leg_parities ket_lp{eo, eo, eooe, eo, eo};
  const fg_ftensor<tensor> bra = fgdp_make_site<tensor>(bra_lp, 0, 41);
  const fg_ftensor<tensor> ket = fgdp_make_site<tensor>(ket_lp, 1, 42);

  const tensor got = fgf::detail::doubled_pipeline(bra, ket);

  // Reference: conj(bra) (x) ket outer product -> joint swaps -> (ket, bra)
  // interleave as a GRADED transpose -> adjacent fuse (plain reshape).
  // Production's six-pair routing, written out independently of
  // reduced.hpp: all ordered label pairs (x, y) with y in {0=left,
  // 3=bottom}.  Each pair contributes a ket(x)-bra(y) and a bra(x)-bra(y)
  // elementwise swap sign; the swaps are diagonal, so their order is
  // irrelevant.
  const auto build_reference = [&](bool with_joint_swaps) {
    fg_ftensor<tensor> doubled =
        fgf::tensordot(fgf::conj(bra), ket, mptensor::Axes(), mptensor::Axes());
    // bra layer at axes 0..4, ket layer at axes 5..9; leg label == axis
    // within its layer.
    if (with_joint_swaps) {
      const int mask_pairs[6][2] = {{1, 0}, {2, 0}, {3, 0},
                                    {0, 3}, {1, 3}, {2, 3}};
      for (const auto& xy : mask_pairs) {
        fgf::apply_swap(doubled, 5 + xy[0], xy[1]);
        fgf::apply_swap(doubled, xy[0], xy[1]);
      }
    }
    mptensor::Axes interleaved;
    for (int ax = 0; ax < 4; ++ax) {
      interleaved.push(5 + ax);
      interleaved.push(ax);
    }
    interleaved.push(9);
    interleaved.push(4);
    const fg_ftensor<tensor> ordered = fgf::transpose(doubled, interleaved);
    mptensor::Shape sh;
    for (int ax = 0; ax < 4; ++ax) {
      sh.push(ordered.shape()[2 * ax] * ordered.shape()[2 * ax + 1]);
    }
    sh.push(ordered.shape()[8]);
    sh.push(ordered.shape()[9]);
    return mptensor::reshape(ordered.t, sh);
  };

  const tensor ref = build_reference(true);
  REQUIRE(got.shape() == ref.shape());
  REQUIRE(got.shape().size() == 6);
  REQUIRE(fg_max_abs_entry(ref) > 0.0);
  fg_check_allclose(got, ref, label);

  // Vacuity guard: on this input the joint swaps must actually change the
  // result, otherwise the comparison above would hold with no signs at all.
  const tensor unswapped = build_reference(false);
  double swap_effect = 0.0;
  for (std::size_t n = 0; n < ref.local_size(); ++n) {
    typename tensor::value_type a, b;
    const mptensor::Index idx = ref.global_index(n);
    ref.get_value(idx, a);
    unswapped.get_value(idx, b);
    swap_effect = std::max(swap_effect, std::abs(a - b));
  }
  REQUIRE(swap_effect > 0.0);
}

}  // namespace

TEST_CASE(
    "fold geometry DP: doubled_pipeline with asymmetric bra/ket layers "
    "matches the primitive reference") {
  fgdp_run_asymmetric_case<tenes::real_tensor>("real");
  fgdp_run_asymmetric_case<tenes::complex_tensor>("complex");
}

// ---- T16: doubled_pipeline_traced -------------------------------------------
//
// Contract: work/fermion/twosite-speedup/CONTRACT-task1.md.
//   C-1  detail::doubled_pipeline_traced(bra, ket) equals
//        contract(detail::doubled_pipeline(bra, ket), Axes(4), Axes(5))
//        elementwise, for parity-even rank-5 layers whose virtual legs may
//        differ in dimension and ledger between the two layers.
//   C-2  build_reduced(Tn) does not move against that same truth formula,
//        and build_reduced_op(Tn) keeps its rank-6 shape.
//   C-4  both layers are parity even on the way in (checked here, not a
//        behavioural requirement on the function).
// The truth side always goes through doubled_pipeline and closes the two
// physical legs here. It must never go through build_reduced: C-2 moves
// build_reduced onto the function under test, so using it as a reference
// would close the loop.
//
// Non-vacuity (contract section 3) is enforced per case rather than by
// construction: ledgers are asserted to carry odd indices (N-1) and to mix
// even with odd (N-2), the built layers are measured for weight in the
// odd-physical sector (N-3, both in the "is populated" and in the stronger
// "actually moves the truth" form), and the truth tensor is asserted
// nonzero (N-4). Complex cases additionally assert a nonzero imaginary
// part, so conjugating the bra layer is not the identity.

namespace {

// N-1: the ledger has at least one odd index.
bool fg_ledger_has_odd(const fgf::parity_vector& p) {
  for (const bool odd : p) {
    if (odd) {
      return true;
    }
  }
  return false;
}

// N-2: the ledger mixes even and odd indices.
bool fg_ledger_mixed(const fgf::parity_vector& p) {
  bool has_even = false;
  bool has_odd = false;
  for (const bool odd : p) {
    if (odd) {
      has_odd = true;
    } else {
      has_even = true;
    }
  }
  return has_even && has_odd;
}

// Largest |element| in the sector where the PHYSICAL index (axis 4) is odd.
// N-3: with no weight there, the physical ledger might as well be all even
// and any sign convention on the physical trace would pass.
template <class tensor>
double fg_max_abs_phys_odd(const fg_ftensor<tensor>& a) {
  const fgf::parity_vector& phys = a.parity[4];
  double m = 0.0;
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    const mptensor::Index idx = a.t.global_index(n);
    if (!phys[idx[4]]) {
      continue;
    }
    typename tensor::value_type v;
    a.t.get_value(idx, v);
    m = std::max(m, std::abs(v));
  }
  return m;
}

// The same layer with that sector deleted; still parity even, since
// dropping elements never creates an odd one. Used for the strong form of
// N-3: the truth must actually depend on the odd-physical sector.
template <class tensor>
fg_ftensor<tensor> fg_drop_phys_odd(const fg_ftensor<tensor>& a) {
  fg_ftensor<tensor> r = a;
  const fgf::parity_vector& phys = r.parity[4];
  for (std::size_t n = 0; n < r.t.local_size(); ++n) {
    const mptensor::Index idx = r.t.global_index(n);
    if (phys[idx[4]]) {
      r.t.set_value(idx, typename tensor::value_type(0));
    }
  }
  return r;
}

inline double fg_abs_imag(double) { return 0.0; }
inline double fg_abs_imag(const std::complex<double>& v) {
  return std::abs(v.imag());
}

template <class tensor>
double fg_max_abs_imag(const tensor& a) {
  double m = 0.0;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    typename tensor::value_type v;
    a.get_value(a.global_index(n), v);
    m = std::max(m, fg_abs_imag(v));
  }
  return m;
}

template <class tensor>
double fg_max_abs_diff(const tensor& a, const tensor& b) {
  REQUIRE(a.shape() == b.shape());
  double m = 0.0;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    const mptensor::Index idx = a.global_index(n);
    typename tensor::value_type va, vb;
    a.get_value(idx, va);
    b.get_value(idx, vb);
    m = std::max(m, std::abs(va - vb));
  }
  return m;
}

// The C-1 / C-2 truth: the existing rank-6 pipeline with its two physical
// legs closed against each other.
template <class tensor>
tensor fg_close_phys(const tensor& open) {
  REQUIRE(open.shape().size() == 6);
  return mptensor::contract(open, mptensor::Axes(4), mptensor::Axes(5));
}

template <class tensor>
tensor fg_pipeline_closed(const fg_ftensor<tensor>& bra,
                          const fg_ftensor<tensor>& ket) {
  return fg_close_phys(fgf::detail::doubled_pipeline(bra, ket));
}

// Mutation defense for N-1: how much the closed tensor moves if the sign
// (-1)^{|s|} of the odd physical indices is dropped on the way in. A case
// where this is 0 is satisfied by any physical sign convention and can
// never fail on a sign mistake, whatever the implementation does.
template <class tensor>
double fg_phys_sign_sensitivity(const tensor& open,
                                const fgf::parity_vector& phys) {
  REQUIRE(open.shape()[4] == phys.size());
  tensor flipped = open;
  std::vector<double> sign(phys.size());
  for (std::size_t i = 0; i < phys.size(); ++i) {
    sign[i] = phys[i] ? -1.0 : 1.0;
  }
  flipped.multiply_vector(sign, 4);
  return fg_max_abs_diff(fg_close_phys(open), fg_close_phys(flipped));
}

struct fg_t16_layers {
  std::string name;
  fgf::leg_parities bra;
  fgf::leg_parities ket;
  //! ket is literally the bra tensor (symmetric double layer).
  bool same;
};

std::vector<fg_t16_layers> fg_t16_layer_cases() {
  const fgf::parity_vector eo = fg_pv("eo");
  const fgf::parity_vector eoo = fg_pv("eoo");
  const fgf::parity_vector eooe = fg_pv("eooe");
  const fgf::parity_vector d2 = ib_phys_parity(2);  // [e,o]
  const fgf::parity_vector d4 = ib_phys_parity(4);  // [e,o,o,e]
  std::vector<fg_t16_layers> cases;
  // bra == ket, the build_reduced shape.
  cases.push_back(
      {"sym d=2 D=2(eo)", {eo, eo, eo, eo, d2}, {eo, eo, eo, eo, d2}, true});
  cases.push_back({"sym d=4 D=3(eoo)",
                   {eoo, eoo, eoo, eoo, d4},
                   {eoo, eoo, eoo, eoo, d4},
                   true});
  // bra != ket with matching ledgers (different site id and seed).
  cases.push_back(
      {"asym d=2 D=2(eo)", {eo, eo, eo, eo, d2}, {eo, eo, eo, eo, d2}, false});
  // bra != ket with per-leg ledgers that differ between the layers.
  cases.push_back({"asym d=4 per-leg mixed ledgers",
                   {eo, eoo, eo, eoo, d4},
                   {eoo, eo, eoo, eo, d4},
                   false});
  // C-3: one ket leg wider than bra's AND on a different ledger (the gate's
  // bundled-k channel fused into a bond leg), horizontal bond.
  cases.push_back({"asym d=2 bundled ket right D=4(eooe)",
                   {eo, eo, eo, eo, d2},
                   {eo, eo, eooe, eo, d2},
                   false});
  // the same on a vertical bond, at d = 4.
  cases.push_back({"asym d=4 bundled ket bottom D=4(eooe)",
                   {eo, eo, eo, eo, d4},
                   {eo, eo, eo, eooe, d4},
                   false});
  // Mirror: the wide leg on the BRA side. C-1 lets either layer carry it
  // ("the virtual dimensions and ledgers may differ from each other"), so a
  // pipeline that quietly assumed ket >= bra would be a contract violation.
  cases.push_back({"asym d=2 bundled bra top D=4(eooe)",
                   {eo, eooe, eo, eo, d2},
                   {eo, eo, eo, eo, d2},
                   false});
  return cases;
}

template <class tensor>
void fg_run_t16_case(const fg_t16_layers& layers, int bra_seed, int ket_seed,
                     bool complex_scalars, const char* type_name) {
  const std::string seeds_desc =
      layers.same ? "seed=" + std::to_string(bra_seed) + " (ket == bra)"
                  : "seeds=(" + std::to_string(bra_seed) + "," +
                        std::to_string(ket_seed) + ")";
  const std::string label =
      std::string("T16 ") + layers.name + " " + type_name + " " + seeds_desc;
  INFO(label);

  // Ledger preconditions: rank 5 layers, matching physical ledgers (C-3),
  // an odd physical index (N-1), mixed virtual ledgers (N-2).
  REQUIRE(layers.bra.size() == 5);
  REQUIRE(layers.ket.size() == 5);
  REQUIRE(layers.bra[4] == layers.ket[4]);
  REQUIRE(fg_ledger_has_odd(layers.bra[4]));
  for (int ax = 0; ax < 4; ++ax) {
    REQUIRE(fg_ledger_mixed(layers.bra[ax]));
    REQUIRE(fg_ledger_mixed(layers.ket[ax]));
  }
  if (layers.same) {
    REQUIRE(layers.bra == layers.ket);
  }

  const fg_ftensor<tensor> bra =
      fgdp_make_site<tensor>(layers.bra, 0, bra_seed);
  const fg_ftensor<tensor> ket =
      layers.same ? bra : fgdp_make_site<tensor>(layers.ket, 1, ket_seed);

  // C-4: the inputs satisfy the precondition the function assumes.
  REQUIRE(fgf::parity_violation(bra) == 0.0);
  REQUIRE(fgf::parity_violation(ket) == 0.0);
  // N-3, populated form.
  REQUIRE(fg_max_abs_phys_odd(bra) > 0.0);
  REQUIRE(fg_max_abs_phys_odd(ket) > 0.0);
  if (complex_scalars) {
    REQUIRE(fg_max_abs_imag(bra.t) > 0.0);
    REQUIRE(fg_max_abs_imag(ket.t) > 0.0);
  }

  const tensor open = fgf::detail::doubled_pipeline(bra, ket);
  const tensor want = fg_close_phys(open);
  REQUIRE(want.shape().size() == 4);
  // N-4.
  REQUIRE(fg_max_abs_entry(want) > 0.0);
  // N-3, strong form: deleting the odd-physical sector must move the truth,
  // otherwise nothing in this case can see the sign under test.
  const tensor want_even_phys =
      fg_pipeline_closed(fg_drop_phys_odd(bra), fg_drop_phys_odd(ket));
  REQUIRE(fg_max_abs_diff(want, want_even_phys) > 0.0);
  // N-1, mutation defense: this case can tell a dropped physical sign apart.
  REQUIRE(fg_phys_sign_sensitivity(open, layers.bra[4]) > 0.0);

  const tensor got = fgf::detail::doubled_pipeline_traced(bra, ket);
  REQUIRE(got.shape().size() == 4);
  fg_check_allclose(got, want, label + " [traced vs closed doubled_pipeline]");
}

template <class tensor>
void fg_run_t16_reduced_case(const std::string& name,
                             const fgf::leg_parities& lp, int seed,
                             bool complex_scalars, const char* type_name) {
  const std::string label = std::string("T16b ") + name + " " + type_name +
                            " seed=" + std::to_string(seed);
  INFO(label);
  REQUIRE(lp.size() == 5);
  REQUIRE(fg_ledger_has_odd(lp[4]));
  for (int ax = 0; ax < 4; ++ax) {
    REQUIRE(fg_ledger_mixed(lp[ax]));
  }

  const fg_ftensor<tensor> Tn = fgdp_make_site<tensor>(lp, 0, seed);
  REQUIRE(fgf::parity_violation(Tn) == 0.0);
  REQUIRE(fg_max_abs_phys_odd(Tn) > 0.0);
  if (complex_scalars) {
    REQUIRE(fg_max_abs_imag(Tn.t) > 0.0);
  }

  // build_reduced_op is not part of this task: it must still be the rank-6
  // pipeline, with the four fused virtual legs and the two physical legs
  // open, element for element.
  const tensor open = fgf::detail::doubled_pipeline(Tn, Tn);
  const tensor op = fgf::build_reduced_op(Tn);
  REQUIRE(op.shape().size() == 6);
  mptensor::Shape expected;
  for (int ax = 0; ax < 4; ++ax) {
    expected.push(lp[ax].size() * lp[ax].size());
  }
  expected.push(lp[4].size());
  expected.push(lp[4].size());
  CHECK(op.shape() == expected);
  fg_check_allclose(op, open, label + " [build_reduced_op vs pipeline]");

  // build_reduced must not move by a single element, judged against the
  // C-1 truth formula written out here rather than taken from the function.
  const tensor want = fg_close_phys(open);
  REQUIRE(fg_max_abs_entry(want) > 0.0);
  const tensor want_even_phys = fg_pipeline_closed(
      fg_drop_phys_odd(Tn), fg_drop_phys_odd(Tn));  // N-3, strong form
  REQUIRE(fg_max_abs_diff(want, want_even_phys) > 0.0);
  // N-1, mutation defense (see fg_phys_sign_sensitivity).
  REQUIRE(fg_phys_sign_sensitivity(open, lp[4]) > 0.0);
  const tensor got = fgf::build_reduced(Tn);
  REQUIRE(got.shape().size() == 4);
  fg_check_allclose(got, want,
                    label + " [build_reduced vs closed doubled_pipeline]");
}

}  // namespace

TEST_CASE(
    "fold geometry T16: doubled_pipeline_traced equals doubled_pipeline with "
    "the two physical legs closed") {
  const std::pair<int, int> seed_pairs[] = {{51, 52}, {67, 71}};
  for (const auto& layers : fg_t16_layer_cases()) {
    for (const auto& seeds : seed_pairs) {
      fg_run_t16_case<tenes::real_tensor>(layers, seeds.first, seeds.second,
                                          false, "real");
      fg_run_t16_case<tenes::complex_tensor>(layers, seeds.first, seeds.second,
                                             true, "complex");
    }
  }
}

TEST_CASE(
    "fold geometry T16b: build_reduced still equals the closed "
    "doubled_pipeline and build_reduced_op keeps its rank-6 shape") {
  const fgf::parity_vector eo = fg_pv("eo");
  const fgf::parity_vector eoo = fg_pv("eoo");
  const fgf::parity_vector d2 = ib_phys_parity(2);
  const fgf::parity_vector d4 = ib_phys_parity(4);
  const std::vector<std::pair<std::string, fgf::leg_parities>> cases = {
      {"d=2 D=2(eo)", {eo, eo, eo, eo, d2}},
      {"d=4 D=3(eoo)", {eoo, eoo, eoo, eoo, d4}},
      {"d=4 per-leg mixed ledgers", {eo, eoo, eo, eoo, d4}},
  };
  for (const auto& c : cases) {
    for (const int seed : {51, 67}) {
      fg_run_t16_reduced_case<tenes::real_tensor>(c.first, c.second, seed,
                                                  false, "real");
      fg_run_t16_reduced_case<tenes::complex_tensor>(c.first, c.second, seed,
                                                     true, "complex");
    }
  }
}

// ---- T16c: the debug parity precondition of doubled_pipeline_traced --------
//
// The derivation of the physical-traced pipeline (DESIGN section 4.4, R-1)
// rests on both layers being parity even, and the function checks that
// premise in debug builds. T16 and T16b only ever hand it clean layers, so
// deleting the whole #ifndef NDEBUG block leaves all of them green: the
// check is production behaviour with nothing guarding it.
//
// One assertion per mutation this defends against:
//   - a contaminated bra layer is rejected AND named as the bra,
//   - a contaminated ket layer AND named as the ket, so dropping either
//     half of the block, or pasting one message over the other, is caught;
//   - the tolerance is RELATIVE to the layer's magnitude: an odd element
//     three orders above the bare 1e-10 constant but three below
//     1e-10*max|.| must pass, so dropping the max_abs() scaling is caught;
//   - the clean pair does not throw, so the case cannot pass by the guard
//     firing unconditionally.

namespace {

inline void fg16c_set(tenes::real_tensor& t, const mptensor::Index& idx,
                      double v) {
  t.set_value(idx, v);
}

inline void fg16c_set(tenes::complex_tensor& t, const mptensor::Index& idx,
                      double v) {
  t.set_value(idx, std::complex<double>(v, 0.0));
}

// Put `value` on a parity-odd index, leaving the ledger alone: the tensor
// stops being parity even by exactly that amount. The index is passed in
// rather than searched for, so every rank contaminates the same element.
template <class tensor>
fg_ftensor<tensor> fg16c_contaminate(const fg_ftensor<tensor>& a,
                                     const mptensor::Index& idx, double value) {
  REQUIRE(fgf::count_odd(a.parity, idx) % 2 == 1);
  fg_ftensor<tensor> bad = a;
  fg16c_set(bad.t, idx, value);
  REQUIRE(fgf::parity_violation(bad) == doctest::Approx(value));
  return bad;
}

template <class tensor>
void fg_run_t16c_case(const char* type_name) {
  const std::string label =
      std::string("T16c doubled_pipeline_traced parity guard ") + type_name;
  INFO(label);

  const fgf::parity_vector eo = fg_pv("eo");
  const fgf::leg_parities lp = {eo, eo, eo, eo, ib_phys_parity(2)};
  const fg_ftensor<tensor> bra = fgdp_make_site<tensor>(lp, 0, 71);
  const fg_ftensor<tensor> ket = fgdp_make_site<tensor>(lp, 1, 72);
  REQUIRE(fgf::parity_violation(bra) == 0.0);
  REQUIRE(fgf::parity_violation(ket) == 0.0);

  // Exactly one odd leg is selected, and the entry is one the generator
  // left at zero, so the injected value IS the resulting violation.
  const mptensor::Index odd_idx = mptensor::Index(1, 0, 0, 0, 0);
  REQUIRE(fgf::count_odd(lp, odd_idx) % 2 == 1);

  CHECK_NOTHROW(fgf::detail::doubled_pipeline_traced(bra, ket));

  const double gross = 1.0;
  const fg_ftensor<tensor> bad_bra = fg16c_contaminate(bra, odd_idx, gross);
  const fg_ftensor<tensor> bad_ket = fg16c_contaminate(ket, odd_idx, gross);

#ifndef NDEBUG
  CHECK_THROWS_WITH_AS(fgf::detail::doubled_pipeline_traced(bad_bra, ket),
                       doctest::Contains("bra layer is not parity even"),
                       std::runtime_error);
  CHECK_THROWS_WITH_AS(fgf::detail::doubled_pipeline_traced(bra, bad_ket),
                       doctest::Contains("ket layer is not parity even"),
                       std::runtime_error);
#else
  CHECK_NOTHROW(fgf::detail::doubled_pipeline_traced(bad_bra, ket));
  CHECK_NOTHROW(fgf::detail::doubled_pipeline_traced(bra, bad_ket));
#endif

  // Scale one layer up so that "1e-10 relative" and "1e-10 absolute" are
  // six orders apart, and contaminate strictly between them.
  fg_ftensor<tensor> big = bra;
  const std::vector<double> scale(lp[0].size(), 1.0e6);
  big.multiply_vector(scale, 0);
  const double tol = 1.0e-10 * std::max(1.0, fgf::max_abs(big));
  const double subtle = 1.0e-3 * tol;
  REQUIRE(subtle > 1.0e-10);
  REQUIRE(subtle < tol);
  const fg_ftensor<tensor> subtle_bra = fg16c_contaminate(big, odd_idx, subtle);
  CHECK_NOTHROW(fgf::detail::doubled_pipeline_traced(subtle_bra, ket));
}

}  // namespace

TEST_CASE(
    "fold geometry T16c: doubled_pipeline_traced checks its parity premise "
    "in debug builds") {
  fg_run_t16c_case<tenes::real_tensor>("real");
  fg_run_t16c_case<tenes::complex_tensor>("complex");
}

// ---- T15: absorbing the pair halves into the CTM environment ----------------
//
// Contract: work/fermion/twosite-speedup/CONTRACT-task3.md.
//   C-1  contract_reduced_pair_halves_density_CTM(C.., eT.., halves) equals
//        the existing blob closure
//        contract_reduced_pair_{horizontal,vertical}_density_CTM(C.., eT..,
//        blob) for the same sites, gate and environment (relative 1e-12).
//   C-2  C-1 for the operator halves (build_reduced_pair_halves) AND for the
//        norm halves (build_reduced_identity_halves).
//   C-3  build_reduced_pair_direct() and build_reduced_identity_pair() still
//        return, elementwise, tensordot(PA, PB, axis_a(), axis_b()) - the
//        truth side is written out here, not taken from the new API - and the
//        operator blob still equals build_reduced_pair_naive().
//   C-4  axis_a() / axis_b() depend on the direction only: they agree between
//        the operator and the norm halves, and across every row of the matrix
//        that shares a direction. The concrete values are an internal
//        convention and are deliberately NOT pinned.
//
// SCOPE LIMIT (contract section 3). The two sides of C-1 share the SAME
// halves object, so C-1 constrains only the CLOSURE - how the environment is
// absorbed around the two halves. It is blind to how the halves themselves
// are built: fold routing, crossing signs, leg order. Any error common to
// both halves moves both sides of C-1 by the same amount and leaves it
// green. What grounds the halves in physical truth is upstream and already
// exists: T9-T13 run the blob path against the single-layer graded
// contraction and the Fock oracle. T15 passing does NOT mean the two-site
// measurement is correct; it means the absorbing closure agrees with the
// blob closure it replaces.
//
// This block is appended after T16/T16b rather than sitting between T13b and
// T14 because it reuses the helpers introduced with T16 (fg_max_abs_diff,
// fg_ledger_mixed, fg_ledger_has_odd).

namespace {

// General, structureless environment entries (contract N-1): every element of
// every corner and edge comes from its own seed, so no two environment
// tensors are equal and none is symmetric under exchanging its two chi legs.
// ib_random_value() takes the rank explicitly, so it serves the rank-2
// corners and the rank-3 edges alike.
inline void fg15_set_env(tenes::real_tensor& t, const mptensor::Index& idx,
                         int seed, std::size_t rank) {
  t.set_value(idx, ib_random_value(seed, idx, rank));
}

inline void fg15_set_env(tenes::complex_tensor& t, const mptensor::Index& idx,
                         int seed, std::size_t rank) {
  t.set_value(idx,
              std::complex<double>(ib_random_value(seed, idx, rank),
                                   ib_random_value(seed + 500, idx, rank)));
}

template <class tensor>
tensor fg15_make_env_tensor(const mptensor::Shape& sh, int seed) {
  tensor t(sh);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    fg15_set_env(t, idx, seed, sh.size());
  }
  return t;
}

// The ten CTM tensors of a two-site window, named by the argument position
// they occupy in the closure calls (corners rank 2 (chi, chi), edges rank 3
// (chi, chi, D^2)).
template <class tensor>
struct fg15_env {
  tensor C1, C2, C3, C4;
  tensor eT1, eT2, eT3, eT4, eT5, eT6;
};

template <class tensor>
fg15_env<tensor> fg15_make_environment(std::size_t chi, std::size_t dd,
                                       int seed) {
  const mptensor::Shape corner(chi, chi);
  const mptensor::Shape edge(chi, chi, dd);
  fg15_env<tensor> e;
  e.C1 = fg15_make_env_tensor<tensor>(corner, seed + 1);
  e.C2 = fg15_make_env_tensor<tensor>(corner, seed + 2);
  e.C3 = fg15_make_env_tensor<tensor>(corner, seed + 3);
  e.C4 = fg15_make_env_tensor<tensor>(corner, seed + 4);
  e.eT1 = fg15_make_env_tensor<tensor>(edge, seed + 5);
  e.eT2 = fg15_make_env_tensor<tensor>(edge, seed + 6);
  e.eT3 = fg15_make_env_tensor<tensor>(edge, seed + 7);
  e.eT4 = fg15_make_env_tensor<tensor>(edge, seed + 8);
  e.eT5 = fg15_make_env_tensor<tensor>(edge, seed + 9);
  e.eT6 = fg15_make_env_tensor<tensor>(edge, seed + 10);
  return e;
}

// N-1: the environment must be general. Two equal environment tensors, or one
// that is symmetric under exchanging its two chi legs, would make some
// swapped wiring produce the very same number, and no test built on this
// environment could see that mistake.
template <class tensor>
void fg15_require_general_environment(const fg15_env<tensor>& e,
                                      const std::string& label) {
  INFO(label << " [N-1 general environment]");
  const tensor* corners[4] = {&e.C1, &e.C2, &e.C3, &e.C4};
  const tensor* edges[6] = {&e.eT1, &e.eT2, &e.eT3, &e.eT4, &e.eT5, &e.eT6};
  for (int i = 0; i < 4; ++i) {
    REQUIRE(fg_max_abs_entry(*corners[i]) > 0.0);
    REQUIRE(fg_max_abs_diff(
                *corners[i],
                mptensor::transpose(*corners[i], mptensor::Axes(1, 0))) > 0.0);
    for (int j = i + 1; j < 4; ++j) {
      REQUIRE(fg_max_abs_diff(*corners[i], *corners[j]) > 0.0);
    }
  }
  for (int i = 0; i < 6; ++i) {
    REQUIRE(fg_max_abs_entry(*edges[i]) > 0.0);
    REQUIRE(fg_max_abs_diff(
                *edges[i],
                mptensor::transpose(*edges[i], mptensor::Axes(1, 0, 2))) > 0.0);
    for (int j = i + 1; j < 6; ++j) {
      REQUIRE(fg_max_abs_diff(*edges[i], *edges[j]) > 0.0);
    }
  }
}

// The reference side: the existing blob closure, dispatched on direction.
template <class tensor>
typename tensor::value_type fg15_close_blob(const fg15_env<tensor>& e, char dir,
                                            const tensor& blob) {
  if (dir == 'h') {
    return fgf::contract_reduced_pair_horizontal_density_CTM(
        e.C1, e.C2, e.C3, e.C4, e.eT1, e.eT2, e.eT3, e.eT4, e.eT5, e.eT6, blob);
  }
  return fgf::contract_reduced_pair_vertical_density_CTM(
      e.C1, e.C2, e.C3, e.C4, e.eT1, e.eT2, e.eT3, e.eT4, e.eT5, e.eT6, blob);
}

// A copy of the environment with two slots exchanged; every corner has the
// same shape as every other corner, and likewise for the edges, so any of
// these swaps is a legal - and wrong - wiring of the same window.
template <class tensor>
fg15_env<tensor> fg15_swapped(const fg15_env<tensor>& e, int which) {
  fg15_env<tensor> p = e;
  switch (which) {
    case 0:
      std::swap(p.C1, p.C2);
      break;
    case 1:
      std::swap(p.eT1, p.eT2);
      break;
    case 2:
      std::swap(p.eT4, p.eT5);
      break;
    default:
      throw std::runtime_error("fg15_swapped: unknown swap");
  }
  return p;
}

const char* fg15_swap_name(int which) {
  switch (which) {
    case 0:
      return "C1<->C2";
    case 1:
      return "eT1<->eT2";
    case 2:
      return "eT4<->eT5";
    default:
      return "?";
  }
}

// N-2, wiring sensitivity. Exchange two environment slots and require the
// REFERENCE value to move by a relative amount far above round-off. Only the
// reference side is perturbed, so nothing here assumes how the absorbing
// closure is implemented. On a case where this fails, the environment cannot
// tell any wiring apart and C-1 would hold for a mis-wired implementation
// too - the case would be vacuous.
template <class tensor>
void fg15_require_wiring_sensitivity(const fg15_env<tensor>& e, char dir,
                                     const tensor& blob,
                                     typename tensor::value_type ref,
                                     const std::string& label) {
  const double floor = 1.0e-8 * std::abs(ref);
  for (int which = 0; which < 3; ++which) {
    const auto moved = fg15_close_blob(fg15_swapped(e, which), dir, blob);
    INFO(label << " [N-2 swap " << std::string(fg15_swap_name(which))
               << "]: ref=" << ref << " swapped=" << moved
               << " |diff|=" << std::abs(moved - ref) << " floor=" << floor);
    REQUIRE(std::abs(moved - ref) > floor);
  }
}

// C-1 judgment: relative 1e-12 against the blob-closure reference, anchored
// on the SCALE of the closure rather than on the compared value alone.
//
// Round-off in a contraction is set by the size of its largest intermediate,
// not by how much the final sum happens to cancel. C-1 and C-2 of one row
// close the same environment against blobs of the same shape, and their
// |got - ref| comes out bit-identical - while |ref| itself differs by up to
// 4700x between them. Measuring the smaller closure against itself therefore
// asks for far more than the arithmetic can deliver: it left 1.5x of head
// room and turned a 2-ulp difference into a compiler-dependent failure (CI,
// 2026-09-05). `scale` is the largest REFERENCE closure of the row, so an
// implementation that returns a huge `got` cannot inflate its own tolerance;
// and because it never exceeds the larger closure, the judgment on that
// larger closure is numerically unchanged. This is the design
// fg_check_scalar() (see above) already uses.
template <class V>
void fg15_check_rel(const std::string& label, V got, V ref, double scale) {
  const double tol =
      1.0e-12 * std::max(std::max(std::abs(got), std::abs(ref)), scale);
  INFO(label << ": got=" << got << " ref=" << ref << " |diff|="
             << std::abs(got - ref) << " tol=" << tol << " scale=" << scale);
  CHECK(std::abs(got - ref) <= tol);
}

// Number of singular values the GRADED svd of the loaded gate carries above
// the noise floor. This is the number of operator channels the gate has to
// bundle into the shared bond, computed here with the same public primitive
// the contract describes the construction in terms of (fgf::svd, as used by
// T7/T8) - never from the halves under test. A rank-1 gate (nn = n_A x n_B
// at d = 2, say) has exactly one channel and can legitimately bundle to
// exactly D^2; see the N-4 note in fg15_run_case().
template <class tensor>
std::size_t fg15_graded_gate_rank(const fg_ftensor<tensor>& gate) {
  fg_ftensor<tensor> u, vt;
  std::vector<double> s;
  REQUIRE(fgf::svd(gate, mptensor::Axes(0, 2), mptensor::Axes(1, 3), u, s,
                   vt) == 0);
  REQUIRE(!s.empty());
  // The graded svd fills its singular values sector by sector, so they are
  // not globally sorted; take the largest explicitly.
  double smax = 0.0;
  for (const double sv : s) {
    smax = std::max(smax, sv);
  }
  REQUIRE(smax > 0.0);
  std::size_t rank = 0;
  for (const double sv : s) {
    if (sv > 1.0e-12 * smax) {
      ++rank;
    }
  }
  return rank;
}

// C-4 accumulator: the axes observed so far per direction (0 = horizontal,
// 1 = vertical), shared by every row of the matrix. It also counts the rows
// that took the strict form of N-4, so the TEST_CASE can insist at least one
// of them ran.
struct fg15_axis_ledger {
  bool seen[2] = {false, false};
  std::size_t a[2] = {0, 0};
  std::size_t b[2] = {0, 0};
  int strict_n4_rows = 0;
};

template <class tensor>
void fg15_run_case(int d, const std::string& vps, char dir,
                   const std::string& op_kind, int op_seed, std::size_t chi,
                   int env_seed, const char* type_name,
                   fg15_axis_ledger& ledger) {
  const std::string label =
      "T15 d=" + std::to_string(d) + " vp=" + vps + " dir=" + dir +
      " op=" + op_kind + " chi=" + std::to_string(chi) + " " + type_name +
      " op_seed=" + std::to_string(op_seed) +
      " env_seed=" + std::to_string(env_seed) + " site_seeds=(61,62)";
  INFO(label);

  const fgf::parity_vector vp = fg_pv(vps);
  const fgf::parity_vector phys = ib_phys_parity(d);
  const std::size_t D = vp.size();
  const std::size_t dd = D * D;
  // Contract section 5: the ledgers must mix even and odd, and chi must
  // differ from D^2 (equal shapes would let an environment leg and a blob leg
  // be exchanged without any shape error).
  REQUIRE(fg_ledger_mixed(vp));
  REQUIRE(fg_ledger_has_odd(phys));
  REQUIRE(chi != dd);

  const fg_ftensor<tensor> TnA = fg_make_bond_site<tensor>(0, vp, phys, 61);
  const fg_ftensor<tensor> TnB = fg_make_bond_site<tensor>(1, vp, phys, 62);
  REQUIRE(fgf::parity_violation(TnA) == 0.0);
  REQUIRE(fgf::parity_violation(TnB) == 0.0);

  const tensor plain = fg_op_plain<tensor>(op_kind, d, op_seed);
  REQUIRE(fg_max_abs_entry(plain) > 0.0);
  // Production loading convention (twosite_obs.cpp): wrap_twosite_gate, i.e.
  // the input-leg swap only.
  const fg_ftensor<tensor> gate = fgf::wrap_twosite_gate(plain, phys, phys);
  const auto dir_e = dir == 'h' ? fgf::reduced_pair_direction::horizontal
                                : fgf::reduced_pair_direction::vertical;

  const fg15_env<tensor> env = fg15_make_environment<tensor>(chi, dd, env_seed);
  fg15_require_general_environment(env, label);

  // ---- reference side: the blob path that already exists.
  const tensor blob_op = fgf::build_reduced_pair_direct(TnA, TnB, gate, dir_e);
  const tensor blob_naive =
      fgf::build_reduced_pair_naive(TnA, TnB, gate, dir_e);
  const tensor blob_id = fgf::build_reduced_identity_pair(TnA, TnB, dir_e);
  REQUIRE(blob_op.shape().size() == 6);
  REQUIRE(blob_id.shape().size() == 6);
  for (std::size_t ax = 0; ax < 6; ++ax) {
    REQUIRE(blob_op.shape()[ax] == dd);
    REQUIRE(blob_id.shape()[ax] == dd);
  }
  // C-3, independent-reference clause: the operator blob must not have moved
  // against the untouched naive builder.
  fg_check_allclose(blob_op, blob_naive,
                    label + " [C-3 direct blob vs naive blob]");

  const auto ref_op = fg15_close_blob(env, dir, blob_op);
  const auto ref_id = fg15_close_blob(env, dir, blob_id);
  // N-3: neither reference value is zero.
  INFO(label << " [N-3 closure values]: op=" << ref_op << " id=" << ref_id);
  REQUIRE(std::abs(ref_op) > 0.0);
  REQUIRE(std::abs(ref_id) > 0.0);
  fg15_require_wiring_sensitivity(env, dir, blob_op, ref_op,
                                  label + " [operator blob]");
  fg15_require_wiring_sensitivity(env, dir, blob_id, ref_id,
                                  label + " [identity blob]");

  // ---- the halves under test.
  const auto halves_op = fgf::build_reduced_pair_halves(TnA, TnB, gate, dir_e);
  const auto halves_id = fgf::build_reduced_identity_halves(TnA, TnB, dir_e);

  // C-4: the shared-bond axes are a function of the direction alone. The
  // concrete values are not contracted, so they are only required to be
  // consistent - between the operator and the norm halves here, and across
  // every row of the matrix through the shared ledger.
  REQUIRE(halves_op.direction == dir_e);
  REQUIRE(halves_id.direction == dir_e);
  REQUIRE(halves_op.axis_a() == halves_id.axis_a());
  REQUIRE(halves_op.axis_b() == halves_id.axis_b());
  REQUIRE(halves_op.axis_a() < 4u);
  REQUIRE(halves_op.axis_b() < 4u);
  const int slot = dir == 'h' ? 0 : 1;
  if (ledger.seen[slot]) {
    INFO(label << " [C-4 axes vs earlier rows of the same direction]");
    REQUIRE(ledger.a[slot] == halves_op.axis_a());
    REQUIRE(ledger.b[slot] == halves_op.axis_b());
  } else {
    ledger.seen[slot] = true;
    ledger.a[slot] = halves_op.axis_a();
    ledger.b[slot] = halves_op.axis_b();
  }

  REQUIRE(halves_op.PA.shape().size() == 4);
  REQUIRE(halves_op.PB.shape().size() == 4);
  REQUIRE(halves_id.PA.shape().size() == 4);
  REQUIRE(halves_id.PB.shape().size() == 4);
  // N-3: the halves themselves are not zero.
  REQUIRE(fg_max_abs_entry(halves_op.PA) > 0.0);
  REQUIRE(fg_max_abs_entry(halves_op.PB) > 0.0);
  REQUIRE(fg_max_abs_entry(halves_id.PA) > 0.0);
  REQUIRE(fg_max_abs_entry(halves_id.PB) > 0.0);
  // N-4 (contract revision, 2026-08-29): the operator channel really is
  // bundled into the shared bond. Everywhere, the bundled dimension is a
  // whole number of virtual pairs and at least one of them; and on a gate
  // whose graded SVD carries more than one channel it is strictly fatter than
  // a bare pair.
  //
  // The first draft required "strictly fatter" unconditionally. That would
  // have become a false alarm two tasks later: once exact zero singular
  // values are dropped (a separate task, out of scope here per contract
  // section 7), a rank-1 gate bundles exactly one channel and D^2 is then the
  // correct width. The rank is measured from the gate itself rather than
  // hard-coded per operator name, so the strict rows stay strict for the
  // right reason whichever gates the matrix grows.
  const std::size_t bundled_a = halves_op.PA.shape()[halves_op.axis_a()];
  const std::size_t bundled_b = halves_op.PB.shape()[halves_op.axis_b()];
  const std::size_t gate_rank = fg15_graded_gate_rank(gate);
  INFO(label << " [N-4 shared bond]: bundled=" << bundled_a << " D^2=" << dd
             << " graded gate rank=" << gate_rank);
  REQUIRE(bundled_a == bundled_b);
  REQUIRE(bundled_a % dd == 0);
  REQUIRE(bundled_a >= dd);
  if (gate_rank > 1) {
    REQUIRE(bundled_a > dd);
    ++ledger.strict_n4_rows;
  }
  REQUIRE(halves_id.PA.shape()[halves_id.axis_a()] ==
          halves_id.PB.shape()[halves_id.axis_b()]);

  // C-3: the blob builders are the halves contracted over the shared bond.
  // The truth side is written out here rather than taken from the new API.
  const tensor rebuilt_op =
      fg_dot(halves_op.PA, halves_op.PB, mptensor::Axes(halves_op.axis_a()),
             mptensor::Axes(halves_op.axis_b()));
  fg_check_allclose(rebuilt_op, blob_op,
                    label +
                        " [C-3 tensordot(PA,PB) vs "
                        "build_reduced_pair_direct]");
  const tensor rebuilt_id =
      fg_dot(halves_id.PA, halves_id.PB, mptensor::Axes(halves_id.axis_a()),
             mptensor::Axes(halves_id.axis_b()));
  fg_check_allclose(rebuilt_id, blob_id,
                    label +
                        " [C-3 tensordot(PA,PB) vs "
                        "build_reduced_identity_pair]");

  // The common scale of the two closures of this row, taken from the
  // reference side only (see fg15_check_rel).
  const double closure_scale = std::max(std::abs(ref_op), std::abs(ref_id));

  // C-1 / C-2: the absorbing closure equals the blob closure, for the
  // operator halves and for the norm halves.
  const auto got_op = fgf::contract_reduced_pair_halves_density_CTM(
      env.C1, env.C2, env.C3, env.C4, env.eT1, env.eT2, env.eT3, env.eT4,
      env.eT5, env.eT6, halves_op);
  fg15_check_rel(label + " [C-1 operator halves closure vs blob closure]",
                 got_op, ref_op, closure_scale);
  const auto got_id = fgf::contract_reduced_pair_halves_density_CTM(
      env.C1, env.C2, env.C3, env.C4, env.eT1, env.eT2, env.eT3, env.eT4,
      env.eT5, env.eT6, halves_id);
  fg15_check_rel(label + " [C-2 identity halves closure vs blob closure]",
                 got_id, ref_id, closure_scale);
}

struct fg15_row {
  int d;
  const char* vps;
  char dir;
  const char* op;
  std::size_t chi;
  bool complex_case;
};

// Contract section 5. Both directions; both halves kinds run inside every
// row; d in {2, 4} with the mandated physical ledgers ([e,o] and [e,o,o,e],
// supplied by ib_phys_parity); three mixed virtual ledgers (eo -> D=2,
// eeo/eoo -> D=3); hopping and the diagonal nn as the two mandated gate
// shapes, plus a dense parity-even random gate; real and complex scalars;
// and chi in {3, 5}, never equal to D^2 (4 or 9).
const fg15_row fg15_rows[] = {
    {2, "eo", 'h', "hopping", 3, false},  {2, "eo", 'v', "hopping", 3, false},
    {2, "eo", 'h', "nn", 5, false},       {2, "eo", 'v', "nn", 5, false},
    {2, "eeo", 'h', "hopping", 5, false}, {2, "eeo", 'v', "nn", 5, false},
    {2, "eoo", 'h', "random", 3, false},  {2, "eo", 'h', "random", 3, true},
    {2, "eo", 'v', "hopping", 5, true},   {4, "eo", 'h', "hopping", 3, false},
    {4, "eo", 'v', "hopping", 3, false},  {4, "eo", 'h', "nn", 5, false},
    {4, "eo", 'v', "nn", 5, false},       {4, "eeo", 'h', "hopping", 5, false},
    {4, "eeo", 'v', "nn", 5, false},      {4, "eo", 'h', "hopping", 5, true},
    {4, "eo", 'v', "random", 3, true},
};

}  // namespace

TEST_CASE(
    "fold geometry T15: closing the pair halves by absorbing them into the "
    "CTM environment equals the blob closure") {
  fg15_axis_ledger ledger;
  int row_id = 0;
  for (const auto& row : fg15_rows) {
    ++row_id;
    const int env_seed = 100 + 20 * row_id;
    if (row.complex_case) {
      fg15_run_case<tenes::complex_tensor>(row.d, row.vps, row.dir, row.op, 31,
                                           row.chi, env_seed, "complex",
                                           ledger);
    } else {
      fg15_run_case<tenes::real_tensor>(row.d, row.vps, row.dir, row.op, 31,
                                        row.chi, env_seed, "real", ledger);
    }
  }
  // Both directions must have been exercised (contract section 5), and at
  // least one row must have taken the strict form of N-4 - otherwise the
  // matrix would consist entirely of single-channel gates and nothing here
  // would pin the bundling at all.
  REQUIRE(ledger.seen[0]);
  REQUIRE(ledger.seen[1]);
  REQUIRE(ledger.strict_n4_rows >= 1);
}
