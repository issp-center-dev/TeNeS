// Design-stage probe for the impurity-blob join dressing.
//
// Builds the factorized pair blob with a faithful copy of the intended
// pipeline (graded QR split -> impurity doubling with the k leg threaded
// outside the frozen joint mask -> plain join over bond and k), compares it
// elementwise against the pinned cluster path build_reduced_pair, and
// identifies the missing diagonal sign dressing automatically:
//   - 8 internal-sign candidates: (-1)^{a pk_bond + b pb_bond + c p(k)}
//     applied to impB before contraction,
//   - for each candidate, an affine F2 fit of the residual sign over the
//     12 open-leg component parity bits (+ constant).
// A candidate that fits exactly IS the dressing specification for T3.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "src/fermion/fops.hpp"
#include "src/fermion/reduced.hpp"
#include "src/tensor.hpp"

using tensor = tenes::real_tensor;
using ft = tenes::fermion::ftensor<tensor>;
namespace f = tenes::fermion;
using mptensor::Axes;
using mptensor::Index;
using mptensor::Shape;

namespace {

f::parity_vector phys_parity(int d) {
  if (d == 2) return {false, true};
  return {false, true, true, false};
}

f::parity_vector vp_mixed(int D) {
  f::parity_vector v(D, false);
  for (int i = 1; i < D; i += 2) v[i] = true;
  return v;
}

double rnd_val(int seed, const Index& idx, std::size_t rank) {
  const double w[6] = {3.1, 5.3, 7.9, 11.7, 13.1, 2.3};
  double x = 0.7 * (seed + 1) + 0.013 * seed * seed;
  for (std::size_t ax = 0; ax < rank; ++ax) x += w[ax] * (double)idx[ax];
  return 0.5 * std::sin(x) + 0.25 * std::cos(1.7 * x + 0.3 * seed);
}

ft make_site(int d, const f::parity_vector& vp, int seed) {
  const std::size_t D = vp.size();
  tensor t(Shape(D, D, D, D, d));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const Index idx = t.global_index(n);
    t.set_value(idx, rnd_val(seed, idx, 5));
  }
  return ft{t, {vp, vp, vp, vp, phys_parity(d)}};
}

tensor zeros4(int d) {
  tensor t(Shape(d, d, d, d));
  for (std::size_t n = 0; n < t.local_size(); ++n) t[n] = 0.0;
  return t;
}

tensor hop_plain_d2() {
  tensor t = zeros4(2);
  t.set_value(Index(0, 1, 1, 0), 1.0);
  t.set_value(Index(1, 0, 0, 1), 1.0);
  return t;
}

tensor nn_plain_d2() {
  tensor t = zeros4(2);
  t.set_value(Index(1, 1, 1, 1), 1.0);
  return t;
}

tensor random_even_plain(int d, int seed) {
  const f::parity_vector p = phys_parity(d);
  tensor t = zeros4(d);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const Index idx = t.global_index(n);
    const bool odd = ((p[idx[0]] != p[idx[1]]) != (p[idx[2]] != p[idx[3]]));
    if (!odd) t.set_value(idx, rnd_val(seed, idx, 4));
  }
  return t;
}

// graded QR split, both halves with leg order (in, out, k)
std::pair<ft, ft> split_op(const ft& op12) {
  ft q, r;
  f::qr(op12, Axes(0, 2), Axes(1, 3), q, r);
  ft opB = f::transpose(r, Axes(1, 2, 0));
  return {q, opB};
}

// Impurity doubling: doubled_pipeline with the k leg threaded outside the
// frozen joint mask. Returns plain (L,T,R,B,k) plus the ftensor parity info
// of the unfused boundary components via out-params.
tensor impurity(const ft& Tn, const ft& op_half) {
  ft ket = f::tensordot(Tn, op_half, Axes(4), Axes(0));  // (l,t,r,b,out,k)
  const auto forms =
      f::detail::joint_swap_forms({0, 1, 2, 3}, {5, 6, 7, 8}, {0, 1, 2, 3});
  ft bra = f::conj(Tn);
  f::apply_swap_form(bra, forms.bra);
  ft doubled = f::tensordot(bra, ket, Axes(), Axes());  // bra 0..4, ket 5..10
  Axes inter;
  for (int ax = 0; ax < 4; ++ax) {
    inter.push(5 + ax);
    inter.push(ax);
  }
  inter.push(9);   // ket out
  inter.push(4);   // bra s
  inter.push(10);  // k
  f::transpose_with_swap_form(doubled, forms.cross, inter);
  // (kl,bl, kt,bt, kr,br, kb,bb, out, s_bra, k) -> trace out x s_bra (plain,
  // same as fuse_doubled_cluster / build_reduced)
  tensor tr = mptensor::contract(doubled.t, Axes(8), Axes(9));
  // (kl,bl,kt,bt,kr,br,kb,bb,k) -> fuse boundary pairs column-major
  Shape sh;
  const Shape ts = tr.shape();
  for (int ax = 0; ax < 4; ++ax) sh.push(ts[2 * ax] * ts[2 * ax + 1]);
  sh.push(ts[8]);
  return mptensor::reshape(tr, sh);
}

struct JoinResult {
  tensor blob;
  bool ok;
  double maxdiff;
};

// candidate internal dressing bits: a*pk_bond + b*pb_bond + c*p(k)
tensor join_with_dressing(const tensor& impA, const tensor& impB_in,
                          const f::parity_vector& vp,
                          const f::parity_vector& kp, bool dir_horizontal,
                          const ft& TnA, const ft& TnB, unsigned cand) {
  const bool use_pk = cand & 1u, use_pb = cand & 2u, use_k = cand & 4u;
  tensor impB = impB_in;
  const std::size_t D = vp.size();
  // dress impB's connecting bond leg (horizontal: L = axis 0; vertical:
  // T = axis 1) and k leg (axis 4) with the candidate diagonal signs
  {
    std::vector<double> bond_sign(D * D, 1.0);
    for (std::size_t bra = 0; bra < D; ++bra)
      for (std::size_t ket = 0; ket < D; ++ket) {
        bool flip = false;
        if (use_pk && vp[ket]) flip = !flip;
        if (use_pb && vp[bra]) flip = !flip;
        bond_sign[ket + D * bra] = flip ? -1.0 : 1.0;
      }
    impB.multiply_vector(bond_sign, dir_horizontal ? 0 : 1);
  }
  if (use_k) {
    std::vector<double> ks(kp.size(), 1.0);
    for (std::size_t k = 0; k < kp.size(); ++k)
      if (kp[k]) ks[k] = -1.0;
    impB.multiply_vector(ks, 4);
  }
  tensor blob = dir_horizontal
                    ? mptensor::tensordot(impA, impB, Axes(2, 4), Axes(0, 4))
                    : mptensor::tensordot(impA, impB, Axes(3, 4), Axes(1, 4));
  // existing fused-leg gauge of the cluster path (design section 2.2)
  if (dir_horizontal) {
    f::detail::apply_fused_leg_gauge(blob, TnA.parity[3], 2, true);
    f::detail::apply_fused_leg_gauge(blob, TnB.parity[3], 5, false);
  } else {
    f::detail::apply_fused_leg_gauge(blob, TnA.parity[0], 0, true);
    f::detail::apply_fused_leg_gauge(blob, TnB.parity[0], 3, false);
  }
  return blob;
}

// Affine F2 fit: sign(idx) == (-1)^{ c0 + sum_b x_b bit_b(idx) } over the 12
// open-leg component parity bits. Gaussian elimination over F2; returns the
// solution mask description or "" if inconsistent / residual not +-1.
std::string f2_fit(const tensor& oldb, const tensor& newb,
                   const f::parity_vector& vp, double scale) {
  const std::size_t D = vp.size();
  const double tol = 1e-10 * scale;
  std::vector<std::pair<unsigned, int>> samples;  // (bits, sign)
  for (std::size_t n = 0; n < oldb.local_size(); ++n) {
    const Index idx = oldb.global_index(n);
    double vo, vn;
    oldb.get_value(idx, vo);
    newb.get_value(idx, vn);
    if (std::fabs(vo) < tol && std::fabs(vn) < tol) continue;
    if (std::fabs(vo) < tol || std::fabs(vn) < tol) return "";  // structural
    const double ratio = vo / vn;
    int sign;
    if (std::fabs(ratio - 1.0) < 1e-6)
      sign = 0;
    else if (std::fabs(ratio + 1.0) < 1e-6)
      sign = 1;
    else
      return "";  // not a pure sign pattern
    unsigned bits = 0;
    for (int leg = 0; leg < 6; ++leg) {
      const std::size_t fused = idx[leg];
      const bool pk = vp[fused % D];
      const bool pb = vp[fused / D];
      if (pk) bits |= (1u << (2 * leg));
      if (pb) bits |= (1u << (2 * leg + 1));
    }
    samples.push_back({bits, sign});
  }
  if (samples.empty()) return "trivial (no significant elements)";
  // Gaussian elimination over F2 with 13 unknowns (12 bits + constant)
  std::vector<unsigned> rows;  // bit 0..11 coeffs, bit 12 = constant, bit 13 = rhs
  for (const auto& s : samples) {
    unsigned row = s.first | (1u << 12) | ((unsigned)s.second << 13);
    // eliminate
    for (const unsigned p : rows) {
      const unsigned lead = p & ~(1u << 13);
      unsigned low = lead & (~lead + 1u);
      if (row & low) row ^= p;
    }
    if (row & ((1u << 13) - 1u)) {
      rows.push_back(row);
    } else if (row & (1u << 13)) {
      return "";  // inconsistent
    }
  }
  // back-substitute a particular solution: free vars = 0
  unsigned sol = 0;
  for (auto it = rows.rbegin(); it != rows.rend(); ++it) {
    const unsigned lead_mask = *it & ((1u << 13) - 1u);
    const unsigned low = lead_mask & (~lead_mask + 1u);
    unsigned acc = (*it >> 13) & 1u;
    for (unsigned b = 0; b < 13; ++b)
      if ((lead_mask & (1u << b)) && (1u << b) != low && (sol >> b & 1u))
        acc ^= 1u;
    if (acc) sol |= low;
  }
  static const char* leg_h[6] = {"L_A", "T_A", "B_A", "T_B", "R_B", "B_B"};
  std::string out = "sign = (-1)^{";
  bool first = true;
  for (unsigned b = 0; b < 12; ++b)
    if (sol >> b & 1u) {
      if (!first) out += " + ";
      out += std::string(leg_h[b / 2]) + (b % 2 ? ".bra" : ".ket");
      first = false;
    }
  if (sol >> 12 & 1u) {
    if (!first) out += " + ";
    out += "1";
    first = false;
  }
  if (first) out += "0";
  out += "}";
  return out;
}

void run_case(const char* name, int d, int D, const tensor& op_plain,
              bool source_second, bool dir_horizontal, int seedA, int seedB) {
  const f::parity_vector vp = D == 1 ? f::parity_vector{false} : vp_mixed(D);
  const ft TnA = make_site(d, vp, seedA);
  const ft TnB = make_site(d, vp, seedB);
  ft op = f::wrap_reduced_pair_op(op_plain, phys_parity(d), phys_parity(d));
  if (source_second) op = f::transpose(op, Axes(1, 0, 3, 2));
  const auto dir = dir_horizontal ? f::reduced_pair_direction::horizontal
                                  : f::reduced_pair_direction::vertical;
  const tensor oldb = f::build_reduced_pair(TnA, TnB, op, dir);
  double scale = 0.0;
  for (std::size_t n = 0; n < oldb.local_size(); ++n)
    scale = std::max(scale, std::fabs(oldb[n]));

  const auto halves = split_op(op);
  const tensor impA = impurity(TnA, halves.first);
  const tensor impB = impurity(TnB, halves.second);
  const f::parity_vector kp = halves.first.parity[2];

  std::printf("%-28s d=%d D=%d dir=%c src=%d:", name, d, D,
              dir_horizontal ? 'h' : 'v', source_second ? 2 : 1);
  bool any = false;
  for (unsigned cand = 0; cand < 8; ++cand) {
    const tensor newb = join_with_dressing(impA, impB, vp, kp, dir_horizontal,
                                           TnA, TnB, cand);
    // quick exact check first
    double maxdiff = 0.0;
    for (std::size_t n = 0; n < oldb.local_size(); ++n) {
      const Index idx = oldb.global_index(n);
      double vo, vn;
      oldb.get_value(idx, vo);
      newb.get_value(idx, vn);
      maxdiff = std::max(maxdiff, std::fabs(vo - vn));
    }
    if (maxdiff <= 1e-12 * std::max(1.0, scale)) {
      std::printf("  [cand=%u%s%s%s EXACT]", cand, cand & 1 ? " pkB" : "",
                  cand & 2 ? " pbB" : "", cand & 4 ? " pK" : "");
      any = true;
      continue;
    }
    const std::string fit = f2_fit(oldb, newb, vp, scale);
    if (!fit.empty()) {
      std::printf("  [cand=%u%s%s%s + open-leg %s]", cand,
                  cand & 1 ? " pkB" : "", cand & 2 ? " pbB" : "",
                  cand & 4 ? " pK" : "", fit.c_str());
      any = true;
    }
  }
  if (!any) std::printf("  NO CANDIDATE MATCHES (structural mismatch)");
  std::printf("\n");
}

}  // namespace

int main() {
  // ladder (a): even-sector op, D=1 then D=2
  run_case("nn (even k)", 2, 1, nn_plain_d2(), false, true, 1, 2);
  run_case("nn (even k)", 2, 2, nn_plain_d2(), false, true, 1, 2);
  run_case("nn (even k)", 2, 2, nn_plain_d2(), false, false, 1, 2);
  // ladder (b): odd-sector op, D=1 (bond trivial)
  run_case("hopping (odd k)", 2, 1, hop_plain_d2(), false, true, 3, 4);
  run_case("hopping (odd k)", 2, 1, hop_plain_d2(), false, false, 3, 4);
  // ladder (c): odd-sector op, mixed bond
  run_case("hopping (odd k)", 2, 2, hop_plain_d2(), false, true, 5, 6);
  run_case("hopping (odd k)", 2, 2, hop_plain_d2(), false, false, 5, 6);
  run_case("hopping src2", 2, 2, hop_plain_d2(), true, true, 7, 8);
  run_case("hopping src2", 2, 2, hop_plain_d2(), true, false, 7, 8);
  // mixed-sector random even op (k even and odd together)
  run_case("random even", 2, 2, random_even_plain(2, 21), false, true, 9, 10);
  run_case("random even", 2, 2, random_even_plain(2, 21), false, false, 9, 10);
  run_case("random even src2", 2, 2, random_even_plain(2, 22), true, false, 11,
           12);
  // d=4 spot checks
  run_case("d4 random even", 4, 2, random_even_plain(4, 23), false, true, 13,
           14);
  run_case("d4 random even", 4, 2, random_even_plain(4, 23), false, false, 13,
           14);
  return 0;
}
