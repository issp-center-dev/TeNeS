// Decisive probe for the "lean fuse" identity:
//
//   fuse_doubled_cluster's rank-16 outer product + interleave transpose +
//   fuse + physical trace  ==  direct plain contraction of the physical legs
//   (rank-12 intermediate, blob-sized) with the swap-form terms redistributed:
//     - bra-internal terms       -> mask on the bra pair tensor,
//     - ket-internal terms       -> mask on the ket pair tensor,
//     - cross terms hitting a traced leg -> rewritten onto its delta-trace
//       twin (bra s -> ket s', and (bra s, ket s') -> linear parity mask),
//     - cross open-open terms    -> mask on the rank-12 result.
//
// If this reproduces build_reduced_pair elementwise, the memory fix needs no
// new sign conventions at all.
#include <cmath>
#include <cstdio>
#include <utility>
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

tensor pairing_plain_d2() {
  tensor t = zeros4(2);
  t.set_value(Index(1, 1, 0, 0), 1.0);
  t.set_value(Index(0, 0, 1, 1), 1.0);
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

// elementwise sign mask (-1)^{p_i p_j} for a list of axis pairs, on a plain
// tensor with explicit leg parities (probe-local; small sizes only)
void apply_pair_terms(tensor& t, const f::leg_parities& parity,
                      const std::vector<std::pair<int, int>>& terms) {
  if (terms.empty()) return;
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const Index idx = t.global_index(n);
    int flips = 0;
    for (const auto& pr : terms) {
      if (parity[pr.first][idx[pr.first]] && parity[pr.second][idx[pr.second]])
        ++flips;
    }
    if (flips % 2) t[n] = -t[n];
  }
}

void apply_parity_axes(tensor& t, const f::leg_parities& parity,
                       const std::vector<int>& axes) {
  if (axes.empty()) return;
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const Index idx = t.global_index(n);
    int flips = 0;
    for (const int ax : axes)
      if (parity[ax][idx[ax]]) ++flips;
    if (flips % 2) t[n] = -t[n];
  }
}

// Lean replacement for fuse_doubled_cluster(conj(ket_ab), ket_op, leg_ids).
tensor lean_fuse(const ft& bra_pair, const ft& ket_pair,
                 const std::vector<int>& leg_ids) {
  const std::vector<int> bra_axes = {0, 1, 2, 4, 5, 6};
  std::vector<int> ket_axes;
  for (const int ax : bra_axes) ket_axes.push_back(ax + 8);
  const auto forms = f::detail::joint_swap_forms(bra_axes, ket_axes, leg_ids);

  // the 16-leg interleave the cluster path uses
  Axes interleaved;
  for (std::size_t i = 0; i < 6; ++i) {
    interleaved.push(ket_axes[i]);
    interleaved.push(bra_axes[i]);
  }
  interleaved.push(11);
  interleaved.push(15);
  interleaved.push(3);
  interleaved.push(7);

  // total pairwise sign terms of the cluster path after the outer product:
  // forms.cross plus the Koszul terms of the interleave transpose
  f::SwapForm total = forms.cross;
  {
    const f::SwapForm koszul = f::detail::transpose_sign_form(interleaved);
    for (const auto& pr : koszul.terms()) total.toggle(pr.first, pr.second);
  }

  // classify: axes 0..7 bra, 8..15 ket; traced legs bra{3,7} <-> ket{11,15}
  auto is_bra = [](int ax) { return ax < 8; };
  auto trace_twin = [](int ax) {  // defined only for traced legs
    if (ax == 3) return 11;
    if (ax == 7) return 15;
    if (ax == 11) return 3;
    if (ax == 15) return 7;
    return -1;
  };
  std::vector<std::pair<int, int>> bra_terms, ket_terms, post_terms;
  std::vector<int> ket_parity_axes;  // linear parity masks on ket
  auto to_ket_local = [](int ax) { return ax - 8; };
  // positions of open legs in the rank-12 contraction result
  // tensordot(bra', ket', Axes(3,7), Axes(3,7)):
  //   result = (bra 0,1,2,4,5,6 | ket 0,1,2,4,5,6) in original relative order
  auto post_pos = [&](int ax) {
    static const int bra_pos[8] = {0, 1, 2, -1, 3, 4, 5, -1};
    if (is_bra(ax)) return bra_pos[ax];
    return 6 + bra_pos[ax - 8];
  };
  for (const auto& pr : total.terms()) {
    int a = pr.first, b = pr.second;
    const bool a_traced = (a == 3 || a == 7 || a == 11 || a == 15);
    const bool b_traced = (b == 3 || b == 7 || b == 11 || b == 15);
    if (a_traced && b_traced && trace_twin(a) == b) {
      // p(s)*p(s') with s == s' under the delta trace: a linear parity mask
      ket_parity_axes.push_back(to_ket_local(std::max(a, b)));
      continue;
    }
    // rewrite any traced BRA leg onto its ket twin (delta trace equates them)
    if (a_traced && is_bra(a)) a = trace_twin(a);
    if (b_traced && is_bra(b)) b = trace_twin(b);
    // now traced legs, if any, are ket-side (11 or 15)
    if (is_bra(a) && is_bra(b)) {
      bra_terms.push_back({a, b});
    } else if (!is_bra(a) && !is_bra(b)) {
      ket_terms.push_back({to_ket_local(a), to_ket_local(b)});
    } else {
      // cross open-open (traced legs can no longer appear on the bra side,
      // and a ket-traced leg paired with a bra open leg was rewritten? no:
      // rewrite only moved bra->ket. ket traced x bra open stays cross; move
      // it to the bra twin instead.
      int ka = is_bra(a) ? b : a;  // ket-side axis
      int br = is_bra(a) ? a : b;  // bra-side axis
      if (ka == 11 || ka == 15) {
        bra_terms.push_back({trace_twin(ka), br});
      } else {
        post_terms.push_back({post_pos(br), post_pos(ka)});
      }
    }
  }

  ft bra = bra_pair;
  f::apply_swap_form(bra, forms.bra);
  tensor bt = bra.t;
  apply_pair_terms(bt, bra.parity, bra_terms);
  tensor kt = ket_pair.t;
  apply_pair_terms(kt, ket_pair.parity, ket_terms);
  apply_parity_axes(kt, ket_pair.parity, ket_parity_axes);

  tensor lean = mptensor::tensordot(bt, kt, Axes(3, 7), Axes(3, 7));
  // parities of the rank-12 result for the post masks
  f::leg_parities post_parity;
  for (const int ax : bra_axes) post_parity.push_back(bra.parity[ax]);
  for (const int ax : bra_axes) post_parity.push_back(ket_pair.parity[ax]);
  apply_pair_terms(lean, post_parity, post_terms);

  // interleave (ket_i, bra_i) and fuse, matching the cluster output
  Axes perm;
  for (int i = 0; i < 6; ++i) {
    perm.push(6 + i);
    perm.push(i);
  }
  lean = mptensor::transpose(lean, perm);
  Shape sh;
  const Shape ls = lean.shape();
  for (int i = 0; i < 6; ++i) sh.push(ls[2 * i] * ls[2 * i + 1]);
  return mptensor::reshape(lean, sh);
}

// full lean build_reduced_pair (mirrors reduced.hpp build_reduced_pair with
// lean_fuse in place of fuse_doubled_cluster)
tensor lean_pair(const ft& TnA, const ft& TnB, const ft& op12,
                 f::reduced_pair_direction dir) {
  const ft ket_ab = f::build_pair_state(TnA, TnB, dir);
  std::vector<int> leg_ids;
  if (dir == f::reduced_pair_direction::horizontal)
    leg_ids = {0, 1, 3, 1, 2, 3};
  else
    leg_ids = {0, 1, 2, 0, 2, 3};
  ft ket_op = f::apply_pair_op(ket_ab, op12);
  tensor ret = lean_fuse(f::conj(ket_ab), ket_op, leg_ids);
  if (dir == f::reduced_pair_direction::horizontal) {
    f::detail::apply_fused_leg_gauge(ret, TnA.parity[3], 2, true);
    f::detail::apply_fused_leg_gauge(ret, TnB.parity[3], 5, false);
  } else {
    f::detail::apply_fused_leg_gauge(ret, TnA.parity[0], 0, true);
    f::detail::apply_fused_leg_gauge(ret, TnB.parity[0], 3, false);
  }
  return ret;
}

int g_fail = 0;

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
  const tensor newb = lean_pair(TnA, TnB, op, dir);
  double scale = 0.0, maxdiff = 0.0;
  for (std::size_t n = 0; n < oldb.local_size(); ++n) {
    const Index idx = oldb.global_index(n);
    double vo, vn;
    oldb.get_value(idx, vo);
    newb.get_value(idx, vn);
    scale = std::max(scale, std::fabs(vo));
    maxdiff = std::max(maxdiff, std::fabs(vo - vn));
  }
  const bool ok = maxdiff <= 1e-12 * std::max(1.0, scale);
  if (!ok) ++g_fail;
  std::printf("%-22s d=%d D=%d dir=%c src=%d  maxdiff=%.3e  %s\n", name, d, D,
              dir_horizontal ? 'h' : 'v', source_second ? 2 : 1, maxdiff,
              ok ? "OK" : "MISMATCH");
}

}  // namespace

int main() {
  for (const bool h : {true, false}) {
    run_case("nn", 2, 1, nn_plain_d2(), false, h, 1, 2);
    run_case("hopping", 2, 1, hop_plain_d2(), false, h, 3, 4);
    run_case("nn", 2, 2, nn_plain_d2(), false, h, 5, 6);
    run_case("hopping", 2, 2, hop_plain_d2(), false, h, 7, 8);
    run_case("hopping src2", 2, 2, hop_plain_d2(), true, h, 9, 10);
    run_case("pairing", 2, 2, pairing_plain_d2(), false, h, 11, 12);
    run_case("random even", 2, 2, random_even_plain(2, 21), false, h, 13, 14);
    run_case("random even src2", 2, 2, random_even_plain(2, 22), true, h, 15,
             16);
    run_case("random even", 2, 3, random_even_plain(2, 23), false, h, 17, 18);
    run_case("d4 random even", 4, 2, random_even_plain(4, 24), false, h, 19,
             20);
    run_case("d4 random src2", 4, 2, random_even_plain(4, 25), true, h, 21,
             22);
  }
  std::printf("%s\n", g_fail == 0 ? "LEAN FUSE IDENTITY HOLDS"
                                  : "LEAN FUSE IDENTITY FAILS");
  return g_fail == 0 ? 0 : 1;
}
