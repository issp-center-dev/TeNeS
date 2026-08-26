// Design-stage probe: does the graded QR split of a wrapped two-site
// operator (rows {0,2} = A legs, cols {1,3} = B legs) recombine to the
// original elementwise, and does the internal leg carry the [F..,T..]
// parity ladder?  Uses only existing primitives; no product code changes.
#include <cmath>
#include <cstdio>
#include <vector>

#include "src/fermion/fops.hpp"
#include "src/tensor.hpp"

using tensor = tenes::real_tensor;
using ft = tenes::fermion::ftensor<tensor>;
namespace f = tenes::fermion;

namespace {

// d=2 spinless: parity [even, odd]
f::parity_vector p_d2() { return {false, true}; }
// d=4 electron: parity [0,1,1,0]
f::parity_vector p_d4() { return {false, true, true, false}; }

tensor zeros(std::size_t d) {
  tensor t(mptensor::Shape(d, d, d, d));
  for (std::size_t n = 0; n < t.local_size(); ++n) t[n] = 0.0;
  return t;
}

void set4(tensor& t, std::size_t a, std::size_t b, std::size_t c,
          std::size_t d, double v) {
  t.set_value(mptensor::Index(a, b, c, d), v);
}

// d=2 hopping  c^dag_A c_B + c^dag_B c_A  (in_A,in_B,out_A,out_B)
tensor hop_d2() {
  tensor t = zeros(2);
  set4(t, 0, 1, 1, 0, 1.0);
  set4(t, 1, 0, 0, 1, 1.0);
  return t;
}

// d=2 density n_A n_B
tensor nn_d2() {
  tensor t = zeros(2);
  set4(t, 1, 1, 1, 1, 1.0);
  return t;
}

// random parity-even op for dimension d with parity p
tensor random_even(const f::parity_vector& p, unsigned seed) {
  const std::size_t d = p.size();
  tensor t = zeros(d);
  unsigned s = seed;
  auto rnd = [&s]() {
    s = s * 1664525u + 1013904223u;
    return (double)(s >> 8) / (double)(1u << 24) - 0.5;
  };
  for (std::size_t a = 0; a < d; ++a)
    for (std::size_t b = 0; b < d; ++b)
      for (std::size_t c = 0; c < d; ++c)
        for (std::size_t e = 0; e < d; ++e) {
          const bool odd = (p[a] ^ p[b] ^ p[c] ^ p[e]);
          if (!odd) set4(t, a, b, c, e, rnd());
        }
  return t;
}

double check_roundtrip(const tensor& op, const f::parity_vector& p1,
                       const f::parity_vector& p2, const char* name) {
  ft fop = f::wrap_reduced_pair_op(op, p1, p2);
  // rows = A legs (in_A=0, out_A=2), cols = B legs (in_B=1, out_B=3)
  ft q, r;
  const int info = f::qr(fop, mptensor::Axes(0, 2), mptensor::Axes(1, 3), q, r);
  // q: (in_A, out_A, k)   r: (k, in_B, out_B)
  // recombine over k and restore leg order (in_A,in_B,out_A,out_B)
  ft rec = f::tensordot(q, r, mptensor::Axes(2), mptensor::Axes(0));
  rec = f::transpose(rec, mptensor::Axes(0, 2, 1, 3));
  double maxdiff = 0.0;
  const std::size_t d1 = p1.size(), d2v = p2.size();
  for (std::size_t a = 0; a < d1; ++a)
    for (std::size_t b = 0; b < d2v; ++b)
      for (std::size_t c = 0; c < d1; ++c)
        for (std::size_t e = 0; e < d2v; ++e) {
          double va, vb;
          fop.get_value(mptensor::Index(a, b, c, e), va);
          rec.get_value(mptensor::Index(a, b, c, e), vb);
          maxdiff = std::max(maxdiff, std::fabs(va - vb));
        }
  // parity ladder of k
  const auto& kp = q.parity[2];
  bool ladder = true;
  for (std::size_t i = 1; i < kp.size(); ++i)
    if (kp[i - 1] && !kp[i]) ladder = false;
  std::printf("%-16s k_dim=%zu ladder=%s maxdiff=%.3e info=%d\n", name,
              kp.size(), ladder ? "ok" : "BROKEN", maxdiff, info);
  return maxdiff;
}

}  // namespace

int main() {
  double worst = 0.0;
  worst = std::max(worst, check_roundtrip(hop_d2(), p_d2(), p_d2(), "d2 hop"));
  worst = std::max(worst, check_roundtrip(nn_d2(), p_d2(), p_d2(), "d2 nn"));
  worst = std::max(
      worst, check_roundtrip(random_even(p_d2(), 7u), p_d2(), p_d2(), "d2 rand"));
  worst = std::max(
      worst, check_roundtrip(random_even(p_d4(), 3u), p_d4(), p_d4(), "d4 rand"));
  worst = std::max(worst, check_roundtrip(random_even(p_d4(), 11u), p_d4(),
                                          p_d4(), "d4 rand2"));
  std::printf("worst = %.3e -> %s\n", worst,
              worst < 1e-12 ? "SPLIT PREMISE OK" : "SPLIT PREMISE FAILS");
  return worst < 1e-12 ? 0 : 1;
}
