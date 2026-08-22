// ===== SS: fused fermionic sign-sweep kernel ==============================
//
// Pins tenes::fermion::SwapForm / apply_swap_form / apply_sign_sweep /
// apply_leg_gauges / detail::use_sign_table (src/fermion/sign_sweep.hpp)
// against the task-1 contract:
//   C1-1 toggle algebra of SwapForm
//   C1-2 agreement with a sequence of apply_swap calls (order independent)
//   C1-3 cancellation of a pair given in both orientations -- the shape
//        apply_joint_swaps() in src/fermion/reduced.hpp actually produces
//        ({bra0,bra3} appears as (bra0,bra3) and as (bra3,bra0))
//   C1-4 the empty form is the identity
//   C1-5 mptensor index decode equivalence (make_l2g_map/global_index_l2g_map
//        vs. global_index_fast) -- the premise of the whole optimisation,
//        deliberately independent of the new header
//   C1-6 table / direct / automatic evaluation paths agree
//   C1-7 boundaries of use_sign_table()
//   C1-8 diagonal leg gauges vs. mptensor multiply_vector
//   C1-9 non-regression of the existing apply_swap, including the equal-axis
//        case apply_swap(a, x, x) == apply_parity(a, x) (ruling R6)
//
// Every comparison is exact (==). The swap signs are exact +-1, and the gauge
// factors are compared against the same sequence of multiply_vector calls the
// kernel is required to reproduce (gauges applied one at a time, in list
// order, never pre-multiplied), so no rounding difference can arise. The
// factors in C1-8 are deliberately NOT powers of two, so that a kernel which
// folds them together before the sweep is detected.

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../../src/fermion/sign_sweep.hpp"

namespace fermion_sign_sweep {

using tenes::fermion::ftensor;
using tenes::fermion::leg_parities;
using tenes::fermion::LegGauge;
using tenes::fermion::parity_vector;
using tenes::fermion::SignEval;
using tenes::fermion::SwapForm;

using pair_list = std::vector<std::pair<int, int>>;

// ---------------------------------------------------------------- helpers

// Per-axis parity ledgers that differ from axis to axis: the positions (and
// for dim >= 3 the number) of the odd entries rotate with the axis index, so
// a kernel that reads the ledger of the wrong axis cannot reproduce the
// reference signs.  Every axis of dimension >= 2 carries both parities, so
// no swap in these tests degenerates into a no-op.
inline parity_vector ss_parity(std::size_t axis, std::size_t dim) {
  static const unsigned masks[8] = {0x2u, 0x1u, 0x6u, 0x9u,
                                    0xau, 0x5u, 0xcu, 0x3u};
  const unsigned m = masks[axis % 8];
  parity_vector p(dim);
  for (std::size_t i = 0; i < dim; ++i) {
    p[i] = ((m >> (i % 4)) & 1u) != 0;
  }
  if (dim >= 2) {
    bool all_same = true;
    for (std::size_t i = 1; i < dim; ++i) {
      if (p[i] != p[0]) {
        all_same = false;
        break;
      }
    }
    if (all_same) {
      p[dim - 1] = !p[0];
    }
  }
  return p;
}

inline leg_parities ss_parities(const mptensor::Shape& sh) {
  leg_parities p;
  for (std::size_t ax = 0; ax < sh.size(); ++ax) {
    p.push_back(ss_parity(ax, sh[ax]));
  }
  return p;
}

template <class T>
struct ss_sampler;

template <>
struct ss_sampler<double> {
  static double get(std::mt19937& gen) {
    std::uniform_real_distribution<double> mag(0.5, 1.5);
    std::uniform_int_distribution<int> sgn(0, 1);
    return (sgn(gen) == 0 ? -1.0 : 1.0) * mag(gen);
  }
};

template <>
struct ss_sampler<std::complex<double>> {
  static std::complex<double> get(std::mt19937& gen) {
    const double re = ss_sampler<double>::get(gen);
    const double im = ss_sampler<double>::get(gen);
    return std::complex<double>(re, im);
  }
};

// Every element has magnitude >= 0.5, so a wrongly (un)flipped sign is
// always visible; a tensor with zeros could hide one.
template <class tensor>
tensor ss_random_tensor(const mptensor::Shape& sh, unsigned seed) {
  using value_type = typename tensor::value_type;
  tensor t(sh);
  std::mt19937 gen(seed);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    t.set_value(t.global_index(n), ss_sampler<value_type>::get(gen));
  }
  return t;
}

template <class tensor>
ftensor<tensor> ss_random_ftensor(const mptensor::Shape& sh, unsigned seed) {
  return ftensor<tensor>{ss_random_tensor<tensor>(sh, seed), ss_parities(sh)};
}

template <class tensor>
std::size_t ss_count_diff(const tensor& a, const tensor& b) {
  std::size_t bad = 0;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    const mptensor::Index idx = a.global_index(n);
    typename tensor::value_type va, vb;
    a.get_value(idx, va);
    b.get_value(idx, vb);
    if (!(va == vb)) {
      ++bad;
    }
  }
  return bad;
}

template <class tensor>
ftensor<tensor> ss_sequential_swaps(const ftensor<tensor>& a,
                                    const pair_list& pairs) {
  ftensor<tensor> ret = a;
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    tenes::fermion::apply_swap(ret, pairs[i].first, pairs[i].second);
  }
  return ret;
}

inline SwapForm ss_form(const pair_list& pairs) {
  SwapForm form;
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    form.toggle(pairs[i].first, pairs[i].second);
  }
  return form;
}

// Reference for the gauges: mptensor's own multiply_vector, one factor at a
// time, in the order the gauges are listed.
template <class tensor>
tensor ss_sequential_gauges(const tensor& a,
                            const std::vector<LegGauge>& gauges) {
  tensor ret = a;
  for (std::size_t i = 0; i < gauges.size(); ++i) {
    ret.multiply_vector(gauges[i].factor,
                        static_cast<std::size_t>(gauges[i].axis));
  }
  return ret;
}

// The kernel must apply the gauges one after another, in the order they appear
// in the vector, exactly as the reference multiply_vector calls do -- it may
// NOT fold several factors of one axis into a pre-multiplied product.  The
// sequence of floating point operations is then identical, so `==` is
// legitimate for arbitrary factors.  ss_generic() therefore deliberately
// avoids powers of two: a kernel that pre-multiplies would round differently
// and is caught here.
inline std::vector<double> ss_pm_one(std::size_t dim, unsigned seed) {
  std::vector<double> v(dim);
  for (std::size_t i = 0; i < dim; ++i) {
    v[i] = ((i + seed) % 2 == 0) ? 1.0 : -1.0;
  }
  return v;
}

inline std::vector<double> ss_generic(std::size_t dim, unsigned seed) {
  static const double values[6] = {0.1, 3.7, -1.3, 0.9, -0.7, 1.1};
  std::vector<double> v(dim);
  for (std::size_t i = 0; i < dim; ++i) {
    v[i] = values[(i + seed) % 6];
  }
  return v;
}

}  // namespace fermion_sign_sweep

// ------------------------------------------------------------------ C1-1

TEST_CASE("SS C1-1 SwapForm toggle cancels pairs and normalises them") {
  using fermion_sign_sweep::pair_list;
  using tenes::fermion::SwapForm;

  SUBCASE("a fresh form is empty") {
    SwapForm form;
    CHECK(form.empty());
    CHECK(form.terms() == pair_list{});
  }

  SUBCASE("one toggle stores the pair with first < second") {
    SwapForm form;
    form.toggle(3, 1);
    CHECK_FALSE(form.empty());
    CHECK(form.terms() == pair_list{{1, 3}});
  }

  SUBCASE("toggle(x,y) and toggle(y,x) are the same unordered pair") {
    SwapForm form;
    form.toggle(3, 1);
    form.toggle(1, 3);
    CHECK(form.empty());
    CHECK(form.terms() == pair_list{});
  }

  SUBCASE("an even number of toggles cancels, an odd number does not") {
    for (int repeat = 1; repeat <= 5; ++repeat) {
      SwapForm form;
      for (int k = 0; k < repeat; ++k) {
        // alternate the orientation so the cancellation cannot depend on it
        if (k % 2 == 0) {
          form.toggle(2, 6);
        } else {
          form.toggle(6, 2);
        }
      }
      INFO("repeat = " << repeat);
      if (repeat % 2 == 0) {
        CHECK(form.empty());
        CHECK(form.terms() == pair_list{});
      } else {
        CHECK_FALSE(form.empty());
        CHECK(form.terms() == pair_list{{2, 6}});
      }
    }
  }

  SUBCASE("toggle(x,x) does not change the form") {
    SwapForm empty_form;
    empty_form.toggle(4, 4);
    CHECK(empty_form.empty());
    CHECK(empty_form.terms() == pair_list{});

    SwapForm form;
    form.toggle(4, 0);
    form.toggle(2, 2);
    form.toggle(0, 0);
    form.toggle(7, 7);
    CHECK(form.terms() == pair_list{{0, 4}});
  }

  SUBCASE("terms() is sorted ascending by (first, second)") {
    SwapForm form;
    form.toggle(5, 2);
    form.toggle(1, 0);
    form.toggle(4, 3);
    form.toggle(3, 2);
    form.toggle(0, 6);
    CHECK(form.terms() == pair_list{{0, 1}, {0, 6}, {2, 3}, {2, 5}, {3, 4}});
  }

  SUBCASE("cancellation in the middle of a longer form") {
    SwapForm form;
    form.toggle(0, 2);
    form.toggle(1, 3);
    form.toggle(2, 0);
    form.toggle(4, 1);
    CHECK(form.terms() == pair_list{{1, 3}, {1, 4}});
  }
}

// ------------------------------------------------------------------ C1-2

TEST_CASE_TEMPLATE("SS C1-2 apply_swap_form matches repeated apply_swap",
                   tensor, tenes::real_tensor, tenes::complex_tensor) {
  using namespace fermion_sign_sweep;

  const mptensor::Shape sh(2, 3, 2, 4, 3);
  const ftensor<tensor> a = ss_random_ftensor<tensor>(sh, 20260822u);

  // four distinct unordered pairs
  const pair_list pairs{{0, 2}, {1, 4}, {3, 0}, {2, 4}};
  const ftensor<tensor> expected = ss_sequential_swaps(a, pairs);

  // the reference must actually differ from the input, otherwise the whole
  // comparison below would be vacuous
  REQUIRE(ss_count_diff(expected.t, a.t) > 0);

  ftensor<tensor> got = a;
  tenes::fermion::apply_swap_form(got, ss_form(pairs));
  CHECK(ss_count_diff(got.t, expected.t) == 0);
  CHECK(got.parity == a.parity);

  SUBCASE("the order in which the pairs are toggled does not matter") {
    const pair_list permuted[3] = {
        {{2, 4}, {0, 3}, {4, 1}, {2, 0}},
        {{4, 1}, {2, 0}, {2, 4}, {3, 0}},
        {{3, 0}, {4, 2}, {0, 2}, {1, 4}},
    };
    for (int i = 0; i < 3; ++i) {
      INFO("permutation " << i);
      ftensor<tensor> other = a;
      tenes::fermion::apply_swap_form(other, ss_form(permuted[i]));
      CHECK(ss_count_diff(other.t, expected.t) == 0);
    }
  }

  SUBCASE("a single pair reproduces apply_swap") {
    const pair_list one{{1, 3}};
    ftensor<tensor> single = a;
    tenes::fermion::apply_swap_form(single, ss_form(one));
    CHECK(ss_count_diff(single.t, ss_sequential_swaps(a, one).t) == 0);
  }

  SUBCASE("applying the same form twice restores the input") {
    ftensor<tensor> twice = a;
    tenes::fermion::apply_swap_form(twice, ss_form(pairs));
    tenes::fermion::apply_swap_form(twice, ss_form(pairs));
    CHECK(ss_count_diff(twice.t, a.t) == 0);
  }
}

// ------------------------------------------------------------------ C1-3

TEST_CASE_TEMPLATE("SS C1-3 a pair given in both orientations cancels", tensor,
                   tenes::real_tensor, tenes::complex_tensor) {
  using namespace fermion_sign_sweep;

  const mptensor::Shape sh(3, 2, 4, 2, 3);
  const ftensor<tensor> a = ss_random_ftensor<tensor>(sh, 4711u);

  SUBCASE("(0,2) and (2,0) drop out of the middle of a list") {
    const pair_list with_dup{{0, 2}, {1, 3}, {2, 0}, {4, 1}};
    const pair_list without{{1, 3}, {4, 1}};
    const ftensor<tensor> expected = ss_sequential_swaps(a, without);
    REQUIRE(ss_count_diff(expected.t, a.t) > 0);

    ftensor<tensor> got = a;
    tenes::fermion::apply_swap_form(got, ss_form(with_dup));
    CHECK(ss_count_diff(got.t, expected.t) == 0);
  }

  SUBCASE("the surviving pair list can become empty") {
    const pair_list all_cancel{{0, 2}, {1, 4}, {2, 0}, {4, 1}};
    const SwapForm form = ss_form(all_cancel);
    CHECK(form.empty());
    ftensor<tensor> got = a;
    tenes::fermion::apply_swap_form(got, form);
    CHECK(ss_count_diff(got.t, a.t) == 0);
  }

  SUBCASE("the shape produced by apply_joint_swaps in reduced.hpp") {
    // {b0,b3} appears as (b0,b3) and as (b3,b0) and must cancel; the other
    // two bra pairs survive.  Axes 0..4 stand for b0..b3 plus one ket leg.
    const pair_list produced{{0, 3}, {1, 0}, {2, 3}, {3, 0}, {4, 1}};
    const pair_list surviving{{1, 0}, {2, 3}, {4, 1}};
    const ftensor<tensor> expected = ss_sequential_swaps(a, surviving);
    REQUIRE(ss_count_diff(expected.t, a.t) > 0);

    ftensor<tensor> got = a;
    tenes::fermion::apply_swap_form(got, ss_form(produced));
    CHECK(ss_count_diff(got.t, expected.t) == 0);
  }
}

// ------------------------------------------------------------------ C1-4

TEST_CASE_TEMPLATE("SS C1-4 an empty SwapForm leaves every element alone",
                   tensor, tenes::real_tensor, tenes::complex_tensor) {
  using namespace fermion_sign_sweep;

  const mptensor::Shape sh(2, 3, 4, 2);
  const ftensor<tensor> a = ss_random_ftensor<tensor>(sh, 99u);

  const SwapForm empty_form;
  REQUIRE(empty_form.empty());

  ftensor<tensor> got = a;
  tenes::fermion::apply_swap_form(got, empty_form);
  CHECK(ss_count_diff(got.t, a.t) == 0);
  CHECK(got.parity == a.parity);

  ftensor<tensor> swept = a;
  tenes::fermion::apply_sign_sweep(swept, empty_form, std::vector<LegGauge>{});
  CHECK(ss_count_diff(swept.t, a.t) == 0);
  CHECK(swept.parity == a.parity);

  SUBCASE("a form emptied by cancellation behaves the same") {
    SwapForm cancelled;
    cancelled.toggle(1, 3);
    cancelled.toggle(3, 1);
    REQUIRE(cancelled.empty());
    ftensor<tensor> other = a;
    tenes::fermion::apply_swap_form(other, cancelled);
    CHECK(ss_count_diff(other.t, a.t) == 0);
  }
}

// ------------------------------------------------------------------ C1-5

TEST_CASE_TEMPLATE("SS C1-5 global_index_l2g_map agrees with global_index_fast",
                   tensor, tenes::real_tensor, tenes::complex_tensor) {
  const std::vector<std::vector<std::size_t>> shapes = {
      {3, 5},
      {2, 3, 4},
      {2, 3, 4, 2},
      {2, 3, 2, 4, 3},
      {2, 3, 4, 2, 3, 2},
      {3, 2, 2, 4, 2, 3, 2},
      {2, 3, 2, 2, 3, 2, 2, 2},
      {2, 2, 3, 2, 2, 2, 2, 3, 2, 2},
      {2, 2, 3, 2, 2, 2, 3, 2, 2, 2, 2, 2},
      {2, 2, 2, 3, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2},
      {2, 2, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2},
  };

  for (std::size_t s = 0; s < shapes.size(); ++s) {
    const std::vector<std::size_t>& dims = shapes[s];
    const std::size_t rank = dims.size();
    INFO("rank = " << rank);

    const mptensor::Shape sh(dims);
    tensor t(sh);
    REQUIRE(t.rank() == rank);
    REQUIRE(t.local_size() > 0);

    t.make_l2g_map();

    mptensor::Index fast;
    fast.resize(rank);
    std::vector<std::size_t> mapped(rank, 0);

    std::size_t bad = 0;
    for (std::size_t n = 0; n < t.local_size(); ++n) {
      t.global_index_fast(n, fast);
      t.global_index_l2g_map(n, &mapped[0]);
      for (std::size_t ax = 0; ax < rank; ++ax) {
        if (fast[ax] != mapped[ax]) {
          ++bad;
        }
      }
    }
    CHECK(bad == 0);
  }
}

// ------------------------------------------------------------------ C1-6

TEST_CASE_TEMPLATE("SS C1-6 table, direct and automatic evaluation agree",
                   tensor, tenes::real_tensor, tenes::complex_tensor) {
  using namespace fermion_sign_sweep;

  // shape 1: local_size 144 >= 2^5, so `automatic` may build a table
  // shape 2: local_size 32 < 2^6 (one axis of dimension 1), so `automatic`
  //          has to evaluate directly.  Both branches are exercised.
  const mptensor::Shape shapes[2] = {mptensor::Shape(2, 3, 2, 4, 3),
                                     mptensor::Shape(2, 2, 2, 2, 2, 1)};
  const pair_list pairs[2] = {{{0, 2}, {1, 4}, {3, 0}, {2, 4}},
                              {{0, 5}, {1, 3}, {2, 4}, {4, 0}}};

  for (int s = 0; s < 2; ++s) {
    INFO("shape index " << s);
    const ftensor<tensor> a = ss_random_ftensor<tensor>(shapes[s], 31337u + s);
    const SwapForm form = ss_form(pairs[s]);

    std::vector<LegGauge> gauges;
    gauges.push_back(LegGauge{1, ss_pm_one(shapes[s][1], 0)});
    gauges.push_back(LegGauge{3, ss_generic(shapes[s][3], 1)});

    ftensor<tensor> table = a;
    ftensor<tensor> direct = a;
    ftensor<tensor> automatic = a;
    tenes::fermion::apply_swap_form(table, form, SignEval::table);
    tenes::fermion::apply_swap_form(direct, form, SignEval::direct);
    tenes::fermion::apply_swap_form(automatic, form, SignEval::automatic);
    CHECK(ss_count_diff(table.t, direct.t) == 0);
    CHECK(ss_count_diff(automatic.t, table.t) == 0);
    CHECK(ss_count_diff(automatic.t, direct.t) == 0);

    ftensor<tensor> sweep_table = a;
    ftensor<tensor> sweep_direct = a;
    ftensor<tensor> sweep_auto = a;
    tenes::fermion::apply_sign_sweep(sweep_table, form, gauges,
                                     SignEval::table);
    tenes::fermion::apply_sign_sweep(sweep_direct, form, gauges,
                                     SignEval::direct);
    tenes::fermion::apply_sign_sweep(sweep_auto, form, gauges,
                                     SignEval::automatic);
    CHECK(ss_count_diff(sweep_table.t, sweep_direct.t) == 0);
    CHECK(ss_count_diff(sweep_auto.t, sweep_table.t) == 0);
    CHECK(ss_count_diff(sweep_auto.t, sweep_direct.t) == 0);

    // and the default argument is the automatic path
    ftensor<tensor> defaulted = a;
    tenes::fermion::apply_sign_sweep(defaulted, form, gauges);
    CHECK(ss_count_diff(defaulted.t, sweep_auto.t) == 0);

    ftensor<tensor> defaulted_form = a;
    tenes::fermion::apply_swap_form(defaulted_form, form);
    CHECK(ss_count_diff(defaulted_form.t, automatic.t) == 0);

    // neither path may be a silent no-op
    REQUIRE(ss_count_diff(table.t, a.t) > 0);
    REQUIRE(ss_count_diff(sweep_table.t, a.t) > 0);
  }
}

// ------------------------------------------------------------------ C1-7

TEST_CASE("SS C1-7 use_sign_table boundaries") {
  using tenes::fermion::detail::kMaxTableRank;
  using tenes::fermion::detail::use_sign_table;

  const std::size_t huge = std::numeric_limits<std::size_t>::max();

  SUBCASE("above kMaxTableRank the answer is false for any local size") {
    const std::size_t ranks[6] = {
        kMaxTableRank + 1, kMaxTableRank + 2, 32, 33, 64, 1000};
    const std::size_t sizes[5] = {0, 1, 1000, 1u << 20, huge};
    for (int r = 0; r < 6; ++r) {
      for (int s = 0; s < 5; ++s) {
        INFO("rank = " << ranks[r] << " local_size = " << sizes[s]);
        CHECK_FALSE(use_sign_table(ranks[r], sizes[s]));
      }
    }
  }

  SUBCASE("at or below kMaxTableRank the rule is 2^rank <= local_size") {
    for (std::size_t rank = 0; rank <= kMaxTableRank; ++rank) {
      const std::size_t entries = std::size_t(1) << rank;
      INFO("rank = " << rank);
      CHECK(use_sign_table(rank, entries));
      CHECK(use_sign_table(rank, entries + 1));
      if (entries > 0) {
        CHECK_FALSE(use_sign_table(rank, entries - 1));
      }
      CHECK(use_sign_table(rank, huge));
      CHECK_FALSE(use_sign_table(rank, 0));
    }
  }

  SUBCASE("the boundary right at kMaxTableRank") {
    const std::size_t entries = std::size_t(1) << kMaxTableRank;
    CHECK(use_sign_table(kMaxTableRank, entries));
    CHECK_FALSE(use_sign_table(kMaxTableRank, entries - 1));
    CHECK_FALSE(use_sign_table(kMaxTableRank + 1, entries));
    CHECK_FALSE(use_sign_table(kMaxTableRank + 1, huge));
  }

  SUBCASE("small ranks") {
    CHECK_FALSE(use_sign_table(0, 0));
    CHECK(use_sign_table(0, 1));
    CHECK_FALSE(use_sign_table(1, 1));
    CHECK(use_sign_table(1, 2));
    CHECK_FALSE(use_sign_table(3, 7));
    CHECK(use_sign_table(3, 8));
    CHECK(use_sign_table(3, 9));
    CHECK_FALSE(use_sign_table(16, 65535));
    CHECK(use_sign_table(16, 65536));
  }
}

// ------------------------------------------------------------------ C1-8

TEST_CASE_TEMPLATE("SS C1-8 leg gauges match repeated multiply_vector", tensor,
                   tenes::real_tensor, tenes::complex_tensor) {
  using namespace fermion_sign_sweep;

  const mptensor::Shape sh(2, 3, 4, 2, 3);
  const tensor base = ss_random_tensor<tensor>(sh, 8080u);
  const ftensor<tensor> fbase{base, ss_parities(sh)};

  SUBCASE("plus/minus one factors on three distinct axes") {
    std::vector<LegGauge> gauges;
    gauges.push_back(LegGauge{0, ss_pm_one(sh[0], 0)});
    gauges.push_back(LegGauge{2, ss_pm_one(sh[2], 1)});
    gauges.push_back(LegGauge{4, ss_pm_one(sh[4], 0)});

    const tensor expected = ss_sequential_gauges(base, gauges);
    REQUIRE(ss_count_diff(expected, base) > 0);

    tensor got = base;
    tenes::fermion::apply_leg_gauges(got, gauges);
    CHECK(ss_count_diff(got, expected) == 0);

    ftensor<tensor> swept = fbase;
    tenes::fermion::apply_sign_sweep(swept, SwapForm{}, gauges);
    CHECK(ss_count_diff(swept.t, expected) == 0);
    CHECK(swept.parity == fbase.parity);
  }

  SUBCASE("generic factors that are not powers of two") {
    // 0.1, 3.7, -1.3, ... : the products below are inexact, so this only
    // matches the reference if the kernel multiplies one factor at a time in
    // the order the gauges are listed.
    std::vector<LegGauge> gauges;
    gauges.push_back(LegGauge{1, ss_generic(sh[1], 0)});
    gauges.push_back(LegGauge{3, ss_generic(sh[3], 2)});
    gauges.push_back(LegGauge{4, ss_generic(sh[4], 3)});

    const tensor expected = ss_sequential_gauges(base, gauges);
    REQUIRE(ss_count_diff(expected, base) > 0);

    tensor got = base;
    tenes::fermion::apply_leg_gauges(got, gauges);
    CHECK(ss_count_diff(got, expected) == 0);

    ftensor<tensor> swept = fbase;
    tenes::fermion::apply_sign_sweep(swept, SwapForm{}, gauges);
    CHECK(ss_count_diff(swept.t, expected) == 0);
  }

  SUBCASE("two gauges on the same axis") {
    // Both factors of axis 2 are generic, so folding them into one vector
    // before the sweep changes the last bits and fails this check.
    std::vector<LegGauge> gauges;
    gauges.push_back(LegGauge{2, ss_generic(sh[2], 0)});
    gauges.push_back(LegGauge{2, ss_generic(sh[2], 3)});
    gauges.push_back(LegGauge{0, ss_generic(sh[0], 1)});

    const tensor expected = ss_sequential_gauges(base, gauges);
    REQUIRE(ss_count_diff(expected, base) > 0);

    tensor got = base;
    tenes::fermion::apply_leg_gauges(got, gauges);
    CHECK(ss_count_diff(got, expected) == 0);

    ftensor<tensor> swept = fbase;
    tenes::fermion::apply_sign_sweep(swept, SwapForm{}, gauges);
    CHECK(ss_count_diff(swept.t, expected) == 0);
  }

  SUBCASE("an empty gauge list is the identity") {
    const std::vector<LegGauge> none;
    tensor got = base;
    tenes::fermion::apply_leg_gauges(got, none);
    CHECK(ss_count_diff(got, base) == 0);
  }

  SUBCASE("swaps and gauges combine in one sweep") {
    const pair_list pairs{{0, 2}, {1, 4}, {3, 0}};
    std::vector<LegGauge> gauges;
    gauges.push_back(LegGauge{1, ss_generic(sh[1], 1)});
    gauges.push_back(LegGauge{4, ss_pm_one(sh[4], 0)});
    gauges.push_back(LegGauge{2, ss_generic(sh[2], 4)});

    ftensor<tensor> reference = ss_sequential_swaps(fbase, pairs);
    reference.t = ss_sequential_gauges(reference.t, gauges);
    REQUIRE(ss_count_diff(reference.t, base) > 0);

    ftensor<tensor> got = fbase;
    tenes::fermion::apply_sign_sweep(got, ss_form(pairs), gauges);
    CHECK(ss_count_diff(got.t, reference.t) == 0);
    CHECK(got.parity == fbase.parity);
  }
}

// ------------------------------------------------------------------ C1-9

TEST_CASE_TEMPLATE("SS C1-9 apply_swap still negates odd-odd elements only",
                   tensor, tenes::real_tensor, tenes::complex_tensor) {
  using namespace fermion_sign_sweep;

  const mptensor::Shape sh(2, 3, 4, 2, 3);
  const ftensor<tensor> a = ss_random_ftensor<tensor>(sh, 606u);

  const int axis_pairs[4][2] = {{0, 2}, {1, 4}, {2, 3}, {0, 4}};
  for (int k = 0; k < 4; ++k) {
    const int ax1 = axis_pairs[k][0];
    const int ax2 = axis_pairs[k][1];
    INFO("apply_swap(" << ax1 << ", " << ax2 << ")");

    ftensor<tensor> b = a;
    tenes::fermion::apply_swap(b, ax1, ax2);

    std::size_t bad = 0;
    std::size_t flipped = 0;
    for (std::size_t n = 0; n < a.t.local_size(); ++n) {
      const mptensor::Index idx = a.t.global_index(n);
      typename tensor::value_type va, vb;
      a.t.get_value(idx, va);
      b.t.get_value(idx, vb);
      const bool odd_odd = a.parity[ax1][idx[ax1]] && a.parity[ax2][idx[ax2]];
      const typename tensor::value_type expected = odd_odd ? -va : va;
      if (!(vb == expected)) {
        ++bad;
      }
      if (odd_odd) {
        ++flipped;
      }
    }
    CHECK(bad == 0);
    REQUIRE(flipped > 0);
    CHECK(b.parity == a.parity);

    // swapping the two arguments must not change anything
    ftensor<tensor> c = a;
    tenes::fermion::apply_swap(c, ax2, ax1);
    CHECK(ss_count_diff(c.t, b.t) == 0);

    // applying it twice is the identity
    ftensor<tensor> d = b;
    tenes::fermion::apply_swap(d, ax1, ax2);
    CHECK(ss_count_diff(d.t, a.t) == 0);
  }

  // Ruling R6: apply_swap must stay bit-identical to the current
  // implementation for EVERY input, the equal-axis case included.  With
  // ax1 == ax2 the two parity conditions coincide, so every element whose
  // parity on that axis is odd flips -- exactly what apply_parity does.
  // (SwapForm::toggle(x, x) is a no-op by contract, so apply_swap must NOT
  // route the equal-axis case through it.)
  for (int ax = 0; ax < a.rank(); ++ax) {
    INFO("apply_swap(" << ax << ", " << ax << ")");

    ftensor<tensor> b = a;
    tenes::fermion::apply_swap(b, ax, ax);

    std::size_t bad = 0;
    std::size_t flipped = 0;
    for (std::size_t n = 0; n < a.t.local_size(); ++n) {
      const mptensor::Index idx = a.t.global_index(n);
      typename tensor::value_type va, vb;
      a.t.get_value(idx, va);
      b.t.get_value(idx, vb);
      const bool odd = a.parity[ax][idx[ax]];
      const typename tensor::value_type expected = odd ? -va : va;
      if (!(vb == expected)) {
        ++bad;
      }
      if (odd) {
        ++flipped;
      }
    }
    CHECK(bad == 0);
    REQUIRE(flipped > 0);
    CHECK(b.parity == a.parity);

    // ... which is exactly apply_parity on that axis
    ftensor<tensor> viaparity = a;
    tenes::fermion::apply_parity(viaparity, ax);
    CHECK(ss_count_diff(b.t, viaparity.t) == 0);

    // and it is still an involution
    ftensor<tensor> twice = b;
    tenes::fermion::apply_swap(twice, ax, ax);
    CHECK(ss_count_diff(twice.t, a.t) == 0);
  }
}

// ===== SS C2: graded transpose folded into the fused sign sweep ===========
//
// `ss_reference` below holds VERBATIM COPIES of the pre-rewrite code, taken
// from the tag wip/fermion_20260822:
//
//   transpose_sign            <- src/fermion/ftensor.hpp,
//   detail::transpose_sign transpose_inplace         <-
//   src/fermion/ftensor.hpp, ftensor::transpose
//                                (member turned into a free function; the body
//                                 is character for character the original)
//   transposed                <- src/fermion/fops.hpp, fermion::transpose
//   apply_transpose_sign_mask <- src/fermion/fops.hpp,
//                                detail::apply_transpose_sign_mask
//   tensordot                 <- src/fermion/fops.hpp, fermion::tensordot,
//                                with the mask above substituted
//   svd_values                <- src/fermion/fops.hpp, fermion::svd, truncated
//                                after the two block SVDs (the contract only
//                                asks for the singular values), with the
//                                reference transpose above substituted
//
// Only namespace qualification was added. The task-2 contract requires this:
// the rewrite changes src/, so src/ cannot serve as its own reference.

namespace ss_reference {

using tenes::fermion::ftensor;
using tenes::fermion::leg_parities;
using tenes::fermion::parity_vector;

// verbatim: tenes::fermion::detail::transpose_sign (ftensor.hpp)
inline int transpose_sign(const leg_parities& parity,
                          const mptensor::Index& idx,
                          const mptensor::Axes& axes) {
  int k = 0;
  for (std::size_t x = 0; x < axes.size(); ++x) {
    for (std::size_t y = x + 1; y < axes.size(); ++y) {
      if (axes[x] > axes[y] && parity[axes[y]][idx[axes[y]]] &&
          parity[axes[x]][idx[axes[x]]]) {
        ++k;
      }
    }
  }
  return (k % 2 == 0) ? 1 : -1;
}

// verbatim: tenes::fermion::ftensor<tensor>::transpose (ftensor.hpp)
template <class tensor>
ftensor<tensor>& transpose_inplace(ftensor<tensor>& self,
                                   const mptensor::Axes& axes) {
  mptensor::Index idx;
  idx.resize(self.t.shape().size());
  for (std::size_t n = 0; n < self.t.local_size(); ++n) {
    self.t.global_index_fast(n, idx);
    self.t[n] *= transpose_sign(self.parity, idx, axes);
  }
  self.t.transpose(axes);
  leg_parities next;
  next.reserve(axes.size());
  for (std::size_t i = 0; i < axes.size(); ++i) {
    next.push_back(self.parity[axes[i]]);
  }
  self.parity = next;
  return self;
}

// verbatim: tenes::fermion::transpose (fops.hpp)
template <class tensor>
ftensor<tensor> transposed(const ftensor<tensor>& a,
                           const mptensor::Axes& axes) {
  ftensor<tensor> ret = a;
  transpose_inplace(ret, axes);
  return ret;
}

// verbatim: tenes::fermion::detail::apply_transpose_sign_mask (fops.hpp)
template <class tensor>
ftensor<tensor> apply_transpose_sign_mask(const ftensor<tensor>& a,
                                          const mptensor::Axes& axes) {
  ftensor<tensor> ret = a;
  mptensor::Index idx;
  idx.resize(ret.t.shape().size());
  for (std::size_t n = 0; n < ret.t.local_size(); ++n) {
    ret.t.global_index_fast(n, idx);
    ret.t[n] *= transpose_sign(ret.parity, idx, axes);
  }
  return ret;
}

// verbatim: tenes::fermion::tensordot (fops.hpp), reference mask substituted
template <class tensor>
ftensor<tensor> tensordot(const ftensor<tensor>& a, const ftensor<tensor>& b,
                          const mptensor::Axes& axes_a,
                          const mptensor::Axes& axes_b) {
  namespace fdetail = tenes::fermion::detail;
  fdetail::validate_contracted_parity(a.parity, b.parity, axes_a, axes_b);
  ftensor<tensor> a_masked = apply_transpose_sign_mask(
      a, fdetail::tensordot_left_perm(a.parity.size(), axes_a));
  ftensor<tensor> b_masked = apply_transpose_sign_mask(
      b, fdetail::tensordot_right_perm(b.parity.size(), axes_b));
  ftensor<tensor> ret;
  ret.t = mptensor::tensordot(a_masked.t, b_masked.t, axes_a, axes_b);
  ret.parity = fdetail::free_leg_parities(a.parity, axes_a);
  leg_parities b_free = fdetail::free_leg_parities(b.parity, axes_b);
  ret.parity.insert(ret.parity.end(), b_free.begin(), b_free.end());
  return ret;
}

// verbatim: the singular-value half of tenes::fermion::svd (fops.hpp), with
// the reference transpose substituted for the production one.
template <class tensor>
int svd_values(const ftensor<tensor>& a, const mptensor::Axes& rows,
               const mptensor::Axes& cols, std::vector<double>& s) {
  namespace fdetail = tenes::fermion::detail;
  ftensor<tensor> a_ordered = transposed(a, rows + cols);
  mptensor::Shape row_shape = fdetail::shape_from_axes(a.shape(), rows);
  mptensor::Shape col_shape = fdetail::shape_from_axes(a.shape(), cols);
  const std::size_t drow = fdetail::product_shape(row_shape);
  const std::size_t dcol = fdetail::product_shape(col_shape);
  tensor mat = mptensor::reshape(a_ordered.t, mptensor::Shape(drow, dcol));

  parity_vector row_parity = fdetail::fuse_axes(a.parity, rows);
  parity_vector col_parity = fdetail::fuse_axes(a.parity, cols);
  std::vector<std::size_t> row_perm =
      tenes::fermion::parity_sort_perm(row_parity);
  std::vector<std::size_t> col_perm =
      tenes::fermion::parity_sort_perm(col_parity);
  tensor prow = tenes::fermion::make_perm_matrix<tensor>(row_perm);
  tensor pcol = tenes::fermion::make_perm_matrix<tensor>(col_perm);
  tensor sorted =
      mptensor::tensordot(prow, mat, mptensor::Axes(1), mptensor::Axes(0));
  sorted =
      mptensor::tensordot(sorted, pcol, mptensor::Axes(1), mptensor::Axes(1));

  const std::size_t row_even = fdetail::count_even(row_parity);
  const std::size_t col_even = fdetail::count_even(col_parity);
  fdetail::validate_block_diagonal(sorted, row_even, col_even,
                                   "ss_reference svd");
  const std::size_t row_odd = drow - row_even;
  const std::size_t col_odd = dcol - col_even;
  const std::size_t size_even = std::min(row_even, col_even);
  const std::size_t size_odd = std::min(row_odd, col_odd);

  std::vector<double> s_even, s_odd;
  int info = 0;
  if (size_even > 0) {
    tensor even_block = mptensor::slice(sorted, mptensor::Index(0, 0),
                                        mptensor::Index(row_even, col_even));
    tensor ue, vte;
    info = mptensor::svd(even_block, mptensor::Axes(0), mptensor::Axes(1), ue,
                         s_even, vte);
  }
  if (size_odd > 0) {
    tensor odd_block =
        mptensor::slice(sorted, mptensor::Index(row_even, col_even),
                        mptensor::Index(drow, dcol));
    tensor uo, vto;
    int info_odd = mptensor::svd(odd_block, mptensor::Axes(0),
                                 mptensor::Axes(1), uo, s_odd, vto);
    if (info == 0) {
      info = info_odd;
    }
  }
  s = s_even;
  s.insert(s.end(), s_odd.begin(), s_odd.end());
  return info;
}

}  // namespace ss_reference

// ------------------------------------------------------- C2 test helpers

namespace fermion_sign_sweep {

struct ss_named_axes {
  std::string name;
  mptensor::Axes axes;
};

inline std::string ss_axes_string(const mptensor::Axes& axes) {
  std::ostringstream os;
  os << "[";
  for (std::size_t i = 0; i < axes.size(); ++i) {
    if (i > 0) {
      os << ",";
    }
    os << axes[i];
  }
  os << "]";
  return os.str();
}

inline mptensor::Axes ss_identity_axes(std::size_t rank) {
  std::vector<std::size_t> p(rank);
  for (std::size_t i = 0; i < rank; ++i) {
    p[i] = i;
  }
  return mptensor::Axes(p);
}

// identity, full reverse, two cyclic shifts, two adjacent swaps, the end
// swap (0, rank-1) -- which always exchanges two odd-carrying legs, so the
// inversion count of the odd legs is never trivially zero -- and three
// deterministic random permutations.
inline std::vector<ss_named_axes> ss_axes_catalogue(std::size_t rank) {
  std::vector<ss_named_axes> out;
  std::vector<std::size_t> p(rank);
  for (std::size_t i = 0; i < rank; ++i) {
    p[i] = i;
  }
  out.push_back({"identity", mptensor::Axes(p)});
  {
    std::vector<std::size_t> q = p;
    std::reverse(q.begin(), q.end());
    out.push_back({"reverse", mptensor::Axes(q)});
  }
  for (std::size_t k = 1; k <= 2 && k < rank; ++k) {
    std::vector<std::size_t> q(rank);
    for (std::size_t i = 0; i < rank; ++i) {
      q[i] = (i + k) % rank;
    }
    out.push_back({"cyclic+" + std::to_string(k), mptensor::Axes(q)});
  }
  {
    std::vector<std::size_t> q = p;
    std::swap(q[0], q[1]);
    out.push_back({"swap(0,1)", mptensor::Axes(q)});
  }
  if (rank >= 3) {
    std::vector<std::size_t> q = p;
    std::swap(q[rank - 2], q[rank - 1]);
    out.push_back({"swap(rank-2,rank-1)", mptensor::Axes(q)});
  }
  {
    std::vector<std::size_t> q = p;
    std::swap(q[0], q[rank - 1]);
    out.push_back({"swap(0,rank-1)", mptensor::Axes(q)});
  }
  for (unsigned seed = 0; seed < 3; ++seed) {
    std::vector<std::size_t> q = p;
    std::mt19937 gen(1000u * static_cast<unsigned>(rank) + seed);
    std::shuffle(q.begin(), q.end(), gen);
    out.push_back({"random#" + std::to_string(seed), mptensor::Axes(q)});
  }
  return out;
}

// rank 2 .. 10, no two axes of a shape share the same dimension pattern and
// the dimensions are never all equal.
inline std::vector<std::vector<std::size_t>> ss_c2_shapes() {
  return {{3, 4},
          {2, 3, 4},
          {2, 3, 4, 2},
          {2, 3, 2, 4, 3},
          {2, 3, 4, 2, 3, 2},
          {3, 2, 2, 4, 2, 3, 2},
          {2, 3, 2, 2, 3, 2, 2, 2},
          {2, 2, 3, 2, 2, 2, 3, 2, 2},
          {2, 2, 3, 2, 2, 2, 2, 3, 2, 2}};
}

// How many elements the reference transpose sign turns negative.  Used as a
// hollowness guard: a permutation that never flips a sign proves nothing.
template <class tensor>
std::size_t ss_flipped_count(const tenes::fermion::ftensor<tensor>& a,
                             const mptensor::Axes& axes) {
  std::size_t n_neg = 0;
  mptensor::Index idx;
  idx.resize(a.t.shape().size());
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    a.t.global_index_fast(n, idx);
    if (ss_reference::transpose_sign(a.parity, idx, axes) < 0) {
      ++n_neg;
    }
  }
  return n_neg;
}

// Same, but only counting elements that are actually non-zero: a parity-even
// tensor is mostly zeros, and a sign flipped on a zero is invisible.
template <class tensor>
std::size_t ss_flipped_nonzero(const tenes::fermion::ftensor<tensor>& a,
                               const mptensor::Axes& axes) {
  std::size_t n_neg = 0;
  mptensor::Index idx;
  idx.resize(a.t.shape().size());
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    a.t.global_index_fast(n, idx);
    if (std::abs(a.t[n]) != 0.0 &&
        ss_reference::transpose_sign(a.parity, idx, axes) < 0) {
      ++n_neg;
    }
  }
  return n_neg;
}

// Parity-even random tensor, built the same way as make_r2_tensor() in
// test/fermion/r2_convention.cpp: only elements whose total parity is even
// are filled, the rest stay exactly zero.  fermion::svd / fermion::qr run
// detail::validate_block_diagonal, which throws in debug builds on an input
// that is not parity even, so a plain random tensor cannot be used there.
template <class tensor>
tenes::fermion::ftensor<tensor> ss_random_even_ftensor(
    const mptensor::Shape& sh, unsigned seed) {
  using value_type = typename tensor::value_type;
  const tenes::fermion::leg_parities p = ss_parities(sh);
  tensor t(sh);
  std::mt19937 gen(seed);
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    if (tenes::fermion::count_odd(p, idx) % 2 == 0) {
      t.set_value(idx, ss_sampler<value_type>::get(gen));
    }
  }
  return tenes::fermion::ftensor<tensor>{t, p};
}

template <class tensor>
std::size_t ss_count_nonzero(const tenes::fermion::ftensor<tensor>& a) {
  std::size_t n_nz = 0;
  for (std::size_t n = 0; n < a.t.local_size(); ++n) {
    if (std::abs(a.t[n]) != 0.0) {
      ++n_nz;
    }
  }
  return n_nz;
}

}  // namespace fermion_sign_sweep

// ------------------------------------------------------------------ C2-1

TEST_CASE_TEMPLATE(
    "SS C2-1 transpose_with_swap_form reproduces the reference transpose",
    tensor, tenes::real_tensor, tenes::complex_tensor) {
  using namespace fermion_sign_sweep;

  const std::vector<std::vector<std::size_t>> shapes = ss_c2_shapes();
  for (std::size_t s = 0; s < shapes.size(); ++s) {
    const mptensor::Shape sh(shapes[s]);
    const tenes::fermion::ftensor<tensor> a =
        ss_random_ftensor<tensor>(sh, 5000u + 17u * static_cast<unsigned>(s));

    for (const ss_named_axes& na : ss_axes_catalogue(sh.size())) {
      INFO("rank " << sh.size() << " " << na.name << " "
                   << ss_axes_string(na.axes));

      const tenes::fermion::ftensor<tensor> expected =
          ss_reference::transposed(a, na.axes);

      tenes::fermion::ftensor<tensor> got = a;
      tenes::fermion::transpose_with_swap_form(got, tenes::fermion::SwapForm{},
                                               na.axes);

      REQUIRE(got.t.shape() == expected.t.shape());
      CHECK(ss_count_diff(got.t, expected.t) == 0);
      CHECK(got.parity == expected.parity);
    }

    // The end swap always exchanges two legs that both carry odd entries, so
    // the odd-leg inversion count is genuinely non-zero: without this the
    // whole comparison above could hold with every sign equal to +1.
    std::vector<std::size_t> ends(sh.size());
    for (std::size_t i = 0; i < sh.size(); ++i) {
      ends[i] = i;
    }
    std::swap(ends[0], ends[sh.size() - 1]);
    INFO("rank " << sh.size() << " swap(0,rank-1)");
    REQUIRE(ss_flipped_count(a, mptensor::Axes(ends)) > 0);
  }
}

// ------------------------------------------------------------------ C2-2

TEST_CASE_TEMPLATE(
    "SS C2-2 transpose_with_swap_form composes the swaps with the transpose",
    tensor, tenes::real_tensor, tenes::complex_tensor) {
  using namespace fermion_sign_sweep;

  struct named_form {
    const char* name;
    pair_list pairs;
  };
  const std::vector<named_form> forms = {
      {"single pair", {{0, 1}}},
      {"two disjoint pairs", {{0, 2}, {1, 3}}},
      {"three pairs", {{0, 3}, {1, 2}, {0, 1}}},
      // (0,2)/(2,0) and (1,3)/(3,1) cancel, only (0,3) survives
      {"with duplicate pairs", {{0, 2}, {1, 3}, {2, 0}, {3, 1}, {0, 3}}},
      // everything cancels: the form is empty again
      {"fully cancelling", {{0, 2}, {1, 3}, {2, 0}, {3, 1}}},
  };

  const std::vector<std::vector<std::size_t>> all = ss_c2_shapes();
  // ranks 4, 5, 6 and 8 (every form above only touches axes 0..3)
  const std::size_t picked[4] = {2, 3, 4, 6};

  for (std::size_t k = 0; k < 4; ++k) {
    const mptensor::Shape sh(all[picked[k]]);
    const tenes::fermion::ftensor<tensor> a =
        ss_random_ftensor<tensor>(sh, 7100u + 29u * static_cast<unsigned>(k));

    for (const ss_named_axes& na : ss_axes_catalogue(sh.size())) {
      for (const named_form& nf : forms) {
        INFO("rank " << sh.size() << " " << na.name << " "
                     << ss_axes_string(na.axes) << " form " << nf.name);

        const tenes::fermion::SwapForm form = ss_form(nf.pairs);

        // reference: swap signs first, then the pre-rewrite transpose
        tenes::fermion::ftensor<tensor> swapped = a;
        tenes::fermion::apply_swap_form(swapped, form);
        const tenes::fermion::ftensor<tensor> expected =
            ss_reference::transposed(swapped, na.axes);

        tenes::fermion::ftensor<tensor> got = a;
        tenes::fermion::transpose_with_swap_form(got, form, na.axes);

        REQUIRE(got.t.shape() == expected.t.shape());
        CHECK(ss_count_diff(got.t, expected.t) == 0);
        CHECK(got.parity == expected.parity);
      }
    }

    // the swap part alone must not be a no-op
    tenes::fermion::ftensor<tensor> swapped_only = a;
    tenes::fermion::apply_swap_form(swapped_only, ss_form(forms[1].pairs));
    INFO("rank " << sh.size());
    REQUIRE(ss_count_diff(swapped_only.t, a.t) > 0);
  }
}

// ------------------------------------------------------------------ C2-3

TEST_CASE_TEMPLATE("SS C2-3 transpose_with_swap_form table and direct agree",
                   tensor, tenes::real_tensor, tenes::complex_tensor) {
  using namespace fermion_sign_sweep;

  // shape 0: local_size 144 >= 2^5, so `automatic` may build a table
  // shape 1: local_size 32 < 2^6 (one axis of dimension 1), so `automatic`
  //          has to evaluate directly.  Both sides of the boundary are hit.
  const mptensor::Shape shapes[2] = {mptensor::Shape(2, 3, 2, 4, 3),
                                     mptensor::Shape(2, 2, 2, 2, 2, 1)};
  const pair_list pairs[2] = {{{0, 2}, {1, 4}, {3, 0}},
                              {{0, 5}, {1, 3}, {2, 4}}};

  for (int s = 0; s < 2; ++s) {
    const mptensor::Shape& sh = shapes[s];
    const tenes::fermion::ftensor<tensor> a =
        ss_random_ftensor<tensor>(sh, 8200u + 13u * static_cast<unsigned>(s));
    const tenes::fermion::SwapForm form = ss_form(pairs[s]);

    for (const ss_named_axes& na : ss_axes_catalogue(sh.size())) {
      INFO("shape " << s << " " << na.name << " " << ss_axes_string(na.axes));

      tenes::fermion::ftensor<tensor> table = a;
      tenes::fermion::ftensor<tensor> direct = a;
      tenes::fermion::ftensor<tensor> automatic = a;
      tenes::fermion::transpose_with_swap_form(table, form, na.axes,
                                               tenes::fermion::SignEval::table);
      tenes::fermion::transpose_with_swap_form(
          direct, form, na.axes, tenes::fermion::SignEval::direct);
      tenes::fermion::transpose_with_swap_form(
          automatic, form, na.axes, tenes::fermion::SignEval::automatic);

      CHECK(ss_count_diff(table.t, direct.t) == 0);
      CHECK(ss_count_diff(automatic.t, table.t) == 0);
      CHECK(ss_count_diff(automatic.t, direct.t) == 0);
      CHECK(table.parity == direct.parity);
      CHECK(automatic.parity == direct.parity);

      // and the default argument is the automatic path
      tenes::fermion::ftensor<tensor> defaulted = a;
      tenes::fermion::transpose_with_swap_form(defaulted, form, na.axes);
      CHECK(ss_count_diff(defaulted.t, automatic.t) == 0);

      // both paths must agree with the pre-rewrite reference as well
      tenes::fermion::ftensor<tensor> swapped = a;
      tenes::fermion::apply_swap_form(swapped, form);
      const tenes::fermion::ftensor<tensor> expected =
          ss_reference::transposed(swapped, na.axes);
      CHECK(ss_count_diff(table.t, expected.t) == 0);
      CHECK(ss_count_diff(direct.t, expected.t) == 0);
    }
  }
}

// ------------------------------------------------------------------ C2-4

TEST_CASE_TEMPLATE(
    "SS C2-4 transpose_with_swap_form short-circuits the identity permutation",
    tensor, tenes::real_tensor, tenes::complex_tensor) {
  using namespace fermion_sign_sweep;

  const std::vector<std::vector<std::size_t>> shapes = ss_c2_shapes();
  for (std::size_t s = 0; s < shapes.size(); ++s) {
    const mptensor::Shape sh(shapes[s]);
    const mptensor::Axes id = ss_identity_axes(sh.size());
    const tenes::fermion::ftensor<tensor> a =
        ss_random_ftensor<tensor>(sh, 9300u + 31u * static_cast<unsigned>(s));
    INFO("rank " << sh.size());

    SUBCASE("empty form and identity leave the tensor alone") {
      tenes::fermion::ftensor<tensor> got = a;
      tenes::fermion::transpose_with_swap_form(got, tenes::fermion::SwapForm{},
                                               id);
      CHECK(ss_count_diff(got.t, a.t) == 0);
      CHECK(got.parity == a.parity);
      CHECK(got.t.shape() == a.t.shape());
    }

    SUBCASE("identity applies the swap signs only") {
      pair_list pairs;
      pairs.push_back({0, static_cast<int>(sh.size()) - 1});
      if (sh.size() >= 4) {
        pairs.push_back({1, 2});
      }
      const tenes::fermion::SwapForm form = ss_form(pairs);

      tenes::fermion::ftensor<tensor> expected = a;
      tenes::fermion::apply_swap_form(expected, form);
      REQUIRE(ss_count_diff(expected.t, a.t) > 0);

      tenes::fermion::ftensor<tensor> got = a;
      tenes::fermion::transpose_with_swap_form(got, form, id);
      CHECK(ss_count_diff(got.t, expected.t) == 0);
      CHECK(got.parity == a.parity);
    }
  }
}

// ------------------------------------------------------------------ C2-5

TEST_CASE_TEMPLATE("SS C2-5 ftensor::transpose still matches the reference",
                   tensor, tenes::real_tensor, tenes::complex_tensor) {
  using namespace fermion_sign_sweep;

  const std::vector<std::vector<std::size_t>> shapes = ss_c2_shapes();
  for (std::size_t s = 0; s < shapes.size(); ++s) {
    const mptensor::Shape sh(shapes[s]);
    const tenes::fermion::ftensor<tensor> a =
        ss_random_ftensor<tensor>(sh, 5000u + 17u * static_cast<unsigned>(s));

    for (const ss_named_axes& na : ss_axes_catalogue(sh.size())) {
      INFO("rank " << sh.size() << " " << na.name << " "
                   << ss_axes_string(na.axes));

      const tenes::fermion::ftensor<tensor> expected =
          ss_reference::transposed(a, na.axes);

      tenes::fermion::ftensor<tensor> got = a;
      got.transpose(na.axes);

      REQUIRE(got.t.shape() == expected.t.shape());
      CHECK(ss_count_diff(got.t, expected.t) == 0);
      CHECK(got.parity == expected.parity);

      // the free function in fops.hpp goes through the same member
      const tenes::fermion::ftensor<tensor> viafree =
          tenes::fermion::transpose(a, na.axes);
      CHECK(ss_count_diff(viafree.t, expected.t) == 0);
      CHECK(viafree.parity == expected.parity);
    }
  }
}

// ------------------------------------------------------------------ C2-6

TEST_CASE_TEMPLATE(
    "SS C2-6 apply_transpose_sign_mask short-circuits the identity permutation",
    tensor, tenes::real_tensor, tenes::complex_tensor) {
  using namespace fermion_sign_sweep;

  const std::vector<std::vector<std::size_t>> shapes = ss_c2_shapes();
  for (std::size_t s = 0; s < shapes.size(); ++s) {
    const mptensor::Shape sh(shapes[s]);
    const tenes::fermion::ftensor<tensor> a =
        ss_random_ftensor<tensor>(sh, 6400u + 23u * static_cast<unsigned>(s));

    const mptensor::Axes id = ss_identity_axes(sh.size());
    INFO("rank " << sh.size() << " identity");
    const tenes::fermion::ftensor<tensor> unchanged =
        tenes::fermion::detail::apply_transpose_sign_mask(a, id);
    CHECK(ss_count_diff(unchanged.t, a.t) == 0);
    CHECK(unchanged.parity == a.parity);
    CHECK(unchanged.t.shape() == a.t.shape());

    for (const ss_named_axes& na : ss_axes_catalogue(sh.size())) {
      INFO("rank " << sh.size() << " " << na.name << " "
                   << ss_axes_string(na.axes));
      const tenes::fermion::ftensor<tensor> expected =
          ss_reference::apply_transpose_sign_mask(a, na.axes);
      const tenes::fermion::ftensor<tensor> got =
          tenes::fermion::detail::apply_transpose_sign_mask(a, na.axes);
      // the mask only multiplies signs: neither the shape nor the ledger move
      REQUIRE(got.t.shape() == a.t.shape());
      CHECK(got.parity == a.parity);
      CHECK(ss_count_diff(got.t, expected.t) == 0);
    }

    // the permutations tensordot actually feeds to the mask
    for (std::size_t nc = 0; nc + 1 < sh.size(); ++nc) {
      mptensor::Axes contracted;
      for (std::size_t i = 0; i <= nc; ++i) {
        contracted.push(sh.size() - 1 - i);
      }
      const mptensor::Axes left =
          tenes::fermion::detail::tensordot_left_perm(sh.size(), contracted);
      const mptensor::Axes right =
          tenes::fermion::detail::tensordot_right_perm(sh.size(), contracted);
      INFO("rank " << sh.size() << " tensordot perms for "
                   << ss_axes_string(contracted));
      CHECK(ss_count_diff(
                tenes::fermion::detail::apply_transpose_sign_mask(a, left).t,
                ss_reference::apply_transpose_sign_mask(a, left).t) == 0);
      CHECK(ss_count_diff(
                tenes::fermion::detail::apply_transpose_sign_mask(a, right).t,
                ss_reference::apply_transpose_sign_mask(a, right).t) == 0);
    }
  }
}

// ------------------------------------------------------------------ C2-7

TEST_CASE_TEMPLATE("SS C2-7 fermion transpose, tensordot and svd are unchanged",
                   tensor, tenes::real_tensor, tenes::complex_tensor) {
  using namespace fermion_sign_sweep;

  SUBCASE("fermion::transpose") {
    const mptensor::Shape sh(2, 3, 4, 2, 3);
    const tenes::fermion::ftensor<tensor> a =
        ss_random_ftensor<tensor>(sh, 11100u);
    for (const ss_named_axes& na : ss_axes_catalogue(sh.size())) {
      INFO(na.name << " " << ss_axes_string(na.axes));
      const tenes::fermion::ftensor<tensor> expected =
          ss_reference::transposed(a, na.axes);
      const tenes::fermion::ftensor<tensor> got =
          tenes::fermion::transpose(a, na.axes);
      REQUIRE(got.t.shape() == expected.t.shape());
      CHECK(ss_count_diff(got.t, expected.t) == 0);
      CHECK(got.parity == expected.parity);
    }
  }

  SUBCASE("fermion::tensordot with no contracted axis (outer product)") {
    const tenes::fermion::ftensor<tensor> a =
        ss_random_ftensor<tensor>(mptensor::Shape(2, 3, 4), 12100u);
    const tenes::fermion::ftensor<tensor> b =
        ss_random_ftensor<tensor>(mptensor::Shape(3, 2), 12200u);

    const tenes::fermion::ftensor<tensor> expected =
        ss_reference::tensordot(a, b, mptensor::Axes(), mptensor::Axes());
    const tenes::fermion::ftensor<tensor> got =
        tenes::fermion::tensordot(a, b, mptensor::Axes(), mptensor::Axes());

    REQUIRE(got.t.shape() == expected.t.shape());
    REQUIRE(got.t.shape().size() == 5);
    CHECK(ss_count_diff(got.t, expected.t) == 0);
    CHECK(got.parity == expected.parity);
  }

  SUBCASE("fermion::tensordot with contracted axes") {
    const mptensor::Shape sh(2, 3, 4, 2);
    const tenes::fermion::ftensor<tensor> a =
        ss_random_ftensor<tensor>(sh, 13100u);

    // b shares a's ledger on the legs that get contracted, otherwise
    // validate_contracted_parity rejects the call.
    tenes::fermion::ftensor<tensor> b = ss_random_ftensor<tensor>(sh, 13200u);
    b.parity[0] = a.parity[0];
    b.parity[2] = a.parity[2];
    b.parity[3] = a.parity[3];
    b.parity[1] = ss_parity(13, sh[1]);

    struct named_pair {
      const char* name;
      mptensor::Axes axes_a;
      mptensor::Axes axes_b;
    };
    const std::vector<named_pair> cases = {
        {"one axis", mptensor::Axes(2), mptensor::Axes(2)},
        {"two axes", mptensor::Axes(0, 3), mptensor::Axes(0, 3)},
        {"two axes reversed", mptensor::Axes(3, 0), mptensor::Axes(3, 0)},
        {"three axes", mptensor::Axes(0, 2, 3), mptensor::Axes(0, 2, 3)},
    };

    for (const named_pair& nc : cases) {
      INFO(std::string(nc.name));
      const tenes::fermion::ftensor<tensor> expected =
          ss_reference::tensordot(a, b, nc.axes_a, nc.axes_b);
      const tenes::fermion::ftensor<tensor> got =
          tenes::fermion::tensordot(a, b, nc.axes_a, nc.axes_b);
      REQUIRE(got.t.shape() == expected.t.shape());
      CHECK(ss_count_diff(got.t, expected.t) == 0);
      CHECK(got.parity == expected.parity);

      // the masks tensordot applies must not be trivial for these cases
      REQUIRE(
          ss_flipped_count(a, tenes::fermion::detail::tensordot_left_perm(
                                  a.parity.size(), nc.axes_a)) +
              ss_flipped_count(b, tenes::fermion::detail::tensordot_right_perm(
                                      b.parity.size(), nc.axes_b)) >
          0);
    }
  }

  SUBCASE("fermion::svd singular values") {
    // parity even, otherwise validate_block_diagonal throws in debug builds
    const mptensor::Shape sh(2, 3, 2, 3);
    const tenes::fermion::ftensor<tensor> a =
        ss_random_even_ftensor<tensor>(sh, 14100u);
    REQUIRE(ss_count_nonzero(a) > 0);

    struct named_split {
      const char* name;
      mptensor::Axes rows;
      mptensor::Axes cols;
    };
    const std::vector<named_split> splits = {
        {"rows(0,2) cols(1,3)", mptensor::Axes(0, 2), mptensor::Axes(1, 3)},
        {"rows(2,0) cols(3,1)", mptensor::Axes(2, 0), mptensor::Axes(3, 1)},
        {"rows(1,3) cols(0,2)", mptensor::Axes(1, 3), mptensor::Axes(0, 2)},
        {"rows(3,1) cols(2,0)", mptensor::Axes(3, 1), mptensor::Axes(2, 0)},
        // rows+cols is the identity here: svd then exercises the
        // short-circuit of C2-4/C2-6 rather than a sign-carrying transpose.
        {"rows(0,1) cols(2,3)", mptensor::Axes(0, 1), mptensor::Axes(2, 3)},
    };

    // At least one of the splits must actually flip the sign of a non-zero
    // element, otherwise the whole subcase would hold with every sign +1.
    std::size_t sign_carrying = 0;
    for (const named_split& ns : splits) {
      INFO(std::string(ns.name));
      sign_carrying += ss_flipped_nonzero(a, ns.rows + ns.cols);

      std::vector<double> s_ref;
      const int info_ref = ss_reference::svd_values(a, ns.rows, ns.cols, s_ref);

      tenes::fermion::ftensor<tensor> u, vt;
      std::vector<double> s;
      const int info = tenes::fermion::svd(a, ns.rows, ns.cols, u, s, vt);

      CHECK(info == info_ref);
      REQUIRE(s_ref.size() > 0);
      REQUIRE(s.size() == s_ref.size());
      CHECK(s == s_ref);
    }
    REQUIRE(sign_carrying > 0);
  }
}
