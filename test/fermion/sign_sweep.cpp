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

#include <complex>
#include <cstddef>
#include <limits>
#include <random>
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
