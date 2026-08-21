// ===== MF: mean-field two-site measurement through the pair state ==========
//
// Pins build_pair_state / apply_pair_op / contract_pair_MF
// (src/fermion/reduced.hpp, src/fermion/reduced_measure.hpp) against
//   layer 1: the Fock oracle (test/fermion/fock_oracle.py, `mf_*` lines) on
//            d = 2 patches whose open legs carry both parities and a
//            mean-field weight lambda = {1.0, 0.7};
//   layer 2: the Fock-verified direct f-primitive path at d = 4
//            (r2_expect_two with the input-swapped operator, see R5);
//   layer 3: linearity of the trace over fixed open-leg labels and the
//            commutation of slicing with the pair construction.
// Operators are loaded with wrap_twosite_gate (input swap only): the direct
// trace path needs the same convention adapter as the simple-update gate.

namespace {

// The open variant of make_r2_tensor: every leg is {even, odd}. Same formula
// as fock_oracle.deterministic_tensor(site, parities, seed).
static ft make_mf2_tensor(int site, int seed) {
  const tenes::fermion::parity_vector p{false, true};
  const tenes::fermion::leg_parities par{p, p, p, p, p};
  tenes::real_tensor t(mptensor::Shape(2, 2, 2, 2, 2));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    auto idx = t.global_index(n);
    if (tenes::fermion::count_odd(par, idx) % 2 == 0) {
      const double x = static_cast<double>(
          (site + 2) * (1 + seed + (3 + seed % 5) * idx[0] +
                        (4 + seed % 5) * idx[1] + (5 + seed % 5) * idx[2] +
                        (6 + seed % 5) * idx[3] + (7 + seed % 5) * idx[4]));
      t.set_value(idx, 0.19 * std::sin(x) + 0.13 * std::cos(0.37 * x));
    }
  }
  return ft{t, par};
}

static const std::vector<double> mf_lambda{1.0, 0.7};

// Legs that stay open once the pair is formed (the internal bond A.r-B.l or
// A.b-B.t is left undressed). Leg indices: 0 l, 1 t, 2 r, 3 b.
static std::vector<int> mf_open_legs_a(
    tenes::fermion::reduced_pair_direction direction) {
  if (direction == tenes::fermion::reduced_pair_direction::horizontal) {
    return {0, 1, 3};
  }
  return {0, 1, 2};
}

static std::vector<int> mf_open_legs_b(
    tenes::fermion::reduced_pair_direction direction) {
  if (direction == tenes::fermion::reduced_pair_direction::horizontal) {
    return {1, 2, 3};
  }
  return {0, 2, 3};
}

static ft mf_dress(ft a, const std::vector<int>& legs) {
  for (const int leg : legs) {
    a.t.multiply_vector(mf_lambda, leg);
  }
  return a;
}

// Site 0 is A (left / top), site 1 is B (right / bottom), lambda on the open
// legs only.
static ft mf2_site_a(tenes::fermion::reduced_pair_direction direction,
                     int seed) {
  return mf_dress(make_mf2_tensor(0, seed), mf_open_legs_a(direction));
}

static ft mf2_site_b(tenes::fermion::reduced_pair_direction direction,
                     int seed) {
  return mf_dress(make_mf2_tensor(1, seed), mf_open_legs_b(direction));
}

static ft mf2_wrap(const ft& op) {
  const tenes::fermion::parity_vector p{false, true};
  return tenes::fermion::wrap_twosite_gate(op.t, p, p);
}

static ft mf2_nn_op() {
  ft op{tenes::real_tensor(mptensor::Shape(2, 2, 2, 2)),
        {{false, true}, {false, true}, {false, true}, {false, true}}};
  op.t.set_value(mptensor::Index(1, 1, 1, 1), 1.0);
  return op;
}

// Axes on the pair state: (l_A, t_A, b_A, s_A, t_B, r_B, b_B, s_B) for
// horizontal, (l_A, t_A, r_A, s_A, l_B, r_B, b_B, s_B) for vertical; the
// open legs of A and B land on pair axes {0, 1, 2} and {4, 5, 6}.
static const int mf_pair_open_axes[6] = {0, 1, 2, 4, 5, 6};

struct mf_oracle_ref {
  const char* name;
  tenes::fermion::reduced_pair_direction direction;
  int seed;
  double norm;
  double n0;
  double n1;
  double hop01;
  double pair01;
  double nn01;
};

// venv/bin/python3 test/fermion/fock_oracle.py | grep '^mf_'
static const mf_oracle_ref mf_refs[] = {
    {"mf_horizontal_2site", tenes::fermion::reduced_pair_direction::horizontal,
     0, 2.60686408978207385e-02, 2.50466769038458825e-01,
     6.27313313350885027e-01, -3.73251249983831065e-02, 7.91808258802335907e-02,
     1.53798349354499875e-01},
    {"mf_vertical_2site", tenes::fermion::reduced_pair_direction::vertical, 0,
     2.71047999509600387e-02, 3.15800751792648993e-01, 5.54446879147443394e-01,
     4.51903604013348132e-03, -2.89028494305034822e-03,
     1.76906556259372849e-01},
    {"mf_seed173_horizontal_2site",
     tenes::fermion::reduced_pair_direction::horizontal, 173,
     2.28155787945912847e-02, 4.55938311982168720e-01, 6.10843923486141693e-01,
     1.02835266183819024e-01, 1.01072667937597435e-01, 2.74730680754215972e-01},
    {"mf_seed173_vertical_2site",
     tenes::fermion::reduced_pair_direction::vertical, 173,
     1.90343286632445076e-02, 4.78271666664231543e-01, 4.89389198205553710e-01,
     9.45737037673597195e-03, -9.68639099935675665e-03,
     3.14906351782756222e-01}};

static ft mf4_site(int site, tenes::fermion::reduced_pair_direction direction,
                   int seed) {
  if (direction == tenes::fermion::reduced_pair_direction::horizontal) {
    return make_r4_tensor(2, 1, site, seed);
  }
  return make_r4_tensor(1, 2, site, seed);
}

// The direct f-primitive pair state (R5), independent of build_pair_state.
static ft mf_direct_pair(const ft& a, const ft& b,
                         tenes::fermion::reduced_pair_direction direction) {
  if (direction == tenes::fermion::reduced_pair_direction::horizontal) {
    return tenes::fermion::tensordot(a, b, mptensor::Axes(2),
                                     mptensor::Axes(0));
  }
  return tenes::fermion::tensordot(a, b, mptensor::Axes(3), mptensor::Axes(1));
}

static tenes::complex_tensor mf_complexify(const tenes::real_tensor& a) {
  tenes::complex_tensor t(a.shape());
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    auto idx = a.global_index(n);
    double v;
    a.get_value(idx, v);
    t.set_value(idx, std::complex<double>(v, 0.0));
  }
  return t;
}

static tenes::fermion::ftensor<tenes::complex_tensor> mf_to_complex(
    const ft& a) {
  return tenes::fermion::ftensor<tenes::complex_tensor>{mf_complexify(a.t),
                                                        a.parity};
}

}  // namespace

TEST_CASE("MF layer1 d=2 pair measurement matches the Fock oracle") {
  using tenes::fermion::build_pair_state;
  using tenes::fermion::contract_pair_MF;
  const ft n_a = mf2_wrap(r2_density_a_op());
  const ft n_b = mf2_wrap(r2_density_b_op());
  const ft hop = mf2_wrap(r2_hop_op());
  const ft nn = mf2_wrap(mf2_nn_op());

  for (const auto& ref : mf_refs) {
    const ft a = mf2_site_a(ref.direction, ref.seed);
    const ft b = mf2_site_b(ref.direction, ref.seed);
    const ft pair = build_pair_state(a, b, ref.direction);
    REQUIRE(pair.rank() == 8);

    const double norm = contract_pair_MF(pair);
    const double n0 = contract_pair_MF(pair, n_a) / norm;
    const double n1 = contract_pair_MF(pair, n_b) / norm;
    const double hop01 = contract_pair_MF(pair, hop) / norm;
    const double nn01 = contract_pair_MF(pair, nn) / norm;
    std::cout << std::setprecision(17) << "MF layer1 " << ref.name
              << " norm=" << norm << " n0=" << n0 << " n1=" << n1
              << " hop01=" << hop01 << " nn01=" << nn01 << std::endl;

    CHECK(norm == doctest::Approx(ref.norm).epsilon(1e-12));
    CHECK(n0 == doctest::Approx(ref.n0).epsilon(1e-12));
    CHECK(n1 == doctest::Approx(ref.n1).epsilon(1e-12));
    CHECK(hop01 == doctest::Approx(ref.hop01).epsilon(1e-12));
    CHECK(nn01 == doctest::Approx(ref.nn01).epsilon(1e-12));
  }
}

// The oracle's pair01 is 2 pairing(b, a) / norm, the expectation of the
// ordered-basis matrix |00><11| + |11><00| = c_b c_a + c_a^dag c_b^dag, so
// it is the wrapped r2_pair_op that must reproduce it; the unwrapped op is
// the antisymmetric partner c_b c_a - c_a^dag c_b^dag, which vanishes on a
// real state.
TEST_CASE("MF layer1 d=2 pair01 matches the Fock oracle pairing") {
  using tenes::fermion::build_pair_state;
  using tenes::fermion::contract_pair_MF;
  const ft pair_op = mf2_wrap(r2_pair_op());
  const ft pair_op_plain = r2_pair_op();

  for (const auto& ref : mf_refs) {
    const ft a = mf2_site_a(ref.direction, ref.seed);
    const ft b = mf2_site_b(ref.direction, ref.seed);
    const ft pair = build_pair_state(a, b, ref.direction);
    const double norm = contract_pair_MF(pair);
    const double pair01 = contract_pair_MF(pair, pair_op) / norm;
    const double pair01_plain = contract_pair_MF(pair, pair_op_plain) / norm;
    std::cout << std::setprecision(17) << "MF layer1 " << ref.name
              << " pair01=" << pair01 << " plain=" << pair01_plain << std::endl;
    CHECK(pair01 == doctest::Approx(ref.pair01).epsilon(1e-12));
    // negative control: verbatim loading is a different operator
    CHECK(std::abs(pair01_plain - ref.pair01) > 0.1 * std::abs(ref.pair01));
  }
}

TEST_CASE("MF layer2 d=4 pair measurement matches the direct path") {
  using tenes::fermion::build_pair_state;
  using tenes::fermion::contract_pair_MF;
  using tenes::fermion::reduced_pair_direction;
  for (const int seed : {0, 173}) {
    for (const auto direction : {reduced_pair_direction::horizontal,
                                 reduced_pair_direction::vertical}) {
      const bool horizontal = direction == reduced_pair_direction::horizontal;
      const std::string label =
          std::string(horizontal ? "H" : "V") + " seed=" + std::to_string(seed);
      const ft a = mf4_site(0, direction, seed);
      const ft b = mf4_site(1, direction, seed);
      const ft psi = mf_direct_pair(a, b, direction);
      const ft pair = build_pair_state(a, b, direction);

      // leg order (l_A, t_A, b_A, s_A, t_B, r_B, b_B, s_B) or
      // (l_A, t_A, r_A, s_A, l_B, r_B, b_B, s_B): the physical legs carry
      // r4_phys, the others the bond or boundary parity of their source.
      REQUIRE(pair.rank() == 8);
      CHECK(pair.parity[0] == a.parity[0]);
      CHECK(pair.parity[1] == a.parity[1]);
      CHECK(pair.parity[2] == a.parity[horizontal ? 3 : 2]);
      CHECK(pair.parity[3] == r4_phys);
      CHECK(pair.parity[4] == b.parity[horizontal ? 1 : 0]);
      CHECK(pair.parity[5] == b.parity[2]);
      CHECK(pair.parity[6] == b.parity[3]);
      CHECK(pair.parity[7] == r4_phys);

      const double norm = contract_pair_MF(pair);
      CHECK(norm == doctest::Approx(r2_norm(psi)).epsilon(1e-12));

      const ft hop = r4_wrap(r4_hop_plain(), true, false);
      const ft nn = r4_wrap(r4_nn_plain(), true, false);
      const double hop_ref = r2_expect_two(psi, 3, 7, hop);
      const double nn_ref = r2_expect_two(psi, 3, 7, nn);
      const double hop_mf = contract_pair_MF(pair, hop) / norm;
      const double nn_mf = contract_pair_MF(pair, nn) / norm;
      CHECK(hop_mf == doctest::Approx(hop_ref).epsilon(1e-12));
      CHECK(nn_mf == doctest::Approx(nn_ref).epsilon(1e-12));

      // negative controls: both swaps (blob convention) and no swap
      const double hop_both =
          contract_pair_MF(pair, r4_wrap(r4_hop_plain(), true, true)) / norm;
      const double hop_plain =
          contract_pair_MF(pair, r4_wrap(r4_hop_plain(), false, false)) / norm;
      const double hop_scale = std::max(1.0e-6, std::abs(hop_ref));
      CHECK(std::abs(hop_both - hop_ref) > 0.1 * hop_scale);
      CHECK(std::abs(hop_plain - hop_ref) > 0.1 * hop_scale);

      // measuring from the B end: the graded transpose of the operator is
      // the same bond operator seen from B, the plain transpose is not.
      const ft o_b = tenes::fermion::transpose(hop, mptensor::Axes(1, 0, 3, 2));
      const ft o_b_plain{mptensor::transpose(hop.t, mptensor::Axes(1, 0, 3, 2)),
                         hop.parity};
      const double hop_b = contract_pair_MF(pair, o_b) / norm;
      const double hop_b_plain = contract_pair_MF(pair, o_b_plain) / norm;
      CHECK(hop_b == doctest::Approx(hop_ref).epsilon(1e-12));
      CHECK(std::abs(hop_b_plain - hop_ref) > 0.1 * hop_scale);

      std::cout << std::setprecision(17) << "MF layer2 " << label
                << " norm=" << norm << " hop_ref=" << hop_ref
                << " hop_mf=" << hop_mf << " hop_both=" << hop_both
                << " hop_plain=" << hop_plain << " hop_b=" << hop_b
                << " hop_b_plain=" << hop_b_plain << " nn_ref=" << nn_ref
                << " nn_mf=" << nn_mf << std::endl;
    }
  }
}

TEST_CASE("MF layer3 contract_pair_MF sums over fixed open-leg labels") {
  using tenes::fermion::build_pair_state;
  using tenes::fermion::contract_pair_MF;
  constexpr auto direction = tenes::fermion::reduced_pair_direction::horizontal;
  const ft a = mf2_site_a(direction, 0);
  const ft b = mf2_site_b(direction, 0);
  const ft hop = mf2_wrap(r2_hop_op());
  const std::vector<int> legs_a = mf_open_legs_a(direction);
  const std::vector<int> legs_b = mf_open_legs_b(direction);

  const ft pair = build_pair_state(a, b, direction);
  const double norm = contract_pair_MF(pair);
  const double hop_value = contract_pair_MF(pair, hop);

  double norm_sum = 0.0;
  double hop_sum = 0.0;
  for (int labels = 0; labels < 64; ++labels) {
    ft a_x = a;
    ft b_x = b;
    for (int i = 0; i < 3; ++i) {
      const std::size_t la = (labels >> i) & 1;
      const std::size_t lb = (labels >> (3 + i)) & 1;
      a_x = tenes::fermion::slice(a_x, legs_a[i], la, la + 1);
      b_x = tenes::fermion::slice(b_x, legs_b[i], lb, lb + 1);
    }
    const ft pair_x = build_pair_state(a_x, b_x, direction);
    norm_sum += contract_pair_MF(pair_x);
    hop_sum += contract_pair_MF(pair_x, hop);
  }
  std::cout << std::setprecision(17) << "MF layer3 norm=" << norm
            << " sum=" << norm_sum << " hop=" << hop_value << " sum=" << hop_sum
            << std::endl;
  CHECK(norm_sum == doctest::Approx(norm).epsilon(1e-12));
  CHECK(hop_sum == doctest::Approx(hop_value).epsilon(1e-12));
}

TEST_CASE("MF layer3 slicing open legs commutes with build_pair_state") {
  using tenes::fermion::build_pair_state;
  constexpr auto direction = tenes::fermion::reduced_pair_direction::horizontal;
  const ft a = mf2_site_a(direction, 0);
  const ft b = mf2_site_b(direction, 0);
  const std::vector<int> legs_a = mf_open_legs_a(direction);
  const std::vector<int> legs_b = mf_open_legs_b(direction);

  // labels in the order (l_A, t_A, b_A, t_B, r_B, b_B): every open leg on
  // the odd label, then a mixed pattern that tells the open axes apart
  const std::array<std::array<std::size_t, 6>, 2> patterns = {
      {{1, 1, 1, 1, 1, 1}, {1, 0, 1, 0, 1, 0}}};
  for (const auto& labels : patterns) {
    ft a_x = a;
    ft b_x = b;
    for (int i = 0; i < 3; ++i) {
      a_x = tenes::fermion::slice(a_x, legs_a[i], labels[i], labels[i] + 1);
      b_x = tenes::fermion::slice(b_x, legs_b[i], labels[3 + i],
                                  labels[3 + i] + 1);
    }
    const ft from_sites = build_pair_state(a_x, b_x, direction);

    ft from_pair = build_pair_state(a, b, direction);
    for (int i = 0; i < 6; ++i) {
      from_pair = tenes::fermion::slice(from_pair, mf_pair_open_axes[i],
                                        labels[i], labels[i] + 1);
    }

    REQUIRE(from_sites.rank() == 8);
    REQUIRE(from_sites.shape() == from_pair.shape());
    CHECK(from_sites.parity == from_pair.parity);
    for (int i = 0; i < 6; ++i) {
      // a one-dimensional leg fixed to an odd label stays odd
      CHECK(from_sites.parity[mf_pair_open_axes[i]] ==
            tenes::fermion::parity_vector{labels[i] == 1});
    }
    CHECK(from_sites.parity[3] == tenes::fermion::parity_vector{false, true});
    CHECK(from_sites.parity[7] == tenes::fermion::parity_vector{false, true});
    for (std::size_t n = 0; n < from_sites.t.local_size(); ++n) {
      auto idx = from_sites.t.global_index(n);
      double got;
      double expected;
      from_sites.t.get_value(idx, got);
      from_pair.t.get_value(idx, expected);
      CHECK(got == doctest::Approx(expected).epsilon(1e-12));
    }
  }
}

TEST_CASE("MF complex pair state reproduces the real measurement") {
  using cft = tenes::fermion::ftensor<tenes::complex_tensor>;
  using tenes::fermion::build_pair_state;
  using tenes::fermion::contract_pair_MF;
  constexpr auto direction = tenes::fermion::reduced_pair_direction::horizontal;
  const ft a = mf2_site_a(direction, 0);
  const ft b = mf2_site_b(direction, 0);
  const ft hop = mf2_wrap(r2_hop_op());

  const ft pair = build_pair_state(a, b, direction);
  const double norm = contract_pair_MF(pair);
  const double hop_value = contract_pair_MF(pair, hop);

  const cft a_c = mf_to_complex(a);
  const cft b_c = mf_to_complex(b);
  const tenes::fermion::parity_vector p{false, true};
  const cft hop_c =
      tenes::fermion::wrap_twosite_gate(mf_complexify(r2_hop_op().t), p, p);
  const cft pair_c = build_pair_state(a_c, b_c, direction);
  const std::complex<double> norm_c = contract_pair_MF(pair_c);
  const std::complex<double> hop_c_value = contract_pair_MF(pair_c, hop_c);
  std::cout << std::setprecision(17) << "MF complex norm=" << norm_c
            << " hop=" << hop_c_value << " real norm=" << norm
            << " hop=" << hop_value << std::endl;

  CHECK(norm_c.real() == doctest::Approx(norm).epsilon(1e-12));
  CHECK(std::abs(norm_c.imag()) <= 1e-12);
  CHECK(hop_c_value.real() == doctest::Approx(hop_value).epsilon(1e-12));
  CHECK(std::abs(hop_c_value.imag()) <= 1e-12);
}

TEST_CASE("MF build_pair_state and apply_pair_op reject wrong ranks") {
  using tenes::fermion::apply_pair_op;
  using tenes::fermion::build_pair_state;
  using tenes::fermion::reduced_pair_direction;
  const ft a = mf2_site_a(reduced_pair_direction::horizontal, 0);
  const ft b = mf2_site_b(reduced_pair_direction::horizontal, 0);
  const ft hop = mf2_wrap(r2_hop_op());

  // a four-leg tensor with the physical leg dropped: tensordot over the
  // bond alone would still succeed, so only the rank guard can refuse it
  const ft a4 = tenes::fermion::tensordot(
      a, ft{tenes::real_tensor(mptensor::Shape(2)), {{false, true}}},
      mptensor::Axes(4), mptensor::Axes(0));
  REQUIRE(a4.rank() == 4);
  CHECK_THROWS_AS(build_pair_state(a4, b, reduced_pair_direction::horizontal),
                  std::runtime_error);
  CHECK_THROWS_AS(build_pair_state(a, a4, reduced_pair_direction::vertical),
                  std::runtime_error);

  // apply_pair_op wants an eight-leg pair and a four-leg operator
  CHECK_THROWS_AS(apply_pair_op(a, hop), std::runtime_error);
  CHECK_THROWS_AS(
      apply_pair_op(mf_direct_pair(a, b, reduced_pair_direction::horizontal),
                    r2_number_op()),
      std::runtime_error);
}
