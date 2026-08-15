namespace {

static tenes::fermion::leg_parities r2_parities(int lx, int ly, int site) {
  const int x = site % lx;
  const int y = site / lx;
  tenes::fermion::parity_vector even{false};
  tenes::fermion::parity_vector p{false, true};
  return {x > 0 ? p : even, y > 0 ? p : even, x + 1 < lx ? p : even,
          y + 1 < ly ? p : even, p};
}

static ft make_r2_tensor(int lx, int ly, int site, int seed = 0) {
  const auto p = r2_parities(lx, ly, site);
  tenes::real_tensor t(mptensor::Shape(p[0].size(), p[1].size(), p[2].size(),
                                       p[3].size(), p[4].size()));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    auto idx = t.global_index(n);
    if (tenes::fermion::count_odd(p, idx) % 2 == 0) {
      const double x = static_cast<double>(
          (site + 2) * (1 + seed + (3 + seed % 5) * idx[0] +
                        (4 + seed % 5) * idx[1] + (5 + seed % 5) * idx[2] +
                        (6 + seed % 5) * idx[3] + (7 + seed % 5) * idx[4]));
      t.set_value(idx, 0.19 * std::sin(x) + 0.13 * std::cos(0.37 * x));
    }
  }
  return ft{t, p};
}

static ft r2_apply_one_site(const ft& psi, int axis, const ft& op) {
  ft applied = tenes::fermion::tensordot(psi, op, mptensor::Axes(axis),
                                         mptensor::Axes(0));
  const int rank = psi.rank();
  mptensor::Axes perm;
  for (int i = 0; i < axis; ++i) {
    perm.push(i);
  }
  perm.push(rank - 1);
  for (int i = axis; i < rank - 1; ++i) {
    perm.push(i);
  }
  return tenes::fermion::transpose(applied, perm);
}

static ft r2_apply_two_site(const ft& psi, int axis0, int axis1, const ft& op) {
  ft applied = tenes::fermion::tensordot(psi, op, mptensor::Axes(axis0, axis1),
                                         mptensor::Axes(0, 1));
  const int rank = psi.rank();
  const int out0 = rank - 2;
  const int out1 = rank - 1;
  int free_axis = 0;
  mptensor::Axes perm;
  for (int ax = 0; ax < rank; ++ax) {
    if (ax == axis0) {
      perm.push(out0);
    } else if (ax == axis1) {
      perm.push(out1);
    } else {
      perm.push(free_axis++);
    }
  }
  return tenes::fermion::transpose(applied, perm);
}

static double r2_norm(const ft& psi) {
  mptensor::Axes axes;
  for (int ax = 0; ax < psi.rank(); ++ax) {
    axes.push(ax);
  }
  return tenes::fermion::trace(tenes::fermion::conj(psi), psi, axes, axes);
}

static double r2_expect_one(const ft& psi, int axis, const ft& op) {
  mptensor::Axes axes;
  for (int ax = 0; ax < psi.rank(); ++ax) {
    axes.push(ax);
  }
  return tenes::fermion::trace(tenes::fermion::conj(psi),
                               r2_apply_one_site(psi, axis, op), axes, axes) /
         r2_norm(psi);
}

static double r2_expect_two(const ft& psi, int axis0, int axis1, const ft& op) {
  mptensor::Axes axes;
  for (int ax = 0; ax < psi.rank(); ++ax) {
    axes.push(ax);
  }
  return tenes::fermion::trace(tenes::fermion::conj(psi),
                               r2_apply_two_site(psi, axis0, axis1, op), axes,
                               axes) /
         r2_norm(psi);
}

static ft r2_number_op() {
  ft op{tenes::real_tensor(mptensor::Shape(2, 2)),
        {{false, true}, {false, true}}};
  op.t.set_value(mptensor::Index(1, 1), 1.0);
  return op;
}

static ft r2_hop_op() {
  ft op{tenes::real_tensor(mptensor::Shape(2, 2, 2, 2)),
        {{false, true}, {false, true}, {false, true}, {false, true}}};
  op.t.set_value(mptensor::Index(0, 1, 1, 0), 1.0);
  op.t.set_value(mptensor::Index(1, 0, 0, 1), 1.0);
  return op;
}

static ft r2_pair_op() {
  ft op{tenes::real_tensor(mptensor::Shape(2, 2, 2, 2)),
        {{false, true}, {false, true}, {false, true}, {false, true}}};
  op.t.set_value(mptensor::Index(1, 1, 0, 0), 1.0);
  op.t.set_value(mptensor::Index(0, 0, 1, 1), 1.0);
  return op;
}

static ft r2_horizontal_pair() {
  return tenes::fermion::tensordot(make_r2_tensor(2, 1, 0),
                                   make_r2_tensor(2, 1, 1), mptensor::Axes(2),
                                   mptensor::Axes(0));
}

static ft r2_vertical_pair() {
  return tenes::fermion::tensordot(make_r2_tensor(1, 2, 0),
                                   make_r2_tensor(1, 2, 1), mptensor::Axes(3),
                                   mptensor::Axes(1));
}

static ft r2_plaquette() {
  ft top = tenes::fermion::tensordot(make_r2_tensor(2, 2, 0),
                                     make_r2_tensor(2, 2, 1), mptensor::Axes(2),
                                     mptensor::Axes(0));
  ft bottom = tenes::fermion::tensordot(make_r2_tensor(2, 2, 2),
                                        make_r2_tensor(2, 2, 3),
                                        mptensor::Axes(2), mptensor::Axes(0));
  return tenes::fermion::tensordot(top, bottom, mptensor::Axes(2, 6),
                                   mptensor::Axes(1, 4));
}

static double r3_sum_entries(const tenes::real_tensor& tensor) {
  double ret = 0.0;
  for (std::size_t n = 0; n < tensor.local_size(); ++n) {
    double v;
    tensor.get_value(tensor.global_index(n), v);
    ret += v;
  }
  return ret;
}

static tenes::real_tensor r3_number_op_tensor() {
  tenes::real_tensor op(mptensor::Shape(2, 2));
  op.set_value(mptensor::Index(1, 1), 1.0);
  return op;
}

struct r3_dataset_ref {
  int seed;
  double h_norm;
  double h_n0;
  double h_n1;
  double v_norm;
  double v_n0;
  double v_n1;
  double p_norm;
  std::array<double, 4> p_n;
  std::array<double, 4> p_hop;
};

static const r3_dataset_ref r3_refs[] = {
    {0,
     5.46520500622191588e-04,
     5.35874068543760670e-02,
     5.35874068543760670e-02,
     5.23025921617477092e-04,
     1.10740924627062767e-02,
     1.10740924627062767e-02,
     9.97335335163249894e-06,
     {2.08718041601653991e-02, 9.04495304092564401e-01, 4.35775657274577038e-02,
      9.25148998317175564e-01},
     {2.73390825968589414e-02, 1.94309350032059692e-03, 9.63079101517500154e-03,
      -6.61038902919681842e-02}},
    {173,
     6.99752779917210531e-05,
     9.97583488910690597e-01,
     9.97583488910690597e-01,
     1.10983813879685873e-04,
     9.98476390120925150e-01,
     9.98476390120925150e-01,
     1.76986289332680877e-06,
     {1.61225025287331764e-01, 9.94507681265915267e-01, 1.63971607379076845e-01,
      9.87796732360337271e-01},
     {-8.33515521666127329e-02, 6.24562346176388110e-03,
      1.23072758949854030e-03, 1.50398236202810870e-01}}};

static tenes::real_tensor r3_reduced_site(
    int lx, int ly, int site, bool with_number_op,
    const tenes::fermion::reduced_variant& variant, int seed = 0) {
  if (!with_number_op) {
    return tenes::fermion::build_reduced(make_r2_tensor(lx, ly, site, seed),
                                         variant);
  }
  return mptensor::tensordot(tenes::fermion::build_reduced_op(
                                 make_r2_tensor(lx, ly, site, seed), variant),
                             r3_number_op_tensor(), mptensor::Axes(4, 5),
                             mptensor::Axes(0, 1));
}

static double r3_horizontal_norm(bool op0, bool op1,
                                 const tenes::fermion::reduced_variant& variant,
                                 int seed) {
  return r3_sum_entries(
      mptensor::tensordot(r3_reduced_site(2, 1, 0, op0, variant, seed),
                          r3_reduced_site(2, 1, 1, op1, variant, seed),
                          mptensor::Axes(2), mptensor::Axes(0)));
}

static double r3_vertical_norm(bool op0, bool op1,
                               const tenes::fermion::reduced_variant& variant,
                               int seed) {
  return r3_sum_entries(
      mptensor::tensordot(r3_reduced_site(1, 2, 0, op0, variant, seed),
                          r3_reduced_site(1, 2, 1, op1, variant, seed),
                          mptensor::Axes(3), mptensor::Axes(1)));
}

static double r3_plaquette_norm(int number_site,
                                const tenes::fermion::reduced_variant& variant,
                                int seed) {
  tenes::real_tensor top = mptensor::tensordot(
      r3_reduced_site(2, 2, 0, number_site == 0, variant, seed),
      r3_reduced_site(2, 2, 1, number_site == 1, variant, seed),
      mptensor::Axes(2), mptensor::Axes(0));
  tenes::real_tensor bottom = mptensor::tensordot(
      r3_reduced_site(2, 2, 2, number_site == 2, variant, seed),
      r3_reduced_site(2, 2, 3, number_site == 3, variant, seed),
      mptensor::Axes(2), mptensor::Axes(0));
  return r3_sum_entries(mptensor::tensordot(top, bottom, mptensor::Axes(2, 5),
                                            mptensor::Axes(1, 3)));
}

static double r3_plaquette_norm(int number_site = -1) {
  return r3_plaquette_norm(number_site, tenes::fermion::reduced_variant{}, 0);
}

static double r3_pair_plaquette_norm(
    int a, int b, const tenes::fermion::reduced_variant& variant, int seed) {
  if (a == 0 && b == 1) {
    tenes::real_tensor blob = tenes::fermion::build_reduced_pair(
        make_r2_tensor(2, 2, 0, seed), make_r2_tensor(2, 2, 1, seed),
        r2_hop_op(), tenes::fermion::reduced_pair_direction::horizontal,
        variant);
    tenes::real_tensor bottom =
        mptensor::tensordot(r3_reduced_site(2, 2, 2, false, variant, seed),
                            r3_reduced_site(2, 2, 3, false, variant, seed),
                            mptensor::Axes(2), mptensor::Axes(0));
    return r3_sum_entries(mptensor::tensordot(
        blob, bottom, mptensor::Axes(2, 5), mptensor::Axes(1, 3)));
  }
  if (a == 2 && b == 3) {
    tenes::real_tensor top =
        mptensor::tensordot(r3_reduced_site(2, 2, 0, false, variant, seed),
                            r3_reduced_site(2, 2, 1, false, variant, seed),
                            mptensor::Axes(2), mptensor::Axes(0));
    tenes::real_tensor blob = tenes::fermion::build_reduced_pair(
        make_r2_tensor(2, 2, 2, seed), make_r2_tensor(2, 2, 3, seed),
        r2_hop_op(), tenes::fermion::reduced_pair_direction::horizontal,
        variant);
    return r3_sum_entries(mptensor::tensordot(top, blob, mptensor::Axes(2, 5),
                                              mptensor::Axes(1, 3)));
  }
  if (a == 0 && b == 2) {
    tenes::real_tensor blob = tenes::fermion::build_reduced_pair(
        make_r2_tensor(2, 2, 0, seed), make_r2_tensor(2, 2, 2, seed),
        r2_hop_op(), tenes::fermion::reduced_pair_direction::vertical, variant);
    tenes::real_tensor right =
        mptensor::tensordot(r3_reduced_site(2, 2, 1, false, variant, seed),
                            r3_reduced_site(2, 2, 3, false, variant, seed),
                            mptensor::Axes(3), mptensor::Axes(1));
    return r3_sum_entries(mptensor::tensordot(blob, right, mptensor::Axes(2, 4),
                                              mptensor::Axes(0, 3)));
  }
  tenes::real_tensor left =
      mptensor::tensordot(r3_reduced_site(2, 2, 0, false, variant, seed),
                          r3_reduced_site(2, 2, 2, false, variant, seed),
                          mptensor::Axes(3), mptensor::Axes(1));
  tenes::real_tensor blob = tenes::fermion::build_reduced_pair(
      make_r2_tensor(2, 2, 1, seed), make_r2_tensor(2, 2, 3, seed), r2_hop_op(),
      tenes::fermion::reduced_pair_direction::vertical, variant);
  return r3_sum_entries(mptensor::tensordot(left, blob, mptensor::Axes(2, 4),
                                            mptensor::Axes(0, 3)));
}

static double r3_pair_plaquette_norm(int a, int b) {
  return r3_pair_plaquette_norm(a, b, tenes::fermion::reduced_variant{}, 0);
}

struct r3_eval {
  std::array<double, 16> values;
};

static r3_eval r3_evaluate(const r3_dataset_ref& ref,
                           const tenes::fermion::reduced_variant& variant) {
  r3_eval ret{};
  const double h_norm = r3_horizontal_norm(false, false, variant, ref.seed);
  const double v_norm = r3_vertical_norm(false, false, variant, ref.seed);
  const double p_norm = r3_plaquette_norm(-1, variant, ref.seed);
  ret.values = {h_norm,
                r3_horizontal_norm(true, false, variant, ref.seed) / h_norm,
                r3_horizontal_norm(false, true, variant, ref.seed) / h_norm,
                v_norm,
                r3_vertical_norm(true, false, variant, ref.seed) / v_norm,
                r3_vertical_norm(false, true, variant, ref.seed) / v_norm,
                p_norm,
                r3_plaquette_norm(0, variant, ref.seed) / p_norm,
                r3_plaquette_norm(1, variant, ref.seed) / p_norm,
                r3_plaquette_norm(2, variant, ref.seed) / p_norm,
                r3_plaquette_norm(3, variant, ref.seed) / p_norm,
                r3_pair_plaquette_norm(0, 1, variant, ref.seed) / p_norm,
                r3_pair_plaquette_norm(0, 2, variant, ref.seed) / p_norm,
                r3_pair_plaquette_norm(1, 3, variant, ref.seed) / p_norm,
                r3_pair_plaquette_norm(2, 3, variant, ref.seed) / p_norm,
                0.0};
  return ret;
}

static std::array<double, 16> r3_reference_values(const r3_dataset_ref& ref) {
  return {ref.h_norm,   ref.h_n0,     ref.h_n1,     ref.v_norm,
          ref.v_n0,     ref.v_n1,     ref.p_norm,   ref.p_n[0],
          ref.p_n[1],   ref.p_n[2],   ref.p_n[3],   ref.p_hop[0],
          ref.p_hop[1], ref.p_hop[2], ref.p_hop[3], 0.0};
}

static bool r3_close(double got, double ref) {
  return std::abs(got - ref) <= 1.0e-10 * std::max(1.0, std::abs(ref));
}

static double r3_variant_score(const tenes::fermion::reduced_variant& variant,
                               std::array<double, 16>* first_values = nullptr) {
  double score = 0.0;
  for (std::size_t d = 0; d < 2; ++d) {
    const auto eval = r3_evaluate(r3_refs[d], variant);
    const auto ref = r3_reference_values(r3_refs[d]);
    if (d == 0 && first_values != nullptr) {
      *first_values = eval.values;
    }
    for (std::size_t i = 0; i < ref.size() - 1; ++i) {
      score += std::abs(eval.values[i] - ref[i]) /
               std::max(1.0e-30, std::abs(ref[i]));
    }
  }
  return score;
}

static bool r3_variant_passes(const tenes::fermion::reduced_variant& variant) {
  for (const auto& dataset : r3_refs) {
    const auto eval = r3_evaluate(dataset, variant);
    const auto ref = r3_reference_values(dataset);
    for (std::size_t i = 0; i < ref.size() - 1; ++i) {
      if (!r3_close(eval.values[i], ref[i])) {
        return false;
      }
    }
  }
  return true;
}

static std::string r3_mask_name(unsigned mask) {
  const char* names[] = {"l", "t", "r", "b"};
  std::ostringstream ss;
  ss << "{";
  bool first = true;
  for (int i = 0; i < 4; ++i) {
    if ((mask & (1u << i)) != 0) {
      if (!first) {
        ss << ",";
      }
      ss << names[i];
      first = false;
    }
  }
  ss << "}";
  return ss.str();
}

static std::string r3_variant_name(
    const tenes::fermion::reduced_variant& variant) {
  std::ostringstream ss;
  ss << "S_bra=" << r3_mask_name(variant.bra_mask)
     << " S_ket=" << r3_mask_name(variant.ket_mask)
     << " order=" << (variant.ket_first ? "(x,xbar)" : "(xbar,x)")
     << " phys_twist=" << (variant.phys_twist ? "true" : "false");
  return ss.str();
}

}  // namespace

TEST_CASE("R2 oracle agrees with f-primitives on two-site open patches") {
  const ft n_op = r2_number_op();
  const ft hop_op = r2_hop_op();
  const ft pair_op = r2_pair_op();

  const ft hp = r2_horizontal_pair();
  CHECK(r2_norm(hp) == doctest::Approx(5.46520500622191588e-04).epsilon(1e-12));
  CHECK(r2_expect_one(hp, 3, n_op) ==
        doctest::Approx(5.35874068543760670e-02).epsilon(1e-12));
  CHECK(r2_expect_one(hp, 7, n_op) ==
        doctest::Approx(5.35874068543760670e-02).epsilon(1e-12));
  CHECK(r2_expect_two(hp, 3, 7, hop_op) ==
        doctest::Approx(0.00000000000000000e+00).epsilon(1e-12));
  CHECK(r2_expect_two(hp, 3, 7, pair_op) ==
        doctest::Approx(0.00000000000000000e+00).epsilon(1e-12));

  const ft vp = r2_vertical_pair();
  CHECK(r2_norm(vp) == doctest::Approx(5.23025921617477092e-04).epsilon(1e-12));
  CHECK(r2_expect_one(vp, 3, n_op) ==
        doctest::Approx(1.10740924627062767e-02).epsilon(1e-12));
  CHECK(r2_expect_one(vp, 7, n_op) ==
        doctest::Approx(1.10740924627062767e-02).epsilon(1e-12));
  CHECK(r2_expect_two(vp, 3, 7, hop_op) ==
        doctest::Approx(0.00000000000000000e+00).epsilon(1e-12));
  CHECK(r2_expect_two(vp, 3, 7, pair_op) ==
        doctest::Approx(0.00000000000000000e+00).epsilon(1e-12));
}

TEST_CASE(
    "R2 plaquette loop comparison records current f-primitive convention") {
  const ft psi = r2_plaquette();
  const ft n_op = r2_number_op();
  const ft hop_op = r2_hop_op();
  const ft pair_op = r2_pair_op();

  const double norm = r2_norm(psi);
  const double n0 = r2_expect_one(psi, 2, n_op);
  const double n1 = r2_expect_one(psi, 5, n_op);
  const double n2 = r2_expect_one(psi, 8, n_op);
  const double n3 = r2_expect_one(psi, 11, n_op);
  const double hop01 = r2_expect_two(psi, 2, 5, hop_op);
  const double hop02 = r2_expect_two(psi, 2, 8, hop_op);
  const double hop13 = r2_expect_two(psi, 5, 11, hop_op);
  const double hop23 = r2_expect_two(psi, 8, 11, hop_op);
  const double pair01 = r2_expect_two(psi, 2, 5, pair_op);
  const double pair02 = r2_expect_two(psi, 2, 8, pair_op);

  std::cout << "R2 plaquette fprimitive norm=" << norm << " n=[" << n0 << ", "
            << n1 << ", " << n2 << ", " << n3 << "] hop=[" << hop01 << ", "
            << hop02 << ", " << hop13 << ", " << hop23 << "] pair=[" << pair01
            << ", " << pair02 << "]" << std::endl;

  CHECK(norm == doctest::Approx(9.97335335163249894e-06).epsilon(1e-12));
  CHECK(n0 == doctest::Approx(2.08718041601653991e-02).epsilon(1e-12));
  CHECK(n1 == doctest::Approx(9.04495304092564401e-01).epsilon(1e-12));
  CHECK(n2 == doctest::Approx(4.35775657274577038e-02).epsilon(1e-12));
  CHECK(n3 == doctest::Approx(9.25148998317175564e-01).epsilon(1e-12));
  CHECK(hop01 == doctest::Approx(2.73390825968589414e-02).epsilon(1e-12));
  CHECK(hop02 == doctest::Approx(1.94309350032059692e-03).epsilon(1e-12));
  CHECK(hop13 == doctest::Approx(9.63079101517500154e-03).epsilon(1e-12));
  CHECK(hop23 == doctest::Approx(-6.61038902919681842e-02).epsilon(1e-12));
  CHECK(pair01 == doctest::Approx(0.00000000000000000e+00).epsilon(1e-12));
  CHECK(pair02 == doctest::Approx(0.00000000000000000e+00).epsilon(1e-12));
}

TEST_CASE("R3 brute-force reduced convention search") {
  std::vector<tenes::fermion::reduced_variant> passing;
  struct candidate {
    double score;
    tenes::fermion::reduced_variant variant;
    std::array<double, 16> first_values;
  };
  std::vector<candidate> best;

  for (unsigned bra_mask = 0; bra_mask < 16; ++bra_mask) {
    for (unsigned ket_mask = 0; ket_mask < 16; ++ket_mask) {
      for (bool ket_first : {false, true}) {
        for (bool phys_twist : {false, true}) {
          tenes::fermion::reduced_variant variant{bra_mask, ket_mask, ket_first,
                                                  phys_twist};
          std::array<double, 16> first_values{};
          const double score = r3_variant_score(variant, &first_values);
          best.push_back(candidate{score, variant, first_values});
          if (r3_variant_passes(variant)) {
            passing.push_back(variant);
          }
        }
      }
    }
  }

  std::sort(
      best.begin(), best.end(),
      [](const candidate& a, const candidate& b) { return a.score < b.score; });

  std::cout << std::setprecision(17)
            << "R3 search passing count=" << passing.size() << std::endl;
  for (const auto& variant : passing) {
    std::cout << "R3 search passing " << r3_variant_name(variant) << std::endl;
  }

  if (passing.empty()) {
    const auto ref = r3_reference_values(r3_refs[0]);
    for (std::size_t i = 0; i < std::min<std::size_t>(10, best.size()); ++i) {
      std::cout << "R3 search best[" << i << "] score=" << best[i].score << " "
                << r3_variant_name(best[i].variant) << std::endl;
      static const char* labels[] = {"h_norm", "h_n0",  "h_n1",   "v_norm",
                                     "v_n0",   "v_n1",  "p_norm", "p_n0",
                                     "p_n1",   "p_n2",  "p_n3",   "hop01",
                                     "hop02",  "hop13", "hop23"};
      for (std::size_t j = 0; j < ref.size() - 1; ++j) {
        std::cout << "  " << labels[j] << " got=" << best[i].first_values[j]
                  << " ref=" << ref[j]
                  << " ratio=" << best[i].first_values[j] / ref[j]
                  << " diff=" << best[i].first_values[j] - ref[j] << std::endl;
      }
    }
  }

  CHECK(best.size() == 1024);
}

TEST_CASE("R3 reduced physical-open tensor traces to reduced norm tensor") {
  tenes::real_tensor traced = mptensor::contract(
      tenes::fermion::build_reduced_op(make_r2_tensor(2, 2, 0)),
      mptensor::Axes(4), mptensor::Axes(5));
  tenes::real_tensor reduced =
      tenes::fermion::build_reduced(make_r2_tensor(2, 2, 0));
  REQUIRE(traced.shape() == reduced.shape());
  for (std::size_t n = 0; n < reduced.local_size(); ++n) {
    auto idx = reduced.global_index(n);
    double got;
    double expected;
    traced.get_value(idx, got);
    reduced.get_value(idx, expected);
    CHECK(got == doctest::Approx(expected).epsilon(1e-12));
  }
}

TEST_CASE("R3 default reduced pair blob diagnostic is finite") {
  const double norm = r3_plaquette_norm();
  const double hop01 = r3_pair_plaquette_norm(0, 1) / norm;
  const double hop02 = r3_pair_plaquette_norm(0, 2) / norm;
  const double hop13 = r3_pair_plaquette_norm(1, 3) / norm;
  const double hop23 = r3_pair_plaquette_norm(2, 3) / norm;
  std::cout << std::setprecision(17) << "R3 reduced pair plaquette hop=["
            << hop01 << ", " << hop02 << ", " << hop13 << ", " << hop23 << "]"
            << std::endl;
  CHECK(std::isfinite(hop01));
  CHECK(std::isfinite(hop02));
  CHECK(std::isfinite(hop13));
  CHECK(std::isfinite(hop23));
}
