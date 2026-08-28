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

static ft r2_density_a_op() {
  ft op{tenes::real_tensor(mptensor::Shape(2, 2, 2, 2)),
        {{false, true}, {false, true}, {false, true}, {false, true}}};
  op.t.set_value(mptensor::Index(1, 0, 1, 0), 1.0);
  op.t.set_value(mptensor::Index(1, 1, 1, 1), 1.0);
  return op;
}

static ft r2_density_b_op() {
  ft op{tenes::real_tensor(mptensor::Shape(2, 2, 2, 2)),
        {{false, true}, {false, true}, {false, true}, {false, true}}};
  op.t.set_value(mptensor::Index(0, 1, 0, 1), 1.0);
  op.t.set_value(mptensor::Index(1, 1, 1, 1), 1.0);
  return op;
}

static ft r2_mixed_op() {
  constexpr double mu = 0.7;
  ft op{tenes::real_tensor(mptensor::Shape(2, 2, 2, 2)),
        {{false, true}, {false, true}, {false, true}, {false, true}}};
  op.t.set_value(mptensor::Index(0, 1, 1, 0), -1.0);
  op.t.set_value(mptensor::Index(1, 0, 0, 1), -1.0);
  op.t.set_value(mptensor::Index(0, 1, 0, 1), -0.25 * mu);
  op.t.set_value(mptensor::Index(1, 0, 1, 0), -0.25 * mu);
  op.t.set_value(mptensor::Index(1, 1, 1, 1), -0.5 * mu);
  return op;
}

static ft r2_density_pair_op() {
  ft op{tenes::real_tensor(mptensor::Shape(2, 2, 2, 2)),
        {{false, true}, {false, true}, {false, true}, {false, true}}};
  op.t.set_value(mptensor::Index(0, 1, 0, 1), 1.0);
  op.t.set_value(mptensor::Index(1, 0, 1, 0), 1.0);
  op.t.set_value(mptensor::Index(1, 1, 1, 1), 2.0);
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
  double h_mixed;
  double v_mixed;
  std::array<double, 4> p_mixed;
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
      -6.61038902919681842e-02},
     -1.87555923990316217e-02,
     -3.87593236194719651e-03,
     {-1.89278326541086644e-01, -1.32217332306546406e-02,
      -3.29818543936879482e-01, -1.03423258415842642e-01}},
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
      1.23072758949854030e-03, 1.50398236202810870e-01},
     -3.49154221118741670e-01,
     -3.49466736542323797e-01,
     {-1.18901671480205487e-01, -6.31550341783853886e-02,
      -3.48133999974092723e-01, -3.51957695657208358e-01}}};

static tenes::real_tensor r3_reduced_site(int lx, int ly, int site,
                                          bool with_number_op, int seed = 0) {
  if (!with_number_op) {
    return tenes::fermion::build_reduced(make_r2_tensor(lx, ly, site, seed));
  }
  return mptensor::tensordot(
      tenes::fermion::build_reduced_op(make_r2_tensor(lx, ly, site, seed)),
      r3_number_op_tensor(), mptensor::Axes(4, 5), mptensor::Axes(0, 1));
}

static double r3_horizontal_norm(bool op0, bool op1, int seed) {
  return r3_sum_entries(mptensor::tensordot(
      r3_reduced_site(2, 1, 0, op0, seed), r3_reduced_site(2, 1, 1, op1, seed),
      mptensor::Axes(2), mptensor::Axes(0)));
}

static double r3_vertical_norm(bool op0, bool op1, int seed) {
  return r3_sum_entries(mptensor::tensordot(
      r3_reduced_site(1, 2, 0, op0, seed), r3_reduced_site(1, 2, 1, op1, seed),
      mptensor::Axes(3), mptensor::Axes(1)));
}

static double r3_plaquette_norm(int number_site, int seed) {
  tenes::real_tensor top =
      mptensor::tensordot(r3_reduced_site(2, 2, 0, number_site == 0, seed),
                          r3_reduced_site(2, 2, 1, number_site == 1, seed),
                          mptensor::Axes(2), mptensor::Axes(0));
  tenes::real_tensor bottom =
      mptensor::tensordot(r3_reduced_site(2, 2, 2, number_site == 2, seed),
                          r3_reduced_site(2, 2, 3, number_site == 3, seed),
                          mptensor::Axes(2), mptensor::Axes(0));
  return r3_sum_entries(mptensor::tensordot(top, bottom, mptensor::Axes(2, 5),
                                            mptensor::Axes(1, 3)));
}

// Since 2026-08-28 the blob builders take operators in the single-layer
// convention (wrap_twosite_gate: INPUT swap only). The r2_*_op helpers
// build the plain matrix elements as a verbatim ftensor, so the pair-norm
// helpers below apply that input swap before handing the operator to the
// blob. The reference values are physical and unchanged.
static ft r3_blob_op(const ft& op) {
  ft o = op;
  tenes::fermion::apply_swap(o, 0, 1);
  return o;
}

static double r3_horizontal_pair_norm(const ft& op, int seed) {
  return r3_sum_entries(tenes::fermion::build_reduced_pair_naive(
      make_r2_tensor(2, 1, 0, seed), make_r2_tensor(2, 1, 1, seed),
      r3_blob_op(op), tenes::fermion::reduced_pair_direction::horizontal));
}

static double r3_vertical_pair_norm(const ft& op, int seed) {
  return r3_sum_entries(tenes::fermion::build_reduced_pair_naive(
      make_r2_tensor(1, 2, 0, seed), make_r2_tensor(1, 2, 1, seed),
      r3_blob_op(op), tenes::fermion::reduced_pair_direction::vertical));
}

static double r3_pair_plaquette_norm(int a, int b, const ft& op, int seed) {
  const ft o = r3_blob_op(op);
  if (a == 0 && b == 1) {
    tenes::real_tensor blob = tenes::fermion::build_reduced_pair_naive(
        make_r2_tensor(2, 2, 0, seed), make_r2_tensor(2, 2, 1, seed), o,
        tenes::fermion::reduced_pair_direction::horizontal);
    tenes::real_tensor bottom =
        mptensor::tensordot(r3_reduced_site(2, 2, 2, false, seed),
                            r3_reduced_site(2, 2, 3, false, seed),
                            mptensor::Axes(2), mptensor::Axes(0));
    return r3_sum_entries(mptensor::tensordot(
        blob, bottom, mptensor::Axes(2, 5), mptensor::Axes(1, 3)));
  }
  if (a == 2 && b == 3) {
    tenes::real_tensor top =
        mptensor::tensordot(r3_reduced_site(2, 2, 0, false, seed),
                            r3_reduced_site(2, 2, 1, false, seed),
                            mptensor::Axes(2), mptensor::Axes(0));
    tenes::real_tensor blob = tenes::fermion::build_reduced_pair_naive(
        make_r2_tensor(2, 2, 2, seed), make_r2_tensor(2, 2, 3, seed), o,
        tenes::fermion::reduced_pair_direction::horizontal);
    return r3_sum_entries(mptensor::tensordot(top, blob, mptensor::Axes(2, 5),
                                              mptensor::Axes(1, 3)));
  }
  if (a == 0 && b == 2) {
    tenes::real_tensor blob = tenes::fermion::build_reduced_pair_naive(
        make_r2_tensor(2, 2, 0, seed), make_r2_tensor(2, 2, 2, seed), o,
        tenes::fermion::reduced_pair_direction::vertical);
    tenes::real_tensor right =
        mptensor::tensordot(r3_reduced_site(2, 2, 1, false, seed),
                            r3_reduced_site(2, 2, 3, false, seed),
                            mptensor::Axes(3), mptensor::Axes(1));
    return r3_sum_entries(mptensor::tensordot(blob, right, mptensor::Axes(2, 4),
                                              mptensor::Axes(0, 3)));
  }
  tenes::real_tensor left =
      mptensor::tensordot(r3_reduced_site(2, 2, 0, false, seed),
                          r3_reduced_site(2, 2, 2, false, seed),
                          mptensor::Axes(3), mptensor::Axes(1));
  tenes::real_tensor blob = tenes::fermion::build_reduced_pair_naive(
      make_r2_tensor(2, 2, 1, seed), make_r2_tensor(2, 2, 3, seed), o,
      tenes::fermion::reduced_pair_direction::vertical);
  return r3_sum_entries(mptensor::tensordot(left, blob, mptensor::Axes(2, 4),
                                            mptensor::Axes(0, 3)));
}

static double r3_pair_plaquette_norm(int a, int b, int seed) {
  return r3_pair_plaquette_norm(a, b, r2_hop_op(), seed);
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

TEST_CASE("R3 diagnostic split density pair blob matches site impurities") {
  const ft n_a = r2_density_a_op();
  const ft n_b = r2_density_b_op();

  for (const auto& ref : r3_refs) {
    const double h_norm = r3_horizontal_norm(false, false, ref.seed);
    const double h_a = r3_horizontal_pair_norm(n_a, ref.seed) / h_norm;
    const double h_b = r3_horizontal_pair_norm(n_b, ref.seed) / h_norm;
    std::cout << "R3 split density horizontal seed=" << ref.seed << " A=" << h_a
              << " refA=" << ref.h_n0 << " B=" << h_b << " refB=" << ref.h_n1
              << std::endl;
    CHECK(h_a == doctest::Approx(ref.h_n0).epsilon(1e-12));
    CHECK(h_b == doctest::Approx(ref.h_n1).epsilon(1e-12));

    const double v_norm = r3_vertical_norm(false, false, ref.seed);
    const double v_a = r3_vertical_pair_norm(n_a, ref.seed) / v_norm;
    const double v_b = r3_vertical_pair_norm(n_b, ref.seed) / v_norm;
    std::cout << "R3 split density vertical seed=" << ref.seed << " A=" << v_a
              << " refA=" << ref.v_n0 << " B=" << v_b << " refB=" << ref.v_n1
              << std::endl;
    CHECK(v_a == doctest::Approx(ref.v_n0).epsilon(1e-12));
    CHECK(v_b == doctest::Approx(ref.v_n1).epsilon(1e-12));
  }
}

TEST_CASE("R3 reduced convention agrees with oracle on both datasets") {
  const int bonds[][2] = {{0, 1}, {0, 2}, {1, 3}, {2, 3}};
  const ft density_op = r2_density_pair_op();
  const ft hop_op = r2_hop_op();
  const ft mixed_op = r2_mixed_op();
  for (const auto& ref : r3_refs) {
    const double h_norm = r3_horizontal_norm(false, false, ref.seed);
    CHECK(h_norm == doctest::Approx(ref.h_norm).epsilon(1e-12));
    CHECK(r3_horizontal_norm(true, false, ref.seed) / h_norm ==
          doctest::Approx(ref.h_n0).epsilon(1e-12));
    CHECK(r3_horizontal_norm(false, true, ref.seed) / h_norm ==
          doctest::Approx(ref.h_n1).epsilon(1e-12));
    CHECK(r3_horizontal_pair_norm(hop_op, ref.seed) / h_norm ==
          doctest::Approx(0.0).epsilon(1e-12));
    CHECK(r3_horizontal_pair_norm(density_op, ref.seed) / h_norm ==
          doctest::Approx(ref.h_n0 + ref.h_n1).epsilon(1e-12));
    CHECK(r3_horizontal_pair_norm(mixed_op, ref.seed) / h_norm ==
          doctest::Approx(ref.h_mixed).epsilon(1e-12));

    const double v_norm = r3_vertical_norm(false, false, ref.seed);
    CHECK(v_norm == doctest::Approx(ref.v_norm).epsilon(1e-12));
    CHECK(r3_vertical_norm(true, false, ref.seed) / v_norm ==
          doctest::Approx(ref.v_n0).epsilon(1e-12));
    CHECK(r3_vertical_norm(false, true, ref.seed) / v_norm ==
          doctest::Approx(ref.v_n1).epsilon(1e-12));
    CHECK(r3_vertical_pair_norm(hop_op, ref.seed) / v_norm ==
          doctest::Approx(0.0).epsilon(1e-12));
    CHECK(r3_vertical_pair_norm(density_op, ref.seed) / v_norm ==
          doctest::Approx(ref.v_n0 + ref.v_n1).epsilon(1e-12));
    CHECK(r3_vertical_pair_norm(mixed_op, ref.seed) / v_norm ==
          doctest::Approx(ref.v_mixed).epsilon(1e-12));

    const double p_norm = r3_plaquette_norm(-1, ref.seed);
    CHECK(p_norm == doctest::Approx(ref.p_norm).epsilon(1e-12));
    for (int site = 0; site < 4; ++site) {
      CHECK(r3_plaquette_norm(site, ref.seed) / p_norm ==
            doctest::Approx(ref.p_n[site]).epsilon(1e-12));
    }
    for (int i = 0; i < 4; ++i) {
      const int a = bonds[i][0];
      const int b = bonds[i][1];
      CHECK(r3_pair_plaquette_norm(a, b, ref.seed) / p_norm ==
            doctest::Approx(ref.p_hop[i]).epsilon(1e-12));
      CHECK(r3_pair_plaquette_norm(a, b, density_op, ref.seed) / p_norm ==
            doctest::Approx(ref.p_n[a] + ref.p_n[b]).epsilon(1e-12));
      CHECK(r3_pair_plaquette_norm(a, b, mixed_op, ref.seed) / p_norm ==
            doctest::Approx(ref.p_mixed[i]).epsilon(1e-12));
    }
  }
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

// ===== R5: d = 4 (electron) oracle pinning of the reduced-pair blob ========
//
// Since 2026-08-28 the blob builders take operators in the single-layer
// convention: wrap_twosite_gate, i.e. the INPUT swap only. R5 pins, on
// d = 4 two-site shapes, that a blob loaded this way reproduces the
// Fock-verified direct f-primitive expectation (the reference is the
// direct graded expectation with the same input-swapped operator, verified
// against exact Fock evolution at d = 4 in the plaquette and chain
// diagnostics). d = 4 is what makes the loading convention observable at
// all: the electron hopping has (odd,odd) input and output channels at
// linear order, so both mis-loadings -- verbatim (no swap) and the retired
// blob convention (in+out swaps) -- give visibly different hopping values,
// which the negative checks below keep distinguishable.

namespace {

tenes::fermion::leg_parities r4_parities(int lx, int ly, int site) {
  const int x = site % lx;
  const int y = site / lx;
  tenes::fermion::parity_vector even{false};
  tenes::fermion::parity_vector bond{false, true};
  tenes::fermion::parity_vector phys{false, true, true, false};
  return {x > 0 ? bond : even, y > 0 ? bond : even, x + 1 < lx ? bond : even,
          y + 1 < ly ? bond : even, phys};
}

ft make_r4_tensor(int lx, int ly, int site, int seed = 0) {
  const auto p = r4_parities(lx, ly, site);
  tenes::real_tensor t(mptensor::Shape(p[0].size(), p[1].size(), p[2].size(),
                                       p[3].size(), p[4].size()));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    auto idx = t.global_index(n);
    if (tenes::fermion::count_odd(p, idx) % 2 == 0) {
      const double x = static_cast<double>(
          (site + 2) * (1 + seed + (3 + seed % 5) * idx[0] +
                        (4 + seed % 5) * idx[1] + (5 + seed % 5) * idx[2] +
                        (6 + seed % 5) * idx[3] + (11 + seed % 7) * idx[4]));
      t.set_value(idx, 0.19 * std::sin(x) + 0.13 * std::cos(0.41 * x));
    }
  }
  return ft{t, p};
}

const tenes::fermion::parity_vector r4_phys{false, true, true, false};

ft r4_number_op() {
  ft op{tenes::real_tensor(mptensor::Shape(4, 4)), {r4_phys, r4_phys}};
  op.t.set_value(mptensor::Index(1, 1), 1.0);
  op.t.set_value(mptensor::Index(2, 2), 1.0);
  op.t.set_value(mptensor::Index(3, 3), 2.0);
  return op;
}

// plain rank-4 tensors (no swaps applied)
tenes::real_tensor r4_nn_plain() {
  tenes::real_tensor op(mptensor::Shape(4, 4, 4, 4));
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

tenes::real_tensor r4_hop_plain() {
  const auto h = electron_bond_hamiltonian(1.0, 0.0, 0.0);
  tenes::real_tensor op(mptensor::Shape(4, 4, 4, 4));
  for (int i1 = 0; i1 < 4; ++i1) {
    for (int i2 = 0; i2 < 4; ++i2) {
      for (int o1 = 0; o1 < 4; ++o1) {
        for (int o2 = 0; o2 < 4; ++o2) {
          const double v = h[o1 * 4 + o2][i1 * 4 + i2];
          if (std::abs(v) > 1.0e-15) {
            op.set_value(mptensor::Index(i1, i2, o1, o2), v);
          }
        }
      }
    }
  }
  return op;
}

ft r4_wrap(const tenes::real_tensor& plain, bool swap_in, bool swap_out) {
  ft op{plain, {r4_phys, r4_phys, r4_phys, r4_phys}};
  if (swap_in) {
    tenes::fermion::apply_swap(op, 0, 1);
  }
  if (swap_out) {
    tenes::fermion::apply_swap(op, 2, 3);
  }
  return op;
}

double r4_pair_norm_blob(const ft& op, bool horizontal, int seed) {
  if (horizontal) {
    return r3_sum_entries(tenes::fermion::build_reduced_pair_naive(
        make_r4_tensor(2, 1, 0, seed), make_r4_tensor(2, 1, 1, seed), op,
        tenes::fermion::reduced_pair_direction::horizontal));
  }
  return r3_sum_entries(tenes::fermion::build_reduced_pair_naive(
      make_r4_tensor(1, 2, 0, seed), make_r4_tensor(1, 2, 1, seed), op,
      tenes::fermion::reduced_pair_direction::vertical));
}

double r4_norm_blob(bool horizontal, int seed) {
  auto site = [&](int s) {
    return horizontal ? make_r4_tensor(2, 1, s, seed)
                      : make_r4_tensor(1, 2, s, seed);
  };
  if (horizontal) {
    return r3_sum_entries(
        mptensor::tensordot(tenes::fermion::build_reduced(site(0)),
                            tenes::fermion::build_reduced(site(1)),
                            mptensor::Axes(2), mptensor::Axes(0)));
  }
  return r3_sum_entries(
      mptensor::tensordot(tenes::fermion::build_reduced(site(0)),
                          tenes::fermion::build_reduced(site(1)),
                          mptensor::Axes(3), mptensor::Axes(1)));
}

}  // namespace

TEST_CASE("R5 d=4 reduced pair blob against the Fock-verified direct path") {
  for (const int seed : {0, 173}) {
    for (const bool horizontal : {true, false}) {
      const ft psi =
          horizontal
              ? tenes::fermion::tensordot(make_r4_tensor(2, 1, 0, seed),
                                          make_r4_tensor(2, 1, 1, seed),
                                          mptensor::Axes(2), mptensor::Axes(0))
              : tenes::fermion::tensordot(make_r4_tensor(1, 2, 0, seed),
                                          make_r4_tensor(1, 2, 1, seed),
                                          mptensor::Axes(3), mptensor::Axes(1));
      const double norm = r2_norm(psi);
      const std::string label =
          std::string(horizontal ? "H" : "V") + " seed=" + std::to_string(seed);

      // layer 1: norm
      const double norm_blob = r4_norm_blob(horizontal, seed);
      CHECK(norm_blob / norm == doctest::Approx(1.0).epsilon(1e-10));

      // layer 2: one-site number density via composed one-site ops
      const ft n_op = r4_number_op();
      const double n0_ref = r2_expect_one(psi, 3, n_op);
      const double n1_ref = r2_expect_one(psi, 7, n_op);

      // layer 3: rank-4 diagonal nn (odd-odd diagonal channel)
      const double nn_ref = [&] {
        ft tmp = r2_apply_one_site(psi, 7, n_op);
        tmp = r2_apply_one_site(tmp, 3, n_op);
        mptensor::Axes all;
        for (int ax = 0; ax < psi.rank(); ++ax) {
          all.push(ax);
        }
        return tenes::fermion::trace(tenes::fermion::conj(psi), tmp, all, all) /
               norm;
      }();
      // cross-check the reference route: canonical (in-swapped) two-site
      const double nn_direct =
          r2_expect_two(psi, 3, 7, r4_wrap(r4_nn_plain(), true, false));
      CHECK(nn_direct == doctest::Approx(nn_ref).epsilon(1e-10));

      // layer 4: rank-4 electron hopping ((odd,odd)->(even,even) channel)
      const double hop_ref =
          r2_expect_two(psi, 3, 7, r4_wrap(r4_hop_plain(), true, false));

      // blob candidates: the production loading (input swap only) and the
      // two mis-loadings the d = 4 hopping can tell apart
      const double nn_in =
          r4_pair_norm_blob(r4_wrap(r4_nn_plain(), true, false), horizontal,
                            seed) /
          norm_blob;
      const double hop_in =
          r4_pair_norm_blob(r4_wrap(r4_hop_plain(), true, false), horizontal,
                            seed) /
          norm_blob;
      const double hop_plain =
          r4_pair_norm_blob(r4_wrap(r4_hop_plain(), false, false), horizontal,
                            seed) /
          norm_blob;
      const double hop_both =
          r4_pair_norm_blob(r4_wrap(r4_hop_plain(), true, true), horizontal,
                            seed) /
          norm_blob;

      std::cout << "R5 " << label << " n_ref=(" << n0_ref << ", " << n1_ref
                << ") nn_ref=" << nn_ref << " nn_in=" << nn_in
                << " hop_ref=" << hop_ref << " hop_in=" << hop_in
                << " hop_plain=" << hop_plain << " hop_both=" << hop_both
                << std::endl;

      // The production blob convention (wrap_twosite_gate = input swap
      // only) must reproduce the Fock-verified direct path; verbatim
      // loading and the retired in+out-swap convention are both wrong on
      // the hopping.
      CHECK(nn_in == doctest::Approx(nn_ref).epsilon(1e-10));
      CHECK(hop_in == doctest::Approx(hop_ref).epsilon(1e-10));
      const double hop_scale = std::max(1.0e-6, std::abs(hop_ref));
      CHECK(std::abs(hop_plain - hop_ref) > 0.1 * hop_scale);
      CHECK(std::abs(hop_both - hop_ref) > 0.1 * hop_scale);
    }
  }
}
