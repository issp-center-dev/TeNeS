namespace {

static tenes::fermion::leg_parities r2_parities(int lx, int ly, int site) {
  const int x = site % lx;
  const int y = site / lx;
  tenes::fermion::parity_vector even{false};
  tenes::fermion::parity_vector p{false, true};
  return {x > 0 ? p : even, y > 0 ? p : even, x + 1 < lx ? p : even,
          y + 1 < ly ? p : even, p};
}

static ft make_r2_tensor(int lx, int ly, int site) {
  const auto p = r2_parities(lx, ly, site);
  tenes::real_tensor t(mptensor::Shape(p[0].size(), p[1].size(), p[2].size(),
                                       p[3].size(), p[4].size()));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    auto idx = t.global_index(n);
    if (tenes::fermion::count_odd(p, idx) % 2 == 0) {
      const double x = static_cast<double>(
          (site + 2) *
          (1 + 3 * idx[0] + 4 * idx[1] + 5 * idx[2] + 6 * idx[3] + 7 * idx[4]));
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
