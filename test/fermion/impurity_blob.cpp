// Tests for the impurity factorization of the fermion two-site blob:
// build_reduced_pair_direct against the build_reduced_pair_naive oracle.
// Contract: misc/fermion_twosite_blob/contract.md. The existing
// build_reduced_pair_naive is the oracle; all helpers live in this file's own
// anonymous namespace (this file is included into the test_fermion_layer TU).

namespace {

// ---- basic parity vectors -------------------------------------------------

tenes::fermion::parity_vector ib_phys_parity(int d) {
  if (d == 2) {
    return {false, true};
  }
  if (d == 4) {
    return {false, true, true, false};
  }
  throw std::runtime_error("ib_phys_parity: unsupported physical dimension");
}

tenes::fermion::parity_vector ib_vp_even(int D) {
  return tenes::fermion::parity_vector(D, false);
}

tenes::fermion::parity_vector ib_vp_odd(int D) {
  return tenes::fermion::parity_vector(D, true);
}

tenes::fermion::parity_vector ib_vp_mixed(int D) {
  tenes::fermion::parity_vector v(D, false);
  for (int i = 1; i < D; i += 2) {
    v[i] = true;
  }
  return v;
}

tenes::fermion::parity_vector ib_vp(int D, const std::string& kind) {
  if (kind == "even") {
    return ib_vp_even(D);
  }
  if (kind == "odd") {
    return ib_vp_odd(D);
  }
  if (kind == "mixed") {
    return ib_vp_mixed(D);
  }
  throw std::runtime_error("ib_vp: unknown parity kind");
}

// ---- deterministic pseudo-random fills ------------------------------------

double ib_random_value(int seed, const mptensor::Index& idx, std::size_t rank) {
  const double w[5] = {3.1, 5.3, 7.9, 11.7, 13.1};
  double x = 0.7 * (seed + 1) + 0.013 * seed * seed;
  for (std::size_t ax = 0; ax < rank; ++ax) {
    x += w[ax] * static_cast<double>(idx[ax]);
  }
  return 0.5 * std::sin(x) + 0.25 * std::cos(1.7 * x + 0.3 * seed);
}

// Random dense site tensor, legs (l, t, r, b, s); all four virtual legs share
// the same parity vector. Parity conservation is deliberately NOT imposed
// (contract section 2: the equivalence tests pin with dense random tensors).
ft ib_make_site_tensor(int d, const tenes::fermion::parity_vector& vp,
                       int seed) {
  const tenes::fermion::leg_parities p{vp, vp, vp, vp, ib_phys_parity(d)};
  tenes::real_tensor t(
      mptensor::Shape(vp.size(), vp.size(), vp.size(), vp.size(), d));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    t.set_value(idx, ib_random_value(seed, idx, 5));
  }
  return ft{t, p};
}

// ---- two-site operators (plain, leg order (in_A, in_B, out_A, out_B)) -----

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

// c_B c_A + c^dag_A c^dag_B, element table taken verbatim from the contract.
tenes::real_tensor ib_pairing_plain() {
  tenes::real_tensor op(mptensor::Shape(2, 2, 2, 2));
  op.set_value(mptensor::Index(1, 1, 0, 0), 1.0);
  op.set_value(mptensor::Index(0, 0, 1, 1), 1.0);
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

ft ib_wrap(const tenes::real_tensor& plain, int d) {
  return tenes::fermion::wrap_reduced_pair_op(plain, ib_phys_parity(d),
                                              ib_phys_parity(d));
}

// Wrapped operator; source_second applies the graded transpose used when the
// operator's first leg sits on the window's second site (twosite_obs.cpp).
ft ib_op12(const std::string& kind, int d, int seed, bool source_second) {
  tenes::real_tensor plain;
  if (kind == "identity") {
    plain = ib_identity_plain(d);
  } else if (kind == "hopping") {
    plain = ib_hop_plain(d);
  } else if (kind == "nn") {
    plain = ib_nn_plain(d);
  } else if (kind == "pairing") {
    plain = ib_pairing_plain();
  } else if (kind == "random") {
    plain = ib_random_even_plain(d, seed);
  } else {
    throw std::runtime_error("ib_op12: unknown operator kind");
  }
  ft op = ib_wrap(plain, d);
  if (source_second) {
    op = tenes::fermion::transpose(op, mptensor::Axes(1, 0, 3, 2));
  }
  return op;
}

// ---- elementwise comparisons ----------------------------------------------

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

double ib_max_abs_entry(const tenes::real_tensor& a) {
  double m = 0.0;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    double v;
    a.get_value(a.global_index(n), v);
    m = std::max(m, std::abs(v));
  }
  return m;
}

// |got - want| <= 1e-12 * max(1, max|want|), scanning every element
// (contract section 3; norm- or sum-based comparison is forbidden).
void ib_check_allclose(const tenes::real_tensor& got,
                       const tenes::real_tensor& want,
                       const std::string& label) {
  REQUIRE(got.shape() == want.shape());
  const double tol = 1.0e-12 * std::max(1.0, ib_max_abs_entry(want));
  double max_dev = 0.0;
  mptensor::Index worst;
  for (std::size_t n = 0; n < want.local_size(); ++n) {
    const mptensor::Index idx = want.global_index(n);
    double w, g;
    want.get_value(idx, w);
    got.get_value(idx, g);
    const double dev = std::abs(g - w);
    if (dev > max_dev) {
      max_dev = dev;
      worst = idx;
    }
  }
  INFO(label << ": max |new-old| = " << max_dev << " (tol " << tol
             << ") at index " << ib_index_string(worst));
  CHECK(max_dev <= tol);
}

// ---- T3 equivalence runner ------------------------------------------------

void ib_check_pair_equivalence(const std::string& label, int d,
                               const tenes::fermion::parity_vector& vp,
                               const ft& op12,
                               tenes::fermion::reduced_pair_direction dir,
                               int seedA, int seedB) {
  INFO(label << " seedA=" << seedA << " seedB=" << seedB);
  const ft TnA = ib_make_site_tensor(d, vp, seedA);
  const ft TnB = ib_make_site_tensor(d, vp, seedB);
  const tenes::real_tensor blob_old =
      tenes::fermion::build_reduced_pair_naive(TnA, TnB, op12, dir);
  const tenes::real_tensor blob_new =
      tenes::fermion::build_reduced_pair_direct(TnA, TnB, op12, dir);
  ib_check_allclose(blob_new, blob_old,
                    "factorized blob vs build_reduced_pair_naive");
}

struct ib_sweep_case {
  const char* op;
  int d;
  int D;
  const char* parity_kind;
  char dir;  // 'h' or 'v'
  bool source_second;
  int op_seed;
};

tenes::fermion::reduced_pair_direction ib_dir(char c) {
  return c == 'h' ? tenes::fermion::reduced_pair_direction::horizontal
                  : tenes::fermion::reduced_pair_direction::vertical;
}

std::string ib_sweep_label(const ib_sweep_case& sc) {
  std::string label =
      std::string("op=") + sc.op + " d=" + std::to_string(sc.d) +
      " D=" + std::to_string(sc.D) + " parity=" + sc.parity_kind +
      " dir=" + sc.dir + " source=" + (sc.source_second ? "second" : "first");
  if (std::string(sc.op) == "random") {
    label += " op_seed=" + std::to_string(sc.op_seed);
  }
  return label;
}

void ib_run_sweep_case(const ib_sweep_case& sc, int case_id) {
  ib_check_pair_equivalence(ib_sweep_label(sc), sc.d,
                            ib_vp(sc.D, sc.parity_kind),
                            ib_op12(sc.op, sc.d, sc.op_seed, sc.source_second),
                            ib_dir(sc.dir), 41 + case_id, 71 + case_id);
}

}  // namespace

// ===== T3-1: equivalence against the oracle ================================

TEST_CASE(
    "impurity blob T3-1: mandatory hopping sweep, direction x source x (d,D) "
    "with mixed parity") {
  // Contract section 6: this block must not be thinned; it is the detector
  // for the odd k sector combined with bond-crossing signs.
  const std::pair<int, int> dims[] = {{2, 2}, {2, 3}, {4, 2}};
  int case_id = 0;
  for (const auto& dD : dims) {
    for (const char dir : {'h', 'v'}) {
      for (const bool second : {false, true}) {
        ++case_id;
        const ib_sweep_case sc{"hopping", dD.first, dD.second, "mixed",
                               dir,       second,   0};
        ib_run_sweep_case(sc, case_id);
      }
    }
  }
}

TEST_CASE("impurity blob T3-1: parity and operator sweeps") {
  // Representatives along the remaining axes (fixed seeds); the mandatory
  // combinations live in the dedicated test case above.
  const ib_sweep_case cases[] = {
      // virtual-parity axis, hopping, source first
      {"hopping", 2, 1, "even", 'h', false, 0},
      {"hopping", 4, 1, "even", 'v', false, 0},
      {"hopping", 2, 2, "even", 'h', false, 0},
      {"hopping", 2, 2, "odd", 'v', false, 0},
      {"hopping", 2, 3, "even", 'v', false, 0},
      {"hopping", 2, 3, "odd", 'h', false, 0},
      {"hopping", 4, 2, "even", 'h', false, 0},
      {"hopping", 4, 2, "odd", 'v', false, 0},
      // operator axis at (d=2, D=2, mixed)
      {"identity", 2, 2, "mixed", 'h', false, 0},
      {"identity", 2, 2, "mixed", 'v', false, 0},
      {"nn", 2, 2, "mixed", 'h', false, 0},
      {"nn", 2, 2, "mixed", 'v', false, 0},
      {"pairing", 2, 2, "mixed", 'h', false, 0},
      {"pairing", 2, 2, "mixed", 'v', false, 0},
      {"random", 2, 2, "mixed", 'h', false, 21},
      {"random", 2, 2, "mixed", 'v', false, 21},
      {"random", 2, 2, "mixed", 'h', true, 22},
      {"random", 2, 2, "mixed", 'v', true, 22},
      // operator axis at (d=4, D=2, mixed)
      {"identity", 4, 2, "mixed", 'h', false, 0},
      {"identity", 4, 2, "mixed", 'v', false, 0},
      {"nn", 4, 2, "mixed", 'h', false, 0},
      {"nn", 4, 2, "mixed", 'v', false, 0},
      {"random", 4, 2, "mixed", 'h', false, 21},
      {"random", 4, 2, "mixed", 'v', false, 21},
      {"random", 4, 2, "mixed", 'h', true, 22},
      {"random", 4, 2, "mixed", 'v', true, 22},
      // D=1 operator variety
      {"identity", 2, 1, "even", 'h', false, 0},
      {"random", 2, 1, "even", 'v', false, 24},
      {"nn", 4, 1, "even", 'h', false, 0},
      // D=3 with a transposed generic operator
      {"random", 2, 3, "mixed", 'v', true, 23},
  };
  int case_id = 100;
  for (const auto& sc : cases) {
    ++case_id;
    ib_run_sweep_case(sc, case_id);
  }
}

// ===== T3-2: slow d=4, D=3 (env gated) =====================================

TEST_CASE(
    "impurity blob T3-2: slow d=4 D=3 equivalence "
    "(TENES_RUN_IMPURITY_BLOB_SLOW)") {
  if (std::getenv("TENES_RUN_IMPURITY_BLOB_SLOW") == nullptr) {
    return;
  }
  const ib_sweep_case cases[] = {
      {"hopping", 4, 3, "mixed", 'h', false, 0},
      {"hopping", 4, 3, "mixed", 'v', false, 0},
      {"random", 4, 3, "mixed", 'h', true, 91},
      {"random", 4, 3, "mixed", 'v', true, 91},
  };
  int case_id = 200;
  for (const auto& sc : cases) {
    ++case_id;
    ib_run_sweep_case(sc, case_id);
  }
}

// ===== T3-1 (addendum 2): complex-tensor small equivalence cases ==========

namespace {

using cft = tenes::fermion::ftensor<tenes::complex_tensor>;

std::complex<double> ibc_random_value(int seed, const mptensor::Index& idx,
                                      std::size_t rank) {
  return {ib_random_value(seed, idx, rank),
          ib_random_value(seed + 1000, idx, rank)};
}

// Complex twin of ib_make_site_tensor: dense deterministic pseudo-random
// fill in both the real and the imaginary part (contract addendum 2).
cft ibc_make_site_tensor(int d, const tenes::fermion::parity_vector& vp,
                         int seed) {
  const tenes::fermion::leg_parities p{vp, vp, vp, vp, ib_phys_parity(d)};
  tenes::complex_tensor t(
      mptensor::Shape(vp.size(), vp.size(), vp.size(), vp.size(), d));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const mptensor::Index idx = t.global_index(n);
    t.set_value(idx, ibc_random_value(seed, idx, 5));
  }
  return cft{t, p};
}

// hopping: matrix elements stay real, only the tensor type is complex.
tenes::complex_tensor ibc_hop_plain(int d) {
  const tenes::real_tensor re = ib_hop_plain(d);
  tenes::complex_tensor op(mptensor::Shape(d, d, d, d));
  for (std::size_t n = 0; n < re.local_size(); ++n) {
    const mptensor::Index idx = re.global_index(n);
    double v;
    re.get_value(idx, v);
    op.set_value(idx, std::complex<double>(v, 0.0));
  }
  return op;
}

// random total-parity-even operator with complex elements.
tenes::complex_tensor ibc_random_even_plain(int d, int seed) {
  const tenes::fermion::parity_vector phys = ib_phys_parity(d);
  tenes::complex_tensor op(mptensor::Shape(d, d, d, d));
  for (std::size_t n = 0; n < op.local_size(); ++n) {
    const mptensor::Index idx = op.global_index(n);
    const bool odd =
        ((phys[idx[0]] != phys[idx[1]]) != (phys[idx[2]] != phys[idx[3]]));
    if (!odd) {
      op.set_value(idx, ibc_random_value(seed, idx, 4));
    }
  }
  return op;
}

cft ibc_op12(const std::string& kind, int d, int seed, bool source_second) {
  tenes::complex_tensor plain;
  if (kind == "hopping") {
    plain = ibc_hop_plain(d);
  } else if (kind == "random") {
    plain = ibc_random_even_plain(d, seed);
  } else {
    throw std::runtime_error("ibc_op12: unknown operator kind");
  }
  cft op = tenes::fermion::wrap_reduced_pair_op(plain, ib_phys_parity(d),
                                                ib_phys_parity(d));
  if (source_second) {
    op = tenes::fermion::transpose(op, mptensor::Axes(1, 0, 3, 2));
  }
  return op;
}

// |new - old| <= 1e-12 * max(1, max|old|) on the complex modulus, scanning
// every element (contract addendum 2 / section 3).
void ibc_check_allclose(const tenes::complex_tensor& got,
                        const tenes::complex_tensor& want,
                        const std::string& label) {
  REQUIRE(got.shape() == want.shape());
  double scale = 1.0;
  for (std::size_t n = 0; n < want.local_size(); ++n) {
    std::complex<double> v;
    want.get_value(want.global_index(n), v);
    scale = std::max(scale, std::abs(v));
  }
  const double tol = 1.0e-12 * scale;
  double max_dev = 0.0;
  mptensor::Index worst;
  for (std::size_t n = 0; n < want.local_size(); ++n) {
    const mptensor::Index idx = want.global_index(n);
    std::complex<double> w, g;
    want.get_value(idx, w);
    got.get_value(idx, g);
    const double dev = std::abs(g - w);
    if (dev > max_dev) {
      max_dev = dev;
      worst = idx;
    }
  }
  INFO(label << ": max |new-old| = " << max_dev << " (tol " << tol
             << ") at index " << ib_index_string(worst));
  CHECK(max_dev <= tol);
}

void ibc_check_pair_equivalence(const std::string& label, int d,
                                const tenes::fermion::parity_vector& vp,
                                const cft& op12,
                                tenes::fermion::reduced_pair_direction dir,
                                int seedA, int seedB) {
  INFO(label << " seedA=" << seedA << " seedB=" << seedB);
  const cft TnA = ibc_make_site_tensor(d, vp, seedA);
  const cft TnB = ibc_make_site_tensor(d, vp, seedB);
  const tenes::complex_tensor blob_old =
      tenes::fermion::build_reduced_pair_naive(TnA, TnB, op12, dir);
  const tenes::complex_tensor blob_new =
      tenes::fermion::build_reduced_pair_direct(TnA, TnB, op12, dir);
  ibc_check_allclose(blob_new, blob_old,
                     "factorized blob vs build_reduced_pair_naive (complex)");
}

}  // namespace

TEST_CASE(
    "impurity blob T3-1 complex: factorized pair blob equals oracle at d=2 "
    "D=2 mixed parity") {
  // Contract addendum 2: five small complex_tensor equivalence cases.
  struct ibc_case {
    const char* op;
    char dir;  // 'h' or 'v'
    bool source_second;
    int op_seed;
  };
  const ibc_case cases[] = {
      {"hopping", 'h', false, 0}, {"hopping", 'v', false, 0},
      {"random", 'h', false, 51}, {"random", 'v', false, 51},
      {"random", 'h', true, 52},
  };
  const tenes::fermion::parity_vector vp = ib_vp_mixed(2);
  int case_id = 300;
  for (const auto& c : cases) {
    ++case_id;
    std::string label = std::string("complex op=") + c.op + " d=2 D=2 " +
                        "parity=mixed dir=" + c.dir +
                        " source=" + (c.source_second ? "second" : "first");
    if (std::string(c.op) == "random") {
      label += " op_seed=" + std::to_string(c.op_seed);
    }
    ibc_check_pair_equivalence(label, 2, vp,
                               ibc_op12(c.op, 2, c.op_seed, c.source_second),
                               ib_dir(c.dir), 41 + case_id, 71 + case_id);
  }
}

// T3-3 (dressing pin) was withdrawn by contract addendum 1, which also
// retired the split_pair_op / build_reduced_impurity APIs and their T1/T2
// test cases from this file.

// ===== addendum 3: traced-leg parity mismatch guard ========================

TEST_CASE(
    "impurity blob guard: fuse_doubled_cluster_direct rejects traced-leg "
    "parity mismatch") {
  // Contract addendum 3: the guard is structurally unreachable from
  // production inputs, so it is exercised directly. D=1, d=2, all-zero
  // values; horizontal pair-state leg order (l_A,t_A,b_A,s_A,t_B,r_B,b_B,s_B).
  const tenes::fermion::parity_vector even1{false};
  const tenes::fermion::parity_vector phys{false, true};
  const tenes::fermion::parity_vector phys_flipped{true, false};
  const tenes::real_tensor zero(mptensor::Shape(1, 1, 1, 2, 1, 1, 1, 2));
  const std::vector<int> leg_ids{0, 1, 3, 1, 2, 3};

  const ft bra{zero, {even1, even1, even1, phys, even1, even1, even1, phys}};

  // matched parities: the guard must stay silent (return value not checked)
  const ft ket_ok = bra;
  CHECK_NOTHROW(tenes::fermion::detail::fuse_doubled_cluster_direct(bra, ket_ok,
                                                                    leg_ids));

  // ket trace leg 3 disagrees with the bra twin
  const ft ket_bad3{
      zero, {even1, even1, even1, phys_flipped, even1, even1, even1, phys}};
  CHECK_THROWS_WITH_AS(
      tenes::fermion::detail::fuse_doubled_cluster_direct(bra, ket_bad3,
                                                          leg_ids),
      "fuse_doubled_cluster_direct: traced-leg parity mismatch",
      std::runtime_error);

  // ket trace leg 7 disagrees with the bra twin
  const ft ket_bad7{
      zero, {even1, even1, even1, phys, even1, even1, even1, phys_flipped}};
  CHECK_THROWS_WITH_AS(
      tenes::fermion::detail::fuse_doubled_cluster_direct(bra, ket_bad7,
                                                          leg_ids),
      "fuse_doubled_cluster_direct: traced-leg parity mismatch",
      std::runtime_error);
}
