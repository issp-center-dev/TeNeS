// ===== A bond described from either end must give the same update ==========
//
// A vertical bond between a site s and t = top(s) can be handed to
// iTPS<tensor>::simple_update in two ways:
//
//   (source = s, source_leg = 1)   -- s is named first
//   (source = t, source_leg = 3)   -- t is named first
//
// They are the same bond in the same geometry acting on the same state, so
// they must produce the same result. They do not take the same code path: the
// driver normalises a bond given from the JW-later side, and its condition is
// `source_leg == 0 || source_leg == 1`, so only the first form is rewritten.
//
// Rewriting means expressing the two-site operator with the site roles
// exchanged, which is a graded operation: conjugating by
// S(|a> (x) |b>) = (-1)^{|a||b|} |b> (x) |a> gives
//
//   O'[j1,j2,k1,k2] = (-1)^{p(j1)p(j2) + p(k1)p(k2)} O[j2,j1,k2,k1]
//
// A bond Hamiltonian is physically symmetric between its two ends, so it is a
// fixed point of that map and **the same matrix serves both descriptions** —
// which is why both forms above are fed the identical gate below. The test
// asserts that fixed-point property first, so it cannot silently rest on a
// false premise.
//
// The prefactor is -1 exactly on the (odd,odd) <-> (even,even) channels. For a
// particle-number-conserving spinless model those are empty — (odd,odd) has
// N = 2 and (even,even) has N = 0 — so at d = 2 the graded and the plain
// transpose coincide and either code path is fine; that case is kept as the
// control. For electrons they are the doublon-holon channel,
// (|dn>,|up>) -> (|up dn>,|0>), present at first order in the hopping, so a
// plain transpose corrupts every bond given in the leg-1 form. Since
// tenes_std emits dy = +1 bonds exactly that way, in a 2D run that is every
// vertical bond, i.e. half of them.

namespace {

//! Check that a bond gate is invariant under the graded exchange of its two
//! sites, which is what lets the same matrix describe the bond from either
//! end. This is the premise of the test below, not an incidental check.
static void bo_require_graded_symmetry(
    const tenes::real_tensor &gate,
    const tenes::fermion::parity_vector &parity) {
  const int d = static_cast<int>(parity.size());
  double worst = 0.0;
  for (int j1 = 0; j1 < d; ++j1) {
    for (int j2 = 0; j2 < d; ++j2) {
      for (int k1 = 0; k1 < d; ++k1) {
        for (int k2 = 0; k2 < d; ++k2) {
          double direct = 0.0;
          double swapped = 0.0;
          gate.get_value(mptensor::Index(j1, j2, k1, k2), direct);
          gate.get_value(mptensor::Index(j2, j1, k2, k1), swapped);
          const bool flip =
              (parity[j1] && parity[j2]) != (parity[k1] && parity[k2]);
          worst =
              std::max(worst, std::abs(direct - (flip ? -swapped : swapped)));
        }
      }
    }
  }
  INFO("largest deviation from graded exchange symmetry: " << worst);
  REQUIRE(worst < 1.0e-12);
}

/*!
 * Run @p nsteps of simple update on the vertical bonds of a 2x2 cell and
 * return the bond spectra, keyed by the lower site of each bond.
 *
 * @param[in] from_lower When true each bond is named by its lower site
 *        (source_leg = 1); when false by its upper site (source_leg = 3).
 *        Both name the same four bonds in the same order.
 */
static std::vector<std::vector<double>> bo_bond_spectra(
    int phys_dim, const tenes::fermion::parity_vector &parity,
    const tenes::real_tensor &gate, bool from_lower, int nsteps,
    const std::string &outdir) {
  using tensor = tenes::real_tensor;

  tenes::SquareLattice lattice(2, 2);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    lattice.physical_dims[site] = phys_dim;
    lattice.virtual_dims[site] = {1, 2, 1, 2};
    lattice.initial_dirs[site] = {0.0};
    lattice.noises[site] = 1.0;
  }

  tenes::itps::PEPS_Parameters params;
  params.fermion = true;
  params.phys_parity.assign(lattice.N_UNIT, parity);
  params.print_level = tenes::PrintLevel::none;
  params.outdir = outdir;
  params.CHI = 8;
  params.Max_CTM_Iteration = 10;
  params.CTM_Convergence_Epsilon = 1.0e-8;
  params.Use_RSVD = false;
  params.seed = 11;

  std::vector<tenes::EvolutionOperator<tensor>> updates;
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    // The bond (site, top(site)), named from one end or the other.
    const int source = from_lower ? site : lattice.top(site);
    const int leg = from_lower ? 1 : 3;
    updates.push_back(
        tenes::make_twosite_EvolutionOperator(source, leg, 0, gate));
  }

  tenes::itps::iTPS<tensor> state(
      MPI_COMM_WORLD, params, lattice, updates,
      tenes::EvolutionOperators<tensor>{}, tenes::Operators<tensor>{},
      tenes::Operators<tensor>{}, tenes::Operators<tensor>{},
      tenes::itps::CorrelationParameter{},
      tenes::itps::TransferMatrix_Parameters{});

  for (int step = 0; step < nsteps; ++step) {
    for (const auto &update : updates) {
      state.simple_update(update);
    }
  }

  auto &lambda = tenes::itps::iTPSTestAccessor::lambda_tensor(state);
  std::vector<std::vector<double>> spectra(lattice.N_UNIT);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    spectra[site] = sorted_desc(lambda[site][1]);
  }
  return spectra;
}

static void bo_check_both_ends(int phys_dim,
                               const tenes::fermion::parity_vector &parity,
                               const tenes::real_tensor &gate, int nsteps,
                               const std::string &tag) {
  bo_require_graded_symmetry(gate, parity);
  const auto lower =
      bo_bond_spectra(phys_dim, parity, gate, true, nsteps,
                      "output_test_bond_orientation_" + tag + "_lower");
  const auto upper =
      bo_bond_spectra(phys_dim, parity, gate, false, nsteps,
                      "output_test_bond_orientation_" + tag + "_upper");
  for (std::size_t bond = 0; bond < lower.size(); ++bond) {
    INFO(tag << " bond " << bond
             << " named from the lower site=" << vector_to_string(lower[bond])
             << " from the upper site=" << vector_to_string(upper[bond]));
    CHECK(lambda_relative_diff(lower[bond], upper[bond]) < 1.0e-10);
  }
}

}  // namespace

TEST_CASE("simple update does not depend on which end names a bond") {
  const int nsteps = 20;

  SUBCASE("spinless, d = 2 (control: the graded swap sign cannot show here)") {
    bo_check_both_ends(2, tenes::fermion::parity_vector{false, true},
                       make_free_fermion_gate(0.01), nsteps, "d2");
  }

  SUBCASE("electron, d = 4 (the doublon-holon channel is populated)") {
    bo_check_both_ends(4,
                       tenes::fermion::parity_vector{false, true, true, false},
                       electron_gate(1.0, 0.0, 0.0, 0.01), nsteps, "d4");
  }
}
