/* TeNeS - Massively parallel tensor network solver /
/ Copyright (C) 2019- The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify /
/ it under the terms of the GNU General Public License as published by /
/ the Free Software Foundation, either version 3 of the License, or /
/ (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, /
/ but WITHOUT ANY WARRANTY; without even the implied warranty of /
/ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the /
/ GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License /
/ along with this program. If not, see http://www.gnu.org/licenses/. */

// ===== T6: the fermion-mode guards that need a whole run =====================
//
// Contract: docs/superpowers/specs/2026-09-04-fermion-full-update-contract.md
// section 3.5 (T6) with the 2026-09-05 revision on configuration and
// placement. The two cases that drive the solver end to end (the fast
// full-update fallback, and the ledger consistency after a full update that
// gates bonds from their raster-later end) used to live in test/input.cpp
// with a deliberately short schedule (20 simple-update steps, 20 CTM sweeps).
// That schedule leaves the fold CTM unconverged, its parity leak trips the
// forbidden-block check of build_full_update_environment (1e-8), and the
// cases fail for a reason that has nothing to do with the guards. They are
// now run on a schedule the CTM converges on - 200 simple-update steps,
// iteration_max 100, convergence_epsilon 1e-8, D = 2, chi = 8 - and assert
// as a PREMISE that the run never printed "CTM did not converge". The
// requirements themselves are unchanged: the fallback announces itself with
// a warning naming Full_Use_FastFullUpdate and completes, and the virtual
// parity ledgers agree at both ends of every bond afterwards.
//
// The solver prints the CTM warning on standard output and the fallback
// warning on standard error, both only at print level warn or above; the
// runs are therefore made at PrintLevel::warn with both streams captured.
// The pure load checks of T6 ("accepted", "refused") stay in test/input.cpp.
//
// Included into the test_fermion_layer TU; it uses that TU's
// iTPSTestAccessor (Tn, finfo).

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <tuple>

#include "../../src/iTPS/load_toml.hpp"
#include "../../src/iTPS/main.hpp"
#include "../../src/iTPS/transfer_matrix.hpp"

namespace {

//! Redirect one standard stream into a string for the lifetime of the
//! object. A doctest cannot read the process' own stdout/stderr, so the
//! stream buffer is swapped instead.
class fgd_stream_capture {
 public:
  explicit fgd_stream_capture(std::ostream& stream) : stream_(stream) {
    saved_ = stream_.rdbuf(buffer_.rdbuf());
  }
  ~fgd_stream_capture() { stream_.rdbuf(saved_); }
  fgd_stream_capture(fgd_stream_capture const&) = delete;
  fgd_stream_capture& operator=(fgd_stream_capture const&) = delete;
  std::string str() const { return buffer_.str(); }

 private:
  std::ostream& stream_;
  std::ostringstream buffer_;
  std::streambuf* saved_ = nullptr;
};

//! Contract 3.5 (T6, 2026-09-05 revision): the schedule the CTM converges
//! on. Reused verbatim by both cases below.
constexpr int fgd_simple_steps = 200;
constexpr int fgd_ctm_iteration_max = 100;
constexpr double fgd_ctm_epsilon = 1.0e-8;
constexpr int fgd_D = 2;
constexpr int fgd_chi = 8;
constexpr const char* fgd_ctm_warning = "CTM did not converge";

//! (source_site, source_leg) of eight gates that cover the eight bonds of a
//! 2x2 unit cell exactly once, using every source_leg value.
//!
//! src/SquareLattice.cpp maps leg 0 to (x-1, y), leg 1 to (x, y+1), leg 2 to
//! (x+1, y) and leg 3 to (x, y-1), with site = x + 2 * y. The raster-earlier
//! site of a bond is its left (horizontal) or upper (vertical) end, so legs 2
//! and 3 name a bond from the earlier end and legs 0 and 1 from the later
//! one; the latter are the cases that force the driver to swap the two sites.
constexpr int fgd_full_update_bonds[8][2] = {{0, 2}, {2, 2}, {0, 0}, {2, 0},
                                             {0, 1}, {2, 1}, {1, 3}, {3, 3}};

std::string fgd_cell_toml() {
  return R"(
[tensor]
L_sub = [2, 2]
skew = 0
[[tensor.unitcell]]
index = []
physical_dim = 2
virtual_dim = 2
parity = [0, 1]
noise = 0.01
)";
}

}  // namespace

TEST_CASE(
    "fermion guards T6: fermion mode falls back from the fast full update "
    "with a warning on a converged CTM") {
  using namespace tenes;
  using namespace tenes::itps;

  // Precondition: fastfullupdate really is on by default, so the input
  // below exercises the fallback rather than the plain path.
  {
    auto defaults =
        gen_param(toml::parse_str(R"([parameter])").at("parameter"));
    REQUIRE(defaults.Full_Use_FastFullUpdate == true);
  }

  const std::string input_filename = "fermion_guards_fast_full_update.toml";
  const std::string outdir = "output_fermion_guards_fast_full_update";

  // A parity-even hopping gate; the diagonal keeps the state from collapsing
  // and the off-diagonal block makes the bond genuinely correlated.
  const std::string gate_elements =
      "0 0 0 0  1.0 0.0\n"
      "0 1 0 1  1.0 0.0\n"
      "1 0 1 0  1.0 0.0\n"
      "1 1 1 1  1.0 0.0\n"
      "0 1 1 0  0.1 0.0\n"
      "1 0 0 1  0.1 0.0";

  {
    std::ostringstream ofs;
    ofs << R"(
[parameter]
[parameter.general]
is_real = true
fermion = true
output = ")"
        << outdir << R"("

[parameter.simple_update]
tau = 0.01
num_step = )"
        << fgd_simple_steps << R"(

[parameter.full_update]
tau = 0.01
num_step = 1

[parameter.ctm]
dimension = )"
        << fgd_chi << R"(
convergence_epsilon = )"
        << fgd_ctm_epsilon << R"(
iteration_max = )"
        << fgd_ctm_iteration_max << R"(

[parameter.random]
seed = 11
)" << fgd_cell_toml()
        << R"(
[observable]
[[observable.onesite]]
name = "n"
group = 0
sites = []
dim = 2
elements = """
1 1  1.0 0.0
"""

[[observable.twosite]]
name = "bond_hamiltonian"
group = 0
bonds = """
0 1 0
1 1 0
2 1 0
3 1 0
0 0 1
1 0 1
2 0 1
3 0 1
"""
dim = [2, 2]
elements = """
0 1 1 0  -1.0 0.0
1 0 0 1  -1.0 0.0
"""

[evolution]
)";
    for (int leg : {2, 1}) {
      for (int site = 0; site < 4; ++site) {
        for (const char* kind : {"simple", "full"}) {
          ofs << "[[evolution." << kind << "]]\n"
              << "group = 0\n"
              << "source_site = " << site << "\n"
              << "source_leg = " << leg << "\n"
              << "dimensions = [2, 2, 2, 2]\n"
              << "elements = \"\"\"\n"
              << gate_elements << "\n\"\"\"\n";
        }
      }
    }
    std::ofstream out(input_filename.c_str());
    out << ofs.str();
  }
  // Premise: the written input really carries the converged-CTM schedule
  // (the toml is generated above; read it back through the loader).
  {
    const auto parsed = toml::parse(input_filename);
    const PEPS_Parameters p = gen_param(parsed.at("parameter"));
    REQUIRE(p.fermion == true);
    REQUIRE(p.num_simple_step.size() == 1);
    REQUIRE(p.num_simple_step[0] >= 200);
    REQUIRE(p.Max_CTM_Iteration >= 100);
    REQUIRE(p.CTM_Convergence_Epsilon == fgd_ctm_epsilon);
    REQUIRE(p.CHI == fgd_chi);
    REQUIRE(p.num_full_step.size() == 1);
    REQUIRE(p.num_full_step[0] == 1);
    REQUIRE(p.Full_Use_FastFullUpdate == true);
  }

  std::string captured_out;
  std::string captured_err;
  {
    fgd_stream_capture err(std::cerr);
    fgd_stream_capture out(std::cout);
    CHECK_NOTHROW(tenes::itps::itps_main(input_filename, MPI_COMM_WORLD,
                                         PrintLevel::warn));
    captured_out = out.str();
    captured_err = err.str();
  }
  INFO("captured stdout: " << captured_out);
  INFO("captured stderr: " << captured_err);
  // Premise: every CTM of the run converged, so the full update ran on a
  // parity-clean environment and any failure below is about the guard.
  REQUIRE(captured_out.find(fgd_ctm_warning) == std::string::npos);
  REQUIRE(captured_err.find(fgd_ctm_warning) == std::string::npos);
  // T6: the fallback announces itself on standard error and the run
  // completes (the CHECK_NOTHROW above).
  CHECK(captured_err.find("Full_Use_FastFullUpdate") != std::string::npos);
  // ... and it really ran the full update: the run wrote its density file.
  CHECK(std::filesystem::exists(std::filesystem::path(outdir) / "density.dat"));

  std::remove(input_filename.c_str());
  std::error_code ec;
  std::filesystem::remove_all(outdir, ec);
}

TEST_CASE(
    "fermion guards T6: a fermionic full update leaves the parity ledger "
    "consistent on a converged CTM") {
  using namespace tenes;
  using namespace tenes::itps;
  using ptensor = complex_tensor;
  MPI_Comm comm = MPI_COMM_WORLD;

  auto tensor_toml = toml::parse_str(fgd_cell_toml());
  SquareLattice lattice = gen_lattice(tensor_toml.at("tensor"));

  std::ostringstream param_text;
  param_text << R"(
[parameter]
[parameter.general]
fermion = true
[parameter.simple_update]
tau = 0.01
num_step = )" << fgd_simple_steps
             << R"(
[parameter.full_update]
tau = 0.01
num_step = 1
fastfullupdate = false
[parameter.ctm]
dimension = )"
             << fgd_chi << R"(
convergence_epsilon = )"
             << fgd_ctm_epsilon << R"(
iteration_max = )"
             << fgd_ctm_iteration_max << R"(
[parameter.random]
seed = 11
)";
  auto param_toml = toml::parse_str(param_text.str());
  PEPS_Parameters peps_parameters = gen_param(param_toml.at("parameter"));
  peps_parameters.phys_parity =
      gen_phys_parity(tensor_toml.at("tensor"), lattice);
  // warn: the CTM non-convergence warning (the premise below) is printed at
  // this level and not below it.
  peps_parameters.print_level = PrintLevel::warn;
  peps_parameters.outdir = "output_fermion_guards_full_ledger";
  REQUIRE(peps_parameters.num_simple_step[0] >= 200);
  REQUIRE(peps_parameters.Max_CTM_Iteration >= 100);
  REQUIRE(peps_parameters.CTM_Convergence_Epsilon == fgd_ctm_epsilon);
  REQUIRE(peps_parameters.CHI == fgd_chi);
  REQUIRE(lattice.virtual_dims[0][0] == fgd_D);

  auto gate = [&]() {
    ptensor op(comm, mptensor::Shape(2, 2, 2, 2));
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        op.set_value(mptensor::Index(i, j, i, j), 1.0);
      }
    }
    op.set_value(mptensor::Index(0, 1, 1, 0), 0.1);
    op.set_value(mptensor::Index(1, 0, 0, 1), 0.1);
    return op;
  }();

  EvolutionOperators<ptensor> simple_updates;
  for (int leg : {2, 1}) {
    for (int site = 0; site < 4; ++site) {
      simple_updates.push_back(
          make_twosite_EvolutionOperator<ptensor>(site, leg, 0, gate));
    }
  }
  EvolutionOperators<ptensor> full_updates;
  for (const auto& bond : fgd_full_update_bonds) {
    full_updates.push_back(
        make_twosite_EvolutionOperator<ptensor>(bond[0], bond[1], 0, gate));
  }

  // Preconditions. The bonds must be gated from the raster-later end at
  // least once (source_leg 0 and 1), because that is the path whose ledger
  // bookkeeping this test is about; and the eight gates must name eight
  // different bonds.
  auto bond_name = [&](int site, int leg) {
    const int x = site % lattice.LX;
    const int y = site / lattice.LX;
    if (leg == 2) return std::make_tuple(0, x, y);
    if (leg == 0) {
      return std::make_tuple(0, (x - 1 + lattice.LX) % lattice.LX, y);
    }
    if (leg == 1) return std::make_tuple(1, x, y);
    return std::make_tuple(1, x, (y - 1 + lattice.LY) % lattice.LY);
  };
  std::set<std::tuple<int, int, int>> named;
  bool has_leg0 = false;
  bool has_leg1 = false;
  for (const auto& up : full_updates) {
    named.insert(bond_name(up.source_site, up.source_leg));
    has_leg0 = has_leg0 || up.source_leg == 0;
    has_leg1 = has_leg1 || up.source_leg == 1;
  }
  REQUIRE(has_leg0);
  REQUIRE(has_leg1);
  REQUIRE(named.size() == 8);
  REQUIRE(peps_parameters.num_full_step[0] == 1);
  REQUIRE(peps_parameters.Full_Use_FastFullUpdate == false);

  iTPS<ptensor> state(comm, peps_parameters, lattice, simple_updates,
                      full_updates, Operators<ptensor>{}, Operators<ptensor>{},
                      Operators<ptensor>{}, CorrelationParameter{},
                      TransferMatrix_Parameters{});
  REQUIRE(iTPSTestAccessor::finfo(state).enabled);

  std::string captured_out;
  std::string captured_err;
  {
    fgd_stream_capture err(std::cerr);
    fgd_stream_capture out(std::cout);
    CHECK_NOTHROW(state.optimize());
    captured_out = out.str();
    captured_err = err.str();
  }
  INFO("captured stdout: " << captured_out);
  INFO("captured stderr: " << captured_err);
  // Premise: every CTM of the run converged (see the file comment).
  REQUIRE(captured_out.find(fgd_ctm_warning) == std::string::npos);
  REQUIRE(captured_err.find(fgd_ctm_warning) == std::string::npos);

  // The full update must actually have moved the state: a no-op would keep
  // every ledger trivially consistent and make the checks below hollow. The
  // baseline runs the identical simple update and stops there.
  {
    PEPS_Parameters su_only = peps_parameters;
    su_only.num_full_step.assign(su_only.num_full_step.size(), 0);
    su_only.print_level = PrintLevel::none;
    iTPS<ptensor> baseline(comm, su_only, lattice, simple_updates, full_updates,
                           Operators<ptensor>{}, Operators<ptensor>{},
                           Operators<ptensor>{}, CorrelationParameter{},
                           TransferMatrix_Parameters{});
    baseline.optimize();
    const auto& after = iTPSTestAccessor::Tn(state);
    const auto& before = iTPSTestAccessor::Tn(baseline);
    REQUIRE(after.size() == before.size());
    double moved = 0.0;
    for (std::size_t site = 0; site < after.size(); ++site) {
      moved = std::max(moved, mptensor::max_abs(after[site] - before[site]));
    }
    INFO("largest change of Tn caused by the full update: " << moved);
    REQUIRE(moved > 1e-8);
  }

  const auto& fi = iTPSTestAccessor::finfo(state);
  CHECK_NOTHROW(tenes::fermion::validate_neighbor_consistency(fi, lattice));

  // The per-bond contract (T3-iii) asks for 1e-12 * max_abs on a single
  // bond update; this runs two hundred simple-update steps and a full-update
  // sweep before looking, so the bound is relaxed by two decades. It is
  // still far tighter than anything a leaked odd component would satisfy.
  const auto& tensors = iTPSTestAccessor::Tn(state);
  REQUIRE(tensors.size() == static_cast<std::size_t>(lattice.N_UNIT));
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    auto ft = tenes::fermion::wrap_Tn(tensors[site], fi, site);
    const double scale = tenes::fermion::max_abs(ft);
    REQUIRE(scale > 0.0);
    CHECK(tenes::fermion::parity_violation(ft) <= 1e-10 * scale);
  }

  std::error_code ec;
  std::filesystem::remove_all(peps_parameters.outdir, ec);
}
