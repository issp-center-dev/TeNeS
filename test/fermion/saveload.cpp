// ===== SL: save_tensors / load_tensors with the fermionic parity ledger ====
//
// Pins iTPS<tensor>::save_tensors() / load_tensors() under
// peps_parameters.fermion = true against the save/load contract:
//   layer 1: the virtual parity ledger (finfo.virt) survives the round trip
//            element by element. The ledger used here is deliberately NOT
//            even_first_parity(4) -- with the default split a loader that
//            never reads fermion.dat would pass, so the test would be hollow.
//   layer 2: the mean-field measured values (measure_onesite /
//            measure_twosite, MeanField_Env = true) and lambda_tensor agree
//            between the saving and the loading state.
//   layer 3: every guard V1..V8 of the contract throws, and the front
//            validation (V5) leaves the already-built tensors untouched.
//   layer 3b: validate_fermion_constraints accepts tensor_save / tensor_load.
// References are taken from the values held before saving, from the input
// parameters and from exceptions only -- never from the (not yet existing)
// serialisation helpers.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../../src/exception.hpp"
#include "../../src/util/file.hpp"
#include "../../src/util/string.hpp"
#include "../../src/iTPS/load_toml.hpp"

namespace fermion_saveload {

using tenes::itps::Bond;
using tenes::itps::iTPSTestAccessor;
using lambda_table = std::vector<std::vector<std::vector<double>>>;

// 2x2 cell, D = 4 on every virtual leg, d = 2.
constexpr int sl_LX = 2;
constexpr int sl_LY = 2;
constexpr int sl_N_UNIT = sl_LX * sl_LY;
constexpr int sl_D = 4;
constexpr int sl_pdim = 2;

const tenes::fermion::parity_vector sl_phys{false, true};
// even_first_parity(4) is {0, 0, 1, 1} (two even, two odd); this ledger is
// three even and one odd, so nothing but fermion.dat can produce it.
const tenes::fermion::parity_vector sl_virt{false, false, false, true};
// A second, equally valid ledger (odd first, and again not the default one)
// used to break the agreement between fermion.dat and the stored Tn (V7).
const tenes::fermion::parity_vector sl_virt_other{true, false, false, false};

inline tenes::SquareLattice make_lattice(int vdim = sl_D) {
  tenes::SquareLattice lattice(sl_LX, sl_LY);
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    lattice.physical_dims[site] = sl_pdim;
    lattice.virtual_dims[site] = {vdim, vdim, vdim, vdim};
    lattice.initial_dirs[site] = {0.0};
    lattice.noises[site] = 1.0;
  }
  return lattice;
}

// One parent directory for every case: running this binary outside the build
// tree then leaves a single directory behind instead of one per case.
inline std::string case_dir_name(const std::string &tag) {
  return "output_test_fermion_saveload/" + tag;
}

inline tenes::itps::PEPS_Parameters make_params(bool fermion,
                                                const std::string &tag,
                                                const std::string &save_dir,
                                                const std::string &load_dir) {
  tenes::itps::PEPS_Parameters params;
  params.fermion = fermion;
  if (fermion) {
    params.phys_parity.assign(sl_N_UNIT, sl_phys);
  }
  params.MeanField_Env = true;
  params.print_level = tenes::PrintLevel::none;
  params.outdir = case_dir_name(tag);
  params.CHI = 4;
  params.Use_RSVD = false;
  params.seed = 11;
  params.tensor_save_dir = save_dir;
  params.tensor_load_dir = load_dir;
  return params;
}

inline std::string save_dir_name(const std::string &tag) {
  return case_dir_name(tag) + "/tensors";
}

// The tests create their own directories; a stale one from an earlier run
// would let a case pass on files it did not write.
inline void reset_dir(const std::string &dir) {
  const std::string cmd = "rm -rf '" + dir + "'";
  std::system(cmd.c_str());
}

inline tenes::fermion::leg_parities site_parities(
    const tenes::fermion::parity_vector &virt,
    const tenes::fermion::parity_vector &phys) {
  return {virt, virt, virt, virt, phys};
}

// Only the parity-even elements of the ledger are non-zero.
inline tenes::real_tensor deterministic_tn(
    int site, const tenes::fermion::leg_parities &parity) {
  tenes::real_tensor t(mptensor::Shape(parity[0].size(), parity[1].size(),
                                       parity[2].size(), parity[3].size(),
                                       parity[4].size()));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const auto idx = t.global_index(n);
    if (tenes::fermion::count_odd(parity, idx) % 2 != 0) {
      continue;
    }
    double x = 1.0;
    for (int ax = 0; ax < 5; ++ax) {
      x += static_cast<double>((ax + 3) * idx[ax]);
    }
    x *= site + 2;
    t.set_value(idx, 0.19 * std::sin(x) + 0.13 * std::cos(0.37 * x));
  }
  return t;
}

inline std::vector<tenes::real_tensor> deterministic_state(
    const tenes::fermion::parity_vector &virt) {
  std::vector<tenes::real_tensor> Tn;
  for (int site = 0; site < sl_N_UNIT; ++site) {
    Tn.push_back(deterministic_tn(site, site_parities(virt, sl_phys)));
  }
  return Tn;
}

// A tensor with no parity structure at all, used as a marker: if a failed
// load overwrote Tn, these values would be gone.
inline tenes::real_tensor marker_tn(int site) {
  tenes::real_tensor t(mptensor::Shape(sl_D, sl_D, sl_D, sl_D, sl_pdim));
  for (std::size_t n = 0; n < t.local_size(); ++n) {
    const auto idx = t.global_index(n);
    double v = 100.0 + site;
    for (int ax = 0; ax < 5; ++ax) {
      v += static_cast<double>(idx[ax]) / (ax + 2.0);
    }
    t.set_value(idx, v);
  }
  return t;
}

// Every leg of every site carries its own lambda; the two ends of a bond
// share it, as the simple update leaves them
// (lambda[s][2] == lambda[right(s)][0], lambda[s][1] == lambda[up(s)][3]).
// None of these values survives a six-significant-digit text round trip, so
// the exact comparison of lambda_tensor holds only if the values are written
// with full double precision (see survives_default_precision below).
inline lambda_table bond_lambdas(const tenes::SquareLattice &lattice) {
  const std::vector<std::vector<double>> horizontal{
      {1.0 / 3.0, std::sqrt(2.0) / 2.0, 0.1234567890123, 2.0 / 7.0},
      {1.0 / 7.0, std::sqrt(3.0) / 3.0, 0.9876543210987, 3.0 / 11.0},
      {2.0 / 3.0, std::sqrt(5.0) / 4.0, 0.3141592653589, 5.0 / 13.0},
      {4.0 / 9.0, std::sqrt(7.0) / 5.0, 0.2718281828459, 6.0 / 17.0}};
  const std::vector<std::vector<double>> vertical{
      {5.0 / 6.0, std::sqrt(11.0) / 6.0, 0.4142135623730, 7.0 / 19.0},
      {8.0 / 9.0, std::sqrt(13.0) / 7.0, 0.5772156649015, 8.0 / 23.0},
      {3.0 / 7.0, std::sqrt(17.0) / 8.0, 0.6180339887498, 9.0 / 29.0},
      {5.0 / 11.0, std::sqrt(19.0) / 9.0, 0.7320508075688, 10.0 / 31.0}};
  lambda_table lambda(lattice.N_UNIT, std::vector<std::vector<double>>(4));
  for (int site = 0; site < lattice.N_UNIT; ++site) {
    lambda[site][2] = horizontal[site];
    lambda[lattice.other(site, 1, 0)][0] = horizontal[site];
    lambda[site][1] = vertical[site];
    lambda[lattice.other(site, 0, 1)][3] = vertical[site];
  }
  return lambda;
}

// Real symmetric one-site operator; the diagonal keeps the expectation
// value away from zero even though the off-diagonal part is parity-odd.
inline tenes::real_tensor onesite_op() {
  tenes::real_tensor op(mptensor::Shape(sl_pdim, sl_pdim));
  op.set_value(mptensor::Index(0, 0), 0.2);
  op.set_value(mptensor::Index(0, 1), 0.5);
  op.set_value(mptensor::Index(1, 0), 0.5);
  op.set_value(mptensor::Index(1, 1), -0.3);
  return op;
}

// op[in1, in2, out1, out2]: -(c1^dag c2 + c2^dag c1).
inline tenes::real_tensor hopping_op() {
  tenes::real_tensor op(mptensor::Shape(2, 2, 2, 2));
  op.set_value(mptensor::Index(0, 1, 1, 0), -1.0);
  op.set_value(mptensor::Index(1, 0, 0, 1), -1.0);
  return op;
}

inline tenes::Operators<tenes::real_tensor> make_onesite_operators() {
  tenes::Operators<tenes::real_tensor> ops;
  for (int site = 0; site < sl_N_UNIT; ++site) {
    ops.emplace_back("onesite", 0, site, onesite_op());
  }
  return ops;
}

inline tenes::Operators<tenes::real_tensor> make_twosite_operators() {
  tenes::Operators<tenes::real_tensor> ops;
  for (int site = 0; site < sl_N_UNIT; ++site) {
    ops.emplace_back("twosite", 0, site, 1, 0, hopping_op());
    ops.emplace_back("twosite", 0, site, 0, 1, hopping_op());
  }
  return ops;
}

using state_type = tenes::itps::iTPS<tenes::real_tensor>;

inline state_type make_state(const tenes::SquareLattice &lattice,
                             const tenes::itps::PEPS_Parameters &params,
                             bool with_operators) {
  return state_type(MPI_COMM_WORLD, params, lattice,
                    tenes::EvolutionOperators<tenes::real_tensor>{},
                    tenes::EvolutionOperators<tenes::real_tensor>{},
                    with_operators ? make_onesite_operators()
                                   : tenes::Operators<tenes::real_tensor>{},
                    with_operators ? make_twosite_operators()
                                   : tenes::Operators<tenes::real_tensor>{},
                    tenes::Operators<tenes::real_tensor>{},
                    tenes::itps::CorrelationParameter{},
                    tenes::itps::TransferMatrix_Parameters{});
}

// The message of the tenes::load_error a fermionic load throws, or a marker
// string when it throws something else (or nothing). Used to check that a
// guard did not merely stand in for the one under test.
inline std::string load_error_message(const tenes::SquareLattice &lattice,
                                      const std::string &tag,
                                      const std::string &dir) {
  try {
    make_state(lattice, make_params(true, tag, "", dir), false);
  } catch (const tenes::load_error &e) {
    return e.what();
  } catch (const std::exception &e) {
    return std::string("<not a load_error: ") + e.what() + ">";
  }
  return "<no exception>";
}

inline void inject(state_type &state, const std::vector<tenes::real_tensor> &Tn,
                   const lambda_table &lambda,
                   const tenes::fermion::parity_vector &virt) {
  auto &state_Tn = iTPSTestAccessor::Tn(state);
  REQUIRE(state_Tn.size() == Tn.size());
  for (std::size_t site = 0; site < Tn.size(); ++site) {
    state_Tn[site] = Tn[site];
  }
  iTPSTestAccessor::lambda_tensor(state) = lambda;
  auto &finfo = iTPSTestAccessor::finfo(state);
  REQUIRE(finfo.enabled);
  REQUIRE(finfo.virt.size() == Tn.size());
  for (auto &site_virt : finfo.virt) {
    for (auto &leg_parity : site_virt) {
      leg_parity = virt;
    }
  }
}

inline double max_abs_diff(const tenes::real_tensor &a,
                           const tenes::real_tensor &b) {
  double diff = 0.0;
  for (std::size_t n = 0; n < a.local_size(); ++n) {
    const auto idx = a.global_index(n);
    double va = 0.0;
    double vb = 0.0;
    a.get_value(idx, va);
    b.get_value(idx, vb);
    diff = std::max(diff, std::abs(va - vb));
  }
  return diff;
}

// True if `v` comes back unchanged from a text round trip at the default
// stream precision (six significant digits). The lambda fixture is built
// from values for which this is false, so lambda_tensor can only match
// exactly if save_tensors() writes full double precision.
inline bool survives_default_precision(double v) {
  std::ostringstream oss;
  oss << v;
  return std::stod(oss.str()) == v;
}

inline bool any_lambda_survives_default_precision(const lambda_table &lambda) {
  for (const auto &site : lambda) {
    for (const auto &leg : site) {
      for (const double v : leg) {
        if (survives_default_precision(v)) {
          return true;
        }
      }
    }
  }
  return false;
}

inline std::vector<std::string> read_lines(const std::string &path) {
  std::ifstream ifs(path.c_str());
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(ifs, line)) {
    lines.push_back(line);
  }
  while (!lines.empty() && tenes::util::strip(lines.back()).empty()) {
    lines.pop_back();
  }
  return lines;
}

inline void write_lines(const std::string &path,
                        const std::vector<std::string> &lines) {
  std::ofstream ofs(path.c_str());
  for (const auto &line : lines) {
    ofs << line << "\n";
  }
}

inline void patch_line(const std::string &path, std::size_t index,
                       const std::string &text) {
  auto lines = read_lines(path);
  REQUIRE(lines.size() > index);
  lines[index] = text;
  write_lines(path, lines);
}

inline std::vector<int> parse_ints(const std::string &line) {
  std::vector<int> values;
  for (const auto &token :
       tenes::util::split(tenes::util::drop_comment(line))) {
    values.push_back(std::stoi(token));
  }
  return values;
}

inline std::vector<int> to_ints(const tenes::fermion::parity_vector &p) {
  std::vector<int> values;
  for (std::size_t i = 0; i < p.size(); ++i) {
    values.push_back(p[i] ? 1 : 0);
  }
  return values;
}

// Line layout of fermion.dat: 0 version, 1 N_UNIT, 2 L_sub, 3 skew, then
// five lines per site (physical leg, then virtual legs 0..3).
constexpr std::size_t sl_header_lines = 4;
inline std::size_t phys_line(int site) {
  return sl_header_lines + 5 * static_cast<std::size_t>(site);
}
inline std::size_t virt_line(int site, int leg) {
  return phys_line(site) + 1 + static_cast<std::size_t>(leg);
}

// Saves a fermionic state built from `virt_state` while declaring
// `virt_ledger` in finfo; the two coincide except in the V7 fixture.
inline std::string save_reference(
    const std::string &tag, const tenes::fermion::parity_vector &virt_state,
    const tenes::fermion::parity_vector &virt_ledger) {
  const std::string dir = save_dir_name(tag);
  reset_dir(case_dir_name(tag));
  const auto lattice = make_lattice();
  auto state = make_state(lattice, make_params(true, tag, dir, ""), false);
  inject(state, deterministic_state(virt_state), bond_lambdas(lattice),
         virt_ledger);
  state.save_tensors();
  return dir;
}

inline std::string save_reference(const std::string &tag) {
  return save_reference(tag, sl_virt, sl_virt);
}

}  // namespace fermion_saveload

TEST_CASE("SL save writes the fermionic parity ledger") {
  using namespace fermion_saveload;
  const std::string dir = save_reference("format");
  const std::string path = dir + "/fermion.dat";
  REQUIRE(tenes::util::path_exists(path));

  const auto lines = read_lines(path);
  REQUIRE(lines.size() == sl_header_lines + 5 * sl_N_UNIT);
  CHECK(parse_ints(lines[0]) == std::vector<int>{1});
  CHECK(parse_ints(lines[1]) == std::vector<int>{sl_N_UNIT});
  CHECK(parse_ints(lines[2]) == std::vector<int>{sl_LX, sl_LY});
  CHECK(parse_ints(lines[3]) == std::vector<int>{0});
  for (int site = 0; site < sl_N_UNIT; ++site) {
    INFO("site " << site);
    CHECK(parse_ints(lines[phys_line(site)]) == to_ints(sl_phys));
    for (int leg = 0; leg < 4; ++leg) {
      INFO("leg " << leg);
      CHECK(parse_ints(lines[virt_line(site, leg)]) == to_ints(sl_virt));
    }
  }
}

TEST_CASE("SL boson save writes no parity ledger") {
  using namespace fermion_saveload;
  const std::string tag = "boson_format";
  const std::string dir = save_dir_name(tag);
  reset_dir(case_dir_name(tag));
  const auto lattice = make_lattice();
  make_state(lattice, make_params(false, tag, dir, ""), false).save_tensors();
  REQUIRE(tenes::util::path_exists(dir + "/params.dat"));
  CHECK_FALSE(tenes::util::path_exists(dir + "/fermion.dat"));
}

TEST_CASE("SL layer1 the virtual parity ledger survives the round trip") {
  using namespace fermion_saveload;
  // Anti-hollowness: the saved ledger must differ from the default one that
  // initialize_tensors() installs, otherwise a loader that ignores
  // fermion.dat would pass this test.
  REQUIRE(sl_virt != tenes::fermion::even_first_parity(sl_D));

  const std::string dir = save_reference("ledger");
  const auto lattice = make_lattice();
  auto loaded =
      make_state(lattice, make_params(true, "ledger_load", "", dir), false);

  auto &finfo = iTPSTestAccessor::finfo(loaded);
  REQUIRE(finfo.enabled);
  REQUIRE(finfo.virt.size() == static_cast<std::size_t>(sl_N_UNIT));
  REQUIRE(finfo.phys.size() == static_cast<std::size_t>(sl_N_UNIT));
  for (int site = 0; site < sl_N_UNIT; ++site) {
    INFO("site " << site);
    CHECK(finfo.phys[site] == sl_phys);
    for (int leg = 0; leg < 4; ++leg) {
      INFO("leg " << leg);
      CHECK(finfo.virt[site][leg] == sl_virt);
    }
  }
}

TEST_CASE("SL layer2 measured values survive the round trip") {
  using namespace fermion_saveload;
  const std::string tag = "measure";
  const std::string dir = save_dir_name(tag);
  reset_dir(case_dir_name(tag));

  const auto lattice = make_lattice();
  const auto Tn = deterministic_state(sl_virt);
  const auto lambda = bond_lambdas(lattice);
  // Anti-hollowness for the lambda comparison at the end of this case: not a
  // single entry of the fixture survives a six-significant-digit text round
  // trip.
  REQUIRE_FALSE(any_lambda_survives_default_precision(lambda));

  auto saved = make_state(lattice, make_params(true, tag, dir, ""), true);
  inject(saved, Tn, lambda, sl_virt);
  saved.save_tensors();
  const auto saved_onesite = saved.measure_onesite();
  const auto saved_twosite = saved.measure_twosite();

  auto loaded =
      make_state(lattice, make_params(true, tag + "_load", "", dir), true);
  const auto loaded_onesite = loaded.measure_onesite();
  const auto loaded_twosite = loaded.measure_twosite();

  REQUIRE(saved_onesite.size() == loaded_onesite.size());
  for (std::size_t group = 0; group < saved_onesite.size(); ++group) {
    INFO("onesite group " << group);
    REQUIRE(saved_onesite[group].size() == loaded_onesite[group].size());
    for (std::size_t site = 0; site < saved_onesite[group].size(); ++site) {
      INFO("site " << site);
      CHECK(loaded_onesite[group][site] ==
            doctest::Approx(saved_onesite[group][site]).epsilon(1e-12));
    }
  }

  REQUIRE(saved_twosite.size() == loaded_twosite.size());
  for (std::size_t group = 0; group < saved_twosite.size(); ++group) {
    INFO("twosite group " << group);
    REQUIRE(saved_twosite[group].size() == loaded_twosite[group].size());
    for (const auto &entry : saved_twosite[group]) {
      const Bond &bond = entry.first;
      INFO("bond (" << bond.source_site << ", " << bond.dx << ", " << bond.dy
                    << ")");
      REQUIRE(loaded_twosite[group].count(bond) == 1);
      CHECK(loaded_twosite[group].at(bond) ==
            doctest::Approx(entry.second).epsilon(1e-12));
    }
  }

  // The saved state must not be trivial, otherwise the comparison above
  // would hold for any pair of states.
  double max_value = 0.0;
  for (const auto &group : saved_twosite) {
    for (const auto &entry : group) {
      max_value = std::max(max_value, std::abs(entry.second));
    }
  }
  CHECK(max_value > 1.0e-8);

  CHECK(iTPSTestAccessor::lambda_tensor(loaded) == lambda);
}

TEST_CASE("SL V1 a fermionic load without fermion.dat is an error") {
  using namespace fermion_saveload;
  const std::string tag = "v1";
  const std::string dir = save_dir_name(tag);
  reset_dir(case_dir_name(tag));
  const auto lattice = make_lattice();
  make_state(lattice, make_params(false, tag, dir, ""), false).save_tensors();
  REQUIRE_FALSE(tenes::util::path_exists(dir + "/fermion.dat"));

  CHECK_THROWS_AS(
      make_state(lattice, make_params(true, tag + "_load", "", dir), false),
      tenes::load_error);
}

TEST_CASE("SL V2 an unknown ledger format version is an error") {
  using namespace fermion_saveload;
  const std::string dir = save_reference("v2");
  const std::string path = dir + "/fermion.dat";
  REQUIRE(tenes::util::path_exists(path));
  patch_line(path, 0, "2 # Fermion_Format_Version");

  const auto lattice = make_lattice();
  CHECK_THROWS_AS(
      make_state(lattice, make_params(true, "v2_load", "", dir), false),
      tenes::load_error);
}

TEST_CASE("SL V3 a ledger with a different N_UNIT is an error") {
  using namespace fermion_saveload;
  const std::string dir = save_reference("v3");
  const std::string path = dir + "/fermion.dat";
  REQUIRE(tenes::util::path_exists(path));
  patch_line(path, 1, "5 # N_UNIT");

  const auto lattice = make_lattice();
  CHECK_THROWS_AS(
      make_state(lattice, make_params(true, "v3_load", "", dir), false),
      tenes::load_error);
}

TEST_CASE("SL V4 a ledger disagreeing with the input parity is an error") {
  using namespace fermion_saveload;
  const std::string dir = save_reference("v4");
  const std::string path = dir + "/fermion.dat";
  REQUIRE(tenes::util::path_exists(path));
  // input parity is [0, 1]
  patch_line(path, phys_line(0), "1 0 # parity of the physical leg of Tn[0]");

  const auto lattice = make_lattice();
  CHECK_THROWS_AS(
      make_state(lattice, make_params(true, "v4_load", "", dir), false),
      tenes::load_error);
}

TEST_CASE("SL V5 a ledger disagreeing with virtual_dim is an error") {
  using namespace fermion_saveload;
  const std::string dir = save_reference("v5");
  const std::string path = dir + "/fermion.dat";
  REQUIRE(tenes::util::path_exists(path));

  // The saved ledger has four entries per virtual leg; reading it back with
  // virtual_dim = 2 must be refused.
  const auto narrow_lattice = make_lattice(2);
  CHECK_THROWS_AS(
      make_state(narrow_lattice, make_params(true, "v5_load", "", dir), false),
      tenes::load_error);
}

TEST_CASE("SL V5 a refused load leaves the tensors untouched") {
  using namespace fermion_saveload;
  const std::string dir = save_reference("v5b");
  const std::string path = dir + "/fermion.dat";
  REQUIRE(tenes::util::path_exists(path));

  const auto lattice = make_lattice();
  auto loaded =
      make_state(lattice, make_params(true, "v5b_load", "", dir), false);

  // Overwrite the freshly loaded state with markers, then break the ledger
  // length so that the next load is refused. If the refusal came only after
  // the tensor files had been read, the markers would be gone.
  std::vector<tenes::real_tensor> markers;
  for (int site = 0; site < sl_N_UNIT; ++site) {
    markers.push_back(marker_tn(site));
  }
  const lambda_table marker_lambda(
      sl_N_UNIT,
      std::vector<std::vector<double>>(4, std::vector<double>(sl_D, 3.5)));
  auto &state_Tn = iTPSTestAccessor::Tn(loaded);
  REQUIRE(state_Tn.size() == markers.size());
  for (std::size_t site = 0; site < markers.size(); ++site) {
    state_Tn[site] = markers[site];
  }
  iTPSTestAccessor::lambda_tensor(loaded) = marker_lambda;

  for (int site = 0; site < sl_N_UNIT; ++site) {
    for (int leg = 0; leg < 4; ++leg) {
      patch_line(path, virt_line(site, leg),
                 "0 1 # parity of the virtual leg " + std::to_string(leg) +
                     " of Tn[" + std::to_string(site) + "]");
    }
  }

  CHECK_THROWS_AS(loaded.load_tensors(), tenes::load_error);

  for (int site = 0; site < sl_N_UNIT; ++site) {
    INFO("site " << site);
    const auto shape = state_Tn[site].shape();
    REQUIRE(shape.size() == 5);
    CHECK(shape[0] == static_cast<std::size_t>(sl_D));
    CHECK(shape[1] == static_cast<std::size_t>(sl_D));
    CHECK(shape[2] == static_cast<std::size_t>(sl_D));
    CHECK(shape[3] == static_cast<std::size_t>(sl_D));
    CHECK(shape[4] == static_cast<std::size_t>(sl_pdim));
    CHECK(max_abs_diff(state_Tn[site], markers[site]) == 0.0);
  }
  CHECK(iTPSTestAccessor::lambda_tensor(loaded) == marker_lambda);
}

TEST_CASE("SL V6a a ledger inconsistent between neighbors is an error") {
  using namespace fermion_saveload;
  const std::string dir = save_reference("v6a");
  const std::string path = dir + "/fermion.dat";
  REQUIRE(tenes::util::path_exists(path));
  // Leg 2 of site 0 faces leg 0 of site 1 on a 2x2 cell; changing only one
  // side breaks validate_neighbor_consistency without changing any length.
  patch_line(path, virt_line(0, 2),
             "1 0 0 0 # parity of the virtual leg 2 of Tn[0]");

  const auto lattice = make_lattice();
  // Not a bare CHECK_THROWS: this fixture also violates V7 (the tensors do
  // not respect the patched ledger), so "something threw" would stay green
  // even with the neighbor-consistency check deleted. The message pins which
  // guard fired -- the V7 message ("breaks fermion parity under the loaded
  // parity ledger") does not contain this text.
  CHECK_THROWS_WITH_AS(
      make_state(lattice, make_params(true, "v6a_load", "", dir), false),
      doctest::Contains("inconsistent virtual parity ledger"),
      tenes::load_error);
}

TEST_CASE("SL V6b a ledger with a different L_sub is an error") {
  using namespace fermion_saveload;
  const std::string dir = save_reference("v6b");
  const std::string path = dir + "/fermion.dat";
  REQUIRE(tenes::util::path_exists(path));
  patch_line(path, 2, "3 3 # L_sub");

  const auto lattice = make_lattice();
  CHECK_THROWS_AS(
      make_state(lattice, make_params(true, "v6b_load", "", dir), false),
      tenes::load_error);
}

TEST_CASE("SL V7 tensors violating the restored ledger are an error") {
  using namespace fermion_saveload;
  // Tn is built for sl_virt but the ledger written out is sl_virt_other,
  // so the restored ledger and the stored tensors disagree.
  REQUIRE(sl_virt != sl_virt_other);
  const auto fixture = deterministic_state(sl_virt);
  // The fixture is parity-clean under the ledger it was built with and
  // parity-violating under the one that gets written out.
  CHECK(tenes::fermion::parity_violation(
            tenes::fermion::ftensor<tenes::real_tensor>{
                fixture[0], site_parities(sl_virt, sl_phys)}) == 0.0);
  REQUIRE(tenes::fermion::parity_violation(
              tenes::fermion::ftensor<tenes::real_tensor>{
                  fixture[0], site_parities(sl_virt_other, sl_phys)}) > 1.0e-8);
  const std::string dir = save_reference("v7", sl_virt, sl_virt_other);

  const auto lattice = make_lattice();
  CHECK_THROWS_AS(
      make_state(lattice, make_params(true, "v7_load", "", dir), false),
      tenes::load_error);
}

TEST_CASE("SL V8 a bosonic load of a fermionic save is an error") {
  using namespace fermion_saveload;
  const std::string dir = save_reference("v8");
  REQUIRE(tenes::util::path_exists(dir + "/fermion.dat"));

  const auto lattice = make_lattice();
  CHECK_THROWS_AS(
      make_state(lattice, make_params(false, "v8_load", "", dir), false),
      tenes::load_error);
}

TEST_CASE("SL a malformed fermion.dat is a load error") {
  using namespace fermion_saveload;
  const auto lattice = make_lattice();

  SUBCASE("a parity entry that is not a number") {
    const std::string dir = save_reference("parse_value");
    const std::string path = dir + "/fermion.dat";
    REQUIRE(tenes::util::path_exists(path));
    // Four entries, so a length check (V5) cannot fire first; only the
    // conversion of "x" can fail.
    patch_line(path, virt_line(0, 0),
               "0 1 x 0 # parity of the virtual leg 0 of Tn[0]");
    CHECK_THROWS_WITH_AS(
        make_state(lattice, make_params(true, "parse_value_load", "", dir),
                   false),
        doctest::Contains("fermion.dat"), tenes::load_error);
    // "fermion.dat" on its own is a weak marker (the V7 hint names the file
    // too), so make sure the parity guard is not the one answering here.
    const std::string what =
        load_error_message(lattice, "parse_value_msg", dir);
    INFO("message: " << what);
    CHECK(what.find("breaks fermion parity") == std::string::npos);
  }

  SUBCASE("a header line with no value at all") {
    const std::string dir = save_reference("parse_empty");
    const std::string path = dir + "/fermion.dat";
    REQUIRE(tenes::util::path_exists(path));
    patch_line(path, 1, "# N_UNIT");
    CHECK_THROWS_WITH_AS(
        make_state(lattice, make_params(true, "parse_empty_load", "", dir),
                   false),
        doctest::Contains("fermion.dat"), tenes::load_error);
    const std::string what =
        load_error_message(lattice, "parse_empty_msg", dir);
    INFO("message: " << what);
    CHECK(what.find("breaks fermion parity") == std::string::npos);
  }
}

TEST_CASE("SL fermion accepts tensor save and load directories") {
  using namespace tenes;
  using namespace tenes::itps;
  using ptensor = tenes::real_tensor;

  auto param_toml = toml::parse_str(R"(
[parameter]
[parameter.general]
fermion = true
tensor_save = "save_dir"
tensor_load = "load_dir"
)");
  auto tensor_toml = toml::parse_str(R"(
[tensor]
L_sub = [2, 2]
[[tensor.unitcell]]
index = []
physical_dim = 2
virtual_dim = 2
parity = [0, 1]
)");
  PEPS_Parameters peps_parameters = gen_param(param_toml.at("parameter"));
  SquareLattice lattice = gen_lattice(tensor_toml.at("tensor"));
  peps_parameters.phys_parity =
      gen_phys_parity(tensor_toml.at("tensor"), lattice);
  REQUIRE(peps_parameters.tensor_save_dir == "save_dir");
  REQUIRE(peps_parameters.tensor_load_dir == "load_dir");
  CHECK_NOTHROW(validate_fermion_constraints(
      peps_parameters, lattice, EvolutionOperators<ptensor>{},
      EvolutionOperators<ptensor>{}, Operators<ptensor>{}, Operators<ptensor>{},
      Operators<ptensor>{}, CorrelationParameter{}));
}
