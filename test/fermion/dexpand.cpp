// ===== DE: extend_parity() and the D-expanding fermionic reload ===========
//
// Pins the "tensor_load with a larger virtual_dim" contract
// (docs/superpowers/specs/
//  2026-09-06-fermion-tensor-load-d-expansion-contract.md):
//   layer 1: extend_parity() as a pure function -- the six literal rows of
//            the R1 table, the refusal to shrink, and the invariant that the
//            entries already in the ledger are never moved nor flipped
//            (exhaustive over every ledger of length 1..6 and every target
//            dimension up to 8).
//   layer 2: load_fermion_ledger() through the iTPS constructor -- a D = 2
//            checkpoint read back with virtual_dim = 4 keeps its saved
//            entries, is padded exactly as the R1 rule says, passes
//            validate_neighbor_consistency(), and arrives with the tensors
//            and the Schmidt weights resized, the saved values in front and
//            the new slots zero; read back at the dimension it was saved at,
//            the ledger is unchanged.
//   layer 3: shrinking is still refused (a D = 4 checkpoint read back with
//            virtual_dim = 2 is a tenes::load_error naming the virtual
//            dimension, i.e. it is that guard and not another one firing).
//   layer 4: a checkpoint whose fermion.dat is longer than the tensors
//            params.dat describes is refused before any tensor file is read
//            (the ledger-versus-saved-shape check, with the marker method of
//            SL V5b).
// Every expected ledger is a literal: layer 1 copies the contract table,
// layer 2 uses values derived by hand from the contract rule. Nothing here
// compares the function under test against a second call of itself.
//
// The fixtures of test/fermion/saveload.cpp (namespace fermion_saveload) are
// reused; the tags are prefixed with "dexpand_" so the two files cannot
// collide in the shared output directory.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/exception.hpp"
#include "../../src/util/file.hpp"

namespace fermion_dexpand {

using tenes::fermion::parity_vector;
using tenes::itps::iTPSTestAccessor;
using lambda_table = std::vector<std::vector<std::vector<double>>>;

// Spelling ledgers as 0/1 keeps the R1 table readable side by side with the
// contract, where 0 = even and 1 = odd.
inline parity_vector pv(std::initializer_list<int> bits) {
  parity_vector p;
  p.reserve(bits.size());
  for (const int b : bits) {
    p.push_back(b != 0);
  }
  return p;
}

inline std::string show(const parity_vector &p) {
  std::string s = "{";
  for (std::size_t i = 0; i < p.size(); ++i) {
    s += (i == 0 ? "" : ", ");
    s += (p[i] ? '1' : '0');
  }
  return s + "}";
}

// Schmidt weights of the right length for a D = `dim` checkpoint; the values
// are irrelevant here (this file never compares measured quantities), they
// only have to exist so that save_tensors() can write lambda_*.dat.
inline lambda_table flat_lambda(int n_unit, std::size_t dim) {
  lambda_table lambda(n_unit, std::vector<std::vector<double>>(4));
  for (int site = 0; site < n_unit; ++site) {
    for (int leg = 0; leg < 4; ++leg) {
      for (std::size_t k = 0; k < dim; ++k) {
        lambda[site][leg].push_back(1.0 / (2.0 + site + leg + k));
      }
    }
  }
  return lambda;
}

// Writes a fermionic checkpoint whose every virtual leg carries `virt`; the
// bond dimension of the checkpoint is virt.size().
inline std::string save_checkpoint(const std::string &tag,
                                   const parity_vector &virt) {
  using namespace fermion_saveload;
  const std::string dir = save_dir_name(tag);
  reset_dir(case_dir_name(tag));
  const auto lattice = make_lattice(static_cast<int>(virt.size()));
  auto state = make_state(lattice, make_params(true, tag, dir, ""), false);
  inject(state, deterministic_state(virt), flat_lambda(sl_N_UNIT, virt.size()),
         virt);
  state.save_tensors();
  return dir;
}

}  // namespace fermion_dexpand

// ---------------------------------------------------------------- layer 1

TEST_CASE("DE layer1 extend_parity reproduces the contract table") {
  using namespace fermion_dexpand;
  using tenes::fermion::extend_parity;

  // The six rows of R1, copied literally from the contract.
  CHECK(extend_parity(pv({0, 1}), 4) == pv({0, 1, 0, 1}));
  CHECK(extend_parity(pv({0}), 4) == pv({0, 1, 0, 1}));
  CHECK(extend_parity(pv({0, 0, 0}), 4) == pv({0, 0, 0, 1}));
  CHECK(extend_parity(pv({1, 1}), 4) == pv({1, 1, 0, 0}));
  CHECK(extend_parity(pv({0, 1}), 2) == pv({0, 1}));
  CHECK(extend_parity(pv({0, 1}), 5) == pv({0, 1, 0, 1, 0}));

  // Anti-hollowness: four of the six rows differ from even_first_parity() at
  // the same length, so an implementation that simply reinvented the default
  // ledger (and thereby forgot the entries it must not touch) cannot pass.
  CHECK(pv({0, 1, 0, 1}) != tenes::fermion::even_first_parity(4));
  CHECK(pv({0, 0, 0, 1}) != tenes::fermion::even_first_parity(4));
  CHECK(pv({1, 1, 0, 0}) != tenes::fermion::even_first_parity(4));
  CHECK(pv({0, 1, 0, 1, 0}) != tenes::fermion::even_first_parity(5));
}

TEST_CASE("DE layer1 extend_parity refuses to shrink") {
  using namespace fermion_dexpand;
  using tenes::fermion::extend_parity;

  CHECK_THROWS_AS(extend_parity(pv({0, 1, 0, 1}), 3), std::runtime_error);
  CHECK_THROWS_AS(extend_parity(pv({0, 1}), 1), std::runtime_error);
  CHECK_THROWS_AS(extend_parity(pv({0}), 0), std::runtime_error);
  // The boundary the guard must not swallow: equal length is legal.
  CHECK_NOTHROW(extend_parity(pv({0, 1, 0, 1}), 4));
}

TEST_CASE("DE layer1 extend_parity never touches the existing entries") {
  using namespace fermion_dexpand;
  using tenes::fermion::extend_parity;

  constexpr std::size_t max_len = 6;
  constexpr std::size_t max_dim = 8;
  std::size_t checked = 0;
  for (std::size_t len = 1; len <= max_len; ++len) {
    for (std::size_t bits = 0; bits < (std::size_t(1) << len); ++bits) {
      parity_vector p(len, false);
      for (std::size_t i = 0; i < len; ++i) {
        p[i] = ((bits >> i) & 1u) != 0;
      }
      for (std::size_t new_dim = len; new_dim <= max_dim; ++new_dim) {
        INFO("p = " << show(p) << ", new_dim = " << new_dim);
        const auto q = extend_parity(p, new_dim);
        REQUIRE(q.size() == new_dim);
        for (std::size_t i = 0; i < len; ++i) {
          INFO("entry " << i << " of " << show(q));
          CHECK(q[i] == p[i]);
        }
        ++checked;
      }
    }
  }
  // Guards the loop bounds themselves: 2^1 + ... + 2^6 = 126 ledgers, and
  // ledger of length L is extended to L..8, i.e. 9 - L targets.
  std::size_t expected = 0;
  for (std::size_t len = 1; len <= max_len; ++len) {
    expected += (std::size_t(1) << len) * (max_dim - len + 1);
  }
  CHECK(checked == expected);
}

// ---------------------------------------------------------------- layer 2

namespace fermion_dexpand {

struct ExpandCase {
  const char *tag;
  parity_vector saved;   //!< ledger written into the checkpoint
  int new_dim;           //!< virtual_dim the reloading input asks for
  parity_vector expect;  //!< ledger the load must install, from the R1 rule
};

// Saves a checkpoint whose virtual legs all carry `c.saved`, reloads it with
// virtual_dim = c.new_dim, and checks the ledger and the tensors that come
// back. Each case is its own TEST_CASE below: an exception escaping one of
// them (which is what the unimplemented state does) then cannot hide the
// verdict of the others.
inline void check_reload(const ExpandCase &c) {
  using namespace fermion_saveload;
  INFO("case " << std::string(c.tag));
  // Anti-hollowness: the expected ledger is not the one a fresh run would
  // invent for this leg dimension, so a loader that kept the default (i.e.
  // never read fermion.dat) fails here.
  REQUIRE(c.expect != tenes::fermion::even_first_parity(
                          static_cast<std::size_t>(c.new_dim)));

  const std::string dir = save_checkpoint(c.tag, c.saved);
  REQUIRE(tenes::util::path_exists(dir + "/fermion.dat"));

  const auto wide = make_lattice(c.new_dim);
  const auto load_params =
      make_params(true, std::string(c.tag) + "_load", "", dir);
  auto loaded = make_state(wide, load_params, false);

  auto &finfo = iTPSTestAccessor::finfo(loaded);
  REQUIRE(finfo.enabled);
  REQUIRE(finfo.virt.size() == static_cast<std::size_t>(sl_N_UNIT));
  for (int site = 0; site < sl_N_UNIT; ++site) {
    INFO("site " << site);
    CHECK(finfo.phys[site] == sl_phys);
    for (int leg = 0; leg < 4; ++leg) {
      INFO("leg " << leg << " got " << show(finfo.virt[site][leg]));
      REQUIRE(finfo.virt[site][leg].size() ==
              static_cast<std::size_t>(c.new_dim));
      // The entries that were saved must come back untouched ...
      for (std::size_t i = 0; i < c.saved.size(); ++i) {
        CHECK(finfo.virt[site][leg][i] == c.saved[i]);
      }
      // ... and the padding must follow the R1 rule.
      CHECK(finfo.virt[site][leg] == c.expect);
    }
  }

  // The ledger the loader installed has to describe a consistent lattice.
  CHECK_NOTHROW(tenes::fermion::validate_neighbor_consistency(finfo, wide));

  // The tensors travelled with the ledger: every leg is at the new dimension
  // and the amplitudes the expansion introduced are exactly zero (F1/F3 of
  // the contract, and the reason the parity check survives any padding rule).
  auto &Tn = iTPSTestAccessor::Tn(loaded);
  REQUIRE(Tn.size() == static_cast<std::size_t>(sl_N_UNIT));
  double max_new = 0.0;
  double max_old = 0.0;
  for (int site = 0; site < sl_N_UNIT; ++site) {
    INFO("site " << site);
    const auto shape = Tn[site].shape();
    REQUIRE(shape.size() == 5);
    for (int leg = 0; leg < 4; ++leg) {
      CHECK(shape[leg] == static_cast<std::size_t>(c.new_dim));
    }
    CHECK(shape[4] == static_cast<std::size_t>(sl_pdim));
    for (std::size_t n = 0; n < Tn[site].local_size(); ++n) {
      const auto idx = Tn[site].global_index(n);
      double v = 0.0;
      Tn[site].get_value(idx, v);
      bool is_new = false;
      for (int leg = 0; leg < 4; ++leg) {
        if (static_cast<std::size_t>(idx[leg]) >= c.saved.size()) {
          is_new = true;
        }
      }
      if (is_new) {
        max_new = std::max(max_new, std::abs(v));
      } else {
        max_old = std::max(max_old, std::abs(v));
      }
    }
  }
  CHECK(max_new == 0.0);
  // Without this the zero check above would also hold for an all-zero
  // tensor, i.e. for a load that dropped the checkpoint on the floor.
  CHECK(max_old > 1.0e-8);

  // The Schmidt weights have to make the same trip: load_tensors_v1() reads
  // as many values per leg as the checkpoint holds and then resizes to the
  // new virtual_dim, so the saved weights stay in front and the new slots are
  // zero (the simple update turns a zero weight into a zero inverse through
  // Inverse_lambda_cut, which is what makes the padded state usable).
  // saved_lambda is the fixture save_checkpoint() wrote, i.e. an input of
  // this test, not a second reading of the code under test.
  const auto saved_lambda = flat_lambda(sl_N_UNIT, c.saved.size());
  auto &lambda = iTPSTestAccessor::lambda_tensor(loaded);
  REQUIRE(lambda.size() == static_cast<std::size_t>(sl_N_UNIT));
  for (int site = 0; site < sl_N_UNIT; ++site) {
    INFO("site " << site);
    REQUIRE(lambda[site].size() == 4u);
    for (int leg = 0; leg < 4; ++leg) {
      INFO("leg " << leg);
      REQUIRE(lambda[site][leg].size() == static_cast<std::size_t>(c.new_dim));
      for (std::size_t k = 0; k < c.saved.size(); ++k) {
        // The saved weights come back unchanged ...
        CHECK(lambda[site][leg][k] == saved_lambda[site][leg][k]);
        // ... and they are not zero, so the tail check below has content:
        // a load that zeroed lambda outright would pass it otherwise.
        CHECK(saved_lambda[site][leg][k] != 0.0);
      }
      for (std::size_t k = c.saved.size(); k < lambda[site][leg].size(); ++k) {
        CHECK(lambda[site][leg][k] == 0.0);
      }
    }
  }
}

}  // namespace fermion_dexpand

// The two expected ledgers below are worked out by hand from the R1 rule.
//   {0,1} -> length 3 wants ceil(3/2) = 2 evens and has 1, so an even is
//            appended; length 4 wants 2 and has 2, so an odd is appended.
//   {1,0} -> same counts, hence {1,0,0,1}.
TEST_CASE("DE layer2 a D=2 checkpoint with ledger 01 loads at virtual_dim 4") {
  using namespace fermion_dexpand;
  check_reload({"dexpand_e01", pv({0, 1}), 4, pv({0, 1, 0, 1})});
}

// {1,0} is deliberately not even_first_parity(2) = {0,1}: with the default
// ledger on both ends of the round trip, a loader that never reads
// fermion.dat would pass the previous case.
TEST_CASE("DE layer2 a D=2 checkpoint with ledger 10 loads at virtual_dim 4") {
  using namespace fermion_dexpand;
  check_reload({"dexpand_e10", pv({1, 0}), 4, pv({1, 0, 0, 1})});
}

TEST_CASE("DE layer2 loading at the saved virtual_dim is unchanged") {
  using namespace fermion_dexpand;
  check_reload({"dexpand_same", pv({1, 0}), 2, pv({1, 0})});
}

// ---------------------------------------------------------------- layer 3

TEST_CASE("DE layer3 a smaller virtual_dim is still refused") {
  using namespace fermion_dexpand;
  using namespace fermion_saveload;

  // Not even_first_parity(4) = {0,0,1,1}, so the refusal cannot be an
  // accident of the ledger happening to be the default one.
  const parity_vector wide_ledger = pv({0, 0, 0, 1});
  REQUIRE(wide_ledger != tenes::fermion::even_first_parity(4));
  const std::string dir = save_checkpoint("dexpand_shrink", wide_ledger);
  REQUIRE(tenes::util::path_exists(dir + "/fermion.dat"));

  const auto narrow = make_lattice(2);
  // CHECK_THROWS alone would stay green if some other guard (the parity
  // violation check, say) fired instead. The "ERROR: the virtual dimension of
  // the leg" prefix is matched in full on purpose: load_tensors_v1() also
  // prints "WARNING: virtual dimension of the leg ..." for every leg of a
  // D-changing load, so the bare phrase would not identify the guard (it is
  // what made run4 of FreeFermionSaveLoad hollow, review I-1).
  CHECK_THROWS_WITH_AS(
      make_state(narrow, make_params(true, "dexpand_shrink_load", "", dir),
                 false),
      doctest::Contains("ERROR: the virtual dimension of the leg"),
      tenes::load_error);
}

// ---------------------------------------------------------------- layer 4

TEST_CASE("DE layer4 a ledger longer than the saved tensors is refused") {
  using namespace fermion_dexpand;
  using namespace fermion_saveload;

  // A D = 2 checkpoint, which by itself loads fine at virtual_dim = 4 ...
  const std::string dir = save_checkpoint("dexpand_long_ledger", pv({0, 1}));
  const std::string path = dir + "/fermion.dat";
  REQUIRE(tenes::util::path_exists(path));

  const auto wide = make_lattice(4);
  auto loaded = make_state(
      wide, make_params(true, "dexpand_long_ledger_load", "", dir), false);

  // ... is overwritten with markers, so that a refusal arriving only after
  // the tensor files had been read would leave a trace.
  std::vector<tenes::real_tensor> markers;
  for (int site = 0; site < sl_N_UNIT; ++site) {
    markers.push_back(marker_tn(site));
  }
  const lambda_table marker_lambda(
      sl_N_UNIT,
      std::vector<std::vector<double>>(4, std::vector<double>(4, 3.5)));
  auto &state_Tn = iTPSTestAccessor::Tn(loaded);
  REQUIRE(state_Tn.size() == markers.size());
  for (std::size_t site = 0; site < markers.size(); ++site) {
    state_Tn[site] = markers[site];
  }
  iTPSTestAccessor::lambda_tensor(loaded) = marker_lambda;

  // Now stretch the ledger alone: fermion.dat claims four parities per virtual
  // leg while params.dat still says the tensors were saved at D = 2. The input
  // virtual_dim is 4, so the shrink guard cannot fire and the lengths the
  // ledger and the input agree on are the same -- only the check of the ledger
  // against the *saved* shape can catch this checkpoint.
  for (int site = 0; site < sl_N_UNIT; ++site) {
    for (int leg = 0; leg < 4; ++leg) {
      patch_line(path, virt_line(site, leg),
                 "0 1 0 1 # parity of the virtual leg " + std::to_string(leg) +
                     " of Tn[" + std::to_string(site) + "]");
    }
  }

  // The full ERROR prefix, not just "parity ledger": that phrase also occurs
  // in the neighbor-consistency message, and matching loosely is how run4 of
  // FreeFermionSaveLoad came to assert nothing (review I-1).
  CHECK_THROWS_WITH_AS(
      loaded.load_tensors(),
      doctest::Contains("ERROR: the virtual parity ledger of the leg"),
      tenes::load_error);

  // Nothing was read back: the markers are still there, at their D = 4 shape.
  for (int site = 0; site < sl_N_UNIT; ++site) {
    INFO("site " << site);
    const auto shape = state_Tn[site].shape();
    REQUIRE(shape.size() == 5);
    for (int leg = 0; leg < 4; ++leg) {
      CHECK(shape[leg] == 4u);
    }
    CHECK(shape[4] == static_cast<std::size_t>(sl_pdim));
    CHECK(max_abs_diff(state_Tn[site], markers[site]) == 0.0);
  }
  CHECK(iTPSTestAccessor::lambda_tensor(loaded) == marker_lambda);
}
