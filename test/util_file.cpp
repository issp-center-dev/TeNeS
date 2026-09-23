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

// Unit tests for the filesystem helpers of src/util/file.*.
//
// The one thing pinned here that no other test can reach is what
// util::remove_all() answers when the deletion did not happen. Its callers
// (src/iTPS/saveload_tensors.cpp) branch on the return value: the save reports
// success and drops the .tenes-save-incomplete marker on the strength of it,
// so a "true" that only means "I could not find out" turns into a checkpoint
// that the next load refuses.

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <unistd.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>

#include "../src/mpi.hpp"
#include "../src/util/file.hpp"
#include "test_workdir.hpp"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  doctest::Context context(argc, argv);
  const int res = context.run();
  MPI_Finalize();
  return res;
}

namespace {

namespace fs = std::filesystem;

//! Give back owner read/write/search on everything under path, top down.
//! os.walk-style traversal skips what it cannot list, and a locked directory
//! is exactly what these cases leave behind, so the recursion is by hand.
void unlock_tree(const fs::path &path) {
  std::error_code ec;
  const fs::file_status st = fs::symlink_status(path, ec);
  if (ec || !fs::exists(st)) {
    return;
  }
  if (fs::is_symlink(st)) {
    return;
  }
  fs::permissions(path, fs::perms::owner_all, fs::perm_options::add, ec);
  if (fs::is_directory(st)) {
    fs::directory_iterator it(path, ec);
    if (ec) {
      return;
    }
    const fs::directory_iterator end;
    for (; it != end; it.increment(ec)) {
      if (ec) {
        return;
      }
      unlock_tree(it->path());
    }
  }
}

//! Set a path's permission bits for the lifetime of the object and put the
//! previous ones back however the body leaves - a failing REQUIRE throws, and
//! a directory left at mode 0 would take the rest of the run with it.
class locked_mode {
 public:
  locked_mode(fs::path path, fs::perms mode) : path_(std::move(path)) {
    std::error_code ec;
    previous_ = fs::status(path_, ec).permissions();
    if (ec) {
      previous_ = fs::perms::owner_all;
    }
    fs::permissions(path_, mode, fs::perm_options::replace, ec);
  }
  ~locked_mode() {
    std::error_code ec;
    fs::permissions(path_, previous_, fs::perm_options::replace, ec);
  }
  locked_mode(const locked_mode &) = delete;
  locked_mode &operator=(const locked_mode &) = delete;

 private:
  fs::path path_;
  fs::perms previous_;
};

//! A fresh directory to work in, gone again when the case ends.
class scratch_dir {
 public:
  explicit scratch_dir(const std::string &name) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    path_ =
        fs::path("util_file_scratch") / (name + "_rank" + std::to_string(rank));
    clear();
    std::error_code ec;
    fs::create_directories(path_, ec);
  }
  ~scratch_dir() { clear(); }
  scratch_dir(const scratch_dir &) = delete;
  scratch_dir &operator=(const scratch_dir &) = delete;

  const fs::path &path() const { return path_; }

 private:
  void clear() {
    unlock_tree(path_);
    std::error_code ec;
    fs::remove_all(path_, ec);
    // The shared parent, only once the last case has left it empty.
    fs::remove(path_.parent_path(), ec);
  }
  fs::path path_;
};

void write_file(const fs::path &path, const std::string &contents) {
  std::ofstream f(path);
  f << contents;
}

//! Why the "could not be removed" case cannot be arranged here, or an empty
//! string when it can. Running as root, or on a filesystem that ignores the
//! mode bits, makes an undeletable path impossible to build: the case then
//! proves nothing and says so rather than passing quietly.
std::string permission_skip_reason(const fs::path &scratch) {
  if (::geteuid() == 0) {
    return "running as root: permission bits are not enforced";
  }
  const fs::path probe = scratch / "permission_probe";
  std::error_code ec;
  fs::create_directories(probe, ec);
  if (ec) {
    return "could not create the probe directory: " + ec.message();
  }
  write_file(probe / "inside.txt", "x\n");
  bool readable_while_locked = false;
  {
    locked_mode lock(probe, fs::perms::none);
    std::ifstream f(probe / "inside.txt");
    readable_while_locked = static_cast<bool>(f);
  }
  unlock_tree(probe);
  fs::remove_all(probe, ec);
  if (readable_while_locked) {
    return "this filesystem does not enforce directory permissions";
  }
  return "";
}

}  // namespace

TEST_CASE("remove_all answers whether the path is gone") {
  using tenes::util::remove_all;
  scratch_dir scratch("remove_all");

  SUBCASE("an existing file is removed and reported as removed") {
    const fs::path victim = scratch.path() / "a_file.dat";
    write_file(victim, "content\n");
    REQUIRE(fs::exists(victim));

    CHECK(remove_all(victim.string()));
    CHECK_FALSE(fs::exists(victim));
  }

  SUBCASE("an existing directory tree is removed and reported as removed") {
    const fs::path victim = scratch.path() / "a_tree";
    fs::create_directories(victim / "nested");
    write_file(victim / "nested" / "inside.dat", "content\n");
    REQUIRE(fs::exists(victim / "nested" / "inside.dat"));

    CHECK(remove_all(victim.string()));
    CHECK_FALSE(fs::exists(victim));
  }

  SUBCASE("a path that was never there is reported as gone") {
    const fs::path victim = scratch.path() / "never_existed";
    REQUIRE_FALSE(fs::exists(victim));

    CHECK(remove_all(victim.string()));
  }

  // The contract: a deletion that did not happen must not be reported as one.
  // Taking search permission off the parent makes the path below it
  // undeletable AND unstattable at the same time, which is the combination
  // that "the path is not there, as far as I can tell" gets wrong.
  SUBCASE("a path that could not be removed is reported as not removed") {
    const std::string reason = permission_skip_reason(scratch.path());
    if (!reason.empty()) {
      MESSAGE("SKIP (cannot make an undeletable path here): " << reason);
    } else {
      const fs::path parent = scratch.path() / "locked_parent";
      const fs::path victim = parent / "victim";
      fs::create_directories(victim);
      write_file(victim / "inside.dat", "content\n");
      REQUIRE(fs::exists(victim / "inside.dat"));

      bool answer = false;
      {
        locked_mode lock(parent, fs::perms::none);
        answer = remove_all(victim.string());
      }

      // The premise of the case: nothing was actually deleted. Checked after
      // the permissions are back, because checking it through the locked
      // parent would fail for the same reason remove_all() did.
      REQUIRE(fs::exists(victim));
      REQUIRE(fs::exists(victim / "inside.dat"));

      INFO("remove_all(" << victim.string() << ") returned " << answer
                         << " while the path and its contents are still "
                            "there");
      CHECK_FALSE(answer);
    }
  }
}

TEST_CASE("path_exists says false when the answer cannot be determined") {
  // Not a contract clause: the documented behaviour of path_exists() is what
  // makes remove_all()'s "!path_exists(path)" the wrong question, so it is
  // pinned here to keep the two from being changed apart. A path_exists()
  // that threw or that reported "true, I could not tell" would need
  // remove_all() rewritten either way.
  scratch_dir scratch("path_exists");
  const std::string reason = permission_skip_reason(scratch.path());
  if (!reason.empty()) {
    MESSAGE("SKIP (cannot make an unstattable path here): " << reason);
  } else {
    const fs::path parent = scratch.path() / "locked_parent";
    const fs::path hidden = parent / "hidden.dat";
    fs::create_directories(parent);
    write_file(hidden, "content\n");
    REQUIRE(fs::exists(hidden));

    bool answer = true;
    {
      locked_mode lock(parent, fs::perms::none);
      answer = tenes::util::path_exists(hidden.string());
    }
    REQUIRE(fs::exists(hidden));
    CHECK_FALSE(answer);
  }
}
