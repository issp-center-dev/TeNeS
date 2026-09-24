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

#ifndef TENES_TEST_TEST_WORKDIR_HPP_
#define TENES_TEST_TEST_WORKDIR_HPP_

#include <filesystem>
#include <system_error>

namespace tenes {
namespace test {

// The unit tests use relative paths: they write output directories and result
// files into the current directory, and full_update / simple_update read their
// reference data from ./data. ctest runs them in the build tree, where both
// belong. Running a test binary by hand from somewhere else - typically the
// repository root - used to litter that directory with output_* directories
// and to fail to find the reference data.
//
// Moving to the build tree at start-up makes a direct run behave exactly like
// a ctest run. Under ctest the destination is already the current directory,
// so this is a no-op there.
//
// Include this header in any test that touches the filesystem; the
// initialization below runs before main.
inline bool move_to_test_workdir() {
#ifdef TENES_TEST_WORKDIR
  std::error_code ec;
  std::filesystem::current_path(TENES_TEST_WORKDIR, ec);
#endif
  return true;
}

namespace {
const bool test_workdir_ready = move_to_test_workdir();
}  // namespace

}  // namespace test
}  // namespace tenes

#endif  // TENES_TEST_TEST_WORKDIR_HPP_
