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

#ifndef TENES_SRC_VERSION_HPP_
#define TENES_SRC_VERSION_HPP_

#include <string>

namespace tenes {
//! version number, as written in the top-level CMakeLists.txt
const char *version();
//! full hash of the commit the executable was built from
//! ("" when it is not known)
const char *git_hash();
//! whether the source tree had uncommitted changes
bool git_dirty();
//! "<version> (<commit>)", where <commit> is the first 8 digits of the hash,
//! followed by "-dirty" for a tree with uncommitted changes;
//! "<version>" when the commit is not known
std::string version_string();
}  // namespace tenes

#define TENES_VERSION (::tenes::version())

#endif // TENES_SRC_VERSION_HPP_
