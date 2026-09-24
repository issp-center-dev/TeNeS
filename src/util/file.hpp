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

#ifndef TENES_SRC_UTIL_FILE_HPP_
#define TENES_SRC_UTIL_FILE_HPP_

#include <string>
#include <vector>

namespace tenes::util {

//! Whether path exists. False, rather than an exception, when the answer
//! cannot be determined (no search permission on a parent, say): these are
//! called from rank-0-only blocks that go on to a collective, where an
//! exception on one rank alone hangs the run.
bool path_exists(const std::string& path);
//! Whether path is a directory; false when it cannot be determined.
bool isdir(const std::string& path);

//! Create a directory (and its parents if necessary).
//! @return true if the directory exists when this function returns.
bool mkdir(const std::string& path);

//! Delete a file, or a directory and everything under it.
//! @return true if the path is gone when this function returns.
bool remove_all(const std::string& path);

//! Move a file or directory. Both paths must be on the same filesystem for
//! this to be the cheap metadata operation the checkpoint swap relies on.
//! @return true on success.
bool rename(const std::string& from, const std::string& to);

std::string basename(const std::string& path);

//! Absolute form of path, or path itself when one cannot be formed. Used where
//! a path is recorded for a person to act on later, possibly from a different
//! working directory, so a relative spelling would be ambiguous.
std::string absolute_path(const std::string& path);

//! Full paths of the entries of directory whose names begin with prefix.
//! Empty, or short, if the directory cannot be read in full: a caller that
//! must know every entry was listed has to confirm it some other way, for
//! example with is_empty_directory() after acting on each one.
std::vector<std::string> entries_with_prefix(const std::string& directory,
                                             const std::string& prefix);

//! Whether path is a directory with no entries. False when it has entries,
//! is not a directory, or cannot be read.
bool is_empty_directory(const std::string& path);

}  // end of namespace tenes::util

#endif  // TENES_SRC_UTIL_FILE_HPP_
