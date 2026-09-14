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

#include <filesystem>
#include <string>
#include <system_error>
#include <vector>

#include "file.hpp"

namespace tenes::util {

namespace fs = std::filesystem;

bool path_exists(const std::string& path) {
  std::error_code ec;
  return fs::exists(path, ec);
}

bool isdir(const std::string& path) {
  std::error_code ec;
  return fs::is_directory(path, ec);
}

bool mkdir(const std::string& path) {
  std::error_code ec;
  fs::create_directories(path, ec);
  // isdir(), not fs::is_directory(): the throwing overload is unreachable here
  // because !ec short-circuits, but this file's rule is that nothing in it
  // throws -- its callers sit in rank-0-only blocks ahead of a collective.
  return !ec && isdir(path);
}

bool remove_all(const std::string& path) {
  std::error_code ec;
  fs::remove_all(path, ec);
  return !path_exists(path);
}

bool rename(const std::string& from, const std::string& to) {
  std::error_code ec;
  fs::rename(from, to, ec);
  return !ec;
}

std::string basename(const std::string& path) {
  return fs::path(path).filename().string();
}

std::string absolute_path(const std::string& path) {
  std::error_code ec;
  const fs::path abs = fs::absolute(path, ec);
  if (ec || abs.empty()) {
    return path;
  }
  return abs.lexically_normal().string();
}

bool is_empty_directory(const std::string& path) {
  std::error_code ec;
  if (!fs::is_directory(path, ec) || ec) {
    return false;
  }
  const bool empty = fs::is_empty(path, ec);
  return !ec && empty;
}

std::vector<std::string> entries_with_prefix(const std::string& directory,
                                             const std::string& prefix) {
  std::vector<std::string> ret;
  std::error_code ec;
  fs::directory_iterator it(directory, ec);
  if (ec) {
    return ret;
  }
  // increment(ec) rather than ++it: the range-for form throws when the
  // directory changes under us, and a listing failure here must not take
  // down the caller.
  const fs::directory_iterator end;
  for (; it != end; it.increment(ec)) {
    if (ec) {
      break;
    }
    const std::string name = it->path().filename().string();
    if (name.size() >= prefix.size() &&
        name.compare(0, prefix.size(), prefix) == 0) {
      ret.push_back(it->path().string());
    }
  }
  return ret;
}

}  // end of namespace tenes::util
