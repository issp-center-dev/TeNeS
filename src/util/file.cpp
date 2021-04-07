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

#include <cstddef>    // for size_t
#include <sys/stat.h>  // for stat, mkdir, S_ISDIR
#include <fstream>     // for string, operator<<, stringstream, basic_ostream
#include <string>      // for char_traits
#include <sstream>

#include "string.hpp"

#include "file.hpp"

using std::size_t;

namespace tenes {
namespace util {

bool path_exists(const std::string& path) {
  struct stat status;
  return stat(path.c_str(), &status) == 0;
}

bool isdir(const std::string& path) {
  struct stat status;
  int code = stat(path.c_str(), &status);
  if (code != 0) return false;
  return S_ISDIR(status.st_mode);
}

bool mkdir(const std::string& path) {
  int ret = 0;
  struct stat st;
  if (stat(path.c_str(), &st) != 0) {
    ret = ::mkdir(path.c_str(), 0755);
  }
  return ret == 0;
}

std::string joinpath(std::vector<std::string> const& xs) {
  const size_t n = xs.size();
  if (n == 0) {
    return std::string("");
  }
  std::stringstream ss;
  ss << xs[0];
  for (size_t i = 1; i < n; ++i) {
    ss << "/";
    ss << xs[i];
  }
  return ss.str();
}

std::string basename(const std::string& path) {
  auto xs = util::split(path, "/");
  return xs[xs.size() - 1];
}

}  // end of namespace util
}  // end of namespace tenes
