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

#include "file.hpp"

namespace tenes::util {

namespace fs = std::filesystem;

bool path_exists(const std::string& path) { return fs::exists(path); }

bool isdir(const std::string& path) { return fs::is_directory(path); }

bool mkdir(const std::string& path) {
  std::error_code ec;
  fs::create_directories(path, ec);
  return !ec && fs::is_directory(path);
}

std::string basename(const std::string& path) {
  return fs::path(path).filename().string();
}

}  // end of namespace tenes::util
