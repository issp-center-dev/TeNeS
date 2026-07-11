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

namespace tenes::util {

bool path_exists(const std::string& path);
bool isdir(const std::string& path);

//! Create a directory (and its parents if necessary).
//! @return true if the directory exists when this function returns.
bool mkdir(const std::string& path);

std::string basename(const std::string& path);

}  // end of namespace tenes::util

#endif  // TENES_SRC_UTIL_FILE_HPP_
