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

#ifndef UTIL_FILE_HPP
#define UTIL_FILE_HPP

#include <vector>
#include <string>

namespace tenes {
namespace util{

bool path_exists(const std::string& path);
bool isdir(const std::string& path);
bool mkdir(const std::string& path);
std::string joinpath(std::vector<std::string> const& xs);
std::string basename(const std::string& path);

}  // end of namespace util
}  // end of namespace tenes

#endif  // UTIL_FILE_HPP
