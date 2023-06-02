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

#ifndef TENES_SRC_UTIL_STRING_HPP_
#define TENES_SRC_UTIL_STRING_HPP_

#include <string>
#include <vector>

namespace tenes {
namespace util {

std::vector<std::string> split(std::string const &str,
                               std::string const &delim = " \t");

std::string lstrip(std::string const &str, std::string const &delim = " \t\n");

std::string rstrip(std::string const &str, std::string const &delim = " \t\n");

std::string strip(std::string const &str, std::string const &delim = " \t\n");

std::string drop_comment(std::string const &str);

bool startswith(std::string const &str, std::string const &prefix);
bool endswith(std::string const &str, std::string const &suffix);

}  // namespace util
}  // namespace tenes

#endif  // TENES_SRC_UTIL_STRING_HPP_
