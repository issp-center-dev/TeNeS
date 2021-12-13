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

#ifndef TENES_SRC_UTIL_DATETIME_HPP_
#define TENES_SRC_UTIL_DATETIME_HPP_

#include <ctime>  // IWYU pragma: export
#include <string> // IWYU pragma: export

namespace tenes {
namespace util {

std::string datetime(std::time_t const &t);
std::string datetime();

}  // end of namespace util

}  // end of namespace tenes

#endif // TENES_SRC_UTIL_DATETIME_HPP_
