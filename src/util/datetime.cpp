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

#include "datetime.hpp"

namespace tenes {
namespace util {

std::string datetime(std::time_t const &t) {
  std::tm *lt = std::localtime(&t);
  char dt[64];
  std::strftime(dt, 64, "%FT%T", lt);
  char tz[7];
  std::strftime(tz, 7, "%z", lt);
  if (tz[0] != 'Z' && tz[0] != 'z') {
    for (size_t i = 6; i > 3; --i) {
      tz[i] = tz[i - 1];
    }
    tz[3] = ':';
  }
  std::string ret = dt;
  return ret + tz;
}
std::string datetime() {
  std::time_t t = std::time(nullptr);
  return datetime(t);
}

}  // end of namespace util

}  // end of namespace tenes
