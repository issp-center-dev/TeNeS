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

#ifndef UTIL_STRING_HPP
#define UTIL_STRING_HPP

#include "string.hpp"

#include <string>
#include <vector>

namespace tenes {
namespace util {

std::vector<std::string> split(std::string const &str,
                               std::string const &delim = " \t") {
  std::vector<std::string> words;
  auto index = 0;
  auto last = str.size();
  while (index != std::string::npos) {
    auto next = str.find_first_of(delim, index);
    if(next == std::string::npos){
      words.push_back(str.substr(index, last - index + 1));
      break;
    }
    words.push_back(str.substr(index, next - index + 1));
    index = str.find_first_not_of(delim, next);
  }
  return words;
}

std::string lstrip(std::string const &str, std::string const &delim = " \t\n") {
  std::string ret("");
  auto index = str.find_first_not_of(delim);
  if (index != std::string::npos) {
    ret = str.substr(index, str.size() - index + 1);
  }
  return ret;
}

std::string rstrip(std::string const &str, std::string const &delim = " \t\n") {
  std::string ret("");
  auto index = str.find_last_not_of(delim);
  if (index != std::string::npos) {
    ret = str.substr(0, index + 1);
  }
  return ret;
}

std::string strip(std::string const &str, std::string const &delim = " \t\n") {
  return rstrip(lstrip(str, delim), delim);
}

std::string drop_comment(std::string const &str) {
  auto index = str.find_first_of("#");
  return str.substr(0, index);
}

}  // namespace util
}  // namespace tenes

#endif  // UTIL_STRING_HPP
