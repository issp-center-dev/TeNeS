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

#ifndef TENES_SRC_UTIL_READ_TENSOR_HPP_
#define TENES_SRC_UTIL_READ_TENSOR_HPP_

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <mptensor/mptensor.hpp>

#include "string.hpp"
#include "type_traits.hpp"
#include "../exception.hpp"

namespace tenes {

namespace util {

template <class ptensor>
ptensor read_tensor(std::string const &str, mptensor::Shape dims,
                    double atol = 0.0) {
  using value_type = typename ptensor::value_type;
  ptensor ret(dims);
  const size_t rank = ret.rank();

  static const std::string delim = " \t";
  std::string line;

  std::stringstream ss(str);
  int linenum = 0;
  while (std::getline(ss, line)) {
    mptensor::Index index;
    index.resize(rank);

    line = util::drop_comment(line);
    line = util::strip(line);
    auto fields = util::split(line);
    if (fields.empty()) {
      ++linenum;
      continue;
    }

    if (fields.size() != rank + 2) {
      std::stringstream msg;
      msg << "cannot parse tensor; the number of columns differs from tensor "
             "rank + 2";
      std::stringstream ss2(str);
      int linenum2 = 0;
      while (std::getline(ss2, line)) {
        msg << "\n" << line;
        if (linenum2 == linenum) {
          msg << "\n";
          for (int i = 0; i < line.length(); ++i) {
            msg << "^";
          }
        }
        ++linenum2;
      }
      throw tenes::input_error(msg.str());
    }

    for (size_t i = 0; i < rank; ++i) {
      index[i] = std::stoi(fields[i]);
    }
    double re = std::stod(fields[rank]);
    re = (std::abs(re) >= atol) ? re : 0.0;
    double im = std::stod(fields[rank + 1]);
    im = (std::abs(im) >= atol) ? im : 0.0;
    ret.set_value(index,
                  convert_complex<value_type>(std::complex<double>(re, im)));
    ++linenum;
  }
  return ret;
}

}  // namespace util
}  // namespace tenes

#endif  // TENES_SRC_UTIL_READ_TENSOR_HPP_
