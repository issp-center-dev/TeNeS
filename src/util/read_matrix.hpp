#ifndef UTIL_READ_MATRIX_HPP
#define UTIL_READ_MATRIX_HPP

#include <algorithm>
#include <mptensor.hpp>
#include <sstream>
#include <string>
#include <vector>

#include "string.hpp"

namespace tenes {

namespace util {
template <class ptensor>
ptensor read_matrix(std::string const &str) {
  std::vector<std::vector<double>> A;
  const static std::string delim = " \t";
  std::string line;

  std::stringstream ss(str);
  while (std::getline(ss, line)) {
    auto fields = util::split(line);
    std::vector<double> row;
    for (auto field : fields) {
      row.push_back(std::stof(field));
    }
    A.push_back(row);
  }

  const auto Nrow = A.size();
  const auto Ncol = A[0].size();
  for (auto const &row : A) {
    assert(row.size() == Ncol);
  }

  ptensor ret(mptensor::Shape(Nrow, Ncol));
  for (int i = 0; i < Nrow; ++i) {
    for (int j = 0; j < Ncol; ++j) {
      ret.set_value(mptensor::Index(i, j), A[i][j]);
    }
  }
  return ret;
}

template <class ptensor>
std::vector<ptensor> read_matrix(std::vector<std::string> const &strs) {
  std::vector<ptensor> ret;
  std::transform(strs.begin(), strs.end(), std::back_inserter(ret),
                 [](std::string s) { return read_matrix<ptensor>(s); });
  return ret;
}


template <class ptensor>
ptensor read_tensor(std::string const &str, mptensor::Shape dims){
  ptensor ret(dims);
  const size_t rank = ret.rank();

  const static std::string delim = " \t";
  std::string line;

  std::stringstream ss(str);
  while (std::getline(ss, line)) {
    mptensor::Index index;
    index.resize(rank);
    auto fields = util::split(line);
    assert(fields.size() == rank+2);
    for(size_t i=0; i<rank; ++i){
      index[i] = std::stoi(fields[i]);
    }
    ret.set_value(index, std::stof(fields[rank]));
  }
  return ret;
}

}  // namespace util
}  // namespace tenes

#endif  // UTIL_READ_MATRIX_HPP
