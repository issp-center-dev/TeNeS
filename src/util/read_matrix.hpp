#ifndef UTIL_READ_MATRIX_HPP
#define UTIL_READ_MATRIX_HPP

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

#include <mptensor.hpp>

#include "string.hpp"

namespace util{
template <class ptensor>
ptensor read_matrix(std::string const& str) {
  std::vector<std::vector<double>> A;
  const static std::string delim = " \t";
  std::string line;

  std::stringstream ss(str);
  while(std::getline(ss, line)){
    auto fields = util::split(line);
    std::vector<double> row;
    for(auto field: fields){
      row.push_back(std::stof(field));
    }
    A.push_back(row);
  }

  const auto N = A.size();
  for(auto const& row: A){
    assert(row.size()==N);
  }

  ptensor ret(mptensor::Shape(N,N));
  for(int i=0; i<N; ++i){
    for(int j=0; j<N; ++j){
      ret.set_value(mptensor::Index(i,j), A[i][j]);
    }
  }
  return ret;
}

template <class ptensor>
std::vector<ptensor> read_matrix(std::vector<std::string> const& strs) {
  std::vector<ptensor> ret;
  std::transform(strs.begin(), strs.end(), std::back_inserter(ret), [](std::string s){return read_matrix<ptensor>(s);});
  return ret;
}

} // namespace util

#endif // UTIL_READ_MATRIX_HPP
