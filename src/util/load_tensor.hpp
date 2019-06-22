#ifndef UTIL_LOAD_TENSOR_HPP
#define UTIL_LOAD_TENSOR_HPP

#include <vector>
#include <string>
#include <sstream>

#include <mptensor.hpp>
// #include <toml11/toml.hpp>

#include "string.hpp"

typedef mptensor::Tensor<mptensor::scalapack::Matrix, double> ptensor;

namespace util{
ptensor load_tensor(std::string const& str) {
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

  ptensor ret(mptensor::Shape(4,4));
  for(int i=0; i<N; ++i){
    for(int j=0; j<N; ++j){
      ret.set_value(mptensor::Index(i,j), A[i][j]);
    }
  }
  return ret;
}

} // namespace util

#endif // UTIL_LOAD_TENSOR_HPP
