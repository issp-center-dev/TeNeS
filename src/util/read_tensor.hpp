#ifndef UTIL_READ_TENSOR_HPP
#define UTIL_READ_TENSOR_HPP

#include <algorithm>
#include <mptensor.hpp>
#include <sstream>
#include <string>
#include <vector>

#include "string.hpp"
#include "type_traits.hpp"
#include "../exception.hpp"

namespace tenes {

namespace util {

template <class ptensor>
ptensor read_tensor(std::string const &str, mptensor::Shape dims){
  using value_type = typename ptensor::value_type;
  ptensor ret(dims);
  const size_t rank = ret.rank();

  const static std::string delim = " \t";
  std::string line;

  std::stringstream ss(str);
  int linenum = 0;
  while (std::getline(ss, line)) {
    mptensor::Index index;
    index.resize(rank);
    auto fields = util::split(line);
    if(fields.empty()){
      ++linenum;
      continue;
    }

    if(fields.size() != rank+2){
      std::stringstream msg;
      msg << "cannot parse tensor; the number of columns differs from tensor rank + 2";
      std::stringstream ss2(str);
      int linenum2 = 0;
      while (std::getline(ss2, line)){
        msg << "\n" << line;
        if(linenum2 == linenum){
          msg << "\n";
          for(int i=0; i<line.length(); ++i){
            msg << "^";
          }
        }
        ++linenum2;
      }
      throw tenes::input_error(msg.str());
    }

    for(size_t i=0; i<rank; ++i){
      index[i] = std::stoi(fields[i]);
    }
    ret.set_value(index, convert_complex<value_type>(std::complex<double>(std::stod(fields[rank]), std::stod(fields[rank+1]))));
    ++linenum;
  }
  return ret;
}

}  // namespace util
}  // namespace tenes

#endif  // UTIL_READ_TENSOR_HPP
