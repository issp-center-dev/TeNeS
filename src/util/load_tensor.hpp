#ifndef UTIL_LOAD_TENSOR_HPP
#define UTIL_LOAD_TENSOR_HPP

#include <string>
#include <mptensor.hpp>

namespace util{

using ptensor = mptensor::Tensor<mptensor::scalapack::Matrix, double>;
ptensor load_tensor(std::string const& str);

} // namespace util

#endif // UTIL_LOAD_TENSOR_HPP
