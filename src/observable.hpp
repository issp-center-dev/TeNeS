#ifndef OBSERVABLE_HPP
#define OBSERVABLE_HPP

#include <vector>
#include <mptensor.hpp>
#include <toml11/toml/value.hpp>

using ptensor = mptensor::Tensor<mptensor::scalapack::Matrix, double>;

std::vector<ptensor> load_local_operators(const char* filename);
std::vector<ptensor> load_local_operators(toml::value const& table);

#endif // OBSERVABLE_HPP
