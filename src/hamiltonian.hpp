#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

#include <vector>
#include <mptensor.hpp>
#include <toml11/toml/value.hpp>

using ptensor = mptensor::Tensor<mptensor::scalapack::Matrix, double>;

std::vector<ptensor> load_hamiltonians(const char* filename);
std::vector<ptensor> load_hamiltonians(toml::value const& table);

#endif // HAMILTONIAN_HPP
