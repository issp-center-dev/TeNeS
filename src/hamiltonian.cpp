#include "hamiltonian.hpp"

#include <vector>

#include <mptensor.hpp>
#include <toml11/toml.hpp>

#include "util/load_tensor.hpp"

typedef mptensor::Tensor<mptensor::scalapack::Matrix, double> ptensor;

std::vector<ptensor> load_hamiltonians(toml::Table const& table){
  std::vector<ptensor> hams;
  auto ham_strs = toml::find<std::vector<std::string>>(toml::find(table, "bond"), "hamiltonian");
  for(auto str: ham_strs){
    hams.push_back(util::load_tensor(str));
  }
  return hams;
}
