#include "observable.hpp"

#include <vector>

#include <mptensor.hpp>
#include <toml11/toml.hpp>

#include "util/load_tensor.hpp"

typedef mptensor::Tensor<mptensor::scalapack::Matrix, double> ptensor;

std::vector<ptensor> load_local_operators(toml::Table const& table){
  std::vector<ptensor> lops;
  auto lops_strs = toml::find<std::vector<std::string>>(toml::find(table, "observable"), "local_operator");
  for(auto str: lops_strs){
    lops.push_back(util::load_tensor(str));
  }
  return lops;
}
