#define _USE_MATH_DEFINES
#include <random>
#include <sys/stat.h>
#include <vector>
#include <utility>

#include <cpptoml.h>

#include "Lattice.hpp"
#include "correlation.hpp"
#include "PEPS_Parameters.hpp"
#include "edge.hpp"
#include "tenes.hpp"
#include "util/read_matrix.hpp"

Lattice gen_lattice(decltype(cpptoml::parse_file("")) toml,
                    const char *tablename = "lattice") {
  auto Lsub = toml->get_array_of<int64_t>("L_sub");
  if (!Lsub) {
    std::cerr << "cannot find Lsub in the section [" << tablename << "]"
              << std::endl;
    // ERROR
  }
  return Lattice((*Lsub)[0], (*Lsub)[1]);
}

Edges gen_edges(decltype(cpptoml::parse_file("")) toml, const char *key,
                const char *tablename) {
  auto str = toml->get_as<std::string>(key);
  if (!str) {
    std::cerr << "cannot find " << key << " in the section [" << tablename
              << "]" << std::endl;
    // ERROR
  }
  return make_edges(*str);
}

template <typename tensor>
std::vector<tensor> gen_matrices(decltype(cpptoml::parse_file("")) toml,
                                 const char *key, const char *tablename) {
  auto strs = toml->get_array_of<std::string>(key);
  if (!strs) {
    std::cerr << "cannot find " << key << " in the section [" << tablename
              << "]" << std::endl;
    // ERROR
  }
  return util::read_matrix<tensor>(*strs);
}

CorrelationParameter gen_corparam(decltype(cpptoml::parse_file("")) toml,
    const char *tablename = "correlation") {
  auto r = toml->get_as<int64_t>("r_max");
  if (!r){
    std::cerr << "cannot find r_max in the section [" << tablename
    << "]" << std::endl;
    // ERROR
  }
  int rmax = static_cast<int>(*r);

  auto oplist = toml->get_array_of<cpptoml::array>("operators");
  if (!oplist){
    std::cerr << "cannot find operators in the section [" << tablename
    << "]" << std::endl;
    // ERROR
  }

  std::vector<std::pair<int,int>> ops;
  for(auto op: *oplist){
    auto i = op->get_array_of<int64_t>();
    ops.push_back(std::make_pair(static_cast<int>((*i)[0]), static_cast<int>((*i)[1])));
  }

  return CorrelationParameter{rmax, ops};
}
