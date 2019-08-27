#define _USE_MATH_DEFINES
#include <random>
#include <sys/stat.h>

#include <cpptoml.h>

#include "Lattice.hpp"
#include "PEPS_Parameters.hpp"
#include "edge.hpp"
#include "util/read_matrix.hpp"
#include "tenes.hpp"

Lattice gen_lattice(decltype(cpptoml::parse_file("")) toml, const char *tablename = "lattice"){
  auto Lsub = toml->get_array_of<int64_t>("Lsub");
  if(!Lsub){
    std::cerr << "cannot find Lsub in the section [" << tablename << "]" << std::endl;
    // ERROR
  }
  return Lattice((*Lsub)[0], (*Lsub)[1]);
}

Edges gen_edges(decltype(cpptoml::parse_file("")) toml, const char *key, const char *tablename){
  auto str = toml->get_as<std::string>(key);
  if(!str){
    std::cerr << "cannot find " << key << " in the section [" << tablename << "]" << std::endl;
    // ERROR
  }
  return make_edges(*str);
}

template <typename tensor>
std::vector<tensor> gen_matrices(decltype(cpptoml::parse_file("")) toml, const char *key, const char *tablename){
  auto strs = toml->get_array_of<std::string>(key);
  if(!strs){
    std::cerr << "cannot find " << key << " in the section [" << tablename << "]" << std::endl;
    // ERROR
  }
  return util::read_matrix<tensor>(*strs);
}

