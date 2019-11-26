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


template <typename T>
inline T find_or(decltype(cpptoml::parse_file("")) param, const char *key,
                 T value) {
  return param->get_as<T>(key).value_or(value);
}

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


PEPS_Parameters gen_param(decltype(cpptoml::parse_file("")) param) {
  PEPS_Parameters pparam;

  // Tensor
  auto tensor = param->get_table("tensor");
  if(tensor != nullptr){
    pparam.D = find_or(tensor, "D", 2);
    pparam.CHI = find_or(tensor, "CHI", 4);
  }

  // Simple update
  auto simple = param->get_table("simple_update");
  if(simple != nullptr){
    pparam.num_simple_step = find_or(simple, "num_step", 0);
    pparam.Inverse_lambda_cut = find_or(simple, "inverse_lambda_cutoff", 1e-12);
  }

  // Full update
  auto full = param->get_table("full_update");
  if(full != nullptr){
    pparam.num_full_step = find_or(full, "num_step", 0);
    pparam.Full_Inverse_precision = find_or(full, "inverse_precision", 1e-12);
    pparam.Inverse_projector_cut =
        find_or(full, "inverse_projector_cutoff", 1e-12);
    pparam.Full_Convergence_Epsilon = find_or(full, "convergence_epsilon", 1e-12);
    pparam.Full_max_iteration = find_or(full, "iteration_max", 1000);
    pparam.Full_Gauge_Fix = find_or(full, "gauge_fix", true);
    pparam.Full_Use_FastFullUpdate = find_or(full, "fastfullupdate", true);
  }

  // Environment
  auto ctm = param->get_table("ctm");
  if(ctm != nullptr){
    pparam.Inverse_Env_cut = find_or(ctm, "inverse_projector_cutoff", 1e-12);
    pparam.CTM_Convergence_Epsilon = find_or(ctm, "convergence_epsilon", 1e-10);
    pparam.Max_CTM_Iteration = find_or(ctm, "iteration_max", 100);
    pparam.CTM_Projector_corner = find_or(ctm, "projector_corner", false);
    pparam.Use_RSVD = find_or(ctm, "use_rsvd", false);
    pparam.RSVD_Oversampling_factor = find_or(ctm, "rsvd_oversampling_factor", 2);
  }

  // random
  auto random = param->get_table("random");
  if(random != nullptr){
    pparam.seed = find_or(random, "seed", 11);
  }

  return pparam;
}
