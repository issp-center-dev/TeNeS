/* TeNeS - Massively parallel tensor network solver /
/ Copyright (C) 2019- The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify /
/ it under the terms of the GNU General Public License as published by /
/ the Free Software Foundation, either version 3 of the License, or /
/ (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, /
/ but WITHOUT ANY WARRANTY; without even the implied warranty of /
/ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the /
/ GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License /
/ along with this program. If not, see http://www.gnu.org/licenses/. */

#define _USE_MATH_DEFINES
#include <random>
#include <sys/stat.h>
#include <tuple>
#include <utility>
#include <vector>

#include <typeinfo>
#include <boost/core/demangle.hpp>

#include <cpptoml.h>

#include "Lattice.hpp"
#include "PEPS_Parameters.hpp"
#include "correlation.hpp"
#include "operator.hpp"
#include "exception.hpp"
#include "util/read_tensor.hpp"
#include "util/string.hpp"

namespace tenes {

namespace detail{
std::string msg_cannot_find(std::string key, std::string section=""){
  std::stringstream ss;
  if(section.empty()){
    ss << "cannot find \"" << key << "\"";
  }else{
    ss << "cannot find \"" << key << "\" in a section \"" << section << "\"";
  }
  return ss.str();
}

template <class T>
std::string type_name(){
  return boost::core::demangle(typeid(T).name());
}

template <>
std::string type_name<std::string>(){
  return "string";
}

} // end of namespace detail


template <class T>
inline void load_if(T& dst, decltype(cpptoml::parse_file("")) param, const char *key){
  if(param->contains(key)){
    auto v = param->get_as<T>(key);
    if (v){
      dst = *v;
    }else{
      auto vv = param->get(key);
      std::stringstream ss;
      ss << "\"" << key << "\" requires a value of type " << detail::type_name<T>() << ", but given value is " << *vv;
      throw input_error(ss.str());
    }
  }
}

template <class T>
inline T find_or(decltype(cpptoml::parse_file("")) param, const char *key,
                 T value) {
  T ret = value;
  if(param->contains(key)){
    load_if(ret, param, key);
  }
  return ret;
}

template <class T>
inline T find(decltype(cpptoml::parse_file("")) param, const char *key){
  T ret;
  if(param->contains(key)){
    load_if(ret, param, key);
  }else{
    throw input_error(detail::msg_cannot_find(key));
  }
  return ret;
}


Lattice gen_lattice(decltype(cpptoml::parse_file("")) toml,
                    const char *tablename = "tensor") {
  auto Lsub = toml->get_array_of<int64_t>("L_sub");
  if (!Lsub) {
    throw input_error(detail::msg_cannot_find("L_sub", tablename));
  }

  auto skew = find_or<int>(toml, "skew", 0);

  Lattice lat((*Lsub)[0], (*Lsub)[1], skew);

  auto sites = toml->get_table_array("unitcell");
  if (!sites) {
    throw input_error(detail::msg_cannot_find("unitcell", tablename));
  }
  for (const auto &site : *sites) {
    std::vector<int> indices;
    auto index_int = site->get_as<int>("index");
    auto index_arr = site->get_array_of<int64_t>("index");

    if (index_arr) {
      indices.assign(index_arr->begin(), index_arr->end());
      if (indices.empty()) {
        for (int i = 0; i < lat.N_UNIT; ++i) {
          indices.push_back(i);
        }
      }
    } else if (index_int) {
      indices.push_back(*index_int);
    } else {
      throw input_error(detail::msg_cannot_find("index", tablename));
    }

    for (int index : indices) {
      auto dir = site->get_array_of<double>("initial_state");
      if (dir) {
        lat.initial_dirs[index] = *dir;
      } else {
        lat.initial_dirs[index] = std::vector<double>{0.0};
      }

      lat.noises[index] = find_or(site, "noise", 0.0);

      lat.physical_dims[index] = find<int>(site, "physical_dim");

      auto vdim_int = site->get_as<int>("virtual_dim");
      auto vdim_arr = site->get_array_of<int64_t>("virtual_dim");
      std::vector<int> vdim;

      if (vdim_arr) {
        vdim.assign((*vdim_arr).begin(), (*vdim_arr).end());
      } else if (vdim_int) {
        vdim.assign(4, *vdim_int);
      }else{
        throw input_error(detail::msg_cannot_find("virtual_dim", "tensor.unitcell"));
      }
      for (int i = 0; i < 4; ++i) {
        lat.virtual_dims[index][i] = vdim[i];
      }
    } // end of for indices
  }

  lat.check_dims();

  return lat;
}

CorrelationParameter gen_corparam(decltype(cpptoml::parse_file("")) toml,
                                  const char *tablename = "correlation") {
  int rmax = find<int>(toml, "r_max");

  auto oplist = toml->get_array_of<cpptoml::array>("operators");
  if (!oplist) {
    throw input_error(detail::msg_cannot_find("operators", tablename));
  }

  std::vector<std::tuple<int, int>> ops;
  for (auto op : *oplist) {
    auto i = op->get_array_of<int64_t>();
    ops.emplace_back(static_cast<int>((*i)[0]), static_cast<int>((*i)[1]));
  }

  return CorrelationParameter{rmax, ops};
}

PEPS_Parameters gen_param(decltype(cpptoml::parse_file("")) param) {
  PEPS_Parameters pparam;

  // general
  auto general = param->get_table("general");
  if (general != nullptr) {
    load_if(pparam.is_real, general, "is_real");
    load_if(pparam.iszero_tol, general, "iszero_tol");
    load_if(pparam.to_measure, general, "measure");
    load_if(pparam.outdir, general, "output");
    load_if(pparam.tensor_load_dir, general, "tensor_load");
    load_if(pparam.tensor_save_dir, general, "tensor_save");
  }

  // Simple update
  auto simple = param->get_table("simple_update");
  if (simple != nullptr) {
    load_if(pparam.num_simple_step, simple, "num_step");
    load_if(pparam.Inverse_lambda_cut, simple, "lambda_cutoff");
    // load_if(pparam.Simple_Gauge_Fix, simple, "gauge_fix");
  }

  // Full update
  auto full = param->get_table("full_update");
  if (full != nullptr) {
    load_if(pparam.num_full_step, full, "num_step");
    load_if(pparam.Full_Inverse_precision, full, "inverse_precision");
    load_if(pparam.Full_Convergence_Epsilon, full, "convergence_epsilon");
    load_if(pparam.Inverse_Env_cut, full, "env_cutoff");
    load_if(pparam.Full_max_iteration, full, "iteration_max");
    load_if(pparam.Full_Gauge_Fix, full, "gauge_fix");
    load_if(pparam.Full_Use_FastFullUpdate, full, "fastfullupdate");
  }

  // Environment
  auto ctm = param->get_table("ctm");
  if (ctm != nullptr) {
    load_if(pparam.CHI, ctm, "dimension");
    load_if(pparam.Inverse_projector_cut, ctm, "projector_cutoff");
    load_if(pparam.CTM_Convergence_Epsilon, ctm, "convergence_epsilon");
    load_if(pparam.Max_CTM_Iteration, ctm, "iteration_max");
    load_if(pparam.CTM_Projector_corner, ctm, "projector_corner");
    load_if(pparam.Use_RSVD, ctm, "use_rsvd");
    load_if(pparam.RSVD_Oversampling_factor, ctm, "rsvd_oversampling_factor");
    load_if(pparam.MeanField_Env, ctm, "meanfield_env");

    if (pparam.RSVD_Oversampling_factor < 1.0) {
      std::string msg = "rsvd_oversampling_factor must be >= 1.0";
      throw tenes::input_error(msg);
    }
  }

  // random
  auto random = param->get_table("random");
  if (random != nullptr) {
    load_if(pparam.seed, random, "seed");
  }

  return pparam;
}

std::tuple<int, int, int> read_bond(std::string line) {
  using std::stoi;
  auto words = util::split(util::strip(line));
  if(words.size() < 3){
    throw input_error("A bond should have 3 integers");
  }
  return std::make_tuple(stoi(words[0]), stoi(words[1]), stoi(words[2]));
}
std::vector<std::tuple<int, int, int>> read_bonds(std::string str) {
  std::vector<std::tuple<int, int, int>> ret;
  std::string line;
  std::stringstream ss(str);
  while (std::getline(ss, line)) {
    line = util::strip(util::drop_comment(line));
    if (util::strip(line).empty()) {
      continue;
    }
    auto bond = read_bond(line);
    ret.push_back(bond);
  }
  return ret;
}

template <class tensor>
Operators<tensor> load_operator(decltype(cpptoml::parse_file("")) param,
                                int nsites, int nbody, double atol=0.0,
                                const char* tablename="observable.onesite") {
  assert(nbody == 1 || nbody == 2);

  auto elements = param->get_as<std::string>("elements");
  auto ops = param->get_array_of<int64_t>("ops");
  if(elements && ops){
    std::stringstream ss;
    ss << "Both elements and ops are defined in a section " << tablename;
    throw tenes::input_error(ss.str());
  }
  tensor A;
  std::vector<int> op_ind;
  if(nbody==1 && !elements){
    std::stringstream ss;
    ss << "elements not found in a section " << tablename;
    throw tenes::input_error(ss.str());
  }
  if(elements){
    mptensor::Shape shape;
    auto dim_arr = param->get_array_of<int64_t>("dim");
    auto dim_int = param->get_as<int>("dim");
    if (dim_arr) {
      if(dim_arr->size() != nbody){
        std::stringstream ss;
        ss << "operator is " << nbody << "-sites but dim has " << dim_arr->size() << " integers";
        throw input_error(ss.str());
      }
      for (int d : *dim_arr) {
        shape.push(d);
      }
    } else if (dim_int) {
      for (int i = 0; i < nbody; ++i) {
        shape.push(*dim_int);
      }
    } else{
      throw input_error(detail::msg_cannot_find("dim", tablename));
    }
    for (int i = 0; i < nbody; ++i) {
      shape.push(shape[i]);
    }
    A = util::read_tensor<tensor>(*elements, shape, atol);
  }else if(ops){
    op_ind.assign(ops->begin(), ops->end());
  }else{
    std::stringstream ss;
    ss << "Neither elements nor ops are not defined in a section " << tablename;
    throw tenes::input_error(ss.str());
  }

  auto group = find<int>(param, "group");
  auto name = find_or(param, "name", std::string(""));

  std::vector<Operator<tensor>> ret;

  if (nbody == 1) {
    auto site_int = param->get_as<int>("sites");
    auto site_arr = param->get_array_of<int64_t>("sites");
    std::vector<int> sites;
    if (site_arr) {
      sites.assign(site_arr->begin(), site_arr->end());
      if (sites.empty()) {
        for (int i = 0; i < nsites; ++i) {
          sites.push_back(i);
        }
      }
    } else if (site_int) {
      sites.push_back(*site_int);
    } else {
      throw input_error(detail::msg_cannot_find("sites", tablename));
    }
    for (int s : sites) {
      ret.emplace_back(name, group, s, A);
    }
  } else { // nbody == 2
    auto bonds_str = find<std::string>(param, "bonds");
    auto bonds = read_bonds(bonds_str);
    for (auto bond : bonds) {
      if(elements){
        ret.emplace_back(name, group, std::get<0>(bond), std::get<1>(bond),
                         std::get<2>(bond), A);
      }else{
        ret.emplace_back(name, group, std::get<0>(bond), std::get<1>(bond),
                         std::get<2>(bond), op_ind);
      }
    }
  }
  return ret;
}

template <class tensor>
Operators<tensor> load_operators(decltype(cpptoml::parse_file("")) param,
                                 int nsites, int nbody,
                                 double atol,
                                 std::string const &key
                                 ) {
  Operators<tensor> ret;
  auto tables = param->get_table_array_qualified(key);
  for (const auto &table : *tables) {
    auto obs = load_operator<tensor>(table, nsites, nbody, atol, key.c_str());
    std::copy(obs.begin(), obs.end(), std::back_inserter(ret));
    // std::move(obs.begin(), obs.end(), std::back_inserter(ret));
  }
  return ret;
}

template <class tensor>
NNOperator<tensor> load_nn_operator(decltype(cpptoml::parse_file("")) param, double atol=0.0, const char* tablename="evolution.simple") {
  auto source_site = find<int>(param, "source_site");
  auto source_leg = find<int>(param, "source_leg");
  auto dimensions = param->get_array_of<int64_t>("dimensions");
  if(!dimensions){
    throw input_error(detail::msg_cannot_find("dimensions", tablename));
  }
  auto elements = find<std::string>(param, "elements");
  auto shape = mptensor::Shape();
  for (auto d : *dimensions) {
    shape.push(d);
  }
  if(shape.size() != 4){
    std::stringstream ss;
    ss << tablename << ".dimensions should have 4 integers";
    throw input_error(ss.str());
  }
  tensor A = util::read_tensor<tensor>(elements, shape, atol);
  return NNOperator<tensor>(source_site, source_leg, A);
}

template <class tensor>
NNOperators<tensor> load_updates(decltype(cpptoml::parse_file("")) param,
                                 double atol,
                                 std::string const &key) {
  NNOperators<tensor> ret;
  auto tables = param->get_table_array_qualified(key);
  for (const auto &table : *tables) {
    ret.push_back(load_nn_operator<tensor>(table, atol, key.c_str()));
  }
  return ret;
}
template <class tensor>
NNOperators<tensor>
load_simple_updates(decltype(cpptoml::parse_file("")) param, double atol=0.0) {
  return load_updates<tensor>(param, atol, "evolution.simple");
}
template <class tensor>
NNOperators<tensor> load_full_updates(decltype(cpptoml::parse_file("")) param, double atol=0.0) {
  return load_updates<tensor>(param, atol, "evolution.full");
}

} // end of namespace tenes
