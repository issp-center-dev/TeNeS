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

#include "load_toml.hpp"

#define _USE_MATH_DEFINES
#include <algorithm>
#include <array>
#include <iterator>
#include <string>
#include <tuple>
#include <typeinfo>
#include <vector>

#include <boost/core/demangle.hpp>

#include "../exception.hpp"
#include "../util/read_tensor.hpp"
#include "../util/string.hpp"
#include "../tensor.hpp"

#include "PEPS_Parameters.hpp"

namespace tenes {
namespace itps {

namespace detail {
std::string msg_cannot_find(std::string key, std::string section = "") {
  std::stringstream ss;
  if (section.empty()) {
    ss << "cannot find \"" << key << "\"";
  } else {
    ss << "cannot find \"" << key << "\" in a section \"" << section << "\"";
  }
  return ss.str();
}

template <class T>
std::string type_name() {
  return boost::core::demangle(typeid(T).name());
}

template <>
std::string type_name<std::string>() {
  return "string";
}

template <class T>
struct bittype {
  typedef T type;
};

template <>
struct bittype<short int> {
  typedef int64_t type;
};

template <>
struct bittype<int> {
  typedef int64_t type;
};

template <>
struct bittype<long long int> {
  typedef int64_t type;
};

}  // end of namespace detail

template <class T>
std::vector<T> get_array_of(decltype(cpptoml::parse_file("")) param,
                            const char *key);

template <class T>
std::vector<T> get_array_of(decltype(cpptoml::parse_file("")) param,
                            const char *key, T default_value);

template <class T>
std::vector<T> get_array_of(decltype(cpptoml::parse_file("")) param,
                            const char *key,
                            std::vector<T> const &default_value);

/*! @brief Load a parameter if key is defined.
 *
 * @param dst Destination of the parameter.
 * @param param TOML file object.
 * @param key Key of the parameter.
 *
 * @throw input_error if the parameter is found but the type is not matched.
 */
template <class T>
inline void load_if(T &dst, decltype(cpptoml::parse_file("")) param,
                    const char *key) {
  if (param->contains(key)) {
    auto v = param->get_as<T>(key);
    if (v) {
      dst = *v;
    } else {
      auto vv = param->get(key);
      std::stringstream ss;
      ss << "\"" << key << "\" requires a value of type "
         << detail::type_name<T>() << ", but given value is " << *vv;
      throw input_error(ss.str());
    }
  }
}

/*! @brief Load a parameter if key is defined.
 *
 * @param dst Destination of the parameter.
 * @param param TOML file object.
 * @param key Key of the parameter.
 *
 * @throw input_error if the parameter is found but the type is not matched.
 */
template <class T>
inline void load_if(std::vector<T> &dst,
                    decltype(cpptoml::parse_file("")) param, const char *key) {
  if (param->contains(key)) {
    auto xs = get_array_of<typename detail::bittype<T>::type>(param, key);
    dst.assign(xs.begin(), xs.end());
  }
}

/*! @brief
 * @param param TOML file object.
 * @param key Key of the parameter.
 * @param value Default value.
 * @return Value of the parameter if the parameter is defined, otherwise default
 * is returned.
 * @throw input_error if the parameter is found but the type is not matched.
 */
template <class T>
inline T find_or(decltype(cpptoml::parse_file("")) param, const char *key,
                 T value) {
  T ret = value;
  if (param->contains(key)) {
    load_if(ret, param, key);
  }
  return ret;
}

/*! @brief Find a parameter from TOML object.
 *  @param param TOML file object.
 *  @param key Key of the parameter.
 *  @return Value of the parameter.
 *  @throw input_error if the parameter is not found.
 */
template <class T>
inline T find(decltype(cpptoml::parse_file("")) param, const char *key) {
  T ret;
  if (param->contains(key)) {
    load_if(ret, param, key);
  } else {
    throw input_error(detail::msg_cannot_find(key));
  }
  return ret;
}

/*! @brief Find a list parameter from TOML object.
 *
 * If the parameter is scalar, it is converted to a list with single element.
 *
 * @param param TOML file object.
 * @param key Key of the parameter.
 * @return Value of the parameter.
 * @throw input_error if the parameter is not found.
 */
template <class T>
std::vector<T> get_array_of(decltype(cpptoml::parse_file("")) param,
                            const char *key) {
  if (!param->contains(key)) {
    throw input_error(detail::msg_cannot_find(key));
  }
  auto scalar = param->get_as<T>(key);
  auto arr = param->get_array_of<T>(key);
  std::vector<T> ret;
  if (scalar) {
    ret.push_back(*scalar);
  } else if (arr) {
    ret.assign(arr->begin(), arr->end());
  } else {
    auto vv = param->get(key);
    std::stringstream ss;
    ss << "\"" << key << "\" requires value(s) of type "
       << detail::type_name<T>() << ", but given value is " << *vv;
    throw input_error(ss.str());
  }
  return ret;
}

/*! @brief Find a list parameter from TOML object.
 *
 * If the parameter is scalar, it is converted to a list with single element.
 *
 * @param param TOML file object.
 * @param key Key of the parameter.
 * @param default_value Default value.
 * @return Value of the parameter.
 */
template <class T>
std::vector<T> get_array_of(decltype(cpptoml::parse_file("")) param,
                            const char *key, T default_value) {
  if (!param->contains(key)) {
    return std::vector<T>{default_value};
  } else {
    auto xs = get_array_of<typename detail::bittype<T>::type>(param, key);
    std::vector<T> ret(xs.begin(), xs.end());
    return ret;
  }
}

template <class T>
std::vector<T> get_array_of(decltype(cpptoml::parse_file("")) param,
                            const char *key,
                            std::vector<T> const &default_value) {
  if (!param->contains(key)) {
    std::vector<T> ret = default_value;
    return ret;
  } else {
    auto xs = get_array_of<typename detail::bittype<T>::type>(param, key);
    std::vector<T> ret(xs.begin(), xs.end());
    return ret;
  }
}

template <>
std::vector<int32_t> get_array_of<int32_t>(
    decltype(cpptoml::parse_file("")) param, const char *key) {
  auto xs = get_array_of<int64_t>(param, key);
  std::vector<int32_t> ret(xs.begin(), xs.end());
  return ret;
}

SquareLattice gen_lattice(decltype(cpptoml::parse_file("")) toml,
                          const char *tablename) {
  auto Lsub = toml->get_array_of<int64_t>("L_sub");
  if (!Lsub) {
    throw input_error(detail::msg_cannot_find("L_sub", tablename));
  }

  auto skew = find_or<int>(toml, "skew", 0);

  SquareLattice lat((*Lsub)[0], (*Lsub)[1], skew);

  auto sites = toml->get_table_array("unitcell");
  if (!sites) {
    throw input_error(detail::msg_cannot_find("unitcell", tablename));
  }
  for (const auto &site : *sites) {
    auto indices = get_array_of<int64_t>(site, "index");
    if (indices.empty()) {
      for (int i = 0; i < lat.N_UNIT; ++i) {
        indices.push_back(i);
      }
    }

    for (int index : indices) {
      lat.initial_dirs[index] =
          get_array_of<double>(site, "initial_state", 0.0);

      lat.noises[index] = find_or(site, "noise", 0.0);

      lat.physical_dims[index] = find<int>(site, "physical_dim");

      auto vdim = get_array_of<int64_t>(site, "virtual_dim");
      if (vdim.size() == 1) {
        vdim.resize(4, vdim[0]);
      }
      if (vdim.size() != 4) {
        throw input_error(
            "The size of tensor.unitcell.virtual_dim must be 1 or 4.");
      }

      for (int i = 0; i < 4; ++i) {
        lat.virtual_dims[index][i] = vdim[i];
      }
    }  // end of for indices
  }

  lat.check_dims();

  return lat;
}

CorrelationParameter gen_corparam(decltype(cpptoml::parse_file("")) toml,
                                  const char *tablename) {
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

TransferMatrix_Parameters gen_transfer_matrix_parameter(
    decltype(cpptoml::parse_file("")) toml, const char *tablename) {
  TransferMatrix_Parameters clength;

  load_if(clength.to_calculate, toml, "measure");
  load_if(clength.num_eigvals, toml, "num_eigvals");
  load_if(clength.maxdim_dense_eigensolver, toml, "maxdim_dense_eigensolver");
  load_if(clength.arnoldi_maxdim, toml, "arnoldi_maxdim");
  load_if(clength.arnoldi_restartdim, toml, "arnoldi_restartdim");
  load_if(clength.arnoldi_maxiter, toml, "arnoldi_maxiterations");
  load_if(clength.arnoldi_rtol, toml, "arnoldi_rtol");
  return clength;
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

    std::string mode_str =
        find_or(general, "mode", std::string("ground state"));
    if (util::startswith(mode_str, "ground")) {
      pparam.calcmode = PEPS_Parameters::CalculationMode::ground_state;
    } else if (util::startswith(mode_str, "time")) {
      pparam.calcmode = PEPS_Parameters::CalculationMode::time_evolution;
    } else if (util::startswith(mode_str, "finite")) {
      pparam.calcmode = PEPS_Parameters::CalculationMode::finite_temperature;
    } else {
      throw input_error("Invalid mode: " + mode_str);
    }
    load_if(pparam.measure_interval, general, "measure_interval");
  }

  // Simple update
  auto simple = param->get_table("simple_update");
  if (simple != nullptr) {
    load_if(pparam.num_simple_step, simple, "num_step");
    load_if(pparam.tau_simple_step, simple, "tau");
    load_if(pparam.Inverse_lambda_cut, simple, "lambda_cutoff");
    load_if(pparam.Simple_Gauge_Fix, simple, "gauge_fix");
    load_if(pparam.Simple_Gauge_maxiter, simple, "gauge_maxiter");
    load_if(pparam.Simple_Gauge_Convergence_Epsilon, simple,
            "gauge_convergence_epsilon");
  }

  // Full update
  auto full = param->get_table("full_update");
  if (full != nullptr) {
    load_if(pparam.num_full_step, full, "num_step");
    load_if(pparam.tau_full_step, simple, "tau");
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

  // Contraction

  auto contraction = param->get_table("contraction");
  if (contraction != nullptr) {
    std::string contraction_mode_str =
        find_or(contraction, "optimize", std::string("automatic"));
    if (contraction_mode_str == "automatic") {
      pparam.cpath_opt = PEPS_Parameters::CPathOptimization::automatic;
    } else if (contraction_mode_str == "never") {
      pparam.cpath_opt = PEPS_Parameters::CPathOptimization::never;
    } else if (contraction_mode_str == "always") {
      pparam.cpath_opt = PEPS_Parameters::CPathOptimization::always;
    } else if (contraction_mode_str == "old") {
      pparam.cpath_opt = PEPS_Parameters::CPathOptimization::old;
    } else {
      throw input_error("Invalid contraction optimize: " + contraction_mode_str);
    }

    pparam.contraction_path_file =
        find_or(contraction, "pathfile", std::string(""));
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
  if (words.size() < 3) {
    throw input_error("A bond should have 3 integers");
  }
  return std::make_tuple(stoi(words[0]), stoi(words[1]), stoi(words[2]));
}

std::tuple<int, std::vector<int>, std::vector<int>> read_multisites(
    std::string line) {
  using std::stoi;
  auto words = util::split(util::strip(line));
  if (words.size() < 5) {
    throw input_error("A multisite should have at least 5 integers");
  }
  if (words.size() % 2 != 1) {
    throw input_error("A multisite should have odd number of integers");
  }

  std::vector<int> dx, dy;
  for (int i = 1; i < words.size(); i += 2) {
    dx.push_back(stoi(words[i]));
    dy.push_back(stoi(words[i + 1]));
  }
  return std::make_tuple(stoi(words[0]), dx, dy);
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

std::vector<std::tuple<int, std::vector<int>, std::vector<int>>>
read_multisitesset(std::string str) {
  std::vector<std::tuple<int, std::vector<int>, std::vector<int>>> ret;
  std::string line;
  std::stringstream ss(str);
  while (std::getline(ss, line)) {
    line = util::strip(util::drop_comment(line));
    if (util::strip(line).empty()) {
      continue;
    }
    auto bond = read_multisites(line);
    ret.push_back(bond);
  }
  return ret;
}

template <class tensor>
Operators<tensor> load_operator(decltype(cpptoml::parse_file("")) param,
                                MPI_Comm comm, int nsites, int nbody,
                                double atol, const char *tablename) {
  auto elements = param->get_as<std::string>("elements");
  auto ops = param->get_array_of<int64_t>("ops");
  if (elements && ops) {
    std::stringstream ss;
    ss << "Both elements and ops are defined in a section " << tablename;
    throw tenes::input_error(ss.str());
  }
  tensor A(comm);
  std::vector<int> op_ind;
  if (nbody == 1 && !elements) {
    std::stringstream ss;
    ss << "elements not found in a section " << tablename;
    throw tenes::input_error(ss.str());
  }
  if (elements) {
    if (nbody > 2) {
      std::stringstream ss;
      ss << "observable.multisite does not support elements";
      throw tenes::input_error(ss.str());
    }
    mptensor::Shape shape;
    auto dim_arr = param->get_array_of<int64_t>("dim");
    auto dim_int = param->get_as<int>("dim");
    if (dim_arr) {
      if (dim_arr->size() != nbody) {
        std::stringstream ss;
        ss << "operator is " << nbody << "-sites but dim has "
           << dim_arr->size() << " integers";
        throw input_error(ss.str());
      }
      for (int d : *dim_arr) {
        shape.push(d);
      }
    } else if (dim_int) {
      for (int i = 0; i < nbody; ++i) {
        shape.push(*dim_int);
      }
    } else {
      throw input_error(detail::msg_cannot_find("dim", tablename));
    }
    for (int i = 0; i < nbody; ++i) {
      shape.push(shape[i]);
    }
    A = util::read_tensor<tensor>(*elements, shape, comm, atol);
  } else if (ops) {
    op_ind.assign(ops->begin(), ops->end());
  } else {
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
  } else if (nbody == 2) {
    auto bonds_str = find<std::string>(param, "bonds");
    auto bonds = read_bonds(bonds_str);
    for (auto bond : bonds) {
      if (elements) {
        ret.emplace_back(name, group, std::get<0>(bond), std::get<1>(bond),
                         std::get<2>(bond), A);
      } else {
        ret.emplace_back(name, group, std::get<0>(bond), std::get<1>(bond),
                         std::get<2>(bond), op_ind);
      }
    }
  } else {
    auto ms_str = find<std::string>(param, "multisites");
    auto mss = read_multisitesset(ms_str);
    for (auto ms : mss) {
      if (elements) {
        ret.emplace_back(name, group, std::get<0>(ms), std::get<1>(ms),
                         std::get<2>(ms), A);
      } else {
        ret.emplace_back(name, group, std::get<0>(ms), std::get<1>(ms),
                         std::get<2>(ms), op_ind);
      }
    }
  }
  return ret;
}

template <class tensor>
Operators<tensor> load_operators(decltype(cpptoml::parse_file("")) param,
                                 MPI_Comm comm, int nsites, int nbody,
                                 double atol, std::string const &key) {
  Operators<tensor> ret;
  auto tables = param->get_table_array_qualified(key);
  if (!tables) {
    std::cout << "INFO: " << key << " is not found in the input file. (Skipped)"
              << std::endl;
    return ret;
  }
  for (const auto &table : *tables) {
    auto obs =
        load_operator<tensor>(table, comm, nsites, nbody, atol, key.c_str());
    std::copy(obs.begin(), obs.end(), std::back_inserter(ret));
    // std::move(obs.begin(), obs.end(), std::back_inserter(ret));
  }
  return ret;
}

template <class tensor>
EvolutionOperator<tensor> load_Evolution_operator(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, double atol,
    const char *tablename) {
  auto dimensions = param->get_array_of<int64_t>("dimensions");
  if (!dimensions) {
    throw input_error(detail::msg_cannot_find("dimensions", tablename));
  }
  auto elements = find<std::string>(param, "elements");
  auto shape = mptensor::Shape();
  for (auto d : *dimensions) {
    shape.push(d);
  }
  auto group = find_or<int>(param, "group", 0);

  if (shape.size() == 2) {
    // siteoperator
    auto site = find<int>(param, "site");

    tensor A = util::read_tensor<tensor>(elements, shape, comm, atol);
    return make_onesite_EvolutionOperator<tensor>(site, group, A);
  } else if (shape.size() == 4) {
    // nnoperator
    auto source_site = find<int>(param, "source_site");
    auto source_leg = find<int>(param, "source_leg");

    tensor A = util::read_tensor<tensor>(elements, shape, comm, atol);
    return make_twosite_EvolutionOperator<tensor>(source_site, source_leg,
                                                  group, A);
  } else {
    std::stringstream ss;
    ss << tablename << ".dimensions should have 2 or 4 integers";
    throw input_error(ss.str());
  }
}

template <class tensor>
EvolutionOperators<tensor> load_updates(decltype(cpptoml::parse_file("")) param,
                                        MPI_Comm comm, double atol,
                                        std::string const &key) {
  EvolutionOperators<tensor> ret;
  auto tables = param->get_table_array_qualified(key);
  for (const auto &table : *tables) {
    ret.push_back(
        load_Evolution_operator<tensor>(table, comm, atol, key.c_str()));
  }
  return ret;
}
template <class tensor>
EvolutionOperators<tensor> load_simple_updates(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, double atol) {
  return load_updates<tensor>(param, comm, atol, "evolution.simple");
}
template <class tensor>
EvolutionOperators<tensor> load_full_updates(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, double atol) {
  return load_updates<tensor>(param, comm, atol, "evolution.full");
}

template <class tensor>
std::map<TNC_map_key, TensorNetworkContractor<tensor>> load_contraction_paths(
    std::string const &path_file) {
  std::map<TNC_map_key, TensorNetworkContractor<tensor>> ret;

  auto input_toml = cpptoml::parse_file(path_file);
  auto params_toml = input_toml->get_table("params");
  if (!params_toml) {
    throw input_error(detail::msg_cannot_find("params", path_file));
  }
  const bool is_tpo = find<bool>(params_toml, "is_TPO");
  const bool is_mf = find<bool>(params_toml, "is_mf");

  auto tensor_toml = input_toml->get_table("tensor");
  if (!tensor_toml) {
    throw input_error(detail::msg_cannot_find("tensor", path_file));
  }
  const auto L_sub = get_array_of<int32_t>(tensor_toml, "L_sub");
  if (L_sub.size() != 2) {
    throw input_error("L_sub should have 2 integers");
  }
  const int LX = L_sub[0];
  const int LY = L_sub[1];
  const int N_UNIT = LX * LY;

  const auto shape_type = get_array_of<int32_t>(tensor_toml, "shape_types");
  if (shape_type.size() != N_UNIT) {
    throw input_error("shape_types should have " + std::to_string(N_UNIT) +
                      " integers");
  }
  int num_types_shape = 0;
  for (auto s : shape_type) {
    num_types_shape = std::max(num_types_shape, s);
  }
  ++num_types_shape;

  const auto shape_dims_o =
      tensor_toml->get_array_of<cpptoml::array>("shape_dims");
  if (shape_dims_o->size() != num_types_shape) {
    throw input_error("shape_dims should have " +
                      std::to_string(num_types_shape) + " lists");
  }
  const auto nlegs = is_tpo ? 6 : 5;
  std::vector<std::vector<int32_t>> shape_dims(num_types_shape);
  for (int i = 0; i < num_types_shape; ++i) {
    auto v = (*shape_dims_o)[i]->get_array_of<int64_t>();
    shape_dims[i].assign(v->begin(), v->end());
    if (shape_dims[i].size() != nlegs) {
      throw input_error("shape_dims[" + std::to_string(i) + "] should have " +
                        std::to_string(nlegs) + " integers");
    }
  }

  auto contractions = input_toml->get_table_array_qualified("contraction");
  if (!contractions) {
    throw input_error(detail::msg_cannot_find("contraction", path_file));
  }
  for (const auto &ctable : *contractions) {
    const int nrow = find<int>(ctable, "nrow");
    const int ncol = find<int>(ctable, "ncol");
    const int nsites = nrow * ncol;
    const auto shape_types = get_array_of<int32_t>(ctable, "shape_types");
    if (shape_types.size() != nsites) {
      throw input_error("shape_types should have " + std::to_string(nsites) +
                        " integers");
    }
    TensorNetworkContractor<tensor> tnc(nrow, ncol, is_tpo, is_mf);
    const auto cpath = get_array_of<int32_t>(ctable, "contraction_path");
    tnc.set_path(cpath);

    ret[TNC_map_key(nrow, ncol, shape_types)] = tnc;
  }

  return ret;
}

// template instantiations

template Operators<real_tensor> load_operator(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, int nsites,
    int nbody, double atol, const char *tablename);
template Operators<complex_tensor> load_operator(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, int nsites,
    int nbody, double atol, const char *tablename);

template Operators<real_tensor> load_operators(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, int nsites,
    int nbody, double atol, std::string const &key);

template Operators<complex_tensor> load_operators(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, int nsites,
    int nbody, double atol, std::string const &key);

template EvolutionOperator<real_tensor> load_Evolution_operator(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, double atol,
    const char *tablename);
template EvolutionOperator<complex_tensor> load_Evolution_operator(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, double atol,
    const char *tablename);

template EvolutionOperators<real_tensor> load_updates(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, double atol,
    std::string const &key);
template EvolutionOperators<complex_tensor> load_updates(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, double atol,
    std::string const &key);

template EvolutionOperators<real_tensor> load_simple_updates(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, double atol);
template EvolutionOperators<complex_tensor> load_simple_updates(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, double atol);

template EvolutionOperators<real_tensor> load_full_updates(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, double atol);
template EvolutionOperators<complex_tensor> load_full_updates(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, double atol);

template std::map<TNC_map_key, TensorNetworkContractor<real_tensor>>
load_contraction_paths(std::string const &path_file);
template std::map<TNC_map_key, TensorNetworkContractor<complex_tensor>>
load_contraction_paths(std::string const &path_file);

// end of template instantiations

}  // namespace itps
}  // namespace tenes
