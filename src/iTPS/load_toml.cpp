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
#include <type_traits>
#include <vector>

#include "../exception.hpp"
#include "../util/read_tensor.hpp"
#include "../util/string.hpp"
#include "../tensor.hpp"

#include "PEPS_Parameters.hpp"

namespace tenes::itps {

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
const char *type_name() {
  if constexpr (std::is_same_v<T, bool>) {
    return "boolean";
  } else if constexpr (std::is_integral_v<T>) {
    return "integer";
  } else if constexpr (std::is_floating_point_v<T>) {
    return "float";
  } else {
    return "string";
  }
}

/*! @brief Convert a TOML value into T.
 *
 * As a special case, an integer value may be read as a floating-point
 * number (e.g., tau = 1 is accepted where a float is expected).
 * Any other type mismatch is an error.
 *
 * @throw input_error if the value has a wrong type.
 */
template <class T>
T as_value(const toml::value &v) {
  if constexpr (std::is_same_v<T, bool>) {
    if (v.is_boolean()) {
      return static_cast<T>(v.as_boolean());
    }
  } else if constexpr (std::is_integral_v<T>) {
    if (v.is_integer()) {
      return static_cast<T>(v.as_integer());
    }
  } else {
    if (v.is_floating()) {
      return static_cast<T>(v.as_floating());
    }
    if (v.is_integer()) {
      return static_cast<T>(v.as_integer());
    }
  }
  throw input_error(toml::format_error(
      std::string("[error] a value of type ") + type_name<T>() + " is required",
      v, "given here"));
}

template <>
std::string as_value<std::string>(const toml::value &v) {
  if (v.is_string()) {
    return v.as_string();
  }
  throw input_error(toml::format_error(
      "[error] a value of type string is required", v, "given here"));
}

/*! @brief Convert a TOML array into std::vector<T>.
 *
 * @throw input_error if the value is not an array or an element has a wrong
 * type.
 */
template <class T>
std::vector<T> as_array(const toml::value &v) {
  if (!v.is_array()) {
    throw input_error(toml::format_error(
        std::string("[error] an array of ") + type_name<T>() + " is required",
        v, "given here"));
  }
  std::vector<T> ret;
  ret.reserve(v.as_array().size());
  for (const auto &elem : v.as_array()) {
    ret.push_back(as_value<T>(elem));
  }
  return ret;
}

/*! @brief Follow a dotted key (e.g. "observable.onesite") from root.
 *
 * @return Pointer to the found value, or nullptr if some component is
 * missing.
 */
const toml::value *find_path(const toml::value &root,
                             const std::string &dotted_key) {
  const toml::value *current = &root;
  for (const auto &key : util::split(dotted_key, ".")) {
    if (!current->is_table() || !current->contains(key)) {
      return nullptr;
    }
    current = &current->at(key);
  }
  return current;
}

/*! @brief Get a subtable if defined.
 *
 * @return Pointer to the subtable, or nullptr if key is not defined.
 * @throw input_error if the value is not a table.
 */
const toml::value *opt_table(const toml::value &param, const char *key) {
  if (!param.contains(key)) {
    return nullptr;
  }
  const toml::value &v = param.at(key);
  if (!v.is_table()) {
    throw input_error(toml::format_error(
        std::string("[error] \"") + key + "\" should be a table", v,
        "given here"));
  }
  return &v;
}

template <class T>
T reduce(double re, double im) {
  if constexpr (std::is_floating_point_v<T>) {
    return re;
  } else {
    return T(re, im);
  }
}

}  // end of namespace detail

template <class T>
std::vector<T> get_array_of(const toml::value &param, const char *key);

template <class T>
std::vector<T> get_array_of(const toml::value &param, const char *key,
                            T default_value);

template <class T>
std::vector<T> get_array_of(const toml::value &param, const char *key,
                            std::vector<T> const &default_value);

/*! @brief Load a parameter if key is defined.
 *
 * @param dst Destination of the parameter.
 * @param param TOML table.
 * @param key Key of the parameter.
 *
 * @throw input_error if the parameter is found but the type is not matched.
 */
template <class T>
inline void load_if(T &dst, const toml::value &param, const char *key) {
  if (param.contains(key)) {
    dst = detail::as_value<T>(param.at(key));
  }
}

/*! @brief Load a parameter if key is defined.
 *
 * A scalar value is loaded as a list with a single element.
 *
 * @param dst Destination of the parameter.
 * @param param TOML table.
 * @param key Key of the parameter.
 *
 * @throw input_error if the parameter is found but the type is not matched.
 */
template <class T>
inline void load_if(std::vector<T> &dst, const toml::value &param,
                    const char *key) {
  if (param.contains(key)) {
    dst = get_array_of<T>(param, key);
  }
}

/*! @brief
 * @param param TOML table.
 * @param key Key of the parameter.
 * @param value Default value.
 * @return Value of the parameter if the parameter is defined, otherwise
 * default is returned.
 * @throw input_error if the parameter is found but the type is not matched.
 */
template <class T>
inline T find_or(const toml::value &param, const char *key, T value) {
  T ret = value;
  load_if(ret, param, key);
  return ret;
}

/*! @brief Find a parameter from a TOML table.
 *  @param param TOML table.
 *  @param key Key of the parameter.
 *  @return Value of the parameter.
 *  @throw input_error if the parameter is not found.
 */
template <class T>
inline T find(const toml::value &param, const char *key) {
  if (!param.contains(key)) {
    throw input_error(detail::msg_cannot_find(key));
  }
  return detail::as_value<T>(param.at(key));
}

/*! @brief Find a list parameter from a TOML table.
 *
 * If the parameter is scalar, it is converted to a list with single element.
 *
 * @param param TOML table.
 * @param key Key of the parameter.
 * @return Value of the parameter.
 * @throw input_error if the parameter is not found.
 */
template <class T>
std::vector<T> get_array_of(const toml::value &param, const char *key) {
  if (!param.contains(key)) {
    throw input_error(detail::msg_cannot_find(key));
  }
  const toml::value &v = param.at(key);
  if (v.is_array()) {
    return detail::as_array<T>(v);
  }
  return std::vector<T>{detail::as_value<T>(v)};
}

/*! @brief Find a list parameter from a TOML table.
 *
 * If the parameter is scalar, it is converted to a list with single element.
 *
 * @param param TOML table.
 * @param key Key of the parameter.
 * @param default_value Default value.
 * @return Value of the parameter.
 */
template <class T>
std::vector<T> get_array_of(const toml::value &param, const char *key,
                            T default_value) {
  if (!param.contains(key)) {
    return std::vector<T>{default_value};
  }
  return get_array_of<T>(param, key);
}

template <class T>
std::vector<T> get_array_of(const toml::value &param, const char *key,
                            std::vector<T> const &default_value) {
  if (!param.contains(key)) {
    return default_value;
  }
  return get_array_of<T>(param, key);
}

SquareLattice gen_lattice(const toml::value &toml, const char *tablename) {
  if (!toml.contains("L_sub")) {
    throw input_error(detail::msg_cannot_find("L_sub", tablename));
  }
  const auto Lsub = detail::as_array<int64_t>(toml.at("L_sub"));
  if (Lsub.size() < 2) {
    throw input_error(toml::format_error("[error] L_sub should have 2 integers",
                                         toml.at("L_sub"), "given here"));
  }

  auto skew = find_or<int>(toml, "skew", 0);

  SquareLattice lat(Lsub[0], Lsub[1], skew);

  if (!toml.contains("unitcell")) {
    throw input_error(detail::msg_cannot_find("unitcell", tablename));
  }
  const toml::value &sites = toml.at("unitcell");
  if (!sites.is_array()) {
    throw input_error(toml::format_error(
        "[error] unitcell should be an array of tables", sites, "given here"));
  }
  for (const auto &site : sites.as_array()) {
    if (!site.is_table()) {
      throw input_error(
          toml::format_error("[error] an element of unitcell should be a table",
                             site, "given here"));
    }
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

CorrelationParameter gen_corparam(const toml::value &toml,
                                  const char *tablename) {
  int rmax = find<int>(toml, "r_max");

  if (!toml.contains("operators")) {
    throw input_error(detail::msg_cannot_find("operators", tablename));
  }
  const toml::value &oplist = toml.at("operators");
  if (!oplist.is_array()) {
    throw input_error(toml::format_error(
        "[error] operators should be an array of pairs of integers", oplist,
        "given here"));
  }

  std::vector<std::tuple<int, int>> ops;
  for (const auto &op : oplist.as_array()) {
    auto i = detail::as_array<int64_t>(op);
    if (i.size() < 2) {
      throw input_error(toml::format_error(
          "[error] an element of operators should have 2 integers", op,
          "given here"));
    }
    ops.emplace_back(static_cast<int>(i[0]), static_cast<int>(i[1]));
  }

  return CorrelationParameter{rmax, ops};
}

TransferMatrix_Parameters gen_transfer_matrix_parameter(const toml::value &toml,
                                                        const char *tablename) {
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

PEPS_Parameters gen_param(const toml::value &param) {
  PEPS_Parameters pparam;

  // general
  const toml::value *general = detail::opt_table(param, "general");
  if (general != nullptr) {
    load_if(pparam.is_real, *general, "is_real");
    load_if(pparam.iszero_tol, *general, "iszero_tol");
    load_if(pparam.to_measure, *general, "measure");
    load_if(pparam.outdir, *general, "output");
    load_if(pparam.tensor_load_dir, *general, "tensor_load");
    load_if(pparam.tensor_save_dir, *general, "tensor_save");

    std::string mode_str =
        find_or(*general, "mode", std::string("ground state"));
    if (util::startswith(mode_str, "ground")) {
      pparam.calcmode = PEPS_Parameters::CalculationMode::ground_state;
    } else if (util::startswith(mode_str, "time")) {
      pparam.calcmode = PEPS_Parameters::CalculationMode::time_evolution;
    } else if (util::startswith(mode_str, "finite")) {
      pparam.calcmode = PEPS_Parameters::CalculationMode::finite_temperature;
    } else {
      throw input_error("Invalid mode: " + mode_str);
    }
    load_if(pparam.measure_interval, *general, "measure_interval");
  }

  // Simple update
  const toml::value *simple = detail::opt_table(param, "simple_update");
  if (simple != nullptr) {
    load_if(pparam.num_simple_step, *simple, "num_step");
    load_if(pparam.tau_simple_step, *simple, "tau");
    load_if(pparam.Inverse_lambda_cut, *simple, "lambda_cutoff");
    load_if(pparam.Simple_Gauge_Fix, *simple, "gauge_fix");
    load_if(pparam.Simple_Gauge_maxiter, *simple, "gauge_maxiter");
    load_if(pparam.Simple_Gauge_Convergence_Epsilon, *simple,
            "gauge_convergence_epsilon");
  }

  // Full update
  const toml::value *full = detail::opt_table(param, "full_update");
  if (full != nullptr) {
    load_if(pparam.num_full_step, *full, "num_step");
    load_if(pparam.tau_full_step, *full, "tau");
    load_if(pparam.Full_Inverse_precision, *full, "inverse_precision");
    load_if(pparam.Full_Convergence_Epsilon, *full, "convergence_epsilon");
    load_if(pparam.Inverse_Env_cut, *full, "env_cutoff");
    load_if(pparam.Full_max_iteration, *full, "iteration_max");
    load_if(pparam.Full_Gauge_Fix, *full, "gauge_fix");
    load_if(pparam.Full_Use_FastFullUpdate, *full, "fastfullupdate");
  }

  // Environment
  const toml::value *ctm = detail::opt_table(param, "ctm");
  if (ctm != nullptr) {
    load_if(pparam.CHI, *ctm, "dimension");
    load_if(pparam.Inverse_projector_cut, *ctm, "projector_cutoff");
    load_if(pparam.CTM_Convergence_Epsilon, *ctm, "convergence_epsilon");
    load_if(pparam.Max_CTM_Iteration, *ctm, "iteration_max");
    load_if(pparam.CTM_Projector_corner, *ctm, "projector_corner");
    load_if(pparam.Use_RSVD, *ctm, "use_rsvd");
    load_if(pparam.RSVD_Oversampling_factor, *ctm, "rsvd_oversampling_factor");
    load_if(pparam.MeanField_Env, *ctm, "meanfield_env");

    if (pparam.RSVD_Oversampling_factor < 1.0) {
      std::string msg = "rsvd_oversampling_factor must be >= 1.0";
      throw tenes::input_error(msg);
    }
  }

  // random
  const toml::value *random = detail::opt_table(param, "random");
  if (random != nullptr) {
    load_if(pparam.seed, *random, "seed");
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
  for (std::size_t i = 1; i < words.size(); i += 2) {
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
Operators<tensor> load_operator(const toml::value &param, MPI_Comm comm,
                                int nsites, int nbody, double atol,
                                const char *tablename) {
  if (!param.is_table()) {
    throw input_error(toml::format_error(
        std::string("[error] a table is required in an array \"") + tablename +
            "\"",
        param, "given here"));
  }
  const bool has_elements = param.contains("elements");
  const bool has_ops = param.contains("ops");
  if (has_elements && has_ops) {
    std::stringstream ss;
    ss << "Both elements and ops are defined in a section " << tablename;
    throw tenes::input_error(ss.str());
  }
  tensor A(comm);
  std::vector<int> op_ind;
  if (nbody == 1 && !has_elements) {
    std::stringstream ss;
    ss << "elements not found in a section " << tablename;
    throw tenes::input_error(ss.str());
  }
  if (has_elements) {
    if (nbody > 2) {
      std::stringstream ss;
      ss << "observable.multisite does not support elements";
      throw tenes::input_error(ss.str());
    }
    auto elements = detail::as_value<std::string>(param.at("elements"));
    mptensor::Shape shape;
    if (!param.contains("dim")) {
      throw input_error(detail::msg_cannot_find("dim", tablename));
    }
    const toml::value &dim = param.at("dim");
    if (dim.is_array()) {
      auto dim_arr = detail::as_array<int64_t>(dim);
      if (dim_arr.size() != static_cast<std::size_t>(nbody)) {
        std::stringstream ss;
        ss << "operator is " << nbody << "-sites but dim has " << dim_arr.size()
           << " integers";
        throw input_error(ss.str());
      }
      for (int d : dim_arr) {
        shape.push(d);
      }
    } else {
      const int dim_int = detail::as_value<int>(dim);
      for (int i = 0; i < nbody; ++i) {
        shape.push(dim_int);
      }
    }
    for (int i = 0; i < nbody; ++i) {
      shape.push(shape[i]);
    }
    A = util::read_tensor<tensor>(elements, shape, comm, atol);
  } else if (has_ops) {
    auto ops = detail::as_array<int64_t>(param.at("ops"));
    op_ind.assign(ops.begin(), ops.end());
  } else {
    std::stringstream ss;
    ss << "Neither elements nor ops are not defined in a section " << tablename;
    throw tenes::input_error(ss.str());
  }

  auto group = find<int>(param, "group");
  auto name = find_or(param, "name", std::string(""));
  auto coeff_re = find_or(param, "coeff", 1.0);
  auto coeff_im = find_or(param, "coeff_im", 0.0);
  typename tensor::value_type coeff =
      detail::reduce<typename tensor::value_type>(coeff_re, coeff_im);
  auto is_real = std::is_same_v<typename tensor::value_type, double>;
  if (is_real && coeff_im != 0.0) {
    std::stringstream ss;
    ss << "parameter.general.is_real is true but coeff_im is not zero in a "
          "section "
       << tablename;
    throw tenes::input_error(ss.str());
  }

  std::vector<Operator<tensor>> ret;
  if (nbody == 1) {
    if (!param.contains("sites")) {
      throw input_error(detail::msg_cannot_find("sites", tablename));
    }
    const toml::value &sites_val = param.at("sites");
    std::vector<int> sites;
    if (sites_val.is_array()) {
      auto site_arr = detail::as_array<int64_t>(sites_val);
      sites.assign(site_arr.begin(), site_arr.end());
      if (sites.empty()) {
        for (int i = 0; i < nsites; ++i) {
          sites.push_back(i);
        }
      }
    } else {
      sites.push_back(detail::as_value<int>(sites_val));
    }
    for (int s : sites) {
      ret.emplace_back(name, group, s, A, coeff);
    }
  } else if (nbody == 2) {
    auto bonds_str = find<std::string>(param, "bonds");
    auto bonds = read_bonds(bonds_str);
    for (const auto &[source, dx, dy] : bonds) {
      if (has_elements) {
        ret.emplace_back(name, group, source, dx, dy, A, coeff);
      } else {
        ret.emplace_back(name, group, source, dx, dy, op_ind, coeff);
      }
    }
  } else {
    auto ms_str = find<std::string>(param, "multisites");
    auto mss = read_multisitesset(ms_str);
    for (const auto &[source, dx, dy] : mss) {
      if (has_elements) {
        ret.emplace_back(name, group, source, dx, dy, A, coeff);
      } else {
        ret.emplace_back(name, group, source, dx, dy, op_ind, coeff);
      }
    }
  }
  return ret;
}

template <class tensor>
Operators<tensor> load_operators(const toml::value &param, MPI_Comm comm,
                                 int nsites, int nbody, double atol,
                                 std::string const &key) {
  Operators<tensor> ret;
  const toml::value *tables = detail::find_path(param, key);
  if (tables == nullptr) {
    std::cout << "INFO: " << key << " is not found in the input file. (Skipped)"
              << std::endl;
    return ret;
  }
  if (!tables->is_array()) {
    throw input_error(toml::format_error(
        std::string("[error] \"") + key +
            "\" should be an array of tables ([[" + key + "]])",
        *tables, "given here"));
  }
  for (const auto &table : tables->as_array()) {
    auto obs =
        load_operator<tensor>(table, comm, nsites, nbody, atol, key.c_str());
    std::copy(obs.begin(), obs.end(), std::back_inserter(ret));
    // std::move(obs.begin(), obs.end(), std::back_inserter(ret));
  }
  return ret;
}

template <class tensor>
EvolutionOperator<tensor> load_Evolution_operator(const toml::value &param,
                                                  MPI_Comm comm, double atol,
                                                  const char *tablename) {
  if (!param.is_table()) {
    throw input_error(toml::format_error(
        std::string("[error] a table is required in an array \"") + tablename +
            "\"",
        param, "given here"));
  }
  if (!param.contains("dimensions")) {
    throw input_error(detail::msg_cannot_find("dimensions", tablename));
  }
  const auto dimensions = detail::as_array<int64_t>(param.at("dimensions"));
  auto elements = find<std::string>(param, "elements");
  auto shape = mptensor::Shape();
  for (auto d : dimensions) {
    shape.push(d);
  }
  auto group = find<int>(param, "group");

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
EvolutionOperators<tensor> load_updates(const toml::value &param, MPI_Comm comm,
                                        double atol, std::string const &key) {
  EvolutionOperators<tensor> ret;
  const toml::value *tables = detail::find_path(param, key);
  if (tables == nullptr) {
    return ret;
  }
  if (!tables->is_array()) {
    throw input_error(toml::format_error(
        std::string("[error] \"") + key +
            "\" should be an array of tables ([[" + key + "]])",
        *tables, "given here"));
  }
  for (const auto &table : tables->as_array()) {
    ret.push_back(
        load_Evolution_operator<tensor>(table, comm, atol, key.c_str()));
  }
  return ret;
}
template <class tensor>
EvolutionOperators<tensor> load_simple_updates(const toml::value &param,
                                               MPI_Comm comm, double atol) {
  return load_updates<tensor>(param, comm, atol, "evolution.simple");
}
template <class tensor>
EvolutionOperators<tensor> load_full_updates(const toml::value &param,
                                             MPI_Comm comm, double atol) {
  return load_updates<tensor>(param, comm, atol, "evolution.full");
}

// template instantiations

template Operators<real_tensor> load_operator(const toml::value &param,
                                              MPI_Comm comm, int nsites,
                                              int nbody, double atol,
                                              const char *tablename);
template Operators<complex_tensor> load_operator(const toml::value &param,
                                                 MPI_Comm comm, int nsites,
                                                 int nbody, double atol,
                                                 const char *tablename);

template Operators<real_tensor> load_operators(const toml::value &param,
                                               MPI_Comm comm, int nsites,
                                               int nbody, double atol,
                                               std::string const &key);

template Operators<complex_tensor> load_operators(const toml::value &param,
                                                  MPI_Comm comm, int nsites,
                                                  int nbody, double atol,
                                                  std::string const &key);

template EvolutionOperator<real_tensor> load_Evolution_operator(
    const toml::value &param, MPI_Comm comm, double atol,
    const char *tablename);
template EvolutionOperator<complex_tensor> load_Evolution_operator(
    const toml::value &param, MPI_Comm comm, double atol,
    const char *tablename);

template EvolutionOperators<real_tensor> load_updates(const toml::value &param,
                                                      MPI_Comm comm,
                                                      double atol,
                                                      std::string const &key);
template EvolutionOperators<complex_tensor> load_updates(
    const toml::value &param, MPI_Comm comm, double atol,
    std::string const &key);

template EvolutionOperators<real_tensor> load_simple_updates(
    const toml::value &param, MPI_Comm comm, double atol);
template EvolutionOperators<complex_tensor> load_simple_updates(
    const toml::value &param, MPI_Comm comm, double atol);

template EvolutionOperators<real_tensor> load_full_updates(
    const toml::value &param, MPI_Comm comm, double atol);
template EvolutionOperators<complex_tensor> load_full_updates(
    const toml::value &param, MPI_Comm comm, double atol);

}  // namespace tenes::itps
