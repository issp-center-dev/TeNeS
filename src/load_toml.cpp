#define _USE_MATH_DEFINES
#include <random>
#include <sys/stat.h>
#include <tuple>
#include <utility>
#include <vector>

#include <cpptoml.h>

#include "Lattice.hpp"
#include "PEPS_Parameters.hpp"
#include "correlation.hpp"
#include "operator.hpp"
#include "tenes.hpp"
#include "util/read_matrix.hpp"
#include "util/string.hpp"

namespace tenes {

template <typename T>
inline T find_or(decltype(cpptoml::parse_file("")) param, const char *key,
                 T value) {
  return param->get_as<T>(key).value_or(value);
}

Lattice gen_lattice(decltype(cpptoml::parse_file("")) toml,
                    const char *tablename = "tensor") {
  auto Lsub = toml->get_array_of<int64_t>("L_sub");
  if (!Lsub) {
    std::cerr << "cannot find L_sub in the section [" << tablename << "]"
              << std::endl;
    assert(false);
  }

  Lattice lat((*Lsub)[0], (*Lsub)[1]);

  auto sites = toml->get_table_array("unitcell");
  if (!sites) {
    std::cerr << "cannot find unitcell in the section [" << tablename << "]"
              << std::endl;
    assert(false);
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
      std::cerr << "cannot find \"index\" in the section [[" << tablename
                << ".site]]" << std::endl;
      assert(false);
    }

    for (int index : indices) {
      auto dir = site->get_array_of<double>("initial_state");
      if (dir) {
        lat.initial_dirs[index] = *dir;
      } else {
        lat.initial_dirs[index] = std::vector<double>{0.0};
      }

      auto noise = site->get_as<double>("noise");
      if (noise) {
        lat.noises[index] = *noise;
      } else {
        lat.noises[index] = 0.0;
      }

      auto pdim = site->get_as<int>("physical_dim");
      lat.physical_dims[index] = *pdim;

      auto vdim_int = site->get_as<int>("virtual_dim");
      auto vdim_arr = site->get_array_of<int64_t>("virtual_dim");
      std::vector<int> vdim;

      if (vdim_arr) {
        vdim.assign((*vdim_arr).begin(), (*vdim_arr).end());
      } else if (vdim_int) {
        vdim.assign(4, *vdim_int);
      }
      for (int i = 0; i < 4; ++i) {
        lat.virtual_dims[index][i] = vdim[i];
      }
    } // end of for indices
  }

  return lat;
}

template <typename tensor>
std::vector<tensor> gen_matrices(decltype(cpptoml::parse_file("")) toml,
                                 const char *key, const char *tablename) {
  auto strs = toml->get_array_of<std::string>(key);
  if (!strs) {
    std::cerr << "cannot find " << key << " in the section [" << tablename
              << "]" << std::endl;
    assert(false);
  }
  return util::read_matrix<tensor>(*strs);
}

CorrelationParameter gen_corparam(decltype(cpptoml::parse_file("")) toml,
                                  const char *tablename = "correlation") {
  auto r = toml->get_as<int64_t>("r_max");
  if (!r) {
    std::cerr << "cannot find r_max in the section [" << tablename << "]"
              << std::endl;
    // ERROR
  }
  int rmax = static_cast<int>(*r);

  auto oplist = toml->get_array_of<cpptoml::array>("operators");
  if (!oplist) {
    std::cerr << "cannot find operators in the section [" << tablename << "]"
              << std::endl;
    // ERROR
  }

  std::vector<std::pair<int, int>> ops;
  for (auto op : *oplist) {
    auto i = op->get_array_of<int64_t>();
    ops.push_back(
        std::make_pair(static_cast<int>((*i)[0]), static_cast<int>((*i)[1])));
  }

  return CorrelationParameter{rmax, ops};
}

PEPS_Parameters gen_param(decltype(cpptoml::parse_file("")) param) {
  PEPS_Parameters pparam;

  // Tensor
  auto tensor = param->get_table("tensor");
  if (tensor != nullptr) {
    pparam.tensor_load_dir = find_or(tensor, "load_dir", std::string(""));
    pparam.tensor_save_dir = find_or(tensor, "save_dir", std::string(""));
  }

  // Simple update
  auto simple = param->get_table("simple_update");
  if (simple != nullptr) {
    pparam.num_simple_step = find_or(simple, "num_step", 0);
    pparam.Inverse_lambda_cut = find_or(simple, "lambda_cutoff", 1e-12);
  }

  // Full update
  auto full = param->get_table("full_update");
  if (full != nullptr) {
    pparam.num_full_step = find_or(full, "num_step", 0);
    pparam.Full_Inverse_precision = find_or(full, "inverse_precision", 1e-12);
    pparam.Full_Convergence_Epsilon =
        find_or(full, "convergence_epsilon", 1e-12);
    pparam.Inverse_Env_cut = find_or(full, "env_cutoff", 1e-12);
    pparam.Full_max_iteration = find_or(full, "iteration_max", 1000);
    pparam.Full_Gauge_Fix = find_or(full, "gauge_fix", true);
    pparam.Full_Use_FastFullUpdate = find_or(full, "fastfullupdate", true);
  }

  // Environment
  auto ctm = param->get_table("ctm");
  if (ctm != nullptr) {
    pparam.CHI = find_or(ctm, "dimension", 4);
    pparam.Inverse_projector_cut = find_or(ctm, "projector_cutoff", 1e-12);
    pparam.CTM_Convergence_Epsilon = find_or(ctm, "convergence_epsilon", 1e-10);
    pparam.Max_CTM_Iteration = find_or(ctm, "iteration_max", 100);
    pparam.CTM_Projector_corner = find_or(ctm, "projector_corner", false);
    pparam.Use_RSVD = find_or(ctm, "use_rsvd", false);
    pparam.RSVD_Oversampling_factor =
        find_or(ctm, "rsvd_oversampling_factor", 2.0);
    if (pparam.RSVD_Oversampling_factor < 1.0) {
      std::string msg = "rsvd_oversampling_factor must be >= 1.0";
      throw std::runtime_error(msg);
    }
  }

  // random
  auto random = param->get_table("random");
  if (random != nullptr) {
    pparam.seed = find_or(random, "seed", 11);
  }

  return pparam;
}

std::tuple<int, int, int, int> read_bond(std::string line) {
  using std::stoi;
  auto words = util::split(util::strip(line));
  assert(words.size() == 4);
  return std::make_tuple(stoi(words[0]), stoi(words[1]), stoi(words[2]),
                         stoi(words[3]));
}
std::vector<std::tuple<int, int, int, int>> read_bonds(std::string str) {
  std::vector<std::tuple<int, int, int, int>> ret;
  std::string line;
  std::stringstream ss(str);
  while (std::getline(ss, line)) {
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
                                int nsites, int nbody) {
  assert(nbody == 1 || nbody == 2);

  mptensor::Shape shape;
  auto dim_arr = param->get_array_of<int64_t>("dim");
  auto dim_int = param->get_as<int>("dim");
  if (dim_arr) {
    assert(dim_arr->size() == nbody);
    for (int d : *dim_arr) {
      shape.push(d);
    }
  } else if (dim_int) {
    for (int i = 0; i < nbody; ++i) {
      shape.push(*dim_int);
    }
  } else{
    assert(false);
  }

  for (int i = 0; i < nbody; ++i) {
    shape.push(shape[i]);
  }

  auto elements = param->get_as<std::string>("elements");
  assert(elements);
  tensor A = util::read_tensor<tensor>(*elements, shape);

  auto group = param->get_as<int>("group");
  assert(group);

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
      std::cerr << "cannot find \"sites\"" << std::endl;
      assert(false);
    }
    for (int s : sites) {
      ret.emplace_back(*group, s, A);
    }
  } else {
    auto bonds_str = param->get_as<std::string>("bonds");
    assert(bonds_str);
    auto bonds = read_bonds(*bonds_str);
    for (auto bond : bonds) {
      ret.emplace_back(*group, std::get<0>(bond), std::get<1>(bond),
                       std::get<2>(bond), std::get<3>(bond), A);
    }
  }
  return ret;
}

template <class tensor>
Operators<tensor> load_operators(decltype(cpptoml::parse_file("")) param,
                                 int nsites, int nbody,
                                 std::string const &key) {
  Operators<tensor> ret;
  auto tables = param->get_table_array_qualified(key);
  for (const auto &table : *tables) {
    auto obs = load_operator<tensor>(table, nsites, nbody);
    std::copy(obs.begin(), obs.end(), std::back_inserter(ret));
    // std::move(obs.begin(), obs.end(), std::back_inserter(ret));
  }
  return ret;
}

template <class tensor>
NNOperator<tensor> load_nn_operator(decltype(cpptoml::parse_file("")) param) {
  auto source_site = param->get_as<int>("source_site");
  auto source_leg = param->get_as<int>("source_leg");
  auto dimensions = param->get_array_of<int64_t>("dimensions");
  auto elements = param->get_as<std::string>("elements");
  auto shape = mptensor::Shape();
  for (auto d : *dimensions) {
    shape.push(d);
  }
  tensor A = util::read_tensor<tensor>(*elements, shape);
  return NNOperator<tensor>(*source_site, *source_leg, A);
}

template <class tensor>
NNOperators<tensor> load_updates(decltype(cpptoml::parse_file("")) param,
                                 std::string const &key) {
  NNOperators<tensor> ret;
  auto tables = param->get_table_array_qualified(key);
  for (const auto &table : *tables) {
    ret.push_back(load_nn_operator<tensor>(table));
  }
  return ret;
}
template <class tensor>
NNOperators<tensor>
load_simple_updates(decltype(cpptoml::parse_file("")) param) {
  return load_updates<tensor>(param, "evolution.simple");
}
template <class tensor>
NNOperators<tensor> load_full_updates(decltype(cpptoml::parse_file("")) param) {
  return load_updates<tensor>(param, "evolution.full");
}

} // end of namespace tenes
