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

#ifndef TENES_SRC_ITPS_LOAD_TOML_HPP_
#define TENES_SRC_ITPS_LOAD_TOML_HPP_

#define _USE_MATH_DEFINES
#include <tuple>
#include <vector>
#include <string>

#include <limits>
#include "cpptoml.h"  // IWYU pragma: export

// IWYU pragma: begin_exports
#include "../SquareLattice.hpp"
#include "../operator.hpp"
#include "../mpi.hpp"

#include "transfer_matrix.hpp"
#include "correlation_function.hpp"
#include "PEPS_Parameters.hpp"
#include "contraction_path.hpp"

// IWYU pragma: end_exports

namespace tenes {
namespace itps {

SquareLattice gen_lattice(decltype(cpptoml::parse_file("")) toml,
                          const char *tablename = "tensor");

CorrelationParameter gen_corparam(decltype(cpptoml::parse_file("")) toml,
                                  const char *tablename = "correlation");

TransferMatrix_Parameters gen_transfer_matrix_parameter(
    decltype(cpptoml::parse_file("")) toml,
    const char *tablename = "correlation_length");

PEPS_Parameters gen_param(decltype(cpptoml::parse_file("")) param);

std::tuple<int, int, int> read_bond(std::string line);
std::vector<std::tuple<int, int, int>> read_bonds(std::string str);

template <class tensor>
Operators<tensor> load_operator(decltype(cpptoml::parse_file("")) param,
                                MPI_Comm comm, int nsites, int nbody,
                                double atol = 0.0,
                                const char *tablename = "observable.onesite");

template <class tensor>
Operators<tensor> load_operators(decltype(cpptoml::parse_file("")) param,
                                 MPI_Comm comm, int nsites, int nbody,
                                 double atol, std::string const &key);

template <class tensor>
EvolutionOperator<tensor> load_Evolution_operator(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, double atol = 0.0,
    const char *tablename = "evolution.simple");

template <class tensor>
EvolutionOperators<tensor> load_updates(decltype(cpptoml::parse_file("")) param,
                                        MPI_Comm comm, double atol,
                                        std::string const &key);
template <class tensor>
EvolutionOperators<tensor> load_simple_updates(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, double atol = 0.0);
template <class tensor>
EvolutionOperators<tensor> load_full_updates(
    decltype(cpptoml::parse_file("")) param, MPI_Comm comm, double atol = 0.0);

template <class tensor>
std::map<TNC_map_key, TensorNetworkContractor<tensor>> load_contraction_paths(
    std::string const &path_file);

}  // namespace itps
}  // namespace tenes

#endif  // TENES_SRC_ITPS_LOAD_TOML_HPP_
