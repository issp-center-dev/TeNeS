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

#ifndef TENES_SRC_LOAD_TOML_HPP_
#define TENES_SRC_LOAD_TOML_HPP_

#define _USE_MATH_DEFINES
#include <cpptoml.h>
#include <tuple>
#include <vector>
#include <string>

#include "Lattice.hpp"
#include "PEPS_Parameters.hpp"
#include "correlation.hpp"
#include "correlation_length.hpp"
#include "operator.hpp"

namespace tenes {

Lattice gen_lattice(decltype(cpptoml::parse_file("")) toml,
                    const char *tablename = "tensor");

CorrelationParameter gen_corparam(decltype(cpptoml::parse_file("")) toml,
                                  const char *tablename = "correlation");

CorrelationLengthCalculator_Parameters gen_correlationlength_parameter(
    decltype(cpptoml::parse_file("")) toml,
    const char *tablename = "correlation_length");

PEPS_Parameters gen_param(decltype(cpptoml::parse_file("")) param);

std::tuple<int, int, int> read_bond(std::string line);
std::vector<std::tuple<int, int, int>> read_bonds(std::string str);

template <class tensor>
Operators<tensor> load_operator(decltype(cpptoml::parse_file("")) param,
                                int nsites, int nbody, double atol = 0.0,
                                const char *tablename = "observable.onesite");

template <class tensor>
Operators<tensor> load_operators(decltype(cpptoml::parse_file("")) param,
                                 int nsites, int nbody, double atol,
                                 std::string const &key);

template <class tensor>
NNOperator<tensor> load_nn_operator(decltype(cpptoml::parse_file("")) param,
                                    double atol = 0.0,
                                    const char *tablename = "evolution.simple");

template <class tensor>
NNOperators<tensor> load_updates(decltype(cpptoml::parse_file("")) param,
                                 double atol, std::string const &key);
template <class tensor>
NNOperators<tensor> load_simple_updates(decltype(cpptoml::parse_file("")) param,
                                        double atol = 0.0);
template <class tensor>
NNOperators<tensor> load_full_updates(decltype(cpptoml::parse_file("")) param,
                                      double atol = 0.0);

}  // end of namespace tenes

#endif  // TENES_SRC_LOAD_TOML_HPP_
