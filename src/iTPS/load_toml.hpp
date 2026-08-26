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

/*! @file
 *  @brief Parsing input.toml into the objects the solver is built from.
 *
 *  The gen_* functions read one section each (lattice, parameters,
 *  correlation settings); the load_* functions read the operator tables
 *  (observables and evolution gates).
 */

#ifndef TENES_SRC_ITPS_LOAD_TOML_HPP_
#define TENES_SRC_ITPS_LOAD_TOML_HPP_

#define _USE_MATH_DEFINES
#include <tuple>
#include <vector>
#include <string>

#include <limits>
#include <toml.hpp>  // IWYU pragma: export

// IWYU pragma: begin_exports
#include "../SquareLattice.hpp"
#include "../operator.hpp"
#include "../mpi.hpp"

#include "transfer_matrix.hpp"
#include "correlation_function.hpp"
#include "PEPS_Parameters.hpp"
// IWYU pragma: end_exports

namespace tenes::itps {

//! Read the unit-cell geometry from the given table (default: [tensor]).
SquareLattice gen_lattice(const toml::value &toml,
                          const char *tablename = "tensor");

//! Read the correlation-function settings (default: [correlation]).
CorrelationParameter gen_corparam(const toml::value &toml,
                                  const char *tablename = "correlation");

//! Read the correlation-length / transfer-matrix settings
//! (default: [correlation_length]).
TransferMatrix_Parameters gen_transfer_matrix_parameter(
    const toml::value &toml, const char *tablename = "correlation_length");

//! Read the [parameter] table into a PEPS_Parameters.
PEPS_Parameters gen_param(const toml::value &param);

//! Parse "source_site dx dy" (one bond) from one line.
std::tuple<int, int, int> read_bond(std::string line);
//! Parse a multi-line string of bonds with read_bond().
std::vector<std::tuple<int, int, int>> read_bonds(std::string str);

//! Read one observable table (an array of nbody-site operators).
template <class tensor>
Operators<tensor> load_operator(const toml::value &param, MPI_Comm comm,
                                int nsites, int nbody, double atol = 0.0,
                                const char *tablename = "observable.onesite");

//! Read the observable table named by key via load_operator().
template <class tensor>
Operators<tensor> load_operators(const toml::value &param, MPI_Comm comm,
                                 int nsites, int nbody, double atol,
                                 std::string const &key);

//! Read one evolution (Trotter gate) table entry.
template <class tensor>
EvolutionOperator<tensor> load_Evolution_operator(
    const toml::value &param, MPI_Comm comm, double atol = 0.0,
    const char *tablename = "evolution.simple");

//! Read the evolution-operator table named by key.
template <class tensor>
EvolutionOperators<tensor> load_updates(const toml::value &param, MPI_Comm comm,
                                        double atol, std::string const &key);
//! Read the simple-update gates (evolution.simple).
template <class tensor>
EvolutionOperators<tensor> load_simple_updates(const toml::value &param,
                                               MPI_Comm comm,
                                               double atol = 0.0);
//! Read the full-update gates (evolution.full).
template <class tensor>
EvolutionOperators<tensor> load_full_updates(const toml::value &param,
                                             MPI_Comm comm, double atol = 0.0);

}  // namespace tenes::itps

#endif  // TENES_SRC_ITPS_LOAD_TOML_HPP_
