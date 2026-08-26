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
 *  (observables and evolution gates). validate_fermion_constraints() is
 *  the input gate of fermion mode: everything the fermion layer under
 *  src/fermion/ assumes about its inputs is enforced here, at load time.
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

/*!
 * @brief Read the per-site physical-leg parities for fermion mode.
 *
 * Reads the optional "parity" array of each tensor.unitcell entry: one 0
 * (even) or 1 (odd) per physical index, whose length must equal the
 * site's physical dimension. Sites without a "parity" key get an empty
 * vector — validate_fermion_constraints() rejects that later if fermion
 * mode is on, so bosonic inputs stay valid without the key.
 *
 * @return One parity vector per site of the unit cell.
 */
std::vector<std::vector<bool>> gen_phys_parity(
    const toml::value &toml, const SquareLattice &lattice,
    const char *tablename = "tensor");

/*!
 * @brief Input gate of fermion mode: reject everything the fermion layer
 *        does not support.
 *
 * No-op unless parameter.general.fermion is true. Otherwise throws
 * tenes::input_error for inputs the fermionic code paths cannot handle,
 * so the library code under src/fermion/ can assume its preconditions
 * instead of re-checking them. Rejected are, in this version:
 *
 *   - missing or wrong-length tensor.unitcell.parity metadata,
 *   - unit cells smaller than 2x2 or with skew != 0 (skewed cells were
 *     measured to give wrong fermionic numbers),
 *   - any calculation mode other than ground state, the full update,
 *     Simple_Gauge_Fix, Use_RSVD, and correlation.r_max > 0
 *     (long-range correlation functions),
 *   - multi-site operators, two-site operators in the ops form or acting
 *     beyond nearest neighbors,
 *   - product initial states with weight on an odd-parity basis state,
 *   - parity-odd operators and gates (any tensor element connecting the
 *     even and odd sectors, checked elementwise across all processes).
 */
template <class tensor>
void validate_fermion_constraints(
    const PEPS_Parameters &peps_parameters, const SquareLattice &lattice,
    const EvolutionOperators<tensor> &simple_updates,
    const EvolutionOperators<tensor> &full_updates,
    const Operators<tensor> &onesite_operators,
    const Operators<tensor> &twosite_operators,
    const Operators<tensor> &multisite_operators,
    const CorrelationParameter &corparam);

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

/*! @brief Add identity two-site gates on every bond no gate acts on.
 *
 *  A bond with virtual dimension > 1 that carries no evolution operator is
 *  never touched by the simple update: it keeps its random, non-canonical
 *  initialisation, and the CTM does not converge on the resulting state.
 *  An identity gate is physically a no-op but makes the update truncate and
 *  canonicalise that bond every step.  One gate is added per ungated bond and
 *  per group present among the two-site operators (group 0 if there is none).
 *
 *  @param simple_updates the gates read from the input.
 *  @param lattice the unit-cell geometry listing the bonds to cover.
 *  @param comm the communicator the generated gate tensors are allocated on.
 *              It must be the one the other operators and the state live on;
 *              it is not inferred from @p simple_updates, which may contain
 *              no two-site operator at all.
 *  @return the input operators followed by the added identity gates.
 */
template <class tensor>
EvolutionOperators<tensor> complete_ungated_bonds(
    const EvolutionOperators<tensor> &simple_updates,
    const SquareLattice &lattice, MPI_Comm comm);

}  // namespace tenes::itps

#endif  // TENES_SRC_ITPS_LOAD_TOML_HPP_
