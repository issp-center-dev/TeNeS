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
 *  @brief Operators the solver works with: measured observables
 *         (::tenes::Operator) and Trotter evolution gates
 *         (::tenes::EvolutionOperator).
 */

#ifndef TENES_SRC_OPERATOR_HPP_
#define TENES_SRC_OPERATOR_HPP_

#include <vector>
#include <string>
#include <iostream>

namespace tenes {

/*! @brief One term of a measured observable.
 *
 *  Acts on source_site plus one further site per entry of dx/dy (empty
 *  for a one-site operator). The matrix elements are given either as one
 *  tensor (op) or as a product of one-site operators referenced by index
 *  (ops_indices); terms sharing a group are summed into one reported
 *  observable.
 */
template <class tensor>
struct Operator {
  //! Scalar type of the tensor elements.
  using value_type = typename tensor::value_type;

  std::string name;     //!< display name of the observable
  int group;            //!< observable this term belongs to
  int source_site;      //!< site index within the unit cell
  std::vector<int> dx;  //!< x displacement of each further site
  std::vector<int> dy;  //!< y displacement of each further site
  tensor op;            //!< matrix elements (empty if ops_indices is used)
  //! Indices of the one-site operators whose product forms this term
  //! (empty if op is used).
  std::vector<int> ops_indices;
  typename tensor::value_type coeff;  //!< overall coefficient of the term

  /*!
   * @brief Constructor for a one-site operator.
   *
   * @param[in] name Name of the operator.
   * @param[in] group Group of the operator.
   * @param[in] site Site of the operator.
   * @param[in] op Operator tensor.
   * @param[in] coeff coefficient of operator (default: 1.0).
   */
  Operator(std::string const &name, int group, int site, tensor const &op,
           value_type coeff = static_cast<value_type>(1.0))
      : name(name),
        group(group),
        source_site(site),
        dx(0),
        dy(0),
        op(op),
        coeff(coeff) {
    if (op.rank() != 2) {
      throw std::runtime_error("Operator tensor must be rank 2.");
    }
  }

  /*!
   * @brief Constructor for a two-site operator.
   *
   * @param[in] name Name of the operator.
   * @param[in] group Group of the operator.
   * @param[in] source_site Index of a site.
   * @param[in] dx X displacement of the other site.
   * @param[in] dy Y displacement of the other site.
   * @param[in] op Operator tensor.
   * @param[in] coeff coefficient of operator (default: 1.0).
   */
  Operator(std::string const &name, int group, int source_site, int dx, int dy,
           tensor const &op, value_type coeff = static_cast<value_type>(1.0))
      : name(name),
        group(group),
        source_site(source_site),
        dx(1, dx),
        dy(1, dy),
        op(op),
        coeff(coeff) {
    if (op.rank() != 4) {
      throw std::runtime_error("Operator tensor must be rank 4.");
    }
  }

  /*!
   * @brief Constructor for a two-site operator represented by the product of
   * two one-site operators.
   *
   * @param[in] name Name of the operator.
   * @param[in] group Group of the operator.
   * @param[in] source_site Index of a site.
   * @param[in] dx X displacement of the other site.
   * @param[in] dy Y displacement of the other site.
   * @param[in] ops_indices Onesite operator indices.
   * @param[in] coeff coefficient of operator (default: 1.0).
   */
  Operator(std::string const &name, int group, int source_site, int dx, int dy,
           std::vector<int> const &ops_indices,
           value_type coeff = static_cast<value_type>(1.0))
      : name(name),
        group(group),
        source_site(source_site),
        dx(1, dx),
        dy(1, dy),
        ops_indices(ops_indices),
        coeff(coeff) {
    if (ops_indices.size() != 2) {
      throw std::runtime_error(
          "Operator must be a product of two one-site operators.");
    }
  }

  /*!
   * @brief Constructor for a multi-site operator.
   *
   * @param[in] name Name of the operator.
   * @param[in] group Group of the operator.
   * @param[in] source_site Index of a site.
   * @param[in] dx X displacement of the other sites.
   * @param[in] dy Y displacement of the other sites.
   * @param[in] op Operator tensor.
   * @param[in] coeff coefficient of operator (default: 1.0).
   */
  Operator(std::string const &name, int group, int source_site,
           std::vector<int> const &dx, std::vector<int> const &dy,
           tensor const &op, value_type coeff = static_cast<value_type>(1.0))
      : name(name),
        group(group),
        source_site(source_site),
        dx(dx),
        dy(dy),
        op(op),
        coeff(coeff) {
    if (dx.size() != dy.size()) {
      throw std::runtime_error("dx and dy must have the same size.");
    }
    if (op.rank() != 2 * (dx.size() + 1)) {
      throw std::runtime_error(
          "Operator tensor must be rank 2 * (dx.size()+1).");
    }
  }

  /*!
   * @brief Constructor for a multi-site operator represented by the product of
   * one-site operators.
   *
   * @param[in] name Name of the operator.
   * @param[in] group Group of the operator.
   * @param[in] source_site Index of a site.
   * @param[in] dx X displacement of the other sites.
   * @param[in] dy Y displacement of the other sites.
   * @param[in] ops_indices Onesite operator indices.
   * @param[in] coeff coefficient of operator (default: 1.0).
   */
  Operator(std::string const &name, int group, int source_site,
           std::vector<int> const &dx, std::vector<int> const &dy,
           std::vector<int> const &ops_indices,
           value_type coeff = static_cast<value_type>(1.0))
      : name(name),
        group(group),
        source_site(source_site),
        dx(dx),
        dy(dy),
        ops_indices(ops_indices),
        coeff(coeff) {
    if (dx.size() != dy.size()) {
      throw std::runtime_error("dx and dy must have the same size.");
    }
    if (ops_indices.size() != dx.size() + 1) {
      throw std::runtime_error(
          "Operator must be a product of dx.size()+1 one-site operators.");
    }
  }

  //! True iff the term acts on a single site.
  bool is_onesite() const { return dx.empty(); }
  //! Number of sites the term acts on.
  int nsites() const { return dx.size() + 1; }
};

//! All terms of all observables of one arity (one-, two-, or multi-site).
template <class tensor>
using Operators = std::vector<Operator<tensor>>;

/*! @brief One Trotter gate of the (imaginary or real) time evolution.
 *
 *  A two-site gate acts on the bond leaving source_site through
 *  source_leg (0=left, 1=top, 2=right, 3=bottom); a one-site gate has
 *  source_leg = -1. Gates sharing a group form one Trotter step and are
 *  driven with that group's tau and num_step parameters.
 */
template <class tensor>
struct EvolutionOperator {
  int source_site;  //!< site index within the unit cell
  int source_leg;   //!< bond direction from source_site; -1 for one-site
  int group;        //!< evolution-operator group this gate belongs to
  tensor op;        //!< gate elements: (in, out) or (in1, in2, out1, out2)

  //! Prefer the make_*_EvolutionOperator() factories, which validate.
  EvolutionOperator(int source_site, int source_leg, int group,
                    tensor const &op)
      : source_site(source_site),
        source_leg(source_leg),
        group(group),
        op(op) {}

  //! True iff the gate acts on a single site.
  bool is_onesite() const { return source_leg < 0; }
  //! True iff the gate acts on a bond.
  bool is_twosite() const { return !is_onesite(); }
  //! True iff the bond is horizontal (left/right).
  bool is_horizontal() const { return source_leg % 2 == 0; }
  //! True iff the bond is vertical (top/bottom).
  bool is_vertical() const { return !is_horizontal(); }
};
//! Build a one-site gate, validating the arguments.
template <class tensor>
EvolutionOperator<tensor> make_onesite_EvolutionOperator(int source_site,
                                                         int group,
                                                         tensor const &op) {
  if (source_site < 0) {
    throw std::runtime_error("source_site must be non-negative");
  }
  if (group < 0) {
    throw std::runtime_error("group must be non-negative");
  }
  return EvolutionOperator<tensor>(source_site, -1, group, op);
}
//! Build a two-site (bond) gate, validating the arguments.
template <class tensor>
EvolutionOperator<tensor> make_twosite_EvolutionOperator(int source_site,
                                                         int source_leg,
                                                         int group,
                                                         tensor const &op) {
  if (source_site < 0) {
    throw std::runtime_error("source_site must be non-negative");
  }
  if (source_leg < 0 || source_leg > 3) {
    throw std::runtime_error("source_leg must be 0, 1, 2, or 3");
  }
  if (group < 0) {
    throw std::runtime_error("group must be non-negative");
  }
  return EvolutionOperator<tensor>(source_site, source_leg, group, op);
}

//! All gates of one update kind (simple or full), all groups interleaved.
template <class tensor>
using EvolutionOperators = std::vector<EvolutionOperator<tensor>>;

}  // namespace tenes

#endif  // TENES_SRC_OPERATOR_HPP_
