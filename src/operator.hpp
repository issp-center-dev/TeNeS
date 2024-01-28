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

#ifndef TENES_SRC_OPERATOR_HPP_
#define TENES_SRC_OPERATOR_HPP_

#include <vector>
#include <string>
#include <iostream>

namespace tenes {

template <class tensor>
struct Operator {
  using value_type = typename tensor::value_type;

  std::string name;
  int group;
  int source_site;
  std::vector<int> dx;
  std::vector<int> dy;
  tensor op;
  std::vector<int> ops_indices;
  typename tensor::value_type coeff;

  /*!
   * @brief Constructor for a one-site operator.
   *
   * @param[in] name Name of the operator.
   * @param[in] group Group of the operator.
   * @param[in] site Site of the operator.
   * @param[in] op Operator tensor.
   * @param[in] coeff coefficient of operator (default: 1.0).
   */
  Operator(std::string const &name, int group, int site, tensor const &op, value_type coeff=static_cast<value_type>(1.0))
      : name(name), group(group), source_site(site), dx(0), dy(0), op(op), coeff(coeff) {
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
           tensor const &op, value_type coeff=static_cast<value_type>(1.0))
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
           std::vector<int> const &ops_indices, value_type coeff=static_cast<value_type>(1.0))
      : name(name),
        group(group),
        source_site(source_site),
        dx(1, dx),
        dy(1, dy),
        ops_indices(ops_indices),
        coeff(coeff){
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
           tensor const &op, value_type coeff=static_cast<value_type>(1.0))
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
    if (op.rank() != 2 * (dx.size()+1)) {
      throw std::runtime_error("Operator tensor must be rank 2 * (dx.size()+1).");
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
           std::vector<int> const &ops_indices, value_type coeff=static_cast<value_type>(1.0))
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
    if (ops_indices.size() != dx.size()+1) {
      throw std::runtime_error(
          "Operator must be a product of dx.size()+1 one-site operators.");
    }
  }

  bool is_onesite() const { return dx.empty(); }
  int nsites() const { return dx.size() + 1; }
};

template <class tensor>
using Operators = std::vector<Operator<tensor>>;

template <class tensor>
struct EvolutionOperator {
  int source_site;
  int source_leg;
  int group;
  tensor op;

  EvolutionOperator(int source_site, int source_leg, int group,
                    tensor const &op)
      : source_site(source_site),
        source_leg(source_leg),
        group(group),
        op(op) {}

  bool is_onesite() const { return source_leg < 0; }
  bool is_twosite() const { return !is_onesite(); }
  bool is_horizontal() const { return source_leg % 2 == 0; }
  bool is_vertical() const { return !is_horizontal(); }
};
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

template <class tensor>
using EvolutionOperators = std::vector<EvolutionOperator<tensor>>;

}  // namespace tenes

#endif  // TENES_SRC_OPERATOR_HPP_
