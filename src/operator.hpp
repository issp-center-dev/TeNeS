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

namespace tenes {

template <class tensor>
struct Operator {
  std::string name;
  int group;
  int source_site;
  std::vector<int> dx;
  std::vector<int> dy;
  tensor op;
  std::vector<int> ops_indices;

  // onesite
  Operator(std::string const &name, int group, int site, tensor const &op)
      : name(name), group(group), source_site(site), dx(0), dy(0), op(op) {}

  // twosite
  Operator(std::string const &name, int group, int source_site, int dx, int dy,
           tensor const &op)
      : name(name),
        group(group),
        source_site(source_site),
        dx(1, dx),
        dy(1, dy),
        op(op) {}
  Operator(std::string const &name, int group, int source_site,
           std::vector<int> const &dx, std::vector<int> const &dy,
           tensor const &op)
      : name(name),
        group(group),
        source_site(source_site),
        dx(dx),
        dy(dy),
        op(op) {}
  Operator(std::string const &name, int group, int source_site, int dx, int dy,
           std::vector<int> const &ops_indices)
      : name(name),
        group(group),
        source_site(source_site),
        dx(1, dx),
        dy(1, dy),
        ops_indices(ops_indices) {}
  Operator(std::string const &name, int group, int source_site,
           std::vector<int> const &dx, std::vector<int> const &dy,
           std::vector<int> const &ops_indices)
      : name(name),
        group(group),
        source_site(source_site),
        dx(dx),
        dy(dy),
        ops_indices(ops_indices) {}

  bool is_onesite() const { return dx.empty(); }
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
