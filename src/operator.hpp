#ifndef TENES_OPERATOR_HPP
#define TENES_OPERATOR_HPP

#include <vector>
#include <string>

namespace tenes {

template <class tensor> struct Operator {
  std::string name;
  int group;
  int source_site;
  std::vector<int> dx;
  std::vector<int> dy;
  tensor op;
  std::vector<int> ops_indices;

  // onesite
  Operator(std::string const& name, int group, int site, tensor const &op): name(name), group(group), source_site(site), dx(0), dy(0), op(op){}

  // twosite
  Operator(std::string const& name, int group, int source_site, int dx, int dy, tensor const &op): name(name), group(group), source_site(source_site), dx(1, dx), dy(1, dy), op(op) {}
  Operator(std::string const& name, int group, int source_site, std::vector<int> const &dx, std::vector<int> const &dy, tensor const &op): name(name), group(group), source_site(source_site), dx(dx), dy(dy), op(op) {}
  Operator(std::string const& name, int group, int source_site, int dx, int dy, std::vector<int> const &ops_indices)
      : name(name), group(group), source_site(source_site), dx(1, dx), dy(1, dy), ops_indices(ops_indices) {}
  Operator(std::string const& name, int group, int source_site, std::vector<int> const &dx, std::vector<int> const &dy, std::vector<int> const &ops_indices)
      : name(name), group(group), source_site(source_site), dx(dx), dy(dy), ops_indices(ops_indices) {}

  bool is_onesite() const {return dx.empty();}
};

template <class tensor> using Operators = std::vector<Operator<tensor>>;

template <class tensor> struct NNOperator {
  int source_site;
  int source_leg;
  tensor op;

  NNOperator(int site, int leg, tensor const &op)
      : source_site(site), source_leg(leg), op(op) {}

  bool is_horizontal() const { return source_leg % 2 == 0; }
  bool is_vertical() const { return !is_horizontal(); }
};

template <class tensor> using NNOperators = std::vector<NNOperator<tensor>>;

} // namespace tenes

#endif // TENES_OPERATOR_HPP
