#ifndef TENES_OPERATOR_HPP
#define TENES_OPERATOR_HPP

#include <vector>

namespace tenes {

template <class tensor> struct Operator {
  int source_site;
  int source_leg;
  tensor op;

  Operator(int site_, int leg_, tensor const &op_)
      : source_site(site_), source_leg(leg_), op(op_) {}

  bool is_horizontal() const { return source_leg % 2 == 0; }
  bool is_vertical() const { return !is_horizontal(); }
};

template <class tensor> using Operators = std::vector<Operator<tensor>>;

} // namespace tenes

#endif // TENES_OPERATOR_HPP
