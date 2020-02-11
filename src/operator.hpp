#ifndef TENES_OPERATOR_HPP
#define TENES_OPERATOR_HPP

#include <vector>

namespace tenes {

template <class tensor> struct Operator {
  int group;
  int source_site;
  int target_site;
  int offset_x;
  int offset_y;
  tensor op;

  Operator(int group, int site, tensor const &op)
      : group(group), source_site(site), target_site(-1), offset_x(0),
        offset_y(0), op(op) {}
  Operator(int group, int source_site, int target_site, int offset_x,
           int offset_y, tensor const &op)
      : group(group), source_site(source_site), target_site(target_site),
        offset_x(offset_x), offset_y(offset_y), op(op) {}

  bool is_onsite() const { return target_site < 0; }
  bool is_twobody() const { return !is_onsite(); }
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
