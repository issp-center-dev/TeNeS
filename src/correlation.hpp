#ifndef CORRELATION_HPP
#define CORRELATION_HPP

#include <utility>
#include <vector>

namespace tenes {

struct Correlation {
  int left_index, right_index;
  int offset_x, offset_y;
  int left_op, right_op;
  double real, imag;
};

struct CorrelationParameter {
  int r_max;
  std::vector<std::pair<int, int>> operators;
  CorrelationParameter() : r_max(0) {}
  CorrelationParameter(int r_max, std::vector<std::pair<int, int>> const& ops)
      : r_max(r_max), operators(ops) {}
};

}  // end of namespace tenes

#endif  // CORRELATION_HPP
