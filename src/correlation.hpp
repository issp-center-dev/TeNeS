#ifndef CORRELATION_HPP
#define CORRELATION_HPP

#include <vector>
#include <utility>

struct Correlation {
  int left_index, right_index;
  int offset_x, offset_y;
  int left_op, right_op;
  double real, imag;
};

struct CorrelationParameter{
  int r_max;
  std::vector<std::pair<int,int>> operators;
};

#endif // CORRELATION_HPP
