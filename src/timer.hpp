#ifndef TENES_TIME_HPP
#define TENES_TIME_HPP

#include <chrono>

template <class CT = std::chrono::high_resolution_clock>
class Timer{
public:
  using clock_type = CT;

  Timer():start(clock_type::now()){}
  void reset(){start = clock_type::now();}
  double elapsed() const
  {
    const auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(clock_type::now()-start);
    return elapsed_time.count()*1.0e-9;
  }
private:
  typename clock_type::time_point start;
};

#endif // TENES_TIME_HPP
