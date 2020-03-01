#ifndef TENES_EXCEPTION_HPP
#define TENES_EXCEPTION_HPP

#include <stdexcept>

namespace tenes{

class logic_error : public std::logic_error {
public:
  logic_error(const std::string &what_arg) : std::logic_error(what_arg) {}
  logic_error(const char *what_arg) : std::logic_error(what_arg) {}
};
class unimplemented_error : public tenes::logic_error {
public:
  unimplemented_error(const std::string &what_arg) : tenes::logic_error(what_arg) {}
  unimplemented_error(const char *what_arg) : tenes::logic_error(what_arg) {}
};

class runtime_error : public std::runtime_error {
public:
  runtime_error(const std::string &what_arg) : std::runtime_error(what_arg) {}
  runtime_error(const char *what_arg) : std::runtime_error(what_arg) {}
};
class input_error : public runtime_error {
public:
  input_error(const std::string &what_arg) : runtime_error(what_arg) {}
  input_error(const char *what_arg) : runtime_error(what_arg) {}
};

} // end of namespace tenes

#endif // TENES_EXCEPTION_HPP
