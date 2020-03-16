#ifndef UTIL_STRING_HPP
#define UTIL_STRING_HPP

#include "string.hpp"

#include <string>
#include <vector>

namespace tenes {
namespace util {

std::vector<std::string> split(std::string const &str,
                               std::string const &delim = " \t") {
  std::vector<std::string> words;
  auto index = 0;
  auto last = str.size();
  while (index != std::string::npos) {
    auto next = str.find_first_of(delim, index);
    if(next == std::string::npos){
      words.push_back(str.substr(index, last - index + 1));
      break;
    }
    words.push_back(str.substr(index, next - index + 1));
    index = str.find_first_not_of(delim, next);
  }
  return words;
}

std::string lstrip(std::string const &str, std::string const &delim = " \t\n") {
  std::string ret("");
  auto index = str.find_first_not_of(delim);
  if (index != std::string::npos) {
    ret = str.substr(index, str.size() - index + 1);
  }
  return ret;
}

std::string rstrip(std::string const &str, std::string const &delim = " \t\n") {
  std::string ret("");
  auto index = str.find_last_not_of(delim);
  if (index != std::string::npos) {
    ret = str.substr(0, index + 1);
  }
  return ret;
}

std::string strip(std::string const &str, std::string const &delim = " \t\n") {
  return rstrip(lstrip(str, delim), delim);
}

std::string drop_comment(std::string const &str) {
  auto index = str.find_first_of("#");
  return str.substr(0, index);
}

}  // namespace util
}  // namespace tenes

#endif  // UTIL_STRING_HPP
