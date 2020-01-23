#ifndef UTIL_STRING_HPP
#define UTIL_STRING_HPP

#include <string>
#include <vector>

namespace tenes {
namespace util {

std::vector<std::string> split(std::string const &str,
                               std::string const &delim = " \t");

std::string lstrip(std::string const &str, std::string const &delim = " \t\n");

std::string rstrip(std::string const &str, std::string const &delim = " \t\n");

std::string strip(std::string const &str, std::string const &delim = " \t\n");

std::string drop_comment(std::string const &str);

}  // namespace util
}  // namespace tenes

#endif  // UTIL_STRING_HPP
