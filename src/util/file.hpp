#ifndef UTIL_FILE_HPP
#define UTIL_FILE_HPP

#include <vector>
#include <string>

namespace tenes {
namespace util{

bool path_exists(const std::string& path);
bool isdir(const std::string& path);
bool mkdir(const std::string& path);
std::string joinpath(std::vector<std::string> const& xs);
std::string basename(const std::string& path);

}  // end of namespace util
}  // end of namespace tenes

#endif  // UTIL_FILE_HPP
