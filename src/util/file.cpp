#include <fstream>
#include <string>
#include <sstream>
#include <sys/stat.h>

#include "../exception.hpp"
#include "string.hpp"

#include "file.hpp"

namespace tenes {
namespace util{

bool path_exists(const std::string& path) {
  struct stat status;
  return stat(path.c_str(), &status) == 0;
}

bool isdir(const std::string& path) {
  struct stat status;
  int code = stat(path.c_str(), &status);
  if(code != 0) return false;
  return S_ISDIR(status.st_mode);
}

bool mkdir(const std::string& path){
  int ret = 0;
  struct stat st;
  if(stat(path.c_str(), &st) != 0){
    ret = ::mkdir(path.c_str(), 0755);
  }
  return ret == 0;
}

std::string joinpath(std::vector<std::string> const& xs){
  const size_t n = xs.size();
  if(n==0){ return std::string(""); }
  std::stringstream ss;
  ss << xs[0];
  for(size_t i=1; i<n; ++i){
    ss << "/";
    ss << xs[i];
  }
  return ss.str();
}

std::string basename(const std::string& path){
  auto xs = util::split(path, "/");
  return xs[xs.size()-1];
}

}  // end of namespace util
}  // end of namespace tenes
