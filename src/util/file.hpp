#ifndef UTIL_FILE_HPP
#define UTIL_FILE_HPP

#include <string>
#include <fstream>

inline bool file_exists(const std::string& path)
{
  std::ifstream ifs(path);
  return ifs.is_open();
}

#endif // UTIL_FILE_HPP
