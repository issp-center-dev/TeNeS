#ifndef UTIL_FIND_OR_HPP
#define UTIL_FIND_OR_HPP

#include <toml11/toml.hpp>

namespace util{

template <class T>
T find_or(const toml::table& tab, const toml::key key, T opt)
{
  if(tab.count(key) == 0){
    return opt;
  }else{
    return toml::get<T>(tab.at(key));
  }
}

} // end of namespace util

#endif // UTIL_FIND_OR_HPP
