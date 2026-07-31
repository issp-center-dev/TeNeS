add_library(
  tenes_util STATIC
  ./file.cpp
  ./string.cpp
  ./datetime.cpp
  )

# GCC 8 (libstdc++) and LLVM libc++ < 9 provide std::filesystem in a separate
# library (stdc++fs / c++fs). Checking the compiler id is not enough: icpx and
# clang pick up the system libstdc++, so probe whether a std::filesystem
# program actually links and add the extra library only when needed.
include(CheckCXXSourceCompiles)
set(TENES_FILESYSTEM_TEST_SOURCE [[
#include <filesystem>
int main() {
  return std::filesystem::exists(std::filesystem::path("/")) ? 0 : 1;
}
]])
check_cxx_source_compiles("${TENES_FILESYSTEM_TEST_SOURCE}" TENES_FILESYSTEM_LINKS)
if(NOT TENES_FILESYSTEM_LINKS)
  set(CMAKE_REQUIRED_LIBRARIES stdc++fs)
  check_cxx_source_compiles("${TENES_FILESYSTEM_TEST_SOURCE}" TENES_FILESYSTEM_LINKS_STDCXXFS)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(TENES_FILESYSTEM_LINKS_STDCXXFS)
    target_link_libraries(tenes_util PUBLIC stdc++fs)
  else()
    set(CMAKE_REQUIRED_LIBRARIES c++fs)
    check_cxx_source_compiles("${TENES_FILESYSTEM_TEST_SOURCE}" TENES_FILESYSTEM_LINKS_CXXFS)
    unset(CMAKE_REQUIRED_LIBRARIES)
    if(TENES_FILESYSTEM_LINKS_CXXFS)
      target_link_libraries(tenes_util PUBLIC c++fs)
    else()
      message(FATAL_ERROR
        "std::filesystem could not be compiled or linked, even with stdc++fs/c++fs. "
        "See CMakeFiles/CMakeError.log for details.")
    endif()
  endif()
endif()

target_compile_options(tenes_util PRIVATE -Wall)
target_compile_options(
  tenes_util PRIVATE $<$<CONFIG:Debug>: -Wextra -Wno-unused-parameter
               -Wno-sign-compare -Wpedantic >)
