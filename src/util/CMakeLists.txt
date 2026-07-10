add_library(
  tenes_util STATIC
  ./file.cpp
  ./string.cpp
  ./datetime.cpp
  )

# GCC 8 provides std::filesystem in a separate library
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
  target_link_libraries(tenes_util PUBLIC stdc++fs)
endif()

target_compile_options(tenes_util PRIVATE -Wall)
target_compile_options(
  tenes_util PRIVATE $<$<CONFIG:Debug>: -Wextra -Wno-unused-parameter
               -Wno-sign-compare -Wpedantic >)
