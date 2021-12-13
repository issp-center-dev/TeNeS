add_library(
  tenes_util STATIC
  ./file.cpp
  ./string.cpp
  ./datetime.cpp
  )

target_compile_options(tenes_util PRIVATE -Wall)
target_compile_options(
  tenes_util PRIVATE $<$<CONFIG:Debug>: -Wextra -Wno-unused-parameter
               -Wno-sign-compare -Wpedantic >)
