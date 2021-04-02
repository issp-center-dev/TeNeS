add_subdirectory(core)

add_library(load_toml STATIC load_toml.cpp)
target_link_libraries(load_toml PUBLIC tenes_mpi tensor)
target_include_directories(load_toml PUBLIC ${CPPTOML_INCLUDE_DIR})

add_library(
  iTPS STATIC
  main.cpp
  iTPS.cpp
  PEPS_Parameters.cpp
  tensors.cpp
  measure.cpp
  onesite_obs.cpp
  twosite_obs.cpp
  simple_update.cpp
  full_update.cpp
  optimize.cpp
  correlation_function.cpp
  correlation_length.cpp
  transfer_matrix.cpp)

target_link_libraries(iTPS PUBLIC iTPS_core SquareLattice load_toml tenes_util)
