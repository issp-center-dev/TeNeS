add_subdirectory(core)

add_library(PEPS_Parameters STATIC PEPS_Parameters.cpp)
target_link_libraries(PEPS_Parameters PUBLIC tenes_mpi)

add_library(load_toml STATIC load_toml.cpp)
target_link_libraries(load_toml PUBLIC tenes_mpi tensor SquareLattice PEPS_Parameters)
target_include_directories(load_toml PUBLIC ${CPPTOML_INCLUDE_DIR} ${BOOST_INCLUDE_DIR})

add_library(
  iTPS_impl STATIC
  iTPS.cpp
  tensors.cpp
  measure.cpp
  onesite_obs.cpp
  twosite_obs.cpp
  multisite_obs.cpp
  density.cpp
  simple_update.cpp
  full_update.cpp
  optimize.cpp
  time_evolution.cpp
  finite_temperature.cpp
  contraction_path.cpp
  correlation_function.cpp
  correlation_length.cpp
  transfer_matrix.cpp)

target_link_libraries(iTPS_impl PUBLIC iTPS_core SquareLattice load_toml tenes_util PEPS_Parameters find_contraction_path)
