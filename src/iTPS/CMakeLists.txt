add_subdirectory(core)

add_library(
  iTPS STATIC
  iTPS.cpp
  tensors.cpp
  measure.cpp
  onesite_obs.cpp
  twosite_obs.cpp
  simple_update.cpp
  full_update.cpp
  optimize.cpp
  correlation_function.cpp
  correlation_length.cpp)

target_compile_options(iTPS PRIVATE ${OMP_FLAG})
target_compile_options(iTPS PRIVATE -Wall)
target_compile_options(
  iTPS PRIVATE $<$<CONFIG:Debug>: -Wextra -Wno-unused-parameter
               -Wno-sign-compare -Wpedantic >)

target_include_directories(iTPS PRIVATE ${MPTENSOR_INCLUDE_DIR})
target_include_directories(iTPS PRIVATE ${MPI_CXX_INCLUDE_DIRS})
target_include_directories(iTPS PRIVATE ${OpenMP_CXX_INCLUDE_DIRS})
target_include_directories(iTPS PRIVATE ${CPPTOML_INCLUDE_DIR})
target_include_directories(iTPS PRIVATE ${BUNDLED_DEPS_DIR})

target_link_libraries(iTPS INTERFACE iTPS_core)
target_link_libraries(iTPS INTERFACE tensor)
# target_link_libraries(iTPS INTERFACE mptensor)
# target_link_libraries(iTPS INTERFACE tenes_mpi)
# target_link_libraries(
#   iTPS INTERFACE ${MPI_CXX_LIBRARIES} ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES}
#   ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES})

if(NOT OPENMP_FOUND)
  target_compile_definitions(iTPS PRIVATE -D_NO_OMP)
endif()

if(NOT ENABLE_MPI)
  target_compile_definitions(iTPS PRIVATE -D_NO_MPI)
endif()
