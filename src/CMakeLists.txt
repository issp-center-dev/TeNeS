add_subdirectory(util)

link_directories(${MPTENSOR_LIB_DIR})
if(OpenMP_CXX_LIBRARY_DIRS)
  link_directories(${OpenMP_CXX_LIBRARY_DIRS})
endif()

add_library(tenes_mpi STATIC mpi.cpp)

target_compile_options(tenes_mpi PUBLIC -Wall)
target_compile_options(
  tenes_mpi PUBLIC $<$<CONFIG:Debug>: -Wextra -Wno-unused-parameter
                   -Wno-sign-compare -Wpedantic >)

if(OPENMP_FOUND)
  target_compile_options(tenes_mpi PUBLIC ${OMP_FLAG})
else()
  target_compile_definitions(tenes_mpi PUBLIC -D_NO_OMP)
endif()

if(NOT ENABLE_MPI)
  target_compile_definitions(tenes_mpi PUBLIC -D_NO_MPI)
endif()

target_include_directories(tenes_mpi PUBLIC ${MPI_CXX_INCLUDE_DIRS} ${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(tenes_mpi INTERFACE ${MPI_CXX_LIBRARIES}
                                          ${OpenMP_CXX_LIBRARIES})

add_library(tensor STATIC tensor.cpp arnoldi.cpp)
target_include_directories(tensor PUBLIC ${MPTENSOR_INCLUDE_DIR})
target_link_libraries(tensor INTERFACE ${SCALAPACK_LIBRARIES}
                                       ${LAPACK_LIBRARIES} mptensor)
target_link_libraries(tensor PUBLIC tenes_mpi)

add_library(SquareLattice STATIC SquareLattice.cpp)
target_link_libraries(SquareLattice PUBLIC tenes_mpi)

add_subdirectory(iTPS)

add_executable(tenes main.cpp)
target_link_libraries(tenes PRIVATE iTPS)

install(TARGETS tenes RUNTIME DESTINATION bin)
