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

target_compile_options(tenes_mpi PUBLIC ${TENES_COMPILE_OPTIONS})

target_include_directories(tenes_mpi PUBLIC ${MPI_CXX_INCLUDE_DIRS}
                                            ${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(tenes_mpi INTERFACE ${MPI_CXX_LIBRARIES}
                                          ${OpenMP_CXX_LIBRARIES})

add_library(tensor STATIC tensor.cpp arnoldi.cpp)
target_include_directories(tensor PUBLIC ${MPTENSOR_INCLUDE_DIR})
target_link_libraries(tensor INTERFACE ${SCALAPACK_LIBRARIES}
                                       ${LAPACK_LIBRARIES} mptensor)
target_link_libraries(tensor PUBLIC tenes_mpi)

add_library(SquareLattice STATIC SquareLattice.cpp)
target_link_libraries(SquareLattice PUBLIC tenes_mpi)

add_subdirectory(find_contraction_path)

add_subdirectory(iTPS)

add_library(tenes_static STATIC iTPS/main.cpp initialize_mptensor.cpp)
target_link_libraries(tenes_static PUBLIC iTPS_impl)

add_library(libtenes SHARED iTPS/main.cpp initialize_mptensor.cpp)
set_target_properties(libtenes PROPERTIES LIBRARY_OUTPUT_NAME tenes)
target_link_libraries(libtenes PUBLIC iTPS_impl)

add_executable(tenes main.cpp)
target_link_libraries(tenes PUBLIC tenes_static)

install(TARGETS libtenes LIBRARY DESTINATION lib)
install(TARGETS tenes RUNTIME DESTINATION bin)
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION include/tenes
  FILES_MATCHING
  PATTERN "*.hpp")
