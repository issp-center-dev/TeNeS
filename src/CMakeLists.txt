add_subdirectory(util)

link_directories(${MPTENSOR_LIB_DIR})
if(OpenMP_CXX_LIBRARY_DIRS)
  link_directories(${OpenMP_CXX_LIBRARY_DIRS})
endif()

# version.cpp is rewritten at every build, since the commit changes without
# a new configuration; it is touched only when its content changes
add_custom_target(tenes_version_source
  COMMAND ${CMAKE_COMMAND}
    -DSOURCE_DIR=${TeNeS_SOURCE_DIR}
    -DTENES_VERSION=${TENES_VERSION}
    -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/version.cpp.in
    -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/version.cpp
    -P ${TeNeS_SOURCE_DIR}/config/write_version.cmake
  BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
  VERBATIM)

add_library(tenes_mpi STATIC mpi.cpp timer.cpp
                             ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)
target_include_directories(tenes_mpi PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_dependencies(tenes_mpi tenes_version_source)

target_compile_options(tenes_mpi PUBLIC -Wall)
target_compile_options(
  tenes_mpi PUBLIC $<$<CONFIG:Debug>: -Wextra -Wno-unused-parameter
                   -Wno-sign-compare -Wpedantic >)

target_compile_options(tenes_mpi PUBLIC ${TENES_COMPILE_OPTIONS})

target_include_directories(tenes_mpi PUBLIC ${MPI_CXX_INCLUDE_DIRS}
                                            ${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(tenes_mpi INTERFACE ${MPI_CXX_LIBRARIES}
                                          ${OpenMP_CXX_LIBRARIES})

add_library(tensor STATIC tensor.cpp arnoldi.cpp arpack_solver.cpp)
target_include_directories(tensor PUBLIC ${MPTENSOR_INCLUDE_DIR})
target_link_libraries(tensor INTERFACE ${SCALAPACK_LIBRARIES}
                                       ${LAPACK_LIBRARIES}
                                       ${TENES_ARPACK_LIBRARIES} mptensor)
target_link_libraries(tensor PUBLIC tenes_mpi)

add_library(SquareLattice STATIC SquareLattice.cpp)
target_link_libraries(SquareLattice PUBLIC tenes_mpi)

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
