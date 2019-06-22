include_directories(${MPTENSOR_INCLUDE_DIR})
link_directories(${CMAKE_BINARY_DIR}/deps/mptensor/src)

add_executable(tnsolve
main.cpp
tnsolve.cpp 
PEPS_Parameters.cpp
Parameters.cpp
edge.cpp
hamiltonian.cpp
util/string.cpp
)

target_include_directories(tnsolve PUBLIC ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(tnsolve mptensor)
target_link_libraries(tnsolve ${MPI_CXX_LIBRARIES})
install(TARGETS tnsolve RUNTIME DESTINATION bin)
