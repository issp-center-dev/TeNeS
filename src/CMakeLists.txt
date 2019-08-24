include_directories(${MPTENSOR_INCLUDE_DIR})
link_directories(${MPTENSOR_LIB_DIR})

add_executable(tnsolve
main.cpp
tnsolve.cpp 
PEPS_Parameters.cpp
Lattice.cpp
edge.cpp
util/string.cpp
)
add_sanitizers(tnsolve)

target_include_directories(tnsolve PUBLIC ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(tnsolve mptensor)
target_link_libraries(tnsolve ${MPI_CXX_LIBRARIES} ${SCALAPACK_LIBS})
install(TARGETS tnsolve RUNTIME DESTINATION bin)
