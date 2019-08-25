# include_directories(${MPTENSOR_INCLUDE_DIR})
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

# target_include_directories(tnsolve PUBLIC ${MPI_CXX_INCLUDE_DIRS})
# target_link_libraries(tnsolve mptensor)
# target_link_libraries(tnsolve ${MPI_CXX_LIBRARIES} ${SCALAPACK_LIBS})

target_compile_options(tnsolve PRIVATE ${OMP_FLAG})
target_compile_options(tnsolve PRIVATE -Wall)
target_compile_options(tnsolve PRIVATE $<$<CONFIG:Debug>: -Wextra>)

target_include_directories(tnsolve PRIVATE ${MPTENSOR_INCLUDE_DIR})
target_include_directories(tnsolve PRIVATE ${MPI_CXX_INCLUDE_DIRS})
target_include_directories(tnsolve PRIVATE ${PROJECT_SOURCE_DIR}/deps/cpptoml/include)

target_link_libraries(tnsolve mptensor)
target_link_libraries(tnsolve ${MPI_CXX_LIBRARIES} ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES})

install(TARGETS tnsolve RUNTIME DESTINATION bin)
