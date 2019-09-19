# include_directories(${MPTENSOR_INCLUDE_DIR})
link_directories(${MPTENSOR_LIB_DIR})

add_executable(tenes
main.cpp
main_impl.cpp
tenes.cpp 
PEPS_Parameters.cpp
Lattice.cpp
edge.cpp
util/string.cpp
)
add_sanitizers(tenes)

# target_include_directories(tenes PUBLIC ${MPI_CXX_INCLUDE_DIRS})
# target_link_libraries(tenes mptensor)
# target_link_libraries(tenes ${MPI_CXX_LIBRARIES} ${SCALAPACK_LIBS})

target_compile_options(tenes PRIVATE ${OMP_FLAG})
target_compile_options(tenes PRIVATE -Wall)
target_compile_options(tenes PRIVATE $<$<CONFIG:Debug>: -Wextra>)

target_include_directories(tenes PRIVATE ${MPTENSOR_INCLUDE_DIR})
target_include_directories(tenes PRIVATE ${MPI_CXX_INCLUDE_DIRS})
target_include_directories(tenes PRIVATE ${PROJECT_SOURCE_DIR}/deps/cpptoml/include)
target_include_directories(tenes PRIVATE ${PROJECT_SOURCE_DIR}/deps/CLI11)

target_link_libraries(tenes mptensor)
target_link_libraries(tenes ${MPI_CXX_LIBRARIES} ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES})

install(TARGETS tenes RUNTIME DESTINATION bin)
