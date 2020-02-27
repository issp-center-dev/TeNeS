# include_directories(${MPTENSOR_INCLUDE_DIR})
link_directories(${MPTENSOR_LIB_DIR})

add_executable(tenes
main.cpp
main_impl.cpp
tenes.cpp 
PEPS_Parameters.cpp
Lattice.cpp
util/string.cpp
mpi.cpp
)

if (USE_SANITIZER)
    add_sanitizers(tenes)
endif()

# target_include_directories(tenes PUBLIC ${MPI_CXX_INCLUDE_DIRS})
# target_link_libraries(tenes mptensor)
# target_link_libraries(tenes ${MPI_CXX_LIBRARIES} ${SCALAPACK_LIBS})

if (NOT ENABLE_MPI)
    target_compile_definitions(tenes PRIVATE -D_NO_MPI)
endif()

target_compile_options(tenes PRIVATE ${OMP_FLAG})
target_compile_options(tenes PRIVATE -Wall)
target_compile_options(tenes PRIVATE $<$<CONFIG:Debug>: -Wextra -Wno-unused-parameter -Wno-sign-compare -Wpedantic >)

target_include_directories(tenes PRIVATE ${MPTENSOR_INCLUDE_DIR})
target_include_directories(tenes PRIVATE ${MPI_CXX_INCLUDE_DIRS})
target_include_directories(tenes PRIVATE ${DEPS_DIR}/cpptoml/include)

target_link_libraries(tenes mptensor)
target_link_libraries(tenes ${MPI_CXX_LIBRARIES} ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES})

install(TARGETS tenes RUNTIME DESTINATION bin)
