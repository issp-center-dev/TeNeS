link_directories(${MPTENSOR_LIB_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fulltest.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fulltest.py
               @ONLY)

foreach(basename input simple_update)
    set(testname "test_${basename}")
    add_executable(${testname} "${basename}.cpp")

    target_compile_options(${testname} PRIVATE ${OMP_FLAG})

    target_include_directories(${testname} PRIVATE ${MPTENSOR_INCLUDE_DIR})
    target_include_directories(${testname} PRIVATE ${MPI_CXX_INCLUDE_DIRS})
    target_include_directories(${testname} PRIVATE ${PROJECT_SOURCE_DIR}/src)
    target_include_directories(${testname} PRIVATE ${PROJECT_SOURCE_DIR}/deps/cpptoml/include)

    target_link_libraries(${testname} mptensor)
    target_link_libraries(${testname} ${MPI_CXX_LIBRARIES} ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES})

    add_test(NAME ${testname} COMMAND $<TARGET_FILE:${testname}>)
endforeach()

add_test(NAME fulltest COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/fulltest.py)
