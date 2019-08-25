link_directories(${MPTENSOR_LIB_DIR})

foreach(basename simple_update)
    set(testname "test_${basename}")
    add_executable(${testname} "${basename}.cpp" "${PROJECT_SOURCE_DIR}/src/PEPS_Parameters.cpp")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/${basename}.dat
        ${CMAKE_CURRENT_BINARY_DIR}/data/${basename}.dat
        COPYONLY
        )

    target_include_directories(${testname} PRIVATE ${MPTENSOR_INCLUDE_DIR})
    target_include_directories(${testname} PRIVATE ${MPI_CXX_INCLUDE_DIRS})
    target_include_directories(${testname} PRIVATE ${PROJECT_SOURCE_DIR}/src)
    target_include_directories(${testname} PRIVATE ${PROJECT_SOURCE_DIR}/deps/cpptoml/include)

    target_link_libraries(${testname} mptensor)
    target_link_libraries(${testname} ${MPI_CXX_LIBRARIES} ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES})

    add_test(NAME ${testname} COMMAND $<TARGET_FILE:${testname}>)
endforeach()
