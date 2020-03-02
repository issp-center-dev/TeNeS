link_directories(${MPTENSOR_LIB_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

foreach(basename input simple_update full_update)
    set(testname "test_${basename}")
    add_executable(${testname} "${basename}.cpp")

    target_compile_options(${testname} PRIVATE ${OMP_FLAG})
    if (NOT ENABLE_MPI)
        target_compile_definitions(${testname} PRIVATE -D_NO_MPI)
    endif()

    target_include_directories(${testname} PRIVATE ${MPTENSOR_INCLUDE_DIR})
    target_include_directories(${testname} PRIVATE ${MPI_CXX_INCLUDE_DIRS})
    target_include_directories(${testname} PRIVATE ${PROJECT_SOURCE_DIR}/src)
    target_include_directories(${testname} PRIVATE ${CPPTOML_INCLUDE_DIR})

    target_link_libraries(${testname} mptensor)
    target_link_libraries(${testname} ${MPI_CXX_LIBRARIES} ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${OpenMP_CXX_LIBRARIES})

    if(ENABLE_MPI)
        add_test(NAME ${testname} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${testname}>)
    else()
        add_test(NAME ${testname} COMMAND $<TARGET_FILE:${testname}>)
    endif()
endforeach()

foreach(name simple_mode std_mode)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${name}.py.in
                   ${CMAKE_CURRENT_BINARY_DIR}/${name}.py
                   @ONLY)
    add_test(NAME ${name} COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/${name}.py)
endforeach()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fulltest.py.in ${CMAKE_CURRENT_BINARY_DIR}/fulltest.py @ONLY)

foreach(name AntiferroHeisenberg_real AntiferroHeisenberg_complex J1J2_AFH)
    add_test(NAME ${name} COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/fulltest.py ${name})
endforeach()

