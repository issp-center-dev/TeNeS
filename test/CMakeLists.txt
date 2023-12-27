link_directories(${MPTENSOR_LIB_DIR})
if(OpenMP_CXX_LIBRARY_DIRS)
  link_directories(${OpenMP_CXX_LIBRARY_DIRS})
endif()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

foreach(basename input simple_update full_update)
  set(testname "test_${basename}")
  add_executable(${testname} "${basename}.cpp")

  target_link_libraries(${testname} PUBLIC tenes_static)

  if(ENABLE_MPI)
    add_test(NAME ${testname}
             COMMAND ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 1
                     ${MPIEXEC_POSTFLAGS} $<TARGET_FILE:${testname}>)
  else()
    add_test(NAME ${testname} COMMAND $<TARGET_FILE:${testname}>)
  endif()
endforeach()

foreach(name simple_mode std_mode)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${name}.py.in
                 ${CMAKE_CURRENT_BINARY_DIR}/${name}.py @ONLY)
  add_test(NAME ${name} COMMAND ${TENES_PYTHON_EXECUTABLE}
                                ${CMAKE_CURRENT_BINARY_DIR}/${name}.py)
endforeach()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fulltest.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fulltest.py @ONLY)

add_test(NAME restart COMMAND ${TENES_PYTHON_EXECUTABLE}
                              ${CMAKE_CURRENT_BINARY_DIR}/restart.py)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/restart.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/restart.py @ONLY)

foreach(
  name
  AntiferroHeisenberg_real
  AntiferroHeisenberg_complex
  AntiferroHeisenberg_mf
  RSVD
  J1J2_AFH
  Honeycomb
  Honeycomb_skew
  TE_TFI
  FT_TFI_square
  FT_Kitaev)
  add_test(NAME ${name} COMMAND ${TENES_PYTHON_EXECUTABLE}
                                ${CMAKE_CURRENT_BINARY_DIR}/fulltest.py ${name})
endforeach()
