link_directories(${MPTENSOR_LIB_DIR})
if(OpenMP_CXX_LIBRARY_DIRS)
  link_directories(${OpenMP_CXX_LIBRARY_DIRS})
endif()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(tenes_unit_tests input simple_update full_update tensor_util util
    correlation_length saveload arnoldi timer_registry onesite_rdm)
if(TENES_ARPACK_FOUND)
  list(APPEND tenes_unit_tests arpack)
endif()
foreach(basename IN LISTS tenes_unit_tests)
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

# A syntactically broken input file must lead to a clean error exit,
# not to an abort by an uncaught exception (WILL_FAIL passes on a
# nonzero exit code but still fails if the process dies on a signal).
if(ENABLE_MPI)
  add_test(NAME broken_input
           COMMAND ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 1
                   ${MPIEXEC_POSTFLAGS} $<TARGET_FILE:tenes>
                   data/broken.toml)
else()
  add_test(NAME broken_input COMMAND $<TARGET_FILE:tenes> data/broken.toml)
endif()
set_tests_properties(broken_input PROPERTIES WILL_FAIL TRUE)

add_test(NAME python_unittest
         COMMAND ${TENES_PYTHON_EXECUTABLE} -m pytest
                 ${CMAKE_CURRENT_SOURCE_DIR}/python)

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

foreach(name quasi1d_ctm onesite_rdm_output)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${name}.py.in
                 ${CMAKE_CURRENT_BINARY_DIR}/${name}.py @ONLY)
  add_test(NAME ${name} COMMAND ${TENES_PYTHON_EXECUTABLE}
                                ${CMAKE_CURRENT_BINARY_DIR}/${name}.py)
endforeach()

if(ENABLE_MPI)
  set(BENCH_SMOKE_LAUNCHER
      "${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_POSTFLAGS}")
else()
  set(BENCH_SMOKE_LAUNCHER "")
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/benchmark_smoke.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/benchmark_smoke.py @ONLY)
add_test(NAME benchmark_smoke
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/benchmark_smoke.py)
