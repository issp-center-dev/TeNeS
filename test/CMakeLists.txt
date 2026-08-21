link_directories(${MPTENSOR_LIB_DIR})
if(OpenMP_CXX_LIBRARY_DIRS)
  link_directories(${OpenMP_CXX_LIBRARY_DIRS})
endif()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(tenes_unit_tests input simple_update full_update tensor_util util
    correlation_length saveload arnoldi timer_registry)
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

add_executable(test_fermion_layer test_fermion_layer.cpp)
target_link_libraries(test_fermion_layer PUBLIC tenes_static)
if(ENABLE_MPI)
  add_test(NAME test_fermion_layer
           COMMAND ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 1
                   ${MPIEXEC_POSTFLAGS} $<TARGET_FILE:test_fermion_layer>)
else()
  add_test(NAME test_fermion_layer COMMAND $<TARGET_FILE:test_fermion_layer>)
endif()

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
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fermion)
# The full free-fermion matrix (D = 4 included) takes many hours because the
# two-site reduced-pair measurement scales as D^12; the default quick mode
# runs D = 2 only. TENES_FREE_FERMION_FULL=1 at runtime also switches modes.
option(TENES_FERMION_FULL_TEST
       "Run the free-fermion end-to-end test with the full (slow) D matrix"
       OFF)
if(TENES_FERMION_FULL_TEST)
  set(TENES_FERMION_FULL_TEST_PY "1")
else()
  set(TENES_FERMION_FULL_TEST_PY "0")
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fermion/free_fermion.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion.py @ONLY)

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

add_test(NAME FreeFermion
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion.py)
if(TENES_FERMION_FULL_TEST)
  set_tests_properties(FreeFermion PROPERTIES TIMEOUT 86400)
else()
  set_tests_properties(FreeFermion PROPERTIES TIMEOUT 3600)
endif()

# End-to-end test of the fermionic mean-field environment
# ([parameter.ctm] meanfield_env = true) at half filling, compared against the
# CTM result at D = 2 and the exact free-fermion energy at D = 4.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fermion/free_fermion_mf.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_mf.py @ONLY)
add_test(NAME FreeFermionMF
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_mf.py)
set_tests_properties(FreeFermionMF PROPERTIES TIMEOUT 1800)

# End-to-end test of the full three-stage pipeline (tenes_simple -> tenes_std
# -> tenes) for the fermionic model, checked against the handwritten
# expert-mode sample sample/07_spinless_fermion.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fermion/free_fermion_simple.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_simple.py @ONLY)
add_test(NAME FreeFermionSimple
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_simple.py)
set_tests_properties(FreeFermionSimple PROPERTIES TIMEOUT 600)

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
