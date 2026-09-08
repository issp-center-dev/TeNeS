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
  target_compile_definitions(
    ${testname} PRIVATE TENES_TEST_WORKDIR="${CMAKE_CURRENT_BINARY_DIR}")

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
target_compile_definitions(
  test_fermion_layer PRIVATE TENES_TEST_WORKDIR="${CMAKE_CURRENT_BINARY_DIR}")
if(ENABLE_MPI)
  add_test(NAME test_fermion_layer
           COMMAND ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 1
                   ${MPIEXEC_POSTFLAGS} $<TARGET_FILE:test_fermion_layer>)
else()
  add_test(NAME test_fermion_layer COMMAND $<TARGET_FILE:test_fermion_layer>)
endif()

# The finiteness predicate of the fermion diagnostics, compiled with
# -ffast-math on purpose. Under finite-math-only the compiler may assume no
# operand is Inf or NaN, and a std::isfinite-based check folds to a constant
# true - which is how a real oneAPI run reported "all elements finite" about
# a matrix that had already diverged. Its own executable because the flag
# must not reach any other test. Header-only, so it links nothing.
add_executable(test_fermion_finite_math fermion/finite_under_fast_math.cpp)
target_compile_options(test_fermion_finite_math PRIVATE -ffast-math)
add_test(NAME test_fermion_finite_math
         COMMAND $<TARGET_FILE:test_fermion_finite_math>)

# S-5 of
# docs/superpowers/specs/2026-09-07-fermion-fast-full-update-contract.md: the
# warm start of iTPS::update_CTM(). Its own executable rather than another
# file included into test_fermion_layer.cpp, because it is the only test that
# calls the new signature: a translation unit that does not compile would take
# every other fermion unit test with it.
add_executable(test_fermion_ctm_warm_start fermion/ctm_warm_start.cpp)
target_link_libraries(test_fermion_ctm_warm_start PUBLIC tenes_static)
target_compile_definitions(
  test_fermion_ctm_warm_start
  PRIVATE TENES_TEST_WORKDIR="${CMAKE_CURRENT_BINARY_DIR}")
if(ENABLE_MPI)
  add_test(NAME test_fermion_ctm_warm_start
           COMMAND ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 1
                   ${MPIEXEC_POSTFLAGS}
                   $<TARGET_FILE:test_fermion_ctm_warm_start>)
else()
  add_test(NAME test_fermion_ctm_warm_start
           COMMAND $<TARGET_FILE:test_fermion_ctm_warm_start>)
endif()
set_tests_properties(test_fermion_ctm_warm_start PROPERTIES TIMEOUT 1800)

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

foreach(name quasi1d_ctm onesite_rdm_output)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${name}.py.in
                 ${CMAKE_CURRENT_BINARY_DIR}/${name}.py @ONLY)
  add_test(NAME ${name} COMMAND ${TENES_PYTHON_EXECUTABLE}
                                ${CMAKE_CURRENT_BINARY_DIR}/${name}.py)
endforeach()

# quasi1d_ctm reads the saved tensors back with its own mptensor reader,
# which only understands single-rank files. On more than one rank it
# exits 77 so that the test is skipped instead of reporting a bogus
# reference value or turning red.
set_tests_properties(quasi1d_ctm PROPERTIES SKIP_RETURN_CODE 77)

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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fermion/boson_equivalence_full.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fermion/boson_equivalence_full.py @ONLY)
add_test(NAME BosonEquivalenceFull
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/boson_equivalence_full.py)
set_tests_properties(BosonEquivalenceFull PROPERTIES TIMEOUT 900)

# S-2 of
# docs/superpowers/specs/2026-09-07-fermion-fast-full-update-contract.md: the
# same boson/fermion equivalence with fastfullupdate = true. Same script, one
# argument; the fast comparison is between two different implementations of
# the partial CTM move (bosonic core::*_move on the bare Tn against fermionic
# core::*_move_single on the reduced tensors) and, unlike S-1, does not depend
# on the non-fast fermionic path completing.
add_test(NAME BosonEquivalenceFastFull
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/boson_equivalence_full.py
                 --fast)
set_tests_properties(BosonEquivalenceFastFull PROPERTIES TIMEOUT 900)

# Contract T5 (2026-09-05 revision): D = 3, chi = 12, 1000 simple-update
# steps and 4 full-update sweeps; about a minute in a Release build, and the
# TIMEOUT is well above three times that.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fermion/free_fermion_full.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_full.py @ONLY)
add_test(NAME FreeFermionFull
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_full.py)
set_tests_properties(FreeFermionFull PROPERTIES TIMEOUT 1800)

# Contract T8-iii and T8-v (E2E part): the free-fermion simple update with
# complex tensors (is_real = false) at half filling; the fold CTM has to
# converge without any norm warning, the energy and every one-site RDM trace
# have to be real, and the density has to agree with the real-tensor golden.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fermion/free_fermion_complex.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_complex.py @ONLY)
add_test(NAME FreeFermionComplex
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_complex.py)
set_tests_properties(FreeFermionComplex PROPERTIES TIMEOUT 900)

# End-to-end test of tensor_save / tensor_load in fermion mode: the parity
# ledger written next to the tensors, the reload round trip, a parameter scan
# seeded from a checkpoint, and the two checkpoints that must be refused.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fermion/free_fermion_saveload.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_saveload.py @ONLY)
add_test(NAME FreeFermionSaveLoad
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_saveload.py)
set_tests_properties(FreeFermionSaveLoad PROPERTIES TIMEOUT 1800)

# Contract S-1 and S-2 of
# docs/superpowers/specs/2026-09-06-fermion-saveload-ctm-contract.md: the
# save/load round trip FreeFermionSaveLoad does not cover, namely fermion mode
# with the CTM environment (meanfield_env = false) and a full update, which
# folds the CTM edge tensors and used to make the checkpoint unreadable.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fermion/free_fermion_saveload_ctm.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_saveload_ctm.py @ONLY)
add_test(NAME FreeFermionSaveLoadCTM
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_saveload_ctm.py)
set_tests_properties(FreeFermionSaveLoadCTM PROPERTIES TIMEOUT 1800)

# Contract S-1..S-4 of
# docs/superpowers/specs/2026-09-06-fermion-tensor-load-d-expansion-contract.md:
# a fermionic checkpoint saved at D = 2 has to be loadable by a run asking for
# D = 4, and the expansion has to reach the better fixed point. Two 1000-step
# simple updates (D = 2 then D = 4, chi = 8) take about a minute in a Release
# build with four threads; ctest runs single-threaded, so the TIMEOUT is set
# well above that.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fermion/free_fermion_dexpand.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_dexpand.py @ONLY)
add_test(NAME FreeFermionDExpand
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_dexpand.py)
set_tests_properties(FreeFermionDExpand PROPERTIES TIMEOUT 3600)

# End-to-end test of the full three-stage pipeline (tenes_simple -> tenes_std
# -> tenes) for the fermionic model, checked against the handwritten
# expert-mode sample sample/07_spinless_fermion.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fermion/free_fermion_simple.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_simple.py @ONLY)
add_test(NAME FreeFermionSimple
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_simple.py)
set_tests_properties(FreeFermionSimple PROPERTIES TIMEOUT 600)

# S-1, S-4 and S-6 of
# docs/superpowers/specs/2026-09-07-fermion-fast-full-update-contract.md: the
# free-fermion run with fastfullupdate = false against the same run with
# fastfullupdate = true. The quick condition (D = 2, chi = 8, one full-update
# sweep) takes about 15 s in a Release build; TENES_FERMION_FULL_TEST /
# TENES_FREE_FERMION_FULL switches it to the contract's reference condition
# (D = 4, chi = 16), which the contract measures at 92 s + 22 s.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fermion/fast_full_update_ab.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fermion/fast_full_update_ab.py @ONLY)
add_test(NAME FreeFermionFastFullUpdate
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/fast_full_update_ab.py)
if(TENES_FERMION_FULL_TEST)
  set_tests_properties(FreeFermionFastFullUpdate PROPERTIES TIMEOUT 3600)
else()
  set_tests_properties(FreeFermionFastFullUpdate PROPERTIES TIMEOUT 900)
endif()

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
