macro(openmp)
  set(OMP_FLAG "")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    if("${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS "15.0.0.20140528")
      set(OMP_FLAG "-openmp")
    else()
      set(OMP_FLAG "-qopenmp")
    endif()
    set(OPENMP_FOUND ON)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(_OMP_FLAG "-Xpreprocessor -fopenmp -lomp")
    if(NOT TENES_OpenMP_INCLUDE_DIR)
      set(_OMP_INCLUDE_DIR "/usr/local/include")
    else()
      set(_OMP_INCLUDE_DIR ${TENES_OpenMP_INCLUDE_DIR})
    endif()
    if(NOT TENES_OpenMP_LIBRARY_DIR)
      set(_OMP_LIBRARY_DIR "/usr/local/lib")
    else()
      set(_OMP_LIBRARY_DIR ${TENES_OpenMP_LIBRARY_DIR})
    endif()
    try_compile(_AppleClangOMP
      ${CMAKE_CURRENT_BINARY_DIR}
      SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/config/check_omp.cpp
      CMAKE_FLAGS "-DCOMPILE_DEFINITIONS=${_OMP_FLAG}"
      CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${_OMP_INCLUDE_DIR}"
      CMAKE_FLAGS "-DLINK_DIRECTORIES=${_OMP_LIBRARY_DIR}"
      OUTPUT_VARIABLE LOG)
    if(_AppleClangOMP)
      set(OMP_FLAG -Xpreprocessor -fopenmp)
      message(STATUS "Found OpenMP: ${OMP_FLAG}")
      set(OpenMP_CXX_FLAGS -Xpreprocessor -fopenmp)
      set(OpenMP_CXX_LIBRARIES -lomp)
      set(OpenMP_CXX_INCLUDE_DIRS ${_OMP_INCLUDE_DIR})
      set(OpenMP_CXX_LIBRARY_DIRS ${_OMP_LIBRARY_DIR})
      set(OPENMP_FOUND ON)
    else()
      message(STATUS "Could NOT find OpenMP")
      set(OPENMP_FOUND OFF)
    endif()
  else()
    find_package(OpenMP)
    if(OPENMP_FOUND)
      set(OMP_FLAG ${OpenMP_CXX_FLAGS})
    endif(OPENMP_FOUND)
  endif()

  if(NOT DEFINED OpenMP_CXX_LIBRARIES)
    set(OpenMP_CXX_LIBRARIES ${OMP_FLAG})
  endif()
endmacro()