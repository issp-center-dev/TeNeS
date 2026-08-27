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
    # Apple clang ships without an OpenMP runtime, so libomp has to be
    # installed separately (Homebrew: `brew install libomp`).  Homebrew's
    # mpicxx wraps Apple clang, so MPI builds land here as well.
    #
    # Locate the libomp prefix.  Homebrew is /opt/homebrew on Apple silicon and
    # /usr/local on Intel, so ask `brew` rather than hard-coding either one.
    set(_OMP_PREFIXES)
    if(OpenMP_ROOT)
      list(APPEND _OMP_PREFIXES ${OpenMP_ROOT})
    endif()
    if(TENES_OpenMP_LIBRARY_DIR)
      get_filename_component(_OMP_HINT ${TENES_OpenMP_LIBRARY_DIR} DIRECTORY)
      list(APPEND _OMP_PREFIXES ${_OMP_HINT})
    endif()
    find_program(TENES_BREW_EXECUTABLE brew)
    if(TENES_BREW_EXECUTABLE)
      execute_process(COMMAND ${TENES_BREW_EXECUTABLE} --prefix libomp
                      OUTPUT_VARIABLE _OMP_BREW_PREFIX
                      OUTPUT_STRIP_TRAILING_WHITESPACE
                      ERROR_QUIET)
      if(_OMP_BREW_PREFIX)
        list(APPEND _OMP_PREFIXES ${_OMP_BREW_PREFIX})
      endif()
    endif()
    list(APPEND _OMP_PREFIXES
         /opt/homebrew/opt/libomp /usr/local/opt/libomp /usr/local)

    set(_OMP_PREFIX "")
    foreach(_prefix IN LISTS _OMP_PREFIXES)
      if(EXISTS ${_prefix}/include/omp.h)
        set(_OMP_PREFIX ${_prefix})
        break()
      endif()
    endforeach()

    if(CMAKE_VERSION VERSION_LESS 3.12)
      # Both OpenMP_ROOT (policy CMP0074) and the SHELL: option prefix used
      # below require CMake 3.12, so probe the flags by hand on older CMake.
      # This only feeds TeNeS itself: mptensor runs its own
      # find_package(OpenMP REQUIRED), which without OpenMP_ROOT cannot find a
      # Homebrew libomp and stops the configuration.
      message(STATUS "CMake ${CMAKE_VERSION} predates OpenMP_ROOT; "
                     "CMake >= 3.12 is recommended for Apple clang builds")
      if(TENES_OpenMP_INCLUDE_DIR)
        set(_OMP_INCLUDE_DIR ${TENES_OpenMP_INCLUDE_DIR})
      else()
        set(_OMP_INCLUDE_DIR ${_OMP_PREFIX}/include)
      endif()
      if(TENES_OpenMP_LIBRARY_DIR)
        set(_OMP_LIBRARY_DIR ${TENES_OpenMP_LIBRARY_DIR})
      else()
        set(_OMP_LIBRARY_DIR ${_OMP_PREFIX}/lib)
      endif()
      # try_compile takes a single CMAKE_FLAGS list; repeating the keyword
      # would drop all but the last entry.
      try_compile(_AppleClangOMP
        ${CMAKE_CURRENT_BINARY_DIR}
        SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/config/check_omp.cpp
        CMAKE_FLAGS "-DCOMPILE_DEFINITIONS=-Xpreprocessor -fopenmp -lomp"
                    "-DINCLUDE_DIRECTORIES=${_OMP_INCLUDE_DIR}"
                    "-DLINK_DIRECTORIES=${_OMP_LIBRARY_DIR}"
        OUTPUT_VARIABLE LOG)
      if(_AppleClangOMP)
        set(OMP_FLAG -Xpreprocessor -fopenmp)
        set(OpenMP_CXX_FLAGS -Xpreprocessor -fopenmp)
        set(OpenMP_CXX_LIBRARIES -lomp)
        set(OpenMP_CXX_INCLUDE_DIRS ${_OMP_INCLUDE_DIR})
        set(OpenMP_CXX_LIBRARY_DIRS ${_OMP_LIBRARY_DIR})
        set(OPENMP_FOUND ON)
        message(STATUS "Found OpenMP: ${OMP_FLAG} (libomp in ${_OMP_LIBRARY_DIR})")
      else()
        set(OPENMP_FOUND OFF)
      endif()
    else()
      # CMake's own FindOpenMP knows how to drive Apple clang once it is told
      # where libomp lives, and it is the only detection mptensor sees -- it
      # calls find_package(OpenMP REQUIRED) itself.
      if(_OMP_PREFIX AND NOT OpenMP_ROOT)
        set(OpenMP_ROOT ${_OMP_PREFIX} CACHE PATH "Root directory of libomp")
      endif()
      find_package(OpenMP)
      if(OPENMP_FOUND)
        # Here OpenMP_CXX_FLAGS is one space-separated string ("-Xclang
        # -fopenmp").  It has to reach the compiler as two arguments, and the
        # SHELL: prefix both splits it and keeps the pair together when the
        # same flags also arrive through OpenMP::OpenMP_CXX -- plain lists get
        # de-duplicated down to a dangling "-Xclang".
        set(OMP_FLAG "SHELL:${OpenMP_CXX_FLAGS}")
      endif()
    endif()

    if(NOT OPENMP_FOUND)
      message(STATUS "Could NOT find OpenMP "
                     "(hint: brew install libomp, "
                     "or -DOpenMP_ROOT=<libomp prefix>)")
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
