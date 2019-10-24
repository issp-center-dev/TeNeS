## This file is originally included in mptensor
## https://github.com/smorita/mptensor

# - Try to find ScaLAPACK
# Once done this will define
#
#  SCALAPACK_FOUND        - system has ScaLAPACK
#  SCALAPACK_LIBRARIES     - libraries for ScaLAPACK

if(DEFINED SCALAPACK_FOUND)
  return()
endif(DEFINED SCALAPACK_FOUND)
  
message(STATUS "Checking for ScaLAPACK library")

if(DEFINED SCALAPACK_LIB)
  set(SCALAPACK_FOUND TRUE)
  set(SCALAPACK_LIBRARIES ${SCALAPACK_LIB})
  message(STATUS "ScaLAPACK libraries: ${SCALAPACK_LIBRARIES}")
  return()
endif(DEFINED SCALAPACK_LIB)

unset(_SCALAPACK_LIBRARY)

if(DEFINED BLAS_mkl_core_LIBRARY)

  find_library(_SCALAPACK_LIBRARY
    NAMES mkl_scalapack_lp64
    PATHS $ENV{MKLROOT}/lib/intel64 $ENV{MKLROOT}/lib/em64t
    DOC "The ScaLAPACK library")

else(DEFINED BLAS_mkl_core_LIBRARY)

  # Standard search path
  set(_PATHS "")
  if(SCALAPACK_DIR)
    set(_PATHS ${SCALAPACK_DIR})
  else(SCALAPACK_DIR)
    list(APPEND _PATHS
  	 ${SCALAPACK_ROOT}/${CMAKE_BUILD_TYPE}
	 ${SCALAPACK_ROOT}
  	 $ENV{SCALAPACK_ROOT}/${CMAKE_BUILD_TYPE}
	 $ENV{SCALAPACK_ROOT}
  	 ${ROKKO_SOLVER_ROOT}/scalapack/${CMAKE_BUILD_TYPE}
	 ${ROKKO_SOLVER_ROOT}/scalapack
  	 $ENV{ROKKO_SOLVER_ROOT}/scalapack/${CMAKE_BUILD_TYPE}
	 $ENV{ROKKO_SOLVER_ROOT}/scalapack
	 ${CMAKE_INSTALL_PREFIX}/scalapack/${CMAKE_BUILD_TYPE}
	 ${CMAKE_INSTALL_PREFIX}/${CMAKE_BUILD_TYPE}
	 $ENV{HOME}/rokko/scalapack/${CMAKE_BUILD_TYPE}
	 $ENV{HOME}/rokko/scalapack
	 /opt/rokko/scalapack/${CMAKE_BUILD_TYPE}
	 /opt/rokko/scalapack
	 /opt/rokko/${CMAKE_BUILD_TYPE}
	 /opt/rokko
	 /opt/local /opt
	 )
    list(APPEND _PATHS /usr/lib64/openmpi) # for CentOS
  endif(SCALAPACK_DIR)

  foreach (_PATH ${_PATHS})
    list(APPEND _LIBPATHS "${_PATH}/lib")
  endforeach()

  find_library(_SCALAPACK_LIBRARY
    NAMES scalapack scalapack-openmpi scalapack-mpich
    PATHS ${_LIBPATHS}
    DOC "The ScaLAPACK library")
endif(DEFINED BLAS_mkl_core_LIBRARY)

find_package_handle_standard_args(ScaLAPACK
  FOUND_VAR SCALAPACK_FOUND
  REQUIRED_VARS _SCALAPACK_LIBRARY
  )

unset(_SCALAPACK_LIBRARIES)

if(SCALAPACK_FOUND)
  list(APPEND SCALAPACK_LIBRARIES ${_SCALAPACK_LIBRARY})
  if(DEFINED BLAS_mkl_core_LIBRARY)
    # Check whether SGI MPT is used
    try_compile(_SGI_MPT
      ${CMAKE_CURRENT_BINARY_DIR}
      SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/config/check_sgimpt.cpp
      CMAKE_FLAGS
	"-DINCLUDE_DIRECTORIES:STRING=${MPI_CXX_INCLUDE_DIRS}"
      LINK_LIBRARIES ${MPI_CXX_LIBRARIES}
      OUTPUT_VARIABLE LOG)
    if(_SGI_MPT)
      find_library(_SCALAPACK_BLACS_LIBRARY
	NAMES mkl_blacs_sgimpt_lp64
	PATHS $ENV{MKLROOT}/lib/intel64 $ENV{MKLROOT}/lib/em64t
	DOC "The BLACS library")
      MESSAGE(STATUS "SGI MPT is used")
    else(_SGI_MPT)
       try_compile(_OPENMPI
	 ${CMAKE_CURRENT_BINARY_DIR}
	 SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/config/check_openmpi.cpp
	 CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${MPI_CXX_INCLUDE_DIRS}"
	 OUTPUT_VARIABLE LOG)
      if(_OPENMPI)
	find_library(_SCALAPACK_BLACS_LIBRARY
	  NAMES mkl_blacs_openmpi_lp64
	  PATHS $ENV{MKLROOT}/lib/intel64 $ENV{MKLROOT}/lib/em64t
	  DOC "The BLACS library")
	MESSAGE(STATUS "OpenMPI is used")
      else(_OPENMPI)
	find_library(_SCALAPACK_BLACS_LIBRARY
	  NAMES mkl_blacs_intelmpi_lp64
	  PATHS $ENV{MKLROOT}/lib/intel64 $ENV{MKLROOT}/lib/em64t
	  DOC "The BLACS library")
	MESSAGE(STATUS "Intel MPI/MPICH2/MVAPICH is used")
      endif(_OPENMPI)
    endif(_SGI_MPT)
    if(_SCALAPACK_BLACS_LIBRARY)
      list(APPEND SCALAPACK_LIBRARIES ${_SCALAPACK_BLACS_LIBRARY})
    endif(_SCALAPACK_BLACS_LIBRARY)


  else(DEFINED BLAS_mkl_core_LIBRARY)
    find_library(_BLACS_LIBRARY
      NAMES blacs blacs-openmpi blacs-mpich
      PATHS ${_LIBPATHS}
      DOC "The ScaLAPACK BLACS library")
    if(_BLACS_LIBRARY)
      list(APPEND SCALAPACK_LIBRARIES ${_BLACS_LIBRARY})
    endif(_BLACS_LIBRARY)

    find_library(_BLACSCINIT_LIBRARY
      NAMES blacsCinit blacsCinit-openmpi blacsCinit-mpich
      PATHS ${_LIBPATHS}
      DOC "The ScaLAPACK BLACS Cinit library")
    if(_BLACSCINIT_LIBRARY)
      list(APPEND SCALAPACK_LIBRARIES ${_BLACSCINIT_LIBRARY})
    endif(_BLACSCINIT_LIBRARY)

  endif(DEFINED BLAS_mkl_core_LIBRARY)
endif(SCALAPACK_FOUND)

