cmake_minimum_required(VERSION 3.6...3.14)
message(STATUS "CMake version: " ${CMAKE_VERSION})
project(TeNeS CXX)
set(TENES_VERSION 1.3.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(GIT_SUBMODULE "Check submodules during build" ON)
option(MPTENSOR_ROOT "External directory where mptensor is installed into" OFF)
option(CPPTOML_ROOT "External directory where cpptoml is download into" OFF)
option(TENES_PYTHON_EXECUTABLE "Path to Python interpreter" OFF)
option(Testing "Enable tests" ON)
option(Document "Build docs" OFF)
# option(ENABLE_MPI "MPI" ON)

if(NOT DEFINED ENABLE_MPI)
    if(APPLE)
        set(ENABLE_MPI OFF CACHE BOOL "Use ScaLAPACK" FORCE)
    else()
        set(ENABLE_MPI ON CACHE BOOL "Use ScaLAPACK" FORCE)
    endif()
else()
    set(ENABLE_MPI ${ENABLE_MPI} CACHE BOOL "Use ScaLAPACK" FORCE)
endif()

if(Testing)
    enable_testing()
endif()

# git submodule update
# ref: https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    set(DEPS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps")
    # Update submodules as needed
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()

    foreach(submod mptensor cpptoml)
        if(NOT EXISTS "${DEPS_DIR}/${submod}/CMakeLists.txt")
            message(FATAL_ERROR "The submodule ${submod} were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
        endif()
    endforeach()
endif()
set(DEPS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps")

include(${CMAKE_CURRENT_SOURCE_DIR}/config/openmp.cmake)
openmp()

if(NOT DEFINED OpenMP_CXX_LIBRARIES)
    set(OpenMP_CXX_LIBRARIES "${OMP_FLAG}")
endif()

find_package(LAPACK REQUIRED)

if(ENABLE_MPI)
    if(APPLE)
        message(WARNING "Apple Accelerate BLAS/LAPACK + ScaLAPACK may fail. We STRONGLY recommend disabling MPI")
    endif()

    find_package(MPI REQUIRED)
    if(NOT MPI_CXX_INCLUDE_DIRS)
        if(MPI_CXX_INCLUDE_PATH)
            set(MPI_CXX_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
        else(MPI_CXX_INCLUDE_PATH)
            set(MPI_CXX_INCLUDE_DIRS ${MPI_INCLUDE_PATH})
        endif(MPI_CXX_INCLUDE_PATH)
    endif(NOT MPI_CXX_INCLUDE_DIRS)
    if(NOT MPI_CXX_LIBRARIES)
        set(MPI_CXX_LIBRARIES ${MPI_LIBRARIES})
    endif(NOT MPI_CXX_LIBRARIES)
    if(TENES_MPIEXEC)
        set(MPIEXEC ${TENES_MPIEXEC})
    endif()
endif()

if(NOT TENES_PYTHON_EXECUTABLE)
    set(TENES_PYTHON_EXECUTABLE /usr/bin/env python3)
endif()

if(NOT MPTENSOR_ROOT)
    set(MPTENSOR_BUILD_LIBS_ONLY ON CACHE BOOL "Build only the library in mptensor")
    set(CMAKE_MODULE_PATH "${DEPS_DIR}/mptensor/config" ${CMAKE_MODULE_PATH})
    add_subdirectory("${DEPS_DIR}/mptensor")
    set(MPTENSOR_INCLUDE_DIR ${DEPS_DIR}/mptensor/include)
    set(MPTENSOR_LIB_DIR ${DEPS_DIR}/mptensor/src)
else()
    message(STATUS "To use pre-built mptensor at ${MPTENSOR_ROOT}")
    set(MPTENSOR_INCLUDE_DIR ${MPTENSOR_ROOT}/include)
    set(MPTENSOR_LIB_DIR ${MPTENSOR_ROOT}/lib)
    set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/config" ${CMAKE_MODULE_PATH})
    set(SCALAPACK_LIBS "")
    if(ENABLE_MPI)
        find_package(ScaLAPACK REQUIRED)
        string(REGEX REPLACE ";" " " SCALAPACK_LIBS "${SCALAPACK_LIBRARIES}")
    endif()
    if(NOT EXISTS ${MPTENSOR_LIB_DIR}/libmptensor.a)
        message(FATAL_ERROR "${MPTENSOR_LIB_DIR}/libmptensor.a is not found.")
    endif()
    if(NOT EXISTS ${MPTENSOR_INCLUDE_DIR}/mptensor/mptensor.hpp)
        message(FATAL_ERROR "${MPTENSOR_INCLUDE_DIR}/mptensor/mptensor.hpp is not found.")
    endif()
endif()

if(CPPTOML_ROOT)
    set(CPPTOML_INCLUDE_DIR ${CPPTOML_ROOT}/include)
else()
    set(CPPTOML_INCLUDE_DIR ${DEPS_DIR}/cpptoml/include)
endif()

set(BOOST_INCLUDE_DIR ${DEPS_DIR}/boost/include)

message(STATUS "MPTENSOR_INCLUDE_DIR: ${MPTENSOR_INCLUDE_DIR}")
message(STATUS "MPTENSOR_LIB_DIR: ${MPTENSOR_LIB_DIR}")
message(STATUS "CPPTOML_INCLUDE_DIR: ${CPPTOML_INCLUDE_DIR}")


set(TENES_COMPILE_OPTIONS)
if(OPENMP_FOUND)
    set(TENES_COMPILE_OPTIONS ${OMP_FLAG} ${TENES_COMPILE_OPTIONS})
else()
    set(TENES_COMPILE_OPTIONS -D_NO_OMP ${TENES_COMPILE_OPTIONS})
endif()

if(NOT ENABLE_MPI)
    set(TENES_COMPILE_OPTIONS -D_NO_MPI ${TENES_COMPILE_OPTIONS})
endif()

message(STATUS "TENES_COMPILE_OPTIONS: ${TENES_COMPILE_OPTIONS}")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tool)

if(Testing)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)
endif()

if(Document)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/docs/sphinx)
endif()

# TeNeSConfig.cmake
include(CMakePackageConfigHelpers)
set(TENESCONFIG_INCLUDE_DIRS include)
set(TENESCONFIG_LIBRARY_DIR lib)
configure_package_config_file(config/TeNeSConfig.cmake.in 
    ${CMAKE_CURRENT_BINARY_DIR}/TeNeSConfig.cmake
    INSTALL_DESTINATION share/cmake/tenes
    PATH_VARS TENESCONFIG_INCLUDE_DIRS TENESCONFIG_LIBRARY_DIR
    )
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/TeNeSConfigVersion.cmake
    VERSION ${TENES_VERSION}
    COMPATIBILITY SameMajorVersion
    )
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/TeNeSConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/TeNeSConfigVersion.cmake
    DESTINATION share/cmake/tenes
    )

# sample files
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/sample DESTINATION share/tenes)
