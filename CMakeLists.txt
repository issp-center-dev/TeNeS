cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(TNSolve CXX)
set(TNSolve_VERSION 0.1-Beta)

message(STATUS "CMake version: " ${CMAKE_VERSION})

option(GIT_SUBMODULE "Check submodules during build" ON)
option(MPTENSOR_ROOT "External directory where mptensor is installed into" OFF)

# git submodule update
# ref: https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
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
endif()
foreach(submod mptensor toml11)
    if(NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/${submod}/CMakeLists.txt")
    message(FATAL_ERROR "The submodule ${submod} were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()
endforeach()

include(${PROJECT_SOURCE_DIR}/config/cxx11.cmake)
cxx11()

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/deps/mptensor/config" ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/deps/sanitizers-cmake/cmake" ${CMAKE_MODULE_PATH})
find_package(Sanitizers)

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

if(NOT MPTENSOR_ROOT)
    add_subdirectory("${PROJECT_SOURCE_DIR}/deps/mptensor")
    set(MPTENSOR_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/deps/mptensor/src/)
    set(MPTENSOR_LIB_DIR ${CMAKE_BINARY_DIR}/deps/mptensor/src)
else()
    message(STATUS "To use pre-built mptensor at ${MPTENSOR_ROOT}")
    set(MPTENSOR_INCLUDE_DIR ${MPTENSOR_ROOT}/include)
    set(MPTENSOR_LIB_DIR ${MPTENSOR_ROOT}/lib)
    find_package(ScaLAPACK REQUIRED)
    string(REGEX REPLACE ";" " " SCALAPACK_LIBS "${SCALAPACK_LIBRARIES}")
endif()

if(NOT EXISTS ${MPTENSOR_INCLUDE_DIR}/mptensor.hpp)
    message(FATAL_ERROR "${MPTENSOR_INCLUDE_DIR}/mptensor.hpp is not found.")
endif()
if(NOT EXISTS ${MPTENSOR_LIB_DIR}/libmptensor.a)
    message(FATAL_ERROR "${MPTENSOR_LIB_DIR}/libmptensor.a is not found.")
endif()

include_directories(${PROJECT_SOURCE_DIR}/deps/)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)
