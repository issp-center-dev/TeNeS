cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(TNSolve CXX)
set(TNSolve_VERSION 0.1-Beta)

message(STATUS "CMake version: " ${CMAKE_VERSION})

# git submodule update
# ref: https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
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
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/mptensor/Makefile")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

include(${PROJECT_SOURCE_DIR}/config/cxx11.cmake)
cxx11()

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

add_subdirectory("${PROJECT_SOURCE_DIR}/deps/mptensor")
set(MPTENSOR_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/deps/mptensor/src/)
set(MPTENSOR_LIB ${CMAKE_CURRENT_BINARY_DIR}/deps/mptensor/src/libmptensor.a)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)
