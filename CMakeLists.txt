cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(TNSolve NONE)
set(TNSolve_VERSION 0.1-Beta)

message(STATUS "CMake version: " ${CMAKE_VERSION})

find_package(Git REQUIRED)

execute_process(COMMAND ${GIT_EXECUTABLE} submodule update -i)
