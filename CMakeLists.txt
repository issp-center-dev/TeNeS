cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
message(STATUS "CMake version: " ${CMAKE_VERSION})
project(TeNeS CXX)
set(TeNeS_VERSION 0.1-Beta)

if(${CMAKE_VERSION} VERSION_LESS 3.4)
    # old FindLAPACK.cmake requires C-lang enabled...
    enable_language(C)
endif()

option(GIT_SUBMODULE "Check submodules during build" ON)
option(MPTENSOR_ROOT "External directory where mptensor is installed into" OFF)
option(CPPTOML_ROOT "External directory where cpptoml is download into" OFF)
option(Testing "Enable tests" ON)
option(Document "Build docs" OFF)

if(Testing)
    enable_testing()
endif()


set(USE_SANITIZER OFF)
if (SANITIZE_ADDRESS OR SANITIZE_MEMORY OR SANITIZE_THREAD OR SANITIZE_UNDEFINED)
    set(USE_SANITIZER ON)
endif ()

# git submodule update
# ref: https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    set(DEPS_DIR "${PROJECT_SOURCE_DIR}/deps/submoduled")
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

else()
    set(DEPS_DIR "${PROJECT_SOURCE_DIR}/deps/cloned")
    add_subdirectory("${PROJECT_SOURCE_DIR}/deps/cloned")
endif()


include(${PROJECT_SOURCE_DIR}/config/cxx11.cmake)
cxx11()

if (USE_SANITIZER)
    set(CMAKE_MODULE_PATH "${DEPS_DIR}/sanitizers-cmake/cmake" ${CMAKE_MODULE_PATH})
    find_package(Sanitizers)
endif()

set(OMP_FLAG "")
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    if("${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS "15.0.0.20140528")
        set(OMP_FLAG "-openmp")
    else()
        set(OMP_FLAG "-qopenmp")
    endif()
else()
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(OMP_FLAG "${OpenMP_CXX_FLAGS}")
    endif(OPENMP_FOUND)
endif()

if(NOT DEFINED OpenMP_CXX_LIBRARIES)
    set(OpenMP_CXX_LIBRARIES "${OMP_FLAG}")
endif()

find_package(LAPACK REQUIRED)
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

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    find_package(PythonInterp REQUIRED)
else()
    find_package(Python COMPONENTS Interpreter REQUIRED)
    set(PYTHON_EXECUTABLE ${Python_EXECUTABLE})
endif()

if(NOT MPTENSOR_ROOT)
    set(CMAKE_MODULE_PATH "${DEPS_DIR}/mptensor/config" ${CMAKE_MODULE_PATH})
    add_subdirectory("${DEPS_DIR}/mptensor")
    set(MPTENSOR_INCLUDE_DIR ${DEPS_DIR}/mptensor/src/)
    set(MPTENSOR_LIB_DIR ${DEPS_DIR}/mptensor/src)
else()
    message(STATUS "To use pre-built mptensor at ${MPTENSOR_ROOT}")
    set(MPTENSOR_INCLUDE_DIR ${MPTENSOR_ROOT}/include)
    set(MPTENSOR_LIB_DIR ${MPTENSOR_ROOT}/lib)
    set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/config" ${CMAKE_MODULE_PATH})
    find_package(ScaLAPACK REQUIRED)
    string(REGEX REPLACE ";" " " SCALAPACK_LIBS "${SCALAPACK_LIBRARIES}")
    if(NOT EXISTS ${MPTENSOR_LIB_DIR}/libmptensor.a)
        message(FATAL_ERROR "${MPTENSOR_LIB_DIR}/libmptensor.a is not found.")
    endif()
    if(NOT EXISTS ${MPTENSOR_INCLUDE_DIR}/mptensor.hpp)
        message(FATAL_ERROR "${MPTENSOR_INCLUDE_DIR}/mptensor.hpp is not found.")
    endif()
endif()

if(CPPTOML_ROOT)
    set(CPPTOML_INCLUDE_DIR ${CPPTOML_ROOT}/include)
else()
    set(CPPTOML_INCLUDE_DIR ${DEPS_DIR}/cpptoml/include)
endif()

message(STATUS "MPTENSOR_INCLUDE_DIR: ${MPTENSOR_INCLUDE_DIR}")
message(STATUS "MPTENSOR_LIB_DIR: ${MPTENSOR_LIB_DIR}")
message(STATUS "CPPTOML_INCLUDE_DIR: ${CPPTOML_INCLUDE_DIR}")


add_subdirectory(${PROJECT_SOURCE_DIR}/src)
add_subdirectory(${PROJECT_SOURCE_DIR}/tool)

if(Testing)
    add_subdirectory(${PROJECT_SOURCE_DIR}/test)
endif()

if(Document)
    add_subdirectory(${PROJECT_SOURCE_DIR}/docs/sphinx)
endif()
