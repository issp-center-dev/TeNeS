# Run at build time (cmake -P):
#   -DSOURCE_DIR=... -DTENES_VERSION=... -DINPUT=... -DOUTPUT=... [-DSHEBANG=...]
# Writes OUTPUT from INPUT with @TENES_VERSION@, @TENES_GIT_HASH@ and
# @TENES_GIT_DIRTY@ filled in, and only touches OUTPUT when its content
# changes.

# A script run with -P starts with no policy set. CMake 3 then keeps the old
# behaviour of CMP0053 and expands @VAR@ inside a quoted argument, so a
# pattern written as "@TENES_VERSION@" became the version number itself and
# nothing was replaced. The policies are set here, and the patterns are
# bracket arguments, which no version of CMake evaluates.
cmake_minimum_required(VERSION 3.8...3.14)

include("${CMAKE_CURRENT_LIST_DIR}/git_hash.cmake")
if(TENES_GIT_FAILED AND EXISTS "${OUTPUT}")
  # e.g. "sudo make install": git refuses a repository of another user.
  # Keep what the build wrote.
  return()
endif()
file(READ "${INPUT}" _content)
string(REPLACE [[@TENES_VERSION@]] "${TENES_VERSION}" _content "${_content}")
string(REPLACE [[@TENES_GIT_HASH@]] "${TENES_GIT_HASH}" _content "${_content}")
string(REPLACE [[@TENES_GIT_DIRTY@]] "${TENES_GIT_DIRTY}" _content "${_content}")
string(FIND "${_content}" [[@TENES_]] _left)
if(NOT _left EQUAL -1)
  message(FATAL_ERROR "${INPUT}: a placeholder @TENES_...@ was not filled in")
endif()
if(SHEBANG)
  set(_content "#!${SHEBANG}\n${_content}")
endif()
set(_old "")
if(EXISTS "${OUTPUT}")
  file(READ "${OUTPUT}" _old)
endif()
if(NOT _content STREQUAL _old)
  file(WRITE "${OUTPUT}" "${_content}")
endif()
