# Run at build time (cmake -P):
#   -DSOURCE_DIR=... -DTENES_VERSION=... -DINPUT=... -DOUTPUT=... [-DSHEBANG=...]
# Writes OUTPUT from INPUT with @TENES_VERSION@, @TENES_GIT_HASH@ and
# @TENES_GIT_DIRTY@ filled in, and only touches OUTPUT when its content
# changes.
include("${CMAKE_CURRENT_LIST_DIR}/git_hash.cmake")
if(TENES_GIT_FAILED AND EXISTS "${OUTPUT}")
  # e.g. "sudo make install": git refuses a repository of another user.
  # Keep what the build wrote.
  return()
endif()
file(READ "${INPUT}" _content)
string(REPLACE "@TENES_VERSION@" "${TENES_VERSION}" _content "${_content}")
string(REPLACE "@TENES_GIT_HASH@" "${TENES_GIT_HASH}" _content "${_content}")
string(REPLACE "@TENES_GIT_DIRTY@" "${TENES_GIT_DIRTY}" _content "${_content}")
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
