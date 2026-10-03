# Find the commit the source tree was built from.
#
#   in a git checkout : asked from git
#   in a tarball      : read from config/git_archive.txt, which "git archive"
#                       fills in (export-subst in .gitattributes)
#   otherwise         : not known
#
# input : SOURCE_DIR
# output: TENES_GIT_HASH   the full hash, empty when it is not known
#         TENES_GIT_DIRTY  "true" when the checkout has uncommitted changes
#                          to tracked files, "false" otherwise
#         TENES_GIT_FAILED ON in a git checkout that git could not read
set(TENES_GIT_HASH "")
set(TENES_GIT_DIRTY "false")
set(TENES_GIT_FAILED OFF)
if(EXISTS "${SOURCE_DIR}/.git")
  set(TENES_GIT_FAILED ON)
  find_package(Git QUIET)
  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      WORKING_DIRECTORY "${SOURCE_DIR}"
      OUTPUT_VARIABLE _hash OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET RESULT_VARIABLE _rc)
    if(_rc EQUAL 0 AND _hash MATCHES "^[0-9a-f]+$")
      set(TENES_GIT_HASH "${_hash}")
      set(TENES_GIT_FAILED OFF)
      execute_process(
        COMMAND ${GIT_EXECUTABLE} status --porcelain --untracked-files=no
                --ignore-submodules=untracked
        WORKING_DIRECTORY "${SOURCE_DIR}"
        OUTPUT_VARIABLE _status OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET RESULT_VARIABLE _rc)
      if(_rc EQUAL 0 AND NOT _status STREQUAL "")
        set(TENES_GIT_DIRTY "true")
      endif()
    endif()
  endif()
elseif(EXISTS "${SOURCE_DIR}/config/git_archive.txt")
  file(STRINGS "${SOURCE_DIR}/config/git_archive.txt" _lines LIMIT_COUNT 1)
  # left as "$Format:%H$" when the tree does not come from git archive
  if(_lines MATCHES "^[0-9a-f]+$")
    set(TENES_GIT_HASH "${_lines}")
  endif()
endif()
