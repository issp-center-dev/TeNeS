function(git_clone uri name branch)
    if (NOT IS_DIRECTORY ${name})
        message(STATUS "git clone ${uri}")
        execute_process(COMMAND ${GIT_EXECUTABLE} clone ${uri} ${name}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
    endif ()
    message(STATUS "git checkout ${branch}")
    execute_process(COMMAND ${GIT_EXECUTABLE} checkout ${branch}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${name}
                   )
endfunction()

if (NOT MPTENSOR_ROOT)
    git_clone(https://github.com/smorita/mptensor mptensor v0.3.0)
endif ()

if (NOT CPPTOML_ROOT)
    git_clone(https://github.com/skystrife/cpptoml cpptoml v0.1.1)
endif ()

if (USE_SANITIZER)
    git_clone(https://github.com/arsenm/sanitizers-cmake sanitizers-cmake master)
endif()
