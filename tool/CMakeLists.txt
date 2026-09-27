# "/usr/bin/env python3", the default, is a list of two elements
string(REPLACE ";" " " tenes_python_shebang "${TENES_PYTHON_EXECUTABLE}")

foreach(name tenes_simple tenes_std)
    add_custom_target(${name} ALL
        COMMAND ${CMAKE_COMMAND}
            -DSOURCE_DIR=${TeNeS_SOURCE_DIR}
            -DTENES_VERSION=${TENES_VERSION}
            -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/${name}.py
            -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/${name}
            "-DSHEBANG=${tenes_python_shebang}"
            -P ${TeNeS_SOURCE_DIR}/config/write_version.cmake
        COMMAND chmod +x ${CMAKE_CURRENT_BINARY_DIR}/${name}
        VERBATIM
        )
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${name} DESTINATION bin)
endforeach()
