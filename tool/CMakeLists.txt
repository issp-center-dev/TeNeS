foreach(name tenes_simple tenes_std)
    add_custom_target(${name} ALL
        COMMAND echo '\#!${TENES_PYTHON_EXECUTABLE}'  > ${CMAKE_CURRENT_BINARY_DIR}/${name}
        COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/${name}.py >> ${CMAKE_CURRENT_BINARY_DIR}/${name}
        COMMAND chmod +x ${CMAKE_CURRENT_BINARY_DIR}/${name}
        )
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${name} DESTINATION bin)
endforeach()
