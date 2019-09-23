add_custom_target(tenes_simple ALL
    COMMAND echo '\#!${PYTHON_EXECUTABLE}'  > ${CMAKE_CURRENT_BINARY_DIR}/tenes_simple
    COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/tenes_simple.py >> ${CMAKE_CURRENT_BINARY_DIR}/tenes_simple
    COMMAND chmod +x ${CMAKE_CURRENT_BINARY_DIR}/tenes_simple
    )
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/tenes_simple DESTINATION bin)
