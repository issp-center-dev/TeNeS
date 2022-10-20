set(SPHINX_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/_build)
set(SPHINX_CACHE_DIR ${CMAKE_CURRENT_BINARY_DIR}/_doctrees)
set(SPHINX_HTML_DIR ${CMAKE_CURRENT_BINARY_DIR}/html)
set(SPHINX_PDF_DIR ${CMAKE_CURRENT_BINARY_DIR}/pdf)
set(TENES_DOC_DIR_JP share/tenes/${TENES_VERSION}/doc/ja)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(doc-ja-html ALL
  COMMAND
    env PYTHONPATH="${pythonpath_build}:$ENV{PYTHONPATH}"
    ${SPHINX_EXECUTABLE}
    -b html
    -d ${SPHINX_CACHE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${SPHINX_HTML_DIR}
  )

if(LATEX_FOUND)
  add_custom_target(doc-ja-pdf ALL
    COMMAND
      env PYTHONPATH="${pythonpath_build}:$ENV{PYTHONPATH}"
      ${SPHINX_EXECUTABLE}
      -b latex
      -d ${SPHINX_CACHE_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${SPHINX_PDF_DIR}
    COMMAND
      cd pdf && make
    )
  install(FILES ${SPHINX_PDF_DIR}/TeNeS_ja.pdf
    DESTINATION ${TENES_DOC_DIR_JP}
    )

  add_custom_target(doc-ja DEPENDS doc-ja-html doc-ja-pdf)
else(LATEX_FOUND)
  add_custom_target(doc-ja DEPENDS doc-ja-html)
endif(LATEX_FOUND)

install(DIRECTORY ${SPHINX_HTML_DIR}
  DESTINATION ${TENES_DOC_DIR_JP}
  PATTERN "_sources" EXCLUDE
  PATTERN ".buildinfo" EXCLUDE
  PATTERN "objects.inv" EXCLUDE
  )
