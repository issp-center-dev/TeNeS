add_library(find_contraction_path STATIC finder.cpp)

target_compile_options(find_contraction_path PUBLIC -Wall)
target_compile_options(find_contraction_path PUBLIC $<$<CONFIG:Debug>: -Wextra -Wno-unused-parameter -Wno-sign-compare -Wpedantic >)
target_compile_options(find_contraction_path PUBLIC ${TENES_COMPILE_OPTIONS})
target_include_directories(find_contraction_path PUBLIC ${BOOST_INCLUDE_DIR} ${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(find_contraction_path PUBLIC ${OpenMP_CXX_LIBRARIES})
