add_library(optimize_contraction_order STATIC optimize.cpp)

target_compile_options(optimize_contraction_order PUBLIC -Wall)
target_compile_options(optimize_contraction_order PUBLIC $<$<CONFIG:Debug>: -Wextra -Wno-unused-parameter -Wno-sign-compare -Wpedantic >)
target_compile_options(optimize_contraction_order PUBLIC ${TENES_COMPILE_OPTIONS})
target_include_directories(optimize_contraction_order PUBLIC ${BOOST_INCLUDE_DIR} ${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(optimize_contraction_order PUBLIC ${OpenMP_CXX_LIBRARIES})
