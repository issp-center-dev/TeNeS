add_subdirectory(contract_ctm)
add_subdirectory(contract_mf)

add_library(iTPS_core STATIC ctm.cpp ctm_single.cpp simple_update.cpp full_update.cpp local_gauge.cpp)
target_link_libraries(iTPS_core PUBLIC iTPS_contract_ctm iTPS_contract_mf)
