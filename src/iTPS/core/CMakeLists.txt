add_subdirectory(contract_itps_ctm)
add_subdirectory(contract_itps_mf)
add_subdirectory(contract_density_ctm)

add_library(iTPS_core STATIC ctm.cpp ctm_single.cpp simple_update.cpp
                             full_update.cpp local_gauge.cpp)
target_link_libraries(
  iTPS_core PUBLIC iTPS_contract_itps_ctm iTPS_contract_itps_mf
                   iTPS_contract_density_ctm)
