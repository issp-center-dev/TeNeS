unset(src_files)
set(src_files ctm.cpp correlation.cpp)
foreach(row 1 2 3 4)
  foreach(col 1 2 3 4)
    list(APPEND src_files ctm_${row}x${col}.cpp)
  endforeach()
endforeach()

add_library(iTPS_contract_itps_ctm STATIC ${src_files})
target_link_libraries(iTPS_contract_itps_ctm PUBLIC tensor)
