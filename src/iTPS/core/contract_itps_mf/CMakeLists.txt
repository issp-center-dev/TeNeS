unset(src_files)
set(src_files mf.cpp correlation.cpp)
foreach(row 1 2 3 4)
  foreach(col 1 2 3 4)
    list(APPEND src_files mf_${row}x${col}.cpp)
  endforeach()
endforeach()

add_library(iTPS_contract_itps_mf STATIC ${src_files})

target_link_libraries(iTPS_contract_itps_mf PUBLIC tensor)
