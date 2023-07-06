for nrow in 1 2 3 4; do
  for ncol in 1 2 3 4; do
    julia ./cont_mf.jl --output=${nrow}_${ncol}_mf.dat --vdim=8 --pass_as_vector --tdt_path="./tdt.py" ${nrow} ${ncol}
  done
done
