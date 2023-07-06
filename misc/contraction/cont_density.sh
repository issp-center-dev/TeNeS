for nrow in 1 2 3 4; do
  for ncol in 1 2 3 4; do
    julia ./cont_density.jl --output=${nrow}_${ncol}_den.dat --pass_as_vector --tdt_path="./tdt.py" ${nrow} ${ncol}
  done
done
