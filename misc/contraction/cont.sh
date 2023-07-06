for nrow in 1 2 3; do
  for ncol in 1 2 3; do
    julia ./cont.jl --output=${nrow}_${ncol}.dat --pass_as_vector --tdt_path="./tdt.py" ${nrow} ${ncol}
  done
done
