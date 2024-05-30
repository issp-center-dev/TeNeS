for L in 8 32; do
  for G in 0.0 0.5 0.8 2.0; do
    file_ene=ene-L${L}-G${G}.dat
    file_magx=magx-L${L}-G${G}.dat
    file_spec=spec-L${L}-G${G}.dat
    rm -f $file_ene
    rm -f $file_magx
    rm -f $file_spec

    dla_alg -l lattice_L${L}.toml -h hamiltonian_G${G}.toml

    for T in 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.25 1.5 1.75 2.0 5.0 10.0; do
      cp common_param.in param.in
      beta=`echo "1.0/$T" | bc -l`
      echo "beta = $beta" >> param.in
      srun dla param.in
      output=sample-L${L}-G${G}-T${T}.log
      mv sample.log $output
      grep "ene" $output | awk --assign="T=$T" '{print T,$4,$5}' >> $file_ene
      grep "bmzu" $output | awk --assign="T=$T" '{print T,$4,$5}' >> $file_magx
      grep "spe" $output | awk --assign="T=$T" '{print T,$4,$5}' >> $file_spec
    done
  done
done
