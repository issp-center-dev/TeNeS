for L in 8 32; do
  for G in 0.0 0.5 0.8 2.0; do
    file_ene=ene-L${L}-G${G}.dat
    file_magx=magx-L${L}-G${G}.dat
    file_spec=spec-L${L}-G${G}.dat
    rm -f $file_ene
    rm -f $file_magx
    rm -f $file_spec

    dla_alg -l lattice_L${L}.toml -h hamiltonian_G${G}.toml

    for beta in 0.2 0.4 0.6 0.8 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.5 4.0; do
      cp common_param.in param.in
      echo "beta = $beta" >> param.in
      T=`echo "1.0/$beta" | bc -l`
      srun dla param.in
      output=sample-L${L}-G${G}-beta${beta}.log
      mv sample.log $output
      grep "ene" $output | awk --assign="beta=$beta" '{print beta,$4,$5}' >> $file_ene
      grep "bmzu" $output | awk --assign="beta=$beta" '{print beta,$4,$5}' >> $file_magx
      grep "spe" $output | awk --assign="beta=$beta" '{print beta,$4,$5}' >> $file_spec
    done
  done
done
