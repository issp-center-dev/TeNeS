for name in gapless gapfull; do
  tenes_simple -o std_${name}.toml simple_${name}.toml
  tenes_std -o input_${name}.toml std_${name}.toml
  tenes input_${name}.toml

  awk '$2 == 0 {print $1, $3}' output_${name}/FT_density.dat > energy_${name}.dat
  python3 ./calcspec.py $name
done
