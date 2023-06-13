# obtain ground state before quency
tenes_simple simple.toml
tenes_std std.toml
tenes input.toml

for name in weak strong; do
  tenes_simple -o std_te_${name}.toml simple_te_${name}.toml
  tenes_std -o input_te_${name}.toml std_te_${name}.toml
  tenes input_te_${name}.toml

  awk '$2 == 2 {print $1, $3, $4}' output_te_${name}/TE_density.dat > magnetization_${name}.dat
done
