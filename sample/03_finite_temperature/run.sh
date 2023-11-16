for name in zero weak middle strong; do
  tenes_simple -o std_ft_${name}.toml simple_ft_${name}.toml
  tenes_std -o input_ft_${name}.toml std_ft_${name}.toml
  tenes input_ft_${name}.toml

  awk '$2 == 0 {print $1, $3, $4}' output_ft_${name}/FT_density.dat > energy_${name}.dat
  awk '$2 == 1 {print $1, $3, $4}' output_ft_${name}/FT_density.dat > magnetization_${name}.dat
  awk '$2 == 2 {print $1, $3, $4}' output_ft_${name}/FT_density.dat > magnetization_x_${name}.dat
done

