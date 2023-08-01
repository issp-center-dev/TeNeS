for name in gapless gapfull; do
  HPhi -s stan_${name}.in
  mv output output_${name}
  python3 ./postprocess.py $name
done
