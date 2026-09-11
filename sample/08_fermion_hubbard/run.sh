for D in 2 3 4 5; do
  OMP_NUM_THREADS=1 tenes ./input-opt-D${D}.toml
  OMP_NUM_THREADS=1 tenes ./input-measure-D${D}.toml
done
