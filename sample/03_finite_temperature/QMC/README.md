QMC result for TFI model on a square lattice (32x32) by using ALPS/looper (https://github.com/wistaria/alps-looper).

``` bash
alpspython init_params.py
mpiexec -np 4 loop --mpi params.in.xml
alpspython extract.py
```
