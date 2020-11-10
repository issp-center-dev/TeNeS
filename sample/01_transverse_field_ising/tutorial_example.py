import subprocess

import numpy as np
import toml

MPI_cmd = ""  # e.g., "mpiexec -np 1"

num_hx = 16
min_hx = 0.0
max_hx = 3.0

total = 0
for idx, hx in enumerate(np.linspace(min_hx, max_hx, num=num_hx)):
    print("Caclulation Process: {}/{}".format(idx + 1, num_hx))
    with open("simple.toml") as f:
        dict_toml = toml.load(f)
    dict_toml["parameter"]["general"]["output"] = "output_{}".format(idx)
    dict_toml["model"]["hx"] = float(hx)
    with open("simple_{}.toml".format(idx), "w") as f:
        toml.dump(dict_toml, f)
    cmd = "tenes_simple simple_{}.toml -o std_{}.toml".format(idx, idx)
    subprocess.call(cmd.split())
    cmd = "tenes_std std_{}.toml -o input_{}.toml".format(idx, idx)
    subprocess.call(cmd.split())
    cmd = "{} tenes input_{}.toml".format(MPI_cmd, idx)
    subprocess.call(cmd.split())
