import subprocess

import numpy as np
import toml

MPI_cmd = ""  # e.g., "mpiexec -np 1"

num_hx = 16
min_hx = 0.0
max_hx = 3.0

total = 0
for idx, hx in enumerate(np.linspace(min_hx, max_hx, num=num_hx)):
    print(f"Calculation Process: {idx+1}/{num_hx}")
    with open("simple.toml") as f:
        dict_toml = toml.load(f)
    dict_toml["parameter"]["general"]["output"] = f"output_{idx}"
    dict_toml["model"]["hx"] = float(hx)

    simple_toml = f"simple_{idx}.toml"
    std_toml = f"std_{idx}.toml"
    input_toml = f"input_{idx}.toml"

    with open(simple_toml, "w") as f:
        toml.dump(dict_toml, f)
    cmd = f"tenes_simple {simple_toml} -o {std_toml}"
    subprocess.call(cmd.split())

    cmd = f"tenes_std {std_toml} -o {input_toml}"
    subprocess.call(cmd.split())

    cmd = f"{MPI_cmd} tenes {input_toml}"
    subprocess.call(cmd.split())
