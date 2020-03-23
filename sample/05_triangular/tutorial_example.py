import subprocess

import numpy as np

import toml


num_h = 21
min_h = 0.0
max_h = 5.0

total = 0
for idx, h in enumerate(np.linspace(min_h, max_h, num=num_h)):
    print("Caclulation Process: {}/{}".format(idx+1, num_h))
    with open("basic.toml") as f:
        dict_toml = toml.load(f)
    dict_toml["parameter"]["general"]["output"] = "output_{}".format(idx)
    dict_toml["parameter"]["general"]["tensor_save"] = "tensor_save_{}".format(idx)
    dict_toml["model"]["H"] = float(h)
    with open("simple_{}.toml".format(idx), 'w') as f:
        toml.dump(dict_toml, f)
    cmd = "tenes_simple simple_{}.toml -o std_{}.toml".format(idx, idx)
    subprocess.call(cmd.split())
    cmd = "tenes_std std_{}.toml -o input_{}.toml".format(idx, idx)
    subprocess.call(cmd.split())
    cmd = "tenes input_{}.toml".format(idx)
    subprocess.call(cmd.split())
