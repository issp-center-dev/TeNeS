from os.path import join

import numpy as np

import toml


num_g = 16

for idx, in range(num_g):
    try:
        with open("simple.toml") as f:
            dict_toml = toml.load(f)
        g = dict_toml["model"]["G"]
        with open(join("output_{}".format(idx), "energy.dat")) as f:
            lines = f.readlines()
            ene = lines[0].split('=')[1].strip()
            mag_sz = lines[1].split('=')[1].strip()
            mag_sx = lines[2].split('=')[1].strip()
        print("{} {} {} {}".format(g, ene, mag_sz, mag_sx))
    except:
        continue
