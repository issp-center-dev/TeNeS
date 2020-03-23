from os.path import join

import numpy as np

import toml

num_g = 21

for idx in range(num_g):
    try:
        with open("simple_{}.toml".format(idx)) as f:
            dict_toml = toml.load(f)
        g = dict_toml["model"]["H"]

        with open(join("output_{}".format(idx), "density.dat")) as f:
            lines = f.readlines()
            ene = lines[2].split('=')[1].strip()
            mag_sz = lines[0].split('=')[1].strip()
            mag_sx = lines[1].split('=')[1].strip()
        values=[g, ene, mag_sz]
        print("{} {} {}".format(values[0], values[1], mag_sz))
    except:
        print("Error: idx = {}".format(idx))
        continue
