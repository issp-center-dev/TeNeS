from os.path import join

import numpy as np

import toml


num_hx = 16

for idx in range(num_hx):
    try:
        with open("simple_{}.toml".format(idx)) as f:
            dict_toml = toml.load(f)
        hx = dict_toml["model"]["hx"]
        ene = 0.0
        mag_sz = 0.0
        mag_sx = 0.0
        with open(join("output_{}".format(idx), "density.dat")) as f:
            for line in f:
                words = line.split()
                if words[0] == 'hamiltonian':
                    ene = words[2]
                elif words[0] == 'Sz':
                    mag_sz = words[2]
                elif words[0] == 'Sx':
                    mag_sx = words[2]
        print("{} {} {} {}".format(hx, ene, mag_sz, mag_sx))
    except:
        continue
