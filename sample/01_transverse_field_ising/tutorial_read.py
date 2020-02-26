from os.path import join

import numpy as np


num_g = 21
min_g = 0.0
max_g = 2.0

for idx, g in enumerate(np.linspace(min_g, max_g, num=num_g)):
    try:
        with open(join("output_{}".format(idx), "energy.dat")) as f:
            lines = f.readlines()
            ene = lines[0].split('=')[1].strip()
            mag_sz = lines[1].split('=')[1].strip()
            mag_sx = lines[2].split('=')[1].strip()
        print("{} {} {} {}".format(g, ene, mag_sz, mag_sx))
    except:
        continue
