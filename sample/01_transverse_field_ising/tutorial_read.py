import os

for idx in range(21):
    D = idx * 0.2
    try:
        with open("output_{}/energy.dat".format(idx), "r") as f:
            lines = f.readlines()
            for line in lines:
                ene = line.split()
        with open("output_{}/site_obs.dat".format(idx), "r") as f:
            lines = f.readlines()
            mag_sz = lines[5].split()
            mag_sx = lines[9].split()

        print("{} {} {} {}".format(D, ene[0], mag_sz[2], mag_sx[2]))

    except:
        continue
