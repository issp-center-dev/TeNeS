from os.path import join

import toml

num_hx = 16

print("# $1: h")
print("# $2: ene")
print("# $3: sz")
print("# $4: sx")
print()

for idx in range(num_hx):
    try:
        with open(f"simple_{idx}.toml") as f:
            dict_toml = toml.load(f)
        hx = dict_toml["model"]["hx"]
        ene = 0.0
        mag_sz = 0.0
        mag_sx = 0.0
        with open(join(f"output_{idx}", "density.dat")) as f:
            for line in f:
                words = line.split()
                if words[0] == "Energy":
                    ene = words[2]
                elif words[0] == "Sz":
                    mag_sz = words[2]
                elif words[0] == "Sx":
                    mag_sx = words[2]
        print(f"{hx} {ene} {mag_sz} {mag_sx}")
    except:
        continue
