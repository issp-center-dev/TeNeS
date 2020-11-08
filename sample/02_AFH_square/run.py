import subprocess
from sys import argv
from os.path import join
import numpy as np
import toml

MPI_cmd = ""  # e.g., "mpiexec -np 1"

Ds = [2, 3]
num_step_table = [0, 10]

# Note that D >= 5 takes a few hours (maybe the order of days in a single core)
# Ds = [2, 3, 4, 5, 6]
# num_step_table = [0, 10, 100, 1000]

fmag = open("magnetization.dat", "w")
fene = open("energy.dat", "w")
for D in Ds:
    print("Caclulation Process: D = {}".format(D))
    chi = D * D
    num_pre = 0
    fmag.write("{} ".format(D))
    fene.write("{} ".format(D))
    for num_step in num_step_table:
        ns = num_step - num_pre
        print("Steps: {}".format(num_step))
        with open("basic.toml") as f:
            dict_toml = toml.load(f)
        dict_toml["parameter"]["general"]["output"] = "output_{}_{}".format(D, num_step)
        dict_toml["parameter"]["general"]["tensor_save"] = "tensor_save".format(
            D, num_step
        )
        dict_toml["lattice"]["virtual_dim"] = D
        dict_toml["parameter"]["ctm"]["dimension"] = chi
        dict_toml["parameter"]["simple_update"]["num_step"] = 1000
        dict_toml["parameter"]["full_update"]["num_step"] = ns
        if ns > 0:
            dict_toml["parameter"]["simple_update"]["num_step"] = 0
            dict_toml["parameter"]["general"]["tensor_load"] = "tensor_save".format(
                D, num_pre
            )
        with open("simple_{}_{}.toml".format(D, num_step), "w") as f:
            toml.dump(dict_toml, f)
        cmd = "tenes_simple simple_{}_{}.toml -o std_{}_{}.toml".format(
            D, num_step, D, num_step
        )
        subprocess.call(cmd.split())
        cmd = "tenes_std std_{}_{}.toml -o input_{}_{}.toml".format(
            D, num_step, D, num_step
        )
        subprocess.call(cmd.split())
        cmd = "{} tenes input_{}_{}.toml".format(MPI_cmd, D, num_step)
        subprocess.call(cmd.split())
        with open(join("output_{}_{}".format(D, num_step), "density.dat")) as f:
            for line in f:
                words = line.split()
                if words[0].strip() == "hamiltonian":
                    ene = float(words[2].strip())
        with open(join("output_{}_{}".format(D, num_step), "onesite_obs.dat")) as f:
            for line in f:
                words = line.split()
                if len(words) < 2:
                    continue
                if words[0] == words[1] == "0":
                    mag = float(words[2])
        fene.write("{} ".format(ene))
        fmag.write("{} ".format(mag))
        fene.flush()
        fmag.flush()
        num_pre = num_step
    fene.write("\n")
    fmag.write("\n")
fene.close()
fmag.close()
