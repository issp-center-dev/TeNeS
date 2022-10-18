import subprocess
from sys import argv
from os.path import join
import numpy as np
import toml

MPI_cmd = ""  # e.g., "mpiexec -np 1"

Ds = [2, 3, 4]
num_step_table = [0, 10]

# Note that D >= 5 takes a few hours (maybe the order of days in a single core)
# Ds = [2, 3, 4, 5, 6]
# num_step_table = [0, 10, 100, 1000]

fmag = open("magnetization.dat", "w")
fene = open("energy.dat", "w")

for f in (fmag, fene):
    f.write("# $1: D\n")
    for i, num_step in enumerate(num_step_table, 2):
        f.write(f"# ${i}: {num_step=}\n")
    f.write("\n")

for D in Ds:
    print(f"Calculation Process: D = {D}")
    chi = D * D
    num_pre = 0
    fmag.write(f"{D} ")
    fene.write(f"{D} ")
    for num_step in num_step_table:
        ns = num_step - num_pre
        print(f"Steps: {num_step}")

        simple_toml = f"simple_{D}_{num_step}.toml"
        std_toml = f"std_{D}_{num_step}.toml"
        input_toml = f"input_{D}_{num_step}.toml"
        output_dir = f"output_{D}_{num_step}"

        with open("basic.toml") as f:
            dict_toml = toml.load(f)
        dict_toml["parameter"]["general"]["output"] = output_dir
        dict_toml["parameter"]["general"]["tensor_save"] = "tensor_save"
        dict_toml["lattice"]["virtual_dim"] = D
        dict_toml["parameter"]["ctm"]["dimension"] = chi
        dict_toml["parameter"]["simple_update"]["num_step"] = 1000
        dict_toml["parameter"]["full_update"]["num_step"] = ns
        if ns > 0:
            dict_toml["parameter"]["simple_update"]["num_step"] = 0
            dict_toml["parameter"]["general"]["tensor_load"] = "tensor_save"

        with open(simple_toml, "w") as f:
            toml.dump(dict_toml, f)
        cmd = f"tenes_simple {simple_toml} -o {std_toml}"
        subprocess.call(cmd.split())

        cmd = f"tenes_std {std_toml} -o {input_toml}"
        subprocess.call(cmd.split())

        cmd = f"{MPI_cmd} tenes {input_toml}"
        subprocess.call(cmd.split())

        with open(join(output_dir, "density.dat")) as f:
            for line in f:
                words = line.split()
                if words[0].strip() == "hamiltonian":
                    ene = float(words[2].strip())
        with open(join(output_dir, "onesite_obs.dat")) as f:
            for line in f:
                words = line.split()
                if len(words) < 2:
                    continue
                if words[0] == words[1] == "0":
                    mag = float(words[2])
        fene.write(f"{ene} ")
        fmag.write(f"{mag} ")
        fene.flush()
        fmag.flush()
        num_pre = num_step
    fene.write("\n")
    fmag.write("\n")
fene.close()
fmag.close()
