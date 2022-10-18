import subprocess
from os.path import join
import numpy as np
import toml

MPI_cmd = ""  # e.g., "mpiexec -np 1"

num_h = 21
min_h = 0.0
max_h = 5.0
num_step_table = [100, 200, 500, 1000, 2000]

fmag = open("magnetization.dat", "w")
fene = open("energy.dat", "w")

for f in (fmag, fene):
    f.write("# $1: hz\n")
    for i, num_step in enumerate(num_step_table, 2):
        f.write(f"# ${i}: num_step={num_step}\n")
    f.write("\n")

for idx, h in enumerate(np.linspace(min_h, max_h, num=num_h)):
    print(f"Calculation Process: {idx+1}/{num_h}")
    inum = 0
    num_pre = 0
    fmag.write(f"{h} ")
    fene.write(f"{h} ")
    for num_step in num_step_table:
        ns = num_step - num_pre
        print(f"Steps: {num_step}")
        with open("basic.toml") as f:
            dict_toml = toml.load(f)

        output_dir = f"output_{idx}_{num_step}"

        dict_toml["parameter"]["general"]["output"] = output_dir
        dict_toml["parameter"]["general"]["tensor_save"] = "tensor_save"
        dict_toml["model"]["hz"] = float(h)
        dict_toml["parameter"]["simple_update"]["num_step"] = ns
        if inum > 0:
            dict_toml["parameter"]["general"]["tensor_load"] = "tensor_save"

        simple_toml = f"simple_{idx}_{num_step}.toml"
        std_toml = f"std_{idx}_{num_step}.toml"
        input_toml = f"input_{idx}_{num_step}.toml"

        with open(simple_toml, "w") as f:
            toml.dump(dict_toml, f)
        cmd = f"tenes_simple {simple_toml} -o {std_toml}"
        subprocess.call(cmd.split())

        cmd = f"tenes_std {std_toml} -o {input_toml}"
        subprocess.call(cmd.split())

        cmd = f"{MPI_cmd} tenes {input_toml}"
        subprocess.call(cmd.split())

        ene = 0.0
        mag_sz = 0.0
        with open(join(output_dir, "density.dat")) as f:
            for line in f:
                name, vals = line.split("=")
                if name.strip() == "hamiltonian":
                    re, im = vals.split()
                    ene += float(re)
                elif name.strip() == "Sz":
                    re, im = vals.split()
                    mag_sz += float(re)
        fene.write(f"{ene} ")
        fmag.write(f"{mag_sz} ")
        inum = inum + 1
        num_pre = num_step
    fene.write("\n")
    fmag.write("\n")
fene.close()
fmag.close()
