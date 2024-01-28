from __future__ import annotations
import subprocess
from typing import TextIO
from os.path import join
import numpy as np
import toml

MPI_cmd = ""  # mpiexec -np 4"

# calculate S(q) = \sum <n_i n_j> exp[-i q(r_i - r_j)] or not
# If True, only 3x3 unitcell is available
calculate_sq = True

# if False, skip computing and only gather results into *.dat (e.g., nq.dat)
do_calculation = True

# sweep mu from larger to smaller or not
decrease_mu = True

num_mu = 21
min_mu = -1.0
max_mu = 3.0
num_step_table = [1000, 2000]
D = 2
chi = 4

# Please uncomment the following lines if you hope to obtain more accurate results
# num_mu = 41
# num_step_table = [1500, 3000]
# D = 5
# chi = 25

# file objects for physical quantities
files: dict[str, TextIO] = {}

# N(q) = (1/Nsite) \sum_i <n_i> exp[-i q r_i] with q = (4pi/3, 0)
files["nq"] = open("nq.dat", "w")
if calculate_sq:
    # S(q) = (1/Nsite) \sum_{ij} <n_i n_j> exp[-i q (r_i-r_j)] with q = (4pi/3, 0)
    files["sq"] = open("sq.dat", "w")
# S(q)' = (1/Nsite) \sum_i \sum_dx <n_i n(r_i+dx)> exp[-i q dx] with q = (4pi/3, 0)
files["sq_x"] = open("sq_x.dat", "w")
# N0 = (<B> + <B^dagger>) / 2
files["den0"] = open("offdiag.dat", "w")
# N = (1/Nsite) \sum_i <n_i>
files["den"] = open("density.dat", "w")
# E = <H>
files["ene"] = open("energy.dat", "w")

for f in files.values():
    f.write("# $1: mu\n")
    for i, num_step in enumerate(num_step_table, 2):
        f.write(f"# ${i}: num_step={num_step}\n")
    f.write("\n")


step_mu = -1 if decrease_mu else 1

for idx, mu in enumerate(np.linspace(min_mu, max_mu, num=num_mu)[::step_mu]):
    print(f"Calculation Process: {idx+1}/{num_mu}")
    inum = 0
    num_pre = 0
    for f in files.values():
        f.write("{} ".format(mu))
    for num_step in num_step_table:
        ns = num_step - num_pre
        print(f"Steps: {num_step}")
        output_dir = f"output_{idx}_{num_step}"
        if do_calculation:
            with open("basic.toml") as f:
                dict_toml = toml.load(f)

            L = dict_toml["lattice"]["L"]
            W = dict_toml["lattice"]["W"]
            nsites = L * W
            if calculate_sq and nsites != 9:
                raise RuntimeError("The number of sites is not 9 (3x3)!")

            dict_toml["parameter"]["general"]["output"] = output_dir
            dict_toml["parameter"]["general"]["tensor_save"] = "tensor_save"
            dict_toml["model"]["mu"] = float(mu)
            dict_toml["lattice"]["virtual_dim"] = D
            dict_toml["parameter"]["ctm"]["dimension"] = chi
            dict_toml["parameter"]["simple_update"]["num_step"] = ns

            # load the previous optimized tensor
            if idx > 0 or inum > 0:
                dict_toml["parameter"]["general"]["tensor_load"] = "tensor_save"

            simple_toml = f"simple_{idx}_{num_step}.toml"
            std_toml = f"std_{idx}_{num_step}.toml"
            input_toml = f"input_{idx}_{num_step}.toml"

            with open(simple_toml, "w") as f:
                toml.dump(dict_toml, f)
            cmd = f"tenes_simple {simple_toml} -o {std_toml}"
            subprocess.call(cmd.split())

            cmd = f"tenes_std -o {input_toml} {std_toml}"
            if calculate_sq:
                # twosite observable <n_i n_j> are also calculated
                cmd += " nn_obs.toml"
            subprocess.call(cmd.split())

            cmd = f"{MPI_cmd} tenes {input_toml}"
            subprocess.call(cmd.split())

        with open(join(output_dir, "density.dat")) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                words = line.split("=")
                if words[0].strip() == "N":
                    density = float(words[1].strip().split()[0])
                    files["den"].write(f"{density} ")
                if words[0].strip() == "Energy":
                    ene = float(words[1].strip().split()[0])
                    files["ene"].write(f"{ene} ")

        coords_list = []
        with open("coordinates.dat") as f:
            for line in f:
                words = line.split()
                if len(words) == 0:
                    continue
                if words[0].startswith("#"):
                    continue
                x = float(words[1])
                y = float(words[2])
                coords_list.append([x, y])
        coords = np.array(coords_list)

        nsites = coords.shape[0]
        q = np.array([4 * np.pi / 3, 0])
        nq_coeffs = np.cos(coords @ q)

        n0 = 0.0  # (<B> + <B^dagger>) / 2
        nq = 0.0  # n(q), q = (4pi/a, 0)
        N = np.zeros(nsites)
        with open(join(output_dir, "onesite_obs.dat")) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                words = line.split()
                if len(words) == 0:
                    continue
                if words[0].strip() == "0":  # N
                    i = int(words[1])
                    n = float(words[2])
                    nq += nq_coeffs[i] * n
                    N[i] = n
                if words[0].strip() == "1":  # Bdagger
                    n = float(words[2])
                    n0 += 0.5 * np.abs(n)
                if words[0].strip() == "2":  # B
                    n = float(words[2])
                    n0 += 0.5 * np.abs(n)
        files["den0"].write(f"{n0 / nsites} ")
        files["nq"].write(f"{nq / nsites} ")

        if calculate_sq:
            sq = 0.0  # S(q), q = (4pi/a, 0)
            for i in range(nsites):
                sq += N[i]  # \hat{n}_i^2 = \hat{n}_i for hardcore boson
            with open(join(output_dir, "twosite_obs.dat")) as f:
                for line in f:
                    words = line.split()
                    if len(words) == 0:
                        continue
                    if words[0].strip() == "4":  # <nn>
                        dx = int(words[2])
                        dy = int(words[3])
                        nn = float(words[4])
                        if dx == dy:
                            sq += nn
                        else:
                            sq -= 0.5 * nn
            files["sq"].write(f"{sq / nsites} ")

        sq = 0.0
        for i in range(nsites):
            # C(r=0)
            # \hat{n}_i^2 = \hat{n}_i for hardcore boson
            sq += N[i]
        with open(join(output_dir, "correlation.dat")) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                words = line.split()
                if len(words) == 0:
                    continue
                dx = int(words[3])
                dy = int(words[4])
                nn = float(words[5])

                # Consider only the C(r) along x-axis
                if dy > 0:
                    continue
                if dx % 3 == 0:
                    sq += nn
                else:
                    sq += -0.5 * nn
        files["sq_x"].write(f"{sq / nsites} ")

        inum = inum + 1
        num_pre = num_step
        for f in files.values():
            f.flush()
    for f in files.values():
        f.write("\n")
for f in files.values():
    f.close()
