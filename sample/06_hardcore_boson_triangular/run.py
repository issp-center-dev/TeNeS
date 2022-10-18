import subprocess
from os.path import join
import numpy as np
import toml

MPI_cmd = ""  # mpiexec -np 4"

do_calculation = True
decrease_mu = True # sweep mu from larger to smaller or not

num_h = 21
min_h = -1.0
max_h = 3.0
num_step_table = [1000, 2000]
D = 2
chi = 4
# Please uncomment the following lines if you hope to obtain more accurate results
# num_h = 41
# num_step_table = [1500, 3000]
# D = 5
# chi = 25

fsq = open("sq.dat", "w")
fmag0 = open("offdiag.dat", "w")
fmag = open("density.dat", "w")
fene = open("energy.dat", "w")

for f in (fsq, fmag, fmag0, fene):
    f.write("# $1: mu\n")
    for i, num_step in enumerate(num_step_table, 2):
        f.write(f"# ${i}: num_step={num_step}\n")
    f.write("\n")


step_mu = -1 if decrease_mu else 1

for idx, h in enumerate(np.linspace(min_h, max_h, num=num_h)[::step_mu]):
    print(f"Calculation Process: {idx+1}/{num_h}")
    inum = 0
    num_pre = 0
    fsq.write("{} ".format(h))
    fmag0.write("{} ".format(h))
    fmag.write("{} ".format(h))
    fene.write("{} ".format(h))
    for num_step in num_step_table:
        ns = num_step - num_pre
        print(f"Steps: {num_step}")
        output_dir = f"output_{idx}_{num_step}"
        if do_calculation:
            with open("basic.toml") as f:
                dict_toml = toml.load(f)
            dict_toml["parameter"]["general"]["output"] = output_dir
            dict_toml["parameter"]["general"]["tensor_save"] = "tensor_save"
            dict_toml["model"]["mu"] = float(h)
            dict_toml["lattice"]["virtual_dim"] = D
            dict_toml["parameter"]["ctm"]["dimension"] = chi
            dict_toml["parameter"]["simple_update"]["num_step"] = ns
            if idx > 0 or inum > 0:
                dict_toml["parameter"]["general"]["tensor_load"] = "tensor_save"

            simple_toml = f"simple_{idx}_{num_step}.toml"
            std_toml = f"std_{idx}_{num_step}.toml"
            input_toml = f"input_{idx}_{num_step}.toml"

            with open(simple_toml, "w") as f:
                toml.dump(dict_toml, f)
            cmd = f"tenes_simple {simple_toml} -o {std_toml}"
            subprocess.call(cmd.split())

            with open(std_toml, "a") as fout:
                with open("nn_obs.toml") as fin:
                    for line in fin:
                        fout.write(line)
            cmd = f"tenes_std {std_toml} -o {input_toml}"
            subprocess.call(cmd.split())

            cmd = f"{MPI_cmd} tenes {input_toml}"
            subprocess.call(cmd.split())
        ene = ""
        density = ""
        with open(join(output_dir, "density.dat")) as f:
            for line in f:
                words = line.split("=")
                if words[0].strip() == "N":
                    density = words[1].strip()
                if words[0].strip() == "hamiltonian":
                    ene = words[1].strip()
        fene.write(f"{ene} ")
        fmag.write(f"{density} ")
        fene.flush()
        fmag.flush()

        n0 = 0.0
        N = np.zeros(9)
        with open(join(output_dir, "onesite_obs.dat")) as f:
            for line in f:
                words = line.split()
                if len(words) == 0:
                    continue
                if words[0].strip() == "0":
                    n = float(words[2])
                    N[int(words[1])] = n
                if words[0].strip() == "1":
                    n = float(words[2])
                    n0 += np.abs(n) / 18.0
                if words[0].strip() == "2":
                    n = float(words[2])
                    n0 += np.abs(n) / 18.0
        fmag0.write("{} ".format(n0))
        fmag0.flush()

        sq = 0.0
        for i in range(3):
            sq += N[i]

        with open(join(output_dir, "twosite_obs.dat")) as f:
            for line in f:
                words = line.split()
                if len(words) == 0:
                    continue
                if words[0].strip() == "4":
                    n = float(words[4])
                    if words[2] != words[3]:
                        n *= -0.5
                    sq += n
        sq *= 3.0 / 9 / 9
        fsq.write(f"{sq} ")
        fsq.flush()

        inum = inum + 1
        num_pre = num_step
    fene.write("\n")
    fmag.write("\n")
    fmag0.write("\n")
    fsq.write("\n")
fene.close()
fmag.close()
fmag0.close()
fsq.close()
