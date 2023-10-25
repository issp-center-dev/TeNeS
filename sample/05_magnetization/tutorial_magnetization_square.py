import subprocess
from os.path import join
import toml

MPI_cmd = ""  # e.g., "mpiexec -np 1"

h_table = [
    0.0,
    0.25,
    0.5,
    0.75,
    1.0,
    1.25,
    1.5,
    1.75,
    2.0,
    2.25,
    2.5,
    2.75,
    3.0,
    3.1,
    3.2,
    3.3,
    3.4,
    3.5,
    3.6,
    3.7,
    3.8,
    3.9,
    4.0,
    4.25,
    4.5,
    4.75,
    5.0,
]
num_h = len(h_table)
num_step_table = [100, 200, 500, 1000, 2000]

fene = open("energy_square.dat", "w")
fmag = open("magnetization_square.dat", "w")

for f in (fmag, fene):
    f.write("# $1: hz\n")
    for i, num_step in enumerate(num_step_table, 2):
        f.write(f"# ${i}: num_step={num_step}\n")
    f.write("\n")

for idx in range(num_h):
    h = h_table[idx]
    print(f"Calculation Process: {idx+1}/{num_h}")
    inum = 0
    num_pre = 0
    fene.write(f"{h} ")
    fmag.write(f"{h} ")
    for num_step in num_step_table:
        ns = num_step - num_pre
        print(f"Steps: {num_step}")
        with open("basic_square.toml") as f:
            dict_toml = toml.load(f)

        output_dir = f"output_square_{idx}_{num_step}"
        dict_toml["parameter"]["general"]["output"] = output_dir
        dict_toml["parameter"]["general"]["tensor_save"] = "tensor_save_square"
        dict_toml["model"]["hz"] = float(h)
        dict_toml["parameter"]["simple_update"]["num_step"] = ns
        if inum > 0:
            dict_toml["parameter"]["general"]["tensor_load"] = "tensor_save_square"

        simple_toml = f"simple_square_{idx}_{num_step}.toml"
        std_toml = f"std_square_{idx}_{num_step}.toml"
        input_toml = f"input_square_{idx}_{num_step}.toml"

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
                if name.strip() == "Energy":
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
