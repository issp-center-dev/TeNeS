import subprocess
from os.path import join
import numpy as np
import toml

h_table = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5,
           3.6, 3.7, 3.8, 3.9, 4.0, 4.25, 4.5, 4.75, 5.0]
num_h = len(h_table)
num_step_table = [100, 200, 500, 1000, 2000]

fene = open("energy_square.dat","w")
fmag = open("magnetization_square.dat","w")
for idx in range(num_h):
    h = h_table[idx]
    print("Caclulation Process: {}/{}".format(idx+1, num_h))
    inum = 0
    num_pre = 0
    fene.write("{} ".format(h))
    fmag.write("{} ".format(h))
    for num_step in num_step_table:
        ns = num_step - num_pre
        print("Steps: {}".format(num_step))
        with open("basic_square.toml") as f:
            dict_toml = toml.load(f)
        dict_toml["parameter"]["general"]["output"] = "output_square_{}_{}".format(idx,num_step)
        dict_toml["parameter"]["general"]["tensor_save"] = "tensor_save_square"
        dict_toml["model"]["H"] = float(h)
        dict_toml["parameter"]["simple_update"]["num_step"] = ns
        if inum > 0:
            dict_toml["parameter"]["general"]["tensor_load"] = "tensor_save_square"
        with open("simple_square_{}_{}.toml".format(idx,num_step), 'w') as f:
            toml.dump(dict_toml, f)
        cmd = "tenes_simple simple_square_{}_{}.toml -o std_square_{}_{}.toml".format(idx,num_step,idx,num_step)
        subprocess.call(cmd.split())
        cmd = "tenes_std std_square_{}_{}.toml -o input_square_{}_{}.toml".format(idx,num_step,idx,num_step)
        subprocess.call(cmd.split())
        cmd = "tenes input_square_{}_{}.toml".format(idx,num_step)
        subprocess.call(cmd.split())

        with open(join("output_square_{}_{}".format(idx,num_step), "density.dat")) as f:
            lines = f.readlines()
            ene = lines[2].split('=')[1].strip()
            mag_sz = lines[0].split('=')[1].strip()
        fene.write("{} ".format(ene))
        fmag.write("{} ".format(mag_sz))
        inum = inum + 1
        num_pre = num_step
    fene.write("\n")
    fmag.write("\n")
fene.close()
fmag.close()
