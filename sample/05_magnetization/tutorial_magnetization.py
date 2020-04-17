import subprocess
from os.path import join
import numpy as np
import toml

num_h = 21
min_h = 0.0
max_h = 5.0
num_step_table = [100, 200, 500, 1000, 2000]

fmag = open("magnetization.dat","w")
fene = open("energy.dat","w")
for idx, h in enumerate(np.linspace(min_h, max_h, num=num_h)):
    print("Caclulation Process: {}/{}".format(idx+1, num_h))
    inum = 0
    num_pre = 0
    fmag.write("{} ".format(h))
    fene.write("{} ".format(h))
    for num_step in num_step_table:
        ns = num_step - num_pre
        print("Step number: {}".format(num_step))
        with open("basic.toml") as f:
            dict_toml = toml.load(f)
        dict_toml["parameter"]["general"]["output"] = "output_{}_{}".format(idx,num_step)
        dict_toml["parameter"]["general"]["tensor_save"] = "tensor_save".format(idx,num_step)
        dict_toml["model"]["H"] = float(h)
        dict_toml["parameter"]["simple_update"]["num_step"] = ns
        if inum > 0:
            dict_toml["parameter"]["general"]["tensor_load"] = "tensor_save".format(idx,num_pre)
        with open("simple_{}_{}.toml".format(idx,num_step), 'w') as f:
            toml.dump(dict_toml, f)
        cmd = "tenes_simple simple_{}_{}.toml -o std_{}_{}.toml".format(idx,num_step,idx,num_step)
        subprocess.call(cmd.split())
        cmd = "tenes_std std_{}_{}.toml -o input_{}_{}.toml".format(idx,num_step,idx,num_step)
        subprocess.call(cmd.split())
        cmd = "tenes input_{}_{}.toml".format(idx,num_step)
        subprocess.call(cmd.split())
        with open(join("output_{}_{}".format(idx,num_step), "density.dat")) as f:
            lines = f.readlines()
            mag_sz = lines[0].split('=')[1].strip()
            ene = lines[2].split('=')[1].strip()
        fene.write("{} ".format(ene))
        fmag.write("{} ".format(mag_sz))
        inum = inum + 1
        num_pre = num_step
    fene.write("\n")
    fmag.write("\n")
fene.close()
fmag.close()
