import subprocess

import toml

total = 20
for idx in range(total + 1):
    print("Caclulation Process: {}/{}".format(idx, total))
    dict_toml = toml.load(open("simple.toml"))
    dict_toml["model"]["G"] = 0.2 * idx
    toml.dump(dict_toml, open("simple_{}.toml".format(idx), mode="w"))
    cmd = "tenes_simple simple_{}.toml -o input_{}.toml".format(idx, idx)
    subprocess.call(cmd.split())
    cmd = "tenes input_{}.toml".format(idx)
    subprocess.call(cmd.split())
    cmd = "mv output output_{}".format(idx)
    subprocess.call(cmd.split())
