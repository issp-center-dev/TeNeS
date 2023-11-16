import sys
import numpy as np
from scipy.interpolate import splrep, splev


def spec(name):
    data = np.loadtxt(f"./energy_{name}.dat")
    bs = data[:, 0]
    Es = data[:, 1]
    spl = splrep(bs, Es)
    Cs = -bs * bs * splev(bs, spl, der=1)
    with open(f"./spec_{name}.dat", "w") as f:
        for b, C in zip(bs, Cs):
            f.write(f"{b} {C}\n")


if __name__ == "__main__":
    name = sys.argv[1]
    spec(name)
