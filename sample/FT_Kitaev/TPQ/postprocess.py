import sys
import os
import numpy as np

def foo(name):
    with open(f"stan_{name}.in") as f:
        nrands = -1
        nsteps = -1
        for line in f:
            words = line.split('=')
            if(words[0].strip() == 'NumAve'):
                nrands = int(words[1].strip())
            if(words[0].strip() == 'Lanczos_max'):
                nsteps = int(words[1].strip())

    bs = np.zeros((nsteps, nrands))
    Es = np.zeros((nsteps, nrands))
    E2s = np.zeros((nsteps, nrands))
    N = 0
    for irand in range(nrands):
        A = np.loadtxt(os.path.join(f"output_{name}", f"SS_rand{irand}.dat"))
        bs[:, irand] = A[:, 0]
        Es[:, irand] = A[:, 1]
        E2s[:, irand] = A[:, 2]
        N = A[0,4]
    bs_jk = (np.sum(bs, axis=1).reshape(-1,1) - bs) / (nrands - 1)
    Es_jk = (np.sum(Es, axis=1).reshape(-1,1) - Es) / (nrands - 1)
    E2s_jk = (np.sum(E2s, axis=1).reshape(-1,1) - E2s) / (nrands - 1)
    C_jk = (E2s_jk - Es_jk**2) * (bs_jk**2)

    bs_mean = np.mean(bs_jk, axis=1)
    Es_mean = np.mean(Es_jk, axis=1) / N
    Es_err = np.sqrt((nrands - 1) * Es_jk.var(axis=1, ddof=0)) / N
    C_mean = np.mean(C_jk, axis=1) / N
    C_err = np.sqrt((nrands - 1) * C_jk.var(axis=1, ddof=0)) / N

    with open(f"TPQ_{name}.dat", "w") as f:
        for i in range(nsteps):
            f.write(f"{bs_mean[i]} {Es_mean[i]} {Es_err[i]} {C_mean[i]} {C_err[i]}\n")

if __name__ == '__main__':
    foo(sys.argv[1])
