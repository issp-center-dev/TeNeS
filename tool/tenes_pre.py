from __future__ import print_function

from itertools import chain
import numpy as np
import numpy.linalg as linalg
import toml


def make_evolution_tensor(Ham, tau):
    eigval, eigvec = linalg.eig(Ham)
    return np.dot(np.dot(eigvec, np.diag(np.exp(-tau * eigval))), eigvec.transpose())


def array_to_str(A):
    res = []
    X, Y = A.shape
    for x in range(X):
        row = []
        for y in range(Y):
            row.append(str(A[x, y]))
        res.append(" ".join(row))
    return "\n".join(res)


def gen_spinhamiltonian(param):
    S = param["S"]
    if int(2 * S) != 2 * S:
        msg = "S is neighther integer nor half-integer: {}".format(S)
        raise RuntimeError(msg)
    M = int(2 * S) + 1
    E = np.eye(M)
    Sz = np.zeros((M, M))
    Splus = np.zeros((M, M))
    Sminus = np.zeros((M, M))
    for i in range(M):
        m = S - i
        Sz[i, i] = m
        if i > 0:
            Splus[i, i - 1] = np.sqrt((S - m) * (S + m + 1.0))
            Sminus[i - 1, i] = np.sqrt((S + m) * (S - m + 1.0))
    Sx = 0.5 * (Splus + Sminus)

    z = 4.0
    Jz = param.get("Jz", 1.0)
    Jxy = param.get("Jxy", Jz)
    h = param.get("h", 0.0) / z
    G = param.get("Gamma", 0.0) / z
    D = param.get("D", 0.0) / z

    ham = Jz * np.kron(Sz, Sz)
    ham += 0.5 * Jxy * (np.kron(Splus, Sminus) + np.kron(Sminus, Splus))
    ham -= h * (np.kron(Sz, E) + np.kron(E, Sz))
    ham -= G * (np.kron(Sx, E) + np.kron(E, Sx))
    ham += D * (np.kron(np.dot(Sz, Sz), E) + np.kron(E, np.dot(Sz, Sz)))

    return [Sz, Sx], [ham]


def tenes_pre(param):
    if param["model"]["type"] == "spin":
        lops, hams = gen_spinhamiltonian(param["model"])
    else:
        msg = "Unknown model type: {}".format(param["model"]["type"])
        raise RuntimeError(msg)

    pparam = param.get("parameter", {})
    simple_update = pparam.get("simple_update", {})
    simple_tau = simple_update.get("tau", 0.01)
    simple_evols = [make_evolution_tensor(H, simple_tau) for H in hams]
    full_update = pparam.get("full_update", {})
    full_tau = full_update.get("tau", 0.01)
    full_evols = [make_evolution_tensor(H, full_tau) for H in hams]

    Lsub = param["lattice"]["L_sub"]
    if not (
        isinstance(Lsub, list)
        and len(Lsub) == 2
        and isinstance(Lsub[0], int)
        and isinstance(Lsub[1], int)
    ):
        msg = "Lsub is not a list with two integers"
        raise RuntimeError(msg)

    simpleupdate_str = []
    fullupdate_str = []
    X, Y = Lsub
    N = X * Y
    lattice = np.reshape(range(N), (X, Y))
    for leftx in range(0, X, 2):
        for lefty in range(Y):
            rightx = (leftx + 1) % X
            righty = lefty
            left = lattice[leftx, lefty]
            right = lattice[rightx, righty]
            simpleupdate_str.append("{} {} h 0".format(left, right))
            fullupdate_str.append("{} {} h 1".format(left, right))
    for leftx in range(1, X, 2):
        for lefty in range(Y):
            rightx = (leftx + 1) % X
            righty = lefty
            left = lattice[leftx, lefty]
            right = lattice[rightx, righty]
            simpleupdate_str.append("{} {} h 0".format(left, right))
            fullupdate_str.append("{} {} h 1".format(left, right))
    for bottomy in range(0, Y, 2):
        for bottomx in range(X):
            topx = bottomx
            topy = (bottomy + 1) % Y
            bottom = lattice[bottomx, bottomy]
            top = lattice[topx, topy]
            simpleupdate_str.append("{} {} v 0".format(bottom, top))
            fullupdate_str.append("{} {} v 1".format(bottom, top))
    for bottomy in range(1, Y, 2):
        for bottomx in range(X):
            topx = bottomx
            topy = (bottomy + 1) % Y
            bottom = lattice[bottomx, bottomy]
            top = lattice[topx, topy]
            simpleupdate_str.append("{} {} v 0".format(bottom, top))
            fullupdate_str.append("{} {} v 1".format(bottom, top))

    dict_evolution = {
        "simple_update": "\n".join(simpleupdate_str),
        "full_update": "\n".join(fullupdate_str),
        "matrix": list(map(array_to_str, chain(simple_evols, full_evols))),
    }

    dict_observable = {
        "local_operator": list(map(array_to_str, lops)),
        "hamiltonian": list(map(array_to_str, hams)),
        "hamiltonian_bonds": "\n".join(simpleupdate_str),
    }

    if "correlation" in param:
        corparam = param["correlation"]
        dict_correlation = {}
        dict_correlation["r_max"] = corparam["r_max"]
        if "operators" in corparam:
            dict_correlation["operators"] = corparam["operators"]
        else:
            dict_correlation["operators"] = [[0, 0], [0, 1], [1, 1]]

    res = {}
    if "parameter" in param:
        res["parameter"] = param["parameter"]
    res["lattice"] = param["lattice"]
    res["evolution"] = dict_evolution
    res["observable"] = dict_observable
    if "correlation" in param:
        res["correlation"] = dict_correlation

    return res


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate input files for TeNeS", add_help=True
    )

    parser.add_argument(
        "input",
        help="Input TOML file",
    )

    parser.add_argument(
        "-o", "--output", dest="output", default="input.toml", help="Output TOML file"
    )

    args = parser.parse_args()
    if args.input == args.output:
        print("The names of input and output are the same")
        sys.exit(1)
    res = tenes_pre(toml.load(args.input))

    with open(args.output, 'w') as f:
        toml.dump(res, f)
