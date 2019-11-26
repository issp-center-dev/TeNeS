from __future__ import print_function

from collections import namedtuple
from itertools import product

import numpy as np
import numpy.linalg as linalg

import toml

TeNeSInput = namedtuple("TeNeSInput", "param lattice model lop ham")


# https://stackoverflow.com/a/42913743
def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)


def make_evolution_tensor(Ham, tau):
    eigval, eigvec = linalg.eigh(Ham)
    return np.dot(np.dot(eigvec, np.diag(np.exp(-tau * eigval))), eigvec.transpose())


def array_to_strarr(A):
    res = []
    X, Y = A.shape
    for x in range(X):
        row = []
        for y in range(Y):
            row.append(str(A[x, y]))
        res.append(" ".join(row))
    return res


def index2coord(index, X):
    return index % X, index // X


def coord2index(x, y, X):
    return x + y * X


def getparam(param, name, index, default):
    P = param.get(name, default)
    if isinstance(P, list):
        if index > len(P):
            msg = "parameter {} has too few elements.".format(name)
            raise RuntimeError(msg)
        return P[index]
    return P


def make_Lsub(Lsub):
    if not isinstance(Lsub, list):
        return [Lsub, Lsub]
    else:
        if not (
            len(Lsub) == 2 and isinstance(Lsub[0], int) and isinstance(Lsub[1], int)
        ):
            msg = "Lsub is neither an integer nor a list with two integers"
            raise RuntimeError(msg)
        return Lsub


Bond = namedtuple("Bond", "source target dir type")


def dumpbond(bond):
    return "{} {} {} {}".format(bond.source, bond.target, bond.dir, bond.type)


class Lattice(object):
    def __init__(self):
        self.type = ""
        self.Lsub = [0, 0]
        self.z = 0
        self.bondtypes = 0
        self.bonds = []


class SquareLattice(Lattice):
    def __init__(self, param=None):
        super(SquareLattice, self).__init__()
        self.type = "square lattice"
        self.Lsub = [0, 0]
        self.z = 4
        self.bondtypes = 2
        if param is not None:
            self.load(param)

    def load(self, param):
        self.Lsub = make_Lsub(param["L_sub"])
        X, Y = self.Lsub

        self.bonds = []
        for leftx in range(0, X):
            for lefty in range(0, Y):
                rightx = (leftx + 1) % X
                righty = lefty
                left = coord2index(leftx, lefty, X)
                right = coord2index(rightx, righty, X)
                self.bonds.append(Bond(left, right, "h", 0))

        self.vertical_bonds = []
        for bottomy in range(0, Y):
            for bottomx in range(X):
                topx = bottomx
                topy = (bottomy + 1) % Y
                bottom = coord2index(bottomx, bottomy, X)
                top = coord2index(topx, topy, X)
                self.bonds.append(Bond(bottom, top, "v", 1))


class HoneycombLattice(Lattice):
    def __init__(self, param=None):
        super(HoneycombLattice, self).__init__()
        self.type = "honeycomb lattice"
        self.z = 3
        self.bondtypes = 3
        if param is not None:
            self.load(param)

    def load(self, param):
        self.Lsub = make_Lsub(param["L_sub"])
        if not all(map(lambda x: x % 2 == 0, self.Lsub)):
            msg = "All elements of Lsub must be even for a honeycomb lattice."
            raise RuntimeError(msg)
        X, Y = self.Lsub

        sublat_a = []
        sublat_b = []
        for x, y in product(range(X), range(Y)):
            index = x + y * X
            if (x + y) % 2 == 0:
                sublat_a.append(index)
            else:
                sublat_b.append(index)

        self.bonds = []

        for index in sublat_a:
            x, y = index2coord(index, X)
            rightx = (x + 1) % X
            rightindex = coord2index(rightx, y, X)
            self.bonds.append(Bond(index, rightindex, "h", 0))  # x-bond

        for index in sublat_b:
            x, y = index2coord(index, X)
            rightx = (x + 1) % X
            rightindex = coord2index(rightx, y, X)
            self.bonds.append(Bond(index, rightindex, "h", 1))  # y-bond

        for index in sublat_a:
            x, y = index2coord(index, X)
            topy = (y + 1) % Y
            topindex = coord2index(x, topy, X)
            self.bonds.append(Bond(index, topindex, "v", 2))  # z-bond


class SpinModel:
    def __init__(self):
        pass

    def localoperators(self, param):
        """
        Generates spin operators, Sz and Sx

        Parameters
        ----------
        param : dict
            parameter

        Returns
        -------
        Sz : darray
            Sz operator
        Sx : darray
            Sx operator
        Splus : darray
            Splus operator
        Sminus : darray
            Sminus operator

        """

        S = param.get("S", 0.5)

        if int(2 * S) != 2 * S:
            msg = "S is neighther integer nor half-integer: {}".format(S)
            raise RuntimeError(msg)
        M = int(2 * S) + 1
        Sz = np.zeros((M, M))
        Splus = np.zeros((M, M))
        Sminus = np.zeros((M, M))
        for i in range(M):
            m = S - i
            Sz[i, i] = m
            if i > 0:
                Splus[i - 1, i] = np.sqrt((S - m) * (S + m + 1.0))
            if i < M - 1:
                Sminus[i + 1, i] = np.sqrt((S + m) * (S - m + 1.0))
        Sx = 0.5 * (Splus + Sminus)

        return Sz, Sx, Splus, Sminus

    def bondhamiltonian(self, param, bondtype, z):
        """
        Geneates bond Hamiltonian of spin system

        Parameters
        ----------
        param : Dict
            parameter
        bondtype : Integer
            bond type index
        z     : Integer
            Coordinate number

        Returns
        -------
        hamiltonian : ndarray
            bond Hamiltonian

        """

        Sz, Sx, Splus, Sminus = self.localoperators(param)
        Sy = 0.5 * (Splus - Sminus)
        E = np.eye(Sz.shape[0])

        Jx = getparam(param, "Jx", bondtype, 1.0)
        Jy = getparam(param, "Jy", bondtype, 1.0)
        Jz = getparam(param, "Jz", bondtype, 1.0)
        BQ = getparam(param, "BQ", bondtype, 0.0)

        h = param.get("h", 0.0)
        G = param.get("G", 0.0)
        D = param.get("D", 0.0)

        ham = Jz * np.kron(Sz, Sz)
        ham += Jx * np.kron(Sx, Sx)
        ham -= Jy * np.kron(Sy, Sy)
        # ham += 0.5 * Jxy * (np.kron(Splus, Sminus) + np.kron(Sminus, Splus))
        ham += BQ * (
            np.kron(Sz, Sz) + 0.5 * np.kron(Splus, Sminus) + np.kron(Sminus, Splus)
        )
        ham -= (h / z) * (np.kron(Sz, E) + np.kron(E, Sz))
        ham -= (G / z) * (np.kron(Sx, E) + np.kron(E, Sx))
        ham += (D / z) * (np.kron(np.dot(Sz, Sz), E) + np.kron(E, np.dot(Sz, Sz)))

        return ham


def tenes_simple(param):
    """
    Generates TeNeS input for spin system on a square lattice

    Parameters
    ----------
    param : dict
        parameter

    Returns
    -------
    res : dict
        Dictionary file describing TeNeS input
    """

    model_param = param["model"]

    if model_param["type"] == "spin":
        model = SpinModel()
        lops = model.localoperators(model_param)[0:2]
    else:
        msg = "Unknown model type: {}".format(model_param["type"])
        raise RuntimeError(msg)

    lattice_param = param["lattice"]
    if lattice_param["type"] == "square lattice":
        lattice = SquareLattice(lattice_param)
    elif lattice_param["type"] == "honeycomb lattice":
        lattice = HoneycombLattice(lattice_param)
    else:
        msg = "Unknown lattice type: {}".format(lattice_param["type"])
        raise RuntimeError(msg)

    hams = [
        model.bondhamiltonian(model_param, bt, lattice.z)
        for bt in range(lattice.bondtypes)
    ]

    for bt in range(lattice.bondtypes):
        if not check_symmetric(hams[bt]):
            msg = "Bond Hamiltonian {} is not symmetric".format(bt)
            raise RuntimeError(msg)

    dict_observable = {"local_operator": lops, "hamiltonian": hams}

    res = {}
    if "parameter" in param:
        res["parameter"] = param["parameter"]
    res["lattice"] = lattice
    res["model"] = model
    res["observable"] = dict_observable
    if "correlation" in param:
        corparam = param["correlation"]
        dict_correlation = {}
        dict_correlation["r_max"] = corparam["r_max"]
        if "operators" in corparam:
            dict_correlation["operators"] = corparam["operators"]
        else:
            dict_correlation["operators"] = [[0, 0], [0, 1], [1, 1]]
        res["correlation"] = dict_correlation

    return res


def dump(param):
    ret = []

    ret.append("[parameter]")
    pparam = param["parameter"]
    for name in ("tensor", "simple_update", "full_update", "ctm", "random"):
        if name in pparam:
            ret.append("[parameter.{}]".format(name))
            for k, v in pparam[name].items():
                if isinstance(v, str):
                    ret.append('{} = "{}"'.format(k, v))
                else:
                    ret.append("{} = {}".format(k, v))

    ret.append("")
    ret.append("[lattice]")
    lattice = param["lattice"]
    ret.append('type = "{}"'.format(lattice.type))
    ret.append("L_sub = {}".format(lattice.Lsub))

    ret.append("")
    ret.append("[observable]")
    obs = param["observable"]
    lops = obs["local_operator"]
    hams = obs["hamiltonian"]
    ret.append("local_operator = [")
    for lop in lops:
        ret.append('"""')
        for line in array_to_strarr(lop):
            ret.append(line)
        ret.append('""",')
    ret.append("]")

    ret.append("hamiltonian = [")
    for ham in hams:
        ret.append('"""')
        for line in array_to_strarr(ham):
            ret.append(line)
        ret.append('""",')
    ret.append("]")

    ret.append('hamiltonian_bonds = """')
    for bond in lattice.bonds:
        ret.append(dumpbond(bond))
    ret.append('"""')

    ret.append("")
    ret.append("[evolution]")
    ret.append('simple_update = """')
    for bond in lattice.bonds:
        ret.append(dumpbond(bond))
    ret.append('"""')

    nsimple = lattice.bondtypes
    ret.append('full_update = """')
    for bond in lattice.bonds:
        ret.append(
            "{} {} {} {}".format(
                bond.source, bond.target, bond.dir, bond.type + nsimple
            )
        )
    ret.append('"""')

    ret.append("matrix = [")
    simple_update = pparam.get("simple_update", {})
    simple_tau = simple_update.get("tau", 0.01)
    for ham in hams:
        ret.append('"""')
        for line in array_to_strarr(make_evolution_tensor(ham, simple_tau)):
            ret.append(line)
        ret.append('""",')

    full_update = pparam.get("full_update", {})
    full_tau = full_update.get("tau", 0.01)
    for ham in hams:
        ret.append('"""')
        for line in array_to_strarr(make_evolution_tensor(ham, full_tau)):
            ret.append(line)
        ret.append('""",')

    ret.append("]")

    if "correlation" in param:
        corparam = param["correlation"]
        ret.append("")
        ret.append("[correlation]")
        ret.append("r_max = {}".format(corparam["r_max"]))
        ret.append("operators = {}".format(corparam["operators"]))

    return "\n".join(ret)


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(
        description="Simple input generator for TeNeS", add_help=True
    )

    parser.add_argument("input", help="Input TOML file")

    parser.add_argument(
        "-o", "--output", dest="output", default="input.toml", help="Output TOML file"
    )

    args = parser.parse_args()
    if args.input == args.output:
        print("The names of input and output are the same")
        sys.exit(1)
    res = tenes_simple(toml.load(args.input))

    with open(args.output, "w") as f:
        f.write(dump(res))
        f.write("\n")
