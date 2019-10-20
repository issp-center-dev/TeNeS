from __future__ import print_function

from itertools import chain, product
from collections import namedtuple
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


Bond = namedtuple("Bond", "source target type")


class SquareLattice:
    def __init__(self, param=None):
        self.Lsub = [0, 0]
        self.z = 4
        self.bondtypes = 2
        if param is not None:
            self.load(param)

    def load(self, param):
        self.L_sub = make_Lsub(param["L_sub"])
        X, Y = self.L_sub

        self.horizontal_bonds = []
        for leftx in range(0, X):
            for lefty in range(0, Y):
                rightx = (leftx + 1) % X
                righty = lefty
                left = coord2index(leftx, lefty, X)
                right = coord2index(rightx, righty, X)
                self.horizontal_bonds.append(Bond(left, right, 0))

        self.vertical_bonds = []
        for bottomy in range(0, Y):
            for bottomx in range(X):
                topx = bottomx
                topy = (bottomy + 1) % Y
                bottom = coord2index(bottomx, bottomy, X)
                top = coord2index(topx, topy, X)
                self.vertical_bonds.append(Bond(bottom, top, 1))


class HoneycombLattice:
    def __init__(self, param=None):
        self.z = 3
        self.bondtypes = 3
        if param is not None:
            self.load(param)

    def load(self, param):
        self.L_sub = make_Lsub(param["L_sub"])
        if not all(map(lambda x: x % 2 == 0, self.L_sub)):
            msg = "All elements of Lsub must be even for a honeycomb lattice."
            raise RuntimeError(msg)
        X, Y = self.L_sub

        sublat_a = []
        sublat_b = []
        for x, y in product(range(X), range(Y)):
            index = x + y * X
            if (x + y) % 2 == 0:
                sublat_a.append(index)
            else:
                sublat_b.append(index)

        self.horizontal_bonds = []
        self.vertical_bonds = []

        for index in sublat_a:
            x, y = index2coord(index, X)
            rightx = (x + 1) % X
            topy = (y + 1) % Y
            rightindex = coord2index(rightx, y, X)
            topindex = coord2index(x, topy, X)
            self.horizontal_bonds.append(Bond(index, rightindex, 0))  # x-bond
            self.vertical_bonds.append(Bond(index, topindex, 2))  # z-bond

        for index in sublat_b:
            x, y = index2coord(index, X)
            rightx = (x + 1) % X
            rightindex = coord2index(rightx, y, X)
            self.horizontal_bonds.append(Bond(index, rightindex, 1))  # y-bond


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
                Splus[i, i - 1] = np.sqrt((S - m) * (S + m + 1.0))
                Sminus[i - 1, i] = np.sqrt((S + m) * (S - m + 1.0))
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

        S = param.get("S", 0.5)
        Sz, Sx, Splus, Sminus = self.localoperators(param)
        E = np.eye(Sz.shape[0])

        Jz = getparam(param, "Jz", bondtype, 1.0)
        Jxy = getparam(param, "Jxy", bondtype, 1.0)
        BQ = getparam(param, "BQ", bondtype, 0.0)

        h = param.get("h", 0.0)
        G = param.get("G", 0.0)
        D = param.get("D", 0.0)

        ham = Jz * np.kron(Sz, Sz)
        ham += 0.5 * Jxy * (np.kron(Splus, Sminus) + np.kron(Sminus, Splus))
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
        model.bondhamiltonian(model_param, bt, lattice.z) for bt in range(lattice.bondtypes)
    ]

    pparam = param.get("parameter", {})
    simple_update = pparam.get("simple_update", {})
    simple_tau = simple_update.get("tau", 0.01)
    simple_evols = [make_evolution_tensor(H, simple_tau) for H in hams]
    full_update = pparam.get("full_update", {})
    full_tau = full_update.get("tau", 0.01)
    full_evols = [make_evolution_tensor(H, full_tau) for H in hams]

    nsimple_evols = len(simple_evols)

    simpleupdate_str = []
    fullupdate_str = []

    for bond in lattice.horizontal_bonds:
        simpleupdate_str.append(
            "{} {} h {}".format(bond.source, bond.target, bond.type)
        )
        fullupdate_str.append(
            "{} {} h {}".format(bond.source, bond.target, bond.type + nsimple_evols)
        )

    for bond in lattice.vertical_bonds:
        simpleupdate_str.append(
            "{} {} v {}".format(bond.source, bond.target, bond.type)
        )
        fullupdate_str.append(
            "{} {} v {}".format(bond.source, bond.target, bond.type + nsimple_evols)
        )

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
        toml.dump(res, f)
