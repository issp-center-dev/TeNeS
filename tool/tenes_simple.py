import copy
import re

from collections import namedtuple
from itertools import product
from typing import Any, Dict, Iterable, List

import numpy as np

import toml

TeNeSInput = namedtuple("TeNeSInput", "param tensor ham obs")


def index2coord(index: int, X: int):
    return index % X, index // X


def coord2index(x: int, y: int, X: int):
    return x + y * X


def dump_op(op: np.ndarray) -> Iterable[str]:
    it = np.nditer(op, flags=["multi_index"], op_flags=["readonly"], order="F")
    while not it.finished:
        index = it.multi_index
        v = op[index]
        if np.abs(v) != 0.0:
            ret = " ".join(map(str, index))
            ret += " {} {}".format(np.real(v), np.imag(v))
            yield ret
        it.iternext()


def getparam(param: dict, name: str, index: int, default):
    P = param.get(name, default)
    if isinstance(P, list):
        if index > len(P):
            msg = "parameter {} has too few elements.".format(name)
            raise RuntimeError(msg)
        return P[index]
    return P


Bond = namedtuple("Bond", "source dx dy")

Hamiltonian = namedtuple("Hamiltonian", "elements bonds")


def dumpbond(bond: Bond) -> str:
    return "{} {} {}".format(bond.source, bond.dx, bond.dy)


class Lattice(object):
    def __init__(self, param: Dict[str, Any]):
        self.type = ""
        self.z = 0
        self.skew = 0
        self.L = param["L"]
        self.W = param.get("W", self.L)
        self.vdims = [[param["virtual_dim"]] * 4]
        self.sublattice = [[]]
        self.bonds = [[[] for j in range(3)] for i in range(3)]
        self.initial_states = param.get("initial", "random")
        self.noise = param.get("noise", 1e-2)

    def to_dict(self, physdim: int) -> Dict[str, Any]:
        ret = {}
        ret["L_sub"] = [self.L, self.W]
        ret["skew"] = self.skew

        ret["unitcell"] = []
        for i, sublat in enumerate(self.sublattice):
            unitcell = {}
            unitcell["index"] = sublat
            unitcell["physical_dim"] = physdim
            unitcell["virtual_dim"] = self.vdims[i]
            ret["unitcell"].append(unitcell)

        return ret


class SquareLattice(Lattice):
    def __init__(self, param: Dict[str, Any]):
        super().__init__(param)
        self.type = "square lattice"
        self.z = 4
        self.skew = 0
        L, W = self.L, self.W
        if W == 1:
            self.skew = 1
        assert L > 1

        if self.initial_states == "ferro":
            self.sublattice = [[]]
        elif self.initial_states == "antiferro":
            self.sublattice = [[], []]
            self.vdims.append(copy.copy(self.vdims[0]))

        for source in range(L * W):
            if self.initial_states == "antiferro":
                x, y = index2coord(source, L)
                if (x + y) % 2 == 0:
                    self.sublattice[0].append(source)
                else:
                    self.sublattice[1].append(source)

            # 1st neighbors
            self.bonds[0][0].append(Bond(source, 1, 0))
            self.bonds[0][1].append(Bond(source, 0, 1))

            # 2nd neighbors
            self.bonds[1][0].append(Bond(source, 1, 1))
            self.bonds[1][1].append(Bond(source, -1, 1))

            # 2nd neighbors
            self.bonds[2][0].append(Bond(source, 2, 0))
            self.bonds[2][1].append(Bond(source, 0, 2))


class HoneycombLattice(Lattice):
    def __init__(self, param: Dict[str, Any]):
        super().__init__(param)
        self.type = "honeycomb lattice"
        self.z = 3
        self.bondtypes = 3

        L, W = self.L, self.W
        if W % 2 != 0:
            self.skew = 1
        assert L > 1

        self.vdims.append(copy.copy(self.vdims[0]))
        self.vdims[0][0] = 1
        self.vdims[1][2] = 1

        self.sublattice.append([])

        for source in range(L * W):
            x, y = index2coord(source, L)
            if (x + y) % 2 == 0:
                # sublattice A
                self.sublattice[0].append(source)

                # 1st neighbors
                self.bonds[0][0].append(Bond(source, 1, 0))
                self.bonds[0][1].append(Bond(source, 0, 1))

                # 2nd neighbors
                self.bonds[1][0].append(Bond(source, 0, 2))
                self.bonds[1][1].append(Bond(source, 1, 1))
                self.bonds[1][2].append(Bond(source, -1, 1))

                # 3rd neighbors
                self.bonds[2][2].append(Bond(source, 1, 2))
            else:
                # sublattice B
                self.sublattice[1].append(source)

                # 1st neighbors
                self.bonds[0][2].append(Bond(source, 0, 1))

                # 2nd neighbors
                self.bonds[1][0].append(Bond(source, 0, 2))
                self.bonds[1][1].append(Bond(source, 1, 1))
                self.bonds[1][2].append(Bond(source, -1, 1))

                # 3rd neighbors
                self.bonds[2][0].append(Bond(source, 1, 0))
                self.bonds[2][1].append(Bond(source, -1, 2))


class TriangularLattice(Lattice):
    def __init__(self, param: Dict[str, Any]):
        super().__init__(param)
        self.type = "triangular lattice"
        self.z = 6

        L, W = self.L, self.W
        assert L > 1 and W > 1

        if self.initial_states == "ferro":
            self.sublattice = [[]]
        elif self.initial_states == "antiferro":
            self.sublattice = [[], []]
            self.vdims.append(copy.copy(self.vdims[0]))
            nhops = np.ones((L + 1, W + 1), dtype=np.int) * 100000000
            nhops[0, 0] = 0

        for source in range(L * W):
            if self.initial_states == "antiferro":
                x, y = index2coord(source, L)
                nhop = nhops[x, y]
                nhops[x + 1, y] = min(nhop + 1, nhops[x + 1, y])
                nhops[x, y + 1] = min(nhop + 1, nhops[x, y + 1])
                nhops[x + 1, y + 1] = min(nhop + 2, nhops[x + 1, y + 1])
                if nhop % 3 == 0:
                    self.sublattice[0].append(source)
                else:
                    self.sublattice[1].append(source)

            # 1st neighbors
            self.bonds[0][0].append(Bond(source, 1, 0))
            self.bonds[0][1].append(Bond(source, 0, 1))
            self.bonds[0][2].append(Bond(source, 1, 1))

            # 2nd neighbors
            self.bonds[1][0].append(Bond(source, 1, 2))
            self.bonds[1][1].append(Bond(source, 2, 1))
            self.bonds[1][2].append(Bond(source, 1, -1))

            # 3rd neighbors
            self.bonds[2][0].append(Bond(source, 2, 0))
            self.bonds[2][1].append(Bond(source, 0, 2))
            self.bonds[2][2].append(Bond(source, 2, 2))


class Model(object):
    def __init__(self):
        self.N = 0
        self.onesite_ops = []
        self.params = [[]]
        self.ham_list = [[]]
        pass

    def onesite_observables_as_dict(self) -> List[Dict[str, Any]]:
        ret = []

        for i, op in enumerate(self.onesite_ops):
            dic = {}
            dic["group"] = i
            dic["sites"] = []
            dic["dim"] = self.N
            dic["elements"] = "\n".join(dump_op(op))
            dic["is_real"] = np.all(np.isreal(op))
            ret.append(dic)
        return ret

    def bondhamiltonian(self, typ: int, z: int = 1) -> np.ndarray:
        n, t = self.ham_list[typ][0]
        return self.model_hamiltonian(z, **self.params[n][t])

    def _sort_ham_groups(self):
        num_nn = len(self.params)
        num_typ = len(self.params[0])
        ham_groups = list(range(num_nn * num_typ))
        for i, (n, typ) in enumerate(product(range(num_nn), range(num_typ))):
            lhs = self.params[n][typ]

            vs = list(lhs.values())
            if np.count_nonzero(vs) == 0:
                ham_groups[i] = -1

            for i2, (n2, typ2) in enumerate(product(range(num_nn), range(num_typ))):
                if i2 <= i:
                    continue
                if ham_groups[i2] != i2:
                    continue
                rhs = self.params[n2][typ2]
                if lhs == rhs:
                    ham_groups[i2] = i
        self.ham_list = []
        for t in np.unique(ham_groups):
            if t < 0:
                continue
            hl = []
            for i, (n, typ) in enumerate(product(range(num_nn), range(num_typ))):
                if ham_groups[i] == t:
                    hl.append((n, typ))
            self.ham_list.append(hl)


class SpinModel(Model):
    def __init__(self, param: Dict[str, Any]):
        super().__init__()

        S = param.get("S", 0.5)
        assert int(2 * S) == 2 * S

        self.N = int(2 * S) + 1
        M = self.N
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
        Sy = 0.5j * (Sminus - Splus)

        self.onesite_ops = [Sz, Sx, Sy]
        self.read_params(param)
        super()._sort_ham_groups()

    def model_hamiltonian(self, z: int, **args) -> np.ndarray:
        """
        Geneates bond Hamiltonian of spin system

        Parameters
        ----------
        z: int
            coordinate number
        Jx: float
            SxSx interaction
        Jy: float
            SySy interaction
        Jz: float
            SzSz interaction
        B: float
            Bilinear-Biquadratic inteaction
        H: float
            Longitudinal magnetic field
        G: float
            Transverse magnetic field
        D: float
            Onsite spin anisotropy

        Returns
        -------
        ham: np.ndarray
            bond Hamiltonian

        """

        Jz = args.get("Jz", 0.0)
        Jx = args.get("Jx", 0.0)
        Jy = args.get("Jy", 0.0)
        B = args.get("B", 0.0)
        h = args.get("H", 0.0) / z
        g = args.get("G", 0.0) / z
        D = args.get("D", 0.0) / z

        E = np.eye(self.N)
        Sz, Sx, Sy = self.onesite_ops
        ham = np.zeros([self.N] * 4, dtype=np.complex)
        SS = np.zeros([self.N] * 4)
        it = np.nditer(ham, flags=["multi_index"], order="F")
        while not it.finished:
            in1, in2, out1, out2 = it.multi_index
            val = 0.0
            val += Jz * Sz[in1, out1] * Sz[in2, out2]
            val += Jx * Sx[in1, out1] * Sx[in2, out2]
            val += Jy * np.real(Sy[in1, out1] * Sy[in2, out2])
            val -= h * (Sz[in1, out1] * E[in2, out2] + E[in1, out1] * Sz[in2, out2])
            val -= D * (
                Sz[in1, out1] ** 2 * E[in2, out2] + E[in1, out1] * Sz[in2, out2] ** 2
            )
            val -= g * (Sx[in1, out1] * E[in2, out2] + E[in1, out1] * Sx[in2, out2])

            SS[in1, in2, out1, out2] = (
                Sz[in1, out1] * Sz[in2, out2]
                + Sx[in1, out1] * Sx[in2, out2]
                + np.real(Sy[in1, out1] * Sy[in2, out2])
            )

            ham[in1, in2, out1, out2] = val
            it.iternext()
        ham += B * SS ** 2
        return ham

    def read_params(self, modelparam: Dict[str, Any]):
        ret = [
            [{}, {}, {}],  # 1st neighbors
            [{}, {}, {}],  # 2nd neighbors
            [{}, {}, {}],  # 3rd neighbors
        ]

        def update(types, names, n, key):
            for typ in types:
                for name in names:
                    if name in ret[n][typ]:
                        msg = "{} is defined twice".format(key)
                        raise RuntimeError(msg)
                    ret[n][typ][name] = modelparam.get(key, 0.0)

        repat_J = re.compile("^J([012]?)('{0,2})([xyz]?)$")
        repat_B = re.compile("^B([012]?)('{0,2})$")

        for key in modelparam.keys():
            if key.startswith("J"):
                ma = repat_J.match(key)
                if not ma:
                    msg = "Unknown keyname {}".format(key)
                    raise RuntimeError(msg)

                gr = ma.groups()
                types = [int(gr[0])] if gr[0] else [0, 1, 2]
                n = len(gr[1])
                names = ["J" + gr[2]] if gr[2] else ["Jx", "Jy", "Jz"]

                update(types, names, n, key)

            if key.startswith("b"):
                ma = repat_B.match(key)
                if not ma:
                    msg = "Unknown keyname {}".format(key)
                    raise RuntimeError(msg)
                gr = ma.groups()
                types = [int(gr[0])] if gr[0] else [0, 1, 2]
                n = len(gr[1])
                update(types, ["B"], n, key)
        for name in ("H", "G", "D"):
            update(range(3), [name], 0, name)
        self.params = ret


def make_lattice(param: Dict[str, Any]) -> Lattice:
    """
    Parameters
    ----------
    param : Dict[str, Any]
        parameter

    Returns
    -------
    lattice: Lattice

    Raises
    ------
    RuntimeError
        Raises RuntimeError when given lattice name is not registerered
    """

    latparam = param["lattice"]
    latname = latparam["type"].strip('""')
    if latname.startswith("square"):
        lattice = SquareLattice(latparam)
    elif latname.startswith("honeycomb"):
        lattice = HoneycombLattice(latparam)
    elif latname.startswith("triangular"):
        lattice = TriangularLattice(latparam)
    else:
        msg = "Unknown lattice: {}".format(latname)
        raise RuntimeError(msg)
    return lattice


def make_model(param: Dict[str, Any]) -> Model:
    """
    Parameters
    ----------
    param : Dict[str, Any]
        parameter

    Returns
    -------
    model: Model

    Raises
    ------
    RuntimeError
        Raises RuntimeError when given model name is not registerered
    """
    modelparam = param["model"]
    if modelparam["type"] == "spin":
        model = SpinModel(modelparam)
    else:
        msg = "Unknown model type: {}".format(modelparam["type"])
        raise RuntimeError(msg)
    return model


def hamiltonians(lattice: Lattice, model: Model) -> List[Hamiltonian]:
    ret = []
    for i, hl_list in enumerate(model.ham_list):
        elem = model.bondhamiltonian(i, lattice.z)
        bonds = []
        for n, typ in hl_list:
            bonds += lattice.bonds[n][typ]
        ret.append(Hamiltonian(elem, bonds))
    return ret


def tenes_simple(param: Dict[str, Any]) -> str:
    """
    Parameters
    ----------
    param : Dict[str, Any]
        parameter
    """

    lattice = make_lattice(param)
    model = make_model(param)
    hams = hamiltonians(lattice, model)

    ret = []
    ret.append("[parameter]")
    pparam = param["parameter"]
    for name in ("general", "simple_update", "full_update", "ctm", "random"):
        if name in pparam:
            ret.append("[parameter.{}]".format(name))
            for k, v in pparam[name].items():
                if isinstance(v, str):
                    ret.append('{} = "{}"'.format(k, v))
                elif isinstance(v, bool):
                    ret.append('{} = {}'.format(k, "true" if v else "false"))
                else:
                    ret.append("{} = {}".format(k, v))
    ret.append("")

    ret.append("[tensor]")
    ret.append('type = "{}"'.format(lattice.type))
    ret.append("L_sub = [{}, {}]".format(lattice.L, lattice.W))
    ret.append("skew = {}".format(lattice.skew))

    ret.append("")
    st = [[0.0] * model.N]
    st[0][0] = 1.0
    st.append(st[0][-1::-1])
    for i, (sl, vdims) in enumerate(zip(lattice.sublattice, lattice.vdims)):
        ret.append("[[tensor.unitcell]]")
        ret.append("index = {}".format(sl))
        ret.append("physical_dim = {}".format(model.N))
        ret.append("virtual_dim = {}".format(vdims))
        ret.append("noise = {}".format(lattice.noise))
        if lattice.initial_states == "random":
            ret.append("initial_state = [0.0]")
        else:
            ret.append("initial_state = {}".format(st[i]))
        ret.append("")

    for ham in hams:
        ret.append("[[hamiltonian]]")
        ret.append("dims = {}".format([model.N, model.N]))
        ret.append('bonds = """')
        for bond in ham.bonds:
            ret.append(dumpbond(bond))
        ret.append('"""')
        ret.append('elements = """')
        for line in dump_op(ham.elements):
            ret.append(line)
        ret.append('"""')
        ret.append("")

    ret.append("[observable]")
    lops = model.onesite_observables_as_dict()
    is_complex = True
    if "general" in pparam:
        is_complex = not pparam["general"].get("is_real", False)
    groups = []
    for lop in lops:
        if is_complex or lop["is_real"]:
            ret.append("[[observable.onesite]]")
            ret.append("group = {}".format(lop["group"]))
            groups.append(lop["group"])
            ret.append("sites = {}".format(lop["sites"]))
            ret.append("dim = {}".format(lop["dim"]))
            ret.append('elements = """')
            for line in lop["elements"].split("\n"):
                ret.append(line)
            ret.append('"""')

    groups = list(set(groups))
    groups.sort()

    if "correlation" in param:
        corparam = param["correlation"]
        ret.append("")
        ret.append("[correlation]")
        ret.append("r_max = {}".format(corparam["r_max"]))
        if not "operators" in corparam:
            corparam["operators"] = []
            for g in groups:
                corparam["operators"].append([g,g])

        ret.append("operators = [")
        for ops in corparam["operators"]:
            ret.append("  {},".format(ops))
        ret.append("]")

    return "\n".join(ret)


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(
        description="Simple input generator for TeNeS", add_help=True
    )

    parser.add_argument("input", help="Input TOML file")

    parser.add_argument(
        "-o", "--output", dest="output", default="std.toml", help="Output TOML file"
    )

    args = parser.parse_args()
    if args.input == args.output:
        print("The names of input and output are the same")
        sys.exit(1)
    res = tenes_simple(toml.load(args.input))

    with open(args.output, "w") as f:
        f.write(res)
        f.write("\n")
