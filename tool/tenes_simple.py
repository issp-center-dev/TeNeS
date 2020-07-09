# TeNeS - Massively parallel tensor network solver
# Copyright (C) 2019- The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses

import re

from collections import namedtuple
from itertools import chain, product
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

from scipy.linalg import expm, norm

import toml

TeNeSInput = namedtuple("TeNeSInput", "param tensor ham obs")


def float_to_str(v: float) -> str:
    if np.isnan(v):
        return "nan"
    ret = ""
    if v < 0.0:
        ret += "-"
        v *= -1
    if np.isinf(v):
        ret += "inf"
        return ret
    if v == 0.0:
        return "0.0"
    exponent = int(np.floor(np.log10(v)))
    if -4 <= exponent <= 4:
        ret += str(v)
    else:
        ret += "{}e{}".format(v * 10 ** (-exponent), exponent)
    return ret


def lower_dict(d: dict) -> dict:
    ks = list(d.keys())
    for k in ks:
        if isinstance(d[k], dict):
            d[k] = lower_dict(d[k])
        if k.islower():
            continue
        d[k.lower()] = d[k]
        d.pop(k)
    return d


def index2coord(index: int, X: int) -> Tuple[int, int]:
    return index % X, index // X


def coord2index(x: int, y: int, X: int) -> int:
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


class SubLattice:
    def __init__(self, vdim: List[int], is_vacancy: bool = False):
        self.sites = []
        self.vdim = vdim
        self.is_vacancy = is_vacancy

    def add_site(self, site: int):
        self.sites.append(site)

    def to_dict(self, physdim: int) -> Dict[str, Any]:
        ret = {}
        ret["index"] = self.sites
        if self.is_vacancy:
            ret["physical_dim"] = 1
        else:
            ret["physical_dim"] = physdim
        ret["virtual_dim"] = self.vdim
        return ret


class Lattice(object):
    def __init__(self, param: Dict[str, Any]):
        self.type = ""
        self.z = 0
        self.skew = 0
        self.L = param["l"]
        self.W = param.get("w", self.L)
        self.vdim = param["virtual_dim"]
        self.sublattice = []
        self.bonds = [[[] for j in range(3)] for i in range(3)]
        self.initial_states = param.get("initial", "random")
        self.noise = param.get("noise", 1e-2)
        self.coords = []
        self.latticevector = np.eye(2)

    def to_dict(self, physdim: int) -> Dict[str, Any]:
        ret = {}
        ret["L_sub"] = [self.L, self.W]
        ret["skew"] = self.skew
        ret["unitcell"] = [sub.to_dict(physdim) for sub in self.sublattice]
        return ret

    def cartesian_coordinate(self, x: int, y: int) -> np.ndarray:
        return np.array([x, y])

    def write_coordinates(self, f):
        f.write("# coord_version = 1\n")
        f.write("# name = {}\n".format(self.type))
        f.write(
            "# a0 = {} {}\n".format(self.latticevector[0, 0], self.latticevector[0, 1])
        )
        f.write(
            "# a1 = {} {}\n".format(self.latticevector[1, 0], self.latticevector[1, 1])
        )
        f.write("# $1: index\n")
        f.write("# $2: x\n")
        f.write("# $3: y\n")
        f.write("\n")
        for i, c in enumerate(self.coords):
            if c is not None:
                # f.write("{} {} {}\n".format(i, c[0], c[1]))
                x, y = index2coord(i, self.L)
                C = self.cartesian_coordinate(x, y)
                f.write("{} {} {}\n".format(i, C[0], C[1]))

    def write_bonds(self, f, nnlevel):
        for i, bonds in enumerate(self.bonds[nnlevel]):
            for b in bonds:
                source = b.source
                x, y = index2coord(source, self.L)
                c = self.cartesian_coordinate(x, y)
                f.write("{} {} {}\n".format(c[0], c[1], i))
                c = self.cartesian_coordinate(x + b.dx, y + b.dy)
                f.write("{} {} {}\n".format(c[0], c[1], i))
                f.write("\n\n")

    def numsites(self):
        return self.L * self.W

    def valid_sites(self):
        ret = []
        for sl in self.sublattice:
            if not sl.is_vacancy:
                ret += sl.sites
        ret.sort()
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

        self.latticevector = np.diag([L, W])

        if self.initial_states == "ferro":
            self.sublattice = [SubLattice([self.vdim] * 4)]
        elif self.initial_states == "antiferro":
            self.sublattice = [SubLattice([self.vdim] * 4), SubLattice([self.vdim] * 4)]
        elif self.initial_states == "random":
            self.sublattice = [SubLattice([self.vdim] * 4)]

        for source in range(L * W):
            x, y = index2coord(source, L)
            if self.initial_states == "antiferro":
                if (x + y) % 2 == 0:
                    self.sublattice[0].add_site(source)
                else:
                    self.sublattice[1].add_site(source)

            self.coords.append(np.array([x, y]))

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

        self.L *= 2

        L, W = self.L, self.W
        self.skew = W % L

        vdim = self.vdim
        self.sublattice.append(SubLattice([vdim, 1, vdim, vdim]))
        self.sublattice.append(SubLattice([vdim, vdim, vdim, 1]))

        self.coords = [np.zeros(2) for _ in range(L * W)]

        NX = L // 2
        NY = W

        self.latticevector = np.array([[np.sqrt(3.0), 0.0], [np.sqrt(3.0) / 2, 1.5]])
        self.latticevector *= np.array([[NX], [NY]])
        a0 = np.array([np.sqrt(3.0), 0.0])
        a1 = np.array([np.sqrt(3.0) / 2, 1.5])
        other = (a0 + a1) / 3.0

        for y in range(NY):
            for X in range(NX):
                c = a0 * X + a1 * y

                #
                # sublattice A
                #
                x = (2 * X + y) % L
                index = coord2index(x, y, L)
                self.coords[index] = c
                self.sublattice[0].add_site(index)

                # 1st neighbors
                self.bonds[0][0].append(Bond(index, 1, 0))

                # 2nd neighbors
                self.bonds[1][0].append(Bond(index, -1, 1))
                self.bonds[1][1].append(Bond(index, 1, 1))
                self.bonds[1][2].append(Bond(index, 2, 0))

                # 3rd neighbors
                self.bonds[2][1].append(Bond(index, 2, -1))
                self.bonds[2][2].append(Bond(index, 0, 1))

                #
                # sublattice B
                #
                x = (2 * X + y + 1) % L
                index = coord2index(x, y, L)
                self.coords[index] = c + other
                self.sublattice[1].add_site(index)

                # 1st neighbors
                self.bonds[0][1].append(Bond(index, 1, 0))
                self.bonds[0][2].append(Bond(index, 0, 1))

                # 2nd neighbors
                self.bonds[1][0].append(Bond(index, -1, 1))
                self.bonds[1][1].append(Bond(index, 1, 1))
                self.bonds[1][2].append(Bond(index, 2, 0))

                # 3rd neighbors
                self.bonds[2][0].append(Bond(index, 2, 1))

    def cartesian_coordinate(self, x: int, y: int) -> np.ndarray:
        a0 = np.array([np.sqrt(3.0), 0.0])
        a1 = np.array([np.sqrt(3.0) / 2, 1.5])
        other = (a0 + a1) / 3.0
        X = x - y
        if X % 2 == 0:
            return a0 * (X // 2) + a1 * y
        else:
            return a0 * (X // 2) + a1 * y + other


class TriangularLattice(Lattice):
    def __init__(self, param: Dict[str, Any]):
        super().__init__(param)
        self.type = "triangular lattice"
        self.z = 6

        L, W = self.L, self.W
        assert L > 1 and W > 1

        self.latticevector = np.array([[1.0, 0.0], [0.5, np.sqrt(3.0) / 2]])
        self.latticevector *= np.array([[L], [W]])

        a0 = np.array([1.0, 0.0])
        a1 = np.array([0.5, np.sqrt(3.0) / 2])

        vdim = self.vdim

        if self.initial_states == "ferro":
            self.sublattice.append(SubLattice([vdim] * 4))
        elif self.initial_states == "antiferro":
            self.sublattice.append(SubLattice([vdim] * 4))
            self.sublattice.append(SubLattice([vdim] * 4))
            self.sublattice.append(SubLattice([vdim] * 4))
            nhops = np.ones((L + 1, W + 1), dtype=np.int) * 100000000
            nhops[0, 0] = 0
        elif self.initial_states == "random":
            self.sublattice.append(SubLattice([vdim] * 4))

        for source in range(L * W):
            x, y = index2coord(source, L)
            if self.initial_states == "antiferro":
                nhop = nhops[x, y]
                nhops[x + 1, y] = min(nhop + 1, nhops[x + 1, y])
                nhops[x, y + 1] = min(nhop + 2, nhops[x, y + 1])
                nhops[x - 1, y + 1] = min(nhop + 1, nhops[x - 1, y + 1])
                if nhop % 3 == 0:
                    self.sublattice[0].add_site(source)
                elif nhop % 3 == 1:
                    self.sublattice[1].add_site(source)
                else:
                    self.sublattice[2].add_site(source)

            self.coords.append(a0 * x + a1 * y)

            # 1st neighbors
            self.bonds[0][0].append(Bond(source, 1, 0))
            self.bonds[0][1].append(Bond(source, 0, 1))
            self.bonds[0][2].append(Bond(source, -1, 1))

            # 2nd neighbors
            self.bonds[1][0].append(Bond(source, -1, 2))
            self.bonds[1][1].append(Bond(source, -2, 1))
            self.bonds[1][2].append(Bond(source, 1, 1))

            # 3rd neighbors
            self.bonds[2][0].append(Bond(source, 2, 0))
            self.bonds[2][1].append(Bond(source, 0, 2))
            self.bonds[2][2].append(Bond(source, -2, 2))

    def cartesian_coordinate(self, x: int, y: int) -> np.ndarray:
        a0 = np.array([1.0, 0.0])
        a1 = np.array([0.5, np.sqrt(3.0) / 2])
        return a0 * x + a1 * y


class KagomeLattice(Lattice):
    def __init__(self, param: Dict[str, Any]):
        super().__init__(param)
        self.type = "kagome lattice"
        self.z = 4

        self.L *= 2
        self.W *= 2

        L, W = self.L, self.W

        self.latticevector = np.array([[1.0, 0.0], [0.5, np.sqrt(3.0) / 2]])
        self.latticevector *= np.array([[L], [W]])

        a0 = np.array([1.0, 0.0])
        a1 = np.array([0.5, np.sqrt(3.0) / 2])

        vd = self.vdim
        self.sublattice.append(SubLattice([vd, vd, vd, vd], is_vacancy=False))
        self.sublattice.append(SubLattice([vd, 1, vd, 1], is_vacancy=False))
        self.sublattice.append(SubLattice([1, vd, 1, vd], is_vacancy=False))
        self.sublattice.append(SubLattice([1, 1, 1, 1], is_vacancy=True))

        for index in range(L * W):
            x, y = index2coord(index, L)

            if x % 2 == 0 and y % 2 == 0:
                #
                # sublattice A
                #
                self.sublattice[0].add_site(index)
                self.coords.append(a0 * x + a1 * y)

                # 1st neighbors
                self.bonds[0][0].append(Bond(index, 1, 0))
                self.bonds[0][0].append(Bond(index, 0, 1))

                # 2nd neighbors
                self.bonds[1][0].append(Bond(index, -1, 2))
                self.bonds[1][0].append(Bond(index, -2, 1))

                # 3rd neighbors
                self.bonds[2][0].append(Bond(index, 2, 0))
                self.bonds[2][0].append(Bond(index, 0, 2))
                self.bonds[2][1].append(Bond(index, -2, 2))
            elif x % 2 == 1 and y % 2 == 0:
                #
                # sublattice B
                #
                self.sublattice[1].add_site(index)
                self.coords.append(a0 * x + a1 * y)

                # 1st neighbors
                self.bonds[0][1].append(Bond(index, 1, 0))
                self.bonds[0][0].append(Bond(index, -1, 1))

                # 2nd neighbors
                self.bonds[1][0].append(Bond(index, 1, 1))
                self.bonds[1][0].append(Bond(index, -1, 2))

                # 3rd neighbors
                self.bonds[2][0].append(Bond(index, 2, 0))
                self.bonds[2][0].append(Bond(index, -2, 2))
                self.bonds[2][1].append(Bond(index, 0, 2))
            elif x % 2 == 0 and y % 2 == 1:
                #
                # sublattice C
                #
                self.sublattice[2].add_site(index)
                self.coords.append(a0 * x + a1 * y)

                # 1st neighbors
                self.bonds[0][1].append(Bond(index, 0, 1))
                self.bonds[0][1].append(Bond(index, -1, 1))

                # 2nd neighbors
                self.bonds[1][0].append(Bond(index, 1, 1))
                self.bonds[1][0].append(Bond(index, -2, 1))

                # 3rd neighbors
                self.bonds[2][0].append(Bond(index, 0, 2))
                self.bonds[2][0].append(Bond(index, -2, 2))
                self.bonds[2][1].append(Bond(index, 2, 0))
            else:
                #
                # sublattice D (vacancy)
                #
                self.sublattice[3].add_site(index)
                self.coords.append(None)

    def cartesian_coordinate(self, x: int, y: int) -> np.ndarray:
        a0 = np.array([1.0, 0.0])
        a1 = np.array([0.5, np.sqrt(3.0) / 2])
        return a0 * x + a1 * y


class Model(object):
    def __init__(self):
        self.N = 0
        self.onesite_ops = []
        self.params = [[]]
        self.ham_list = [[]]
        pass

    def onesite_observables_as_dict(self) -> List[Dict[str, Any]]:
        ret = []

        for i, (name, op) in enumerate(zip(self.onesite_ops_name, self.onesite_ops)):
            dic = {}
            dic["group"] = i
            dic["name"] = name
            dic["sites"] = []
            dic["dim"] = self.N
            dic["elements"] = "\n".join(dump_op(op))
            dic["is_real"] = np.all(np.isreal(op))
            ret.append(dic)
        return ret

    def twosite_observables_as_dict(self) -> List[Dict[str, Any]]:
        ret = []

        for i, (name, op) in enumerate(
            zip(self.twosite_ops_name, self.twosite_ops), start=1
        ):
            dic = {}
            dic["group"] = i
            dic["name"] = name
            dic["dim"] = [self.N, self.N]
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


def Sz(S: float) -> np.ndarray:
    N = int(2 * S) + 1
    ret = np.zeros((N, N))
    for i in range(N):
        m = S - i
        ret[i, i] = m
    return ret


def Splus(S: float) -> np.ndarray:
    N = int(2 * S) + 1
    ret = np.zeros((N, N))
    for i in range(1, N):
        m = S - i
        ret[i - 1, i] = np.sqrt((S - m) * (S + m + 1.0))
    return ret


def Sminus(S: float) -> np.ndarray:
    N = int(2 * S) + 1
    ret = np.zeros((N, N))
    for i in range(N - 1):
        m = S - i
        ret[i + 1, i] = np.sqrt((S + m) * (S - m + 1.0))
    return ret


def Sx(S: float) -> np.ndarray:
    return 0.5 * (Splus(S) + Sminus(S))


def Sy(S: float) -> np.ndarray:
    return -0.5j * (Splus(S) - Sminus(S))


class SpinModel(Model):
    def __init__(self, param: Dict[str, Any]):
        super().__init__()

        self.S = param.get("s", 0.5)
        S = self.S
        assert int(2 * S) == 2 * S

        self.N = int(2 * S) + 1

        self.onesite_ops = [Sz(S), Sx(S), Sy(S)]
        self.onesite_ops_name = ["Sz", "Sx", "Sy"]

        self.twosite_ops = []
        self.twosite_ops_name = []
        for i, name in enumerate(self.onesite_ops_name):
            self.twosite_ops.append((i, i))
            self.twosite_ops_name.append("{}{}".format(name, name))
        self.read_params(param)
        super()._sort_ham_groups()

    def initial_states(self, num_sublattice: int) -> np.ndarray:
        def coherent_state(theta, phi):
            zeta = np.tan(0.5 * theta) * np.complex(np.cos(phi), np.sin(phi))
            sm = Sminus(self.S)
            U = expm(zeta * sm)
            ret = np.zeros(self.N)
            ret[0] = 1.0
            ret = np.dot(U, ret)
            return ret / norm(ret)

        ret = np.zeros((num_sublattice, self.N))
        if num_sublattice == 1:
            ret[0, 0] = 1.0
        elif num_sublattice == 2:
            ret[0, 0] = 1.0
            ret[1, 1] = 1.0
        elif num_sublattice == 3:
            ret[0, 0] = 1.0
            v = coherent_state(2 * np.pi / 3, 0.0)
            ret[1, :] = v.real
            v = coherent_state(2 * np.pi / 3, np.pi)
            ret[2, :] = v.real
        return ret

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
        Hx: float
            Magnetic field along Sx
        Hy: float
            Magnetic field along Sy
        Hz: float
            Magnetic field along Sz
            (longitudinal field)
        D: float
            Onsite spin anisotropy

        Returns
        -------
        ham: np.ndarray
            bond Hamiltonian
        """

        Jz = args.get("jz", 0.0)
        Jx = args.get("jx", 0.0)
        Jy = args.get("jy", 0.0)
        B = args.get("b", 0.0)

        hx = args.get("hx", 0.0) / z
        hy = args.get("hy", 0.0) / z
        hz = args.get("hz", 0.0) / z

        D = args.get("d", 0.0) / z

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
            val -= hz * (Sz[in1, out1] * E[in2, out2] + E[in1, out1] * Sz[in2, out2])
            val -= hx * (Sx[in1, out1] * E[in2, out2] + E[in1, out1] * Sx[in2, out2])
            val -= hy * (Sy[in1, out1] * E[in2, out2] + E[in1, out1] * Sy[in2, out2])
            val -= D * (
                Sz[in1, out1] ** 2 * E[in2, out2] + E[in1, out1] * Sz[in2, out2] ** 2
            )

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

        repat_J = re.compile("^j([012]?)('{0,2})([xyz]?)$")
        repat_B = re.compile("^b([012]?)('{0,2})$")

        for key in modelparam.keys():
            if key.startswith("j"):
                ma = repat_J.match(key)
                if not ma:
                    msg = "Unknown keyname {}".format(key)
                    raise RuntimeError(msg)

                gr = ma.groups()
                types = [int(gr[0])] if gr[0] else [0, 1, 2]
                n = len(gr[1])
                names = ["j" + gr[2]] if gr[2] else ["jx", "jy", "jz"]

                update(types, names, n, key)

            if key.startswith("b"):
                ma = repat_B.match(key)
                if not ma:
                    msg = "Unknown keyname {}".format(key)
                    raise RuntimeError(msg)
                gr = ma.groups()
                types = [int(gr[0])] if gr[0] else [0, 1, 2]
                n = len(gr[1])
                update(types, ["b"], n, key)

        if "h" in modelparam:
            print(
                'WARNING: "h" for the longitudial field is deprecated, and will be removed in future.'
            )
            print('         Use "hz" instead.')
            if "hz" in modelparam:
                print('ERROR: Both "hz" and "h" are specified. Use "hz".')
                sys.exit(1)
            update(range(3), ["hz"], 0, "h")
        else:
            update(range(3), ["hz"], 0, "hz")

        if "g" in modelparam:
            print(
                'WARNING: "g" for the transverse field is deprecated, and will be removed in future.'
            )
            print('         Use "hx" instead.')
            if "hx" in modelparam:
                print('ERROR: Both "hx" and "g" are specified. Use "hx".')
                sys.exit(1)
            update(range(3), ["hx"], 0, "g")
        else:
            update(range(3), ["hx"], 0, "hx")

        update(range(3), ["hy"], 0, "hy")
        update(range(3), ["d"], 0, "d")
        self.params = ret


def bose_creator(nmax: int) -> np.ndarray:
    N = nmax + 1
    ret = np.zeros((N, N))
    for i in range(N - 1):
        ret[i, i + 1] = np.sqrt(i + 1)
    return ret


def bose_annihilator(nmax: int) -> np.ndarray:
    N = nmax + 1
    ret = np.zeros((N, N))
    for i in range(N - 1):
        ret[i + 1, i] = np.sqrt(i + 1)
    return ret


def bose_number(nmax: int) -> np.ndarray:
    N = nmax + 1
    ret = np.zeros((N, N))
    for i in range(N):
        ret[i, i] = i
    return ret


class BoseHubbardModel(Model):
    def __init__(self, param: Dict[str, Any]):
        super().__init__()

        self.nmax = param.get("nmax", 1)
        nmax = self.nmax
        assert self.nmax > 0
        self.N = self.nmax + 1

        self.onesite_ops = [
            bose_number(nmax),
            bose_creator(nmax),
            bose_annihilator(nmax),
        ]
        self.onesite_ops_name = ["N", "Bdagger", "B"]

        self.twosite_ops = []
        self.twosite_ops_name = []
        for i, j in [(0, 0), (1, 2), (2, 1)]:
            self.twosite_ops.append((i, j))
            self.twosite_ops_name.append(
                "{}{}".format(self.onesite_ops_name[i], self.onesite_ops[j])
            )
        self.read_params(param)
        super()._sort_ham_groups()

    def initial_states(self, num_sublattice: int) -> np.ndarray:
        ret = np.zeros((num_sublattice, self.N))
        ret[0, self.N-1] = 1.0
        for i in range(1, num_sublattice):
            ret[i, 0] = 1.0
        return ret

    def model_hamiltonian(self, z: int, **args) -> np.ndarray:
        """
        Generate bond Hamiltonian of bosonic system
        
        Parameters
        ----------
        z: int 
            coordinate number
        t: float
            hopping constant, -t bi^d bj
        mu: float
            chemical potential, -mu ni
        U: float
            onsite repulsion, U/2 ni(ni-1)
        V: float
            offsite repulsion, V ninj

        Returns
        -------
        ham: np.ndarray
            bond Hamiltonian
        """

        t = args.get("t", 0.0)
        mu = args.get("mu", 0.0) / z
        U = args.get("u", 0.0) / z
        V = args.get("v", 0.0)

        E = np.eye(self.N)
        N, C, A = self.onesite_ops
        ham = np.zeros([self.N] * 4, dtype=np.complex)
        it = np.nditer(ham, flags=["multi_index"], order="F")
        while not it.finished:
            in1, in2, out1, out2 = it.multi_index
            val = 0.0
            val -= t * A[in1, out1] * C[in2, out2]
            val -= t * C[in1, out1] * A[in2, out2]
            val -= mu * N[in1, out1] * E[in2, out2]
            val -= mu * E[in1, out1] * N[in2, out2]
            val += 0.5 * U * N[in1, out1] * (N[in1, out1] - 1) * E[in2, out2]
            val += 0.5 * U * N[in2, out2] * (N[in2, out2] - 1) * E[in1, out1]
            val += V * N[in1, out1] * N[in2, out2]

            ham[in1, in2, out1, out2] = val
            it.iternext()
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

        repat = {
            "t": re.compile("^t([012]?)('{0,2})$"),
            "v": re.compile("^v([012]?)('{0,2})$"),
        }

        for key in modelparam.keys():
            if key == "type":
                continue
            for prefix in ("t", "v"):
                if key.startswith(prefix):
                    ma = repat[prefix].match(key)
                    if not ma:
                        msg = "Unknown keyname {}".format(key)
                        raise RuntimeError(msg)
                    gr = ma.groups()
                    types = [int(gr[0])] if gr[0] else [0, 1, 2]
                    n = len(gr[1])
                    update(types, [prefix], n, key)
        for typ in ("mu", "u"):
            update(range(3), [typ], 0, typ)
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
    elif latname.startswith("kagome"):
        lattice = KagomeLattice(latparam)
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
    elif modelparam["type"] == "boson":
        model = BoseHubbardModel(modelparam)
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

    param = lower_dict(param)
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
                    ret.append("{} = {}".format(k, "true" if v else "false"))
                elif isinstance(v, float):
                    ret.append("{} = {}".format(k, float_to_str(v)))
                else:
                    ret.append("{} = {}".format(k, v))
    ret.append("")

    ret.append("[tensor]")
    ret.append('type = "{}"'.format(lattice.type))
    ret.append("L_sub = [{}, {}]".format(lattice.L, lattice.W))
    ret.append("skew = {}".format(lattice.skew))

    ret.append("")

    num_sublattice = 0
    for sl in lattice.sublattice:
        if not sl.is_vacancy:
            num_sublattice += 1
    st = model.initial_states(num_sublattice)
    for i, sl in enumerate(lattice.sublattice):
        ret.append("[[tensor.unitcell]]")
        ret.append("virtual_dim = {}".format(sl.vdim))
        ret.append("index = {}".format(sl.sites))
        if sl.is_vacancy:
            ret.append("physical_dim = {}".format(1))
            ret.append("initial_state = [1.0]")
        else:
            ret.append("physical_dim = {}".format(model.N))
            if lattice.initial_states == "random":
                ret.append("initial_state = [0.0]")
            else:
                v = ", ".join(map(str, st[i, :]))
                ret.append("initial_state = [{}]".format(v))
        ret.append("noise = {}".format(lattice.noise))
        ret.append("")

    for ham in hams:
        ret.append("[[hamiltonian]]")
        ret.append("dim = {}".format([model.N, model.N]))
        ret.append('bonds = """')
        for bond in ham.bonds:
            ret.append(dumpbond(bond))
        ret.append('"""')
        ret.append('elements = """')
        for line in dump_op(ham.elements):
            ret.append(line)
        ret.append('"""')
        ret.append("")

    vsites = lattice.valid_sites()
    if len(vsites) == lattice.numsites():
        vsites = []
    ret.append("[observable]")
    lops = model.onesite_observables_as_dict()
    is_complex = True
    if "general" in pparam:
        is_complex = not pparam["general"].get("is_real", False)
    groups = []
    for lop in lops:
        if is_complex or lop["is_real"]:
            ret.append("[[observable.onesite]]")
            ret.append('name = "{}"'.format(lop["name"]))
            ret.append("group = {}".format(lop["group"]))
            groups.append(lop["group"])
            ret.append("sites = {}".format(vsites))
            ret.append("dim = {}".format(lop["dim"]))
            ret.append('elements = """')
            for line in lop["elements"].split("\n"):
                ret.append(line)
            ret.append('"""')
            ret.append("")
    groups = list(set(groups))
    groups.sort()

    for ham in hams:
        ret.append("[[observable.twosite]]")
        ret.append('name = "hamiltonian"')
        ret.append("group = 0")
        ret.append("dim = {}".format([model.N, model.N]))
        ret.append('bonds = """')
        for bond in ham.bonds:
            ret.append(dumpbond(bond))
        ret.append('"""')
        ret.append('elements = """')
        for line in dump_op(ham.elements):
            ret.append(line)
        ret.append('"""')
        ret.append("")

    # for i in range(len(model.onesite_ops)):
    #     oo = model.onesite_ops[i]
    #     if not is_complex:
    #         if not np.all(np.isreal(oo)):
    #             v = np.einsum("ij,kl -> ikjl", oo, oo)
    #             if not np.all(np.isreal(v)):
    #                 continue
    #     ret.append("[[observable.twosite]]")
    #     ret.append('name = "{name}{name}"'.format(name=model.onesite_ops_name[i]))
    #     ret.append("group = {}".format(i + 1))
    #     ret.append("dim = {}".format([model.N] * 2))
    #     ret.append('bonds = """')
    #     for bond in chain(*lattice.bonds[0]):
    #         ret.append(dumpbond(bond))
    #     ret.append('"""')
    #     if is_complex or np.all(np.isreal(oo)):
    #         ret.append("ops = {}".format([i] * 2))
    #     else:
    #         v = np.einsum("ij,kl -> ikjl", oo, oo)
    #         if is_complex or np.all(np.isreal(v)):
    #             ret.append('elements = """')
    #             for line in dump_op(v):
    #                 ret.append(line)
    #             ret.append('"""')
    #     ret.append("")
    k = 1
    for i, j in model.twosite_ops:
        oi = model.onesite_ops[i]
        oj = model.onesite_ops[j]
        v = np.einsum("ij,kl -> ikjl", oi, oj)
        if not (is_complex or np.all(np.isreal(v))):
            continue
        ret.append("[[observable.twosite]]")
        ret.append('name = "{}{}"'.format(model.onesite_ops_name[i], model.onesite_ops_name[j]))
        ret.append("group = {}".format(k))
        k += 1
        ret.append("dim = {}".format([model.N] * 2))
        ret.append('bonds = """')
        for bond in chain(*lattice.bonds[0]):
            ret.append(dumpbond(bond))
        ret.append('"""')
        if is_complex or (np.all(np.isreal(oi)) and np.all(np.isreal(oj))):
            ret.append("ops = {}".format([i, j]))
        else:
            v = np.einsum("ij,kl -> ikjl", oi, oj)
            ret.append('elements = """')
            for line in dump_op(v):
                ret.append(line)
            ret.append('"""')
        ret.append("")

    if "correlation" in param:
        corparam = param["correlation"]
        ret.append("[correlation]")
        ret.append("r_max = {}".format(corparam["r_max"]))
        if not "operators" in corparam:
            corparam["operators"] = []
            for g in groups:
                corparam["operators"].append([g, g])

        ret.append("operators = [")
        for ops in corparam["operators"]:
            ret.append("  {},".format(ops))
        ret.append("]")

    return "\n".join(ret), lattice


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
    parser.add_argument(
        "-c",
        "--coordinatefile",
        dest="coords",
        default="coordinates.dat",
        help="Site Information file",
    )
    parser.add_argument(
        "-b", "--bondfile", dest="bondfile", default="", help="Bond Information file",
    )
    parser.add_argument(
        "-v", "--version", dest="version", action="version", version="1.1.0"
    )

    args = parser.parse_args()

    if args.input == args.output:
        print("The names of input and output are the same")
        sys.exit(1)
    res, lattice = tenes_simple(toml.load(args.input))

    with open(args.output, "w") as f:
        f.write(res)
        f.write("\n")

    with open(args.coords, "w") as f:
        lattice.write_coordinates(f)

    if args.bondfile:
        for nn in range(3):
            with open("{}-{}.dat".format(args.bondfile, nn), "w") as f:
                lattice.write_bonds(f, nn)
