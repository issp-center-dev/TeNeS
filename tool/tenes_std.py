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

from collections import namedtuple
from itertools import product
from typing import TextIO, List, Tuple, Optional, cast, MutableMapping, Dict, Any, Union

import numpy as np

import scipy.sparse as sparse


def drop_comment(line: str) -> str:
    last = line.find("#")
    if last < 0:
        return line[:]
    else:
        return line[:last]


def lower_dict(d: MutableMapping) -> Dict[str, Any]:
    ks = list(d.keys())
    for k in ks:
        if isinstance(d[k], dict):
            d[k] = lower_dict(d[k])
        if k.islower():
            continue
        d[k.lower()] = d[k]
        d.pop(k)
    return cast(Dict[str, Any], d)


def all_positive(xs, v=0) -> bool:
    return all(map(lambda x: x > v, xs))


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


def value_to_str(v) -> str:
    if isinstance(v, str):
        return "'{}'".format(v)
    elif isinstance(v, bool):
        return "true" if v else "false"
    elif isinstance(v, float):
        return float_to_str(v)
    else:
        return "{}".format(v)


class Bond:
    source_site: int
    dx: int
    dy: int

    def __init__(self, source_site: int, dx: int, dy: int):
        self.source_site = source_site
        self.dx = dx
        self.dy = dy

    def is_site(self) -> bool:
        return self.dx == self.dy == 0

    def is_bond(self) -> bool:
        return not self.is_site()


def parse_bond(line: str) -> Bond:
    line = drop_comment(line)
    if not line:
        return None
    words = line.split()
    source_site = int(words[0])
    dx = int(words[1])
    dy = int(words[2])
    return Bond(source_site, dx, dy)


def load_tensor(elements_str: str, dims: List[int], atol: float = 1e-15) -> np.ndarray:
    A_re = np.zeros(dims)
    A_im = np.zeros(dims)
    ndim = A_re.ndim
    for line in elements_str.strip().splitlines():
        line = drop_comment(line)
        words = line.strip().split()
        if len(words) != len(dims) + 2:
            msg = "number of columns should be {} but {} given".format(
                len(dims) + 2, len(words)
            )
            raise RuntimeError(msg)

        index = [int(x) for x in words[0:ndim]]

        re = float(words[ndim])
        im = float(words[ndim + 1])
        A_re[tuple(index)] = re
        A_im[tuple(index)] = im

    is_real = np.all(np.isclose(A_im, 0.0, rtol=0.0, atol=atol))
    if is_real:
        return A_re
    else:
        return A_re + A_im * 1j


def is_hermite(A: np.ndarray) -> bool:
    nsite = len(A.shape) // 2
    input_dirs = A.shape[:nsite]
    output_dirs = A.shape[nsite:]
    if input_dirs != output_dirs:
        return False
    idir = int(np.prod(input_dirs))
    odir = int(np.prod(output_dirs))
    Amat = A.reshape((idir, odir))
    return bool(np.all(Amat.conjugate().transpose() == Amat))


class LocalTensor:
    """
    Local tensors

    Attributes
    ----------
    phys_dim : int
        Dimension of physical bond.
    virtual_dim : List[int]
        Dimenstions of four virtual bonds.
        Four bonds are ordered as [-x, +y, +x, -y].
    """

    phys_dim: int
    virtual_dim: List[int]

    def __init__(self, tensor_dict: dict = None):
        if tensor_dict is not None:
            self.load_dict(tensor_dict)

    def load_dict(self, tensor_dict: dict):
        """
        Loads information from tensor_dict

        Parameters
        ----------
        tensor_dict : dict
            Parameters
        """
        self.phys_dim = tensor_dict["physical_dim"]

        virtual_dim = tensor_dict["virtual_dim"]
        if isinstance(virtual_dim, int):
            self.virtual_dim = [virtual_dim] * 4
        else:
            self.virtual_dim = virtual_dim
        self.check()

    def check(self):
        """
        Checks if self stores valid information.

        Raises
        ------
        RuntimeError
            Raises RuntimeError if any field is invalid.
        """
        if not (isinstance(self.phys_dim, int) and self.phys_dim > 0):
            msg = "ERROR: physical_dim should be a positive integer"
            raise RuntimeError(msg)

        if not (
            isinstance(self.virtual_dim, list)
            and len(self.virtual_dim) == 4
            and all_positive(self.virtual_dim)
        ):
            msg = "ERROR: virtual_dim must be a positive integer or a list with 4 positive integers"
            raise RuntimeError(msg)


class Unitcell:
    """
    Unitcell

    Attributes
    ----------
    L : List[int]
        Lengths of unit cell, [Lx, Ly]
    sites : LocalTensor
        Sites in a unit cell
    """

    L: List[int]
    sites: List[LocalTensor]

    def __init__(self, lat_dict: dict = None):
        self.L = [1] * 2
        self.sites = cast(List[LocalTensor], [None])
        if lat_dict is not None:
            self.load_dict(lat_dict)

    def load_dict(self, lat_dict: dict):
        L = lat_dict["l_sub"]
        if isinstance(L, int):
            self.L = [L, L]
        else:
            self.L = L

        self.skew = lat_dict.get("skew", 0)

        N = self.L[0] * self.L[1]
        self.sites = cast(List[LocalTensor], [None] * N)
        for site in lat_dict["unitcell"]:
            index = site["index"]
            if isinstance(index, int):
                index = [index]
            if isinstance(index, list) and len(index) == 0:
                index = range(N)
            for i in index:
                self.sites[i] = LocalTensor(site)
        self.check()

    def numsites(self) -> int:
        return len(self.sites)

    def index2coord(self, index: int) -> Tuple[int, int]:
        d, m = divmod(index, self.L[0])
        return m, d

    def coord2index(self, x: int, y: int) -> int:
        x, y, _, _ = self.coord2supercoord(x, y)
        return x + self.L[0] * y

    def coord2supercoord(self, x: int, y: int) -> Tuple[int, int, int, int]:
        offset_y = y // self.L[1]
        x -= offset_y * self.skew
        offset_x = x // self.L[0]
        x -= offset_x * self.L[0]
        y -= offset_y * self.L[1]
        return x, y, offset_x, offset_y

    def wan2coord(self, site: int, ox: int, oy: int) -> Tuple[int, int]:
        x, y = self.index2coord(site)
        x += oy * self.skew
        x += ox * self.L[0]
        y += oy * self.L[1]
        return x, y

    def source_site(self, bond: Bond) -> int:
        return bond.source_site

    def source_coord(self, bond: Bond) -> Tuple[int, int]:
        return self.index2coord(self.source_site(bond))

    def target_site(self, bond: Bond) -> int:
        x, y = self.index2coord(bond.source_site)
        x += bond.dx
        y += bond.dy
        x, y, _, _ = self.coord2supercoord(x, y)
        return self.coord2index(x, y)

    def target_coord(self, bond: Bond) -> Tuple[int, int]:
        x, y = self.index2coord(bond.source_site)
        x += bond.dx
        y += bond.dy
        return x, y

    def target_offset(self, bond: Bond) -> Tuple[int, int]:
        x, y = self.target_coord(bond)
        _, _, ox, oy = self.coord2supercoord(x, y)
        return ox, oy

    def bond_displacement(self, bond: Bond) -> Tuple[int, int]:
        return bond.dx, bond.dy

    def bond_direction(self, bond: Bond) -> int:
        dx, dy = self.bond_displacement(bond)
        nhop = abs(dx) + abs(dy)
        assert nhop == 1, "{} {}".format(dx, dy)
        if dx == 1:
            direction = 2
        elif dx == -1:
            direction = 0
        elif dy == 1:
            direction = 1
        else:
            direction = 3
        return direction

    def bond_dim(self, bond: Bond) -> int:
        direction = self.bond_direction(bond)
        return self.sites[bond.source_site].virtual_dim[direction]

    def make_bond(self, source: int, direction: int) -> Bond:
        dx = 0
        dy = 0
        if direction == 0:
            dx = -1
        elif direction == 1:
            dy = 1
        elif direction == 2:
            dx = 1
        elif direction == 3:
            dy = -1
        else:
            raise RuntimeError()
        return Bond(source, dx, dy)

    def check(self) -> None:
        """
        Check if self stores valid information

        Raises
        ------
        RuntimeError
            Raises RuntimeError if any information is invalid.
        """
        if not (isinstance(self.L, list) and len(self.L) == 2 and all_positive(self.L)):
            msg = "ERROR: L_sub should be a positive integer or a list with 2 positive integers"
            raise RuntimeError(msg)

        failed = False
        for i, site in enumerate(self.sites):
            if site is None:
                print("ERROR: site {} is not defined".format(i))
                failed = True
        for i, isite in enumerate(self.sites):
            for idir in (1, 2):
                idim = isite.virtual_dim[idir]
                j = self.neighbor(i, idir)
                jdir = (idir + 2) % 4
                jdim = self.sites[j].virtual_dim[jdir]
                if idim != jdim:
                    print(
                        "ERROR: The dimension of bond {idir} of site {i} and that of {jdir} of {j} are mismatch.".format(
                            idir=idir, i=i, jdir=jdir, j=j
                        )
                    )
                    failed = True
        if failed:
            msg = "ERROR: some sites have problems"
            raise RuntimeError(msg)

    def neighbor(self, index: int, direction: int) -> int:
        x, y = self.index2coord(index)
        if direction == 0:
            X, Y, _, _ = self.coord2supercoord(x - 1, y)
        elif direction == 1:
            X, Y, _, _ = self.coord2supercoord(x, y + 1)
        elif direction == 2:
            X, Y, _, _ = self.coord2supercoord(x + 1, y)
        else:  # elif direction == 3:
            X, Y, _, _ = self.coord2supercoord(x, y - 1)
        return self.coord2index(X, Y)


class LatticeGraph:
    unitcell: Unitcell
    offset_x_min: int
    offset_y_min: int
    offset_x_max: int
    offset_y_max: int
    offset_Lx: int
    offset_Ly: int
    N: int
    path_pred: np.ndarray
    path_weight: np.ndarray

    def __init__(
        self,
        unitcell: Unitcell,
        offset_x_min: int,
        offset_y_min: int,
        offset_x_max: int,
        offset_y_max: int,
    ):
        self.unitcell = unitcell
        self.offset_x_min = offset_x_min
        self.offset_y_min = offset_y_min
        self.offset_x_max = offset_x_max
        self.offset_y_max = offset_y_max
        self.offset_Lx = offset_x_max - offset_x_min + 1
        self.offset_Ly = offset_y_max - offset_y_min + 1
        self.N = self.unitcell.numsites() * self.offset_Lx * self.offset_Ly
        A = np.zeros((self.N, self.N))
        for i in range(self.unitcell.numsites()):
            bonds = [unitcell.make_bond(i, direction) for direction in (0, 1, 2, 3)]
            for ioy, iox in product(
                range(offset_y_min, offset_y_max + 1),
                range(offset_x_min, offset_x_max + 1),
            ):
                for bond in bonds:
                    x, y = self.unitcell.target_coord(bond)
                    x, y, ox, oy = self.unitcell.coord2supercoord(x, y)
                    D = unitcell.bond_dim(bond)
                    if D == 1:
                        continue
                    if not offset_x_min <= ox + iox <= offset_x_max:
                        continue
                    if not offset_y_min <= oy + ioy <= offset_y_max:
                        continue
                    source_site = self.graph_site(
                        self.unitcell.source_site(bond), iox, ioy
                    )
                    target_site = self.graph_site(
                        self.unitcell.target_site(bond), iox + ox, ioy + oy
                    )

                    source_x, source_y = self.unitcell.wan2coord(
                        self.unitcell.source_site(bond), iox, ioy
                    )
                    target_x, target_y = self.unitcell.wan2coord(
                        self.unitcell.target_site(bond), iox, ioy
                    )
                    bond_x = (source_x + target_x) * 0.5
                    bond_y = (source_y + target_y) * 0.5
                    epsilon = 1.0e-6

                    weight = 1.0 / D / D - epsilon * (100 * bond_x - bond_y)
                    if ox == oy == 0:
                        A[source_site, target_site] = weight
                    else:
                        A[source_site, target_site] = 2 * weight
        self.path_weight, self.path_pred = sparse.csgraph.shortest_path(
            sparse.csr_matrix(A), return_predecessors=True
        )

    def graph_site(self, localsite: int, offset_x: int = 0, offset_y: int = 0) -> int:
        assert self.offset_x_min <= offset_x <= self.offset_x_max
        assert self.offset_y_min <= offset_y <= self.offset_y_max
        offset_x -= self.offset_x_min
        offset_y -= self.offset_y_min
        unitcell_index = offset_x + offset_y * self.offset_Lx
        return localsite + self.unitcell.numsites() * unitcell_index

    def graph_coords(self, index: int) -> Tuple[int, int, int]:
        assert 0 <= index <= self.N

        index, localsite = divmod(index, self.unitcell.numsites())
        offset_y, offset_x = divmod(index, self.offset_Lx)
        return localsite, offset_x, offset_y

    def make_path(self, bond: Bond) -> List[Bond]:
        ox, oy = self.unitcell.target_offset(bond)
        target_index = self.graph_site(self.unitcell.target_site(bond), ox, oy)
        path_pred = self.path_pred[self.graph_site(bond.source_site)]
        targets = [target_index]
        target_index = path_pred[target_index]
        while target_index >= 0:
            targets.append(target_index)
            target_index = path_pred[target_index]

        bonds: List[Bond] = []
        source_site, source_offset_x, source_offset_y = self.graph_coords(targets.pop())
        while len(targets) > 0:
            target_site, target_offset_x, target_offset_y = self.graph_coords(
                targets.pop()
            )
            ox = target_offset_x - source_offset_x
            oy = target_offset_y - source_offset_y
            x, y = self.unitcell.index2coord(source_site)
            X, Y = self.unitcell.wan2coord(target_site, ox, oy)
            b = Bond(source_site, X - x, Y - y)
            bonds.append(b)
            source_site = target_site
            source_offset_x = target_offset_x
            source_offset_y = target_offset_y
        return bonds


class SiteOperator:
    site: int
    elements: np.ndarray

    def __init__(self, site: int, elements: np.ndarray):
        self.site = site
        self.elements = elements

    def to_toml_strs(self, unitcell) -> List[str]:
        ret = []
        ret.append("site = {}".format(self.site))
        ret.append("dimensions = {}".format(list(self.elements.shape)))
        it = np.nditer(
            self.elements, flags=["multi_index"], op_flags=["readonly"], order="F"
        )
        ret.append('elements = """')
        while not it.finished:
            i = it.multi_index
            v = self.elements[i]
            if np.abs(v) == 0.0:
                it.iternext()
                continue
            line = ""
            for j in i:
                line += str(j) + " "
            line += " {} {}".format(np.real(v), np.imag(v))
            ret.append(line)
            it.iternext()
        ret.append('"""')
        return ret


class NNOperator:
    bond: Bond
    elements: Optional[np.ndarray]
    ops: Optional[List[int]]

    def __init__(
        self, bond: Bond, *, elements: np.ndarray = None, ops: List[int] = None
    ):
        self.bond = bond
        if elements is not None:
            self.elements = elements
            self.ops = None
        if ops is not None:
            self.ops = ops
            self.elements = None

    def to_toml_strs(self, unitcell: Unitcell) -> List[str]:
        ret = []
        ret.append("source_site = {}".format(self.bond.source_site))
        ret.append("source_leg = {}".format(unitcell.bond_direction(self.bond)))
        if self.elements is not None:
            ret.append("dimensions = {}".format(list(self.elements.shape)))
            it = np.nditer(
                self.elements, flags=["multi_index"], op_flags=["readonly"], order="F"
            )
            ret.append('elements = """')
            while not it.finished:
                i = it.multi_index
                v = self.elements[i]
                if np.abs(v) == 0.0:
                    it.iternext()
                    continue
                line = ""
                for j in i:
                    line += str(j) + " "
                line += " {} {}".format(np.real(v), np.imag(v))
                ret.append(line)
                it.iternext()
            ret.append('"""')
        else:
            ret.append("ops = {}".format(self.ops))
        return ret


Operator = Union[SiteOperator, NNOperator]


class OnesiteObservable:
    group: int
    elements: np.ndarray
    sites: List[int]
    name: str

    def __init__(self, group: int, elements: np.ndarray, sites: List[int], name: str):
        self.group = group
        assert elements.ndim == 2
        assert elements.shape[0] == elements.shape[1]
        self.elements = elements
        self.sites = sites
        self.name = name

    def to_toml_strs(self) -> List[str]:
        ret = []
        ret.append('name = "{}"'.format(self.name))
        ret.append("group = {}".format(self.group))
        ret.append("sites = {}".format(self.sites))
        dim = self.elements.shape[0]
        ret.append("dim = {}".format(dim))
        it = np.nditer(
            self.elements, flags=["multi_index"], op_flags=["readonly"], order="F"
        )
        ret.append('elements = """')
        while not it.finished:
            i = it.multi_index
            v = self.elements[i]
            if np.abs(v) == 0.0:
                it.iternext()
                continue
            line = ""
            for j in i:
                line += str(j) + " "
            line += " {} {}".format(np.real(v), np.imag(v))
            ret.append(line)
            it.iternext()
        ret.append('"""')
        return ret


class TwositeObservable:
    group: int
    bonds: List[Bond]
    elements: Optional[np.ndarray]
    ops: Optional[List[int]]
    name: str

    def __init__(
        self,
        group: int,
        bonds: List[Bond],
        *,
        elements: np.ndarray = None,
        ops: List[int] = None,
        name: str = ""
    ):
        self.group = group
        if elements is not None:
            dim = elements.ndim
            nbody = dim // 2
            assert nbody == 2
            assert dim == 2 * nbody
            for i in range(nbody):
                assert elements.shape[i] == elements.shape[i + nbody]
            self.elements = elements
            self.ops = None
        else:
            self.elements = None
            self.ops = ops
        self.bonds = bonds
        self.name = name

    def to_toml_strs(self) -> List[str]:
        ret = []
        ret.append('name = "{}"'.format(self.name))
        ret.append("group = {}".format(self.group))
        ret.append('bonds = """')
        for b in self.bonds:
            ret.append("{} {} {}".format(b.source_site, b.dx, b.dy))
        ret.append('"""')
        if self.elements is not None:
            dims = list(self.elements.shape[0:2])
            ret.append("dim = {}".format(dims))
            it = np.nditer(
                self.elements, flags=["multi_index"], op_flags=["readonly"], order="F"
            )
            ret.append('elements = """')
            while not it.finished:
                i = it.multi_index
                v = self.elements[i]
                if np.abs(v) == 0.0:
                    it.iternext()
                    continue
                line = ""
                for j in i:
                    line += str(j) + " "
                line += " {} {}".format(np.real(v), np.imag(v))
                ret.append(line)
                it.iternext()
            ret.append('"""')
        else:
            ret.append("ops = {}".format(self.ops))
        return ret

    def to_twosite_operators(self) -> List[NNOperator]:
        if self.elements is not None:
            return [NNOperator(bond, elements=self.elements) for bond in self.bonds]
        else:
            return [NNOperator(bond, ops=self.ops) for bond in self.bonds]


def make_evolution_onesite(
    hamiltonian: SiteOperator,
    graph: LatticeGraph,
    tau: float,
    result_cutoff: float = 1e-15,
) -> List[SiteOperator]:
    D, V = np.linalg.eigh(hamiltonian.elements)
    evo = np.einsum("il, l, jl -> ij", V, np.exp(-tau * D), V)
    return [SiteOperator(hamiltonian.site, evo)]


def make_evolution_twosite(
    hamiltonian: NNOperator,
    graph: LatticeGraph,
    tau: float,
    result_cutoff: float = 1e-15,
) -> List[NNOperator]:
    if hamiltonian.elements is None:
        raise NotImplementedError(
            "make_evolution for NNOperator with two onesite operators is not yet implemented."
        )
    dims = hamiltonian.elements.shape[0:2]
    H = hamiltonian.elements.reshape((dims[0] * dims[1], dims[0] * dims[1]))
    D, V = np.linalg.eigh(H)
    evo = np.reshape(
        np.dot(V, np.dot(np.diag(np.exp(-tau * D)), V.transpose())),
        hamiltonian.elements.shape,
    )
    bonds = graph.make_path(hamiltonian.bond)
    nhops = len(bonds)
    if nhops == 1:
        return [NNOperator(hamiltonian.bond, elements=evo)]

    mdofs = []
    unitcell = graph.unitcell
    for bond in bonds[1:]:
        mdofs.append(unitcell.sites[bond.source_site].phys_dim)
    nsites = nhops + 1
    nmids = nsites - 2
    index = [0]
    index += list(range(4, 4 + nmids))
    index += [1, 2]
    index += list(range(4 + nmids, 2 * nsites))
    index += [3]
    I = np.eye(np.prod(mdofs)).reshape(mdofs + mdofs)
    A = np.einsum("ijkl,...->ijkl...", evo, I).transpose(index)

    dofs = [dims[0]] + mdofs + [dims[1]]

    ret = []
    for bond in bonds[0:-1]:
        ndim = A.ndim
        nsites = ndim // 2
        index = [0, 1, nsites]
        index += list(range(2, nsites))
        index += list(range(nsites + 1, ndim))
        A = A.transpose(index)
        C = A.reshape((A.shape[0] * A.shape[1] * A.shape[nsites], -1))
        U, S, Vt = np.linalg.svd(C, full_matrices=False)
        U = np.dot(U, np.diag(S))
        B = U.reshape((A.shape[0], A.shape[1], dofs[0], -1))
        A = Vt.reshape([B.shape[3]] + dofs[2:] + dofs[1:])
        ret.append(NNOperator(bond, elements=B))
        dofs.pop(0)
    ret.append(NNOperator(bonds[-1], elements=A))
    return ret


def make_evolution(
    hamiltonian: Operator,
    graph: LatticeGraph,
    tau: float,
    result_cutoff: float = 1e-15,
) -> Union[List[SiteOperator], List[NNOperator]]:
    if isinstance(hamiltonian, SiteOperator):
        return make_evolution_onesite(hamiltonian, graph, tau, result_cutoff)
    else:
        return make_evolution_twosite(hamiltonian, graph, tau, result_cutoff)


class Model:
    param: Dict[str, Dict[str, Any]]
    parameter: Dict[str, Any]
    simple_tau: float
    full_tau: float
    correlation: Dict[str, Any]
    clength: Dict[str, Any]
    unitcell: Unitcell
    hamiltonians: List[Operator]
    graph: LatticeGraph
    onesites: List[OnesiteObservable]
    twobodies: List[TwositeObservable]
    simple_updates: List[Operator]
    full_updates: List[Operator]

    def __init__(self, param: MutableMapping, atol: float = 1e-15):
        param = lower_dict(param)
        self.param = param
        self.parameter = param["parameter"]
        self.simple_tau = self.parameter["simple_update"].get("tau", 0.01)
        self.full_tau = self.parameter["full_update"].get("tau", 0.01)
        self.correlation = param.get("correlation", None)
        self.clength = param.get("correlation_length", None)

        self.unitcell = Unitcell(param["tensor"])
        offset_x_min = offset_x_max = offset_y_min = offset_y_max = 0
        self.hamiltonians = []
        ham_as_onesite_obs = []
        ham_as_twosite_obs = []
        for ham in param["hamiltonian"]:
            if "sites" in ham and "bonds" in ham:
                raise RuntimeError("Both sites and bonds are defined in hamiltonian")
            dims: List[int] = ham["dim"]
            if isinstance(dims, int):
                dims = [dims]
            if len(dims) == 1:
                assert (
                    isinstance(dims[0], int) and dims[0] > 0
                ), "hamiltonian.dims should be list with one or two positive integer(s) but {} is given".format(
                    dims
                )
                elements = load_tensor(ham["elements"], dims + dims, atol=atol)
                assert is_hermite(elements)
                sites: List[int] = ham["sites"]
                if len(sites) == 0:
                    sites = list(range(self.unitcell.numsites()))
                for site in sites:
                    op = SiteOperator(site, elements)
                    self.hamiltonians.append(op)
                ham_as_onesite_obs.append(
                    OnesiteObservable(
                        0, sites=sites, elements=elements, name="hamiltonian"
                    )
                )
            elif len(dims) == 2:
                assert (
                    isinstance(dims[0], int)
                    and isinstance(dims[1], int)
                    and dims[0] > 0
                    and dims[1] > 0
                ), "hamiltonian.dims should be list with one or two positive integer(s) but {} is given".format(
                    dims
                )
                elements = load_tensor(ham["elements"], dims + dims, atol=atol)
                assert is_hermite(elements)
                bonds: List[Bond] = []
                for line in ham["bonds"].strip().splitlines():
                    b = parse_bond(line)
                    if b is None:
                        continue
                    ox, oy = self.unitcell.target_offset(b)
                    offset_x_min = min(ox, offset_x_min)
                    offset_x_max = max(ox, offset_x_max)
                    offset_y_min = min(oy, offset_y_min)
                    offset_y_max = max(oy, offset_y_max)
                    assert (
                        self.unitcell.sites[self.unitcell.source_site(b)].phys_dim
                        == elements.shape[0]
                    )
                    assert (
                        self.unitcell.sites[self.unitcell.target_site(b)].phys_dim
                        == elements.shape[1]
                    )
                    bonds.append(b)
                    op = NNOperator(b, elements=elements)
                    self.hamiltonians.append(op)
                ham_as_twosite_obs.append(
                    TwositeObservable(0, bonds, elements=elements, name="hamiltonian")
                )
            else:
                raise RuntimeError("dims should be a list with two elements")
        self.graph = LatticeGraph(
            self.unitcell, offset_x_min, offset_y_min, offset_x_max, offset_y_max
        )

        self.onesites = []
        self.twobodies = []
        observable = param.get("observable", {})
        has_zero_onesite = False
        for onesite in observable.get("onesite", []):
            name = onesite["name"]
            group = onesite["group"]
            if group == 0:
                has_zero_onesite = True
            sites = onesite["sites"]
            dim = onesite["dim"]
            elements = load_tensor(onesite["elements"], [dim, dim], atol=atol)
            one_obs = OnesiteObservable(
                group=group, elements=elements, sites=sites, name=name
            )
            self.onesites.append(one_obs)
        if not has_zero_onesite:
            for ham in ham_as_onesite_obs:
                self.onesites.append(ham)
        has_zero_twosite = False
        for twosite in observable.get("twosite", []):
            name = twosite["name"]
            group = twosite["group"]
            if group == 0:
                has_zero_twosite = True
            bonds = [
                parse_bond(line)
                for line in twosite["bonds"].strip().splitlines()
                if parse_bond(line) is not None
            ]
            if "elements" in twosite:
                dim = twosite["dim"]
                elements = load_tensor(twosite["elements"], dim + dim, atol=atol)
                two_obs = TwositeObservable(group, bonds, elements=elements, name=name)
            else:
                two_obs = TwositeObservable(group, bonds, ops=twosite["ops"], name=name)
            self.twobodies.append(two_obs)
        if not has_zero_twosite:
            for ham in ham_as_twosite_obs:
                self.twobodies.append(ham)

        self.simple_updates = []
        self.full_updates = []

        for ham in self.hamiltonians:
            for evo in make_evolution(ham, self.graph, self.simple_tau):
                self.simple_updates.append(evo)
            for evo in make_evolution(ham, self.graph, self.full_tau):
                self.full_updates.append(evo)

    def to_toml(self, f: TextIO):
        # parameter
        f.write("[parameter]\n")
        for tablename, table in self.parameter.items():
            f.write("[parameter.{}]\n".format(tablename))
            for k, v in table.items():
                f.write("{} = {}\n".format(k, value_to_str(v)))
        f.write("\n")

        # tensor
        f.write("[tensor]\n")
        f.write("L_sub = {}\n".format(self.param["tensor"]["l_sub"]))
        if "skew" in self.param["tensor"]:
            f.write("skew = {}\n".format(self.param["tensor"]["skew"]))
        for ucell in self.param["tensor"]["unitcell"]:
            f.write("[[tensor.unitcell]]\n")
            f.write("index = {}\n".format(ucell["index"]))
            f.write("physical_dim = {}\n".format(ucell["physical_dim"]))
            f.write("virtual_dim = {}\n".format(ucell["virtual_dim"]))
            if "initial_state" in ucell:
                f.write("initial_state = {}\n".format(ucell["initial_state"]))
            if "noise" in ucell:
                f.write("noise = {}\n".format(ucell["noise"]))
        f.write("\n")

        # observable
        f.write("[observable]\n")
        for onesite in self.onesites:
            f.write("[[observable.onesite]]\n")
            for line in onesite.to_toml_strs():
                f.write(line + "\n")
        f.write("\n")

        for twosite in self.twobodies:
            f.write("[[observable.twosite]]\n")
            for line in twosite.to_toml_strs():
                f.write(line + "\n")
        f.write("\n")

        # correlation
        if self.correlation is not None:
            f.write("[correlation]\n")
            f.write("r_max = {}\n".format(self.correlation["r_max"]))
            f.write("operators = [\n")
            for ops in self.correlation["operators"]:
                f.write("  {},\n".format(ops))
            f.write("]\n")
            f.write("\n")

        # correlation length
        if self.clength is not None:
            f.write("[correlation_length]\n")
            for k, v in self.clength.items():
                f.write("{} = {}\n".format(k, value_to_str(v)))

        f.write("[evolution]\n")
        for update in self.simple_updates:
            f.write("[[evolution.simple]]\n")
            for line in update.to_toml_strs(self.unitcell):
                f.write(line + "\n")
        for update in self.full_updates:
            f.write("[[evolution.full]]\n")
            for line in update.to_toml_strs(self.unitcell):
                f.write(line + "\n")


if __name__ == "__main__":
    import argparse
    import sys

    import toml

    parser = argparse.ArgumentParser(
        description="Input converter for TeNeS", add_help=True
    )

    parser.add_argument("input", help="Input TOML file")

    parser.add_argument(
        "-o", "--output", dest="output", default="input.toml", help="Output TOML file"
    )
    parser.add_argument(
        "-v", "--version", dest="version", action="version", version="1.3.4"
    )

    args = parser.parse_args()

    if args.input == args.output:
        print("The names of input and output are the same")
        sys.exit(1)

    param = toml.load(args.input)
    model = Model(param)

    with open(args.output, "w") as f:
        model.to_toml(f)
