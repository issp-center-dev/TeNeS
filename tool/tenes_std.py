from collections import namedtuple
from itertools import product
from typing import IO, List, Tuple

import numpy as np

import scipy.sparse as sparse


def all_positive(xs, v=0):
    return all(map(lambda x: x > v, xs))


Bond = namedtuple("Bond", ("source_site", "target_site", "offset_x", "offset_y"))


def parse_bond(line: str) -> Bond:
    words = line.split()
    source_site = int(words[0])
    target_site = int(words[1])
    offset_x = int(words[2])
    offset_y = int(words[3])
    return Bond(source_site, target_site, offset_x, offset_y)


def encode_bond(bond: Bond) -> str:
    return "{} {} {} {}".format(
        bond.source_site, bond.target_site, bond.offset_x, bond.offset_y
    )


def load_tensor(elements_str: str, dims: List[int], atol: float = 1e-15) -> np.ndarray:
    A_re = np.zeros(dims)
    A_im = np.zeros(dims)
    ndim = A_re.ndim
    for line in elements_str.strip().splitlines():
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
    Amat = A.reshape((np.prod(A.shape[:nsite]), np.prod(A.shape[nsite:])))
    return np.all(Amat.conjugate().transpose() == Amat)


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

    def __init__(self, tensor_dict: dict = None):
        self.phys_dim = None
        self.virtual_dim = [None] * 4
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

    def __init__(self, lat_dict: dict = None):
        self.L = [1] * 2
        self.sites = [None]
        if lat_dict is not None:
            self.load_dict(lat_dict)

    def load_dict(self, lat_dict: dict):
        L = lat_dict["L_sub"]
        if isinstance(L, int):
            self.L = [L, L]
        else:
            self.L = L

        self.skew = lat_dict.get("skew", 0)

        N = self.L[0] * self.L[1]
        self.sites = [None] * N
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

    def bond_displacement(self, bond: Bond) -> Tuple:
        x, y = self.index2coord(bond.source_site)
        X, Y = self.index2coord(bond.target_site)
        X += (self.L[0]) * bond.offset_x + self.skew * bond.offset_y
        Y += self.L[1] * bond.offset_y
        dx = X - x
        dy = Y - y
        return dx, dy

    def bond_direction(self, bond: Bond) -> int:
        dx, dy = self.bond_displacement(bond)
        nhop = abs(dx) + abs(dy)
        assert nhop == 1
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
        assert 0 <= direction < 4
        x, y = self.index2coord(source)
        X, Y = x, y
        offset_x = 0
        offset_y = 0
        if direction == 0:
            X = x - 1
            if X < 0:
                X = self.L[0] - 1
                offset_x = -1
        elif direction == 1:
            Y = y + 1
            if Y == self.L[1]:
                Y = 0
                offset_y = 1
                X = x - self.skew
                if X < 0:
                    X += self.L[0]
                    offset_x = -1
                elif X >= self.L[0]:
                    X -= self.L[0]
                    offset_x = 1
        elif direction == 2:
            X = x + 1
            if X == self.L[0]:
                X = 0
                offset_x = 1
        else:
            Y = y - 1
            if Y < 0:
                Y = self.L[1] - 1
                offset_y = -1
                X = x + self.skew
                if X < 0:
                    X += self.L[0]
                    offset_x = -1
                elif X >= self.L[0]:
                    X -= self.L[0]
                    offset_x = 1

        return Bond(source, self.coord2index(X, Y), offset_x, offset_y)

    def check(self):
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

    def index2coord(self, index: int) -> Tuple:
        return index % self.L[0], index // self.L[0]

    def coord2index(self, x: int, y: int) -> int:
        return x + self.L[0] * y

    def neighbor(self, index: int, direction: int) -> int:
        x, y = self.index2coord(index)
        if direction == 0:
            x = (x - 1 + self.L[0]) % self.L[0]
        elif direction == 1:
            if self.skew == 0 or y < self.L[1]:
                y = (y + 1) % self.L[1]
            else:
                y = 0
                x = (x + self.skew + self.L[0]) % self.L[0]
        elif direction == 2:
            x = (x + 1) % self.L[0]
        elif direction == 4:
            if self.skew == 0 or y == 0:
                y = (y - 1 + self.L[1]) % self.L[1]
            else:
                y = 0
                x = (x - self.skew + self.L[0]) % self.L[0]
        else:
            msg = "ERROR: given direction is {}, but must be 0, 1, 2, or 3".format(
                direction
            )
            raise RuntimeError(msg)
        return self.coord2index(x, y)


class LatticeGraph:
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
                    D = unitcell.bond_dim(bond)
                    if D == 1:
                        continue
                    if not offset_x_min <= bond.offset_x + iox <= offset_x_max:
                        continue
                    if not offset_y_min <= bond.offset_y + ioy <= offset_y_max:
                        continue
                    source_site = self.graph_site(bond.source_site, iox, ioy)
                    target_site = self.graph_site(
                        bond.target_site, iox + bond.offset_x, ioy + bond.offset_y
                    )
                    weight = 1.0 / D / D
                    if bond.offset_x == bond.offset_y == 0:
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

    def graph_coords(self, index: int):
        assert 0 <= index <= self.N

        index, localsite = divmod(index, self.unitcell.numsites())
        offset_y, offset_x = divmod(index, self.offset_Lx)
        return localsite, offset_x, offset_y

    def make_path(self, bond: Bond) -> List[Bond]:
        target_index = self.graph_site(bond.target_site, bond.offset_x, bond.offset_y)
        path_pred = self.path_pred[self.graph_site(bond.source_site)]
        targets = [target_index]
        target_index = path_pred[target_index]
        while target_index >= 0:
            targets.append(target_index)
            target_index = path_pred[target_index]

        bonds = []
        source_site, source_offset_x, source_offset_y = self.graph_coords(targets.pop())
        while len(targets) > 0:
            target_site, target_offset_x, target_offset_y = self.graph_coords(
                targets.pop()
            )
            dx = target_offset_x - source_offset_x
            dy = target_offset_y - source_offset_y
            b = Bond(source_site, target_site, dx, dy)
            bonds.append(b)
            source_site = target_site
            source_offset_x = target_offset_x
            source_offset_y = target_offset_y
        return bonds


class TwoBodyOperator:
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


class OnsiteObservable:
    def __init__(self, group: int, elements: np.ndarray, sites: List[int]):
        self.group = group
        assert elements.ndim == 2
        assert elements.shape[0] == elements.shape[1]
        self.elements = elements
        self.sites = sites

    def to_toml_strs(self) -> List[str]:
        ret = []
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


class TwoBodyObservables:
    def __init__(
        self,
        group: int,
        bonds: List[Bond],
        *,
        elements: np.ndarray = None,
        ops: List[int] = None
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

    def to_toml_strs(self) -> List[str]:
        ret = []
        ret.append("group = {}".format(self.group))
        ret.append('bonds = """')
        for b in self.bonds:
            ret.append(
                "{} {} {} {}".format(
                    b.source_site, b.target_site, b.offset_x, b.offset_y
                )
            )
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

    def to_twosite_operators(self) -> List[TwoBodyOperator]:
        if self.elements is not None:
            return [
                TwoBodyOperator(bond, elements=self.elements) for bond in self.bonds
            ]
        else:
            return [TwoBodyOperator(bond, ops=self.ops) for bond in self.bonds]


def make_evolution(
    hamiltonian: TwoBodyOperator,
    graph: LatticeGraph,
    tau: float,
    result_cutoff: float = 1e-15,
) -> List[TwoBodyOperator]:
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
        return [TwoBodyOperator(hamiltonian.bond, elements=evo)]

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
        ret.append(TwoBodyOperator(bond, elements=B))
        dofs.pop(0)
    ret.append(TwoBodyOperator(bonds[-1], elements=A))
    return ret


class Model:
    def __init__(self, param: dict, atol: float = 1e-15):
        self.param = param
        self.parameter = param["parameter"]
        self.simple_tau = self.parameter["simple_update"].get("tau", 0.01)
        self.full_tau = self.parameter["full_update"].get("tau", 0.01)
        self.correlation = param.get("correlation", None)

        self.unitcell = Unitcell(param["tensor"])
        offset_x_min = offset_x_max = offset_y_min = offset_y_max = 0
        self.energy_obs = []
        self.hamiltonians = []
        for ham in param["hamiltonian"]:
            dims = ham["dims"]
            assert (
                len(dims) == 2
                and isinstance(dims[0], int)
                and isinstance(dims[1], int)
                and dims[0] > 0
                and dims[1] > 0
            ), "hamiltonian.dims should be list with two positive integers but {} is given".format(
                dims
            )
            elements = load_tensor(ham["elements"], dims + dims, atol=atol)
            assert is_hermite(elements)
            bonds = []
            for line in ham["bonds"].strip().splitlines():
                b = parse_bond(line)
                offset_x_min = min(b.offset_x, offset_x_min)
                offset_x_max = max(b.offset_x, offset_x_max)
                offset_y_min = min(b.offset_y, offset_y_min)
                offset_y_max = max(b.offset_y, offset_y_max)
                assert self.unitcell.sites[b.source_site].phys_dim == elements.shape[0]
                assert self.unitcell.sites[b.target_site].phys_dim == elements.shape[1]
                bonds.append(b)
                op = TwoBodyOperator(b, elements=elements)
                self.hamiltonians.append(op)
            enes = TwoBodyObservables(0, bonds, elements=elements)
            self.energy_obs.append(enes)
        self.graph = LatticeGraph(
            self.unitcell, offset_x_min, offset_y_min, offset_x_max, offset_y_max
        )

        self.onesites = []
        self.twobodies = []
        if "observable" in param:
            observable = param["observable"]
            for onesite in observable.get("onesite", []):
                group = onesite["group"]
                sites = onesite["sites"]
                dim = onesite["dim"]
                elements = load_tensor(onesite["elements"], [dim, dim], atol=atol)
                obs = OnsiteObservable(group, elements, sites)
                self.onesites.append(obs)

            for twosite in observable.get("twosite", []):
                group = twosite["group"]
                bonds = [
                    parse_bond(line) for line in twosite["bonds"].strip().splitlines()
                ]
                if "elements" in twosite:
                    dim = twosite["dim"]
                    elements = load_tensor(twosite["elements"], dim + dim, atol=atol)
                    obs = TwoBodyObservables(group, bonds, elements=elements)
                else:
                    obs = TwoBodyObservables(group, bonds, ops=twosite["ops"])
                self.twobodies.append(obs)

        self.simple_updates = []
        self.full_updates = []

        for ham in self.hamiltonians:
            for evo in make_evolution(ham, self.graph, self.simple_tau):
                self.simple_updates.append(evo)
            for evo in make_evolution(ham, self.graph, self.full_tau):
                self.full_updates.append(evo)

    def to_toml(self, f: IO, *, no_hamiltonian_observe=False):
        # parameter
        f.write("[parameter]\n")
        for tablename, table in self.parameter.items():
            f.write("[parameter.{}]\n".format(tablename))
            for k, v in table.items():
                if isinstance(v, str):
                    f.write("{} = '{}'\n".format(k, v))
                elif isinstance(v, bool):
                    f.write("{} = {}\n".format(k, "true" if v else "false"))
                else:
                    f.write("{} = {}\n".format(k, v))
        f.write("\n")

        # tensor
        f.write("[tensor]\n")
        f.write("L_sub = {}\n".format(self.param["tensor"]["L_sub"]))
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

        if not no_hamiltonian_observe:
            for obs in self.energy_obs:
                f.write("[[observable.twosite]]\n")
                for line in obs.to_toml_strs():
                    f.write(line + "\n")

        for twosite in self.twobodies:
            if twosite.group == 0 and not no_hamiltonian_observe:
                continue
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
        "--no-hamiltonian-observe",
        action="store_true",
        help="Do not overwrite twosite observable with group=0 by bond Hamiltonian",
    )

    args = parser.parse_args()
    if args.input == args.output:
        print("The names of input and output are the same")
        sys.exit(1)

    param = toml.load(args.input)
    model = Model(param)

    with open(args.output, "w") as f:
        model.to_toml(f, no_hamiltonian_observe=args.no_hamiltonian_observe)

