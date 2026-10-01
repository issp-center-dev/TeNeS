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

"""Regenerate the bosonic baseline of test_fermion_chain.py (contract item 4).

Writes three bosonic std.toml inputs (bosonic_*.toml, each with bonds that
tenes_std decomposes into two or three nearest-neighbour gates) and the
input.toml that tool/tenes_std.py renders from them (bosonic_*.input.toml).

    cd test/python/data/fermion_chain && python3 make_bosonic_baseline.py

The stored files were made on 2026-09-30 with tool/tenes_std.py at
1a823e62, before task T1 of
docs/superpowers/plans/2026-09-30-fermion-longrange-hamiltonian.md.  They
are the reference that T1 must not change: rerun this only to change the
bosonic output on purpose.  ``--outdir DIR`` writes to DIR instead, to
compare a fresh rendering with the stored files without touching them.
"""

import argparse
import io
import os
import sys

import numpy as np
import toml

HERE = os.path.dirname(os.path.abspath(__file__))


def load_tool():
    sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "tool"))
    import tenes_std

    return tenes_std


def spin_ops(S2):
    """S = S2 / 2 spin matrices (Sz, S+, S-) in the |m = S, S-1, ...> basis."""
    d = S2 + 1
    S = 0.5 * S2
    m = np.array([S - i for i in range(d)])
    Sz = np.diag(m)
    Sp = np.zeros((d, d))
    for i in range(1, d):
        Sp[i - 1, i] = np.sqrt(S * (S + 1) - m[i] * (m[i] + 1))
    return Sz, Sp, Sp.T


def bond_elements(M, d):
    """Matrix M[(o1 o2), (i1 i2)] -> 'i1 i2 o1 o2 re im' lines."""
    A = M.reshape(d, d, d, d).transpose(2, 3, 0, 1)
    lines = []
    for idx in np.ndindex(*A.shape):
        v = A[idx]
        if abs(v) > 1e-15:
            lines.append(
                "{} {} {} {} {!r} {!r}".format(
                    *idx, float(np.real(v)), float(np.imag(v))
                )
            )
    return "\n".join(lines) + "\n"


def site_elements(M):
    lines = []
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[j, i]
            if abs(v) > 1e-15:
                lines.append(
                    "{} {} {!r} {!r}".format(i, j, float(np.real(v)), float(np.imag(v)))
                )
    return "\n".join(lines) + "\n"


def heisenberg(S2, J=1.0):
    Sz, Sp, Sm = spin_ops(S2)
    return J * (np.kron(Sz, Sz) + 0.5 * (np.kron(Sp, Sm) + np.kron(Sm, Sp)))


def inputs():
    ret = {}

    # spin 1/2 J1-J2-J3 on a 2x2 cell: (1,1)/(1,-1) and (2,0)/(0,2) bonds
    # are two-hop chains
    h = bond_elements(heisenberg(1), 2)
    ret["bosonic_heisenberg_2x2"] = {
        "parameter": {
            "general": {"is_real": True},
            "simple_update": {"num_step": 10, "tau": [0.05, 0.01]},
            "full_update": {"num_step": 0, "tau": 0.01},
        },
        "tensor": {
            "L_sub": [2, 2],
            "unitcell": [{"index": [], "physical_dim": 2, "virtual_dim": 3}],
        },
        "hamiltonian": [
            {
                "dim": [2, 2],
                "bonds": "0 1 0\n1 1 0\n2 1 0\n3 1 0\n0 0 1\n1 0 1\n2 0 1\n3 0 1\n",
                "elements": h,
            },
            {
                "dim": [2, 2],
                "bonds": "0 1 1\n1 1 1\n0 1 -1\n1 1 -1\n",
                "elements": bond_elements(0.4 * heisenberg(1), 2),
            },
            {
                "dim": [2, 2],
                "bonds": "0 2 0\n3 0 2\n",
                "elements": bond_elements(0.2 * heisenberg(1), 2),
            },
            {
                "dim": [2],
                "sites": [],
                "elements": site_elements(-0.3 * spin_ops(1)[0]),
            },
        ],
    }

    # spin 1/2 on a 3x3 cell with non-uniform bond dimensions (the path
    # weights prefer wide bonds) and three-hop bonds
    h1 = heisenberg(1)
    unitcell = []
    for i in range(9):
        vd = [3, 3, 3, 3]
        if i in (0, 4):
            vd = [2, 3, 2, 3]
        unitcell.append({"index": [i], "physical_dim": 2, "virtual_dim": vd})
    # make the horizontal bonds consistent: site i's +x leg meets site
    # (i+1)'s -x leg
    for i in range(9):
        x, y = i % 3, i // 3
        right = ((x + 1) % 3) + 3 * y
        unitcell[right]["virtual_dim"][0] = unitcell[i]["virtual_dim"][2]
        up = x + 3 * ((y + 1) % 3)
        unitcell[up]["virtual_dim"][3] = unitcell[i]["virtual_dim"][1]
    ret["bosonic_mixed_d_3x3"] = {
        "parameter": {
            "general": {"is_real": True},
            "simple_update": {"num_step": 10, "tau": 0.02},
            "full_update": {"num_step": 0, "tau": 0.01},
        },
        "tensor": {"L_sub": [3, 3], "unitcell": unitcell},
        "hamiltonian": [
            {
                "dim": [2, 2],
                "bonds": "".join("{} 1 0\n{} 0 1\n".format(i, i) for i in range(9)),
                "elements": bond_elements(h1, 2),
            },
            {
                "dim": [2, 2],
                "bonds": "0 2 1\n4 -2 1\n8 1 2\n2 0 -2\n",
                "elements": bond_elements(0.3 * h1, 2),
            },
        ],
    }

    # complex Hamiltonian (Dzyaloshinskii-Moriya z term) in real-time mode
    Sz, Sp, Sm = spin_ops(1)
    dm = 0.5j * (np.kron(Sp, Sm) - np.kron(Sm, Sp))
    ret["bosonic_complex_te_2x2"] = {
        "parameter": {
            "general": {"mode": "time evolution", "is_real": False},
            "simple_update": {"num_step": 10, "tau": 0.01},
            "full_update": {"num_step": 0, "tau": 0.01},
        },
        "tensor": {
            "L_sub": [2, 2],
            "unitcell": [{"index": [], "physical_dim": 2, "virtual_dim": 2}],
        },
        "hamiltonian": [
            {
                "dim": [2, 2],
                "bonds": "0 1 0\n0 0 1\n3 -1 0\n3 0 -1\n",
                "elements": bond_elements(heisenberg(1) + 0.3 * dm, 2),
            },
            {
                "dim": [2, 2],
                "bonds": "0 -1 1\n1 0 -2\n",
                "elements": bond_elements(0.5 * dm + 0.2 * heisenberg(1), 2),
            },
        ],
    }
    return ret


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--outdir", default=HERE)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    tool = load_tool()
    for name, param in inputs().items():
        std_path = os.path.join(args.outdir, name + ".toml")
        with open(std_path, "w") as f:
            toml.dump(param, f)
        model = tool.Model(toml.load(std_path))
        buf = io.StringIO()
        model.to_toml(buf)
        with open(os.path.join(args.outdir, name + ".input.toml"), "w") as f:
            f.write(buf.getvalue())
        print("wrote", name)


if __name__ == "__main__":
    main()
