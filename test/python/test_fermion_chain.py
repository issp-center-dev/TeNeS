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

"""Long-range fermion Hamiltonian bonds in tenes_std (task T1).

Behaviour contract: "振る舞い契約書(T1)" in
docs/superpowers/plans/2026-09-30-fermion-longrange-hamiltonian.md, design in
docs/superpowers/specs/2026-09-30-fermion-longrange-hamiltonian-design.md.

How the reference is built (and why it is independent of tenes_std)
-------------------------------------------------------------------
A long-range bond (s, dx, dy) is decomposed along the path
s = x_0, x_1, ..., x_L = t that ``LatticeGraph.make_path`` picks.  The
reference for the whole gate chain is the operator O_st = exp(-tau H_st),
with identity on the intermediate sites, written as a plain matrix in the
Fock space whose modes are ordered along the PATH:

    |n_0 n_1 ... n_L> = (c^dag_{x_0})^{n_0} (c^dag_{x_1})^{n_1} ... |0>

(for Hubbard every site carries the two modes up, dn in this order, and the
local index is i = n_up + 2 n_dn).  H_st is assembled here from
Jordan-Wigner creation/annihilation matrices (``creation`` below), and O_st
is its matrix exponential (scipy.linalg.expm).  The Jordan-Wigner strings
that c_t picks up when it passes the intermediate modes come out of the
Fock-space algebra by themselves; no sign formula of the design (section
4.1) or of tenes_std is used anywhere in the reference.  The two-site input
handed to tenes_std is the same construction with L = 1, i.e. the matrix
element definition of design section 3 (source-first ordered basis).

The gate chain is composed as plain Kronecker embeddings in this path
order: gate k acts on path positions (k, k+1), and the fat chi leg a gate
leaves on position k+1 is treated as the local space of that site until the
next gate consumes it.  This is legitimate because every gate is parity
even once each chi index is given its parity (contract item 2 checks
exactly that), and an even operator on two adjacent positions of an ordered
Fock space is the bare Kronecker embedding: the strings from the modes on
its left cancel in pairs and the modes on its right are untouched.
"""

import copy
import io
import os
import re
import sys
import time
import tracemalloc

import numpy as np
import pytest
import scipy.linalg
import toml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tool"))

import tenes_std  # noqa: E402

DATA_DIR = os.path.join(HERE, "data", "fermion_chain")

ATOL = 1e-12


def popcount(x):
    # int.bit_count() needs Python 3.10; CI runs 3.9.
    return bin(x).count("1")


# ---------------------------------------------------------------------------
# Independent Fock-space reference
# ---------------------------------------------------------------------------


def creation(mode, nmodes):
    """c^dag_mode in the ordered Fock basis, numpy convention mat[out, in].

    Basis state g has bit m equal to n_m, and
    |g> = (c^dag_0)^{n_0} (c^dag_1)^{n_1} ... |0>.  Bringing c^dag_mode to
    its place passes every occupied mode k < mode, one minus sign each.
    """
    dim = 1 << nmodes
    mat = np.zeros((dim, dim))
    below = (1 << mode) - 1
    for g in range(dim):
        if (g >> mode) & 1:
            continue
        sign = -1.0 if popcount(g & below) % 2 else 1.0
        mat[g | (1 << mode), g] = sign
    return mat


class SiteKind:
    """A fermion site: nspin modes; local index i is the occupation bit
    pattern occ[i] (bit s = mode s), by default occ[i] = i."""

    def __init__(self, name, nspin, parity, occ=None):
        self.name = name
        self.nspin = nspin
        self.d = 1 << nspin
        self.occ = list(range(self.d)) if occ is None else list(occ)
        assert sorted(self.occ) == list(range(self.d))
        self.parity = list(parity)
        # the declared parity table must be the particle-number parity of
        # the occupation the local index stands for, or the fixture itself
        # is wrong
        assert self.parity == [popcount(b) % 2 for b in self.occ]


SPINLESS = SiteKind("spinless", 1, [0, 1])
# tool/tenes_simple.py HubbardModel: i = n_up + 2 n_dn, |up dn> = c^dag_up
# c^dag_dn |0>, parity [0, 1, 1, 0]
HUBBARD = SiteKind("hubbard", 2, [0, 1, 1, 0])
# spinless with the local labels swapped: index 0 is the occupied state,
# parity table [1, 0]
SPINLESS_FLIPPED = SiteKind("flipped", 1, [1, 0], occ=[1, 0])


class PathSpace:
    """Fock space of L + 1 sites in path order (source first, target last).

    `kinds` is one SiteKind per path position (or a single SiteKind with
    `nsites`, for a path of identical sites).  Site k owns the modes
    off_k ... off_k + nspin_k - 1, in path order.
    """

    def __init__(self, kinds, nsites=None):
        if isinstance(kinds, SiteKind):
            kinds = [kinds] * nsites
        self.kinds = list(kinds)
        self.nsites = len(self.kinds)
        self.offsets = []
        nm = 0
        for kind in self.kinds:
            self.offsets.append(nm)
            nm += kind.nspin
        self.nmodes = nm
        self.cdag = [
            [creation(off + s, nm) for s in range(kind.nspin)]
            for off, kind in zip(self.offsets, self.kinds)
        ]
        # every matrix is real, so c = (c^dag)^T
        self.c = [[m.T.copy() for m in row] for row in self.cdag]
        self.eye = np.eye(1 << nm)

    def n(self, k, s=None):
        if s is not None:
            return self.cdag[k][s] @ self.c[k][s]
        return sum(self.cdag[k][s] @ self.c[k][s] for s in range(len(self.cdag[k])))

    def majorana(self, k, s):
        return self.cdag[k][s] + self.c[k][s]

    def to_tensor(self, mat):
        """mat[g_out, g_in] -> T[i_0, ..., i_L, o_0, ..., o_L].

        The Fock state of the local indices (i_0, ..., i_L) is
        g = sum_k occ_k[i_k] << off_k.
        """
        dims = [kind.d for kind in self.kinds]
        L = self.nsites
        gvec = np.array(
            [
                sum(
                    kind.occ[i] << off
                    for i, kind, off in zip(idx, self.kinds, self.offsets)
                )
                for idx in np.ndindex(*dims)
            ]
        )
        T = mat[np.ix_(gvec, gvec)].reshape(dims + dims)  # [o..., i...]
        return T.transpose(list(range(L, 2 * L)) + list(range(L)))


# Each operator is a two-site Hamiltonian between path position 0 (source)
# and the last position (target), built from creation/annihilation matrices.
T_HOP = 1.0
V_NN = 0.8


def h_spinless_hop(sp):
    s, t = 0, sp.nsites - 1
    return -T_HOP * (sp.cdag[s][0] @ sp.c[t][0] + sp.cdag[t][0] @ sp.c[s][0])


def h_spinless_nn(sp):
    s, t = 0, sp.nsites - 1
    return V_NN * sp.n(s) @ sp.n(t)


def h_hubbard_hop_v(sp):
    s, t = 0, sp.nsites - 1
    h = V_NN * sp.n(s) @ sp.n(t)
    for spin in range(2):
        h = h - T_HOP * (
            sp.cdag[s][spin] @ sp.c[t][spin] + sp.cdag[t][spin] @ sp.c[s][spin]
        )
    return h


def h_spinless_product(sp):
    # a sum of one-site terms: exp(-tau H) is a product, operator Schmidt
    # rank 1 across (s | t), even channel only
    s, t = 0, sp.nsites - 1
    return -0.7 * sp.n(s) + 0.4 * sp.n(t)


def h_hubbard_product(sp):
    s, t = 0, sp.nsites - 1
    return (
        2.0 * sp.n(s, 0) @ sp.n(s, 1)
        - 0.7 * sp.n(s)
        + 1.5 * sp.n(t, 0) @ sp.n(t, 1)
        + 0.4 * sp.n(t)
    )


def h_majorana(sp):
    # H = i gamma_s gamma_t (gamma = c + c^dag of the first mode) is
    # Hermitian with H^2 = 1, so exp(-tau H) at tau = i pi / 2 is exactly
    # -i H = gamma_s gamma_t: operator Schmidt rank 1 in the ODD channel.
    s, t = 0, sp.nsites - 1
    return 1j * sp.majorana(s, 0) @ sp.majorana(t, 0)


# name -> (site kind, Hamiltonian builder, tau)
OPERATORS = {
    "spinless_hop": (SPINLESS, h_spinless_hop, 0.1),
    "spinless_nn": (SPINLESS, h_spinless_nn, 0.1),
    "hubbard_hopV": (HUBBARD, h_hubbard_hop_v, 0.1),
}

# Extra operators built for contract item 3 (chi == d); see
# test_some_gate_has_chi_equal_to_d_with_a_non_physical_parity_table.
CONSTRUCTED_OPERATORS = {
    "spinless_product": (SPINLESS, h_spinless_product, 0.1),
    "hubbard_product": (HUBBARD, h_hubbard_product, 0.1),
    "spinless_majorana": (SPINLESS, h_majorana, 0.5j * np.pi),
    "hubbard_majorana": (HUBBARD, h_majorana, 0.5j * np.pi),
}

ALL_OPERATORS = dict(OPERATORS)
ALL_OPERATORS.update(CONSTRUCTED_OPERATORS)

# design section 4.3
DISPLACEMENTS = [
    (2, 0),
    (0, 2),
    (-2, 0),
    (0, -2),
    (1, 1),
    (1, -1),
    (-1, 1),
    (-1, -1),
    (2, 1),
    (-2, 1),
]

# cell name -> (L_sub, source site, displacements)
# 4x4: every site on a path of at most three hops is a different unit-cell
#      site.
# 2x2: the same unit-cell site comes back on the path: (2, 0) and (0, -2)
#      end on the source site itself, and in (2, 1) the second intermediate
#      site is the source site again.
CELLS = {
    "4x4": ([4, 4], 5, DISPLACEMENTS),
    "2x2": ([2, 2], 0, [(2, 0), (0, -2), (1, -1), (-1, 1), (2, 1), (-2, 1)]),
}

MAIN_CASES = [
    (cell, disp, op) for cell in CELLS for disp in CELLS[cell][2] for op in OPERATORS
]
CONSTRUCTED_CASES = [
    ("4x4", disp, op) for disp in [(2, 0), (-1, 1)] for op in CONSTRUCTED_OPERATORS
]


def case_id(case):
    cell, (dx, dy), op = case
    return "{}-({},{})-{}".format(cell, dx, dy, op)


def make_unitcell(kind, l_sub, virtual_dim=2):
    return tenes_std.Unitcell(
        {
            "l_sub": l_sub,
            "unitcell": [
                {
                    "index": [],
                    "physical_dim": kind.d,
                    "virtual_dim": virtual_dim,
                    "parity": kind.parity,
                }
            ],
        }
    )


def make_graph(unitcell, bonds):
    """The graph tenes_std.Model would build for these Hamiltonian bonds."""
    ox_min = oy_min = ox_max = oy_max = 0
    for b in bonds:
        ox, oy = unitcell.target_offset(b)
        ox_min, ox_max = min(ox, ox_min), max(ox, ox_max)
        oy_min, oy_max = min(oy, oy_min), max(oy, oy_max)
    return tenes_std.LatticeGraph(unitcell, ox_min, oy_min, ox_max, oy_max)


def two_site_hamiltonian(kind, builder, target_kind=None):
    """H[i_s, i_t, o_s, o_t] = <o_s o_t| H |i_s i_t> (design section 3)."""
    kinds = [kind, kind if target_kind is None else target_kind]
    sp = PathSpace(kinds)
    return sp.to_tensor(builder(sp))


def path_reference(kind, builder, tau, nsites):
    """O_st (x) 1_mid in the path-ordered Fock basis, as
    T[i_0..i_L, o_0..o_L] = <o| exp(-tau H_st) |i>."""
    sp = PathSpace(kind, nsites)  # kind: a SiteKind, or one per position
    return sp.to_tensor(scipy.linalg.expm(-tau * builder(sp)))


def unsigned_embedding(kind, builder, tau, nsites):
    """evo[i_s, i_t, o_s, o_t] * prod_m delta(i_m, o_m), with no sign at
    all: what the bosonic decomposition composes to."""
    evo = path_reference(kind, builder, tau, 2)
    d = kind.d
    nmid = nsites - 2
    T = evo
    for _ in range(nmid):
        T = np.multiply.outer(T, np.eye(d))
    # axes now: i_s, i_t, o_s, o_t, (i_m, o_m) * nmid
    perm = [0] + [4 + 2 * m for m in range(nmid)] + [1, 2]
    perm += [5 + 2 * m for m in range(nmid)] + [3]
    return T.transpose(perm)


def compose(gates, nsites, d_end):
    """Compose the gate chain in path order.

    X[i_0..i_L, c_0..c_L] starts as the identity; gate k contracts its
    (in1, in2) with the current legs c_k, c_{k+1} and puts (out1, out2) in
    their place.  Returns X, or raises AssertionError on a dimension
    mismatch between consecutive gates.
    """
    L = nsites
    d = d_end
    X = np.eye(d**L).reshape([d] * (2 * L))
    for k, G in enumerate(gates):
        a, b = L + k, L + k + 1
        assert G.ndim == 4, "gate {} has rank {}".format(k, G.ndim)
        assert (
            X.shape[a] == G.shape[0] and X.shape[b] == G.shape[1]
        ), "gate {} input legs {} do not match the current legs {}".format(
            k, G.shape[:2], (X.shape[a], X.shape[b])
        )
        Y = np.tensordot(X, G, axes=([a, b], [0, 1]))
        X = np.moveaxis(Y, [Y.ndim - 2, Y.ndim - 1], [a, b])
    return X


class Case:
    """One (cell, displacement, operator) case, with the path tenes_std's
    own graph picks and the gates produced with fermion=True."""

    def __init__(self, cell, disp, opname):
        l_sub, source, _ = CELLS[cell]
        self.kind, self.builder, self.tau = ALL_OPERATORS[opname]
        self.unitcell = make_unitcell(self.kind, l_sub)
        self.bond = tenes_std.Bond(source, disp[0], disp[1])
        self.graph = make_graph(self.unitcell, [self.bond])
        self.path = self.graph.make_path(self.bond)
        self.nsites = len(self.path) + 1
        uc = self.unitcell
        # unit-cell site at each path position
        self.sites = [int(b.source_site) for b in self.path]
        self.sites.append(int(uc.target_site(self.path[-1])))
        self.legs = [uc.bond_direction(b) for b in self.path]
        self.H2 = two_site_hamiltonian(self.kind, self.builder)
        self.hamiltonian = tenes_std.NNOperator(self.bond, elements=self.H2)

    def describe(self):
        steps = ", ".join(
            "({}: {:+d},{:+d} leg {})".format(int(b.source_site), b.dx, b.dy, leg)
            for b, leg in zip(self.path, self.legs)
        )
        return "path of {} hops [{}] through unit-cell sites {}".format(
            len(self.path), steps, self.sites
        )

    def phys_parity(self, position):
        return list(self.unitcell.sites[self.sites[position]].parity)

    def fermion_gates(self):
        return tenes_std.make_evolution_twosite(
            self.hamiltonian, self.graph, self.tau, fermion=True
        )


def infer_parity_chain(case, gates):
    """Walk the chain with the solver's inference rule (design 6.1).

    Returns (tables, problems): tables[k] is the out2 parity table of gate
    k (a list of 0/1) for the non-final gates; problems lists every
    violation of contract item 2 as a string.
    """
    problems = []
    tables = []
    current = case.phys_parity(0)  # parity table of in1 of gate 0
    for k, G in enumerate(gates):
        last = k == len(gates) - 1
        p_in1 = current
        p_in2 = case.phys_parity(k + 1)
        p_out1 = case.phys_parity(k)
        if G.shape[0] != len(p_in1):
            problems.append(
                "gate {}: in1 has {} indices but the leg it consumes has {}".format(
                    k, G.shape[0], len(p_in1)
                )
            )
            return tables, problems
        if G.shape[1] != len(p_in2) or G.shape[2] != len(p_out1):
            problems.append(
                "gate {}: physical legs in2/out1 have shape {}, expected {}".format(
                    k, G.shape[1:3], (len(p_in2), len(p_out1))
                )
            )
            return tables, problems
        nz = np.argwhere(G != 0)
        # every in1 index must be used (a chi index no gate reads is dead)
        used_in1 = set(int(a) for a in nz[:, 0])
        dead_in1 = sorted(set(range(G.shape[0])) - used_in1)
        if dead_in1:
            problems.append(
                "gate {}: in1 indices {} have no nonzero element".format(k, dead_in1)
            )
        if last:
            p_out2 = case.phys_parity(k + 1)
            if G.shape[3] != len(p_out2):
                problems.append(
                    "last gate: out2 has {} indices, target has {}".format(
                        G.shape[3], len(p_out2)
                    )
                )
                return tables, problems
            odd = [
                tuple(int(x) for x in idx)
                for idx in nz
                if (p_in1[idx[0]] + p_in2[idx[1]] + p_out1[idx[2]] + p_out2[idx[3]]) % 2
            ]
            if odd:
                problems.append(
                    "last gate: {} nonzero elements are parity odd, e.g. {}".format(
                        len(odd), odd[:3]
                    )
                )
        else:
            seen = [set() for _ in range(G.shape[3])]
            for a, b, c, e in nz:
                seen[e].add((p_in1[a] + p_in2[b] + p_out1[c]) % 2)
            table = []
            for e, s in enumerate(seen):
                if len(s) == 0:
                    problems.append(
                        "gate {}: out2 index {} has no nonzero element".format(k, e)
                    )
                    table.append(0)
                elif len(s) == 2:
                    problems.append(
                        "gate {}: out2 index {} has mixed parity".format(k, e)
                    )
                    table.append(0)
                else:
                    table.append(s.pop())
            tables.append(table)
            current = table
    return tables, problems


# ---------------------------------------------------------------------------
# Contract item 1: exact composition
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("case", MAIN_CASES + CONSTRUCTED_CASES, ids=case_id)
def test_chain_composes_to_the_path_ordered_operator(case):
    c = Case(*case)
    assert len(c.path) >= 2, c.describe()
    gates = c.fermion_gates()
    info = c.describe()

    assert len(gates) == len(c.path), info
    for g, b in zip(gates, c.path):
        # the chain follows the path, the near end of each hop is the source
        assert (int(g.bond.source_site), int(g.bond.dx), int(g.bond.dy)) == (
            int(b.source_site),
            int(b.dx),
            int(b.dy),
        ), info
        assert g.elements is not None, info

    got = compose([g.elements for g in gates], c.nsites, c.kind.d)
    ref = path_reference(c.kind, c.builder, c.tau, c.nsites)
    assert got.shape == ref.shape, info
    err = np.max(np.abs(got - ref))
    assert err <= ATOL, "{}: max |chain - reference| = {:.3e}".format(info, err)


@pytest.mark.parametrize("case", MAIN_CASES, ids=case_id)
def test_reference_is_not_the_unsigned_embedding_for_odd_channels(case):
    """Self-check of the reference: where the operator has odd channels and
    an intermediate site can be odd, the correct chain is NOT the bosonic
    (unsigned) one, by a margin far above ATOL.  For the purely even
    operator n_s n_t the two coincide (documented here so that nobody reads
    the nn cases as sign tests of the odd channel)."""
    c = Case(*case)
    ref = path_reference(c.kind, c.builder, c.tau, c.nsites)
    plain = unsigned_embedding(c.kind, c.builder, c.tau, c.nsites)
    diff = np.max(np.abs(ref - plain))
    if case[2] == "spinless_nn":
        assert diff <= ATOL, c.describe()
    else:
        assert diff > 1e-3, "{}: diff {:.3e}".format(c.describe(), diff)


def test_case_list_covers_left_and_up_steps_and_revisited_sites():
    """Contract item 1: at least one case steps with source_leg 0 (towards
    -x) and one with source_leg 1 (towards +y); the list also has to keep
    a path whose target is its source's unit-cell site and one whose
    intermediate site is its source's unit-cell site."""
    legs = set()
    target_is_source = []
    mid_revisits = []
    lines = []
    for case in MAIN_CASES:
        c = Case(*case)
        legs.update(c.legs)
        lines.append("{}: {}".format(case_id(case), c.describe()))
        if c.sites[-1] == c.sites[0]:
            target_is_source.append(case_id(case))
        if any(s in (c.sites[0], c.sites[-1]) for s in c.sites[1:-1]):
            mid_revisits.append(case_id(case))
    summary = "\n".join(lines)
    assert 0 in legs, "no step with source_leg 0:\n" + summary
    assert 1 in legs, "no step with source_leg 1:\n" + summary
    assert target_is_source, summary
    assert mid_revisits, summary


def test_make_evolution_dispatches_the_fermion_flag():
    """make_evolution(..., fermion=True) must reach the fermion branch."""
    for case in [("4x4", (-1, 1), "spinless_hop"), ("2x2", (2, 1), "hubbard_hopV")]:
        c = Case(*case)
        gates = tenes_std.make_evolution(c.hamiltonian, c.graph, c.tau, fermion=True)
        got = compose([g.elements for g in gates], c.nsites, c.kind.d)
        ref = path_reference(c.kind, c.builder, c.tau, c.nsites)
        err = np.max(np.abs(got - ref))
        assert err <= ATOL, "{}: {:.3e}".format(c.describe(), err)


# ---------------------------------------------------------------------------
# Contract item 2: every chi index has one parity
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("case", MAIN_CASES + CONSTRUCTED_CASES, ids=case_id)
def test_chi_parity_is_unique_and_chain_is_even(case):
    c = Case(*case)
    gates = c.fermion_gates()
    _, problems = infer_parity_chain(c, [g.elements for g in gates])
    assert not problems, "{}:\n  {}".format(c.describe(), "\n  ".join(problems))


# ---------------------------------------------------------------------------
# Contract item 3: chi == d with a parity table that is not the physical one
# ---------------------------------------------------------------------------


def test_some_gate_has_chi_equal_to_d_with_a_non_physical_parity_table():
    """Contract item 3.

    None of the contract-1 operators can do it: the fat leg a gate leaves
    on an intermediate site m has chi = d_m * r, r being the operator
    Schmidt rank of O_st across (s | t) (the delta on m carries d_m
    independent copies), so chi == d_m needs r == 1, i.e. an O_st that is a
    product.  CONSTRUCTED_OPERATORS supplies such products: even ones
    (exp of one-site terms; the chi table is then the physical table in
    whatever order the implementation lists its blocks) and odd ones
    (gamma_s gamma_t, whose chi table is the physical table flipped).  At
    least one of them has to give a chi table that differs position by
    position from the physical table of the site it sits on.
    """
    found = []
    seen = []
    for case in MAIN_CASES + CONSTRUCTED_CASES:
        c = Case(*case)
        gates = [g.elements for g in c.fermion_gates()]
        tables, problems = infer_parity_chain(c, gates)
        assert not problems, "{}: {}".format(case_id(case), problems)
        for k, table in enumerate(tables):
            phys = c.phys_parity(k + 1)
            seen.append(
                "{} gate {}: chi = {}, d = {}, table {} vs physical {}".format(
                    case_id(case), k, len(table), len(phys), table, phys
                )
            )
            if len(table) == len(phys) and table != phys:
                found.append(seen[-1])
    assert found, "no gate with chi == d and a non-physical table:\n" + "\n".join(seen)


# ---------------------------------------------------------------------------
# Contract items 1 and 2 on paths whose sites differ in kind
#
# With one kind of site on the whole path, every position has the same
# physical dimension and parity table, so taking the ledger of the wrong
# position (source for an intermediate site, one position off, source and
# target swapped) goes unnoticed.  Here the unit cell mixes spinless
# (d = 2, [0, 1]), Hubbard (d = 4, [0, 1, 1, 0]) and a spinless site with
# swapped labels (d = 2, [1, 0]: local index 0 is the occupied state).
# ---------------------------------------------------------------------------

S_, H_, F_ = SPINLESS, HUBBARD, SPINLESS_FLIPPED


def _checker(even, odd):
    return lambda x, y: even if (x + y) % 2 == 0 else odd


def _three_kinds(x, y):
    return [S_, H_, F_][(x + y) % 3]


def _two_by_two(x, y):
    return {(0, 0): H_, (1, 0): S_, (0, 1): F_, (1, 1): S_}[(x, y)]


MIXED_DISPLACEMENTS = [(2, 0), (-2, 0), (0, 2), (1, -1), (-1, 1), (2, 1), (-2, 1)]

# cell name -> (L_sub, source site, kind of the site at (x, y), displacements)
# 4x4-SH: spinless source, Hubbard nearest neighbours: S-H-S, S-H-S-H
# 4x4-HS: the other way round: H-S-H, H-S-H-S
# 4x4-SHF: kind by (x + y) mod 3, source (1, 1) is flipped spinless:
#          e.g. (2, 0) is F-S-H, (2, 1) is F-S-H-F
# 2x2-mixed: sites 0..3 are H, S, F, S; (2, 0) is H-S-H back on the source
#          site, (2, 1) is H-S-H-F with the source site again in the middle
MIXED_CELLS = {
    "4x4-SH": ([4, 4], 5, _checker(S_, H_), MIXED_DISPLACEMENTS),
    "4x4-HS": ([4, 4], 5, _checker(H_, S_), MIXED_DISPLACEMENTS),
    "4x4-SHF": ([4, 4], 5, _three_kinds, DISPLACEMENTS),
    "2x2-mixed": (
        [2, 2],
        0,
        _two_by_two,
        [(2, 0), (-2, 0), (0, -2), (1, 1), (2, 1), (-2, 1)],
    ),
}


def h_hop0_nn(sp):
    return h_spinless_hop(sp) + h_spinless_nn(sp)


# Mode 0 of each site (the up mode of a Hubbard site) hops; the density
# interaction couples the total occupations.
MIXED_OPERATORS = {
    "hop0": (h_spinless_hop, 0.1),
    "nn": (h_spinless_nn, 0.1),
    "hop0_nn": (h_hop0_nn, 0.1),
}

MIXED_CASES = [
    (cell, disp, op)
    for cell in MIXED_CELLS
    for disp in MIXED_CELLS[cell][3]
    for op in MIXED_OPERATORS
]


class MixedCase(Case):
    """Like Case, but every unit-cell site has its own SiteKind."""

    def __init__(self, cell, disp, opname):
        l_sub, source, kind_of, _ = MIXED_CELLS[cell]
        self.builder, self.tau = MIXED_OPERATORS[opname]
        Lx, Ly = l_sub
        self.site_kinds = [kind_of(i % Lx, i // Lx) for i in range(Lx * Ly)]
        self.unitcell = tenes_std.Unitcell(
            {
                "l_sub": l_sub,
                "unitcell": [
                    {
                        "index": [i],
                        "physical_dim": k.d,
                        "virtual_dim": 2,
                        "parity": k.parity,
                    }
                    for i, k in enumerate(self.site_kinds)
                ],
            }
        )
        self.bond = tenes_std.Bond(source, disp[0], disp[1])
        self.graph = make_graph(self.unitcell, [self.bond])
        self.path = self.graph.make_path(self.bond)
        self.nsites = len(self.path) + 1
        uc = self.unitcell
        self.sites = [int(b.source_site) for b in self.path]
        self.sites.append(int(uc.target_site(self.path[-1])))
        self.legs = [uc.bond_direction(b) for b in self.path]
        self.kinds = [self.site_kinds[i] for i in self.sites]
        self.H2 = two_site_hamiltonian(self.kinds[0], self.builder, self.kinds[-1])
        self.hamiltonian = tenes_std.NNOperator(self.bond, elements=self.H2)

    def describe(self):
        return "{}; kinds along the path {}".format(
            Case.describe(self), [k.name for k in self.kinds]
        )


def unsigned_embedding_mixed(c):
    """evo (x) 1_mid with no sign, for per-position dimensions."""
    evo = path_reference(c.kinds[:: len(c.kinds) - 1], c.builder, c.tau, 2)
    T = evo
    mids = c.kinds[1:-1]
    for k in mids:
        T = np.multiply.outer(T, np.eye(k.d))
    nmid = len(mids)
    perm = [0] + [4 + 2 * m for m in range(nmid)] + [1, 2]
    perm += [5 + 2 * m for m in range(nmid)] + [3]
    return T.transpose(perm)


@pytest.mark.parametrize("case", MIXED_CASES, ids=case_id)
def test_mixed_chain_composes_to_the_path_ordered_operator(case):
    c = MixedCase(*case)
    info = c.describe()
    gates = c.fermion_gates()
    assert len(gates) == len(c.path), info
    for g, b in zip(gates, c.path):
        assert (int(g.bond.source_site), int(g.bond.dx), int(g.bond.dy)) == (
            int(b.source_site),
            int(b.dx),
            int(b.dy),
        ), info
    ref = path_reference(c.kinds, c.builder, c.tau, c.nsites)
    got = compose_chain([g.elements for g in gates])
    assert got.shape == ref.shape, "{}: {} vs {}".format(info, got.shape, ref.shape)
    err = np.max(np.abs(got - ref))
    assert err <= ATOL, "{}: max |chain - reference| = {:.3e}".format(info, err)


@pytest.mark.parametrize("case", MIXED_CASES, ids=case_id)
def test_mixed_chi_parity_is_unique_and_chain_is_even(case):
    c = MixedCase(*case)
    gates = c.fermion_gates()
    _, problems = infer_parity_chain(c, [g.elements for g in gates])
    assert not problems, "{}:\n  {}".format(c.describe(), "\n  ".join(problems))


@pytest.mark.parametrize("case", MIXED_CASES, ids=case_id)
def test_mixed_reference_is_not_the_unsigned_embedding(case):
    """Self-check of the mixed reference, as for the uniform cases: the odd
    channel of mode-0 hopping sees the string of every intermediate site;
    the pure density interaction does not."""
    c = MixedCase(*case)
    ref = path_reference(c.kinds, c.builder, c.tau, c.nsites)
    plain = unsigned_embedding_mixed(c)
    diff = np.max(np.abs(ref - plain))
    if case[2] == "nn":
        assert diff <= ATOL, c.describe()
    else:
        assert diff > 1e-3, "{}: diff {:.3e}".format(c.describe(), diff)


def test_mixed_case_list_covers_the_kind_patterns():
    """The mixed list keeps: S-H-S and H-S-H two-hop paths, a path whose
    source, intermediate and target kinds all differ, a three-hop path
    whose two intermediate sites differ in kind, a flipped site at the
    source, in the middle and at the target, steps with source_leg 0 and
    1, and unit-cell sites that come back on the path."""
    names = []
    legs = set()
    revisit = False
    lines = []
    for case in MIXED_CASES:
        c = MixedCase(*case)
        k = [x.name for x in c.kinds]
        names.append(k)
        legs.update(c.legs)
        lines.append("{}: {}".format(case_id(case), c.describe()))
        if len(set(c.sites)) < len(c.sites):
            revisit = True
    summary = "\n".join(lines)
    assert ["spinless", "hubbard", "spinless"] in names, summary
    assert ["hubbard", "spinless", "hubbard"] in names, summary
    assert any(len(k) == 3 and len(set(k)) == 3 for k in names), summary
    assert any(len(k) == 4 and k[1] != k[2] for k in names), summary
    assert any(k[0] == "flipped" for k in names), summary
    assert any("flipped" in k[1:-1] for k in names), summary
    assert any(k[-1] == "flipped" for k in names), summary
    assert {0, 1} <= legs, summary
    assert revisit, summary


# ---------------------------------------------------------------------------
# Contract item 4: one hop and bosonic inputs are unchanged
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("disp", [(1, 0), (0, 1), (-1, 0), (0, -1)])
@pytest.mark.parametrize("opname", list(OPERATORS))
def test_one_hop_fermion_gate_equals_the_bosonic_one(disp, opname):
    kind, builder, tau = OPERATORS[opname]
    uc = make_unitcell(kind, [4, 4])
    bond = tenes_std.Bond(5, disp[0], disp[1])
    graph = make_graph(uc, [bond])
    ham = tenes_std.NNOperator(bond, elements=two_site_hamiltonian(kind, builder))
    for func in (tenes_std.make_evolution_twosite, tenes_std.make_evolution):
        fer = func(ham, graph, tau, fermion=True)
        bos = func(ham, graph, tau, fermion=False)
        assert len(fer) == len(bos) == 1
        assert (fer[0].bond.source_site, fer[0].bond.dx, fer[0].bond.dy) == (
            bos[0].bond.source_site,
            bos[0].bond.dx,
            bos[0].bond.dy,
        )
        assert fer[0].elements.dtype == bos[0].elements.dtype
        assert np.array_equal(fer[0].elements, bos[0].elements)


BOSONIC_INPUTS = [
    "bosonic_heisenberg_2x2.toml",
    "bosonic_mixed_d_3x3.toml",
    "bosonic_complex_te_2x2.toml",
]


def _render(input_path):
    model = tenes_std.Model(toml.load(input_path))
    buf = io.StringIO()
    model.to_toml(buf)
    return buf.getvalue()


def _stored(input_path):
    with open(input_path[: -len(".toml")] + ".input.toml") as f:
        return f.read()


def split_at_evolution(text):
    """(everything before the [evolution] table, the [evolution] table)."""
    lines = text.splitlines(keepends=True)
    heads = [i for i, line in enumerate(lines) if line.strip() == "[evolution]"]
    assert len(heads) == 1, heads
    return "".join(lines[: heads[0]]), "".join(lines[heads[0] :])


class CellGeometry:
    """Unit cell of an input.toml: site -> physical dim, bond targets."""

    def __init__(self, text):
        tensor = toml.loads(text)["tensor"]
        L = tensor["L_sub"]
        self.Lx, self.Ly = (L, L) if isinstance(L, int) else L
        self.skew = tensor.get("skew", 0)
        n = self.Lx * self.Ly
        self.phys = [None] * n
        for uc in tensor["unitcell"]:
            index = uc["index"]
            if isinstance(index, int):
                index = [index]
            if len(index) == 0:
                index = range(n)
            for i in index:
                self.phys[i] = uc["physical_dim"]

    def target(self, source, leg):
        """Nearest neighbour of `source` across leg 0 (-x), 1 (+y), 2 (+x),
        3 (-y), through the periodic boundary with skew."""
        dx, dy = {0: (-1, 0), 1: (0, 1), 2: (1, 0), 3: (0, -1)}[leg]
        x, y = source % self.Lx + dx, source // self.Lx + dy
        oy = y // self.Ly
        x -= oy * self.skew
        return x % self.Lx + self.Lx * (y % self.Ly)


def dense_elements(entry):
    A = np.zeros(entry["dimensions"], dtype=complex)
    for line in entry["elements"].strip().splitlines():
        words = line.split()
        if not words:
            continue
        idx = tuple(int(w) for w in words[:-2])
        A[idx] = float(words[-2]) + 1j * float(words[-1])
    return A


def gate_key(entry):
    """Everything of an [[evolution.*]] entry except the matrix elements."""
    return {k: v for k, v in entry.items() if k != "elements"}


def segment_gates(entries, geom):
    """Split an [[evolution.*]] list into one-site gates, nearest-neighbour
    gates and gate chains of long-range bonds, by the output order alone.

    A two-site gate whose in1 is the physical leg of its source starts a
    unit; the unit goes on while the out2 leg is not the physical leg of
    the target (a fat chi leg), and the next gate must then start at that
    target and consume that leg.  (A chi that happens to equal the physical
    dimension would end a chain early; test_bosonic_baseline_segmentation_
    matches_the_hamiltonian rules that out for the stored data.)
    """
    units = []
    chain = None
    for n, e in enumerate(entries):
        if "site" in e:
            assert chain is None, "one-site gate {} inside a gate chain".format(n)
            units.append([e])
            continue
        dims = e["dimensions"]
        s = e["source_site"]
        t = geom.target(s, e["source_leg"])
        assert dims[1] == geom.phys[t] and dims[2] == geom.phys[s], (n, dims)
        if chain is None:
            assert dims[0] == geom.phys[s], "gate {} starts with a fat leg".format(n)
            chain = [e]
        else:
            prev = chain[-1]
            assert s == geom.target(prev["source_site"], prev["source_leg"]), n
            assert dims[0] == prev["dimensions"][3], n
            chain.append(e)
        if dims[3] == geom.phys[t]:
            units.append(chain)
            chain = None
    assert chain is None, "the last gate chain is not closed"
    return units


def compose_chain(gates):
    """Contract a gate chain back into the multi-site operator
    X[i_0..i_L, o_0..o_L]; independent of the SVD gauge of the chi legs."""
    dims = [gates[0].shape[0]] + [G.shape[1] for G in gates]
    L = len(dims)
    X = np.eye(int(np.prod(dims))).reshape(dims + dims)
    for k, G in enumerate(gates):
        a, b = L + k, L + k + 1
        Y = np.tensordot(X, G, axes=([a, b], [0, 1]))
        X = np.moveaxis(Y, [Y.ndim - 2, Y.ndim - 1], [a, b])
    return X


@pytest.mark.parametrize("name", BOSONIC_INPUTS)
def test_bosonic_output_matches_the_stored_baseline(name):
    """Contract item 4, bosonic part.

    test/python/data/fermion_chain/bosonic_*.input.toml were rendered by
    tool/tenes_std.py before T1 (make_bosonic_baseline.py in the same
    directory).  The comparison does not depend on the LAPACK build:
      * everything before [evolution] (parameter, tensor, observables)
        must be the same text;
      * every gate must have the same group, site / source_site /
        source_leg and dimensions, in the same order;
      * one-site and nearest-neighbour gates must agree element by element
        (atol 1e-12);
      * the gates of a long-range bond's chain are contracted back into the
        multi-site operator before comparing, since the SVD may pick other
        signs or rotations of the chi legs.
    """
    path = os.path.join(DATA_DIR, name)
    stored = _stored(path)
    new = _render(path)
    head_stored, evo_stored = split_at_evolution(stored)
    head_new, evo_new = split_at_evolution(new)
    assert head_new == head_stored

    geom = CellGeometry(stored)
    evo_stored = toml.loads(evo_stored)["evolution"]
    evo_new = toml.loads(evo_new)["evolution"]
    assert sorted(evo_new) == sorted(evo_stored)
    nchains = 0
    for kind in evo_stored:
        old_list, new_list = evo_stored[kind], evo_new[kind]
        assert [gate_key(e) for e in new_list] == [gate_key(e) for e in old_list]
        old_units = segment_gates(old_list, geom)
        new_units = segment_gates(new_list, geom)
        assert [len(u) for u in new_units] == [len(u) for u in old_units]
        for uo, un in zip(old_units, new_units):
            where = "{} {}".format(kind, gate_key(uo[0]))
            if len(uo) == 1:
                A, B = dense_elements(uo[0]), dense_elements(un[0])
            else:
                nchains += 1
                A = compose_chain([dense_elements(e) for e in uo])
                B = compose_chain([dense_elements(e) for e in un])
            err = np.max(np.abs(A - B))
            assert err <= ATOL, "{}: {:.3e}".format(where, err)
    assert nchains > 0


@pytest.mark.parametrize("name", BOSONIC_INPUTS)
def test_bosonic_baseline_segmentation_matches_the_hamiltonian(name):
    """The chains segment_gates finds in the stored data are exactly one
    per Hamiltonian bond, as long as its make_path (one gate per one-site
    term and per nearest-neighbour bond)."""
    path = os.path.join(DATA_DIR, name)
    stored = _stored(path)
    geom = CellGeometry(stored)
    evo = toml.loads(split_at_evolution(stored)[1])["evolution"]
    model = tenes_std.Model(toml.load(path))
    per_group = []
    for ham in model.hamiltonians:
        if isinstance(ham, tenes_std.SiteOperator):
            per_group.append(1)
        else:
            per_group.append(len(model.graph.make_path(ham.bond)))
    for kind, taus in (("simple", model.simple_tau), ("full", model.full_tau)):
        units = segment_gates(evo[kind], geom)
        assert [len(u) for u in units] == per_group * len(taus)


# Bosonic long-range bonds whose path crosses sites of different physical
# dimension.  L_sub = [4, 1] with skew 1 (no site is its own neighbour);
# (dx, dy) = (2, 0) and (3, 0) from site 0 run 0 -> 1 -> 2 (-> 3).
BOSONIC_MIXED_DIM_CASES = {
    "2/3/2 two hops": ([2, 3, 2, 5], (2, 0), False),
    "2/3/4 two hops": ([2, 3, 4, 5], (2, 0), False),
    "2/3/3/2 three hops": ([2, 3, 3, 2], (3, 0), False),
    "3/2/4/2 three hops": ([3, 2, 4, 2], (3, 0), False),
    "2/3/3/2 three hops, complex": ([2, 3, 3, 2], (3, 0), True),
}


@pytest.mark.parametrize("name", list(BOSONIC_MIXED_DIM_CASES))
def test_bosonic_chain_over_sites_of_different_dimensions(name):
    """A bosonic (fermion=False) long-range gate chain over sites of
    different physical dimensions composes to evo (x) 1_mid, with no sign.

    evo is exp(-tau H) of a random Hermitian H, computed here with
    scipy.linalg.expm on H[(i_s i_t), (o_s o_t)]."""
    dims, disp, is_complex = BOSONIC_MIXED_DIM_CASES[name]
    uc = tenes_std.Unitcell(
        {
            "l_sub": [4, 1],
            "skew": 1,
            "unitcell": [
                {"index": [i], "physical_dim": d, "virtual_dim": 2}
                for i, d in enumerate(dims)
            ],
        }
    )
    bond = tenes_std.Bond(0, disp[0], disp[1])
    graph = make_graph(uc, [bond])
    path = graph.make_path(bond)
    sites = [int(b.source_site) for b in path] + [int(uc.target_site(path[-1]))]
    pdims = [dims[i] for i in sites]
    info = "path through sites {} with dimensions {}".format(sites, pdims)
    assert len(path) == abs(disp[0]), info
    assert len(set(pdims)) > 1, info

    ds, dt = pdims[0], pdims[-1]
    rng = np.random.default_rng(7)
    M = rng.normal(size=(ds * dt, ds * dt))
    if is_complex:
        M = M + 1j * rng.normal(size=(ds * dt, ds * dt))
    Hmat = 0.5 * (M + M.conj().T)
    ham = tenes_std.NNOperator(bond, elements=Hmat.reshape(ds, dt, ds, dt))
    tau = 0.1
    gates = tenes_std.make_evolution_twosite(ham, graph, tau, fermion=False)
    assert len(gates) == len(path), info

    evo = scipy.linalg.expm(-tau * Hmat).reshape(ds, dt, ds, dt)
    T = evo
    for d in pdims[1:-1]:
        T = np.multiply.outer(T, np.eye(d))
    nmid = len(pdims) - 2
    perm = [0] + [4 + 2 * m for m in range(nmid)] + [1, 2]
    perm += [5 + 2 * m for m in range(nmid)] + [3]
    ref = T.transpose(perm)

    got = compose_chain([g.elements for g in gates])
    assert got.shape == ref.shape, info
    err = np.max(np.abs(got - ref))
    assert err <= ATOL, "{}: {:.3e}".format(info, err)


def test_bosonic_baseline_inputs_contain_long_range_bonds():
    """The baseline comparison must exercise the SVD decomposition."""
    hops = []
    for name in BOSONIC_INPUTS:
        model = tenes_std.Model(toml.load(os.path.join(DATA_DIR, name)))
        assert not model.parameter.get("general", {}).get("fermion", False)
        for ham in model.hamiltonians:
            if isinstance(ham, tenes_std.NNOperator):
                hops.append(len(model.graph.make_path(ham.bond)))
    assert max(hops) >= 3, hops
    assert 2 in hops, hops


# ---------------------------------------------------------------------------
# Contract items 4-6 through Model: the fermion flag reaches the decomposition
# ---------------------------------------------------------------------------


def elements_string(A):
    lines = []
    for idx in np.ndindex(*A.shape):
        v = A[idx]
        if v != 0:
            lines.append(
                " ".join(str(i) for i in idx)
                + " {!r} {!r}".format(float(np.real(v)), float(np.imag(v)))
            )
    return "\n".join(lines)


def fermion_model_input(kind, builder, l_sub, bonds, parameter=None):
    """A fermion-mode std.toml dict with one two-site Hamiltonian entry."""
    H2 = two_site_hamiltonian(kind, builder)
    param = {
        "parameter": {
            "general": {"fermion": True},
            "simple_update": {"tau": 0.1, "num_step": 10},
        },
        "tensor": {
            "l_sub": l_sub,
            "unitcell": [
                {
                    "index": [],
                    "physical_dim": kind.d,
                    "virtual_dim": 2,
                    "parity": kind.parity,
                }
            ],
        },
        "hamiltonian": [
            {
                "dim": [kind.d, kind.d],
                "bonds": "".join("{} {} {}\n".format(*b) for b in bonds),
                "elements": elements_string(H2),
            }
        ],
    }
    if parameter:
        for table, entries in parameter.items():
            param["parameter"].setdefault(table, {}).update(entries)
    return param


@pytest.mark.parametrize(
    "l_sub, bond, opname",
    [
        ([4, 4], (5, -2, 1), "spinless_hop"),
        ([4, 4], (5, 0, 2), "hubbard_hopV"),
        ([2, 2], (0, 2, 0), "spinless_hop"),
        ([2, 2], (1, -1, 1), "hubbard_hopV"),
    ],
)
def test_model_passes_the_fermion_flag(l_sub, bond, opname):
    """Model with parameter.general.fermion = true decomposes a long-range
    bond with the fermion rule: its simple-update gates compose to the
    path-ordered reference."""
    kind, builder, tau = OPERATORS[opname]
    param = fermion_model_input(kind, builder, l_sub, [bond])
    model = tenes_std.Model(param)
    b = tenes_std.Bond(*bond)
    path = model.graph.make_path(b)
    gates = [g for g in model.simple_updates if g.group == 0]
    assert len(gates) == len(path) >= 2
    got = compose([g.elements for g in gates], len(path) + 1, kind.d)
    ref = path_reference(kind, builder, tau, len(path) + 1)
    err = np.max(np.abs(got - ref))
    assert err <= ATOL, "bond {} path {}: {:.3e}".format(
        bond, [(int(p.source_site), p.dx, p.dy) for p in path], err
    )


def test_model_nearest_neighbour_fermion_gates_are_the_bosonic_ones():
    kind, builder, tau = OPERATORS["hubbard_hopV"]
    param = fermion_model_input(kind, builder, [2, 2], [(0, 1, 0), (0, 0, 1)])
    model = tenes_std.Model(param)
    uc = model.unitcell
    ham = tenes_std.NNOperator(
        tenes_std.Bond(0, 1, 0), elements=two_site_hamiltonian(kind, builder)
    )
    expect = tenes_std.make_evolution_twosite(ham, model.graph, tau)
    gates = [g for g in model.simple_updates if g.group == 0]
    assert len(gates) == 2
    assert uc.bond_direction(gates[0].bond) == 2
    assert np.array_equal(gates[0].elements, expect[0].elements)


# ---------------------------------------------------------------------------
# Contract item 5: input checks
# ---------------------------------------------------------------------------


def spinless_long_range_input(parameter=None, l_sub=(4, 4), bonds=((5, 2, 0),)):
    kind, builder, _ = OPERATORS["spinless_hop"]
    return fermion_model_input(kind, builder, list(l_sub), list(bonds), parameter)


class TestFermionLongRangeInputChecks:
    def test_long_range_bond_is_accepted_without_full_update(self):
        model = tenes_std.Model(spinless_long_range_input())
        assert len([g for g in model.simple_updates if g.group == 0]) == 2

    def test_long_range_bond_is_accepted_with_zero_full_update_steps(self):
        param = spinless_long_range_input({"full_update": {"num_step": 0}})
        tenes_std.Model(param)

    def test_long_range_bond_is_accepted_with_zero_full_update_step_list(self):
        param = spinless_long_range_input(
            {"full_update": {"num_step": [0, 0], "tau": [0.1, 0.05]}}
        )
        tenes_std.Model(param)

    def test_long_range_bond_with_full_update_is_rejected(self):
        param = spinless_long_range_input({"full_update": {"num_step": 1}})
        with pytest.raises(RuntimeError, match="(?i)full"):
            tenes_std.Model(param)

    def test_long_range_bond_with_a_positive_full_update_step_in_a_list(self):
        param = spinless_long_range_input(
            {"full_update": {"num_step": [0, 3], "tau": [0.1, 0.05]}}
        )
        with pytest.raises(RuntimeError, match="(?i)full"):
            tenes_std.Model(param)

    def test_nearest_neighbour_fermion_bonds_keep_the_full_update(self):
        param = spinless_long_range_input(
            {"full_update": {"num_step": 5}}, bonds=((5, 1, 0), (5, 0, 1))
        )
        model = tenes_std.Model(param)
        assert len(model.full_updates) == 2

    def test_bosonic_long_range_bond_keeps_the_full_update(self):
        param = spinless_long_range_input({"full_update": {"num_step": 5}})
        del param["parameter"]["general"]["fermion"]
        for site in param["tensor"]["unitcell"]:
            del site["parity"]
        model = tenes_std.Model(param)
        assert len(model.full_updates) == 2

    def test_missing_parity_is_still_rejected_next_to_a_long_range_bond(self):
        param = spinless_long_range_input()
        param["tensor"]["unitcell"] = [
            {"index": [i], "physical_dim": 2, "virtual_dim": 2, "parity": [0, 1]}
            for i in range(16)
            if i != 6
        ] + [{"index": [6], "physical_dim": 2, "virtual_dim": 2}]
        with pytest.raises(RuntimeError, match="parity"):
            tenes_std.Model(param)

    def test_self_neighbour_cell_is_still_rejected_next_to_a_long_range_bond(self):
        # L_sub = [2, 1] with skew 0: every site is its own vertical neighbour
        param = spinless_long_range_input(l_sub=(2, 1), bonds=((0, 2, 0),))
        with pytest.raises(RuntimeError, match="own nearest neighbour"):
            tenes_std.Model(param)

    def test_multisite_is_still_rejected_next_to_a_long_range_bond(self):
        param = spinless_long_range_input()
        param["observable"] = {
            "multisite": [
                {
                    "name": "three",
                    "group": 0,
                    "multisites": "0 1 0 1 1\n",
                    "ops": [0, 0, 0],
                }
            ]
        }
        with pytest.raises(RuntimeError, match="multisite"):
            tenes_std.Model(param)


# ---------------------------------------------------------------------------
# Design section 1: a fermion Hamiltonian bond outside the 4x4 measurement
# window (|dx| > 3 or |dy| > 3) cannot have its energy measured, so it is
# refused -- also when the input defines its own group-0 two-site observable
# and tenes_std therefore does not add the Hamiltonian as one (then nothing
# downstream would refuse it).
# ---------------------------------------------------------------------------


def with_user_twosite_observable(param):
    """Give the input its own group-0 two-site observable (n_s n_t on a
    nearest-neighbour bond), so the Hamiltonian is not added as one."""
    param["observable"] = {
        "twosite": [
            {
                "name": "nn_user",
                "group": 0,
                "bonds": "5 1 0\n",
                "dim": [2, 2],
                "elements": "1 1 1 1 1.0 0.0",
            }
        ]
    }
    return param


def window_input(bond, user_observable):
    param = spinless_long_range_input(bonds=(bond,))
    if user_observable:
        with_user_twosite_observable(param)
    return param


def displacement_pattern(dx, dy):
    return r"\(\s*{}\s*,\s*{}\s*\)".format(dx, dy)


class TestFermionMeasurementWindow:
    @pytest.mark.parametrize("user_observable", [True, False], ids=["user", "auto"])
    @pytest.mark.parametrize("dx, dy", [(4, 0), (0, -4), (4, 1)])
    def test_bond_outside_the_window_is_rejected(self, dx, dy, user_observable):
        param = window_input((5, dx, dy), user_observable)
        with pytest.raises(RuntimeError) as excinfo:
            tenes_std.Model(param)
        msg = str(excinfo.value)
        assert "4x4" in msg, msg
        assert re.search(displacement_pattern(dx, dy), msg), msg
        assert re.search(r"\b5\b", msg), msg  # the source site

    @pytest.mark.parametrize("user_observable", [True, False], ids=["user", "auto"])
    @pytest.mark.parametrize("dx, dy", [(3, 0), (-3, 3)])
    def test_bond_on_the_window_edge_is_accepted(self, dx, dy, user_observable):
        param = window_input((5, dx, dy), user_observable)
        model = tenes_std.Model(param)
        path = model.graph.make_path(tenes_std.Bond(5, dx, dy))
        assert len([g for g in model.simple_updates if g.group == 0]) == len(path)
        names = [t.name for t in model.twobodies]
        if user_observable:
            # the fixture really suppresses the automatic energy observable
            assert names == ["nn_user"], names
        else:
            assert "bond_hamiltonian" in names, names

    @pytest.mark.parametrize("user_observable", [True, False], ids=["user", "auto"])
    def test_bosonic_bond_outside_the_window_is_accepted(self, user_observable):
        param = window_input((5, 4, 0), user_observable)
        del param["parameter"]["general"]["fermion"]
        for site in param["tensor"]["unitcell"]:
            del site["parity"]
        model = tenes_std.Model(param)
        path = model.graph.make_path(tenes_std.Bond(5, 4, 0))
        assert len(path) == 4
        assert len([g for g in model.simple_updates if g.group == 0]) == 4


# ---------------------------------------------------------------------------
# PR #118 review, item 1: a fermion-mode Hamiltonian bond term must be parity
# even.  The block SVD of the long-range decomposition only keeps the
# parity-diagonal blocks, so an odd or mixed term would otherwise be
# replaced silently by the even part of exp(-tau H); a nearest-neighbour odd
# term would reach tenes unchecked.  Elements no larger than tenes_std's
# atol (1e-15, the default of Model) do not count.
# ---------------------------------------------------------------------------


def h_majorana_source(sp):
    # gamma_s = c_s + c^dag_s: Hermitian and parity odd
    return sp.majorana(0, 0)


def h_nn_plus_majorana(sp):
    # parity mixed: an even and an odd part
    return h_spinless_nn(sp) + 0.3 * sp.majorana(0, 0)


def h_hop_plus(epsilon):
    def builder(sp):
        return h_spinless_hop(sp) + epsilon * sp.majorana(0, 0)

    return builder


def parity_input(builder, bond, user_observable=True, fermion=True):
    param = fermion_model_input(SPINLESS, builder, [4, 4], [bond])
    if user_observable:
        # tenes_std then does not add the Hamiltonian as an observable,
        # which is the input that used to run through silently
        with_user_twosite_observable(param)
    if not fermion:
        del param["parameter"]["general"]["fermion"]
        for site in param["tensor"]["unitcell"]:
            del site["parity"]
    return param


class _ModelPath:
    """Path of one bond of a Model, shaped for infer_parity_chain."""

    def __init__(self, model, bond):
        self.unitcell = model.unitcell
        self.path = model.graph.make_path(tenes_std.Bond(*bond))
        self.sites = [int(b.source_site) for b in self.path]
        self.sites.append(int(self.unitcell.target_site(self.path[-1])))

    def phys_parity(self, position):
        return list(self.unitcell.sites[self.sites[position]].parity)


class TestFermionHamiltonianTermParity:
    @pytest.mark.parametrize("user_observable", [True, False], ids=["user", "auto"])
    @pytest.mark.parametrize(
        "bond",
        [(5, 2, 0), (5, 2, 1), (5, 1, 0)],
        ids=["two-hop", "three-hop", "nearest-neighbour"],
    )
    @pytest.mark.parametrize(
        "builder",
        [h_majorana_source, h_nn_plus_majorana],
        ids=["odd", "mixed"],
    )
    def test_odd_or_mixed_bond_term_is_rejected(self, builder, bond, user_observable):
        param = parity_input(builder, bond, user_observable)
        with pytest.raises(RuntimeError) as excinfo:
            tenes_std.Model(param)
        msg = str(excinfo.value)
        assert re.search("(?i)parity", msg), msg
        assert re.search(displacement_pattern(bond[1], bond[2]), msg), msg

    @pytest.mark.parametrize("bond", [(5, 2, 0), (5, 2, 1), (5, 1, 0)])
    def test_even_bond_term_is_accepted(self, bond):
        model = tenes_std.Model(parity_input(h_hop_plus(0.0), bond))
        path = model.graph.make_path(tenes_std.Bond(*bond))
        assert len([g for g in model.simple_updates if g.group == 0]) == len(path)

    @pytest.mark.parametrize("bond", [(5, 2, 0), (5, 1, 0)])
    def test_odd_bond_term_is_accepted_without_fermion_mode(self, bond):
        model = tenes_std.Model(
            parity_input(h_majorana_source, bond, user_observable=True, fermion=False)
        )
        path = model.graph.make_path(tenes_std.Bond(*bond))
        assert len([g for g in model.simple_updates if g.group == 0]) == len(path)

    @pytest.mark.parametrize("bond", [(5, 2, 1), (5, 1, 0)])
    def test_odd_residue_below_atol_is_accepted(self, bond):
        # 5e-16 < atol = 1e-15: rounding noise, not an odd term.  The gates
        # that come out must still be exactly even (the solver counts every
        # nonzero element).
        model = tenes_std.Model(parity_input(h_hop_plus(5e-16), bond))
        gates = [g.elements for g in model.simple_updates if g.group == 0]
        # for one hop infer_parity_chain checks the single gate against the
        # physical ledgers of both sites
        _, problems = infer_parity_chain(_ModelPath(model, bond), gates)
        assert not problems, problems

    @pytest.mark.parametrize("bond", [(5, 2, 1), (5, 1, 0)])
    def test_odd_part_above_atol_is_rejected(self, bond):
        with pytest.raises(RuntimeError, match="(?i)parity"):
            tenes_std.Model(parity_input(h_hop_plus(1e-9), bond))


# ---------------------------------------------------------------------------
# PR #118 review, item 2: the fermion chain must not build the dense
# operator over the whole path (d^(2 nsites) elements: 4^14 for a Hubbard
# bond of six hops).  The gates themselves need chi <= r * d_m, r being the
# operator Schmidt rank of exp(-tau H) across (source | target).
# ---------------------------------------------------------------------------

LONG_CHAIN_PEAK_BYTES = 64 * 1024 * 1024
LONG_CHAIN_SECONDS = 10.0


def operator_schmidt_rank(evo, tol=1e-12):
    ds, dt = evo.shape[0], evo.shape[1]
    M = evo.transpose(0, 2, 1, 3).reshape(ds * ds, dt * dt)
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s > tol * s[0]))


def test_long_hubbard_chain_stays_small():
    """Hubbard (d = 4) bonds of five hops ((3, 2)) and six hops ((3, 3), 4x4
    window edge) from site 5 of a 4x4 cell: peak numpy allocation (traced
    by tracemalloc) and wall time stay small, every chi is at most r * d_m,
    and every chi index has one parity (contract item 2).

    The five-hop bond comes first on purpose: building the dense path
    operator for it already needs 4^12 * 8 B = 128 MiB, so an implementation
    that still does that fails here before the six-hop bond would ask for
    2 GiB.
    """
    for disp in [(3, 2), (3, 3)]:
        c = Case("4x4", disp, "hubbard_hopV")
        assert len(c.path) == abs(disp[0]) + abs(disp[1]), c.describe()
        tracemalloc.start()
        t0 = time.perf_counter()
        try:
            gates = [g.elements for g in c.fermion_gates()]
            elapsed = time.perf_counter() - t0
            _, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        info = "{}: peak {:.1f} MiB, {:.2f} s".format(
            c.describe(), peak / 2**20, elapsed
        )
        assert peak <= LONG_CHAIN_PEAK_BYTES, info
        assert elapsed <= LONG_CHAIN_SECONDS, info

        evo = path_reference(c.kind, c.builder, c.tau, 2)
        r = operator_schmidt_rank(evo)
        assert len(gates) == len(c.path), info
        for k, G in enumerate(gates[:-1]):
            d_m = len(c.phys_parity(k + 1))
            assert G.shape[3] <= r * d_m, "{}: gate {} chi {} > r d = {} x {}".format(
                info, k, G.shape[3], r, d_m
            )
        _, problems = infer_parity_chain(c, gates)
        assert not problems, "{}:\n  {}".format(info, "\n  ".join(problems))


def test_six_hop_spinless_chain_composes_to_the_path_ordered_operator():
    """The six-hop (3, 3) bond with d = 2: the dense reference is 2^7 x 2^7,
    cheap enough for the exact composition check of contract item 1."""
    c = Case("4x4", (3, 3), "spinless_hop")
    assert len(c.path) == 6, c.describe()
    gates = [g.elements for g in c.fermion_gates()]
    got = compose_chain(gates)
    ref = path_reference(c.kind, c.builder, c.tau, c.nsites)
    err = np.max(np.abs(got - ref))
    assert err <= ATOL, "{}: {:.3e}".format(c.describe(), err)
    _, problems = infer_parity_chain(c, gates)
    assert not problems, problems
