#!/usr/bin/env python3
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
# along with this program. If not, see http://www.gnu.org/licenses/.

"""Constants of test/fermion/longrange_gate.cpp (task T2).

Behaviour contract: "振る舞い契約書(T2)" in
docs/superpowers/plans/2026-09-30-fermion-longrange-hamiltonian.md; design in
docs/superpowers/specs/2026-09-30-fermion-longrange-hamiltonian-design.md.

Run from anywhere:

    python3 test/fermion/gen_longrange_gate.py --write

rewrites the block between the "BEGIN GENERATED" / "END GENERATED" markers of
longrange_gate.cpp; without --write it prints the block to stdout.

What is generated, and where each number comes from
---------------------------------------------------
* The gate chains: ``tenes_std.make_evolution_twosite(..., fermion=True)``
  (tool/tenes_std.py, task T1) applied to a two-site Hamiltonian written here.
  The chains are the INPUT of the solver under test; nothing in the reference
  below is derived from them.
* The expected parity ledger of every gate leg (contract items 3 and 4): read
  off the nonzero elements of the chain with the rule of design section 6.1,
  p(out2) = p(in1) + p(in2) + p(out1), walking the chain from the physical
  tables of the unit cell.  Mixed or empty out2 indices abort the generator.
* The reference overlaps (contract items 1 and 2).  The unit-cell tensors
  before the update are the deterministic formula of fock_oracle.py
  (``det`` below; longrange_gate.cpp evaluates the same formula), nonzero only
  on two labels (one even, one odd) of every virtual leg.  An open patch of
  them is turned into a Fock-space state by ``SparseOracle``, a sparse
  re-implementation of fock_oracle.Oracle that also handles several physical
  modes per site (Hubbard) and a relabelled spinless site; it is checked
  against fock_oracle.Oracle itself on d = 2 patches (``self_check``).  The
  long-range operator is applied as a Fock-space operator,

      O_st = sum evo[i_s, i_t, o_s, o_t] M_{o_s}(s) M_{o_t}(t) P_0(s, t)
                                         M_{i_t}(t)^dag M_{i_s}(s)^dag,

  with M_i the ordered creation monomial of local state i and P_0 the
  projector on the empty s and t.  That is design section 3 read literally
  (source-first ordered two-site basis); the Jordan-Wigner string of every
  other mode comes out of the Fock algebra by itself.  evo = expm(-tau H) is
  computed here with scipy, not taken from tenes_std.  No sign formula of the
  design (section 4.1) or of tenes_std is used for the reference.
* The overlaps <phi_k|psi> with a few fixed dense bra states phi_k (formula
  ``bra`` below, also evaluated in longrange_gate.cpp), for psi before the
  update and for O_st psi.  The test compares them up to one overall scalar.
* For contract item 2 (a unit-cell site that comes back on the path, 2x2
  cell) only the chains are generated; the test compares the 2x2 run with
  the same gates on a supercell (see longrange_gate.cpp, T2-2).

Before writing, the generator checks its own reference: SparseOracle
against fock_oracle.Oracle, and for every case that O_st psi differs from
psi and (for operators with an odd channel) from the state the bosonic
decomposition would give (the Jordan-Wigner string along the path dropped),
both by more than 1e-4, so that no reference is blind to what it is meant
to catch.
"""

import argparse
import cmath
import math
import os
import sys
from itertools import product

import numpy as np
import scipy.linalg

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tool"))
sys.path.insert(0, HERE)

import fock_oracle  # noqa: E402
import tenes_std  # noqa: E402

LEGS = fock_oracle.LEGS  # ("l", "t", "r", "b")
OP_ORDER = fock_oracle.OP_ORDER  # ("b", "r", "t", "l")
NBRA = 6


def popcount(x):
    # int.bit_count() needs Python 3.10; CI runs 3.9.
    return bin(x).count("1")


# ---------------------------------------------------------------------------
# Site kinds
# ---------------------------------------------------------------------------


class Kind:
    """A fermionic site: nspin modes; local state i is the occupation bit
    pattern occ[i] (bit k = mode k), created as c^dag_{k0} c^dag_{k1} ... |0>
    with k0 < k1 < ...  (tool/tenes_simple.py HubbardModel: i = n_up + 2 n_dn,
    |up dn> = c^dag_up c^dag_dn |0>)."""

    def __init__(self, code, name, nspin, occ):
        self.code = code
        self.name = name
        self.nspin = nspin
        self.occ = list(occ)
        self.d = len(self.occ)
        self.parity = [popcount(b) % 2 for b in self.occ]


SPINLESS = Kind(0, "S", 1, [0, 1])
# spinless with the local labels swapped: index 0 is the occupied state
FLIPPED = Kind(1, "F", 1, [1, 0])
HUBBARD = Kind(2, "H", 2, [0, 1, 2, 3])
KINDS = [SPINLESS, FLIPPED, HUBBARD]


# ---------------------------------------------------------------------------
# Two-site Hamiltonians (dense Fock space of the two sites, source first)
# ---------------------------------------------------------------------------


def creation(mode, nmodes):
    dim = 1 << nmodes
    mat = np.zeros((dim, dim))
    below = (1 << mode) - 1
    for g in range(dim):
        if (g >> mode) & 1:
            continue
        sign = -1.0 if popcount(g & below) % 2 else 1.0
        mat[g | (1 << mode), g] = sign
    return mat


def two_site_operator(ks, kt, builder):
    """H[i_s, i_t, o_s, o_t] = <o_s o_t| H |i_s i_t> in the source-first
    ordered basis |i_s i_t> = M_{i_s}(s) M_{i_t}(t) |0> (design section 3)."""
    nm = ks.nspin + kt.nspin
    cd = [
        [creation(k, nm) for k in range(ks.nspin)],
        [creation(ks.nspin + k, nm) for k in range(kt.nspin)],
    ]
    c = [[m.T.copy() for m in row] for row in cd]
    mat = builder(cd, c)
    gvec = np.array(
        [ks.occ[i] | (kt.occ[j] << ks.nspin) for i in range(ks.d) for j in range(kt.d)]
    )
    T = mat[np.ix_(gvec, gvec)].reshape([ks.d, kt.d, ks.d, kt.d])  # [o, i]
    return T.transpose(2, 3, 0, 1)


def _n(cd, c, k):
    return sum(cd[k][s] @ c[k][s] for s in range(len(cd[k])))


def h_hop(cd, c):
    """-t (c^dag_s c_t + h.c.) on mode 0 of each site (t = 1)."""
    return -1.0 * (cd[0][0] @ c[1][0] + cd[1][0] @ c[0][0])


def h_nn(cd, c):
    """V n_s n_t (V = 0.8): even channel only."""
    return 0.8 * _n(cd, c, 0) @ _n(cd, c, 1)


def h_hop_nn(cd, c):
    return h_hop(cd, c) + h_nn(cd, c)


def h_hubbard(cd, c):
    """Hubbard hopping of both spins plus the density interaction V n_s n_t
    (the operator of the T1 tests)."""
    h = 0.8 * _n(cd, c, 0) @ _n(cd, c, 1)
    for s in range(2):
        h = h - (cd[0][s] @ c[1][s] + cd[1][s] @ c[0][s])
    return h


def h_majorana(cd, c):
    """H = i gamma_s gamma_t (gamma = c + c^dag of mode 0): H^2 = 1, so
    exp(-tau H) at tau = i pi / 2 is exactly gamma_s gamma_t, a real product
    operator in the ODD channel.  Its chain has chi = d on the intermediate
    site with the physical parity table flipped (T1 contract item 3)."""
    gs = cd[0][0] + c[0][0]
    gt = cd[1][0] + c[1][0]
    return 1j * gs @ gt


def h_product(cd, c):
    """A sum of one-site terms: exp(-tau H) is a product operator, even
    channel only, chi = d with the physical parity table."""
    return -0.7 * _n(cd, c, 0) + 0.4 * _n(cd, c, 1)


BUILDERS = {
    "hop": h_hop,
    "nn": h_nn,
    "hopnn": h_hop_nn,
    "hubbard": h_hubbard,
    "majorana": h_majorana,
    "product": h_product,
}


def evolution(ks, kt, builder, tau):
    """evo = expm(-tau H) as [i_s, i_t, o_s, o_t] (scipy, independent of
    tenes_std's eigh route)."""
    H = two_site_operator(ks, kt, builder)
    d2 = ks.d * kt.d
    # matrix convention mat[out, in]
    mat = H.transpose(2, 3, 0, 1).reshape(d2, d2)
    E = scipy.linalg.expm(-tau * mat)
    return H, E.reshape(ks.d, kt.d, ks.d, kt.d).transpose(2, 3, 0, 1)


# ---------------------------------------------------------------------------
# Sparse Fock-space oracle (fock_oracle.Oracle, generalised)
# ---------------------------------------------------------------------------


def _create(vec, mode):
    out = {}
    bit = 1 << mode
    below = bit - 1
    for state, amp in vec.items():
        if state & bit:
            continue
        sign = -1.0 if popcount(state & below) % 2 else 1.0
        key = state | bit
        out[key] = out.get(key, 0.0) + sign * amp
    return out


def _annihilate(vec, mode):
    out = {}
    bit = 1 << mode
    below = bit - 1
    for state, amp in vec.items():
        if not state & bit:
            continue
        sign = -1.0 if popcount(state & below) % 2 else 1.0
        key = state ^ bit
        out[key] = out.get(key, 0.0) + sign * amp
    return out


def _axpy(out, vec, coeff):
    for state, amp in vec.items():
        out[state] = out.get(state, 0.0) + coeff * amp


class SparseOracle:
    """fock_oracle.Oracle on sparse states, with a list of internal bonds
    given explicitly and several physical modes per site.

    Mode layout: the physical modes of site 0, of site 1, ... (a Hubbard site
    owns two consecutive modes, up then down), then two auxiliary modes per
    internal bond in bond order, exactly as fock_oracle.Oracle allocates them.
    Every step (bond creators in bond order, then the site projectors in site
    order, virtual legs annihilated in the order b, r, t, l) is that of
    fock_oracle.Oracle.
    """

    def __init__(self, kinds, bonds, tensors, leg_parities):
        self.kinds = kinds
        self.bonds = bonds
        self.tensors = tensors
        self.leg_parities = leg_parities
        self.base = []
        nm = 0
        for kind in kinds:
            self.base.append(nm)
            nm += kind.nspin
        self.nphys = nm
        self.mode = {}
        for a, aleg, b, bleg in bonds:
            self.mode[(a, aleg)] = nm
            nm += 1
            self.mode[(b, bleg)] = nm
            nm += 1
        self.nmode = nm

    def create_local(self, vec, site, i):
        """M_i(site) vec, M_i = c^dag_{k0} c^dag_{k1} ... (k0 < k1)."""
        kind = self.kinds[site]
        ks = [k for k in range(kind.nspin) if (kind.occ[i] >> k) & 1]
        for k in reversed(ks):
            vec = _create(vec, self.base[site] + k)
        return vec

    def annihilate_local(self, vec, site, i):
        """M_i(site)^dag vec = ... c_{k1} c_{k0} vec."""
        kind = self.kinds[site]
        ks = [k for k in range(kind.nspin) if (kind.occ[i] >> k) & 1]
        for k in ks:
            vec = _annihilate(vec, self.base[site] + k)
        return vec

    def state(self):
        vec = {0: 1.0}
        for a, aleg, b, bleg in self.bonds:
            pa = self.leg_parities[a][LEGS.index(aleg)]
            pb = self.leg_parities[b][LEGS.index(bleg)]
            out = {}
            for label in range(min(len(pa), len(pb))):
                if pa[label] != pb[label]:
                    continue
                term = vec
                if pa[label]:
                    term = _create(term, self.mode[(a, aleg)])
                    term = _create(term, self.mode[(b, bleg)])
                _axpy(out, term, 1.0)
            vec = out
        for site in range(len(self.kinds)):
            vec = self.project(vec, site)
        return vec

    def project(self, vec, site):
        parity = self.leg_parities[site]
        tensor = self.tensors[site]
        out = {}
        for idx in product(*[range(len(p)) for p in parity]):
            coeff = tensor[idx]
            if coeff == 0.0:
                continue
            term = vec
            for leg in OP_ORDER:
                label = idx[LEGS.index(leg)]
                if not parity[LEGS.index(leg)][label]:
                    continue
                key = (site, leg)
                if key not in self.mode:
                    term = {}
                    break
                term = _annihilate(term, self.mode[key])
            if not term:
                continue
            term = self.create_local(term, site, idx[4])
            _axpy(out, term, coeff)
        return out

    def physical_state(self):
        vec = self.state()
        mask = (1 << self.nphys) - 1
        return {s: a for s, a in vec.items() if (s & ~mask) == 0 and a != 0.0}

    def apply_twosite(self, vec, s, t, evo, string_sites=(), string_sign=False):
        """O_st vec with O_st = sum evo M_o(s) M_o(t) P_0 M_i(t)^dag M_i(s)^dag.

        With string_sign = True, each term is multiplied by
        (-1)^{p_c * N(string_sites)} (p_c the channel parity): the operator a
        chain composes to when the Jordan-Wigner string along the path is
        dropped (the bosonic decomposition).  Only used to show that the
        reference is not blind to that error."""
        ks = self.kinds[s]
        kt = self.kinds[t]
        stmask = 0
        for site in (s, t):
            for k in range(self.kinds[site].nspin):
                stmask |= 1 << (self.base[site] + k)
        strmask = 0
        for site in string_sites:
            for k in range(self.kinds[site].nspin):
                strmask |= 1 << (self.base[site] + k)
        out = {}
        for i_s in range(ks.d):
            w1 = self.annihilate_local(vec, s, i_s)
            if not w1:
                continue
            for i_t in range(kt.d):
                w = self.annihilate_local(w1, t, i_t)
                w = {k: v for k, v in w.items() if (k & stmask) == 0}
                if not w:
                    continue
                for o_s in range(ks.d):
                    for o_t in range(kt.d):
                        coeff = evo[i_s, i_t, o_s, o_t]
                        if coeff == 0.0:
                            continue
                        term = self.create_local(w, t, o_t)
                        term = self.create_local(term, s, o_s)
                        if string_sign and (kt.parity[i_t] + kt.parity[o_t]) % 2:
                            term = {
                                k: (-v if popcount(k & strmask) % 2 else v)
                                for k, v in term.items()
                            }
                        _axpy(out, term, coeff)
        return out

    def bits(self, local):
        g = 0
        for site, i in enumerate(local):
            g |= self.kinds[site].occ[i] << self.base[site]
        return g


def self_check():
    """SparseOracle against fock_oracle.Oracle on d = 2 patches."""
    for lx, ly, seed in [(2, 1, 0), (1, 2, 3), (2, 2, 5), (3, 1, 7)]:
        patch, tensors, leg_parities = fock_oracle.make_case(
            lx, ly, [False, True], seed
        )
        ref = fock_oracle.Oracle(patch, tensors, leg_parities).physical_state()
        sp = SparseOracle(
            [SPINLESS] * patch.nsite, patch.internal_bonds(), tensors, leg_parities
        )
        got = sp.physical_state()
        # both put the physical modes first, one per site
        dense = np.zeros(1 << patch.nsite)
        for s, a in got.items():
            dense[s] = a
        want = np.zeros(1 << patch.nsite)
        for s in range(len(ref)):
            if ref[s] != 0.0:
                assert s < (1 << patch.nsite), "fock_oracle state has aux bits"
                want[s] = ref[s]
        err = np.max(np.abs(dense - want))
        assert err < 1e-14, "SparseOracle differs from fock_oracle: {}".format(err)


# ---------------------------------------------------------------------------
# Deterministic tensors and bras (mirrored in longrange_gate.cpp)
# ---------------------------------------------------------------------------


def det(site, seed, idx):
    """fock_oracle.deterministic_tensor's formula for one element."""
    x = (site + 2) * (1 + seed + sum((ax + 3 + seed % 5) * idx[ax] for ax in range(5)))
    return 0.19 * math.sin(x) + 0.13 * math.cos(0.37 * x)


def bra(k, local):
    """phi_k(i_0, ..., i_{N-1}) for the local indices of the patch sites in
    raster order."""
    x = 0.4 + 1.1 * k
    for p, i in enumerate(local):
        x += (0.23 + 0.07 * p + 0.03 * k) * (i + 1) * (p + 1)
    return math.cos(x)


def even_first(dim):
    neven = (dim + 1) // 2
    return [0] * neven + [1] * (dim - neven)


def support_labels(dim):
    """The two labels a leg of this dimension carries nonzeros on: the first
    even and the first odd index of its even-first ledger."""
    if dim == 1:
        return [0]
    return [0, (dim + 1) // 2]


# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------


class Cell:
    def __init__(self, lx, ly, kinds):
        self.lx = lx
        self.ly = ly
        self.kinds = kinds  # per unit-cell site

    def site(self, X, Y):
        return (X % self.lx) + self.lx * (Y % self.ly)

    def unitcell(self):
        return tenes_std.Unitcell(
            {
                "l_sub": [self.lx, self.ly],
                "unitcell": [
                    {
                        "index": [i],
                        "physical_dim": k.d,
                        "virtual_dim": 2,
                        "parity": k.parity,
                    }
                    for i, k in enumerate(self.kinds)
                ],
            }
        )


def make_graph(uc, bonds):
    """The graph tenes_std.Model would build for these Hamiltonian bonds."""
    ox_min = oy_min = ox_max = oy_max = 0
    for b in bonds:
        ox, oy = uc.target_offset(b)
        ox_min, ox_max = min(ox, ox_min), max(ox, ox_max)
        oy_min, oy_max = min(oy, oy_min), max(oy, oy_max)
    return tenes_std.LatticeGraph(uc, ox_min, oy_min, ox_max, oy_max)


DIRS = {0: (-1, 0), 1: (0, 1), 2: (1, 0), 3: (0, -1)}


def chain_for(cell, source, disp, opname, tau):
    """Gates of tenes_std for the bond (source, dx, dy), their positions and
    the expected ledgers."""
    uc = cell.unitcell()
    bond = tenes_std.Bond(source, disp[0], disp[1])
    graph = make_graph(uc, [bond])
    target = uc.target_site(bond)
    ks = cell.kinds[source]
    kt = cell.kinds[target]
    H2 = two_site_operator(ks, kt, BUILDERS[opname])
    gates = tenes_std.make_evolution_twosite(
        tenes_std.NNOperator(bond, elements=H2), graph, tau, fermion=True
    )
    path = graph.make_path(bond)
    assert len(gates) == len(path)
    # global coordinates of the path positions, source at (0, 0)
    pos = [(0, 0)]
    sites = [source]
    legs = []
    X0, Y0 = source % cell.lx, source // cell.lx
    for g in gates:
        leg = uc.bond_direction(g.bond)
        legs.append(leg)
        assert int(g.bond.source_site) == sites[-1]
        dx, dy = DIRS[leg]
        pos.append((pos[-1][0] + dx, pos[-1][1] + dy))
        sites.append(cell.site(X0 + pos[-1][0], Y0 + pos[-1][1]))
    assert pos[-1] == tuple(disp), (pos, disp)
    assert sites[-1] == target
    # expected ledgers, walking the chain (design section 6.1)
    current = {s: list(cell.kinds[s].parity) for s in range(len(cell.kinds))}
    ledgers = []
    for k, g in enumerate(gates):
        a, b = sites[k], sites[k + 1]
        G = np.asarray(g.elements)
        p_in1 = current[a]
        p_in2 = current[b]
        p_out1 = list(cell.kinds[a].parity)
        assert G.shape[0] == len(p_in1) and G.shape[1] == len(p_in2)
        assert G.shape[2] == len(p_out1)
        seen = [set() for _ in range(G.shape[3])]
        for i1, i2, o1, o2 in np.argwhere(G != 0):
            seen[o2].add((p_in1[i1] + p_in2[i2] + p_out1[o1]) % 2)
        assert all(len(s) == 1 for s in seen), "out2 not uniquely graded"
        p_out2 = [s.pop() for s in seen]
        ledgers.append((list(p_in1), list(p_in2), p_out1, p_out2))
        current[a] = p_out1
        current[b] = p_out2
    for s in range(len(cell.kinds)):
        assert current[s] == list(cell.kinds[s].parity), "chain does not close"
    _, evo = evolution(ks, kt, BUILDERS[opname], tau)
    return dict(
        gates=gates,
        legs=legs,
        pos=pos,
        sites=sites,
        ledgers=ledgers,
        evo=evo,
        target=target,
    )


def vdims_for(cell, pos, sites, dc):
    """Virtual dimensions: dc on the unit-cell bonds of the path (one value,
    or one per hop), 2 elsewhere."""
    n = len(cell.kinds)
    vd = [[2, 2, 2, 2] for _ in range(n)]
    if isinstance(dc, int):
        dc = [dc] * (len(pos) - 1)
    for k in range(len(pos) - 1):
        (x0, y0), (x1, y1) = pos[k], pos[k + 1]
        leg = {v: key for key, v in DIRS.items()}[(x1 - x0, y1 - y0)]
        a, b = sites[k], sites[k + 1]
        vd[a][leg] = dc[k]
        vd[b][(leg + 2) % 4] = dc[k]
    return vd


def cell_tensor_value(site, seed, small, big_parities, is_complex):
    """Element of the initial unit-cell tensor at the support index `small`
    (0/1 per virtual leg, physical index), or 0 when it is parity odd."""
    odd = sum(big_parities[ax][small[ax]] for ax in range(5)) % 2
    if odd:
        return 0.0
    re = det(site, seed, small)
    if not is_complex:
        return re
    return complex(re, det(site, seed + 1000, small))


def reference(cell, chain, vd, seed, is_complex):
    """Patch, oracle state before, O_st applied, and the bra overlaps."""
    xs = [p[0] for p in chain["pos"]]
    ys = [p[1] for p in chain["pos"]]
    xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)
    ncol = xmax - xmin + 1
    nrow = ymax - ymin + 1
    X0 = chain["sites"][0] % cell.lx
    Y0 = chain["sites"][0] // cell.lx
    patch = fock_oracle.Patch(ncol, nrow)
    psites = []  # unit-cell site at each patch position (raster)
    for row in range(nrow):
        for col in range(ncol):
            X = X0 + xmin + col
            Y = Y0 + ymax - row
            psites.append(cell.site(X, Y))
    assert len(set(psites)) == len(psites), "a unit-cell site repeats in the patch"
    # patch position of each path position
    ppos = [(ymax - y) * ncol + (x - xmin) for x, y in chain["pos"]]
    kinds = [cell.kinds[s] for s in psites]
    internal = set()
    for a, aleg, b, bleg in patch.internal_bonds():
        internal.add((a, aleg))
        internal.add((b, bleg))
    tensors = []
    leg_parities = []
    for p, s in enumerate(psites):
        kind = cell.kinds[s]
        big = []
        for leg in range(4):
            dim = vd[s][leg]
            lab = support_labels(dim)
            big.append([even_first(dim)[l] for l in lab])
        big.append(kind.parity)
        lp = []
        for leg in range(4):
            if (p, LEGS[leg]) in internal:
                lp.append([bool(x) for x in big[leg]])
            else:
                lp.append([False])
        lp.append([bool(x) for x in kind.parity])
        shape = tuple(len(x) for x in lp)
        T = np.zeros(shape, dtype=complex if is_complex else float)
        for idx in product(*[range(n) for n in shape]):
            T[idx] = cell_tensor_value(s, seed, idx, big, is_complex)
        tensors.append(T)
        leg_parities.append(lp)
    oracle = SparseOracle(kinds, patch.internal_bonds(), tensors, leg_parities)
    psi = oracle.physical_state()
    s_pos, t_pos = ppos[0], ppos[-1]
    after = oracle.apply_twosite(psi, s_pos, t_pos, chain["evo"])
    nostring = oracle.apply_twosite(
        psi, s_pos, t_pos, chain["evo"], string_sites=ppos[1:-1], string_sign=True
    )
    dims = [k.d for k in kinds]
    locals_ = list(product(*[range(d) for d in dims]))

    def overlaps(vec):
        out = []
        for k in range(NBRA):
            acc = 0.0
            for loc in locals_:
                amp = vec.get(oracle.bits(loc), 0.0)
                if amp != 0.0:
                    acc += bra(k, loc) * amp
            out.append(acc)
        return np.array(out, dtype=complex)

    ob = overlaps(psi)
    oa = overlaps(after)
    on = overlaps(nostring)
    return dict(
        nrow=nrow,
        ncol=ncol,
        psites=psites,
        ppos=ppos,
        before=ob,
        after=oa,
        nostring=on,
    )


def residual(got, want):
    """|got - c want| / |got| with the best complex scalar c."""
    c = np.vdot(want, got) / np.vdot(want, want)
    return float(np.linalg.norm(got - c * want) / np.linalg.norm(got))


# ---------------------------------------------------------------------------
# Case list
# ---------------------------------------------------------------------------

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

UNIFORM_S = Cell(4, 4, [SPINLESS] * 16)
UNIFORM_H = Cell(4, 4, [HUBBARD] * 16)


def _checker_sf(x, y):
    return SPINLESS if (x + y) % 2 == 0 else FLIPPED


def _three(x, y):
    return [SPINLESS, HUBBARD, FLIPPED][(x + y) % 3]


def cell_from(fn, lx=4, ly=4):
    return Cell(lx, ly, [fn(i % lx, i // lx) for i in range(lx * ly)])


# S at (1,1), (0,1) ...; F on the other colour: the intermediate site of
# every two-hop path from site 5 is F while source and target are S.
MIXED_SF = cell_from(_checker_sf)
# kind by (x + y) mod 3; source 5 = (1, 1) is F.
MIXED_SHF = cell_from(_three)
# Hubbard ends, spinless middle: site 5 = (1, 1) and 11 = (3, 2) Hubbard,
# the rest spinless.  The (-2, 1) path 5 -> 9 -> 8 -> 11 is H-S-S-H: a
# three-hop chain of the Hubbard operator whose middle gate stays small
# (chi = 32 on the spinless sites instead of 64 x 64 elements).
HSSH = Cell(4, 4, [HUBBARD if i in (5, 11) else SPINLESS for i in range(16)])

TAU = 0.1
TAU_C = complex(0.1, 0.05)
MAJ_TAU = 0.5j * math.pi

# (name, cell, source, disp, operator, tau, is_complex, dc, seed)
SU_CASES = []
for dd in DISPLACEMENTS:
    for op in ("hop", "nn"):
        SU_CASES.append(
            ("d2 {} ({},{})".format(op, *dd), UNIFORM_S, 5, dd, op, TAU, False, 16, 3)
        )
SU_CASES += [
    ("d2 complex hopnn (-2,1)", UNIFORM_S, 5, (-2, 1), "hopnn", TAU_C, True, 16, 7),
    ("d2 complex hopnn (0,-2)", UNIFORM_S, 5, (0, -2), "hopnn", TAU_C, True, 16, 7),
    ("d2 complex hopnn (-1,1)", UNIFORM_S, 5, (-1, 1), "hopnn", TAU_C, True, 16, 9),
    ("mixed S/F hopnn (-2,0)", MIXED_SF, 5, (-2, 0), "hopnn", TAU, False, 16, 11),
    ("mixed S/F hopnn (1,1)", MIXED_SF, 5, (1, 1), "hopnn", TAU, False, 16, 11),
    ("mixed S/H/F hopnn (-2,1)", MIXED_SHF, 5, (-2, 1), "hopnn", TAU, False, 32, 13),
    ("mixed S/H/F hopnn (0,2)", MIXED_SHF, 5, (0, 2), "hopnn", TAU, False, 32, 13),
    # chi = d: the intermediate site 9 is F ([1, 0]); the odd product
    # operator gamma_s gamma_t leaves an even-first chi ledger [0, 1] there.
    (
        "chi=d majorana S/F (-1,1)",
        MIXED_SF,
        5,
        (-1, 1),
        "majorana",
        MAJ_TAU,
        False,
        16,
        5,
    ),
    # hops up then left (source_leg 1, then 0) through a Hubbard site whose
    # physical leg becomes chi = 64
    ("d4 hubbard (-1,1)", UNIFORM_H, 5, (-1, 1), "hubbard", TAU, False, 32, 17),
    # three hops 5 -> 9 -> 8 -> 11 (source_leg 1, 0, 0), Hubbard at both ends;
    # an end bond can reach rank 8 x 4 = 32 (the other legs of its d = 4 end
    # times d), the middle bond more (2 x 16 through the chi legs, and the
    # driver's bound is larger), so it gets 64
    (
        "d4 hubbard H-S-S-H (-2,1)",
        HSSH,
        5,
        (-2, 1),
        "hubbard",
        TAU,
        False,
        [32, 64, 32],
        19,
    ),
]


# ---------------------------------------------------------------------------
# Contract item 2: a unit-cell site that comes back on the path (2x2 cell)
#
# Only the chains are generated here (see longrange_gate.cpp, T2-2, for why
# no Fock reference can be given and what the test compares instead).
# (2, 0): 0 -> 1 -> 0, the target is the source site;
# (3, 0): 0 -> 1 -> 0 -> 1, site 0 is the source and the second intermediate
#         site, site 1 the first intermediate site and the target;
# (2, 1): 0 -> 1 -> 0 -> 2, site 0 is the source and the second intermediate
#         site, the last hop goes up.
# ---------------------------------------------------------------------------

CELL_2X2 = Cell(2, 2, [SPINLESS] * 4)
CELL_2X2_MIXED = Cell(2, 2, [SPINLESS, FLIPPED, SPINLESS, FLIPPED])

# name, cell, displacement, virtual dimension of the path bonds, seed.  The
# comparison of T2-2 does not need an untruncated update (the 2x2 run and the
# supercell run truncate alike), so the bonds stay small.
REVISIT_CASES = [
    ("revisit (2,0) in 2x2: target = source site", CELL_2X2, (2, 0), 8, 23),
    ("revisit (3,0) in 2x2: source site again in the middle", CELL_2X2, (3, 0), 8, 29),
    (
        "revisit (2,1) in 2x2 S/F: source site again in the middle, then up",
        CELL_2X2_MIXED,
        (2, 1),
        8,
        31,
    ),
]


def fmt(v):
    return repr(float(v))


def cpp_bool_string(p):
    return '"' + "".join("1" if x else "0" for x in p) + '"'


def emit_entries(arrays, array_name, G):
    """A static array of the nonzero elements of gate G (every element the
    tool writes, down to the round-off of its SVDs, as the solver would read
    them from the TOML file)."""
    items = []
    for idx in np.argwhere(G != 0):
        v = G[tuple(idx)]
        ij = ",".join(str(int(i)) for i in idx)
        if np.imag(v) == 0.0:
            items.append("{{{{{}}},{}}}".format(ij, fmt(np.real(v))))
        else:
            items.append(
                "{{{{{}}},{},{}}}".format(ij, fmt(np.real(v)), fmt(np.imag(v)))
            )
    arrays.append("const lg_entry {}[] = {{".format(array_name))
    for i in range(0, len(items), 3):
        arrays.append("    " + ",".join(items[i : i + 3]) + ",")
    arrays.append("};")


def emit_chain(arrays, out, tag, chain, indent):
    """The gates of a chain as an initializer of std::vector<lg_gate_data>."""
    out.append(indent + "{")
    for k, (g, leg, led) in enumerate(
        zip(chain["gates"], chain["legs"], chain["ledgers"])
    ):
        G = np.asarray(g.elements)
        array_name = "lg_e_{}_{}".format(tag, k)
        emit_entries(arrays, array_name, G)
        out.append(
            indent
            + "    {{{}, {}, {{{}}}, {{{}}}, {}, sizeof({}) / sizeof(lg_entry)}},".format(
                int(g.bond.source_site),
                leg,
                ", ".join(str(x) for x in G.shape),
                ", ".join(cpp_bool_string(p) for p in led),
                array_name,
                array_name,
            )
        )
    out.append(indent + "},")


def emit_case(
    arrays, out, tag, name, cell, source, disp, op, tau, is_complex, dc, seed
):
    chain = chain_for(cell, source, disp, op, tau)
    vd = vdims_for(cell, chain["pos"], chain["sites"], dc)
    ref = reference(cell, chain, vd, seed, is_complex)
    r_before = residual(ref["after"], ref["before"])
    r_nostring = residual(ref["after"], ref["nostring"])
    out.append(
        "  // {}: path {} through unit-cell sites {}".format(
            name, " ".join("leg{}".format(l) for l in chain["legs"]), chain["sites"]
        )
    )
    out.append(
        "  //   residual(after vs before) = {:.3e}, "
        "residual(after vs string dropped) = {:.3e}".format(r_before, r_nostring)
    )
    out.append("  {")
    out.append(
        '      "{}", {}, {}, {}, {}, {}, {}, {},'.format(
            name,
            cell.lx,
            cell.ly,
            source,
            disp[0],
            disp[1],
            "true" if is_complex else "false",
            seed,
        )
    )
    out.append("      {" + ", ".join(str(k.code) for k in cell.kinds) + "},")
    out.append(
        "      {"
        + ", ".join("{" + ", ".join(str(x) for x in v) + "}" for v in vd)
        + "},"
    )
    emit_chain(arrays, out, tag, chain, "      ")
    out.append(
        "      {}, {}, {{{}}}, {{{}}},".format(
            ref["nrow"],
            ref["ncol"],
            ", ".join(str(s) for s in ref["psites"]),
            ", ".join(str(p) for p in ref["ppos"]),
        )
    )
    for key in ("before", "after"):
        out.append(
            "      {"
            + ", ".join("{{{}, {}}}".format(fmt(v.real), fmt(v.imag)) for v in ref[key])
            + "},"
        )
    out.append("  },")
    return r_before, r_nostring


def emit_revisit(arrays, out, tag, name, cell, disp, dc, seed):
    chain = chain_for(cell, 0, disp, "hopnn", TAU)
    vd = vdims_for(cell, chain["pos"], chain["sites"], dc)
    assert len(set(chain["sites"])) < len(chain["sites"]), "no site comes back"
    out.append(
        "  // {}: path {} through unit-cell sites {}".format(
            name, " ".join("leg{}".format(l) for l in chain["legs"]), chain["sites"]
        )
    )
    out.append("  {")
    out.append(
        '      "{}", {}, {}, 0, {}, {}, {},'.format(
            name, cell.lx, cell.ly, disp[0], disp[1], seed
        )
    )
    out.append("      {" + ", ".join(str(k.code) for k in cell.kinds) + "},")
    out.append(
        "      {"
        + ", ".join("{" + ", ".join(str(x) for x in v) + "}" for v in vd)
        + "},"
    )
    emit_chain(arrays, out, tag, chain, "      ")
    out.append("  },")


def generate():
    self_check()
    arrays = []
    out = []
    out.append("const std::vector<lg_case_data> lg_su_cases = {")
    for n, case in enumerate(SU_CASES):
        r_before, r_nostring = emit_case(arrays, out, "su{}".format(n), *case)
        name, op = case[0], case[4]
        # the reference must see the operator at all
        assert r_before > 1e-4, (name, r_before)
        # and, where the operator has an odd channel, the path string
        if op in ("hop", "hopnn", "hubbard", "majorana"):
            assert r_nostring > 1e-4, (name, r_nostring)
    out.append("};")
    out.append("const std::vector<lg_revisit_data> lg_revisit_cases = {")
    for n, case in enumerate(REVISIT_CASES):
        emit_revisit(arrays, out, "rv{}".format(n), *case)
    out.append("};")
    head = [
        "// Generated by test/fermion/gen_longrange_gate.py -- do not edit by",
        "// hand; rerun the generator with --write instead.",
        "const int lg_nbra = {};".format(NBRA),
    ]
    return "\n".join(head + arrays + out)


BEGIN = "// ---- BEGIN GENERATED (gen_longrange_gate.py) ----"
END = "// ---- END GENERATED (gen_longrange_gate.py) ----"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    block = generate()
    if not args.write:
        print(block)
        return
    path = os.path.join(HERE, "longrange_gate.cpp")
    with open(path) as f:
        text = f.read()
    head, rest = text.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    with open(path, "w") as f:
        f.write(head + BEGIN + "\n" + block + "\n" + END + tail)


if __name__ == "__main__":
    main()
