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

"""Constants of test/fermion/unit_d1.cpp (task T4).

Behaviour contract: "振る舞い契約書(T4)" in
docs/superpowers/plans/2026-09-30-fermion-longrange-hamiltonian.md; design
section 7 of docs/superpowers/specs/2026-09-30-fermion-longrange-hamiltonian-
design.md.

Run from anywhere:

    python3 test/fermion/gen_unit_d1.py --write

rewrites the block between the "BEGIN GENERATED" / "END GENERATED" markers of
unit_d1.cpp; without --write it prints the block to stdout.

The unit cells
--------------
The honeycomb and kagome lattices of tool/tenes_simple.py, embedded in the
square lattice as HoneycombLattice / KagomeLattice do it (written out below,
and checked against those classes by ``check_geometry``):

* honeycomb, l x w: an L = 2 l by W = w cell with skew W mod L; sublattice A
  at x = (2 X + y) mod L has virtual dims (l, t, r, b) = (D, 1, D, D), B has
  (D, D, D, 1).  The D = 1 bond is A's top / B's bottom.
* kagome, l x w: a 2 l by 2 w cell; A (even, even) (D, D, D, D), B (odd,
  even) (D, 1, D, 1), C (even, odd) (1, D, 1, D), and the vacancy V (odd,
  odd): physical dimension 1, parity [0], virtual dims (1, 1, 1, 1).

D = 1 legs carry the ledger [even]; a vacancy's physical leg [even].

What is generated, and where each number comes from
---------------------------------------------------
* Contract items 1 and 2 (the simple update).  The gate chains are
  tenes_std's output for a std.toml written here (``Cell.std_toml``: the
  cell above, fermion = true, parity tables, Hamiltonian bonds of
  tenes_simple's lattice up to third neighbours), read by tenes_std.Model;
  the gates of one bond are ``tenes_std.make_evolution`` of that Model's
  Hamiltonian term on that Model's path graph.  They are the INPUT of the
  solver under test.  The reference is the Fock-space construction of
  test/fermion/gen_longrange_gate.py (task T2), copied here as
  ``SparseOracle``: the unit-cell tensors before the update (deterministic
  formula of fock_oracle.py, nonzero on one even and one odd label of every
  leg of dimension > 1, the only label of a D = 1 leg), an open patch around
  the path turned into a Fock state, and O_st applied as a Fock operator,

      O_st = sum evo[i_s, i_t, o_s, o_t] M_{o_s}(s) M_{o_t}(t) P_0(s, t)
                                         M_{i_t}(t)^dag M_{i_s}(s)^dag,

  with evo = expm(-tau H) computed here with scipy from the Fock matrix of
  H.  Nothing of tenes_std's sign formula nor of the simple update enters
  the reference.  The test compares <phi_k|psi'> for six fixed bras with
  <phi_k|O_st psi> up to one overall scalar.
* Contract items 3 and 4 (measurements).  Unit-cell tensors from the same
  deterministic formula on every label (D = 2 legs [e, o], D = 1 legs [e],
  physical [e, o] or [e] for a vacancy).  For every measured window (and
  every correlation chain) the reference is the Fock oracle of the window
  as an open patch:
    - CTM: the environment of the test is CHI = 1, corners 1, and every edge
      diagonal in its fused (ket, bra) leg with entries mu_x^2
      (``ctm_weight``, per site and leg), so the window is the open patch
      with every perimeter leg summed over its labels with weight mu_x^2;
    - MF: the same with lambda (``lambda_value``, per bond).
      Both are the construction of fock_oracle.mf_sum (dangling labels,
      tensors dressed with the weight), which test/fermion/mf_measure.cpp
      pins to the solver's nearest-neighbour mean-field path.
    (An environment fixed to label 0 would not do: a 1 x 2 window is then
    a|00> + b|11>, whose hopping vanishes identically.)
    The overlaps <phi_k|psi> of every window with its perimeter at label 0
    are emitted too; the test's [truth] calibration pins its tensor setup
    with them.
  ``SparseOracle`` is checked against fock_oracle.Oracle itself (patches with
  D = 1 legs, a vacancy, dangling labels) by ``self_check``.

Before writing, the generator checks its own references (see ``generate``):
every window value is either signal (|value| > 1e-4) or a structural zero
(|value| < 1e-15: no D > 1 path joins source and target inside the window);
for every chain O_st psi differs from psi and, where the chain has more than
one hop, from the state the bosonic decomposition would give (the
Jordan-Wigner string along the path dropped), both by more than 1e-4.
Correlation chains that run through a vacancy or along a D = 1 bond are
structural zeros for the same reason (a one-wide chain has no other path).
"""

import argparse
import math
import os
import sys
from itertools import product

import numpy as np
import scipy.linalg
import toml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tool"))
sys.path.insert(0, HERE)

import fock_oracle  # noqa: E402
import tenes_std  # noqa: E402

LEGS = fock_oracle.LEGS  # ("l", "t", "r", "b")
OP_ORDER = fock_oracle.OP_ORDER  # ("b", "r", "t", "l")
NBRA = 6
DIRS = {0: (-1, 0), 1: (0, 1), 2: (1, 0), 3: (0, -1)}
LEG_OF = {v: k for k, v in DIRS.items()}


def popcount(x):
    # int.bit_count() needs Python 3.10; CI runs 3.9.
    return bin(x).count("1")


# ---------------------------------------------------------------------------
# Site kinds
# ---------------------------------------------------------------------------


class Kind:
    """nspin modes; local state i is the occupation bit pattern occ[i] (bit
    k = mode k), created as c^dag_{k0} c^dag_{k1} ... |0> with k0 < k1
    (tool/tenes_simple.py HubbardModel: i = n_up + 2 n_dn)."""

    def __init__(self, code, name, nspin, occ):
        self.code = code
        self.name = name
        self.nspin = nspin
        self.occ = list(occ)
        self.d = len(self.occ)
        self.parity = [popcount(b) % 2 for b in self.occ]


SPINLESS = Kind(0, "S", 1, [0, 1])
HUBBARD = Kind(1, "H", 2, [0, 1, 2, 3])
VACANCY = Kind(2, "V", 0, [0])


# ---------------------------------------------------------------------------
# Unit cells: the honeycomb and kagome embeddings of tenes_simple
# ---------------------------------------------------------------------------


class Cell:
    """A unit cell: size, skew, per-site kind and embedding virtual dims
    (l, t, r, b), and the Hamiltonian bonds tenes_simple would list."""

    def __init__(self, name, lx, ly, skew, kinds, vdims, bonds):
        self.name = name
        self.lx = lx
        self.ly = ly
        self.skew = skew
        self.kinds = kinds
        self.vdims = vdims
        self.bonds = bonds  # list of (source, dx, dy)

    @property
    def n(self):
        return self.lx * self.ly

    def site(self, X, Y):
        """Site at global coordinates (X, Y); T(x, y) = T(x + skew, y + ly)."""
        oy = Y // self.ly
        X -= oy * self.skew
        Y -= oy * self.ly
        return (X % self.lx) + self.lx * Y

    def xy(self, s):
        return s % self.lx, s // self.lx

    def neighbor(self, s, leg):
        x, y = self.xy(s)
        dx, dy = DIRS[leg]
        return self.site(x + dx, y + dy)

    def bond_id(self, s, leg):
        """Canonical id of the unit-cell bond on leg `leg` of site s."""
        if leg in (1, 2):
            return (s, leg)
        return (self.neighbor(s, leg), (leg + 2) % 4)

    def std_toml(self, elements, tau, extra_bonds=()):
        """std.toml text for tenes_std: fermion mode, every Hamiltonian bond
        of the cell with the same two-site elements (dense [i_s, i_t, o_s,
        o_t], written with full precision)."""
        lines = [
            "[parameter]",
            "[parameter.general]",
            "is_real = true",
            "fermion = true",
            "[parameter.simple_update]",
            "tau = {}".format(repr(tau)),
            "num_step = 1",
            "[parameter.full_update]",
            "num_step = 0",
            "",
            "[tensor]",
            "L_sub = [{}, {}]".format(self.lx, self.ly),
            "skew = {}".format(self.skew),
            "",
        ]
        for s in range(self.n):
            k = self.kinds[s]
            lines += [
                "[[tensor.unitcell]]",
                "index = [{}]".format(s),
                "physical_dim = {}".format(k.d),
                "virtual_dim = [{}]".format(", ".join(str(v) for v in self.vdims[s])),
                "parity = [{}]".format(", ".join(str(p) for p in k.parity)),
                "",
            ]
        d = elements.shape[0]
        lines += ["[[hamiltonian]]", "dim = [{}, {}]".format(d, d), 'bonds = """']
        lines += ["{} {} {}".format(*b) for b in self.bonds]
        lines += ["{} {} {}".format(*b) for b in extra_bonds if b not in self.bonds]
        lines += ['"""', 'elements = """']
        for idx in product(*[range(n) for n in elements.shape]):
            v = elements[idx]
            if v != 0.0:
                lines.append("{} {} {} {} {} 0.0".format(*idx, repr(float(np.real(v)))))
        lines += ['"""', ""]
        return "\n".join(lines)


def honeycomb(l, w, vd, kind):
    """tenes_simple.HoneycombLattice with l, w and virtual_dim vd."""
    L, W = 2 * l, w
    skew = W % L
    kinds = [None] * (L * W)
    vdims = [None] * (L * W)
    bonds = []
    for y in range(W):
        for X in range(l):
            a = (2 * X + y) % L + L * y
            b = (2 * X + y + 1) % L + L * y
            kinds[a] = kind
            kinds[b] = kind
            vdims[a] = [vd, 1, vd, vd]
            vdims[b] = [vd, vd, vd, 1]
            bonds += [(a, 1, 0), (b, 1, 0), (b, 0, 1)]
            for s in (a, b):
                bonds += [(s, -1, 1), (s, 1, 1), (s, 2, 0)]
            bonds += [(a, 2, -1), (a, 0, 1), (b, 2, 1)]
    return "honeycomb {}x{}".format(l, w), L, W, skew, kinds, vdims, bonds


def kagome(l, w, vd, kind):
    """tenes_simple.KagomeLattice with l, w and virtual_dim vd."""
    L, W = 2 * l, 2 * w
    kinds = []
    vdims = []
    bonds = []
    for s in range(L * W):
        x, y = s % L, s // L
        if x % 2 == 0 and y % 2 == 0:
            kinds.append(kind)
            vdims.append([vd, vd, vd, vd])
            bonds += [(s, 1, 0), (s, 0, 1), (s, -1, 2), (s, -2, 1)]
            bonds += [(s, 2, 0), (s, 0, 2), (s, -2, 2)]
        elif x % 2 == 1 and y % 2 == 0:
            kinds.append(kind)
            vdims.append([vd, 1, vd, 1])
            bonds += [(s, 1, 0), (s, -1, 1), (s, 1, 1), (s, -1, 2)]
            bonds += [(s, 2, 0), (s, -2, 2), (s, 0, 2)]
        elif x % 2 == 0 and y % 2 == 1:
            kinds.append(kind)
            vdims.append([1, vd, 1, vd])
            bonds += [(s, 0, 1), (s, -1, 1), (s, 1, 1), (s, -2, 1)]
            bonds += [(s, 0, 2), (s, -2, 2), (s, 2, 0)]
        else:
            kinds.append(VACANCY)
            vdims.append([1, 1, 1, 1])
    return "kagome {}x{}".format(l, w), L, W, 0, kinds, vdims, bonds


def make_cell(builder, *args):
    return Cell(*builder(*args))


def check_geometry():
    """The embeddings above are those of tool/tenes_simple.py."""
    import tenes_simple

    for builder, cls, l, w in [
        (honeycomb, tenes_simple.HoneycombLattice, 2, 4),
        (honeycomb, tenes_simple.HoneycombLattice, 2, 2),
        (kagome, tenes_simple.KagomeLattice, 2, 2),
    ]:
        cell = make_cell(builder, l, w, 2, SPINLESS)
        lat = cls({"l": l, "w": w, "virtual_dim": 2})
        assert (cell.lx, cell.ly, cell.skew) == (lat.L, lat.W, lat.skew), cell.name
        for sub in lat.sublattice:
            for s in sub.sites:
                assert cell.vdims[s] == list(sub.vdim), (cell.name, s)
                assert (cell.kinds[s] is VACANCY) == sub.is_vacancy, (cell.name, s)
        mine = sorted(cell.bonds)
        theirs = sorted(
            (b.source, b.dx, b.dy) for level in lat.bonds for typ in level for b in typ
        )
        assert mine == theirs, cell.name


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


def h_hopnn(cd, c):
    """-t (c^dag_s c_t + h.c.) + V n_s n_t, t = 1, V = 0.8 (spinless), or
    the hopping of both spins plus V n_s n_t (Hubbard)."""
    h = 0.8 * _n(cd, c, 0) @ _n(cd, c, 1)
    for s in range(len(cd[0])):
        h = h - (cd[0][s] @ c[1][s] + cd[1][s] @ c[0][s])
    return h


def evolution(ks, kt, tau):
    """evo = expm(-tau H) as [i_s, i_t, o_s, o_t] (scipy, independent of
    tenes_std's eigh route)."""
    H = two_site_operator(ks, kt, h_hopnn)
    d2 = ks.d * kt.d
    mat = H.transpose(2, 3, 0, 1).reshape(d2, d2)  # mat[out, in]
    E = scipy.linalg.expm(-tau * mat)
    return E.reshape(ks.d, kt.d, ks.d, kt.d).transpose(2, 3, 0, 1)


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


def _dot(a, b):
    """<a|b> (a conjugated)."""
    if len(a) > len(b):
        return np.conj(sum(np.conj(v) * a.get(k, 0.0) for k, v in b.items()))
    return sum(np.conj(v) * b.get(k, 0.0) for k, v in a.items())


class SparseOracle:
    """fock_oracle.Oracle on sparse states (dicts), with several physical
    modes per site (Hubbard) or none (vacancy), and the dangling open legs
    of fock_oracle.Oracle (mean-field sums).

    Mode layout, exactly as fock_oracle.Oracle allocates it: the physical
    modes of site 0, of site 1, ... (a Hubbard site owns two, up then down; a
    vacancy none; fock_oracle gives every site one mode, which only renumbers
    modes that are never occupied), then two auxiliary modes per internal
    bond in patch.internal_bonds() order, then one per open leg in
    patch.open_legs() order when dangling labels are given.  The state is
    built as fock_oracle.Oracle.state() does: dangling modes of odd labels
    created in open-leg order, the bond creators in bond order, then the site
    projectors in site order, each annihilating its virtual legs in the order
    b, r, t, l and creating its physical state.
    """

    def __init__(self, kinds, patch, tensors, leg_parities, dangling=None):
        self.kinds = kinds
        self.patch = patch
        self.bonds = patch.internal_bonds()
        self.tensors = tensors
        self.leg_parities = leg_parities
        self.base = []
        nm = 0
        for kind in kinds:
            self.base.append(nm)
            nm += kind.nspin
        self.nphys = nm
        self.mode = {}
        for a, aleg, b, bleg in self.bonds:
            self.mode[(a, aleg)] = nm
            nm += 1
            self.mode[(b, bleg)] = nm
            nm += 1
        self.dangling = None
        self.open_legs = []
        if dangling is not None:
            self.open_legs = patch.open_legs()
            for key in self.open_legs:
                self.mode[key] = nm
                nm += 1
            self.dangling = {key: int(dangling[key]) for key in self.open_legs}
        self.nmode = nm
        self._psi = None

    def create_local(self, vec, site, i):
        kind = self.kinds[site]
        ks = [k for k in range(kind.nspin) if (kind.occ[i] >> k) & 1]
        for k in reversed(ks):
            vec = _create(vec, self.base[site] + k)
        return vec

    def annihilate_local(self, vec, site, i):
        kind = self.kinds[site]
        ks = [k for k in range(kind.nspin) if (kind.occ[i] >> k) & 1]
        for k in ks:
            vec = _annihilate(vec, self.base[site] + k)
        return vec

    def state(self):
        vec = {0: 1.0}
        for site, leg in self.open_legs:
            p = self.leg_parities[site][LEGS.index(leg)]
            if p[self.dangling[(site, leg)]]:
                vec = _create(vec, self.mode[(site, leg)])
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
        fixed = []
        if self.dangling is not None:
            for leg in LEGS:
                if (site, leg) in self.dangling:
                    fixed.append((LEGS.index(leg), self.dangling[(site, leg)]))
        out = {}
        for idx in product(*[range(len(p)) for p in parity]):
            coeff = tensor[idx]
            if coeff == 0.0:
                continue
            if any(idx[axis] != label for axis, label in fixed):
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
        if self._psi is None:
            vec = self.state()
            mask = (1 << self.nphys) - 1
            self._psi = {s: a for s, a in vec.items() if (s & ~mask) == 0 and a != 0.0}
        return self._psi

    def norm(self):
        psi = self.physical_state()
        return float(np.real(_dot(psi, psi)))

    def one_body(self, i, j):
        """<psi| c^dag_i c_j |psi> (mode 0 of the sites)."""
        psi = self.physical_state()
        ket = _annihilate(psi, self.base[j])
        ket = _create(ket, self.base[i])
        return _dot(psi, ket)

    def apply_twosite(self, vec, s, t, evo, string_sites=(), string_sign=False):
        """O_st vec, O_st = sum evo M_o(s) M_o(t) P_0 M_i(t)^dag M_i(s)^dag.

        With string_sign = True each term is multiplied by
        (-1)^{p_c N(string_sites)} (p_c the channel parity): the operator a
        chain composes to when the Jordan-Wigner string along the path is
        dropped (the bosonic decomposition), only used to show that the
        reference is not blind to it."""
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


def _dense_from_fock(ref, nphys):
    want = np.zeros(1 << nphys)
    for s in range(len(ref)):
        if ref[s] != 0.0:
            assert s < (1 << nphys), "fock_oracle state has aux bits"
            want[s] = ref[s]
    return want


def _check_case(lx, ly, vac, d1):
    """Kinds, tensors and ledgers of a self-check patch: d = 2 sites, the
    vacancy `vac` (or None), and D = 1 legs (site, leg) (a D = 1 leg makes
    both ends of an internal bond D = 1)."""
    patch = fock_oracle.Patch(lx, ly)
    internal = {}
    for a, aleg, b, bleg in patch.internal_bonds():
        internal[(a, LEGS.index(aleg))] = (b, LEGS.index(bleg))
        internal[(b, LEGS.index(bleg))] = (a, LEGS.index(aleg))
    leg_parities = []
    kinds = []
    for s in range(patch.nsite):
        lp = []
        for leg in range(4):
            if (s, leg) in d1 or s == vac:
                lp.append([False])
            else:
                lp.append([False, True])
        lp.append([False] if s == vac else [False, True])
        leg_parities.append(lp)
        kinds.append(VACANCY if s == vac else SPINLESS)
    for (s, leg), (t, tleg) in internal.items():
        if len(leg_parities[s][leg]) == 1:
            leg_parities[t][tleg] = [False]
    tensors = [
        fock_oracle.deterministic_tensor(s, leg_parities[s], 3)
        for s in range(patch.nsite)
    ]
    return patch, internal, kinds, tensors, leg_parities


def self_check():
    """SparseOracle against fock_oracle.Oracle (unchanged) on d = 2 patches
    with D = 1 legs and a vacancy: every open leg fixed to label 0 (the CTM
    windows), and every assignment of dangling labels (the mean-field sums;
    small patches only, fock_oracle is dense)."""
    closed_cases = [
        (2, 1, None, []),
        (2, 2, 3, [(3, 0), (2, 2), (3, 1), (1, 3)]),
        (3, 1, None, [(0, 2), (1, 0)]),
        (1, 3, 1, [(1, 1), (0, 3), (1, 3), (2, 1)]),
        (2, 2, None, [(0, 1), (2, 3)]),
    ]
    for lx, ly, vac, d1 in closed_cases:
        patch, internal, kinds, tensors, leg_parities = _check_case(lx, ly, vac, d1)
        closed = [
            [
                lp if (s, leg) in internal else [False]
                for leg, lp in enumerate(leg_parities[s][:4])
            ]
            + [leg_parities[s][4]]
            for s in range(patch.nsite)
        ]
        ctens = [
            tensors[s][
                tuple(
                    slice(None) if (s, leg) in internal else slice(0, 1)
                    for leg in range(4)
                )
                + (slice(None),)
            ]
            for s in range(patch.nsite)
        ]
        ref = fock_oracle.Oracle(patch, ctens, closed)
        sp = SparseOracle(kinds, patch, ctens, closed)
        dense = np.zeros(1 << patch.nsite)
        for st, a in sp.physical_state().items():
            # SparseOracle gives a vacancy no mode: re-insert it
            g = 0
            k = 0
            for site in range(patch.nsite):
                if kinds[site] is VACANCY:
                    continue
                if (st >> k) & 1:
                    g |= 1 << site
                k += 1
            dense[g] = a
        want = _dense_from_fock(ref.physical_state(), patch.nsite)
        err = np.max(np.abs(dense - want))
        scale = np.max(np.abs(want))
        assert scale > 1e-5, ("empty self-check state", lx, ly)
        assert err <= 1e-13 * scale, ("SparseOracle vs fock_oracle", lx, ly, err)
    dangling_cases = [
        (2, 1, None, []),
        (2, 1, None, [(0, 2), (0, 1)]),
        (1, 2, 1, []),
        (1, 2, None, [(0, 3), (1, 0)]),
    ]
    for lx, ly, vac, d1 in dangling_cases:
        patch, internal, kinds, tensors, leg_parities = _check_case(lx, ly, vac, d1)
        open_legs = patch.open_legs()
        ranges = [range(len(leg_parities[s][LEGS.index(leg)])) for s, leg in open_legs]
        pos = [k for k in range(patch.nsite) if kinds[k] is not VACANCY]
        nconf = 0
        for labels in product(*ranges):
            x = dict(zip(open_legs, labels))
            ref = fock_oracle.Oracle(patch, tensors, leg_parities, dangling_labels=x)
            sp = SparseOracle(kinds, patch, tensors, leg_parities, dangling=x)
            assert abs(ref.norm() - sp.norm()) < 1e-14
            for i in pos:
                for j in pos:
                    a = ref.one_body(i, j)
                    b = np.real(sp.one_body(i, j))
                    assert abs(a - b) < 1e-14, ("dangling", lx, ly, x, i, j, a, b)
            nconf += 1
        assert nconf > 1


# ---------------------------------------------------------------------------
# Deterministic tensors and bras (mirrored in unit_d1.cpp)
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
    """The labels a leg carries nonzeros on in the simple-update cases: the
    first even and the first odd index of its even-first ledger (the only
    label of a D = 1 leg)."""
    if dim == 1:
        return [0]
    return [0, (dim + 1) // 2]


def ctm_weight(cell, s, leg):
    """CTM edge weights of the measurement cases (mirrored in unit_d1.cpp):
    the edge tensor that closes leg `leg` of site s (eTl, eTt, eTr, eTb for
    leg 0, 1, 2, 3) is diagonal in its (ket, bra) pair with entries mu_x^2.
    Per (site, leg), not per bond: the two ends of a bond are closed by
    different edge tensors."""
    D = cell.vdims[s][leg]
    if D == 1:
        return [0.8]
    return [1.0, 0.45 + 0.029 * s + 0.07 * leg]


def lambda_value(cell, s, leg):
    """Mean-field weights of the measurement cases (mirrored in unit_d1.cpp):
    one vector per unit-cell bond, the same at both ends."""
    b, bleg = cell.bond_id(s, leg)
    D = cell.vdims[s][leg]
    if D == 1:
        return [0.9]
    return [1.0, 0.35 + 0.037 * b + (0.11 if bleg == 1 else 0.0)]


# ---------------------------------------------------------------------------
# Contract items 1 and 2: gate chains through the simple update
# ---------------------------------------------------------------------------


class StdModel:
    """tenes_std.Model of a cell's std.toml (kept for its Hamiltonian terms
    and its path graph)."""

    def __init__(self, cell, kind, extra_bonds=()):
        H = two_site_operator(kind, kind, h_hopnn)
        self.text = cell.std_toml(H, 0.1, extra_bonds)
        self.model = tenes_std.Model(toml.loads(self.text))
        self.cell = cell

    def gates(self, source, disp, tau):
        for ham in self.model.hamiltonians:
            if not isinstance(ham, tenes_std.NNOperator):
                continue
            b = ham.bond
            if (b.source_site, b.dx, b.dy) == (source, disp[0], disp[1]):
                return tenes_std.make_evolution(
                    ham, self.model.graph, tau, fermion=True
                )
        raise KeyError((source, disp))


def chain_for(std, source, disp, tau):
    cell = std.cell
    gates = std.gates(source, disp, tau)
    pos = [(0, 0)]
    sites = [source]
    legs = []
    X0, Y0 = cell.xy(source)
    for g in gates:
        leg = std.model.unitcell.bond_direction(g.bond)
        legs.append(leg)
        assert int(g.bond.source_site) == sites[-1]
        dx, dy = DIRS[leg]
        pos.append((pos[-1][0] + dx, pos[-1][1] + dy))
        sites.append(cell.site(X0 + pos[-1][0], Y0 + pos[-1][1]))
        # a Hamiltonian path never crosses a D = 1 bond
        assert cell.vdims[sites[-2]][leg] > 1
    assert pos[-1] == tuple(disp), (pos, disp)
    target = sites[-1]
    # expected ledgers, walking the chain (design section 6.1)
    current = {s: list(cell.kinds[s].parity) for s in range(cell.n)}
    ledgers = []
    for k, g in enumerate(gates):
        a, b = sites[k], sites[k + 1]
        G = np.asarray(g.elements)
        p_in1 = current[a]
        p_in2 = current[b]
        p_out1 = list(cell.kinds[a].parity)
        assert G.shape[:3] == (len(p_in1), len(p_in2), len(p_out1))
        seen = [set() for _ in range(G.shape[3])]
        for i1, i2, o1, o2 in np.argwhere(G != 0):
            seen[o2].add((p_in1[i1] + p_in2[i2] + p_out1[o1]) % 2)
        assert all(len(s) == 1 for s in seen), "out2 not uniquely graded"
        p_out2 = [s.pop() for s in seen]
        ledgers.append((list(p_in1), list(p_in2), p_out1, p_out2))
        current[a] = p_out1
        current[b] = p_out2
    for s in range(cell.n):
        assert current[s] == list(cell.kinds[s].parity), "chain does not close"
    evo = evolution(cell.kinds[source], cell.kinds[target], tau)
    return dict(gates=gates, legs=legs, pos=pos, sites=sites, ledgers=ledgers, evo=evo)


def vdims_for(cell, chain, dc):
    """Test virtual dims: the embedding dims (D = 1 legs stay 1, the others
    2), dc on the unit-cell bonds of the path."""
    vd = [[1 if v == 1 else 2 for v in cell.vdims[s]] for s in range(cell.n)]
    for k, leg in enumerate(chain["legs"]):
        a, b = chain["sites"][k], chain["sites"][k + 1]
        vd[a][leg] = dc
        vd[b][(leg + 2) % 4] = dc
    return vd


def cell_tensor_value(site, seed, small, big_parities, is_complex):
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
    X0, Y0 = cell.xy(chain["sites"][0])
    patch = fock_oracle.Patch(ncol, nrow)
    psites = []  # unit-cell site at each patch position (raster, row 0 on top)
    for row in range(nrow):
        for col in range(ncol):
            psites.append(cell.site(X0 + xmin + col, Y0 + ymax - row))
    assert len(set(psites)) == len(psites), "a unit-cell site repeats in the patch"
    ppos = [(ymax - y) * ncol + (x - xmin) for x, y in chain["pos"]]
    # The patch state after the chain is O_st psi only if no gate of the
    # chain (or of a translated copy) acts on a patch bond other than the
    # path bonds, and each path bond once.
    gate_bonds = [cell.bond_id(s, l) for s, l in zip(chain["sites"], chain["legs"])]
    assert len(set(gate_bonds)) == len(gate_bonds)
    path_pairs = set()
    for k in range(len(ppos) - 1):
        path_pairs.add((ppos[k], ppos[k + 1]))
        path_pairs.add((ppos[k + 1], ppos[k]))
    for p, s in enumerate(psites):
        for leg in range(4):
            dx, dy = DIRS[leg]
            q = None
            r, c = p // ncol, p % ncol
            rr, cc = r - dy, c + dx
            if 0 <= rr < nrow and 0 <= cc < ncol:
                q = rr * ncol + cc
            if q is not None and (p, q) in path_pairs:
                continue
            assert cell.bond_id(s, leg) not in gate_bonds, (
                "a gate bond is a non-path patch bond",
                s,
                leg,
            )
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
            big.append([even_first(dim)[lab] for lab in support_labels(dim)])
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
    oracle = SparseOracle(kinds, patch, tensors, leg_parities)
    psi = oracle.physical_state()
    s_pos, t_pos = ppos[0], ppos[-1]
    after = oracle.apply_twosite(psi, s_pos, t_pos, chain["evo"])
    nostring = oracle.apply_twosite(
        psi, s_pos, t_pos, chain["evo"], string_sites=ppos[1:-1], string_sign=True
    )
    locals_ = list(product(*[range(k.d) for k in kinds]))

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

    return dict(
        nrow=nrow,
        ncol=ncol,
        psites=psites,
        ppos=ppos,
        before=overlaps(psi),
        after=overlaps(after),
        nostring=overlaps(nostring),
    )


def residual(got, want):
    c = np.vdot(want, got) / np.vdot(want, want)
    return float(np.linalg.norm(got - c * want) / np.linalg.norm(got))


# ---------------------------------------------------------------------------
# Contract items 3 and 4: measurements
# ---------------------------------------------------------------------------


def measure_tensor(cell, s, seed):
    """Unit-cell tensor of the measurement cases: det on every parity-even
    element (embedding dims, ledgers [e, o] / [e], physical [e, o] / [e])."""
    lp = [[bool(x) for x in even_first(v)] for v in cell.vdims[s]]
    lp.append([bool(x) for x in cell.kinds[s].parity])
    shape = tuple(len(x) for x in lp)
    T = np.zeros(shape)
    for idx in product(*[range(n) for n in shape]):
        if sum(lp[ax][idx[ax]] for ax in range(5)) % 2 == 0:
            T[idx] = det(s, seed, idx)
    return T, lp


def window_sites(cell, source, cols, rows, scol, srow):
    """Unit-cell site of each window cell (row 0 on top), source at (srow,
    scol): cell (row, col) = lattice.other(source, col - scol, srow - row)."""
    X0, Y0 = cell.xy(source)
    return [
        cell.site(X0 + col - scol, Y0 + srow - row)
        for row in range(rows)
        for col in range(cols)
    ]


def window_oracles(cell, psites, ncol, nrow, seed, weight):
    """Yield a SparseOracle for every assignment of labels to the perimeter
    legs of a window (fock_oracle's dangling labels), each perimeter leg of
    site s dressed with weight(cell, s, leg) (so that a label contributes
    its weight squared, as in fock_oracle.mf_sum); weight = None fixes every
    perimeter leg to label 0 instead (one oracle)."""
    patch = fock_oracle.Patch(ncol, nrow)
    internal = set()
    for a, aleg, b, bleg in patch.internal_bonds():
        internal.add((a, aleg))
        internal.add((b, bleg))
    kinds = [cell.kinds[s] for s in psites]
    tensors = []
    lps = []
    for p, s in enumerate(psites):
        T, lp = measure_tensor(cell, s, seed)
        if weight is not None:
            for leg in range(4):
                if (p, LEGS[leg]) in internal:
                    continue
                lam = np.asarray(weight(cell, s, leg))
                shape = [1] * 5
                shape[leg] = len(lam)
                T = T * lam.reshape(shape)
        else:
            sl = []
            for leg in range(4):
                if (p, LEGS[leg]) in internal:
                    sl.append(slice(None))
                else:
                    sl.append(slice(0, 1))
                    lp[leg] = [False]
            T = T[tuple(sl) + (slice(None),)]
        tensors.append(T)
        lps.append(lp)
    if weight is None:
        yield SparseOracle(kinds, patch, tensors, lps)
        return
    open_legs = patch.open_legs()
    ranges = [range(len(lps[s][LEGS.index(leg)])) for s, leg in open_legs]
    for labels in product(*ranges):
        yield SparseOracle(
            kinds, patch, tensors, lps, dangling=dict(zip(open_legs, labels))
        )


def window_value(cell, psites, ncol, nrow, sp, tp, seed, weight):
    """<c^dag_s c_t> of the window (patch positions sp, tp) and a scale that
    does not cancel (sum of |terms| / norm)."""
    num = 0.0
    den = 0.0
    absnum = 0.0
    for orc in window_oracles(cell, psites, ncol, nrow, seed, weight):
        assert orc.kinds[sp].nspin == 1 and orc.kinds[tp].nspin == 1
        v = np.real(orc.one_body(sp, tp))
        num += v
        absnum += abs(v)
        den += orc.norm()
    return num / den, absnum / den, den


def window_overlaps(cell, psites, ncol, nrow, seed):
    """<phi_k|psi> of the CTM window (perimeter fixed to label 0)."""
    orc = next(window_oracles(cell, psites, ncol, nrow, seed, None))
    psi = orc.physical_state()
    locals_ = list(product(*[range(k.d) for k in orc.kinds]))
    out = []
    for k in range(NBRA):
        acc = 0.0
        for loc in locals_:
            amp = psi.get(orc.bits(loc), 0.0)
            if amp != 0.0:
                acc += bra(k, loc) * amp
        out.append(acc)
    return out


def pair_window(cell, source, dx, dy):
    ncol = abs(dx) + 1
    nrow = abs(dy) + 1
    scol = 0 if dx >= 0 else ncol - 1
    srow = nrow - 1 if dy >= 0 else 0
    tcol = ncol - 1 - scol
    trow = nrow - 1 - srow
    psites = window_sites(cell, source, ncol, nrow, scol, srow)
    assert len(set(psites)) == len(psites)
    return psites, ncol, nrow, srow * ncol + scol, trow * ncol + tcol


def x_first_crossing(cell, source, dx, dy):
    """What the x_first relay path from source to target crosses: the D = 1
    bonds (as (site, leg)) and the vacancies strictly between the ends."""
    X, Y = cell.xy(source)
    s = source
    d1 = []
    vac = []
    steps = [(1 if dx > 0 else -1, 0)] * abs(dx) + [(0, 1 if dy > 0 else -1)] * abs(dy)
    for k, (ddx, ddy) in enumerate(steps):
        leg = LEG_OF[(ddx, ddy)]
        if cell.vdims[s][leg] == 1:
            d1.append((s, leg))
        X += ddx
        Y += ddy
        s = cell.site(X, Y)
        if k + 1 < len(steps) and cell.kinds[s] is VACANCY:
            vac.append(s)
    return d1, vac


def chain_value(cell, left, r, vertical, seed, weight, left_dagger):
    """Correlation <A_left B_right> of the solver's chain, right = left
    moved r times right (up when vertical); A = c^dag, B = c (left_dagger)
    or A = c, B = c^dag = -<c^dag_right c_left>."""
    X, Y = cell.xy(left)
    if vertical:
        psites = [cell.site(X, Y + r - row) for row in range(r + 1)]
        ncol, nrow = 1, r + 1
        lp, rp = r, 0
    else:
        psites = [cell.site(X + col, Y) for col in range(r + 1)]
        ncol, nrow = r + 1, 1
        lp, rp = 0, r
    if left_dagger:
        v, sc, _ = window_value(cell, psites, ncol, nrow, lp, rp, seed, weight)
        return v, sc
    v, sc, _ = window_value(cell, psites, ncol, nrow, rp, lp, seed, weight)
    return -v, sc


# ---------------------------------------------------------------------------
# Case lists
# ---------------------------------------------------------------------------

HC44 = make_cell(honeycomb, 2, 4, 2, SPINLESS)  # 4x4, skew 0
HC42 = make_cell(honeycomb, 2, 2, 2, SPINLESS)  # 4x2, skew 2
KG44 = make_cell(kagome, 2, 2, 2, SPINLESS)  # 4x4 with vacancies
KG44H = make_cell(kagome, 2, 2, 2, HUBBARD)  # the same with Hubbard sites
CELLS = [HC44, HC42, KG44, KG44H]
CELL_KIND = {id(HC44): SPINLESS, id(HC42): SPINLESS, id(KG44): SPINLESS}
CELL_KIND[id(KG44H)] = HUBBARD

TAU = 0.1
TAU_C = complex(0.1, 0.05)

# (name, cell, source, disp, tau, is_complex, dc, seed)
# Contract item 1: one nearest-neighbour gate next to D = 1 legs.
NN_CASES = [
    # honeycomb: every site has a D = 1 leg (A top, B bottom): both ends
    ("nn honeycomb A->B (1,0)", HC44, 0, (1, 0), TAU, False, 16, 3),
    ("nn honeycomb B->A (1,0)", HC44, 1, (1, 0), TAU, False, 16, 5),
    ("nn honeycomb B->A (0,1)", HC44, 1, (0, 1), TAU, False, 16, 7),
    # source behind the target in raster order (source_leg 0 and 3)
    ("nn honeycomb A->B (-1,0)", HC44, 5, (-1, 0), TAU, False, 16, 9),
    ("nn honeycomb A->B (0,-1)", HC44, 5, (0, -1), TAU, False, 16, 11),
    # kagome: A has no D = 1 leg, B / C have two: one end only
    ("nn kagome A->B (1,0)", KG44, 0, (1, 0), TAU, False, 16, 13),
    ("nn kagome A->C (0,1)", KG44, 0, (0, 1), TAU, False, 16, 15),
    ("nn kagome B->A (1,0)", KG44, 1, (1, 0), TAU, False, 16, 17),
    ("nn kagome C->A (0,1)", KG44, 4, (0, 1), TAU, False, 16, 19),
    ("nn kagome B->A (-1,0)", KG44, 11, (-1, 0), TAU, False, 16, 21),
    # Hubbard (d = 4) next to D = 1 legs
    ("nn hubbard kagome B->A (1,0)", KG44H, 1, (1, 0), TAU, False, 16, 23),
]
# Contract item 2: chains along D > 1 bonds next to D = 1 legs / vacancies.
CHAIN_CASES = [
    ("chain honeycomb A NNN (1,1)", HC44, 0, (1, 1), TAU, False, 16, 31),
    ("chain honeycomb A NNN (-1,1)", HC44, 0, (-1, 1), TAU, False, 16, 33),
    ("chain honeycomb A NNN (2,0)", HC44, 5, (2, 0), TAU, False, 16, 35),
    ("chain honeycomb B NNN (-1,1)", HC44, 1, (-1, 1), TAU, False, 16, 37),
    ("chain honeycomb B NNN (1,1)", HC44, 6, (1, 1), TAU, False, 16, 39),
    ("chain honeycomb A third (0,1)", HC44, 0, (0, 1), TAU, False, 16, 41),
    ("chain honeycomb skew A NNN (1,1)", HC42, 5, (1, 1), TAU, False, 16, 43),
    ("chain honeycomb skew B NNN (-1,1)", HC42, 4, (-1, 1), TAU, False, 16, 45),
    ("chain kagome NN B->C (-1,1)", KG44, 1, (-1, 1), TAU, False, 16, 47),
    ("chain kagome NN C->B (-1,1)", KG44, 4, (-1, 1), TAU, False, 16, 49),
    ("chain kagome NN B->C (-1,1) far", KG44, 11, (-1, 1), TAU, False, 16, 51),
    ("chain kagome NNN B (1,1)", KG44, 1, (1, 1), TAU, False, 16, 53),
    ("chain complex kagome NN C->B (-1,1)", KG44, 14, (-1, 1), TAU_C, True, 16, 55),
    ("chain hubbard kagome NN B->C (-1,1)", KG44H, 1, (-1, 1), TAU, False, 32, 57),
]
SU_CASES = NN_CASES + CHAIN_CASES

MEASURE_SEED = 7
# (cell, sources, displacements); every (source, disp) window is measured
MEASURE_SETS = [
    (
        HC44,
        [0, 1, 5, 6],
        [(1, 0), (0, 1), (-1, 0), (0, -1), (1, 1), (-1, 1), (2, 0), (1, -1), (2, 1)],
    ),
    (
        KG44,
        [0, 1, 4, 10, 11, 14],
        [(1, 0), (0, 1), (-1, 0), (0, -1), (-1, 1), (1, -1), (1, 1)],
    ),
]
CORR_RMAX = 3


def fmt(v):
    return repr(float(v))


def cpp_bool_string(p):
    return '"' + "".join("1" if x else "0" for x in p) + '"'


def emit_entries(arrays, array_name, G):
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
    arrays.append("const u1_entry {}[] = {{".format(array_name))
    for i in range(0, len(items), 3):
        arrays.append("    " + ",".join(items[i : i + 3]) + ",")
    arrays.append("};")


def cell_index(cell):
    return CELLS.index(cell)


def emit_cells(out):
    out.append("const std::vector<u1_cell_data> u1_cells = {")
    for cell in CELLS:
        out.append(
            '  {{"{}", {}, {}, {}, {{{}}},'.format(
                cell.name + (" hubbard" if cell is KG44H else ""),
                cell.lx,
                cell.ly,
                cell.skew,
                ", ".join(str(k.code) for k in cell.kinds),
            )
        )
        out.append(
            "   {"
            + ", ".join("{" + ", ".join(str(x) for x in v) + "}" for v in cell.vdims)
            + "}},"
        )
    out.append("};")


def emit_su_case(arrays, out, tag, stds, case):
    name, cell, source, disp, tau, is_complex, dc, seed = case
    chain = chain_for(stds[id(cell)], source, disp, tau)
    vd = vdims_for(cell, chain, dc)
    ref = reference(cell, chain, vd, seed, is_complex)
    r_before = residual(ref["after"], ref["before"])
    r_nostring = residual(ref["after"], ref["nostring"])
    assert r_before > 1e-4, (name, r_before)
    if len(chain["gates"]) > 1:
        assert r_nostring > 1e-4, (name, r_nostring)
    out.append(
        "  // {}: path {} through unit-cell sites {}".format(
            name, " ".join("leg{}".format(x) for x in chain["legs"]), chain["sites"]
        )
    )
    out.append(
        "  //   residual(after vs before) = {:.3e}, "
        "residual(after vs string dropped) = {:.3e}".format(r_before, r_nostring)
    )
    out.append("  {")
    out.append(
        '      "{}", {}, {}, {}, {}, {}, {},'.format(
            name,
            cell_index(cell),
            source,
            disp[0],
            disp[1],
            "true" if is_complex else "false",
            seed,
        )
    )
    out.append(
        "      {"
        + ", ".join("{" + ", ".join(str(x) for x in v) + "}" for v in vd)
        + "},"
    )
    out.append("      {")
    for k, (g, leg, led) in enumerate(
        zip(chain["gates"], chain["legs"], chain["ledgers"])
    ):
        G = np.asarray(g.elements)
        array_name = "u1_e_{}_{}".format(tag, k)
        emit_entries(arrays, array_name, G)
        out.append(
            "          {{{}, {}, {{{}}}, {{{}}}, {}, sizeof({}) / sizeof(u1_entry)}},".format(
                int(g.bond.source_site),
                leg,
                ", ".join(str(x) for x in G.shape),
                ", ".join(cpp_bool_string(p) for p in led),
                array_name,
                array_name,
            )
        )
    out.append("      },")
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


MIN_SIGNAL = 1e-4


def emit_measure(out, stds):
    out.append("const int u1_measure_seed = {};".format(MEASURE_SEED))
    out.append("const int u1_corr_rmax = {};".format(CORR_RMAX))
    out.append("const std::vector<u1_measure_set> u1_measure_sets = {")
    for cell, sources, disps in MEASURE_SETS:
        out.append("  {{{}, {{".format(cell_index(cell)))
        for s in sources:
            for dx, dy in disps:
                psites, ncol, nrow, sp, tp = pair_window(cell, s, dx, dy)
                if cell.kinds[psites[tp]] is VACANCY:
                    continue
                ctm, ctm_sc, _ = window_value(
                    cell, psites, ncol, nrow, sp, tp, MEASURE_SEED, ctm_weight
                )
                mf, mf_sc, _ = window_value(
                    cell, psites, ncol, nrow, sp, tp, MEASURE_SEED, lambda_value
                )
                # Either signal, or a structural zero (source and target not
                # connected by any D > 1 path inside the window).
                for v in (ctm, mf):
                    assert abs(v) > MIN_SIGNAL or abs(v) < 1e-15, (
                        cell.name,
                        s,
                        dx,
                        dy,
                        v,
                    )
                ov = window_overlaps(cell, psites, ncol, nrow, MEASURE_SEED)
                d1, vac = x_first_crossing(cell, s, dx, dy)
                out.append(
                    "      // source {} ({},{}): window sites {}; x_first crosses "
                    "D = 1 bonds {} and vacancies {}".format(s, dx, dy, psites, d1, vac)
                )
                out.append(
                    "      {{{}, {}, {}, {}, {}, {}, {}, {}, {}, {{{}}}}},".format(
                        s,
                        dx,
                        dy,
                        fmt(ctm),
                        fmt(ctm_sc),
                        fmt(mf),
                        fmt(mf_sc),
                        len(d1),
                        len(vac),
                        ", ".join(fmt(v) for v in ov),
                    )
                )
        out.append("  }, {")
        for left in range(cell.n):
            if cell.kinds[left] is VACANCY:
                continue
            for vertical in (False, True):
                X, Y = cell.xy(left)
                for r in range(1, CORR_RMAX + 1):
                    right = cell.site(X, Y + r) if vertical else cell.site(X + r, Y)
                    if cell.kinds[right] is VACANCY:
                        continue
                    for lop, rop, dag in ((0, 1, True), (1, 0, False)):
                        ctm, ctm_sc = chain_value(
                            cell, left, r, vertical, MEASURE_SEED, ctm_weight, dag
                        )
                        mf, mf_sc = chain_value(
                            cell, left, r, vertical, MEASURE_SEED, lambda_value, dag
                        )
                        out.append(
                            "      {{{}, {}, {}, {}, {}, {}, {}, {}, {}}},".format(
                                left,
                                0 if vertical else r,
                                r if vertical else 0,
                                lop,
                                rop,
                                fmt(ctm),
                                fmt(ctm_sc),
                                fmt(mf),
                                fmt(mf_sc),
                            )
                        )
        out.append("  }},")
    out.append("};")


def generate():
    check_geometry()
    self_check()
    stds = {}
    for cell in CELLS:
        extra = [(c[2], c[3][0], c[3][1]) for c in SU_CASES if c[1] is cell]
        stds[id(cell)] = StdModel(cell, CELL_KIND[id(cell)], extra)
    arrays = []
    out = []
    emit_cells(out)
    out.append("const std::vector<u1_su_case> u1_su_cases = {")
    for n, case in enumerate(SU_CASES):
        emit_su_case(arrays, out, "su{}".format(n), stds, case)
    out.append("};")
    emit_measure(out, stds)
    head = [
        "// Generated by test/fermion/gen_unit_d1.py -- do not edit by hand;",
        "// rerun the generator with --write instead.",
        "const int u1_nbra = {};".format(NBRA),
    ]
    return "\n".join(head + arrays + out)


BEGIN = "// ---- BEGIN GENERATED (gen_unit_d1.py) ----"
END = "// ---- END GENERATED (gen_unit_d1.py) ----"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument(
        "--std-toml",
        metavar="DIR",
        help="also write the std.toml of every cell into DIR (for inspection)",
    )
    args = parser.parse_args()
    if args.std_toml:
        os.makedirs(args.std_toml, exist_ok=True)
        for cell in CELLS:
            text = cell.std_toml(
                two_site_operator(CELL_KIND[id(cell)], CELL_KIND[id(cell)], h_hopnn),
                0.1,
            )
            fname = cell.name.replace(" ", "_") + ("_hubbard" if cell is KG44H else "")
            with open(os.path.join(args.std_toml, fname + ".std.toml"), "w") as f:
                f.write(text)
    block = generate()
    if not args.write:
        print(block)
        return
    path = os.path.join(HERE, "unit_d1.cpp")
    with open(path) as f:
        text = f.read()
    head, rest = text.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    with open(path, "w") as f:
        f.write(head + BEGIN + "\n" + block + "\n" + END + tail)


if __name__ == "__main__":
    main()
