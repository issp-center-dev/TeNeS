#!/usr/bin/env python3
# TeNeS - Massively parallel tensor network solver
# Copyright (C) 2019- The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np

LEGS = ("l", "t", "r", "b")
OP_ORDER = ("b", "r", "t", "l")


def label_occ(parity, label):
    if len(parity) == 1:
        return int(parity[label])
    if len(parity) == 2:
        return int(parity[label])
    if len(parity) == 4:
        return int(parity[label])
    raise ValueError(f"unsupported parity-vector length: {len(parity)}")


def deterministic_tensor(site, parities):
    shape = tuple(len(p) for p in parities)
    out = np.zeros(shape, dtype=np.float64)
    for idx in product(*[range(n) for n in shape]):
        if sum(int(parities[ax][idx[ax]]) for ax in range(len(shape))) % 2 == 0:
            x = (site + 2) * (1 + sum((ax + 3) * idx[ax] for ax in range(5)))
            out[idx] = 0.19 * np.sin(x) + 0.13 * np.cos(0.37 * x)
    return out


@dataclass(frozen=True)
class Patch:
    lx: int
    ly: int

    @property
    def nsite(self):
        return self.lx * self.ly

    def site(self, x, y):
        return x + self.lx * y

    def internal_bonds(self):
        bonds = []
        for y in range(self.ly):
            for x in range(self.lx):
                s = self.site(x, y)
                if x + 1 < self.lx:
                    bonds.append((s, "r", self.site(x + 1, y), "l"))
                if y + 1 < self.ly:
                    bonds.append((s, "b", self.site(x, y + 1), "t"))
        return bonds


class Oracle:
    def __init__(self, patch, tensors, leg_parities):
        self.patch = patch
        self.tensors = tensors
        self.leg_parities = leg_parities
        self.mode = {}
        next_mode = 0
        for site in range(patch.nsite):
            self.mode[("f", site)] = next_mode
            next_mode += 1
        for bond_id, (a, aleg, b, bleg) in enumerate(patch.internal_bonds()):
            self.mode[("a", a, aleg)] = next_mode
            next_mode += 1
            self.mode[("a", b, bleg)] = next_mode
            next_mode += 1
        self.nmode = next_mode

    def apply_annihilate(self, vec, mode):
        out = np.zeros_like(vec)
        bit = 1 << mode
        for state, amp in enumerate(vec):
            if amp == 0.0 or (state & bit) == 0:
                continue
            sign = -1.0 if ((state & (bit - 1)).bit_count() % 2) else 1.0
            out[state ^ bit] += sign * amp
        return out

    def apply_create(self, vec, mode):
        out = np.zeros_like(vec)
        bit = 1 << mode
        for state, amp in enumerate(vec):
            if amp == 0.0 or (state & bit) != 0:
                continue
            sign = -1.0 if ((state & (bit - 1)).bit_count() % 2) else 1.0
            out[state | bit] += sign * amp
        return out

    def apply_physical_projector(self, vec, site):
        parity = self.leg_parities[site]
        tensor = self.tensors[site]
        out = np.zeros_like(vec)
        for idx in product(*[range(len(p)) for p in parity]):
            coeff = tensor[idx]
            if coeff == 0.0:
                continue
            term = vec
            for leg in OP_ORDER:
                label = idx[LEGS.index(leg)]
                occ = label_occ(parity[LEGS.index(leg)], label)
                if occ == 0:
                    continue
                key = ("a", site, leg)
                if key not in self.mode:
                    term = np.zeros_like(vec)
                    break
                term = self.apply_annihilate(term, self.mode[key])
            if idx[4] == 1:
                term = self.apply_create(term, self.mode[("f", site)])
            out += coeff * term
        return out

    def apply_bond_creator(self, vec, bond):
        a, aleg, b, bleg = bond
        pa = self.leg_parities[a][LEGS.index(aleg)]
        pb = self.leg_parities[b][LEGS.index(bleg)]
        out = np.zeros_like(vec)
        for label in range(min(len(pa), len(pb))):
            occ_a = label_occ(pa, label)
            occ_b = label_occ(pb, label)
            if occ_a != occ_b:
                continue
            term = vec
            if occ_a:
                term = self.apply_create(term, self.mode[("a", a, aleg)])
                term = self.apply_create(term, self.mode[("a", b, bleg)])
            out += term
        return out

    def state(self):
        vec = np.zeros(1 << self.nmode, dtype=np.float64)
        vec[0] = 1.0
        for bond in self.patch.internal_bonds():
            vec = self.apply_bond_creator(vec, bond)
        for site in range(self.patch.nsite):
            vec = self.apply_physical_projector(vec, site)
        return vec

    def physical_state(self):
        vec = self.state()
        aux_mask = 0
        for key, mode in self.mode.items():
            if key[0] == "a":
                aux_mask |= 1 << mode
        if aux_mask == 0:
            return vec
        out = np.zeros_like(vec)
        for state, amp in enumerate(vec):
            if (state & aux_mask) == 0:
                out[state] = amp
        return out

    def norm(self):
        psi = self.physical_state()
        return float(np.vdot(psi, psi).real)

    def one_body(self, i, j):
        psi = self.physical_state()
        ket = self.apply_annihilate(psi, self.mode[("f", j)])
        ket = self.apply_create(ket, self.mode[("f", i)])
        return float(np.vdot(psi, ket).real)

    def pairing(self, i, j):
        psi = self.physical_state()
        ket = self.apply_annihilate(psi, self.mode[("f", j)])
        ket = self.apply_annihilate(ket, self.mode[("f", i)])
        return float(np.vdot(psi, ket).real)

    def observables(self):
        norm = self.norm()
        n = [self.one_body(i, i) / norm for i in range(self.patch.nsite)]
        bonds = []
        for a, _, b, _ in self.patch.internal_bonds():
            hop = (self.one_body(a, b) + self.one_body(b, a)) / norm
            pair = (self.pairing(a, b) + self.pairing(b, a)) / norm
            bonds.append((a, b, hop, pair))
        return norm, n, bonds


def plain_boson_norm(patch, tensors):
    if patch.nsite == 1:
        return float(np.vdot(tensors[0], tensors[0]).real)
    # Boundary labels are forced to zero by the all-even oracle construction.
    total = 0.0
    for phys in product((0, 1), repeat=patch.nsite):
        amp = 0.0
        ranges = [range(tensors[0].shape[0]) for _ in patch.internal_bonds()]
        for bond_labels in product(*ranges):
            labels = [[0, 0, 0, 0, phys[s]] for s in range(patch.nsite)]
            for blabel, (a, aleg, b, bleg) in zip(bond_labels, patch.internal_bonds()):
                labels[a][LEGS.index(aleg)] = blabel
                labels[b][LEGS.index(bleg)] = blabel
            prod_amp = 1.0
            for site in range(patch.nsite):
                prod_amp *= tensors[site][tuple(labels[site])]
            amp += prod_amp
        total += amp * amp
    return float(total)


def make_case(lx, ly, parity):
    patch = Patch(lx, ly)
    even = [False]
    leg_parities = []
    for y in range(ly):
        for x in range(lx):
            leg_parities.append(
                [
                    parity if x > 0 else even,
                    parity if y > 0 else even,
                    parity if x + 1 < lx else even,
                    parity if y + 1 < ly else even,
                    [False, True],
                ]
            )
    tensors = [
        deterministic_tensor(site, leg_parities[site]) for site in range(patch.nsite)
    ]
    return patch, tensors, leg_parities


def print_case(name, lx, ly, parity):
    patch, tensors, leg_parities = make_case(lx, ly, parity)
    oracle = Oracle(patch, tensors, leg_parities)
    norm, n, bonds = oracle.observables()
    print(f"{name}.norm = {norm:.17e}")
    for i, value in enumerate(n):
        print(f"{name}.n{i} = {value:.17e}")
    for a, b, hop, pair in bonds:
        print(f"{name}.hop{a}{b} = {hop:.17e}")
        print(f"{name}.pair{a}{b} = {pair:.17e}")


def self_check():
    patch = Patch(2, 1)
    parity = [False]
    leg_parities = [
        [parity, parity, parity, parity, [False, True]] for _ in range(patch.nsite)
    ]
    tensors = [
        deterministic_tensor(site, leg_parities[site]) for site in range(patch.nsite)
    ]
    got = Oracle(patch, tensors, leg_parities).norm()
    expected = plain_boson_norm(patch, tensors)
    print(f"self_check.oracle_norm = {got:.17e}")
    print(f"self_check.boson_norm = {expected:.17e}")
    print(f"self_check.diff = {abs(got - expected):.17e}")
    if not np.isclose(got, expected, rtol=0.0, atol=1.0e-14):
        raise SystemExit(1)


def apply_physical_annihilate(bits, mode):
    if bits[mode] == 0:
        return None, 0.0
    sign = -1.0 if sum(bits[:mode]) % 2 else 1.0
    out = list(bits)
    out[mode] = 0
    return tuple(out), sign


def apply_physical_create(bits, mode):
    if bits[mode] != 0:
        return None, 0.0
    sign = -1.0 if sum(bits[:mode]) % 2 else 1.0
    out = list(bits)
    out[mode] = 1
    return tuple(out), sign


def task9_chain_values():
    weights = np.array(
        [
            [[1.10, 0.70], [0.0, 0.0]],
            [[0.90, -0.40], [1.30, 0.80]],
            [[1.20, 0.50], [-0.60, 1.10]],
            [[0.75, -0.35], [1.40, 0.65]],
        ],
        dtype=np.float64,
    )
    amplitudes = {}
    for occ in product((0, 1), repeat=4):
        amp = 0.0
        for b1, b2, b3 in product((0, 1), repeat=3):
            bonds = (0, b1, b2, b3, 0)
            term = 1.0
            for site, n in enumerate(occ):
                left = bonds[site]
                right = bonds[site + 1]
                if right != (left ^ n):
                    term = 0.0
                    break
                term *= weights[site, left, n]
            amp += term
        amplitudes[occ] = amp

    norm = sum(amp * amp for amp in amplitudes.values())
    densities = [
        sum(amp * amp * occ[site] for occ, amp in amplitudes.items()) / norm
        for site in range(4)
    ]
    hop01 = 0.0
    for occ, amp in amplitudes.items():
        for ops in (
            (("ann", 1), ("cre", 0)),
            (("ann", 0), ("cre", 1)),
        ):
            state = occ
            sign = 1.0
            for kind, mode in ops:
                if kind == "ann":
                    state, sgn = apply_physical_annihilate(state, mode)
                else:
                    state, sgn = apply_physical_create(state, mode)
                if state is None:
                    break
                sign *= sgn
            if state is not None:
                hop01 += amplitudes[state] * sign * amp
    return norm, densities, hop01 / norm


def task9_check():
    norm, densities, hop01 = task9_chain_values()
    refs = [
        2.0353391950000002,
        0.47988049112374104,
        0.2202854252015719,
        0.40879752477817338,
        0.14348843338616096,
        -0.32763393032383481,
    ]
    got = [norm, *densities, hop01]
    print("task9_chain.values = " + " ".join(f"{value:.17e}" for value in got))
    print(
        "task9_chain.max_abs_diff = "
        f"{max(abs(a - b) for a, b in zip(got, refs)):.17e}"
    )
    if not np.allclose(got, refs, rtol=0.0, atol=1.0e-14):
        raise SystemExit(1)


def main():
    self_check()
    task9_check()
    print_case("horizontal_2site", 2, 1, [False, True])
    print_case("vertical_2site", 1, 2, [False, True])
    print_case("plaquette_2x2", 2, 2, [False, True])


if __name__ == "__main__":
    main()
