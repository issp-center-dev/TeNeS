#!/usr/bin/env python3

import itertools

import numpy as np

PHYS_PARITY = np.array([0, 1], dtype=int)
VIRT_PARITY = np.array([0, 1], dtype=int)

# weights[site, left_virtual_parity, physical_occupation]
WEIGHTS = np.array(
    [
        [[1.10, 0.70], [0.0, 0.0]],
        [[0.90, -0.40], [1.30, 0.80]],
        [[1.20, 0.50], [-0.60, 1.10]],
        [[0.75, -0.35], [1.40, 0.65]],
    ],
    dtype=float,
)


def make_tensor(site):
    left_dim = 1 if site == 0 else 2
    right_dim = 1 if site == 3 else 2
    tensor = np.zeros((left_dim, 2, right_dim), dtype=float)
    left_values = [0] if site == 0 else [0, 1]
    right_values = [0] if site == 3 else [0, 1]
    for li, left in enumerate(left_values):
        for occ in [0, 1]:
            right = left ^ occ
            if right in right_values:
                ri = right_values.index(right)
                tensor[li, occ, ri] = WEIGHTS[site, left, occ]
    return tensor


def state_vector():
    tensors = [make_tensor(site) for site in range(4)]
    psi = np.zeros(16, dtype=float)
    for occs in itertools.product([0, 1], repeat=4):
        amp = 0.0
        for b0 in [0, 1]:
            for b1 in [0, 1]:
                for b2 in [0, 1]:
                    amp += (
                        tensors[0][0, occs[0], b0]
                        * tensors[1][b0, occs[1], b1]
                        * tensors[2][b1, occs[2], b2]
                        * tensors[3][b2, occs[3], 0]
                    )
        state = sum(occs[site] << site for site in range(4))
        psi[state] = amp
    return psi


def popcount(state):
    return sum((state >> i) & 1 for i in range(4))


def annihilate(state, site):
    if ((state >> site) & 1) == 0:
        return None, 0.0
    sign = -1.0 if popcount(state & ((1 << site) - 1)) % 2 else 1.0
    return state & ~(1 << site), sign


def create(state, site):
    if ((state >> site) & 1) == 1:
        return None, 0.0
    sign = -1.0 if popcount(state & ((1 << site) - 1)) % 2 else 1.0
    return state | (1 << site), sign


def hopping_matrix(site_a, site_b):
    mat = np.zeros((16, 16), dtype=float)
    for state in range(16):
        mid, s1 = annihilate(state, site_b)
        if mid is not None:
            out, s2 = create(mid, site_a)
            if out is not None:
                mat[out, state] += s1 * s2
        mid, s1 = annihilate(state, site_a)
        if mid is not None:
            out, s2 = create(mid, site_b)
            if out is not None:
                mat[out, state] += s1 * s2
    return mat


def main():
    psi = state_vector()
    norm = float(np.dot(psi, psi))
    occ = []
    for site in range(4):
        values = np.array([(state >> site) & 1 for state in range(16)], dtype=float)
        occ.append(float(np.dot(psi * values, psi) / norm))
    hop01 = float(psi @ hopping_matrix(0, 1) @ psi / norm)

    print(f"norm = {norm:.17g}")
    for site, value in enumerate(occ):
        print(f"n{site} = {value:.17g}")
    print(f"hop01 = {hop01:.17g}")


if __name__ == "__main__":
    main()
