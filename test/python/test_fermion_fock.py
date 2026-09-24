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

import os
import sys

import numpy as np
import pytest

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tool")
)

import tenes_simple


def cd(m, M):
    return tenes_simple.fock_cop(True, m, M)


def c(m, M):
    return tenes_simple.fock_cop(False, m, M)


def test_fermion_modes():
    assert tenes_simple.fermion_modes(1) == 2
    assert tenes_simple.fermion_modes(2) == 4


def test_anticommutation_relations():
    M = 4
    dim = 1 << M
    for m in range(M):
        for n in range(M):
            expected = np.eye(dim) if m == n else np.zeros((dim, dim))
            assert np.allclose(cd(m, M) @ c(n, M) + c(n, M) @ cd(m, M), expected)
            assert np.allclose(c(m, M) @ c(n, M) + c(n, M) @ c(m, M), 0.0)
            assert np.allclose(cd(m, M) @ cd(n, M) + cd(n, M) @ cd(m, M), 0.0)


def test_local_index_to_occupation_spinless():
    assert tenes_simple.local_index_to_occupation(0, 1) == [0]
    assert tenes_simple.local_index_to_occupation(1, 1) == [1]


def test_local_index_to_occupation_spinful():
    # |0>, |up>, |dn>, |up dn>  ->  i = n_up + 2 n_dn
    assert tenes_simple.local_index_to_occupation(0, 2) == [0, 0]
    assert tenes_simple.local_index_to_occupation(1, 2) == [1, 0]
    assert tenes_simple.local_index_to_occupation(2, 2) == [0, 1]
    assert tenes_simple.local_index_to_occupation(3, 2) == [1, 1]


def test_bond_matrix_leg_order_is_in_in_out_out():
    # c^dag_{1} c_{2} is not hermitian, so a transposed leg order would show up.
    M = 2
    t = tenes_simple.bond_matrix(cd(0, M) @ c(1, M), 1)
    # <1, 0| c^dag_1 c_2 |0, 1> = +1
    assert t[0, 1, 1, 0] == pytest.approx(1.0)
    # the adjoint element must NOT be here
    assert t[1, 0, 0, 1] == pytest.approx(0.0)


def test_spinless_hopping_matches_the_handwritten_sample():
    # sample/07_spinless_fermion/input.toml, [[observable.twosite]] bond_hamiltonian
    M = 2
    h = -(cd(0, M) @ c(1, M) + cd(1, M) @ c(0, M))
    t = tenes_simple.bond_matrix(h, 1)
    expected = np.zeros((2, 2, 2, 2))
    expected[0, 1, 1, 0] = -1.0
    expected[1, 0, 0, 1] = -1.0
    assert np.allclose(t, expected)


def test_known_sign_for_spinful_hopping():
    # <up dn, 0| c^dag_{1 up} c_{2 up} |dn, up> = -1
    M = 4
    t = tenes_simple.bond_matrix(cd(0, M) @ c(2, M), 2)
    assert t[2, 1, 3, 0] == pytest.approx(-1.0)


def test_doublon_is_cdag_up_cdag_dn_on_the_vacuum():
    # |up dn> = c^dag_up c^dag_dn |0>, so the amplitude is +1 (not -1)
    vac = np.zeros(4)
    vac[0] = 1.0
    state = cd(0, 2) @ (cd(1, 2) @ vac)
    assert state[3] == pytest.approx(1.0)


def test_bond_hamiltonian_is_hermitian_and_parity_conserving():
    M = 4
    parity = [0, 1, 1, 0]
    h = np.zeros((1 << M, 1 << M))
    for s in range(2):
        m1, m2 = s, 2 + s
        h = h - (cd(m1, M) @ c(m2, M) + cd(m2, M) @ c(m1, M))
    assert np.allclose(h, h.conj().T)
    t = tenes_simple.bond_matrix(h, 2)
    for i1, i2, o1, o2 in np.ndindex(t.shape):
        if (parity[i1] ^ parity[i2]) != (parity[o1] ^ parity[o2]):
            assert t[i1, i2, o1, o2] == 0.0


def test_zero_parameters_give_a_zero_bond_matrix():
    M = 2
    h = 0.0 * (cd(0, M) @ c(1, M))
    assert np.allclose(tenes_simple.bond_matrix(h, 1), 0.0)


def test_onesite_matrix_transposes_to_in_out():
    cdag_up = tenes_simple.fock_cop(True, 0, 2)  # mat[out, in]
    m = tenes_simple.onesite_matrix(cdag_up)  # op[in, out]
    assert cdag_up[1, 0] == pytest.approx(1.0)  # <up| c^dag_up |0> = 1
    assert m[0, 1] == pytest.approx(1.0)


def test_onesite_number_operators():
    n_up = tenes_simple.onesite_matrix(cd(0, 2) @ c(0, 2))
    n_dn = tenes_simple.onesite_matrix(cd(1, 2) @ c(1, 2))
    assert np.allclose(np.diag(n_up), [0.0, 1.0, 0.0, 1.0])
    assert np.allclose(np.diag(n_dn), [0.0, 0.0, 1.0, 1.0])
