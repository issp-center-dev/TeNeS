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

"""Self-checks for the mean-field (open-leg) extension of the Fock oracle."""

import os
import sys

import numpy as np
import pytest

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fermion"),
)

import fock_oracle  # noqa: E402

FT = [False, True]
LAM = [1.0, 0.7]
TOL = 1.0e-12

# R2 reference values printed by ``print_case("horizontal_2site", 2, 1, FT)``.
R2_NORM = 5.46520500622191588e-04
R2_N = 5.35874068543760670e-02


def observable_vector(oracle):
    return np.array(
        [
            oracle.norm(),
            oracle.one_body(0, 0),
            oracle.one_body(1, 1),
            oracle.one_body(0, 1),
            oracle.pairing(0, 1),
            oracle.density_density(0, 1),
        ]
    )


def parse_printed(capsys, name, lx, ly, parity, seed, lam):
    fock_oracle.print_mf_case(name, lx, ly, parity, seed, lam)
    out = {}
    for line in capsys.readouterr().out.splitlines():
        key, value = line.split(" = ")
        out[key] = float(value)
    return out


# 1. Omitting dangling_labels reproduces the existing R2 reference values.
def test_default_oracle_reproduces_r2_reference():
    patch, tensors, leg_parities = fock_oracle.make_case(2, 1, FT)
    oracle = fock_oracle.Oracle(patch, tensors, leg_parities)
    norm, n, bonds = oracle.observables()
    assert abs(norm - R2_NORM) < TOL
    assert abs(n[0] - R2_N) < TOL
    assert abs(n[1] - R2_N) < TOL
    assert len(bonds) == 1
    assert oracle.dangling_labels is None
    assert oracle.nmode == 4


# 2. All open legs one-dimensional even: mf_sum has one assignment and
#    agrees with observables().
def test_mf_sum_single_assignment_matches_observables():
    patch, tensors, leg_parities = fock_oracle.make_case(2, 1, FT)
    norm, n, bonds = fock_oracle.Oracle(patch, tensors, leg_parities).observables()
    a, b, hop, pair = bonds[0]

    calls = []

    def fn(oracle):
        calls.append(oracle)
        return observable_vector(oracle)

    total = fock_oracle.mf_sum(patch, tensors, leg_parities, fn)
    assert len(calls) == 1
    assert calls[0].dangling_labels == {key: 0 for key in patch.open_legs()}
    assert abs(total[0] - norm) < TOL
    assert abs(total[1] / total[0] - n[0]) < TOL
    assert abs(total[2] / total[0] - n[1]) < TOL
    mf_hop = fock_oracle.mf_sum(
        patch,
        tensors,
        leg_parities,
        lambda o: o.one_body(a, b) + o.one_body(b, a),
    )
    mf_pair = fock_oracle.mf_sum(
        patch,
        tensors,
        leg_parities,
        lambda o: o.pairing(a, b) + o.pairing(b, a),
    )
    assert abs(mf_hop / total[0] - hop) < TOL
    assert abs(mf_pair / total[0] - pair) < TOL


# 4. Reversing the creation order of odd open-leg auxiliary modes leaves every
#    mf_sum observable unchanged.
def test_reversed_dangling_creation_order_is_invariant(monkeypatch):
    patch, tensors, leg_parities = fock_oracle.make_case_open(2, 1, FT)
    original = fock_oracle.Oracle.dangling_creation_order
    forward = fock_oracle.mf_sum(patch, tensors, leg_parities, observable_vector)

    # An assignment with two odd open legs on site 0: reversing two creations
    # must flip the sign of the raw state, which proves the hook is in use.
    labels = {key: 0 for key in patch.open_legs()}
    labels[(0, "l")] = 1
    labels[(0, "t")] = 1
    state_forward = fock_oracle.Oracle(
        patch, tensors, leg_parities, dangling_labels=labels
    ).state()
    assert np.linalg.norm(state_forward) > 0.0

    monkeypatch.setattr(
        fock_oracle.Oracle,
        "dangling_creation_order",
        lambda self: list(reversed(original(self))),
    )
    reversed_oracle = fock_oracle.Oracle(
        patch, tensors, leg_parities, dangling_labels=labels
    )
    assert reversed_oracle.dangling_creation_order() == [(0, "t"), (0, "l")]
    assert np.allclose(reversed_oracle.state(), -state_forward, rtol=0.0, atol=TOL)

    backward = fock_oracle.mf_sum(patch, tensors, leg_parities, observable_vector)
    assert np.allclose(forward, backward, rtol=0.0, atol=TOL)
    assert forward[0] > 0.0


# 5. The seed-0 MF datasets have a detectable hopping amplitude.
@pytest.mark.parametrize(
    "name,lx,ly",
    [("mf_horizontal_2site", 2, 1), ("mf_vertical_2site", 1, 2)],
)
def test_mf_datasets_have_nonzero_hopping(capsys, name, lx, ly):
    printed = parse_printed(capsys, name, lx, ly, FT, 0, LAM)
    assert set(printed) == {
        f"{name}.norm",
        f"{name}.n0",
        f"{name}.n1",
        f"{name}.hop01",
        f"{name}.pair01",
        f"{name}.nn01",
    }
    norm = printed[f"{name}.norm"]
    assert norm > 0.0
    # hop01 is printed already divided by norm.
    assert abs(printed[f"{name}.hop01"]) > 1.0e-3
    assert 0.0 <= printed[f"{name}.n0"] <= 1.0
    assert 0.0 <= printed[f"{name}.n1"] <= 1.0

    # The printed values are the mf_sum ratios of the contract, with lam applied
    # along every open leg of every site.
    patch, tensors, leg_parities = fock_oracle.make_case_open(lx, ly, FT, 0)
    lam = np.asarray(LAM)
    weighted = [t.copy() for t in tensors]
    for site, leg in patch.open_legs():
        shape = [1] * 5
        shape[fock_oracle.LEGS.index(leg)] = 2
        weighted[site] = weighted[site] * lam.reshape(shape)
    total = fock_oracle.mf_sum(patch, weighted, leg_parities, observable_vector)
    hop = fock_oracle.mf_sum(
        patch, weighted, leg_parities, lambda o: o.one_body(0, 1) + o.one_body(1, 0)
    )
    pair = fock_oracle.mf_sum(
        patch, weighted, leg_parities, lambda o: o.pairing(0, 1) + o.pairing(1, 0)
    )
    assert abs(printed[f"{name}.norm"] - total[0]) < TOL
    assert abs(printed[f"{name}.n0"] - total[1] / total[0]) < TOL
    assert abs(printed[f"{name}.n1"] - total[2] / total[0]) < TOL
    assert abs(printed[f"{name}.hop01"] - hop / total[0]) < TOL
    assert abs(printed[f"{name}.pair01"] - pair / total[0]) < TOL
    assert abs(printed[f"{name}.nn01"] - total[5] / total[0]) < TOL


# 3. Boson limit on a single site (virtual [False, False]): every leg is open,
#    so the label filter of the fixed open legs must reproduce the plain numpy
#    sum even when two labels share the same (even) parity.  Internal bonds are
#    deliberately excluded: the Fock encoding cannot tell two equal-parity
#    labels apart on a bond (see the Oracle docstring).
def test_boson_limit_single_site_all_even():
    patch, tensors, leg_parities = fock_oracle.make_case_open(1, 1, [False, False])
    got = fock_oracle.mf_sum(patch, tensors, leg_parities, lambda o: o.norm())
    expected = float(np.sum(tensors[0] ** 2))
    assert expected > 0.0
    assert abs(got - expected) < TOL


# Supplementary: the single-site MF reference values (lam applied along each
# open leg) against plain numpy.
@pytest.mark.parametrize("name,seed", [("mf_single_site", 0), ("mf_seed173", 173)])
def test_single_site_mf_case_matches_plain_numpy(capsys, name, seed):
    patch, tensors, leg_parities = fock_oracle.make_case_open(1, 1, FT, seed)
    lam = np.asarray(LAM)
    weighted = tensors[0]
    for axis in range(4):
        shape = [1] * 5
        shape[axis] = 2
        weighted = weighted * lam.reshape(shape)
    expected_norm = float(np.sum(weighted**2))
    expected_n0 = float(np.sum(weighted[..., 1] ** 2)) / expected_norm
    printed = parse_printed(capsys, name, 1, 1, FT, seed, LAM)
    assert set(printed) == {f"{name}.norm", f"{name}.n0"}
    assert abs(printed[f"{name}.norm"] - expected_norm) < TOL
    assert abs(printed[f"{name}.n0"] - expected_n0) < TOL


# Supplementary: density_density against the Fock amplitudes directly.
def test_density_density_matches_fock_amplitudes():
    patch, tensors, leg_parities = fock_oracle.make_case_open(2, 1, FT)
    labels = {key: 0 for key in patch.open_legs()}
    labels[(0, "l")] = 1
    labels[(1, "b")] = 1
    oracle = fock_oracle.Oracle(patch, tensors, leg_parities, dangling_labels=labels)
    psi = oracle.physical_state()
    both = sum(amp * amp for state, amp in enumerate(psi) if (state & 0b11) == 0b11)
    assert both > 0.0
    assert abs(oracle.density_density(0, 1) - both) < TOL
    assert abs(oracle.density_density(1, 0) - both) < TOL
    for i in range(2):
        assert abs(oracle.density_density(i, i) - oracle.one_body(i, i)) < TOL


# Extension 1: one auxiliary mode per open leg, numbered after the physical and
# internal-bond modes in sorted (site, leg_name) order.
def test_dangling_modes_follow_sorted_open_legs():
    patch, tensors, leg_parities = fock_oracle.make_case_open(2, 1, FT)
    open_legs = [(0, "b"), (0, "l"), (0, "t"), (1, "b"), (1, "r"), (1, "t")]
    assert patch.open_legs() == open_legs
    labels = {key: 0 for key in open_legs}
    oracle = fock_oracle.Oracle(patch, tensors, leg_parities, dangling_labels=labels)
    assert oracle.nmode == 2 + 2 + 6
    assert [oracle.mode[("a",) + key] for key in open_legs] == list(range(4, 10))
    assert oracle.dangling_creation_order() == []
    labels[(1, "t")] = 1
    labels[(0, "b")] = 1
    oracle = fock_oracle.Oracle(patch, tensors, leg_parities, dangling_labels=labels)
    assert oracle.dangling_creation_order() == [(0, "b"), (1, "t")]


# 6. dangling_labels must cover every open leg.
def test_missing_dangling_label_raises():
    patch, tensors, leg_parities = fock_oracle.make_case_open(2, 1, FT)
    labels = {key: 0 for key in patch.open_legs()}
    assert len(labels) == 6
    fock_oracle.Oracle(patch, tensors, leg_parities, dangling_labels=labels)
    incomplete = dict(labels)
    del incomplete[(1, "r")]
    with pytest.raises(ValueError):
        fock_oracle.Oracle(patch, tensors, leg_parities, dangling_labels=incomplete)
