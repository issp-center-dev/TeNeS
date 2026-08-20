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

import copy
import io
import os
import sys

import pytest
import toml

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tool")
)

import numpy as np

import tenes_std


def minimal_std_input():
    """A minimal valid standard-mode input without a [parameter] section."""
    bonds = "\n".join(
        "{} {} {}".format(source, dx, dy)
        for source in range(4)
        for dx, dy in ((1, 0), (0, 1))
    )
    return {
        "tensor": {
            "l_sub": [2, 2],
            "unitcell": [{"index": [], "physical_dim": 2, "virtual_dim": 2}],
        },
        "hamiltonian": [
            {
                "dim": [2, 2],
                "bonds": bonds,
                "elements": "0 0 0 0 1.0 0.0\n1 1 1 1 -1.0 0.0",
            }
        ],
    }


class TestParseBond:
    def test_valid_line(self):
        bond = tenes_std.parse_bond("0 1 -1")
        assert bond.source_site == 0
        assert bond.dx == 1
        assert bond.dy == -1

    def test_comment_line_returns_none(self):
        assert tenes_std.parse_bond("# comment") is None
        assert tenes_std.parse_bond("") is None


class TestIsHermite:
    def test_hermitian(self):
        A = np.array([[1.0, 1.0j], [-1.0j, 2.0]])
        assert tenes_std.is_hermite(A)

    def test_not_hermitian(self):
        A = np.array([[0.0, 1.0], [0.0, 0.0]])
        assert not tenes_std.is_hermite(A)

    def test_tolerates_rounding_error(self):
        A = np.array([[1.0, 0.5], [0.5 + 1e-16, 2.0]])
        assert tenes_std.is_hermite(A)

    def test_nonhermitian_hamiltonian_raises(self):
        param = minimal_std_input()
        param["hamiltonian"][0]["elements"] = "0 0 1 1 1.0 0.0"
        with pytest.raises(RuntimeError):
            tenes_std.Model(param)


class TestModel:
    def test_missing_parameter_section_is_allowed(self):
        model = tenes_std.Model(minimal_std_input())
        assert model.simple_tau == [0.01]
        assert model.full_tau == [0.01]

    def test_missing_tensor_section_raises(self):
        param = minimal_std_input()
        del param["tensor"]
        with pytest.raises(RuntimeError):
            tenes_std.Model(param)

    def test_missing_hamiltonian_section_raises(self):
        param = minimal_std_input()
        del param["hamiltonian"]
        with pytest.raises(RuntimeError):
            tenes_std.Model(param)


class TestMergeInputDict:
    def test_known_subsections_are_merged(self):
        d1 = {"parameter": {"general": {"is_real": True}}}
        d2 = {"parameter": {"simple_update": {"num_step": 100}}}
        tenes_std.merge_input_dict(d1, d2)
        assert d1["parameter"]["general"]["is_real"] is True
        assert d1["parameter"]["simple_update"]["num_step"] == 100

    def test_unknown_subsections_are_kept(self):
        d1 = {"parameter": {"general": {"is_real": True}}}
        d2 = {"parameter": {"tensor": {"save_dir": "ckpt"}}}
        tenes_std.merge_input_dict(d1, d2)
        assert d1["parameter"]["tensor"]["save_dir"] == "ckpt"

    def test_conflicting_keys_raise(self):
        d1 = {"parameter": {"general": {"is_real": True}}}
        d2 = {"parameter": {"general": {"is_real": False}}}
        with pytest.raises(RuntimeError):
            tenes_std.merge_input_dict(d1, d2)


class TestUnitcell:
    def test_valid_unitcell(self):
        unitcell = tenes_std.Unitcell(
            {
                "l_sub": [2, 1],
                "unitcell": [
                    {"index": [], "physical_dim": 2, "virtual_dim": 2},
                ],
            }
        )
        assert unitcell.numsites() == 2

    def test_missing_site_raises_runtime_error(self):
        with pytest.raises(RuntimeError):
            tenes_std.Unitcell(
                {
                    "l_sub": [2, 1],
                    "unitcell": [
                        {"index": [0], "physical_dim": 2, "virtual_dim": 2},
                    ],
                }
            )


# ---------------------------------------------------------------------------
# Fermion-mode fixtures (task-5-contract.md, clauses C1-C5)
#
# `minimal_fermion_std_input` is a single-site, nearest-neighbour-only,
# `fermion = true` std.toml-shaped dict: the smallest input that a correct
# implementation must ACCEPT. Individual tests take a deep copy and change
# exactly one thing to trigger exactly one rejection.
# ---------------------------------------------------------------------------


def minimal_fermion_std_input():
    """A minimal valid fermion-mode standard-mode input.

    One site, L_sub = [1, 1], a single nearest-neighbour bond term
    (dx, dy) = (1, 0), which has make_path length 1 on this lattice.
    """
    return {
        "parameter": {"general": {"fermion": True}},
        "tensor": {
            "l_sub": [1, 1],
            "unitcell": [
                {
                    "index": [],
                    "physical_dim": 2,
                    "virtual_dim": 2,
                    "parity": [0, 1],
                },
            ],
        },
        "hamiltonian": [
            {
                "dim": [2, 2],
                "bonds": "0 1 0\n",
                "elements": "0 0 0 0 1.0 0.0\n1 1 1 1 -1.0 0.0",
            }
        ],
    }


def two_site_fermion_input(parities):
    """A two-site fermion-mode input, with only a one-site Hamiltonian term
    (no bonds), so it isolates the "every unitcell needs parity" check (C3a)
    from the bond-distance check (C3c).

    `parities` maps a subset of {0, 1} to a parity list; a site not present
    in the mapping is emitted without a `parity` key at all.
    """
    unitcell = []
    for index in (0, 1):
        site = {"index": [index], "physical_dim": 2, "virtual_dim": 2}
        if index in parities:
            site["parity"] = parities[index]
        unitcell.append(site)
    return {
        "parameter": {"general": {"fermion": True}},
        "tensor": {"l_sub": [2, 1], "unitcell": unitcell},
        "hamiltonian": [
            {
                "dim": [2],
                "sites": [],
                "elements": "0 0 1.0 0.0\n1 1 -1.0 0.0",
            }
        ],
    }


class TestParityRoundTrip:
    """C1: parity must survive the round trip."""

    def test_parity_is_available_after_parsing(self):
        # LocalTensor.__init__ currently reads only physical_dim and
        # virtual_dim; parity is dropped at parse time (contract C1).
        lt = tenes_std.LocalTensor(
            {"physical_dim": 2, "virtual_dim": 2, "parity": [0, 1]}
        )
        assert list(lt.parity) == [0, 1]

    def test_parity_round_trips_through_to_toml(self):
        model = tenes_std.Model(minimal_fermion_std_input())
        buf = io.StringIO()
        model.to_toml(buf)
        parsed = toml.loads(buf.getvalue())
        unitcells = parsed["tensor"]["unitcell"]
        assert len(unitcells) == 1
        assert unitcells[0]["parity"] == [0, 1]

    def test_local_tensor_without_parity_parses_without_error(self):
        # A unitcell without parity (every bosonic input) must not gain
        # one, and must not break.
        lt = tenes_std.LocalTensor({"physical_dim": 2, "virtual_dim": 2})
        assert getattr(lt, "parity", None) is None


class TestParityValidation:
    """C2: parity must be validated."""

    def test_wrong_length_parity_raises(self):
        with pytest.raises(RuntimeError):
            tenes_std.LocalTensor(
                {"physical_dim": 2, "virtual_dim": 2, "parity": [0, 1, 0]}
            )

    def test_non_binary_parity_raises(self):
        with pytest.raises(RuntimeError):
            tenes_std.LocalTensor(
                {"physical_dim": 2, "virtual_dim": 2, "parity": [0, 2]}
            )

    def test_valid_parity_is_accepted(self):
        # Sanity check: a correctly-shaped, 0/1-valued parity must not be
        # rejected by the C2 validation.
        tenes_std.LocalTensor({"physical_dim": 2, "virtual_dim": 2, "parity": [0, 1]})


class TestFermionModeValidation:
    """C3: fermion-mode input validation, before gate generation."""

    def test_missing_parity_on_any_unitcell_is_rejected(self):
        # C3a: a unitcell with no parity at all, under fermion = true.
        param = two_site_fermion_input({0: [0, 1]})  # site 1 has none
        with pytest.raises(RuntimeError):
            tenes_std.Model(param)

    def test_ops_form_twosite_observable_is_rejected(self):
        # C3b: fermion mode requires explicit `elements`, not `ops = [i, j]`.
        param = copy.deepcopy(minimal_fermion_std_input())
        param["observable"] = {
            "twosite": [
                {
                    "name": "hopping",
                    "group": 1,
                    "bonds": "0 1 0\n",
                    "ops": [0, 1],
                }
            ]
        }
        with pytest.raises(RuntimeError):
            tenes_std.Model(param)

    def test_multihop_bond_is_rejected(self):
        # C3c: a bond term whose graph.make_path length is not 1 must be
        # rejected before gates are built, because make_evolution_twosite's
        # decomposition places an unsigned identity on the intermediate
        # site (wrong for fermions) while still emitting nearest-neighbour
        # gates that no downstream guard would catch.
        param = copy.deepcopy(minimal_fermion_std_input())
        param["hamiltonian"][0]["bonds"] = "0 2 0\n"  # make_path length 2
        with pytest.raises(RuntimeError):
            tenes_std.Model(param)

    def test_nearest_neighbor_bond_is_accepted_in_same_fermionic_config(self):
        # Pins the DISTINCTION C3c cares about: in the exact same
        # fermionic configuration that rejects a 2-hop bond above, a
        # 1-hop bond must be accepted. A test that only shows rejection
        # would not prove the guard discriminates rather than rejecting
        # everything in fermion mode.
        param = copy.deepcopy(minimal_fermion_std_input())
        param["hamiltonian"][0]["bonds"] = "0 1 0\n"  # make_path length 1
        model = tenes_std.Model(param)  # must not raise
        assert len(model.simple_updates) == 1


class TestBosonicInputsUnaffected:
    """C4: bosonic inputs must be entirely unaffected.

    These fixtures carry no `fermion = true` and no `parity`. Every test
    here would go red if the implementer applied a fermion-only check
    (missing-parity, ops-form, or multihop-bond) unconditionally instead
    of gating it on `parameter.general.fermion`.
    """

    def test_long_distance_bond_still_decomposes(self):
        # Same lattice/bond as TestFermionModeValidation's rejected 2-hop
        # bond (make_path length 2), but with no `fermion` key and no
        # `parity` anywhere: this must still go through
        # make_evolution_twosite's SVD decomposition into a chain of two
        # nearest-neighbour gates, exactly as it does today.
        param = {
            "tensor": {
                "l_sub": [1, 1],
                "unitcell": [{"index": [], "physical_dim": 2, "virtual_dim": 2}],
            },
            "hamiltonian": [
                {
                    "dim": [2, 2],
                    "bonds": "0 2 0\n",
                    "elements": "0 0 0 0 1.0 0.0\n1 1 1 1 -1.0 0.0",
                }
            ],
        }
        model = tenes_std.Model(param)  # must not raise
        assert len(model.simple_updates) == 2

    def test_ops_form_twosite_observable_still_works(self):
        param = minimal_std_input()
        param["observable"] = {
            "twosite": [
                {
                    "name": "sxsx",
                    "group": 1,
                    "bonds": "0 1 0\n",
                    "ops": [0, 0],
                }
            ]
        }
        model = tenes_std.Model(param)
        named = [t for t in model.twobodies if t.name == "sxsx"]
        assert len(named) == 1
        assert named[0].ops == [0, 0]
        assert named[0].elements is None

    def test_missing_parity_does_not_raise_without_fermion_flag(self):
        param = two_site_fermion_input({0: [0, 1]})  # site 1 has none
        del param["parameter"]["general"]["fermion"]
        tenes_std.Model(param)  # must not raise

    def test_output_has_no_parity_field(self):
        model = tenes_std.Model(minimal_std_input())
        buf = io.StringIO()
        model.to_toml(buf)
        assert "parity" not in buf.getvalue()


class TestFermionErrorMessageQuality:
    """C5: error-message quality.

    No message may contain the internal milestone labels "M1"/"M2", and a
    rejection must name the offending input rather than print a single
    canned string regardless of what triggered it.
    """

    def _violations(self):
        missing_parity = two_site_fermion_input({0: [0, 1]})  # site 1 has none

        ops_form = copy.deepcopy(minimal_fermion_std_input())
        ops_form["observable"] = {
            "twosite": [
                {"name": "hopping", "group": 1, "bonds": "0 1 0\n", "ops": [0, 1]}
            ]
        }

        multihop = copy.deepcopy(minimal_fermion_std_input())
        multihop["hamiltonian"][0]["bonds"] = "0 2 0\n"

        bad_parity_length = copy.deepcopy(minimal_fermion_std_input())
        bad_parity_length["tensor"]["unitcell"][0]["parity"] = [0, 1, 0]

        return [missing_parity, ops_form, multihop, bad_parity_length]

    def test_no_milestone_labels_in_rejection_messages(self):
        for param in self._violations():
            with pytest.raises(RuntimeError) as excinfo:
                tenes_std.Model(param)
            message = str(excinfo.value)
            assert "M1" not in message
            assert "M2" not in message

    def test_missing_parity_message_differs_by_offending_site(self):
        param_site1 = two_site_fermion_input({0: [0, 1]})  # site 1 has none
        param_site0 = two_site_fermion_input({1: [0, 1]})  # site 0 has none

        with pytest.raises(RuntimeError) as e1:
            tenes_std.Model(param_site1)
        with pytest.raises(RuntimeError) as e0:
            tenes_std.Model(param_site0)

        assert str(e1.value) != str(e0.value)

    def test_multihop_bond_message_differs_by_offending_bond(self):
        param_2hop = copy.deepcopy(minimal_fermion_std_input())
        param_2hop["hamiltonian"][0]["bonds"] = "0 2 0\n"  # make_path length 2

        param_3hop = copy.deepcopy(minimal_fermion_std_input())
        param_3hop["hamiltonian"][0]["bonds"] = "0 3 0\n"  # make_path length 3

        with pytest.raises(RuntimeError) as e2:
            tenes_std.Model(param_2hop)
        with pytest.raises(RuntimeError) as e3:
            tenes_std.Model(param_3hop)

        assert str(e2.value) != str(e3.value)

    def test_ops_form_observable_message_differs_by_offending_observable(self):
        param_a = copy.deepcopy(minimal_fermion_std_input())
        param_a["observable"] = {
            "twosite": [
                {"name": "alpha", "group": 1, "bonds": "0 1 0\n", "ops": [0, 1]}
            ]
        }

        param_b = copy.deepcopy(minimal_fermion_std_input())
        param_b["observable"] = {
            "twosite": [{"name": "beta", "group": 2, "bonds": "0 0 1\n", "ops": [1, 0]}]
        }

        with pytest.raises(RuntimeError) as ea:
            tenes_std.Model(param_a)
        with pytest.raises(RuntimeError) as eb:
            tenes_std.Model(param_b)

        assert str(ea.value) != str(eb.value)
