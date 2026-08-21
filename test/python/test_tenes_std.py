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
import re
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

    A 2x2 unit cell -- the smallest shape fermion mode allows under C2
    (task-11-contract.md: both L_sub entries must be >= 2) -- with one
    shared site definition broadcast to all four positions via
    index = [], and a single nearest-neighbour bond term (dx, dy) = (1, 0),
    which has make_path length 1 on this lattice (unchanged from the
    former 1x1 fixture: the taxicab hop count between two lattice sites
    does not depend on the unit-cell tiling).
    """
    return {
        "parameter": {"general": {"fermion": True}},
        "tensor": {
            "l_sub": [2, 2],
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
    """A four-site (2x2, the minimum fermion-valid shape per C2) fermion-mode
    input, with only a one-site Hamiltonian term (no bonds), so it isolates
    the "every unitcell needs parity" check (C3a) from the bond-distance
    check (C3c).

    `parities` maps a subset of {0, 1} to a parity list for sites 0 and 1;
    a site not present in the mapping is emitted without a `parity` key at
    all, exactly as before. Sites 2 and 3 exist only to satisfy the
    L_sub >= [2, 2] requirement and always carry a valid parity, so they
    never trigger this check themselves.
    """
    unitcell = []
    for index in range(4):
        site = {"index": [index], "physical_dim": 2, "virtual_dim": 2}
        if index in (0, 1):
            if index in parities:
                site["parity"] = parities[index]
        else:
            site["parity"] = [0, 1]
        unitcell.append(site)
    return {
        "parameter": {"general": {"fermion": True}},
        "tensor": {"l_sub": [2, 2], "unitcell": unitcell},
        "hamiltonian": [
            {
                "dim": [2],
                "sites": [],
                "elements": "0 0 1.0 0.0\n1 1 -1.0 0.0",
            }
        ],
    }


def fermion_input_with_missing_parity(missing_index, num_sites=4):
    """A fermion-mode input with `num_sites` unitcell sites (L_sub =
    [2, num_sites // 2], the smallest fermion-valid (>= 2 in both entries,
    per C2) shape that holds `num_sites` sites), all carrying
    `parity = [0, 1]` except `missing_index`, which has none. Only a
    one-site Hamiltonian term is used, so this isolates the "every unitcell
    needs parity" check (C3a).

    `num_sites` defaults to 4 so that the offending index can be chosen
    from {2, 3}: values that cannot collide with any digit already present
    in the boilerplate part of the rejection message (which happens to
    talk about a "0/1" parity entry), keeping the content-pinning checks
    in TestFermionErrorMessageQuality unambiguous.
    """
    assert num_sites % 2 == 0 and num_sites // 2 >= 2
    unitcell = []
    for index in range(num_sites):
        site = {"index": [index], "physical_dim": 2, "virtual_dim": 2}
        if index != missing_index:
            site["parity"] = [0, 1]
        unitcell.append(site)
    return {
        "parameter": {"general": {"fermion": True}},
        "tensor": {"l_sub": [2, num_sites // 2], "unitcell": unitcell},
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
        # Strengthened per code review: it is not enough for the two
        # messages to differ by *something* (a mutant that appended
        # id(self) to an otherwise-generic string would still pass that).
        # Each message must contain ITS OWN offending site index. Sites 2
        # and 3 are used (out of a 4-site unitcell) because neither digit
        # can appear by coincidence in the message's fixed boilerplate
        # (which mentions a "0/1" parity entry, i.e. only the digits 0
        # and 1) -- so a hard-coded message naming the wrong site is
        # guaranteed to fail the corresponding assertion below.
        param_site2 = fermion_input_with_missing_parity(2, num_sites=4)
        param_site3 = fermion_input_with_missing_parity(3, num_sites=4)

        with pytest.raises(RuntimeError) as e2:
            tenes_std.Model(param_site2)
        with pytest.raises(RuntimeError) as e3:
            tenes_std.Model(param_site3)

        msg2, msg3 = str(e2.value), str(e3.value)
        assert re.search(r"\b2\b", msg2), msg2
        assert re.search(r"\b3\b", msg3), msg3
        # Discrimination: a message that hard-codes one site index cannot
        # simultaneously name the other offending site.
        assert not re.search(r"\b3\b", msg2), msg2
        assert not re.search(r"\b2\b", msg3), msg3
        assert msg2 != msg3

    def test_multihop_bond_message_differs_by_offending_bond(self):
        # Strengthened per code review, same rationale as above: pin the
        # actual bond identifiers (source_site, dx, and the resulting hop
        # count) rather than merely observing that two messages differ.
        # dx = 5 and dx = 6 are used (rather than small values like 1-3)
        # so neither digit can coincide with source_site (0) or dy (0) in
        # the message.
        param_5hop = copy.deepcopy(minimal_fermion_std_input())
        param_5hop["hamiltonian"][0]["bonds"] = "0 5 0\n"  # make_path length 5

        param_6hop = copy.deepcopy(minimal_fermion_std_input())
        param_6hop["hamiltonian"][0]["bonds"] = "0 6 0\n"  # make_path length 6

        with pytest.raises(RuntimeError) as e5:
            tenes_std.Model(param_5hop)
        with pytest.raises(RuntimeError) as e6:
            tenes_std.Model(param_6hop)

        msg5, msg6 = str(e5.value), str(e6.value)
        assert re.search(r"\b5\b", msg5), msg5
        assert re.search(r"\b6\b", msg6), msg6
        # Discrimination: a message that hard-codes one bond's
        # displacement/hop-count cannot simultaneously name the other.
        assert not re.search(r"\b6\b", msg5), msg5
        assert not re.search(r"\b5\b", msg6), msg6
        assert msg5 != msg6

    def test_ops_form_observable_message_differs_by_offending_observable(self):
        # Strengthened per code review: pin the observable's own name in
        # its own message. The names are long and distinctive on purpose,
        # so an accidental substring collision with unrelated message text
        # is not a realistic concern, and a message that hard-codes one
        # observable's name is guaranteed to fail for the other.
        param_a = copy.deepcopy(minimal_fermion_std_input())
        param_a["observable"] = {
            "twosite": [
                {
                    "name": "alpha_observable_marker",
                    "group": 1,
                    "bonds": "0 1 0\n",
                    "ops": [0, 1],
                }
            ]
        }

        param_b = copy.deepcopy(minimal_fermion_std_input())
        param_b["observable"] = {
            "twosite": [
                {
                    "name": "beta_observable_marker",
                    "group": 2,
                    "bonds": "0 1 0\n",
                    "ops": [1, 0],
                }
            ]
        }

        with pytest.raises(RuntimeError) as ea:
            tenes_std.Model(param_a)
        with pytest.raises(RuntimeError) as eb:
            tenes_std.Model(param_b)

        msg_a, msg_b = str(ea.value), str(eb.value)
        assert "alpha_observable_marker" in msg_a
        assert "beta_observable_marker" in msg_b
        # Discrimination: a message that hard-codes one observable's name
        # cannot simultaneously name the other.
        assert "beta_observable_marker" not in msg_a
        assert "alpha_observable_marker" not in msg_b
        assert msg_a != msg_b


# ---------------------------------------------------------------------------
# C2 (task-11-contract.md): fermion mode requires BOTH tensor unit-cell
# dimensions to be >= 2.
#
# L_sub = [2, 1] (or [1, 2]) has one entry equal to 1, so with skew = 0 a
# site becomes its own neighbour along that direction; the C++ constructor
# computes skew % LX = 0 in that case, so the pre-existing skew guard can
# never fire. LatticeGraph happily builds a shortest-path graph for this
# geometry -- no crash, no warning -- so it must be rejected up front here.
# ---------------------------------------------------------------------------


class TestFermionUnitCellDimensionGuard:
    def test_l_sub_2_1_is_rejected(self):
        param = copy.deepcopy(minimal_fermion_std_input())
        param["tensor"]["l_sub"] = [2, 1]
        with pytest.raises(RuntimeError) as excinfo:
            tenes_std.Model(param)
        message = str(excinfo.value)
        # Identifying content: name the unit-cell dimensions, not just
        # "raises something".
        assert re.search(r"\bL_sub\b", message, re.I) or re.search(
            r"\bdimension", message, re.I
        )
        assert re.search(r"\b2\b", message)
        assert re.search(r"\b1\b", message)

    def test_l_sub_1_2_is_rejected(self):
        param = copy.deepcopy(minimal_fermion_std_input())
        param["tensor"]["l_sub"] = [1, 2]
        with pytest.raises(RuntimeError) as excinfo:
            tenes_std.Model(param)
        message = str(excinfo.value)
        assert re.search(r"\bL_sub\b", message, re.I) or re.search(
            r"\bdimension", message, re.I
        )
        assert re.search(r"\b1\b", message)
        assert re.search(r"\b2\b", message)

    def test_l_sub_2_1_message_does_not_mention_the_internal_milestone(self):
        param = copy.deepcopy(minimal_fermion_std_input())
        param["tensor"]["l_sub"] = [2, 1]
        with pytest.raises(RuntimeError) as excinfo:
            tenes_std.Model(param)
        message = str(excinfo.value)
        assert "M1" not in message
        assert "M2" not in message

    def test_l_sub_2_2_is_accepted(self):
        # Regression net (i): the exact same fermionic configuration
        # (parity everywhere, a single nearest-neighbour bond) as the
        # rejected cases above, but with both dimensions >= 2, must not
        # raise -- this is minimal_fermion_std_input() itself.
        model = tenes_std.Model(minimal_fermion_std_input())
        assert model.unitcell.L == [2, 2]

    def test_bosonic_l_sub_2_1_still_constructs(self):
        # Regression net (ii): the bug -- and this new guard -- is
        # fermion-specific. A bosonic (no `fermion` key) 2x1 cell, the
        # shape the benchmark harness uses, must be entirely unaffected.
        param = {
            "tensor": {
                "l_sub": [2, 1],
                "unitcell": [
                    {"index": [], "physical_dim": 2, "virtual_dim": 2},
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
        model = tenes_std.Model(param)  # must not raise
        assert model.unitcell.L == [2, 1]


# ---------------------------------------------------------------------------
# C2/C3 (task-4b-contract.md): fermion = true with [tensor] skew != 0 must
# be rejected.
#
# work/skew-validation/FINDINGS.md measured this at the tenes_simple layer
# (a 20.6% energy shift, density drifted off half filling for the standard
# skew = 1 two-site cell); tenes_std must refuse the same combination when
# `skew` arrives as an explicit std.toml field, before evolution operators
# are built. Where the parsed value lives: `Model.unitcell.skew`
# (`Unitcell.load_dict` sets it from `tensor.skew`, defaulting to 0) --
# tests read it from there rather than re-parsing the raw dict.
# ---------------------------------------------------------------------------


class TestFermionSkewGuard:
    def test_fermion_mode_with_nonzero_skew_is_rejected(self):
        param = copy.deepcopy(minimal_fermion_std_input())
        param["tensor"]["skew"] = 1
        with pytest.raises(RuntimeError) as excinfo:
            tenes_std.Model(param)
        message = str(excinfo.value)
        assert re.search(r"\bskew\b", message, re.I)
        assert re.search(r"\b1\b", message)

    def test_fermion_mode_with_nonzero_skew_message_is_fermion_specific(self):
        # C2 requirements mirror C1: the message must read as a measured
        # correctness bug (wrong numbers), not as a generic missing-scope
        # complaint -- which in this file's own established phrasing reads
        # "Fermion mode does not support ...". The skew guard must not
        # reuse that exact framing.
        param = copy.deepcopy(minimal_fermion_std_input())
        param["tensor"]["skew"] = 1
        with pytest.raises(RuntimeError) as excinfo:
            tenes_std.Model(param)
        message = str(excinfo.value)
        assert re.search(r"\bfermion", message, re.I)
        assert re.search(r"measured|wrong|incorrect|known limitation", message, re.I)
        assert "does not support" not in message

    def test_skew_guard_message_does_not_mention_the_internal_milestone(self):
        param = copy.deepcopy(minimal_fermion_std_input())
        param["tensor"]["skew"] = 1
        with pytest.raises(RuntimeError) as excinfo:
            tenes_std.Model(param)
        message = str(excinfo.value)
        assert "M1" not in message
        assert "M2" not in message

    def test_skew_guard_message_names_the_offending_skew_value(self):
        # Identifying content, not merely "the two messages differ" (a
        # content-free string differing only by e.g. id(self) would still
        # pass a bare inequality check): each message must contain its OWN
        # offending skew value. 7 and 9 are used so neither digit can
        # coincide with any digit already in the fixed boilerplate.
        param_a = copy.deepcopy(minimal_fermion_std_input())
        param_a["tensor"]["skew"] = 7
        param_b = copy.deepcopy(minimal_fermion_std_input())
        param_b["tensor"]["skew"] = 9

        with pytest.raises(RuntimeError) as ea:
            tenes_std.Model(param_a)
        with pytest.raises(RuntimeError) as eb:
            tenes_std.Model(param_b)

        msg_a, msg_b = str(ea.value), str(eb.value)
        assert re.search(r"\b7\b", msg_a), msg_a
        assert re.search(r"\b9\b", msg_b), msg_b
        # Discrimination: a message that hard-codes one skew value cannot
        # simultaneously name the other.
        assert not re.search(r"\b9\b", msg_a), msg_a
        assert not re.search(r"\b7\b", msg_b), msg_b
        assert msg_a != msg_b

    def test_fermion_mode_with_default_skew_is_still_accepted(self):
        # Regression net (i): fermion = true with no explicit `skew` key
        # (skew defaults to 0 in Unitcell.load_dict) must keep working
        # exactly as today -- this is minimal_fermion_std_input() itself.
        model = tenes_std.Model(minimal_fermion_std_input())
        assert model.unitcell.skew == 0

    def test_fermion_mode_with_explicit_zero_skew_is_accepted(self):
        # Regression net (i), the explicit form: skew = 0 written out
        # must not be mistaken for "skew present" by the guard.
        param = copy.deepcopy(minimal_fermion_std_input())
        param["tensor"]["skew"] = 0
        model = tenes_std.Model(param)
        assert model.unitcell.skew == 0

    def test_bosonic_input_with_nonzero_skew_is_still_accepted(self):
        # Regression net (ii): the bug is fermion-specific. A bosonic
        # input (no `fermion` key at all) with skew != 0 must be entirely
        # untouched by this guard.
        param = minimal_std_input()
        param["tensor"]["skew"] = 1
        model = tenes_std.Model(param)
        assert model.unitcell.skew == 1
