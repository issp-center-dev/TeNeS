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

    A 2x2 unit cell -- no site is its own nearest neighbour, so the
    unit-cell guard (section 2.2 of
    docs/superpowers/specs/2026-09-11-fermion-skew-guard-contract.md) lets
    it through -- with one shared site definition broadcast to all four
    positions via index = [], and a single nearest-neighbour bond term
    (dx, dy) = (1, 0), which has make_path length 1 on this lattice
    (unchanged from the former 1x1 fixture: the taxicab hop count between
    two lattice sites does not depend on the unit-cell tiling).
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
    """A four-site (2x2, free of self-neighbour sites) fermion-mode
    input, with only a one-site Hamiltonian term (no bonds), so it isolates
    the "every unitcell needs parity" check (C3a) from the bond-distance
    check (C3c).

    `parities` maps a subset of {0, 1} to a parity list for sites 0 and 1;
    a site not present in the mapping is emitted without a `parity` key at
    all, exactly as before. Sites 2 and 3 exist only to keep the cell free
    of self-neighbour sites and always carry a valid parity, so they
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
    [2, num_sites // 2], a cell with at least two sites in each direction,
    so no site is its own neighbour, that holds `num_sites` sites), all
    carrying `parity = [0, 1]` except `missing_index`, which has none. Only a
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
# docs/superpowers/specs/2026-09-11-fermion-skew-guard-contract.md section
# 2.2: fermion mode refuses a [tensor] cell if and only if some site is its
# own nearest neighbour through the periodic + skew boundary. With T(x, y) =
# T(x + skew, y + L_sub[1]) that is exactly L_sub[0] == 1 (horizontal
# self-neighbour, whatever the skew), or L_sub[1] == 1 with skew = 0 mod
# L_sub[0] (vertical self-neighbour).
#
# This replaces two guards whose premise is gone:
#   * "skew != 0 is refused": based on a 2026-08-20 measurement that is
#     RETRACTED. It predated the CTM folding fix 3bef24a4, compared a
#     one-row cell with a 2x2 cell, and read an ansatz restriction as a sign
#     error: [2, 1] skew 1 is exactly the [2, 2] skew-0 calculation
#     restricted to T0 = T3, T1 = T2, and the simple update stays in that
#     symmetric subspace. The note
#     docs/superpowers/notes/2026-09-11-fermion-skew-revisit.md has the
#     details, and shows a skewed cell and its unfolded skew-0 equivalent
#     agree bit for bit in the simple update and to finite-chi precision in
#     the CTM measurement; test/fermion/skew_unfold.cpp keeps it that way.
#   * "both L_sub entries >= 2": right about the hazard (a one-wide cell can
#     make a site its own neighbour, and the simple update then writes that
#     site twice per bond), but broader than it. L_sub = [2, 1] with skew = 1
#     - what tenes_simple builds for a square lattice with W = 1 - has no
#     self-neighbour: every bond of site 0 goes to site 1.
#
# Where the parsed values live: Model.unitcell.L and Model.unitcell.skew
# (Unitcell.load_dict reads tensor.l_sub and tensor.skew, skew defaulting to
# 0). The skew is the raw input value, possibly negative or >= L_sub[0].
# ---------------------------------------------------------------------------


def fermion_cell_input(l_sub, skew=None):
    """minimal_fermion_std_input() on the given cell, with one horizontal and
    one vertical nearest-neighbour bond from site 0 (the vertical one wraps
    through the skewed boundary on a one-row cell). `skew` None leaves the
    key out."""
    param = copy.deepcopy(minimal_fermion_std_input())
    param["tensor"]["l_sub"] = list(l_sub)
    if skew is not None:
        param["tensor"]["skew"] = skew
    param["hamiltonian"][0]["bonds"] = "0 1 0\n0 0 1\n"
    return param


def names_numbers_after(message, word, numbers):
    """True iff `message` contains `word` followed, after characters that are
    neither digits nor minus signs, by `numbers` in order, each ending at a
    non-digit (case-insensitive). "L_sub = [3, 1]", "L_sub=[3,1]" and
    "l_sub (3 x 1)" all name L_sub 3 1; "skew = -2" names skew -2."""
    pattern = re.escape(word)
    for i, n in enumerate(numbers):
        pattern += ("[^0-9-]*" if i == 0 else "[^0-9-]+") + re.escape(str(n))
    pattern += "(?![0-9])"
    return re.search(pattern, message, re.I) is not None


# (L_sub, skew) cells with a self-neighbour site. LX == 1: [1, 1] and [1, 2]
# with any skew, [1, 3] with 0. LY == 1 with skew = 0 mod LX, including
# non-zero multiples of LX and a negative one.
SELF_NEIGHBOUR_CELLS = [
    ([1, 1], 0),
    ([1, 1], 1),
    ([1, 2], 0),
    ([1, 2], 5),
    ([1, 3], 0),
    ([2, 1], 0),
    ([3, 1], 0),
    ([2, 1], 2),
    ([3, 1], 3),
    ([2, 1], -2),
    ([4, 1], 8),
]

# (L_sub, skew) cells without one, refused before this change: skewed cells
# (skew != 0) and one-row cells whose skew is not a multiple of LX. [2, 1]
# skew 1 is what tenes_simple builds for W = 1; [2, 2] skew 7 and [3, 1]
# skew 4 have |skew| >= LX; [2, 1] skew -1 and [3, 1] skew -1 are negative.
NEWLY_ACCEPTED_CELLS = [
    ([2, 2], 1),
    ([3, 2], 1),
    ([3, 3], 2),
    ([2, 3], 1),
    ([2, 1], 1),
    ([2, 1], -1),
    ([3, 1], 1),
    ([3, 1], 2),
    ([3, 1], -1),
    ([4, 1], 2),
    ([2, 2], 7),
    ([3, 1], 4),
]


def cell_id(cell):
    l_sub, skew = cell
    return "{}x{}-skew{}".format(l_sub[0], l_sub[1], skew)


class TestFermionSelfNeighbourCellGuard:
    @pytest.mark.parametrize(
        "l_sub, skew",
        SELF_NEIGHBOUR_CELLS,
        ids=[cell_id(c) for c in SELF_NEIGHBOUR_CELLS],
    )
    def test_cell_with_a_self_neighbour_site_is_rejected(self, l_sub, skew):
        with pytest.raises(RuntimeError) as excinfo:
            tenes_std.Model(fermion_cell_input(l_sub, skew))
        message = str(excinfo.value)
        # The cell as given: both L_sub values and the skew, each after its
        # own name, so a message that hard-codes another cell fails.
        assert names_numbers_after(message, "L_sub", l_sub), message
        assert names_numbers_after(message, "skew", [skew]), message
        # The cause: a site would be its own nearest neighbour.
        assert re.search(r"\bown\b", message, re.I), message
        assert re.search(r"neighbou?r", message, re.I), message
        assert "M1" not in message
        assert "M2" not in message

    def test_rejection_without_a_skew_key_names_the_default_skew(self):
        # skew is optional and defaults to 0; the message names that value.
        with pytest.raises(RuntimeError) as excinfo:
            tenes_std.Model(fermion_cell_input([3, 1]))
        message = str(excinfo.value)
        assert names_numbers_after(message, "L_sub", [3, 1]), message
        assert names_numbers_after(message, "skew", [0]), message
        assert re.search(r"\bown\b", message, re.I), message

    @pytest.mark.parametrize(
        "l_sub, skew",
        NEWLY_ACCEPTED_CELLS,
        ids=[cell_id(c) for c in NEWLY_ACCEPTED_CELLS],
    )
    def test_cell_without_a_self_neighbour_site_is_accepted(self, l_sub, skew):
        model = tenes_std.Model(fermion_cell_input(l_sub, skew))
        assert model.unitcell.L == l_sub
        assert model.unitcell.skew == skew
        # The emitted input.toml carries the cell and the skew unchanged,
        # and the fermion metadata the solver needs.
        buf = io.StringIO()
        model.to_toml(buf)
        emitted = toml.loads(buf.getvalue())
        assert emitted["tensor"]["L_sub"] == l_sub
        assert emitted["tensor"]["skew"] == skew
        assert emitted["parameter"]["general"]["fermion"] is True
        assert emitted["tensor"]["unitcell"]
        for ucell in emitted["tensor"]["unitcell"]:
            assert ucell["parity"] == [0, 1]
        # Both bonds, the vertical one included, became evolution gates.
        assert len(model.simple_updates) == 2

    def test_unskewed_2x2_cell_is_still_accepted(self):
        # Regression net: minimal_fermion_std_input() itself, no skew key.
        model = tenes_std.Model(minimal_fermion_std_input())
        assert model.unitcell.L == [2, 2]
        assert model.unitcell.skew == 0

    def test_explicit_zero_skew_on_a_2x2_cell_is_accepted(self):
        param = copy.deepcopy(minimal_fermion_std_input())
        param["tensor"]["skew"] = 0
        model = tenes_std.Model(param)
        assert model.unitcell.skew == 0

    def test_other_fermion_guards_still_run_on_a_skewed_cell(self):
        # Lifting the cell-shape refusal must not take the other fermion
        # checks with it: a missing parity on an accepted skewed cell is
        # still refused, for that reason.
        param = fermion_cell_input([2, 1], 1)
        param["tensor"]["unitcell"] = [
            {"index": [0], "physical_dim": 2, "virtual_dim": 2, "parity": [0, 1]},
            {"index": [1], "physical_dim": 2, "virtual_dim": 2},
        ]
        with pytest.raises(RuntimeError, match="parity"):
            tenes_std.Model(param)

    def test_bosonic_one_row_cells_are_untouched(self):
        # Regression net: the guard is fermion-specific. Bosonic (no
        # `fermion` key) cells, self-neighbour ones included ([2, 1] skew 0
        # is the shape the benchmark harness uses), construct as before.
        for l_sub, skew in [([2, 1], 0), ([2, 1], 1), ([2, 2], 1)]:
            param = minimal_std_input()
            param["tensor"]["l_sub"] = l_sub
            param["tensor"]["skew"] = skew
            param["hamiltonian"][0]["bonds"] = "0 1 0\n"
            model = tenes_std.Model(param)  # must not raise
            assert model.unitcell.L == l_sub
            assert model.unitcell.skew == skew
