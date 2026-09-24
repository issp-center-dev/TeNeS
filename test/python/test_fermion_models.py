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

import io
import os
import re
import sys

import numpy as np
import pytest
import toml

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tool")
)

import tenes_simple
import tenes_std


def spinless_param(model_extra=None, lattice_extra=None):
    model = {"type": "spinless fermion", "t": 1.0}
    model.update(model_extra or {})
    lattice = {"type": "square lattice", "L": 2, "W": 2, "virtual_dim": 2}
    lattice.update(lattice_extra or {})
    return {"parameter": {"general": {}}, "lattice": lattice, "model": model}


def std_toml(param):
    text, _ = tenes_simple.tenes_simple(param)
    return toml.loads(text)


def hubbard_param(model_extra=None, lattice_extra=None):
    model = {"type": "hubbard", "t": 1.0}
    model.update(model_extra or {})
    lattice = {"type": "square lattice", "L": 2, "W": 2, "virtual_dim": 2}
    lattice.update(lattice_extra or {})
    return {"parameter": {"general": {}}, "lattice": lattice, "model": model}


class TestSpinlessFermionModel:
    def test_is_selected_by_type(self):
        model = tenes_simple.make_model(spinless_param())
        assert isinstance(model, tenes_simple.SpinlessFermionModel)

    def test_physical_dimension_and_parity(self):
        model = tenes_simple.make_model(spinless_param())
        assert model.N == 2
        assert model.parity == [0, 1]
        assert model.is_fermion is True

    def test_bond_hamiltonian_matches_the_handwritten_sample(self):
        # t = 1, v = 0, mu = 0  ->  sample/07_spinless_fermion/input.toml
        model = tenes_simple.make_model(spinless_param())
        h = model.bondhamiltonian(0, 0, z=4)
        expected = np.zeros((2, 2, 2, 2))
        expected[0, 1, 1, 0] = -1.0
        expected[1, 0, 0, 1] = -1.0
        assert np.allclose(h, expected)

    def test_chemical_potential_is_split_over_the_bonds(self):
        model = tenes_simple.make_model(spinless_param({"mu": 2.0}))
        h = model.bondhamiltonian(0, 0, z=4)
        # -mu/z * (n1 + n2), z = 4  ->  -0.5 on each occupied site
        assert h[1, 0, 1, 0] == pytest.approx(-0.5)
        assert h[0, 1, 0, 1] == pytest.approx(-0.5)
        assert h[1, 1, 1, 1] == pytest.approx(-1.0)

    def test_nearest_neighbour_repulsion(self):
        model = tenes_simple.make_model(spinless_param({"v": 3.0}))
        h = model.bondhamiltonian(0, 0, z=4)
        assert h[1, 1, 1, 1] == pytest.approx(3.0)

    def test_onesite_observable_is_the_density(self):
        model = tenes_simple.make_model(spinless_param())
        assert model.onesite_ops_name == ["n"]
        assert np.allclose(model.onesite_ops[0], np.diag([0.0, 1.0]))

    def test_hopping_is_an_explicit_rank4_observable(self):
        model = tenes_simple.make_model(spinless_param())
        names = [name for name, _ in model.twosite_ops_explicit]
        assert "hopping" in names
        op = dict(model.twosite_ops_explicit)["hopping"]
        assert op.shape == (2, 2, 2, 2)
        assert op[0, 1, 1, 0] == pytest.approx(1.0)
        assert op[1, 0, 0, 1] == pytest.approx(1.0)


class TestSpinlessFermionSchema:
    def test_fermion_flag_is_injected(self):
        assert std_toml(spinless_param())["parameter"]["general"]["fermion"] is True

    def test_explicit_fermion_false_is_rejected(self):
        param = spinless_param()
        param["parameter"]["general"]["fermion"] = False
        with pytest.raises(RuntimeError, match="fermion"):
            tenes_simple.tenes_simple(param)

    def test_parity_is_emitted_for_every_unitcell(self):
        parsed = std_toml(spinless_param())
        for ucell in parsed["tensor"]["unitcell"]:
            assert ucell["parity"] == [0, 1]

    def test_bosonic_models_do_not_emit_parity_or_fermion(self):
        param = {
            "parameter": {"general": {}},
            "lattice": {"type": "square lattice", "L": 2, "W": 2, "virtual_dim": 2},
            "model": {"type": "spin"},
        }
        parsed = std_toml(param)
        assert "fermion" not in parsed["parameter"]["general"]
        for ucell in parsed["tensor"]["unitcell"]:
            assert "parity" not in ucell

    def test_no_twosite_observable_uses_the_ops_form(self):
        text, _ = tenes_simple.tenes_simple(spinless_param())
        assert "ops = " not in text

    def test_vacuum_initial_state(self):
        parsed = std_toml(spinless_param(lattice_extra={"initial": "vacuum"}))
        for ucell in parsed["tensor"]["unitcell"]:
            assert np.allclose(ucell["initial_state"], [1.0, 0.0])

    def test_random_initial_state_stays_scalar(self):
        parsed = std_toml(spinless_param(lattice_extra={"initial": "random"}))
        for ucell in parsed["tensor"]["unitcell"]:
            assert ucell["initial_state"] == [0.0]

    @pytest.mark.parametrize("mode", ["ferro", "antiferro", "full", "cdw"])
    def test_unsupported_initial_states_are_rejected(self, mode):
        with pytest.raises(RuntimeError):
            tenes_simple.tenes_simple(spinless_param(lattice_extra={"initial": mode}))


class TestFermionScopeGuards:
    @pytest.mark.parametrize(
        "latname", ["honeycomb lattice", "triangular lattice", "kagome lattice"]
    )
    def test_non_square_lattices_are_rejected(self, latname):
        param = spinless_param(lattice_extra={"type": latname})
        with pytest.raises(RuntimeError, match="square"):
            tenes_simple.tenes_simple(param)

    def test_square_lattice_is_accepted(self):
        tenes_simple.tenes_simple(spinless_param())

    # In read_params the digit is the BOND TYPE and the number of primes is
    # the NEIGHBOUR LEVEL, so t1 / t2 are still nearest neighbour and only the
    # primed keys go beyond it.
    @pytest.mark.parametrize("key", ["t'", "t''", "v'", "v''"])
    def test_beyond_nearest_neighbour_parameters_are_rejected(self, key):
        param = spinless_param({key: 0.5})
        with pytest.raises(RuntimeError, match="nearest"):
            tenes_simple.tenes_simple(param)

    def test_zero_valued_far_neighbour_parameters_are_accepted(self):
        tenes_simple.tenes_simple(spinless_param({"t'": 0.0}))

    def test_bond_type_variants_of_the_first_neighbour_are_accepted(self):
        # t0 is bond type 0 at the FIRST neighbour level, so the scope guard
        # must not mistake it for a beyond-nearest-neighbour term.
        param = spinless_param()
        param["model"] = {"type": "spinless fermion", "t0": 1.0}
        tenes_simple.tenes_simple(param)

    def test_duplicate_bond_type_specification_is_rejected(self):
        # same rule as BoseHubbardModel.read_params
        param = spinless_param()
        param["model"] = {"type": "spinless fermion", "t": 1.0, "t0": 1.0}
        with pytest.raises(RuntimeError, match="defined twice"):
            tenes_simple.tenes_simple(param)

    def test_correlation_is_rejected(self):
        param = spinless_param()
        param["correlation"] = {"r_max": 5, "operators": [[0, 0]]}
        with pytest.raises(RuntimeError, match="correlation"):
            tenes_simple.tenes_simple(param)

    def test_correlation_length_is_rejected(self):
        param = spinless_param()
        param["correlation_length"] = {"measure": True}
        with pytest.raises(RuntimeError, match="correlation_length"):
            tenes_simple.tenes_simple(param)

    def test_bosonic_models_are_untouched_by_the_guards(self):
        param = {
            "parameter": {"general": {}},
            "lattice": {"type": "kagome lattice", "L": 2, "W": 2, "virtual_dim": 2},
            "model": {"type": "spin", "j": 1.0, "j'": 0.5},
            "correlation": {"r_max": 3, "operators": [[0, 0]]},
        }
        tenes_simple.tenes_simple(param)

    def test_the_message_does_not_mention_the_internal_milestone(self):
        param = spinless_param(lattice_extra={"type": "honeycomb lattice"})
        with pytest.raises(RuntimeError) as excinfo:
            tenes_simple.tenes_simple(param)
        assert "M1" not in str(excinfo.value)


# ---------------------------------------------------------------------------
# docs/superpowers/specs/2026-09-11-fermion-skew-guard-contract.md section
# 2.3: fermionic models on the square lattice with W = 1 are accepted.
#
# SquareLattice realises W = 1 as L_sub = [L, 1] with skew = 1. That cell
# used to be refused on the strength of a 2026-08-20 measurement ("a 20.6%
# energy shift" against a skew = 0 control). That measurement is RETRACTED:
# it predated the CTM folding fix 3bef24a4, compared a one-row cell with a
# 2x2 cell, and read an ansatz restriction as a sign error: the [2, 1]
# skew 1 cell is exactly the [2, 2] skew-0 calculation restricted to
# T0 = T3, T1 = T2, and the simple update stays in that symmetric subspace
# (the bosonic XY ferromagnet does the same; details in
# docs/superpowers/notes/2026-09-11-fermion-skew-revisit.md). The same note
# shows the [2, 1] skew 1 cell and its unfolded [2, 2] skew 0 equivalent
# agree bit for bit in the simple update, to 1e-13 in the CTM measurement
# and to 4e-15 in the full update; test/fermion/skew_unfold.cpp keeps it
# that way. The cell has no site that is its own neighbour (every bond of
# site 0 goes to site 1), so it also passes the unit-cell guard of tenes_std
# and of the solver.
#
# tenes_simple cannot produce a self-neighbour fermion cell at all (the
# square lattice asserts L > 1, and W = 1 always gets skew = 1), so it has no
# cell-shape refusal left. Unchanged: non-square lattices are refused for
# fermionic models (TestFermionScopeGuards), W >= 2 gives skew = 0, and
# bosonic models are untouched.
# ---------------------------------------------------------------------------


def assert_one_row_skew_1_std(text, lattice, L):
    """The std.toml of a W = 1 square lattice: L_sub = [L, 1], skew = 1,
    both in the returned lattice, in the text, and in the parsed tensor
    section."""
    assert lattice.skew == 1
    assert lattice.W == 1
    assert re.search(r"^skew = 1$", text, re.M), text
    tensor = toml.loads(text)["tensor"]
    assert tensor["L_sub"] == [L, 1]
    assert tensor["skew"] == 1


def assert_tenes_std_accepts(text, L, parity):
    """Feed a std.toml to tenes_std: it must build the model and emit an
    input.toml in fermion mode that keeps the cell, the skew and the parity
    metadata of every unitcell."""
    model = tenes_std.Model(toml.loads(text))
    assert model.parameter["general"]["fermion"] is True
    assert model.unitcell.L == [L, 1]
    assert model.unitcell.skew == 1
    assert all(site.parity == parity for site in model.unitcell.sites)
    # One gate per nearest-neighbour bond of the cell: L horizontal and L
    # vertical ones, the vertical ones wrapping through the skewed boundary.
    assert len(model.simple_updates) == 2 * L
    buf = io.StringIO()
    model.to_toml(buf)
    emitted = toml.loads(buf.getvalue())
    assert emitted["parameter"]["general"]["fermion"] is True
    assert emitted["tensor"]["L_sub"] == [L, 1]
    assert emitted["tensor"]["skew"] == 1
    assert emitted["tensor"]["unitcell"]
    for ucell in emitted["tensor"]["unitcell"]:
        assert ucell["parity"] == parity


class TestFermionOneRowSkewedCell:
    @pytest.mark.parametrize("L", [2, 3])
    def test_w1_gives_a_one_row_skew_1_cell(self, L):
        text, lattice = tenes_simple.tenes_simple(
            spinless_param(lattice_extra={"L": L, "W": 1})
        )
        assert_one_row_skew_1_std(text, lattice, L)

    @pytest.mark.parametrize("L", [2, 3])
    def test_w1_std_toml_is_accepted_by_tenes_std(self, L):
        text, _ = tenes_simple.tenes_simple(
            spinless_param(lattice_extra={"L": L, "W": 1})
        )
        assert_tenes_std_accepts(text, L, [0, 1])

    def test_no_square_lattice_cell_is_refused_for_its_shape(self):
        # No cell-shape refusal is left in tenes_simple: every L >= 2 and W
        # it can be given produces a std.toml, with skew = 1 exactly for
        # W = 1.
        for L in (2, 3, 4):
            for W in (1, 2, 3):
                text, lattice = tenes_simple.tenes_simple(
                    spinless_param(lattice_extra={"L": L, "W": W})
                )
                assert lattice.skew == (1 if W == 1 else 0), (L, W)
                assert toml.loads(text)["tensor"]["L_sub"] == [L, W]

    def test_fermionic_model_with_square_cell_no_skew_is_accepted(self):
        # Regression net: fermionic + skew = 0 (the default 2x2 cell).
        text, lattice = tenes_simple.tenes_simple(spinless_param())
        assert lattice.skew == 0

    def test_fermionic_model_with_wide_cell_no_skew_is_accepted(self):
        # Regression net, a second shape: any W != 1 keeps skew = 0.
        param = spinless_param(lattice_extra={"L": 3, "W": 3})
        text, lattice = tenes_simple.tenes_simple(param)
        assert lattice.skew == 0

    def test_bosonic_model_with_skewed_cell_is_still_accepted(self):
        # Regression net: W = 1 is standard practice for spins/bosons and
        # stays untouched.
        param = {
            "parameter": {"general": {}},
            "lattice": {"type": "square lattice", "L": 2, "W": 1, "virtual_dim": 2},
            "model": {"type": "spin", "j": 1.0},
        }
        text, lattice = tenes_simple.tenes_simple(param)
        assert lattice.skew == 1
        assert "parity" not in text
        assert "fermion" not in text


# ---------------------------------------------------------------------------
# task-6-contract.md: the gate must expand back to the bond Hamiltonian.
#
# tenes_simple builds h (through the Fock builder); tenes_std computes the
# imaginary-time gate expm(-tau h). Nothing so far checks the two agree end
# to end through the real pipeline (toml round trip included), only that
# each half separately matches a handwritten/schema expectation. This pins
#
#     (gate - identity) / (-tau) -> h   as tau -> 0
# ---------------------------------------------------------------------------


class TestSpinlessFermionGateReducesToBondHamiltonian:
    def test_gate_reduces_to_bond_hamiltonian_as_tau_to_zero(self):
        # C1: t, V, mu all nonzero, so a dropped term cannot hide.
        param = spinless_param({"t": 1.0, "v": 0.7, "mu": 0.3})
        text, _ = tenes_simple.tenes_simple(param)
        std_param = toml.loads(text)

        # tau small enough that O(tau) truncation error is far below the
        # assert tolerance below.
        tau = 1e-6
        std_param.setdefault("parameter", {})
        std_param["parameter"]["simple_update"] = {"tau": [tau]}

        # C3: real pipeline. model.hamiltonians[0] is parsed back from the
        # generated std.toml's [[hamiltonian]] block (dump_op / load_tensor
        # round trip), not recomputed via model.bondhamiltonian. The
        # uniform square lattice carries the same gate on every bond, so
        # simple_updates[0] is representative and, since it comes from the
        # same source term, corresponds to hamiltonians[0].
        model = tenes_std.Model(std_param)
        ham = model.hamiltonians[0]
        evo = model.simple_updates[0]
        assert isinstance(ham, tenes_std.NNOperator)
        assert isinstance(evo, tenes_std.NNOperator)

        h = ham.elements
        gate = evo.elements
        d = h.shape[0]

        # C2: leg-order-honest comparison. Both h and the gate are
        # op[in1, in2, out1, out2]; build the identity with that same leg
        # placement by explicit construction (not by reshaping through a
        # (d*d, d*d) matrix), so a transposed convention in either tool
        # would show up as a mismatch here rather than being absorbed.
        identity = np.zeros((d, d, d, d))
        for i1 in range(d):
            for i2 in range(d):
                identity[i1, i2, i1, i2] = 1.0

        approx_h = (gate - identity) / (-tau)
        assert np.allclose(approx_h, h, atol=1e-4)


# ---------------------------------------------------------------------------
# task-8-contract.md: the fermionic Hubbard model (d = 4).
#
# At d = 2 (the spinless model above) the reduced-pair loading conventions
# and the intra-site mode order both degenerate: there is only one mode per
# site, and every operator conserves particle number, so a wrong
# intra-site-order or a per-site-only sign rule is invisible. d = 4 is where
# these conventions have visible content, so every matrix element pinned
# below is re-derived independently from the stated conventions:
#
#   * local basis |0>, |up>, |dn>, |up dn>, index i = n_up + 2 n_dn
#   * intra-site creation order FIXED: |up dn> = c^dag_up c^dag_dn |0>
#   * bond modes ordered (site1 up, site1 dn, site2 up, site2 dn)
#   * fock_cop's JW sign is (-1)^{sum of occupied modes strictly below the
#     acted-on mode}, read off the occupation bit string of the *global*
#     two-site state (bit m <-> mode m)
#
# fock_cop itself is not re-derived here (it is pinned by its own
# first-principles tests elsewhere); what is new and untested is the
# HubbardModel class and the lattice/schema wiring around it.
# ---------------------------------------------------------------------------


class TestHubbardModel:
    def test_is_selected_by_type(self):
        model = tenes_simple.make_model(hubbard_param())
        assert isinstance(model, tenes_simple.HubbardModel)

    def test_boson_type_still_selects_the_bose_hubbard_model(self):
        # C1: "hubbard" must not shadow the pre-existing Bose-Hubbard model,
        # which stays reachable through type = "boson".
        param = {
            "parameter": {"general": {}},
            "lattice": {"type": "square lattice", "L": 2, "W": 2, "virtual_dim": 2},
            "model": {"type": "boson", "t": 1.0},
        }
        model = tenes_simple.make_model(param)
        assert isinstance(model, tenes_simple.BoseHubbardModel)
        assert not isinstance(model, tenes_simple.HubbardModel)

    def test_physical_dimension_and_parity(self):
        model = tenes_simple.make_model(hubbard_param())
        assert model.N == 4
        assert model.parity == [0, 1, 1, 0]
        assert model.is_fermion is True

    def test_hopping_carries_the_jordan_wigner_sign(self):
        # Derivation (mode order site1-up=0, site1-dn=1, site2-up=2,
        # site2-dn=3; global occupation bit string g = i1 + 4*i2):
        #
        #   |in1, in2> = |dn, up> -> g = 2 + 4*1 = 6 = 0b0110
        #   c_{2,up} (mode 2) on g=6: mode 2 is occupied, one occupied mode
        #     (mode 1, site1-dn) lies strictly below it -> sign = -1,
        #     result state g=2 = 0b0010 = |dn, 0>
        #   c^dag_{1,up} (mode 0) on g=2: mode 0 is empty, nothing below it
        #     -> sign = +1, result state g=3 = 0b0011 = |up dn, 0>
        #
        # so <up dn, 0| c^dag_{1,up} c_{2,up} |dn, up> = -1, and the -t
        # prefactor in H_bond flips it to +t on h[in1=2, in2=1, out1=3,
        # out2=0].  t = 1.0 is the fixture default.
        model = tenes_simple.make_model(hubbard_param())
        h = model.bondhamiltonian(0, 0, z=4)
        assert h[2, 1, 3, 0] == pytest.approx(1.0)

    def test_zeeman_field_diagonal_element_with_hopping_also_present(self):
        # Second hand-derived element, with h (Zeeman) present and t also
        # nonzero (to prove the diagonal is untouched by hopping): the state
        # |up, up> (i1=i2=1) has Sz_1 = Sz_2 = +1/2, so H contributes
        # -(h/z)(Sz_1+Sz_2) = -(h/z)*1 = -0.5 for h=2, z=4.  U, V, mu are
        # zero, and the -t hopping term cannot map |up, up> back to
        # |up, up>: moving a fermion between two singly-occupied "up" modes
        # always changes at least one site's occupation number, so it never
        # contributes to a diagonal entry (true for any t).
        model = tenes_simple.make_model(hubbard_param({"h": 2.0}))
        h = model.bondhamiltonian(0, 0, z=4)
        assert h[1, 1, 1, 1] == pytest.approx(-0.5)

    def test_combined_u_v_h_diagonal_element(self):
        # Third hand-derived element, mixing U, V and h (with t present but
        # inert on the diagonal, as above).  State |up dn, up> (i1=3, i2=1):
        #   U term: (1/z)[U*(n_up n_dn)_1 + U*(n_up n_dn)_2]
        #         = (1/4)[8*(1*1) + 8*(0*0)] = 2.0
        #   V term: V * n1 * n2 = 1.0 * 2 * 1 = 2.0   (n1 = 1+1=2, n2 = 1+0=1)
        #   h term: -(h/z)(Sz_1+Sz_2) = -(2/4)*(0 + 0.5) = -0.25
        #     (Sz_1 = 0.5*(1-1) = 0 for the doubly occupied site,
        #      Sz_2 = 0.5*(1-0) = 0.5 for the singly "up" site)
        #   mu term: 0 (mu not set)
        #   total: 2.0 + 2.0 - 0.25 = 3.75
        model = tenes_simple.make_model(hubbard_param({"u": 8.0, "v": 1.0, "h": 2.0}))
        h = model.bondhamiltonian(0, 0, z=4)
        assert h[3, 1, 3, 1] == pytest.approx(3.75)

    def test_hubbard_u_appears_on_doubly_occupied_sites(self):
        # U/z on each site, z = 4 -> site1 doubly occupied (i1=3), site2
        # empty (i2=0); then both sites doubly occupied.
        model = tenes_simple.make_model(hubbard_param({"u": 8.0}))
        h = model.bondhamiltonian(0, 0, z=4)
        assert h[3, 0, 3, 0] == pytest.approx(2.0)
        assert h[3, 3, 3, 3] == pytest.approx(4.0)

    def test_bond_hamiltonian_is_parity_conserving(self):
        # C6/robustness: all five couplings nonzero simultaneously, so a
        # parity-breaking cross term in any one of them cannot hide behind
        # another being zero.
        model = tenes_simple.make_model(
            hubbard_param({"t": 1.0, "u": 4.0, "v": 1.0, "mu": 2.0, "h": 0.7})
        )
        h = model.bondhamiltonian(0, 0, z=4)
        parity = model.parity
        for i1, i2, o1, o2 in np.ndindex(h.shape):
            if (parity[i1] ^ parity[i2]) != (parity[o1] ^ parity[o2]):
                assert h[i1, i2, o1, o2] == 0.0

    def test_bond_hamiltonian_is_hermitian(self):
        model = tenes_simple.make_model(
            hubbard_param({"t": 1.0, "u": 4.0, "v": 1.0, "mu": 2.0, "h": 0.7})
        )
        h = model.bondhamiltonian(0, 0, z=4)
        # rows are (in1, in2), columns are (out1, out2)
        m = h.reshape(16, 16)
        assert np.allclose(m, m.conj().T)

    def test_onesite_observables(self):
        model = tenes_simple.make_model(hubbard_param())
        assert model.onesite_ops_name == [
            "n",
            "n_up",
            "n_dn",
            "Sz",
            "doublon",
            "holon",
        ]
        ops = dict(zip(model.onesite_ops_name, model.onesite_ops))
        assert np.allclose(np.diag(ops["n"]), [0.0, 1.0, 1.0, 2.0])
        assert np.allclose(np.diag(ops["n_up"]), [0.0, 1.0, 0.0, 1.0])
        assert np.allclose(np.diag(ops["n_dn"]), [0.0, 0.0, 1.0, 1.0])
        assert np.allclose(np.diag(ops["Sz"]), [0.0, 0.5, -0.5, 0.0])
        assert np.allclose(np.diag(ops["doublon"]), [0.0, 0.0, 0.0, 1.0])
        assert np.allclose(np.diag(ops["holon"]), [1.0, 0.0, 0.0, 0.0])

    def test_every_onesite_observable_is_parity_even(self):
        model = tenes_simple.make_model(hubbard_param())
        parity = model.parity
        for op in model.onesite_ops:
            for i, o in np.ndindex(op.shape):
                if parity[i] != parity[o]:
                    assert op[i, o] == 0.0

    def test_nn_and_szsz_are_index_pair_products(self):
        # C3: "nn" and "SzSz" are two-site observables built as index-pair
        # products of the one-site ops (n is index 0, Sz is index 3 in the
        # onesite_ops_name order pinned above), not new operators.
        model = tenes_simple.make_model(hubbard_param())
        assert model.twosite_ops_name == ["nn", "SzSz"]
        assert model.twosite_ops == [(0, 0), (3, 3)]

    def test_hopping_is_an_explicit_rank4_observable(self):
        model = tenes_simple.make_model(hubbard_param())
        names = [name for name, _ in model.twosite_ops_explicit]
        assert "hopping" in names
        op = dict(model.twosite_ops_explicit)["hopping"]
        assert op.shape == (4, 4, 4, 4)
        # Same Jordan-Wigner computation as the C2 derivation above, but
        # without the -t prefactor (this is the bare operator sum_s
        # (c^dag_{1s} c_{2s} + h.c.), used to *measure* the hopping, not to
        # build H): <up dn, 0| c^dag_{1,up} c_{2,up} |dn, up> = -1, and no
        # other term in the sum connects |dn, up> to |up dn, 0>, so the
        # sign is NOT flipped the way it is in bondhamiltonian.
        assert op[2, 1, 3, 0] == pytest.approx(-1.0)


class TestHubbardSchema:
    def test_fermion_flag_is_injected(self):
        assert std_toml(hubbard_param())["parameter"]["general"]["fermion"] is True

    def test_parity_is_emitted_for_every_unitcell(self):
        parsed = std_toml(hubbard_param())
        for ucell in parsed["tensor"]["unitcell"]:
            assert ucell["parity"] == [0, 1, 1, 0]

    def test_no_twosite_observable_uses_the_ops_form(self):
        # C3: the hopping observable is inherently a rank-4 tensor (not an
        # outer product of one-site operators), and is_fermion = True routes
        # even the index-pair products (nn, SzSz) through the explicit
        # "elements" form instead of "ops = [...]".
        text, _ = tenes_simple.tenes_simple(hubbard_param())
        assert "ops = " not in text

    def test_bond_hamiltonian_observable_is_emitted_automatically(self):
        # C3: bond_hamiltonian is emitted by tenes_simple's own pipeline
        # (hamiltonians() + the [[observable.twosite]] block), not
        # constructed by the model itself.
        text, _ = tenes_simple.tenes_simple(hubbard_param())
        assert 'name = "bond_hamiltonian"' in text

    def test_random_initial_state_stays_scalar(self):
        parsed = std_toml(hubbard_param(lattice_extra={"initial": "random"}))
        for ucell in parsed["tensor"]["unitcell"]:
            assert ucell["initial_state"] == [0.0]

    def test_vacuum_initial_state(self):
        parsed = std_toml(hubbard_param(lattice_extra={"initial": "vacuum"}))
        for ucell in parsed["tensor"]["unitcell"]:
            assert np.allclose(ucell["initial_state"], [1.0, 0.0, 0.0, 0.0])

    def test_full_initial_state(self):
        parsed = std_toml(hubbard_param(lattice_extra={"initial": "full"}))
        for ucell in parsed["tensor"]["unitcell"]:
            assert np.allclose(ucell["initial_state"], [0.0, 0.0, 0.0, 1.0])

    def test_cdw_initial_state_alternates(self):
        # C4: cdw is the checkerboard |0> / |up dn> state, which needs TWO
        # sublattices (the contract's mandated SquareLattice change; see
        # TestSquareLatticeCdw below for that change pinned in isolation).
        parsed = std_toml(hubbard_param(lattice_extra={"initial": "cdw"}))
        states = [u["initial_state"] for u in parsed["tensor"]["unitcell"]]
        assert len(states) == 2
        assert np.allclose(states[0], [1.0, 0.0, 0.0, 0.0])
        assert np.allclose(states[1], [0.0, 0.0, 0.0, 1.0])

    @pytest.mark.parametrize("mode", ["ferro", "antiferro", "neel", "bogus"])
    def test_unsupported_initial_states_are_rejected(self, mode):
        # Match on the mode name itself so this cannot be satisfied by an
        # unrelated RuntimeError (e.g. "Unknown model type: hubbard" from a
        # missing dispatch branch) that happens to also be a RuntimeError.
        with pytest.raises(RuntimeError, match=re.escape(mode)):
            tenes_simple.tenes_simple(hubbard_param(lattice_extra={"initial": mode}))

    def test_unsupported_initial_state_message_gives_the_parity_reason(self):
        # C4: Neel-like odd-parity product states are impossible in this
        # scheme (the state vector sits on the even virtual index), and the
        # rejection message must say so, mirroring the spinless model's
        # message pattern (see SpinlessFermionModel.initial_state_vectors).
        with pytest.raises(RuntimeError) as excinfo:
            tenes_simple.tenes_simple(hubbard_param(lattice_extra={"initial": "neel"}))
        message = str(excinfo.value)
        assert re.search(r"\bparity\b", message, re.I)
        assert "M1" not in message and "M2" not in message

    def test_boson_model_with_cdw_initial_state_is_still_rejected(self):
        # C4: the base-class allowlist (random/ferro/antiferro only) must
        # keep rejecting "cdw" for a bosonic model; the SquareLattice/
        # HubbardModel changes for cdw must not leak into BoseHubbardModel.
        param = {
            "parameter": {"general": {}},
            "lattice": {
                "type": "square lattice",
                "L": 2,
                "W": 2,
                "virtual_dim": 2,
                "initial": "cdw",
            },
            "model": {"type": "boson", "t": 1.0},
        }
        with pytest.raises(RuntimeError, match="cdw"):
            tenes_simple.tenes_simple(param)


# ---------------------------------------------------------------------------
# C4's "KNOWN PLAN DEFECT": SquareLattice.__init__ currently has no branch
# for initial = "cdw" at all, so lattice.sublattice comes out empty (not an
# error -- SquareLattice itself never validates `initial`). This is tested
# directly against SquareLattice, independently of HubbardModel, because it
# fails differently from the rest of this file's tests: not with
# "Unknown model type: hubbard" (SquareLattice does not look at the model),
# but with an empty sublattice list.
# ---------------------------------------------------------------------------


class TestSquareLatticeCdw:
    def test_cdw_creates_two_sublattices(self):
        lattice = tenes_simple.SquareLattice(
            {"l": 2, "w": 2, "virtual_dim": 2, "initial": "cdw"}
        )
        assert len(lattice.sublattice) == 2

    def test_cdw_sublattices_are_checkerboarded(self):
        # Same partition rule as the existing "antiferro" branch: site
        # (x, y) goes to sublattice 0 when (x + y) is even, sublattice 1
        # otherwise.
        lattice = tenes_simple.SquareLattice(
            {"l": 2, "w": 2, "virtual_dim": 2, "initial": "cdw"}
        )
        sites0 = set(lattice.sublattice[0].sites)
        sites1 = set(lattice.sublattice[1].sites)
        assert sites0 | sites1 == set(range(4))
        assert sites0 & sites1 == set()
        for idx in sites0:
            x, y = tenes_simple.index2coord(idx, 2)
            assert (x + y) % 2 == 0
        for idx in sites1:
            x, y = tenes_simple.index2coord(idx, 2)
            assert (x + y) % 2 == 1


class TestHubbardScopeGuards:
    @pytest.mark.parametrize(
        "latname", ["honeycomb lattice", "triangular lattice", "kagome lattice"]
    )
    def test_non_square_lattices_are_rejected(self, latname):
        # C6, and the is_fermion-wiring proof the contract asks for: a
        # triangular lattice is rejected specifically *because*
        # HubbardModel.is_fermion is True and the guard is generic across
        # fermionic models, not because of anything spinless-specific.
        param = hubbard_param(lattice_extra={"type": latname})
        with pytest.raises(RuntimeError, match="square"):
            tenes_simple.tenes_simple(param)

    def test_square_lattice_is_accepted(self):
        tenes_simple.tenes_simple(hubbard_param())

    @pytest.mark.parametrize("key", ["t'", "t''", "v'", "v''"])
    def test_beyond_nearest_neighbour_parameters_are_rejected(self, key):
        param = hubbard_param({key: 0.5})
        with pytest.raises(RuntimeError, match="nearest"):
            tenes_simple.tenes_simple(param)

    def test_zero_valued_far_neighbour_parameters_are_accepted(self):
        tenes_simple.tenes_simple(hubbard_param({"t'": 0.0}))

    def test_bond_type_variants_of_the_first_neighbour_are_accepted(self):
        param = hubbard_param()
        param["model"] = {"type": "hubbard", "t0": 1.0}
        tenes_simple.tenes_simple(param)

    def test_duplicate_bond_type_specification_is_rejected(self):
        param = hubbard_param()
        param["model"] = {"type": "hubbard", "t": 1.0, "t0": 1.0}
        with pytest.raises(RuntimeError, match="defined twice"):
            tenes_simple.tenes_simple(param)

    def test_correlation_is_rejected(self):
        param = hubbard_param()
        param["correlation"] = {"r_max": 5, "operators": [[0, 0]]}
        with pytest.raises(RuntimeError, match="correlation"):
            tenes_simple.tenes_simple(param)

    def test_correlation_length_is_rejected(self):
        param = hubbard_param()
        param["correlation_length"] = {"measure": True}
        with pytest.raises(RuntimeError, match="correlation_length"):
            tenes_simple.tenes_simple(param)

    def test_the_message_does_not_mention_the_internal_milestone(self):
        param = hubbard_param(lattice_extra={"type": "honeycomb lattice"})
        with pytest.raises(RuntimeError) as excinfo:
            tenes_simple.tenes_simple(param)
        message = str(excinfo.value)
        assert "M1" not in message and "M2" not in message


class TestHubbardOneRowSkewedCell:
    # Section 2.3 of the contract cited above TestFermionOneRowSkewedCell,
    # for the second fermionic model: W = 1 is accepted (it used to be
    # refused on the retracted 2026-08-20 measurement; see that comment).
    def test_w1_cell_is_accepted_with_skew_1(self):
        text, lattice = tenes_simple.tenes_simple(hubbard_param(lattice_extra={"W": 1}))
        assert_one_row_skew_1_std(text, lattice, 2)
        assert_tenes_std_accepts(text, 2, [0, 1, 1, 0])

    def test_square_cell_no_skew_is_accepted(self):
        text, lattice = tenes_simple.tenes_simple(hubbard_param())
        assert lattice.skew == 0


# ---------------------------------------------------------------------------
# task-8-contract.md C5: the Hubbard gate must expand back to the bond
# Hamiltonian through the real pipeline, the same identity construction and
# tolerance reasoning as TestSpinlessFermionGateReducesToBondHamiltonian
# above (task-6-contract.md), generalised to a shared helper and exercised
# with all five Hubbard couplings nonzero so a dropped term cannot hide.
# ---------------------------------------------------------------------------


def assert_gate_expands_to_hamiltonian(param, tau=1e-6, atol=1e-4):
    text, _ = tenes_simple.tenes_simple(param)
    std_param = toml.loads(text)
    std_param.setdefault("parameter", {})
    std_param["parameter"]["simple_update"] = {"tau": [tau]}

    model = tenes_std.Model(std_param)
    ham = model.hamiltonians[0]
    evo = model.simple_updates[0]
    assert isinstance(ham, tenes_std.NNOperator)
    assert isinstance(evo, tenes_std.NNOperator)

    h = ham.elements
    gate = evo.elements
    d = h.shape[0]

    identity = np.zeros((d, d, d, d))
    for i1 in range(d):
        for i2 in range(d):
            identity[i1, i2, i1, i2] = 1.0

    approx_h = (gate - identity) / (-tau)
    assert np.allclose(approx_h, h, atol=atol)


def test_hubbard_gate_expands_to_the_hamiltonian():
    assert_gate_expands_to_hamiltonian(
        hubbard_param({"t": 1.0, "u": 4.0, "v": 0.5, "mu": 1.0, "h": 0.3})
    )


# ---------------------------------------------------------------------------
# task-11-contract.md C1: use_onesite_hamiltonian must not delete the
# Hubbard model's U/mu/h terms.
#
# Today HubbardModel.model_bondhamiltonian guards its onsite block on
# `use_onesite_hamiltonian`, and model_sitehamiltonian always returns zeros
# -- so tenes_simple(param, use_onesite_hamiltonian=True) with type =
# "hubbard" silently emits a bond Hamiltonian with U, mu and h dropped and
# no site Hamiltonian to carry them: wrong physics through the documented
# --use-site-hamiltonian CLI flag, with no error. The required fix folds
# the onsite terms into the bond Hamiltonian unconditionally (matching the
# spinless sibling, which already has no such conditional) and rejects the
# flag outright for fermionic models instead of silently corrupting the
# Hamiltonian.
# ---------------------------------------------------------------------------


class TestHubbardOnesiteTermsAlwaysInBondHamiltonian:
    def test_bond_hamiltonian_is_identical_for_both_flag_values(self):
        # The sharpest form of C1(a): the two flag values must produce the
        # SAME bond Hamiltonian element-for-element, not merely "both
        # nonzero somewhere".
        model = tenes_simple.make_model(
            hubbard_param({"t": 1.0, "u": 8.0, "v": 0.5, "mu": 1.0, "h": 0.3})
        )
        h_true = model.bondhamiltonian(0, 0, z=4, use_onesite_hamiltonian=True)
        h_false = model.bondhamiltonian(0, 0, z=4, use_onesite_hamiltonian=False)
        assert np.allclose(h_true, h_false)

    def test_u_carrying_element_is_present_regardless_of_the_flag(self):
        # Derivation (same as TestHubbardModel.test_hubbard_u_appears_on_
        # doubly_occupied_sites): u=8, z=4 -> U/z * doublon1 = 2.0 on the
        # site1-doubly-occupied (i1=3), site2-empty (i2=0) diagonal element.
        model = tenes_simple.make_model(hubbard_param({"u": 8.0}))
        h_true = model.bondhamiltonian(0, 0, z=4, use_onesite_hamiltonian=True)
        h_false = model.bondhamiltonian(0, 0, z=4, use_onesite_hamiltonian=False)
        assert h_true[3, 0, 3, 0] == pytest.approx(2.0)
        assert h_false[3, 0, 3, 0] == pytest.approx(2.0)


class TestFermionicModelsRejectUseOnesiteHamiltonian:
    def test_hubbard_rejects_use_onesite_hamiltonian(self):
        param = hubbard_param({"u": 8.0})
        with pytest.raises(RuntimeError) as excinfo:
            tenes_simple.tenes_simple(param, use_onesite_hamiltonian=True)
        message = str(excinfo.value)
        assert re.search(r"\bfermion", message, re.I)
        assert re.search(r"onesite_hamiltonian|one.?site", message, re.I)

    def test_spinless_rejects_use_onesite_hamiltonian(self):
        param = spinless_param({"mu": 1.0})
        with pytest.raises(RuntimeError) as excinfo:
            tenes_simple.tenes_simple(param, use_onesite_hamiltonian=True)
        message = str(excinfo.value)
        assert re.search(r"\bfermion", message, re.I)
        assert re.search(r"onesite_hamiltonian|one.?site", message, re.I)

    def test_message_does_not_mention_the_internal_milestone(self):
        param = hubbard_param({"u": 8.0})
        with pytest.raises(RuntimeError) as excinfo:
            tenes_simple.tenes_simple(param, use_onesite_hamiltonian=True)
        message = str(excinfo.value)
        assert "M1" not in message and "M2" not in message

    def test_spin_model_with_use_onesite_hamiltonian_still_works(self):
        # Regression net: C1's fix is gated on model.is_fermion, so a
        # bosonic/spin model must keep the flag exactly as today.
        param = {
            "parameter": {"general": {}},
            "lattice": {"type": "square lattice", "L": 2, "W": 2, "virtual_dim": 2},
            "model": {"type": "spin", "j": 1.0, "hz": 0.5},
        }
        text, lattice = tenes_simple.tenes_simple(param, use_onesite_hamiltonian=True)
        assert isinstance(text, str) and len(text) > 0


# ---------------------------------------------------------------------------
# task-11-contract.md C3: fermion = true with a non-fermionic model type
# must be rejected up front.
#
# Today [parameter.general] fermion = true with type = "spin"/"boson" is
# passed straight through into std.toml (no parity is emitted, since the
# model is not fermionic), and tenes_std then tells the user to ADD parity
# -- a confusing dead end. _check_fermion_scope must instead reject the
# combination itself, before its `is_fermion` early return.
# ---------------------------------------------------------------------------


class TestFermionFlagWithNonFermionicModelIsRejected:
    def test_spin_with_fermion_flag_is_rejected(self):
        param = {
            "parameter": {"general": {"fermion": True}},
            "lattice": {"type": "square lattice", "L": 2, "W": 2, "virtual_dim": 2},
            "model": {"type": "spin", "j": 1.0},
        }
        with pytest.raises(RuntimeError) as excinfo:
            tenes_simple.tenes_simple(param)
        message = str(excinfo.value)
        # Identifying content: the flag and/or the (non-fermionic) model
        # type must be named, not a generic "invalid input" string.
        assert re.search(r"\bfermion\b", message, re.I)
        assert "spin" in message.lower()

    def test_boson_with_fermion_flag_is_rejected(self):
        param = {
            "parameter": {"general": {"fermion": True}},
            "lattice": {"type": "square lattice", "L": 2, "W": 2, "virtual_dim": 2},
            "model": {"type": "boson", "t": 1.0},
        }
        with pytest.raises(RuntimeError) as excinfo:
            tenes_simple.tenes_simple(param)
        message = str(excinfo.value)
        assert re.search(r"\bfermion\b", message, re.I)
        assert "boson" in message.lower()

    def test_message_does_not_mention_the_internal_milestone(self):
        param = {
            "parameter": {"general": {"fermion": True}},
            "lattice": {"type": "square lattice", "L": 2, "W": 2, "virtual_dim": 2},
            "model": {"type": "spin", "j": 1.0},
        }
        with pytest.raises(RuntimeError) as excinfo:
            tenes_simple.tenes_simple(param)
        message = str(excinfo.value)
        assert "M1" not in message and "M2" not in message

    def test_fermionic_model_with_fermion_flag_still_works(self):
        # A fermionic model type with the flag explicitly set to True (the
        # matching, non-conflicting case) must keep working.
        param = spinless_param()
        param["parameter"]["general"]["fermion"] = True
        text, lattice = tenes_simple.tenes_simple(param)
        assert isinstance(text, str) and len(text) > 0

    def test_boson_without_the_flag_is_untouched(self):
        # Regression net: the ordinary bosonic path (no fermion key at
        # all) must not be touched by this new up-front check.
        param = {
            "parameter": {"general": {}},
            "lattice": {"type": "square lattice", "L": 2, "W": 2, "virtual_dim": 2},
            "model": {"type": "boson", "t": 1.0},
        }
        parsed = std_toml(param)
        assert "fermion" not in parsed["parameter"].get("general", {})


class TestHubbardGateIsExactlyParityEven:
    """expm(-tau h) of a parity-even h is parity-even, but eigh-based
    reconstruction leaves O(1e-16) noise in the odd blocks for the 16x16
    Hubbard gate.  The solver's fermion guard is a strict zero test, so
    that noise is a hard INPUT ERROR.  tenes_std must emit exact zeros."""

    def test_emitted_hubbard_gate_has_no_parity_odd_element(self):
        param = hubbard_param({"t": 1.0, "u": 4.0, "mu": 0.0})
        param["parameter"] = {
            "general": {"is_real": True},
            "simple_update": {"tau": 0.01, "num_step": 1},
            "ctm": {"dimension": 4},
        }
        text, _ = tenes_simple.tenes_simple(param)
        model = tenes_std.Model(toml.loads(text))
        parity = [0, 1, 1, 0]
        for evo in model.simple_updates:
            gate = evo.elements
            odd = []
            for i1, i2, o1, o2 in np.ndindex(gate.shape):
                if (parity[i1] ^ parity[i2]) != (parity[o1] ^ parity[o2]):
                    if gate[i1, i2, o1, o2] != 0.0:
                        odd.append(((i1, i2, o1, o2), gate[i1, i2, o1, o2]))
            assert odd == [], "parity-odd gate elements survive: %s" % odd[:4]


# ---------------------------------------------------------------------------
# The transverse spin correlations SxSx and SySy of the Hubbard model.
#
# An iPEPS breaks SU(2), so SzSz alone is not a third of <S_i . S_j>; to
# compare the nearest-neighbour spin correlation with the literature
# (Qin, Shi, Zhang, PRB 96, 075156 (2017)) the Hubbard model also offers
#
#     SxSx = S^x_1 S^x_2,   SySy = S^y_1 S^y_2,
#     S^a_j = (1/2) sum_{s,s'} c^dag_{j s} sigma^a_{s s'} c_{j s'},
#
# as explicit rank-4 observables appended AFTER "hopping" (so the group
# numbers of the existing observables do not move), in the layout of
# "hopping": op[in1, in2, out1, out2] = <out1 out2| O |in1 in2>.
#
# The expected operators are built here from the Pauli matrices alone, not
# from tenes_simple's Fock helpers. S^a_j is a fermion bilinear of site j,
# hence parity even, and its Jordan-Wigner strings cancel: in the ordered
# two-site basis |i1 i2> = (site-1 creators)(site-2 creators)|0> the product
# S^a_1 S^b_2 is the Kronecker product of the one-site matrices,
#
#     <o1 o2| S^a_1 S^b_2 |i1 i2> = S^a[o1, i1] S^b[o2, i2].
#
# On one site S^a annihilates |0> and |up dn> and acts on the doublet
# |up>, |dn> (local indices 1, 2) as sigma^a / 2.
# TestHubbardSpinReference checks this construction against an explicit
# Jordan-Wigner one, written in this file as well.
# ---------------------------------------------------------------------------

HUBBARD_UP = 1
HUBBARD_DN = 2
HUBBARD_DOUBLET = (HUBBARD_UP, HUBBARD_DN)

# mat[out, in], both indices in the order (up, dn)
PAULI = {
    "x": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    "y": np.array([[0.0, -1.0j], [1.0j, 0.0]]),
    "z": np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
}

TRANSVERSE = [("SxSx", "x"), ("SySy", "y")]


def hubbard_spin(axis):
    """S^axis on one Hubbard site as a 4x4 matrix mat[out, in]."""
    s = np.zeros((4, 4), dtype=complex)
    for a, out in enumerate(HUBBARD_DOUBLET):
        for b, inn in enumerate(HUBBARD_DOUBLET):
            s[out, inn] = 0.5 * PAULI[axis][a, b]
    return s


def product_bond_op(left, right):
    """op[in1, in2, out1, out2] = left[out1, in1] * right[out2, in2]."""
    return np.einsum("ai,bj->ijab", left, right)


def expected_spin_correlation(axis):
    """S^axis_1 S^axis_2 in the layout of the "hopping" observable."""
    op = product_bond_op(hubbard_spin(axis), hubbard_spin(axis))
    assert np.all(op.imag == 0.0)
    return op.real


def bond_op_as_matrix(op):
    """op[in1, in2, out1, out2] -> mat[out, in] on the 16-dimensional pair
    space, pair index 4 * i1 + i2 (the index of np.kron(site1, site2))."""
    d = op.shape[0]
    return op.reshape(d * d, d * d).T


def parse_elements(elements, shape):
    """An "indices... re im" elements block -> dense complex array."""
    op = np.zeros(shape, dtype=complex)
    n = len(shape)
    for line in elements.strip().splitlines():
        words = line.split()
        assert len(words) == n + 2, line
        index = tuple(int(w) for w in words[:n])
        op[index] = float(words[n]) + 1j * float(words[n + 1])
    return op


def bond_lines(bonds):
    return [line.split() for line in bonds.strip().splitlines()]


def twosite_entries(parsed):
    """observable.twosite entries of a parsed std.toml/input.toml, by name."""
    entries = {}
    for entry in parsed["observable"]["twosite"]:
        assert entry["name"] not in entries, entry["name"]
        entries[entry["name"]] = entry
    return entries


def hubbard_explicit_twosite(name):
    model = tenes_simple.make_model(hubbard_param())
    names = [n for n, _ in model.twosite_ops_explicit]
    assert name in names, "{} is not an explicit two-site observable: {}".format(
        name, names
    )
    return dict(model.twosite_ops_explicit)[name]


def emitted_hubbard_twosite(name, param=None):
    entries = twosite_entries(std_toml(param or hubbard_param()))
    assert name in entries, "{} is not emitted; observable.twosite has {}".format(
        name, sorted(entries)
    )
    return entries[name]


class TestHubbardSpinReference:
    """Checks the reference of this section, not tenes_simple."""

    def test_reference_is_the_second_quantized_definition(self):
        # Modes (site 1 up, site 1 dn, site 2 up, site 2 dn) = bits 0..3 of
        # the global occupation index g, so g = i1 + 4 * i2 with the local
        # index i = n_up + 2 n_dn, and
        #   |n_0 n_1 n_2 n_3> = (c^dag_0)^n_0 ... (c^dag_3)^n_3 |0>,
        # so c_m carries (-1)^(number of occupied modes below m).
        nmodes = 4
        dim = 1 << nmodes
        c = []
        for m in range(nmodes):
            cm = np.zeros((dim, dim))
            for g in range(dim):
                if (g >> m) & 1:
                    below = bin(g & ((1 << m) - 1)).count("1")
                    cm[g ^ (1 << m), g] = (-1.0) ** below
            c.append(cm)

        def spin(site, axis):
            s = np.zeros((dim, dim), dtype=complex)
            for a in range(2):
                for b in range(2):
                    cd = c[2 * site + a].T
                    s = s + 0.5 * PAULI[axis][a, b] * (cd @ c[2 * site + b])
            return s

        for axis in ("x", "y", "z"):
            fock = spin(0, axis) @ spin(1, axis)
            op = np.zeros((4, 4, 4, 4), dtype=complex)
            for i1, i2, o1, o2 in np.ndindex(op.shape):
                op[i1, i2, o1, o2] = fock[o1 + 4 * o2, i1 + 4 * i2]
            assert np.allclose(op, expected_spin_correlation(axis)), axis

    def test_reference_sz_is_the_model_sz(self):
        # ties HUBBARD_UP / HUBBARD_DN to the model's local basis
        model = tenes_simple.make_model(hubbard_param())
        sz = model.onesite_ops[model.onesite_ops_name.index("Sz")]
        assert np.allclose(hubbard_spin("z"), sz)


class TestHubbardTransverseSpinObservables:
    def test_hopping_stays_the_first_explicit_observable(self):
        # The emitted group numbers follow this order, so anything inserted
        # before "hopping" renumbers it.
        model = tenes_simple.make_model(hubbard_param())
        assert model.twosite_ops_explicit[0][0] == "hopping"

    def test_sxsx_and_sysy_follow_hopping(self):
        model = tenes_simple.make_model(hubbard_param())
        names = [name for name, _ in model.twosite_ops_explicit]
        assert names[0] == "hopping"
        assert sorted(names[1:]) == ["SxSx", "SySy"], names

    @pytest.mark.parametrize("name", ["SxSx", "SySy"])
    def test_is_a_real_hermitian_parity_even_rank4_operator(self, name):
        op = hubbard_explicit_twosite(name)
        assert op.shape == (4, 4, 4, 4)
        # exactly real: with is_real = true tenes_simple drops an explicit
        # observable that has any nonzero imaginary part
        assert np.all(np.isreal(op))
        m = op.reshape(16, 16)
        assert np.allclose(m, m.conj().T)
        parity = [0, 1, 1, 0]
        for i1, i2, o1, o2 in np.ndindex(op.shape):
            if (parity[i1] ^ parity[i2]) != (parity[o1] ^ parity[o2]):
                assert op[i1, i2, o1, o2] == 0.0, (i1, i2, o1, o2)

    def test_sxsx_named_elements(self):
        up, dn = HUBBARD_UP, HUBBARD_DN
        op = hubbard_explicit_twosite("SxSx")
        # op[in1, in2, out1, out2] = <out1 out2| SxSx |in1 in2>
        assert op[up, dn, dn, up] == pytest.approx(0.25)  # <dn,up|SxSx|up,dn>
        assert op[up, up, dn, dn] == pytest.approx(0.25)  # <dn,dn|SxSx|up,up>
        assert op[dn, up, up, dn] == pytest.approx(0.25)  # <up,dn|SxSx|dn,up>
        # S^x flips the spin on both sites, so nothing is diagonal; with the
        # two input legs swapped the 1/4 of <up,dn|SxSx|dn,up> lands here.
        assert op[up, dn, up, dn] == pytest.approx(0.0)  # <up,dn|SxSx|up,dn>

    def test_sysy_named_elements(self):
        # <dn|S^y|up> = i/2 and <up|S^y|dn> = -i/2, so a double flip of equal
        # spins gets (+-i/2)^2 = -1/4 and one of opposite spins +1/4.
        up, dn = HUBBARD_UP, HUBBARD_DN
        op = hubbard_explicit_twosite("SySy")
        assert op[up, up, dn, dn] == pytest.approx(-0.25)  # <dn,dn|SySy|up,up>
        assert op[dn, dn, up, up] == pytest.approx(-0.25)  # <up,up|SySy|dn,dn>
        assert op[up, dn, dn, up] == pytest.approx(0.25)  # <dn,up|SySy|up,dn>
        assert op[up, dn, up, dn] == pytest.approx(0.0)  # <up,dn|SySy|up,dn>

    @pytest.mark.parametrize("name", ["SxSx", "SySy"])
    def test_nothing_acts_on_an_empty_or_doubly_occupied_site(self, name):
        op = hubbard_explicit_twosite(name)
        for index in np.ndindex(op.shape):
            if any(i not in HUBBARD_DOUBLET for i in index):
                assert op[index] == 0.0, index

    @pytest.mark.parametrize("name, axis", TRANSVERSE)
    def test_is_the_reference_operator(self, name, axis):
        op = hubbard_explicit_twosite(name)
        assert np.allclose(op, expected_spin_correlation(axis), rtol=0.0, atol=1e-14)


class TestHubbardHeisenbergBond:
    """SxSx + SySy + SzSz as emitted (the existing SzSz included) is the
    spin-1/2 Heisenberg coupling S_1 . S_2 of the two singly occupied sites.
    The spectrum and the zero pattern do not use the reference above, and the
    commutator uses only its one-site spin matrices, so these also catch a
    slip that the reference and the implementation share."""

    SINGLY = [4 * i1 + i2 for i1 in HUBBARD_DOUBLET for i2 in HUBBARD_DOUBLET]

    def heisenberg(self):
        total = np.zeros((16, 16), dtype=complex)
        for name in ("SxSx", "SySy", "SzSz"):
            entry = emitted_hubbard_twosite(name)
            total += bond_op_as_matrix(parse_elements(entry["elements"], (4,) * 4))
        return total

    def test_spectrum_on_the_singly_occupied_states(self):
        h = self.heisenberg()
        block = h[np.ix_(self.SINGLY, self.SINGLY)]
        assert np.allclose(block, block.conj().T)
        assert np.allclose(np.linalg.eigvalsh(block), [-0.75, 0.25, 0.25, 0.25])
        # the -3/4 state is the singlet (|up,dn> - |dn,up>) / sqrt(2)
        singlet = np.zeros(16)
        singlet[4 * HUBBARD_UP + HUBBARD_DN] = 1.0 / np.sqrt(2.0)
        singlet[4 * HUBBARD_DN + HUBBARD_UP] = -1.0 / np.sqrt(2.0)
        assert np.allclose(h @ singlet, -0.75 * singlet)

    def test_zero_outside_the_singly_occupied_states(self):
        h = self.heisenberg()
        others = [g for g in range(16) if g not in self.SINGLY]
        assert np.allclose(h[others, :], 0.0)
        assert np.allclose(h[:, others], 0.0)

    @pytest.mark.parametrize("axis", ["x", "y", "z"])
    def test_commutes_with_the_total_spin(self, axis):
        h = self.heisenberg()
        s = hubbard_spin(axis)
        total = np.kron(s, np.eye(4)) + np.kron(np.eye(4), s)
        assert np.allclose(h @ total - total @ h, 0.0)


class TestHubbardTransverseSpinEmission:
    @pytest.mark.parametrize("name, axis", TRANSVERSE)
    def test_is_emitted_with_explicit_elements(self, name, axis):
        entry = emitted_hubbard_twosite(name)
        assert "ops" not in entry
        assert entry["dim"] == [4, 4]
        op = parse_elements(entry["elements"], (4,) * 4)
        assert np.allclose(op, expected_spin_correlation(axis), rtol=0.0, atol=1e-14)

    @pytest.mark.parametrize("name", ["SxSx", "SySy"])
    def test_is_emitted_on_the_bonds_of_the_other_twosite_observables(self, name):
        bonds = bond_lines(emitted_hubbard_twosite(name)["bonds"])
        for other in ("nn", "SzSz", "hopping"):
            assert bonds == bond_lines(emitted_hubbard_twosite(other)["bonds"])

    def test_existing_twosite_group_numbers_are_unchanged(self):
        entries = twosite_entries(std_toml(hubbard_param()))
        groups = {
            name: entry["group"]
            for name, entry in entries.items()
            if name not in ("SxSx", "SySy")
        }
        assert groups == {"bond_hamiltonian": 0, "nn": 1, "SzSz": 2, "hopping": 3}

    def test_new_observables_take_the_groups_after_hopping(self):
        groups = {emitted_hubbard_twosite(n)["group"] for n in ("SxSx", "SySy")}
        assert groups == {4, 5}

    def test_onesite_observables_and_groups_are_unchanged(self):
        onesite = std_toml(hubbard_param())["observable"]["onesite"]
        assert [(o["name"], o["group"]) for o in onesite] == [
            ("n", 0),
            ("n_up", 1),
            ("n_dn", 2),
            ("Sz", 3),
            ("doublon", 4),
            ("holon", 5),
        ]

    @pytest.mark.parametrize("name, axis", TRANSVERSE)
    def test_is_emitted_when_is_real_is_set(self, name, axis):
        # is_real = true drops every observable that is not exactly real
        param = hubbard_param()
        param["parameter"]["general"]["is_real"] = True
        entry = emitted_hubbard_twosite(name, param)
        op = parse_elements(entry["elements"], (4,) * 4)
        assert np.allclose(op, expected_spin_correlation(axis), rtol=0.0, atol=1e-14)

    def test_tenes_std_writes_both_into_the_input_toml(self):
        text, _ = tenes_simple.tenes_simple(hubbard_param())
        model = tenes_std.Model(toml.loads(text))
        buf = io.StringIO()
        model.to_toml(buf)
        entries = twosite_entries(toml.loads(buf.getvalue()))
        std_entries = twosite_entries(toml.loads(text))
        for name, axis in TRANSVERSE:
            assert name in entries, sorted(entries)
            entry = entries[name]
            assert "ops" not in entry
            assert entry["group"] == std_entries[name]["group"]
            assert bond_lines(entry["bonds"]) == bond_lines(entries["nn"]["bonds"])
            op = parse_elements(entry["elements"], (4,) * 4)
            assert np.allclose(op, expected_spin_correlation(axis), atol=1e-14)


class TestOtherModelsKeepTheirObservables:
    # The SxSx / SySy of the Hubbard model must not leak into the other
    # models: their observable lists, group numbers and forms are pinned.

    def test_spinless_fermion(self):
        model = tenes_simple.make_model(spinless_param())
        assert [name for name, _ in model.twosite_ops_explicit] == ["hopping"]
        parsed = std_toml(spinless_param())
        assert [(o["name"], o["group"]) for o in parsed["observable"]["onesite"]] == [
            ("n", 0)
        ]
        assert [
            (o["name"], o["group"], "ops" in o) for o in parsed["observable"]["twosite"]
        ] == [("bond_hamiltonian", 0, False), ("nn", 1, False), ("hopping", 2, False)]

    @pytest.mark.parametrize(
        "model_param, onesite, twosite",
        [
            (
                {"type": "spin", "j": 1.0},
                ["Sz", "Sx", "Sy"],
                [
                    ("bond_hamiltonian", 0, False),
                    ("SzSz", 1, True),
                    ("SxSx", 2, True),
                    ("SySy", 3, True),
                ],
            ),
            (
                {"type": "boson", "t": 1.0},
                ["N", "Bdagger", "B"],
                [
                    ("bond_hamiltonian", 0, False),
                    ("NN", 1, True),
                    ("BdaggerB", 2, True),
                    ("BBdagger", 3, True),
                ],
            ),
        ],
        ids=["spin", "boson"],
    )
    def test_bosonic_models(self, model_param, onesite, twosite):
        def param():
            return {
                "parameter": {"general": {}},
                "lattice": {"type": "square lattice", "L": 2, "W": 2, "virtual_dim": 2},
                "model": dict(model_param),
            }

        model = tenes_simple.make_model(param())
        assert model.twosite_ops_explicit == []
        parsed = std_toml(param())
        assert [o["name"] for o in parsed["observable"]["onesite"]] == onesite
        assert [
            (o["name"], o["group"], "ops" in o) for o in parsed["observable"]["twosite"]
        ] == twosite
