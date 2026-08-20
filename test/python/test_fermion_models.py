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
# C1 (task-4b-contract.md): a fermionic model on a skewed cell must be
# rejected.
#
# work/skew-validation/FINDINGS.md measured this, it is not a precaution:
# free spinless fermions on the standard two-site cell L_sub = [2, 1] (which
# SquareLattice realises as skew = 1 when W == 1) give a 20.6% energy shift
# and a density drifted off half filling relative to the flat, skew = 0
# control -- with no error, no warning, a plausible-looking wrong number.
# Until the C++ sign bug is fixed, tenes_simple must refuse to emit a
# std.toml for this combination.
# ---------------------------------------------------------------------------


class TestFermionSkewGuard:
    def test_fermionic_model_with_skewed_cell_is_rejected(self):
        # W = 1 is the standard two-site cell; SquareLattice sets skew = 1
        # for it. This must raise before any std.toml text is produced.
        param = spinless_param(lattice_extra={"W": 1})
        with pytest.raises(RuntimeError) as excinfo:
            tenes_simple.tenes_simple(param)
        message = str(excinfo.value)
        # Identifying content, not just "raises something": name the
        # offending cell (W = 1) or the skew value it produces (both are
        # literally 1 here), and steer the user to a wider cell (W >= 2).
        assert re.search(r"\bskew\b", message, re.I) or re.search(
            r"\bW\s*=\s*1\b", message
        )
        assert re.search(r"\b2\b", message)

    def test_fermionic_model_with_skewed_cell_message_is_fermion_specific(self):
        # C1 requires the message to say this is a MEASURED limitation of
        # the fermion implementation (wrong numbers), not a generic
        # "unsupported configuration" the way e.g. the non-square-lattice
        # guard reads. Pin that framing without pinning exact prose.
        param = spinless_param(lattice_extra={"W": 1})
        with pytest.raises(RuntimeError) as excinfo:
            tenes_simple.tenes_simple(param)
        message = str(excinfo.value)
        assert re.search(r"\bfermion", message, re.I)
        assert re.search(r"measured|wrong|incorrect|known limitation", message, re.I)
        # The existing missing-feature guards in this same function say
        # exactly this; the skew bug must not be described the same way,
        # since it is a correctness bug, not missing scope.
        assert "is not available for fermionic models" not in message

    def test_skew_guard_message_does_not_mention_the_internal_milestone(self):
        param = spinless_param(lattice_extra={"W": 1})
        with pytest.raises(RuntimeError) as excinfo:
            tenes_simple.tenes_simple(param)
        message = str(excinfo.value)
        assert "M1" not in message
        assert "M2" not in message

    def test_fermionic_model_with_square_cell_no_skew_is_accepted(self):
        # Regression net (i): fermionic + skew = 0 (the default 2x2 cell)
        # must keep working exactly as today.
        text, lattice = tenes_simple.tenes_simple(spinless_param())
        assert lattice.skew == 0

    def test_fermionic_model_with_wide_cell_no_skew_is_accepted(self):
        # Regression net (i), a second shape: any W != 1 keeps skew = 0.
        param = spinless_param(lattice_extra={"L": 3, "W": 3})
        text, lattice = tenes_simple.tenes_simple(param)
        assert lattice.skew == 0

    def test_bosonic_model_with_skewed_cell_is_still_accepted(self):
        # Regression net (ii): the bug is fermion-specific. W = 1 is
        # standard practice for spins/bosons and must be untouched.
        param = {
            "parameter": {"general": {}},
            "lattice": {"type": "square lattice", "L": 2, "W": 1, "virtual_dim": 2},
            "model": {"type": "spin", "j": 1.0},
        }
        text, lattice = tenes_simple.tenes_simple(param)
        assert lattice.skew == 1


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


class TestHubbardSkewGuard:
    def test_skewed_cell_is_rejected(self):
        param = hubbard_param(lattice_extra={"W": 1})
        with pytest.raises(RuntimeError) as excinfo:
            tenes_simple.tenes_simple(param)
        message = str(excinfo.value)
        assert re.search(r"\bskew\b", message, re.I) or re.search(
            r"\bW\s*=\s*1\b", message
        )
        assert "M1" not in message and "M2" not in message

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
