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
import toml

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tool")
)

import tenes_simple


def spinless_param(model_extra=None, lattice_extra=None):
    model = {"type": "spinless fermion", "t": 1.0}
    model.update(model_extra or {})
    lattice = {"type": "square lattice", "L": 2, "W": 2, "virtual_dim": 2}
    lattice.update(lattice_extra or {})
    return {"parameter": {"general": {}}, "lattice": lattice, "model": model}


def std_toml(param):
    text, _ = tenes_simple.tenes_simple(param)
    return toml.loads(text)


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
