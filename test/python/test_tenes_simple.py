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

import pytest

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tool")
)

import tenes_simple


class TestBoseHubbardModel:
    def test_site_hamiltonian_has_onsite_repulsion(self):
        model = tenes_simple.BoseHubbardModel({"type": "boson", "nmax": 2, "u": 3.0})
        ham = model.sitehamiltonian()
        # H = -mu n + U/2 n(n-1); for n = 2: U/2 * 2 * 1 = U
        assert ham[2, 2].real == pytest.approx(3.0)
        assert ham[2, 2].imag == pytest.approx(0.0)
        assert ham[1, 1].real == pytest.approx(0.0)

    def test_twosite_ops_names(self):
        model = tenes_simple.BoseHubbardModel({"type": "boson"})
        assert model.twosite_ops_name == ["NN", "BdaggerB", "BBdagger"]

    def test_site_hamiltonian_has_chemical_potential(self):
        model = tenes_simple.BoseHubbardModel({"type": "boson", "nmax": 2, "mu": 2.0})
        ham = model.sitehamiltonian()
        assert ham[1, 1].real == pytest.approx(-2.0)
        assert ham[2, 2].real == pytest.approx(-4.0)
