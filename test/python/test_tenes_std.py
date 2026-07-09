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
