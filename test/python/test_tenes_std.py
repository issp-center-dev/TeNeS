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

import tenes_std


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
