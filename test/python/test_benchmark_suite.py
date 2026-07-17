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
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmark"),
)

from benchlib import suite


def test_expand_sweep_product():
    got = suite.expand_sweep({"D": [2, 4], "chi_ratio": [1, 2]})
    assert len(got) == 4
    assert {"D": 2, "chi_ratio": 1} in got
    assert {"D": 4, "chi_ratio": 2} in got


def test_expand_sweep_with_list_values():
    got = suite.expand_sweep({"Lsub": [[1, 1], [2, 2]]})
    assert got == [{"Lsub": [1, 1]}, {"Lsub": [2, 2]}]


def test_derive_params_lsub_and_chi():
    got = suite.derive_params({"Lsub": [2, 3], "D": 4, "chi_ratio": 2})
    assert got["L"] == 2
    assert got["W"] == 3
    assert got["chi"] == 32
    assert "Lsub" not in got
    assert "chi_ratio" not in got


def test_format_name():
    name = suite.format_name(
        "afh_${L}x${W}_D${D}_chi${chi}", {"L": 2, "W": 2, "D": 4, "chi": 32}
    )
    assert name == "afh_2x2_D4_chi32"


def test_load_suite_expands_template_cases(tmp_path):
    (tmp_path / "tpl.toml").write_text("D = ${D}\nsteps = ${steps}\n")
    suite_toml = tmp_path / "suite.toml"
    suite_toml.write_text("""
[suite]
name = "s"
repeat = 2

[[case]]
name = "c_D${D}"
template = "tpl.toml"
params = { steps = 10 }
sweep = { D = [2, 4] }
""")
    s = suite.load_suite(suite_toml)
    assert s.name == "s"
    assert s.repeat == 2
    assert [c.name for c in s.cases] == ["c_D2", "c_D4"]
    assert s.cases[0].kind == "template"
    assert s.cases[0].params == {"steps": 10, "D": 2}
    assert s.cases[0].source == tmp_path / "tpl.toml"


def test_load_suite_input_case(tmp_path):
    (tmp_path / "in.toml").write_text("")
    suite_toml = tmp_path / "suite.toml"
    suite_toml.write_text("""
[suite]
name = "s"

[[case]]
name = "existing"
input = "in.toml"
""")
    s = suite.load_suite(suite_toml)
    assert s.repeat == 1
    assert s.cases[0].kind == "input"
    assert s.cases[0].source == tmp_path / "in.toml"


def test_load_suite_rejects_duplicate_names(tmp_path):
    (tmp_path / "tpl.toml").write_text("")
    suite_toml = tmp_path / "suite.toml"
    suite_toml.write_text("""
[suite]
name = "s"

[[case]]
name = "same"
template = "tpl.toml"

[[case]]
name = "same"
template = "tpl.toml"
""")
    with pytest.raises(ValueError):
        suite.load_suite(suite_toml)


def test_render_template():
    assert suite.render_template("L = ${L}\n", {"L": 2}) == "L = 2\n"
