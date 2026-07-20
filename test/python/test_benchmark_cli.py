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

import argparse
import os
import sys

import pytest

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmark"),
)

import bench


def _make_stub(path, body):
    path.write_text("#!/bin/sh\n" + body + "\n")
    path.chmod(0o755)


@pytest.fixture
def stub_tools(tmp_path):
    tool_dir = tmp_path / "tool"
    tool_dir.mkdir()
    _make_stub(tool_dir / "tenes_simple", "touch std.toml")
    _make_stub(tool_dir / "tenes_std", "touch input.toml")
    tenes_dir = tmp_path / "src"
    tenes_dir.mkdir()
    _make_stub(
        tenes_dir / "tenes",
        'mkdir -p output && echo "{}" > output/timers.json'
        ' && echo "0.5" > output/density.dat',
    )
    return tenes_dir, tool_dir


def _write_suite(tmp_path):
    (tmp_path / "tpl.toml").write_text("")
    suite_toml = tmp_path / "suite.toml"
    suite_toml.write_text(
        '[suite]\nname = "s"\n\n[[case]]\nname = "c"\ntemplate = "tpl.toml"\n'
    )
    return suite_toml


def _args(tmp_path, suite_toml, tenes_dir, tool_dir, force=False):
    return argparse.Namespace(
        suite=[str(suite_toml)],
        label="lbl",
        tenes_dir=str(tenes_dir),
        tool_dir=str(tool_dir),
        results_dir=str(tmp_path / "results"),
        launcher=None,
        repeat=None,
        force=force,
    )


def test_run_refuses_existing_label_without_force(tmp_path, stub_tools):
    tenes_dir, tool_dir = stub_tools
    suite_toml = _write_suite(tmp_path)
    (tmp_path / "results" / "lbl").mkdir(parents=True)
    with pytest.raises(SystemExit, match="already exists"):
        bench.cmd_run(_args(tmp_path, suite_toml, tenes_dir, tool_dir))


def test_run_force_replaces_existing_label(tmp_path, stub_tools):
    tenes_dir, tool_dir = stub_tools
    suite_toml = _write_suite(tmp_path)
    label_dir = tmp_path / "results" / "lbl"
    label_dir.mkdir(parents=True)
    stale = label_dir / "stale_case_marker"
    stale.write_text("old")
    bench.cmd_run(_args(tmp_path, suite_toml, tenes_dir, tool_dir, force=True))
    assert not stale.exists()
    assert (label_dir / "meta.json").exists()
    assert (label_dir / "c" / "run_0" / "timers.json").exists()
