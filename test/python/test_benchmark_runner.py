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

import json
import os
import sys

import pytest

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmark"),
)

from benchlib import meta, runner, suite


def _make_stub(path, body):
    path.write_text("#!/bin/sh\n" + body + "\n")
    path.chmod(0o755)


@pytest.fixture
def stub_tools(tmp_path):
    """tenes_simple/tenes_std/tenes の代わりになるスタブ実行ファイル群"""
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


def test_run_context_resolves_relative_dirs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    ctx = runner.RunContext(
        tenes_dir="rel/src", tool_dir="rel/tool", results_dir="rel/res"
    )
    assert ctx.tenes_dir.is_absolute()
    assert ctx.tool_dir.is_absolute()
    assert ctx.results_dir.is_absolute()


def test_collect_meta_in_git_repo(tmp_path):
    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
    got = meta.collect_meta(tmp_path / "no_such_tenes", repo_dir=repo)
    assert len(got["git_commit"]) == 40
    assert isinstance(got["git_dirty"], bool)
    assert got["tenes_version"] is None  # 存在しないバイナリは None
    assert got["hostname"]


def test_prepare_input_renders_template(tmp_path, stub_tools):
    tenes_dir, tool_dir = stub_tools
    tpl = tmp_path / "tpl.toml"
    tpl.write_text("D = ${D}\n")
    case = suite.Case(name="c", kind="template", source=tpl, params={"D": 2})
    ctx = runner.RunContext(
        tenes_dir=tenes_dir, tool_dir=tool_dir, results_dir=tmp_path / "results"
    )
    workdir = tmp_path / "results" / "c" / "work"
    input_toml = runner.prepare_input(case, workdir, ctx)
    assert (workdir / "simple.toml").read_text() == "D = 2\n"
    assert input_toml.exists()


def test_run_case_stores_outputs_per_run(tmp_path, stub_tools):
    tenes_dir, tool_dir = stub_tools
    tpl = tmp_path / "tpl.toml"
    tpl.write_text("D = ${D}\n")
    case = suite.Case(name="c_D2", kind="template", source=tpl, params={"D": 2})
    ctx = runner.RunContext(
        tenes_dir=tenes_dir, tool_dir=tool_dir, results_dir=tmp_path / "results"
    )
    runner.run_case(case, repeat=2, ctx=ctx)
    for i in range(2):
        rundir = tmp_path / "results" / "c_D2" / "run_{}".format(i)
        assert (rundir / "timers.json").exists()
        assert (rundir / "density.dat").exists()


def test_run_suite_writes_meta(tmp_path, stub_tools):
    tenes_dir, tool_dir = stub_tools
    tpl = tmp_path / "tpl.toml"
    tpl.write_text("")
    s = suite.Suite(
        name="s",
        repeat=1,
        cases=[suite.Case(name="c", kind="template", source=tpl, params={})],
    )
    ctx = runner.RunContext(
        tenes_dir=tenes_dir, tool_dir=tool_dir, results_dir=tmp_path / "results"
    )
    runner.run_suite(s, ctx)
    got = json.loads((tmp_path / "results" / "meta.json").read_text())
    assert got["suite"] == "s"
    assert "git_commit" in got
