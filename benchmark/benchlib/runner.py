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
import shutil
import subprocess
from pathlib import Path

from . import meta as meta_mod
from . import suite as suite_mod


class RunContext:
    def __init__(self, tenes_dir, tool_dir, results_dir, launcher=None, repeat=None):
        # Resolve to absolute paths: commands run with cwd=workdir, so
        # relative paths given on the command line would not survive.
        self.tenes_dir = Path(tenes_dir).resolve()
        self.tool_dir = Path(tool_dir).resolve()
        self.results_dir = Path(results_dir).resolve()
        self.launcher = launcher
        self.repeat = repeat  # None -> suite default


def _call(cmd, cwd):
    cmd = [str(c) for c in cmd]
    print("  $ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def prepare_input(case, workdir, ctx):
    """Generate input.toml in workdir via tenes_simple -> tenes_std."""
    workdir.mkdir(parents=True, exist_ok=True)
    simple_toml = workdir / "simple.toml"
    if case.kind == "template":
        simple_toml.write_text(
            suite_mod.render_template(case.source.read_text(), case.params)
        )
    else:
        shutil.copy(case.source, simple_toml)
    _call([ctx.tool_dir / "tenes_simple", "simple.toml"], cwd=workdir)
    _call([ctx.tool_dir / "tenes_std", "std.toml"], cwd=workdir)
    return workdir / "input.toml"


def run_case(case, repeat, ctx):
    """Run one case `repeat` times and store outputs under run_<i>/."""
    casedir = ctx.results_dir / case.name
    workdir = casedir / "work"
    prepare_input(case, workdir, ctx)
    tenes_cmd = []
    if ctx.launcher:
        tenes_cmd += ctx.launcher.split()
    tenes_cmd += [ctx.tenes_dir / "tenes", "input.toml"]
    for i in range(repeat):
        _call(tenes_cmd, cwd=workdir)
        rundir = casedir / "run_{}".format(i)
        rundir.mkdir(parents=True, exist_ok=True)
        outdir = workdir / "output"
        for f in sorted(outdir.glob("*.dat")) + sorted(outdir.glob("*.json")):
            shutil.copy(f, rundir / f.name)


def run_suite(s, ctx):
    ctx.results_dir.mkdir(parents=True, exist_ok=True)
    repo_dir = Path(__file__).resolve().parents[2]
    result = meta_mod.collect_meta(
        ctx.tenes_dir / "tenes", repo_dir=repo_dir, launcher=ctx.launcher
    )
    result["suite"] = s.name
    (ctx.results_dir / "meta.json").write_text(json.dumps(result, indent=2))
    repeat = ctx.repeat if ctx.repeat is not None else s.repeat
    for case in s.cases:
        print("== case: {} (repeat={}) ==".format(case.name, repeat), flush=True)
        run_case(case, repeat, ctx)
