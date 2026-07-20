#!/usr/bin/env python3
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

"""TeNeS benchmark harness.

Run a suite and store results:
    bench.py run suites/contraction.toml --label baseline

Compare two result sets:
    bench.py compare results/baseline results/cotengra
"""

import argparse
import sys
from pathlib import Path

BENCH_DIR = Path(__file__).resolve().parent
REPO_DIR = BENCH_DIR.parent
sys.path.insert(0, str(BENCH_DIR))

from benchlib import compare as compare_mod
from benchlib import runner, suite


def cmd_run(args):
    suites = [suite.load_suite(path) for path in args.suite]
    results_dir = Path(args.results_dir) / args.label
    if results_dir.exists():
        sys.exit("error: {} already exists; use a new label".format(results_dir))
    ctx = runner.RunContext(
        tenes_dir=args.tenes_dir,
        tool_dir=args.tool_dir,
        results_dir=results_dir,
        launcher=args.launcher,
        repeat=args.repeat,
    )
    try:
        runner.run_suites(suites, ctx)
    except ValueError as e:
        sys.exit("error: {}".format(e))
    print("results saved to {}".format(results_dir))


def cmd_compare(args):
    report = compare_mod.compare_results(
        args.dir_a, args.dir_b, rtol=args.rtol, atol=args.atol
    )
    if args.output:
        Path(args.output).write_text(report)
        print("report saved to {}".format(args.output))
    else:
        print(report)


def main():
    parser = argparse.ArgumentParser(description="TeNeS benchmark harness")
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser("run", help="run one or more benchmark suites")
    p_run.add_argument(
        "suite", nargs="+", help="path(s) to suite TOML files (run sequentially)"
    )
    p_run.add_argument("--label", required=True, help="identifier of this run")
    p_run.add_argument("--tenes-dir", default=str(REPO_DIR / "build" / "src"))
    p_run.add_argument("--tool-dir", default=str(REPO_DIR / "build" / "tool"))
    p_run.add_argument("--results-dir", default=str(BENCH_DIR / "results"))
    p_run.add_argument("--launcher", default=None, help='e.g. "mpirun -np 4"')
    p_run.add_argument("--repeat", type=int, default=None, help="override suite repeat")
    p_run.set_defaults(func=cmd_run)

    p_cmp = sub.add_parser("compare", help="compare two result directories")
    p_cmp.add_argument("dir_a")
    p_cmp.add_argument("dir_b")
    p_cmp.add_argument("--output", default=None, help="write report to a file")
    p_cmp.add_argument("--rtol", type=float, default=1e-3)
    p_cmp.add_argument("--atol", type=float, default=1e-4)
    p_cmp.set_defaults(func=cmd_compare)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
