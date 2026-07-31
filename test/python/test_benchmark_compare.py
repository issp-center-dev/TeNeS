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

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmark"),
)

from benchlib import compare


def _write_label(root, label, sums, count=10):
    """テスト用のラベルディレクトリ(1ケース・2反復)を作る"""
    labeldir = root / label
    for i, s in enumerate(sums):
        rundir = labeldir / "case1" / "run_{}".format(i)
        rundir.mkdir(parents=True)
        timers = {
            "total": {
                "count": 1,
                "sum": 10 * s,
                "max_rank": 10 * s,
                "min_rank": 10 * s,
            },
            "phase/environment": {
                "count": 1,
                "sum": 5 * s,
                "max_rank": 5 * s,
                "min_rank": 5 * s,
            },
            "contract/itps_ctm/2x2": {
                "count": count,
                "sum": s,
                "max_rank": s,
                "min_rank": s,
            },
        }
        (rundir / "timers.json").write_text(json.dumps({"meta": {}, "timers": timers}))
        (rundir / "density.dat").write_text("1.0\n")
    (labeldir / "meta.json").write_text(
        json.dumps({"git_commit": label + "0" * 34, "hostname": "h"})
    )
    return labeldir


def test_compare_timers_ratio_and_overlap():
    agg_a = {
        "x": {"count": 1, "count_varies": False, "median": 2.0, "min": 1.9, "max": 2.1}
    }
    agg_b = {
        "x": {"count": 1, "count_varies": False, "median": 1.0, "min": 0.9, "max": 1.1}
    }
    rows = compare.compare_timers(agg_a, agg_b)
    assert rows[0]["ratio"] == 0.5
    assert rows[0]["overlap"] is False
    assert rows[0]["count_mismatch"] is False


def test_compare_timers_one_side_only():
    rows = compare.compare_timers(
        {
            "only_a": {
                "count": 1,
                "count_varies": False,
                "median": 1.0,
                "min": 1.0,
                "max": 1.0,
            }
        },
        {},
    )
    assert rows[0]["a"] is not None
    assert rows[0]["b"] is None


def test_compare_results_end_to_end(tmp_path):
    dir_a = _write_label(tmp_path, "base", [1.0, 1.2])
    dir_b = _write_label(tmp_path, "new", [0.5, 0.6])
    report = compare.compare_results(dir_a, dir_b)
    assert "# Benchmark comparison" in report
    assert "total" in report
    assert "## contract" in report
    assert "case1" in report
    assert "WARNING" not in report


def test_compare_results_detects_observable_mismatch(tmp_path):
    dir_a = _write_label(tmp_path, "base", [1.0])
    dir_b = _write_label(tmp_path, "new", [1.0])
    rundir = dir_b / "case1" / "run_0"
    (rundir / "density.dat").write_text("2.0\n")
    report = compare.compare_results(dir_a, dir_b)
    assert "WARNING" in report
    assert "density.dat" in report


def test_compare_results_flags_count_mismatch(tmp_path):
    dir_a = _write_label(tmp_path, "base", [1.0], count=10)
    dir_b = _write_label(tmp_path, "new", [1.0], count=20)
    report = compare.compare_results(dir_a, dir_b)
    assert "count differs" in report


def test_format_row_one_side_only_column_order():
    stat = {
        "count": 1,
        "count_varies": False,
        "median": 42.0,
        "min": 40.0,
        "max": 44.0,
    }
    rows = compare.compare_timers({}, {"x": stat})
    cells = [
        c.strip() for c in compare._format_row("case1", rows[0]).strip("|").split("|")
    ]
    assert cells[2] == "(absent)"
    assert cells[3].startswith("42")
    assert cells[6] == "B only"

    rows = compare.compare_timers({"x": stat}, {})
    cells = [
        c.strip() for c in compare._format_row("case1", rows[0]).strip("|").split("|")
    ]
    assert cells[2].startswith("42")
    assert cells[3] == "(absent)"
    assert cells[6] == "A only"


def _write_ab_label(root, cases):
    """テスト用のラベルディレクトリ(ケース名 -> 反復ごとのスケール係数)を作る"""
    labeldir = root / "lbl"
    for case, sums in cases.items():
        for i, s in enumerate(sums):
            rundir = labeldir / case / "run_{}".format(i)
            rundir.mkdir(parents=True)
            timers = {
                "total": {
                    "count": 1,
                    "sum": 10 * s,
                    "max_rank": 10 * s,
                    "min_rank": 10 * s,
                },
                "measure/correlation_length": {
                    "count": 1,
                    "sum": s,
                    "max_rank": s,
                    "min_rank": s,
                },
            }
            (rundir / "timers.json").write_text(
                json.dumps({"meta": {}, "timers": timers})
            )
    (labeldir / "meta.json").write_text(json.dumps({"hostname": "h", "suite": "s"}))
    return labeldir


def test_show_results_pairs_ab_cases(tmp_path):
    labeldir = _write_ab_label(
        tmp_path,
        {"cl_builtin_D3": [1.0, 1.0], "cl_arpack_D3": [0.5, 0.5]},
    )
    report = compare.show_results(labeldir)
    assert "cl_builtin_D3 vs cl_arpack_D3" in report
    assert "measure/correlation_length" in report
    assert "0.500" in report  # arpack/builtin の比
    assert "hostname" in report


def test_show_results_unpaired_case_gets_plain_table(tmp_path):
    labeldir = _write_ab_label(tmp_path, {"solo_case": [2.0]})
    report = compare.show_results(labeldir)
    assert "## solo_case" in report
    assert "vs" not in report
    assert "total" in report


def test_show_results_ab_disabled(tmp_path):
    labeldir = _write_ab_label(tmp_path, {"cl_builtin": [1.0], "cl_arpack": [1.0]})
    report = compare.show_results(labeldir, ab=("", ""))
    assert "vs" not in report
    assert "## cl_arpack" in report
    assert "## cl_builtin" in report
