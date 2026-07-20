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

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmark"),
)

from benchlib import obscheck, stats


def _t(count, s):
    return {"count": count, "sum": s}


def test_aggregate_runs_median_min_max():
    runs = [
        {"a": _t(10, 3.0)},
        {"a": _t(10, 1.0)},
        {"a": _t(10, 2.0)},
    ]
    got = stats.aggregate_runs(runs)
    assert got["a"]["median"] == 2.0
    assert got["a"]["min"] == 1.0
    assert got["a"]["max"] == 3.0
    assert got["a"]["count"] == 10
    assert got["a"]["count_varies"] is False


def test_aggregate_runs_uses_slowest_rank_time():
    """MPI 実行では rank 0 の sum ではなく最遅 rank (max_rank) を統計に使う"""
    runs = [
        {"a": {"count": 10, "sum": 1.0, "max_rank": 3.0, "min_rank": 0.5}},
        {"a": {"count": 10, "sum": 1.2, "max_rank": 2.0, "min_rank": 0.4}},
    ]
    got = stats.aggregate_runs(runs)
    assert got["a"]["median"] == 2.5
    assert got["a"]["min"] == 2.0
    assert got["a"]["max"] == 3.0


def test_aggregate_runs_missing_name_and_count_change():
    runs = [
        {"a": _t(10, 1.0), "b": _t(1, 5.0)},
        {"a": _t(20, 2.0)},
    ]
    got = stats.aggregate_runs(runs)
    assert got["a"]["count_varies"] is True
    assert got["b"]["median"] == 5.0


def test_extract_numbers_skips_comments():
    text = "# header 123\n1.0 2.5e-3\n-4 hello5\n"
    assert obscheck.extract_numbers(text) == [1.0, 2.5e-3, -4.0, 5.0]


def test_compare_dat_dirs(tmp_path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir()
    b.mkdir()
    (a / "density.dat").write_text("1.0 2.0\n")
    (b / "density.dat").write_text("1.0001 2.0\n")
    (a / "time.dat").write_text("999\n")
    (b / "time.dat").write_text("1\n")
    assert obscheck.compare_dat_dirs(a, b) == []

    (b / "density.dat").write_text("1.5 2.0\n")
    warnings = obscheck.compare_dat_dirs(a, b)
    assert len(warnings) == 1
    assert "density.dat" in warnings[0]

    (b / "extra.dat").write_text("0\n")
    warnings = obscheck.compare_dat_dirs(a, b)
    assert any("extra.dat" in w for w in warnings)
