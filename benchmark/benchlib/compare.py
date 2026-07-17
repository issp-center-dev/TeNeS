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
from pathlib import Path

from . import obscheck, stats

_META_KEYS = (
    "git_commit",
    "git_dirty",
    "tenes_version",
    "hostname",
    "date",
    "omp_num_threads",
    "launcher",
    "suite",
)


def load_label_dir(label_dir):
    """Load meta.json and per-case aggregated timers from a label directory."""
    label_dir = Path(label_dir)
    meta_path = label_dir / "meta.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    cases = {}
    for casedir in sorted(p for p in label_dir.iterdir() if p.is_dir()):
        runs = []
        for rundir in sorted(casedir.glob("run_*")):
            timers_path = rundir / "timers.json"
            if timers_path.exists():
                runs.append(json.loads(timers_path.read_text())["timers"])
        if runs:
            cases[casedir.name] = {
                "agg": stats.aggregate_runs(runs),
                "run0": casedir / "run_0",
            }
    return meta, cases


def compare_timers(agg_a, agg_b):
    """Compare two aggregate dicts name-by-name."""
    rows = []
    for name in sorted(set(agg_a) | set(agg_b)):
        a = agg_a.get(name)
        b = agg_b.get(name)
        row = {"name": name, "a": a, "b": b}
        if a is not None and b is not None:
            row["ratio"] = b["median"] / a["median"] if a["median"] > 0 else None
            row["overlap"] = a["min"] <= b["max"] and b["min"] <= a["max"]
            row["count_mismatch"] = a["count"] != b["count"]
        rows.append(row)
    return rows


def _group_key(name):
    """Group timers: total and phase/* into 'summary', others by prefix."""
    if name == "total" or name.startswith("phase/"):
        return "summary"
    return name.split("/", 1)[0] if "/" in name else "other"


def _stat_cell(s):
    return "{:.4g} [{:.4g}, {:.4g}]".format(s["median"], s["min"], s["max"])


def _format_row(case, row):
    a, b = row["a"], row["b"]
    if a is None:
        return "| {} | {} | {} | (absent) | n/a |  | B only |".format(
            case, row["name"], _stat_cell(b)
        )
    if b is None:
        return "| {} | {} | {} | (absent) | n/a |  | A only |".format(
            case, row["name"], _stat_cell(a)
        )
    ratio = "{:.3f}".format(row["ratio"]) if row["ratio"] is not None else "n/a"
    notes = []
    if row["overlap"]:
        notes.append("no sig. diff")
    if row["count_mismatch"]:
        notes.append("count differs")
    if a.get("count_varies") or b.get("count_varies"):
        notes.append("count varies across runs")
    return "| {} | {} | {} | {} | {} | {}/{} | {} |".format(
        case,
        row["name"],
        _stat_cell(a),
        _stat_cell(b),
        ratio,
        a["count"],
        b["count"],
        ", ".join(notes),
    )


_TABLE_HEADER = (
    "| case | timer | A: median [min, max] | B: median [min, max]"
    " | ratio B/A | count A/B | note |\n|---|---|---|---|---|---|---|"
)


def compare_results(dir_a, dir_b, rtol=1e-3, atol=1e-4):
    """Render a Markdown report comparing two label directories."""
    meta_a, cases_a = load_label_dir(dir_a)
    meta_b, cases_b = load_label_dir(dir_b)
    lines = []
    lines.append(
        "# Benchmark comparison: {} (A) vs {} (B)".format(
            Path(dir_a).name, Path(dir_b).name
        )
    )
    lines.append("")
    lines.append("| meta | A | B |")
    lines.append("|---|---|---|")
    for key in _META_KEYS:
        lines.append("| {} | {} | {} |".format(key, meta_a.get(key), meta_b.get(key)))
    lines.append("")

    common = sorted(set(cases_a) & set(cases_b))
    only = sorted(set(cases_a) ^ set(cases_b))
    if only:
        lines.append("**WARNING: cases on one side only:** " + ", ".join(only))
        lines.append("")

    obs_warnings = []
    for case in common:
        for w in obscheck.compare_dat_dirs(
            cases_a[case]["run0"], cases_b[case]["run0"], rtol=rtol, atol=atol
        ):
            obs_warnings.append("{}: {}".format(case, w))
    if obs_warnings:
        lines.append("## WARNING: observable mismatches")
        lines.append("")
        for w in obs_warnings:
            lines.append("- " + w)
        lines.append("")

    groups = {}
    for case in common:
        rows = compare_timers(cases_a[case]["agg"], cases_b[case]["agg"])
        for row in rows:
            groups.setdefault(_group_key(row["name"]), []).append((case, row))

    group_names = sorted(groups, key=lambda g: (g != "summary", g))
    for group in group_names:
        title = "Summary (total / phase)" if group == "summary" else group + "/"
        lines.append("## " + title)
        lines.append("")
        lines.append(_TABLE_HEADER)
        for case, row in groups[group]:
            lines.append(_format_row(case, row))
        lines.append("")
    return "\n".join(lines)
