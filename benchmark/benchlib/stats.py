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

import statistics


def aggregate_runs(run_timers):
    """Aggregate repeated runs of one case.

    run_timers: list of timer dicts (name -> {"count": int, "sum": float, ...}).
    Returns name -> {"count", "count_varies", "median", "min", "max"}.
    Names missing from some runs are aggregated over the runs that have them.

    The statistics use the slowest rank's cumulative time (max_rank) so that
    MPI comparisons reflect the wall-clock bottleneck, not rank 0's share;
    entries without max_rank fall back to sum (identical for single-rank runs).
    """
    names = []
    for timers in run_timers:
        for name in timers:
            if name not in names:
                names.append(name)
    result = {}
    for name in names:
        elapsed = [
            t[name].get("max_rank", t[name]["sum"]) for t in run_timers if name in t
        ]
        counts = [t[name]["count"] for t in run_timers if name in t]
        result[name] = {
            "count": counts[0],
            "count_varies": len(set(counts)) > 1 or len(elapsed) != len(run_timers),
            "median": statistics.median(elapsed),
            "min": min(elapsed),
            "max": max(elapsed),
        }
    return result
