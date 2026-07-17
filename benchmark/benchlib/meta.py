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
import socket
import subprocess
from datetime import datetime, timezone


def _run(cmd, cwd=None):
    return subprocess.run(
        cmd, cwd=cwd, check=True, capture_output=True, text=True
    ).stdout.strip()


def collect_meta(tenes_bin, repo_dir, launcher=None):
    """Collect provenance metadata for a benchmark run."""
    result = {
        "date": datetime.now(timezone.utc).astimezone().isoformat(),
        "hostname": socket.gethostname(),
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
        "launcher": launcher,
    }
    try:
        result["git_commit"] = _run(["git", "rev-parse", "HEAD"], cwd=repo_dir)
        result["git_dirty"] = bool(_run(["git", "status", "--porcelain"], cwd=repo_dir))
    except (subprocess.CalledProcessError, OSError):
        result["git_commit"] = None
        result["git_dirty"] = None
    try:
        result["tenes_version"] = _run([str(tenes_bin), "--version"])
    except (subprocess.CalledProcessError, OSError):
        result["tenes_version"] = None
    return result
