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

import re
from pathlib import Path

import numpy as np

_FLOAT_RE = re.compile(r"[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?")


def extract_numbers(text):
    """Extract numeric tokens, ignoring comment lines starting with '#'."""
    numbers = []
    for line in text.splitlines():
        if line.lstrip().startswith("#"):
            continue
        for m in _FLOAT_RE.finditer(line):
            numbers.append(float(m.group(0)))
    return numbers


def compare_dat_dirs(
    dir_a, dir_b, rtol=1e-3, atol=1e-4, exclude=("time.dat", "parameters.dat")
):
    """Numerically compare common .dat files. Returns a list of warnings."""
    dir_a = Path(dir_a)
    dir_b = Path(dir_b)
    files_a = {p.name for p in dir_a.glob("*.dat")} - set(exclude)
    files_b = {p.name for p in dir_b.glob("*.dat")} - set(exclude)
    warnings = []
    for name in sorted(files_a & files_b):
        num_a = extract_numbers((dir_a / name).read_text())
        num_b = extract_numbers((dir_b / name).read_text())
        if len(num_a) != len(num_b):
            warnings.append(
                "{}: number of values differs ({} vs {})".format(
                    name, len(num_a), len(num_b)
                )
            )
        elif not np.allclose(num_a, num_b, rtol=rtol, atol=atol):
            warnings.append("{}: values differ beyond tolerance".format(name))
    for name in sorted(files_a ^ files_b):
        warnings.append("{}: missing on one side".format(name))
    return warnings
