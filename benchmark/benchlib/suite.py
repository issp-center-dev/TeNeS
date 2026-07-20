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

import itertools
import string
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List


def load_toml(path):
    """Load a TOML file with tomllib (Python >= 3.11) or toml as fallback."""
    try:
        import tomllib
    except ImportError:
        import toml

        return toml.load(str(path))
    with open(path, "rb") as f:
        return tomllib.load(f)


@dataclass
class Case:
    name: str
    kind: str  # "template" or "input"
    source: Path
    params: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Suite:
    name: str
    repeat: int
    cases: List[Case]


def expand_sweep(sweep: Dict[str, List[Any]]) -> List[Dict[str, Any]]:
    """Expand {key: [values...]} into the direct product as a list of dicts."""
    keys = list(sweep.keys())
    values = [sweep[k] for k in keys]
    return [dict(zip(keys, comb)) for comb in itertools.product(*values)]


def derive_params(params: Dict[str, Any]) -> Dict[str, Any]:
    """Resolve derived parameters.

    - Lsub = [L, W]   -> L, W
    - chi_ratio = r   -> chi = r * D * D (requires D)
    """
    ret = dict(params)
    if "Lsub" in ret:
        ret["L"], ret["W"] = ret.pop("Lsub")
    if "chi_ratio" in ret:
        ret["chi"] = ret.pop("chi_ratio") * ret["D"] ** 2
    return ret


def _substitute(text: str, params: Dict[str, Any]) -> str:
    return string.Template(text).substitute({k: str(v) for k, v in params.items()})


def format_name(pattern: str, params: Dict[str, Any]) -> str:
    return _substitute(pattern, params)


def render_template(text: str, params: Dict[str, Any]) -> str:
    return _substitute(text, params)


def load_suite(path) -> Suite:
    path = Path(path)
    data = load_toml(path)
    suite_table = data["suite"]
    cases: List[Case] = []
    for case_table in data.get("case", []):
        if "template" in case_table:
            base = dict(case_table.get("params", {}))
            sweep = case_table.get("sweep", {})
            param_sets = expand_sweep(sweep) if sweep else [{}]
            for sweep_params in param_sets:
                params = derive_params({**base, **sweep_params})
                cases.append(
                    Case(
                        name=format_name(case_table["name"], params),
                        kind="template",
                        source=path.parent / case_table["template"],
                        params=params,
                    )
                )
        elif "input" in case_table:
            cases.append(
                Case(
                    name=case_table["name"],
                    kind="input",
                    source=path.parent / case_table["input"],
                )
            )
        else:
            raise ValueError("each [[case]] must have either 'template' or 'input'")
    names = [c.name for c in cases]
    if len(names) != len(set(names)):
        raise ValueError("case names must be unique after sweep expansion")
    return Suite(
        name=suite_table["name"],
        repeat=int(suite_table.get("repeat", 1)),
        cases=cases,
    )
