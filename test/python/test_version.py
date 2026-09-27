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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/.

"""--version of tenes_simple and tenes_std: the version number written in
the top-level CMakeLists.txt, followed by the commit."""

import os
import re
import subprocess
import sys

import pytest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tool"))

import tenes_simple
import tenes_std

HASH = "0123456789abcdef0123456789abcdef01234567"

tools = pytest.mark.parametrize("tool", [tenes_simple, tenes_std])


def version_in_cmakelists():
    with open(os.path.join(ROOT, "CMakeLists.txt")) as f:
        return re.search(r"set\(TENES_VERSION (\S+)\)", f.read()).group(1)


def as_built(monkeypatch, tool, version, githash, dirty):
    monkeypatch.setattr(tool, "_TENES_VERSION", version)
    monkeypatch.setattr(tool, "_TENES_GIT_HASH", githash)
    monkeypatch.setattr(tool, "_TENES_GIT_DIRTY", dirty)


def source_tree(monkeypatch, tool, tmp_path, cmakelists=None, archive=None):
    """Put the script into tmp_path/tool, as if it were run from there."""
    (tmp_path / "tool").mkdir()
    if cmakelists is not None:
        (tmp_path / "CMakeLists.txt").write_text(cmakelists)
    if archive is not None:
        (tmp_path / "config").mkdir()
        (tmp_path / "config" / "git_archive.txt").write_text(archive)
    monkeypatch.setattr(tool, "__file__", str(tmp_path / "tool" / "script.py"))


@tools
class TestBuiltByCMake:
    def test_commit_is_abbreviated_to_eight_digits(self, monkeypatch, tool):
        as_built(monkeypatch, tool, "3.2.1", HASH, "false")
        assert tool.version_string() == "3.2.1 (01234567)"

    def test_uncommitted_changes_are_marked(self, monkeypatch, tool):
        as_built(monkeypatch, tool, "3.3-dev", HASH, "true")
        assert tool.version_string() == "3.3-dev (01234567-dirty)"

    def test_commit_not_known(self, monkeypatch, tool):
        as_built(monkeypatch, tool, "3.2.1", "", "false")
        assert tool.version_string() == "3.2.1"


@tools
class TestRunFromSourceTree:
    def test_version_number_is_the_one_of_cmakelists(self, tool):
        version = tool.version_string().split()[0]
        assert version == version_in_cmakelists()

    def test_format(self, tool):
        assert re.fullmatch(
            r"\S+( \([0-9a-f]{8}(-dirty)?\))?", tool.version_string()
        )

    def test_tarball(self, monkeypatch, tool, tmp_path):
        source_tree(
            monkeypatch,
            tool,
            tmp_path,
            cmakelists="project(TeNeS CXX)\nset(TENES_VERSION 3.2.1)\n",
            archive=HASH + "\n",
        )
        assert tool.version_string() == "3.2.1 (01234567)"

    def test_copy_that_git_archive_did_not_make(self, monkeypatch, tool, tmp_path):
        source_tree(
            monkeypatch,
            tool,
            tmp_path,
            cmakelists="project(TeNeS CXX)\nset(TENES_VERSION 3.2.1)\n",
            archive="$Format:%H$\n",
        )
        assert tool.version_string() == "3.2.1"

    def test_script_outside_the_source_tree(self, monkeypatch, tool, tmp_path):
        # the directory above is a git repository of something else
        subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
        source_tree(monkeypatch, tool, tmp_path)
        assert tool.version_string() == "unknown"


@pytest.mark.parametrize("script", ["tenes_simple.py", "tenes_std.py"])
def test_version_option(script):
    result = subprocess.run(
        [sys.executable, os.path.join(ROOT, "tool", script), "--version"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert re.fullmatch(
        version_in_cmakelists().replace(".", r"\.") + r"( \([0-9a-f]{8}(-dirty)?\))?\n",
        result.stdout,
    )
