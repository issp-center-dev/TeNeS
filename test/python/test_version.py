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
import shutil
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


def commit_of_the_source_tree():
    """The commit --version must show for this source tree, or None where
    none can be known. Asked from git here, not from the code under test:
    with the commit optional in the expected output, a tool that had lost it
    would pass."""
    if os.path.exists(os.path.join(ROOT, ".git")):
        try:
            return subprocess.run(
                ["git", "-C", ROOT, "rev-parse", "HEAD"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return None
    try:
        with open(os.path.join(ROOT, "config", "git_archive.txt")) as f:
            line = f.readline().strip()
    except OSError:
        return None
    return line if re.fullmatch(r"[0-9a-f]+", line) else None


def expected_of_the_source_tree():
    """A pattern for the whole output of --version, without the newline."""
    version = re.escape(version_in_cmakelists())
    commit = commit_of_the_source_tree()
    if commit is None:
        return version
    return version + r" \(" + commit[:8] + r"(-dirty)?\)"


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

    def test_commit_is_the_one_of_the_source_tree(self, tool):
        assert re.fullmatch(expected_of_the_source_tree(), tool.version_string())

    def fake_git(self, monkeypatch, fail):
        """Make git answer rev-parse with HASH and status with a change,
        except for the subcommands in fail. Returns the calls made."""
        calls = []

        def run(cmd, **kwargs):
            calls.append(cmd)
            subcommand = [w for w in cmd[3:] if not w.startswith("-")][0]
            if subcommand in fail:
                raise subprocess.CalledProcessError(128, cmd)
            out = HASH + "\n" if subcommand == "rev-parse" else " M README.md\n"
            return subprocess.CompletedProcess(cmd, 0, stdout=out, stderr="")

        monkeypatch.setattr(subprocess, "run", run)
        return calls

    def checkout(self, monkeypatch, tool, tmp_path):
        source_tree(
            monkeypatch,
            tool,
            tmp_path,
            cmakelists="project(TeNeS CXX)\nset(TENES_VERSION 3.2.1)\n",
        )
        (tmp_path / ".git").mkdir()

    def test_checkout(self, monkeypatch, tool, tmp_path):
        self.checkout(monkeypatch, tool, tmp_path)
        calls = self.fake_git(monkeypatch, fail=())
        assert tool.version_string() == "3.2.1 (01234567-dirty)"
        # status must not take the lock of the index
        status = [c for c in calls if "status" in c]
        assert len(status) == 1 and "--no-optional-locks" in status[0]

    def test_status_that_fails_keeps_the_commit(self, monkeypatch, tool, tmp_path):
        self.checkout(monkeypatch, tool, tmp_path)
        self.fake_git(monkeypatch, fail=("status",))
        assert tool.version_string() == "3.2.1 (01234567)"

    def test_checkout_that_git_cannot_read(self, monkeypatch, tool, tmp_path):
        self.checkout(monkeypatch, tool, tmp_path)
        self.fake_git(monkeypatch, fail=("rev-parse", "status"))
        assert tool.version_string() == "3.2.1"

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
    assert re.fullmatch(expected_of_the_source_tree() + r"\n", result.stdout)


# config/write_version.cmake, the script that fills the version number and
# the commit into version.cpp and into the built tools.
needs_cmake_and_checkout = pytest.mark.skipif(
    shutil.which("cmake") is None or commit_of_the_source_tree() is None
    or not os.path.exists(os.path.join(ROOT, ".git")),
    reason="needs cmake and a git checkout that git can read",
)


def write_version(tmp_path, template, version, git_executable=None):
    (tmp_path / "in").write_text(template)
    cmd = ["cmake", "-DSOURCE_DIR=" + ROOT, "-DTENES_VERSION=" + version]
    cmd += ["-DINPUT=" + str(tmp_path / "in"), "-DOUTPUT=" + str(tmp_path / "out")]
    if git_executable is not None:
        cmd.append("-DGIT_EXECUTABLE=" + git_executable)
    cmd += ["-P", os.path.join(ROOT, "config", "write_version.cmake")]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return (tmp_path / "out").read_text()


TEMPLATE = "{} @TENES_VERSION@ @TENES_GIT_HASH@ @TENES_GIT_DIRTY@\n"


@needs_cmake_and_checkout
def test_write_version_fills_in_the_commit(tmp_path):
    words = write_version(tmp_path, TEMPLATE.format("first"), "3.2.1").split()
    assert words[:3] == ["first", "3.2.1", commit_of_the_source_tree()]
    assert words[3] in ("true", "false")


@needs_cmake_and_checkout
def test_write_version_without_git_keeps_the_commit_only(tmp_path):
    # A build that cannot ask git fills in the commit of the last build that
    # could, and nothing else of it: the template and the version number are
    # those of this build.
    first = write_version(tmp_path, TEMPLATE.format("first"), "3.2.1").split()
    no_git = str(tmp_path / "no_such_git")
    second = write_version(
        tmp_path, TEMPLATE.format("second"), "3.2.2", git_executable=no_git
    ).split()
    assert second == ["second", "3.2.2"] + first[2:]


@needs_cmake_and_checkout
def test_write_version_without_git_from_the_start(tmp_path):
    no_git = str(tmp_path / "no_such_git")
    out = write_version(
        tmp_path, "@TENES_VERSION@|@TENES_GIT_HASH@|@TENES_GIT_DIRTY@\n", "3.2.1",
        git_executable=no_git,
    )
    assert out == "3.2.1||false\n"
