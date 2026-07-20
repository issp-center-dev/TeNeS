# bench.py run --force 実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `bench.py run --force` で既存のラベルディレクトリを削除してから再実行できるようにする。

**Architecture:** `cmd_run` の既存チェックを分岐(`--force` あり→ `shutil.rmtree` して続行、なし→ エラーメッセージに `or --force` を追記して終了)。CLI 層のみの変更で runner/benchlib は触らない。テストは新規 `test_benchmark_cli.py` で `bench` モジュールを import して `cmd_run` を直接検証。

**Tech Stack:** Python 3 / pytest

**Spec:** `docs/superpowers/specs/2026-07-20-bench-force-overwrite-design.md`

## Global Constraints

- 削除対象は `results-dir/<label>` のみ(`--results-dir` 自体や他のラベルに触れない)。
- `--force` は「丸ごと削除してから実行」。マージはしない。
- Python は black でフォーマットしてからコミット。新規テストファイルの先頭には既存
  `test/python/test_tenes_std.py` と同一の 15 行 `#` コメント GPL ヘッダを付ける。
- pytest は repo ルートから `./venv/bin/python3 -m pytest ...`(bare python は `toml` 不足で
  test_tenes_simple の収集に失敗する)。
- コミットメッセージ末尾に必須フッター:

```
Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_0147iVGzythY7Y3iwtHobQxL
```

---

### Task 1: --force フラグ

**Files:**
- Modify: `benchmark/bench.py`(import 部・`cmd_run`・argparse)
- Create: `test/python/test_benchmark_cli.py`
- Modify: `benchmark/README.md`(Quick start の複数スイート段落の後)

**Interfaces:**
- Consumes: 既存 `bench.cmd_run(args)`(argparse.Namespace: `suite` list / `label` / `tenes_dir` / `tool_dir` / `results_dir` / `launcher` / `repeat`)、`runner.run_suites`
- Produces: `args.force: bool` を追加した `cmd_run`

- [ ] **Step 1: 失敗するテストを書く**

`test/python/test_benchmark_cli.py` を作成(GPL ヘッダは test_tenes_std.py からコピー):

```python
# GPLヘッダ(15行、test/python/test_tenes_std.py と同一)

import argparse
import os
import sys

import pytest

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmark"),
)

import bench


def _make_stub(path, body):
    path.write_text("#!/bin/sh\n" + body + "\n")
    path.chmod(0o755)


@pytest.fixture
def stub_tools(tmp_path):
    tool_dir = tmp_path / "tool"
    tool_dir.mkdir()
    _make_stub(tool_dir / "tenes_simple", "touch std.toml")
    _make_stub(tool_dir / "tenes_std", "touch input.toml")
    tenes_dir = tmp_path / "src"
    tenes_dir.mkdir()
    _make_stub(
        tenes_dir / "tenes",
        'mkdir -p output && echo "{}" > output/timers.json'
        ' && echo "0.5" > output/density.dat',
    )
    return tenes_dir, tool_dir


def _write_suite(tmp_path):
    (tmp_path / "tpl.toml").write_text("")
    suite_toml = tmp_path / "suite.toml"
    suite_toml.write_text(
        '[suite]\nname = "s"\n\n[[case]]\nname = "c"\ntemplate = "tpl.toml"\n'
    )
    return suite_toml


def _args(tmp_path, suite_toml, tenes_dir, tool_dir, force=False):
    return argparse.Namespace(
        suite=[str(suite_toml)],
        label="lbl",
        tenes_dir=str(tenes_dir),
        tool_dir=str(tool_dir),
        results_dir=str(tmp_path / "results"),
        launcher=None,
        repeat=None,
        force=force,
    )


def test_run_refuses_existing_label_without_force(tmp_path, stub_tools):
    tenes_dir, tool_dir = stub_tools
    suite_toml = _write_suite(tmp_path)
    (tmp_path / "results" / "lbl").mkdir(parents=True)
    with pytest.raises(SystemExit, match="already exists"):
        bench.cmd_run(_args(tmp_path, suite_toml, tenes_dir, tool_dir))


def test_run_force_replaces_existing_label(tmp_path, stub_tools):
    tenes_dir, tool_dir = stub_tools
    suite_toml = _write_suite(tmp_path)
    label_dir = tmp_path / "results" / "lbl"
    label_dir.mkdir(parents=True)
    stale = label_dir / "stale_case_marker"
    stale.write_text("old")
    bench.cmd_run(_args(tmp_path, suite_toml, tenes_dir, tool_dir, force=True))
    assert not stale.exists()
    assert (label_dir / "meta.json").exists()
    assert (label_dir / "c" / "run_0" / "timers.json").exists()
```

- [ ] **Step 2: 失敗を確認**

Run: `./venv/bin/python3 -m pytest test/python/test_benchmark_cli.py -v 2>&1 | tail -5`
Expected: `test_run_force_replaces_existing_label` が `AttributeError: 'Namespace' object has no attribute 'force'`… ではなく、`cmd_run` 内で `args.force` 未参照のため **SystemExit("already exists")** になり FAIL。`test_run_refuses_existing_label_without_force` は現行実装でも PASS する(メッセージ互換のため)。1 failed, 1 passed であること。

- [ ] **Step 3: 実装**

`benchmark/bench.py` の import 部に `import shutil` を追加(`import sys` の隣)。
`cmd_run` の既存チェックを次に置き換え:

```python
    results_dir = Path(args.results_dir) / args.label
    if results_dir.exists():
        if args.force:
            shutil.rmtree(results_dir)
        else:
            sys.exit(
                "error: {} already exists; use a new label or --force".format(
                    results_dir
                )
            )
```

argparse の `run` サブコマンドにフラグを追加(`--repeat` の行の後):

```python
    p_run.add_argument(
        "--force",
        action="store_true",
        help="delete an existing label directory before running",
    )
```

- [ ] **Step 4: テストが通ることを確認**

Run: `./venv/bin/python3 -m pytest test/python/test_benchmark_cli.py -v`
Expected: 2 passed

- [ ] **Step 5: README 更新**

`benchmark/README.md` の複数スイート説明(Quick start 内)のコードブロックの後に追加:

```markdown
`run` refuses to overwrite an existing label directory; pass `--force` to
delete it and rerun from scratch.
```

- [ ] **Step 6: E2E 検証と全体回帰**

```bash
python3 benchmark/bench.py run benchmark/suites/smoke.toml --label forcecheck \
    --tenes-dir out-gcc/build/src --tool-dir out-gcc/build/tool
python3 benchmark/bench.py run benchmark/suites/smoke.toml --label forcecheck \
    --tenes-dir out-gcc/build/src --tool-dir out-gcc/build/tool; echo "exit=$?"
python3 benchmark/bench.py run benchmark/suites/smoke.toml --label forcecheck --force \
    --tenes-dir out-gcc/build/src --tool-dir out-gcc/build/tool
rm -rf benchmark/results/forcecheck
./venv/bin/python3 -m pytest test/python -q
```

Expected: 1回目成功 → 2回目 `error: ... already exists; use a new label or --force` で exit=1 → 3回目(--force)成功。pytest 51 passed(49 + 新規2)。

- [ ] **Step 7: フォーマットとコミット**

```bash
python3 -m black benchmark/bench.py test/python/test_benchmark_cli.py
./venv/bin/python3 -m pytest test/python/test_benchmark_cli.py -q
git add benchmark/bench.py test/python/test_benchmark_cli.py benchmark/README.md
git commit -m "Add --force to bench.py run to overwrite an existing label"
```
