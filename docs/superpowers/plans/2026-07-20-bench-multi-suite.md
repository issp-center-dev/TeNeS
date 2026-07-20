# bench.py 複数スイート順次実行 実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `bench.py run` に複数のスイート TOML を渡し、同一ラベルの下で順次実行できるようにする。

**Architecture:** `runner.run_suites(suites, ctx)` を追加(スイート横断の重複名チェック → meta.json 1回書き出し → 順次実行)。既存 `run_suite` は `run_suites([s], ctx)` へ委譲。CLI は `suite` 引数を `nargs="+"` にするだけ。結果レイアウトはフラットのまま(compare 無変更)。

**Tech Stack:** Python 3(既存 benchlib)/ pytest

**Spec:** `docs/superpowers/specs/2026-07-20-bench-multi-suite-design.md`

## Global Constraints

- Python は black(line-length 88)でフォーマットしてからコミット。
- 重複ケース名の `ValueError` は **results ディレクトリ作成・meta.json 書き出しより前**に送出する(ラベルディレクトリを作り残さない — bench.py の「存在したら拒否」を壊さないため)。
- `meta.json` の `suite` キーはスイート名のカンマ結合(例: `"s1,s2"`)。単一スイート時は従来と同値。
- `repeat` は各スイート定義の値を使い、`ctx.repeat` 指定時は全スイートに適用。
- コミットメッセージ末尾に必須フッター:

```
Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_0147iVGzythY7Y3iwtHobQxL
```

- テスト実行: pytest は repo ルートから `python3 -m pytest ...`(bare python で `import toml` が失敗する場合は `./venv/bin/python3`)。設定済みビルドは `out-gcc/build`。

---

### Task 1: runner.run_suites

**Files:**
- Modify: `benchmark/benchlib/runner.py:83-94`(`run_suite` を置き換え)
- Test: `test/python/test_benchmark_runner.py`(3テスト追加)

**Interfaces:**
- Consumes: 既存 `run_case(case, repeat, ctx)`、`meta_mod.collect_meta`、`suite.Suite`(`name`/`repeat`/`cases`)、テストの既存フィクスチャ `stub_tools`(`(tenes_dir, tool_dir)` を返す)
- Produces: `run_suites(suites: list[Suite], ctx: RunContext) -> None`(Task 2 の CLI が呼ぶ)。`run_suite(s, ctx)` は互換維持。

- [ ] **Step 1: 失敗するテストを書く**

`test/python/test_benchmark_runner.py` の末尾に追加(`json`, `pytest`, `runner`, `suite` は import 済み):

```python
def _stub_suite(tmp_path, name, case_names, repeat=1):
    tpl = tmp_path / "tpl_shared.toml"
    if not tpl.exists():
        tpl.write_text("")
    return suite.Suite(
        name=name,
        repeat=repeat,
        cases=[
            suite.Case(name=c, kind="template", source=tpl, params={})
            for c in case_names
        ],
    )


def test_run_suites_merges_cases_under_one_label(tmp_path, stub_tools):
    tenes_dir, tool_dir = stub_tools
    s1 = _stub_suite(tmp_path, "s1", ["a"])
    s2 = _stub_suite(tmp_path, "s2", ["b"])
    ctx = runner.RunContext(
        tenes_dir=tenes_dir, tool_dir=tool_dir, results_dir=tmp_path / "results"
    )
    runner.run_suites([s1, s2], ctx)
    meta_data = json.loads((tmp_path / "results" / "meta.json").read_text())
    assert meta_data["suite"] == "s1,s2"
    assert (tmp_path / "results" / "a" / "run_0" / "timers.json").exists()
    assert (tmp_path / "results" / "b" / "run_0" / "timers.json").exists()


def test_run_suites_honors_per_suite_repeat(tmp_path, stub_tools):
    tenes_dir, tool_dir = stub_tools
    s1 = _stub_suite(tmp_path, "s1", ["a"], repeat=2)
    s2 = _stub_suite(tmp_path, "s2", ["b"], repeat=1)
    ctx = runner.RunContext(
        tenes_dir=tenes_dir, tool_dir=tool_dir, results_dir=tmp_path / "results"
    )
    runner.run_suites([s1, s2], ctx)
    assert (tmp_path / "results" / "a" / "run_1").exists()
    assert not (tmp_path / "results" / "b" / "run_1").exists()


def test_run_suites_rejects_cross_suite_duplicate_names(tmp_path, stub_tools):
    tenes_dir, tool_dir = stub_tools
    s1 = _stub_suite(tmp_path, "s1", ["same"])
    s2 = _stub_suite(tmp_path, "s2", ["same"])
    ctx = runner.RunContext(
        tenes_dir=tenes_dir, tool_dir=tool_dir, results_dir=tmp_path / "results"
    )
    with pytest.raises(ValueError, match="same"):
        runner.run_suites([s1, s2], ctx)
    assert not (tmp_path / "results").exists()
```

- [ ] **Step 2: 失敗を確認**

Run: `python3 -m pytest test/python/test_benchmark_runner.py -v 2>&1 | tail -5`
Expected: 3件が `AttributeError: module 'benchlib.runner' has no attribute 'run_suites'` で FAIL。

- [ ] **Step 3: 実装**

`benchmark/benchlib/runner.py` の `run_suite`(83-94行)を次で置き換える:

```python
def run_suites(suites, ctx):
    """Run multiple suites sequentially under one label directory."""
    names = [case.name for s in suites for case in s.cases]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise ValueError(
            "case names collide across suites: {}".format(", ".join(duplicates))
        )
    ctx.results_dir.mkdir(parents=True, exist_ok=True)
    repo_dir = Path(__file__).resolve().parents[2]
    result = meta_mod.collect_meta(
        ctx.tenes_dir / "tenes", repo_dir=repo_dir, launcher=ctx.launcher
    )
    result["suite"] = ",".join(s.name for s in suites)
    (ctx.results_dir / "meta.json").write_text(json.dumps(result, indent=2))
    for s in suites:
        repeat = ctx.repeat if ctx.repeat is not None else s.repeat
        for case in s.cases:
            print(
                "== case: {} (repeat={}) ==".format(case.name, repeat), flush=True
            )
            run_case(case, repeat, ctx)


def run_suite(s, ctx):
    run_suites([s], ctx)
```

- [ ] **Step 4: テストが通ることを確認**

Run: `python3 -m pytest test/python/test_benchmark_runner.py -v`
Expected: 9 passed(既存6 + 新規3)

- [ ] **Step 5: フォーマットとコミット**

```bash
python3 -m black benchmark/benchlib/runner.py test/python/test_benchmark_runner.py
python3 -m pytest test/python/test_benchmark_runner.py -q
git add benchmark/benchlib/runner.py test/python/test_benchmark_runner.py
git commit -m "Add run_suites for sequential multi-suite execution"
```

---

### Task 2: CLI と README

**Files:**
- Modify: `benchmark/bench.py:39-52`(`cmd_run`)と `benchmark/bench.py:70-71`(argparse)
- Modify: `benchmark/README.md`(Quick start セクション)

**Interfaces:**
- Consumes: Task 1 の `runner.run_suites(suites, ctx)`、既存 `suite.load_suite(path)`

- [ ] **Step 1: bench.py を変更**

`cmd_run`(39行目付近)を次に置き換え:

```python
def cmd_run(args):
    suites = [suite.load_suite(path) for path in args.suite]
    results_dir = Path(args.results_dir) / args.label
    if results_dir.exists():
        sys.exit("error: {} already exists; use a new label".format(results_dir))
    ctx = runner.RunContext(
        tenes_dir=args.tenes_dir,
        tool_dir=args.tool_dir,
        results_dir=results_dir,
        launcher=args.launcher,
        repeat=args.repeat,
    )
    try:
        runner.run_suites(suites, ctx)
    except ValueError as e:
        sys.exit("error: {}".format(e))
    print("results saved to {}".format(results_dir))
```

argparse の suite 引数(71行目)を変更:

```python
    p_run.add_argument(
        "suite", nargs="+", help="path(s) to suite TOML files (run sequentially)"
    )
```

`run` サブコマンドの help 文字列(70行目)を `"run one or more benchmark suites"` に変更。

- [ ] **Step 2: README を更新**

`benchmark/README.md` の Quick start セクション(step 1 のコードブロックの直後)に追加:

```markdown
Multiple suite files can be given in one invocation; they run sequentially
under the same label (case names must be unique across the suites):

```bash
python3 benchmark/bench.py run \
    benchmark/suites/contraction.toml benchmark/suites/e2e.toml \
    --label baseline --tenes-dir build/src --tool-dir build/tool
```
```

(既存の見出し・文体に合わせて位置は微調整してよい。内容は上記のまま。)

- [ ] **Step 3: E2E 検証**

```bash
python3 benchmark/bench.py run benchmark/suites/smoke.toml benchmark/suites/ci.toml \
    --label multicheck --tenes-dir out-gcc/build/src --tool-dir out-gcc/build/tool
python3 -c "
import json
m = json.load(open('benchmark/results/multicheck/meta.json'))
assert m['suite'] == 'smoke,ci', m['suite']
print('meta OK')
"
python3 benchmark/bench.py compare benchmark/results/multicheck benchmark/results/multicheck | head -20
python3 benchmark/bench.py run benchmark/suites/smoke.toml benchmark/suites/smoke.toml \
    --label dupcheck --tenes-dir out-gcc/build/src --tool-dir out-gcc/build/tool; echo "exit=$?"
rm -rf benchmark/results/multicheck benchmark/results/dupcheck
```

Expected: multicheck は3ケース(smoke 1 + ci 2)完走、`meta OK`、self-compare に WARNING なし。dupcheck は traceback ではなく `error: case names collide across suites: smoke_2x1_D2_chi8` で exit=1、`benchmark/results/dupcheck` は作られない(rm は空振りでよい)。

- [ ] **Step 4: 全体回帰**

```bash
python3 -m pytest test/python -q
(cd out-gcc/build && ctest -R benchmark_smoke --output-on-failure)
```

Expected: pytest 49 passed(46 + Task 1 の 3)、smoke ctest PASS(単一スイート経路の互換確認)。

- [ ] **Step 5: フォーマットとコミット**

```bash
python3 -m black benchmark/bench.py
git add benchmark/bench.py benchmark/README.md
git commit -m "Accept multiple suite files in bench.py run"
```
