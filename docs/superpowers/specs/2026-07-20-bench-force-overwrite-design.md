# bench.py run --force(既存ラベル上書き)設計

- 日付: 2026-07-20
- ブランチ: `benchmark_automation`
- ステータス: 承認済み(実装計画待ち)
- 前提: `docs/superpowers/specs/2026-07-17-benchmark-automation-design.md` のハーネス

## 目的

`bench.py run` は既存の `results/<label>` があるとエラー終了する。`--force` フラグで
既存ラベルを上書きして再実行できるようにする。

## 設計

### セマンティクス

`--force` は「ラベルディレクトリを丸ごと削除してから実行」。追記・マージはしない。
理由: スイート構成を変えて再実行した場合、マージでは前回のみに存在したケースが残留し、
compare が新旧混在の結果を比較してしまう。結果セットはラベル単位でアトミックに保つ。

削除対象は `results-dir/<label>` のみ。`--results-dir` 自体や他のラベルには触れない。

### CLI(benchmark/bench.py)

- `run` サブコマンドに `--force`(`action="store_true"`、
  help: "delete an existing label directory before running")を追加。
- `cmd_run` の既存チェックを変更:
  - 存在 + `--force` なし → `sys.exit("error: {} already exists; use a new label or --force")`
  - 存在 + `--force` あり → `shutil.rmtree(results_dir)` してから実行
- `import shutil` を bench.py に追加。

### テスト(新規 test/python/test_benchmark_cli.py)

`benchmark/` を sys.path に足して `bench` を import し、`argparse.Namespace` で
`cmd_run` を直接呼ぶ(スタブ実行ファイルは既存 runner テストのパターンを流用):

1. `--force` なし + 既存ラベルディレクトリ → `SystemExit`、メッセージに "already exists"。
2. `--force` あり + 古いマーカーファイル入りの既存ラベル → 正常終了し、マーカーが消え、
   新しい結果(meta.json とケースディレクトリ)が存在する。

### ドキュメント

`benchmark/README.md` に `--force` の一文を追加(既存ラベルの拒否と上書きの説明)。

## スコープ外

- 部分上書き(特定ケースのみ再実行)
- compare 側の変更
