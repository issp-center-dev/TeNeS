# bench.py 複数スイート順次実行 設計

- 日付: 2026-07-20
- ブランチ: `benchmark_automation`
- ステータス: 承認済み(実装計画待ち)
- 前提: `docs/superpowers/specs/2026-07-17-benchmark-automation-design.md` のハーネス

## 目的

`bench.py run` で複数のスイート定義を1コマンド・同一ラベルで順次実行できるようにする。

```bash
python3 benchmark/bench.py run suites/contraction.toml suites/e2e.toml --label baseline
```

## 設計

### CLI(benchmark/bench.py)

- `run` サブコマンドの `suite` 引数を `nargs="+"` に変更(1個以上、引数順に実行)。
- `cmd_run` は全スイートを先にロードし、`runner.run_suites(suites, ctx)` を呼ぶ。
- 既存の「ラベルディレクトリが存在したら拒否」の挙動は不変。
- `compare` サブコマンドは無変更。

### 結果レイアウト(変更なし=フラット統合)

全スイートのケースを現行と同じ `results/<label>/<case>/run_<i>/` に置く。
`compare` / `load_label_dir` は無変更で動作する。

### runner(benchmark/benchlib/runner.py)

- 新規 `run_suites(suites, ctx)`:
  1. 全スイートのケース名を横断で重複チェックし、衝突があれば **何も実行せず**
     `ValueError`(衝突名を列挙)。
  2. `meta.json` を1回書く。`suite` キーは各スイート名のカンマ結合
     (例: `"contraction,e2e"`)。
  3. スイートを引数順に実行。`repeat` は各スイート定義の値を使う
     (`ctx.repeat` 指定時は全スイートに適用)。
- 既存 `run_suite(s, ctx)` は互換のため残し、`run_suites([s], ctx)` へ委譲する。

### エラー処理

ケース失敗時は現行どおり即停止(`_call` の `check=True` / timers.json 欠損の
`RuntimeError`)。keep-going 機能は追加しない(既存フォローアップ項目のまま)。

### テスト(test/python/test_benchmark_runner.py)

pytest を3件追加(既存のスタブ実行ファイルのパターンを流用):

1. 2スイートを `run_suites` で実行 → 両方のケースディレクトリが同一ラベル下にでき、
   `meta.json` の `suite` が `"s1,s2"` になる。
2. per-suite repeat の尊重 → repeat=2 のスイートは `run_0`/`run_1`、repeat=1 の
   スイートは `run_0` のみ。
3. スイート横断のケース名重複 → `ValueError`、ケースディレクトリが1つも作られない。

### ドキュメント

`benchmark/README.md` の Quick start に複数スイート指定の例を1行追加。

## スコープ外

- ケース失敗時の keep-going(継続実行)
- スイートごとのサブディレクトリ分割(案Bとして検討し却下: compare の変更が
  必要になる割に、比較はケース名マッチングなので実益がない)
