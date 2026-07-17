# パフォーマンス計測自動化(ベンチマークハーネス)設計

- 日付: 2026-07-17
- ブランチ: `benchmark_automation`
- ステータス: 承認済み(実装計画待ち)

## 背景と目的

cotengra による縮約順序最適化など、性能に関わる変更を TeNeS に導入する前提として、
変更前後の性能を定量的に比較(A/B 比較)できる自動化された仕組みを整備する。

- **主目的**: 同一マシン上で 2 つのビルド/設定を走らせて性能差を定量化する A/B 比較。
  CI での継続的リグレッション検出は主目的としない(将来、同じハーネスを縮小スイートで
  GitHub Actions に載せる余地は残す)。
- **計測粒度**: 縮約カーネル単位(`Contract_NxM` 等の呼び出し回数・累積時間)まで。
  ただし縮約カーネルは最初の適用例にすぎず、任意の計測点(CTM move、SVD、フェーズ、
  E2E 全体)に一般化できるデータモデルとする。
- **実行環境**: ローカル mac(`_NO_MPI`)、HPC(MPI 並列)、GitHub Actions(参考値)の
  3 環境すべてで同じ仕組みが動くこと。
- **ケース**: unit cell サイズ × ボンド次元 D × 環境ボンド次元 χ をスイープする専用
  スイートを主軸に、実ワークロード代表として既存の `sample/` ベースのケースも含める。

## 検討した代替案

- **Google Benchmark によるカーネル直叩きマイクロベンチ**: 統計処理と A/B 比較ツールが
  既製で手に入るが、MPI 分散テンソルと相性が悪く実質ローカル専用。エンドツーエンドの
  実テンソル形状を反映しない。→ 不採用(将来、カーネル深掘りが必要になったら
  C++ 計測機構を流用して後付け可能)。
- **ASV / github-action-benchmark 等の継続的トラッキング基盤**: 主目的が A/B 比較なので
  オーバーキル。C++ バイナリ + MPI との統合も困難。→ 不採用。

## 設計

### 1. 計測データモデルと C++ 計測機構

計測対象を「階層的な名前を持つ累積タイマーの集合」として一般化する。集計・比較ツールは
名前でマッチングするだけで、特定の計測点の知識を持たない。

**名前空間の規約**(スラッシュ区切り):

```
total                                  # E2E(プロセス全体)
phase/simple_update                    # 既存 time.dat 相当のフェーズ
phase/environment
phase/observable
contract/itps_ctm/2x2                  # 縮約カーネル単位
contract/density_ctm/3x3
ctm/left_move                          # 将来追加する任意の計測点
```

**実装**:

- `TimerRegistry`(`src/timer.hpp` を拡張): `名前 → {呼び出し回数, 累積秒}` のマップ。
  RAII の `ScopedTimer`(コンストラクタで開始、デストラクタで加算)を計測点に置く。
  既存 `Timer` クラスを内部で利用する。
- **常時有効**: オーバーヘッドは縮約 1 回(ms〜秒オーダー)につき chrono 呼び出し 2 回+
  マップ加算 1 回で無視できるため、フラグ分岐は入れない。
- **スレッド安全性**: 計測点は逐次コンテキストに置き、並列性はその内側(BLAS/mptensor)に
  ある構造を原則とする。実装時に呼び出しコンテキストを確認し、OpenMP 並列領域内から
  呼ばれる計測点が必要になった場合のみ per-thread 集計を検討する。
- **MPI 集計**: mptensor の演算は全 rank 集団実行なので呼び出し回数は全 rank で一致する。
  終了時に `MPI_Reduce` で rank 間の max/min/mean を集計し、rank 0 が出力する
  (負荷不均衡の可視化を兼ねる)。
- **出力**: 実行終了時に `output/timers.json` を書き出す。JSON は手書き生成し外部依存を
  追加しない。既存 `time.dat` は互換性のため残す。

```json
{
  "meta": {"tenes_version": "2.2-dev", "mpi_size": 4, "omp_threads": 8},
  "timers": {
    "total": {"count": 1, "sum": 123.4, "max_rank": 125.0, "min_rank": 122.1},
    "contract/itps_ctm/2x2": {"count": 4800, "sum": 45.6, "max_rank": 46.0, "min_rank": 44.9}
  }
}
```

**初期の計測点**: 既存フェーズ計時のレジストリ統合(`phase/*`)、縮約カーネルの
ディスパッチ点(`contract/*`)、`total`。CTM move 単位や SVD などは必要になった時点で追加。

### 2. ベンチマークハーネス(`benchmark/`)

```
benchmark/
  bench.py                 # CLI エントリポイント(run / compare)
  suites/
    contraction.toml       # 専用スイート: unit cell × D × χ スイープ
    e2e.toml               # 実ワークロード代表(既存 sample/ ベース)
  templates/               # tenes_simple 用 toml テンプレート
  results/                 # 計測結果(.gitignore 対象)
    <label>/
      meta.json            # git commit, コンパイラ, ホスト, 日時, スレッド/rank 数
      <case_name>/
        run_0/ run_1/ ...  # 反復ごとの timers.json + 主要観測量
```

**スイート定義**(TOML。Python 3.11+ は tomllib、それ未満は toml パッケージに
フォールバック):

```toml
[suite]
name = "contraction"
repeat = 3                      # 反復回数(比較時に中央値を取る)

[[case]]
name = "afh_square_{Lsub}_D{D}_chi{chi}"
template = "templates/afh_square.toml"   # tenes_simple 用テンプレート
sweep = { Lsub = [[1,1],[2,2],[3,3],[4,4]], D = [2,4], chi_ratio = [1,2] }

[[case]]
name = "honeycomb_e2e"
input = "../test/data/Honeycomb.toml"   # 既存入力をそのまま使う
```

- テンプレートは `simple.toml` へのパラメータ置換で、ハーネスが
  `tenes_simple → tenes_std → tenes` のパイプラインを流す(fulltest.py と同じ流儀)。
  前段 2 つは計測対象外(所要時間の記録のみ)。
- ベンチ用ケースはステップ数を絞り 1 ケース数十秒〜数分に収める。スイープの直積が
  大きくなりすぎる場合に備え、スイートファイル側で明示列挙も可能にする。

**実行**:

```bash
python3 benchmark/bench.py run suites/contraction.toml \
    --label baseline \
    --tenes-dir build/src \
    --launcher "mpirun -np 4"     # HPC では "srun -n 64" など。省略時は直接実行
```

- `--label` が A/B の識別子(例: `baseline` / `cotengra`)。
- `meta.json` に記録するもの: git コミットハッシュ(`git rev-parse HEAD`)と dirty フラグ
  (`git status --porcelain` が空か否か)、`tenes --version` 出力、ホスト名、
  `OMP_NUM_THREADS`、launcher 文字列、日時。
- 反復は同一ケースを連続実行し、生の `timers.json` をそのまま保存する
  (集計は比較時に行う。生データを捨てない)。
- 各 run の `output/` 以下の観測量ファイル(`*.dat`)を timers.json と併せてコピーして
  保存する(物理量一致チェック用)。
- ジョブスクリプト生成機能は作らない(YAGNI)。バッチジョブ内から `bench.py run` を
  呼べればよい。

### 3. 集計・比較(`bench.py compare`)

```bash
python3 benchmark/bench.py compare results/baseline results/cotengra
```

2 つの結果ディレクトリを取り、**ケース × タイマー名**でマッチングして Markdown の
比較レポートを出力する。

- **統計**: 反復の**中央値**を代表値とし、min–max の幅を併記。速度比は
  `median(B) / median(A)`。min–max の幅が重なる場合は「有意差なし」とマークする。
  検定は入れない(反復 3 回程度では意味がない)。
- **粒度混在への対処**: レポートは名前の階層でセクション分けする。まず `total` と
  `phase/*` のサマリ(E2E 視点)、次に `contract/*` などの詳細(カーネル視点)。
  プレフィックスごとに表を分けるだけなので、将来 `ctm/*` や `svd/*` を追加しても
  自動的に新セクションとして現れる。
- **片側にしか存在しない名前**(計測点追加後の比較、縮約パス構造の変化)はエラーに
  せず「A のみ / B のみ」欄に列挙する。比較ツールが名前の集合に何も仮定しないことが
  汎用性の要。
- **count の比較**: 呼び出し回数が変わっている場合はアルゴリズム変更を意味し、時間比較の
  前提が違うのでレポートで明示する。
- **物理量一致チェック**: 両ラベルの観測量を `np.isclose`(初期値は fulltest.py と同じ
  rtol=1e-3, atol=1e-4)で照合し、不一致があればレポート冒頭に警告を出す。
  「速くなったが答えが変わった」ことを見逃さないための安全弁。

### 4. 環境ごとの運用

- **ローカル mac**: `_NO_MPI` ビルドで直接実行。A/B は 2 つのビルドディレクトリを
  `--tenes-dir` で切り替えるか、コミット切替+再ビルドで運用。ハーネスは
  `OMP_NUM_THREADS` を meta に記録するのみで、マシン状態の制御には関与しない。
- **HPC**: ジョブスクリプト内から `--launcher "srun -n 64"` 等で実行。A/B の 2 ラベルを
  同一ジョブ内で連続実行することを推奨手順として README に明記する(ノード配置差の
  ノイズ回避)。
- **GitHub Actions**: 手動トリガー(`workflow_dispatch`)の専用ワークフローを用意。
  PR ごとの自動実行はしない。縮小スイートで「桁が変わる劣化」の検出用と位置づけ、
  結果 JSON は artifact としてアップロードし比較はローカルで行う。

### 5. テスト

- **C++**: `TimerRegistry` の doctest 単体テスト(加算・回数・JSON 出力形式)。
  既存の test_input 等と同列に ctest 登録。
- **Python**: スイート展開(スイープ直積、名前生成)と compare の統計・マッチング
  ロジックに pytest 単体テストを付け、ctest に組み込む。
- **統合**: 最小ケース 1 個 × 1 反復のスモークテストを ctest に追加(数十秒以内)。

### 6. ドキュメント

`benchmark/README.md` に使い方(A/B 比較の標準手順、HPC での注意)を書く。
Sphinx ドキュメントへの統合は当面しない。

## スコープ外

- PR ごとの自動ベンチマーク実行と履歴トラッキング
- ジョブスクリプト生成
- Google Benchmark によるマイクロベンチ(将来の拡張候補)
- cotengra の導入そのもの(このハーネス完成後の別プロジェクト)
