# 外部 OSS ライブラリによる TeNeS 高度化の調査メモ

- 調査日: 2026-07-17(メモ化: 2026-07-20)
- 前提: 性能評価基盤としてベンチマークハーネスを整備済み(PR #110、develop マージ済み)。
  性能に関わる変更は `benchmark/bench.py run` / `compare` で A/B 比較してから採否を決める
  (使い方は `benchmark/README.md`)。

## インパクト大(研究機能・性能に直結)

### 1. GPU 対応 — NVIDIA cuQuantum (cuTensorNet)

mptensor は CPU(LAPACK/ScaLAPACK)のみで、CTMRG の縮約・SVD が全コスト。
cuTensorNet は縮約順序最適化+GPU 実行を提供し、iPEPS 系では variPEPS などの実績がある。
`src/tensor.hpp` が既にバックエンド切替の抽象点(`scalapack::Matrix` vs `lapack::Matrix`)に
なっているので第3のバックエンドとして追加する設計は可能だが、mptensor API
(transpose/reshape/svd)への依存が深く大工事。着手するならまず単一 GPU・`_NO_MPI` 相当から。
規模的に独立プロジェクト級。

### 2. 縮約順序最適化 — cotengra / opt_einsum

`misc/contraction/netcon.py` + `tdt.py` で netcon アルゴリズムを自前実装し、
`src/iTPS/core/*/ctm_NxM.cpp` のコード生成に使っている。
[cotengra](https://github.com/jcmgray/cotengra) はハイパーグラフ分割ベースで大規模
ネットワーク(4x4 超のクラスタや multisite 観測量)にもスケールし、スライシングによる
メモリ削減候補も出せる。**オフラインのコード生成ツールなので本体に依存が入らず低リスク**。
効果は `contract/<backend>/<N>x<M>` タイマーで直接測定できる。着手優先度: 最上位。

### 3. Arnoldi の置き換え — ARPACK-NG(Spectra は不向き)

`src/arnoldi.cpp` の自前 implicit restart Arnoldi は、収束判定やシフト戦略で実績ある
ライブラリに劣る可能性がある。ARPACK-NG は reverse communication インターフェースなので、
現在の `std::function` による matvec 抽象(MPI 分散テンソル対応)をそのまま接続できる。
Spectra はヘッダオンリーで導入は楽だが Eigen ベクトル前提のため分散テンソルと相性が悪い。
相関長計算の信頼性・速度向上が見込める。

### 4. 対称性保存テンソル(長期ロードマップ)

U(1)/Z2 対称性のブロックスパーステンソルは同じ計算コストでボンド次元を数倍にできる。
ITensor(C++/Julia)、Cytnx、TensorKit.jl などがあるが、mptensor 全置換になるため
「高度化要素」ではなくロードマップ案件。mptensor 側に対称性を足すか、バックエンド抽象を
広げるかの設計判断が先。

## 実用性向上(中コスト)

### 5. チェックポイント・出力の HDF5 化 — HighFive

save/load は mptensor の生バイナリ(MPI ランク数に依存)、観測値は `.dat` テキスト。
HighFive(ヘッダオンリー)+並列 HDF5 にすると、ランク数を変えた restart、メタデータ同梱、
Python からの直接読み出し(h5py/pandas)が可能になる。既存フォーマットとの互換層が必要。

### 6. Python バインディング — nanobind

`libtenes`(共有ライブラリ)があるので、nanobind で 3 段パイプラインを Python から
一気通貫で回せるようにできる。パラメータスキャンやノートブック利用に有効。
**注意**: `TimerRegistry` はプロセス単位のシングルトンで、同一プロセス内で
`tenes_itps_main` を複数回呼ぶと計時が二重加算される(PR #110 の最終レビューで指摘済み)。
バインディング実装時にはレジストリのリセット機構が必要。

### 7. Python ツールの近代化

- `toml` パッケージ(メンテ停滞)→ 読み込みは `tomllib`(3.11+ 標準)/`tomli`、
  書き出しは `tomli-w`。benchmark/ ハーネスは既に tomllib+フォールバック方式。
- 入力検証に pydantic を入れると `simple.toml` のエラーメッセージが大幅に親切になる
  (C++ 側の toml11 エラー出力強化と同じ方向性)。

## QoL(低コスト)

- **{fmt}** — 出力整形。C++17 なので `std::format` は不可、fmt はヘッダオンリーモードあり。
- **spdlog** — `src/printlevel.hpp` の自前ログレベル管理の置き換え。
- **CLI11** — `src/main.cpp` の引数解析。

## 推奨着手順

効果/コスト比で **cotengra(2)→ ARPACK-NG(3)→ HDF5(5)**。GPU(1)は別枠プロジェクト。
いずれも導入前後で `benchmark/suites/contraction.toml`(カーネル粒度)と `e2e.toml`
(実ワークロード)の A/B 比較を取ること。

## 関連資料

- ベンチマークハーネス設計: `docs/superpowers/specs/2026-07-17-benchmark-automation-design.md`
- ハーネス使用法: `benchmark/README.md`
- 制約メモ: square lattice の unit cell は `tenes_simple` の制約で最小 2x1(1x1 不可)
