# フェルミオン2サイト blob の省メモリ化 — 補助資料

`src/fermion/reduced.hpp` の `fuse_doubled_cluster_direct` /
`build_reduced_pair_direct`（コミット `37a66c82`）の設計記録と、
設計段階で使った検証プローブ。実装のコメントからここを参照している。

## 背景（何を直したか）

フェルミオン模式の 2 サイト観測量（CTM 経路）は、blob を作るのに
bra 層と ket 層の **rank-16 全外積**（D¹²d⁴ 要素）を実体化していた。
Hubbard d=4, D=4 で 1 コピー 34 GB、転置と再分配で 4〜5 コピー同時に生きて
140〜180 GB に達し、実機（PBS、np=4、mem=240gb）で cgroup OOM kill となった。

採用した解が **direct fuse**: 外積 + 符号付き転置 + 物理脚の δ trace は、
符号項を掛ける先ごとに仕分ければ「物理脚の直接縮約（rank-12 = blob と同サイズ）」
と同値になる、という恒等式に基づく。新しい符号規約は導入していない。
詳細は `design.md` §3（改訂2）。

実測（2×2 Hubbard、χ=40、serial）:

| | 旧 | 新 |
|---|---|---|
| D=3 2サイト測定 | 521 s / peak RSS 4.85 GB | 24.9 s / 1.0 GB |
| D=4 2サイト測定 | 実行不能（実機で OOM） | 200 s / 4.4 GB |

## ファイル

| | |
|---|---|
| `design.md` | 設計書。改訂2 が実装された方式（改訂1 の当初案は §3 冒頭の経緯参照） |
| `contract.md` | `test/fermion/impurity_blob.cpp` の振る舞い契約書（追補3 まで） |
| `probe_direct.cpp` | **direct fuse の恒等式の検証**。実装より短く読める参照実装 |
| `probe_join.cpp` | 当初案（演算子の graded QR 分解 + サイト局所二重化）が対角 gauge で閉じないことを示した実験。ピボットの根拠 |
| `probe_split.cpp` | graded QR 分割の厳密再結合と k パリティ序列の検証。当初案は破棄したが、分割自体は正しく、blob を作らないストリーミング化（将来課題）の部品になる |

プローブはいずれも `main()` を持つ単体プログラムで、既存 API だけを使い
製品コードを変更しない。合否を標準出力に出し、終了コードで返す。

## プローブのビルドと実行

mptensor の静的ライブラリが要る（TeNeS を一度ビルドしておく）。
このディレクトリを CWD にして:

```sh
BUILD=../../out-gcc-release/build      # 自分のビルドディレクトリに読み替える
g++ -std=c++17 -O2 -D_NO_MPI -I../.. -I../../deps/mptensor/include \
    probe_direct.cpp -o probe_direct \
    $BUILD/deps/mptensor/src/libmptensor.a -framework Accelerate   # macOS
./probe_direct
```

Linux では `-framework Accelerate` を `-llapack -lblas`（または使用中の BLAS）
に置き換える。macOS + Homebrew GCC では `SDKROOT` の指定が必要なことがある:

```sh
SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk g++-16 ...
```

`probe_direct` は 22 ケース（方向 {h,v} × 演算子 {nn, hopping, pairing,
ランダム偶, source 第2サイト} × d {2,4} × D {1,2,3}）で新旧経路の
elementwise 一致を確認し、`DIRECT FUSE IDENTITY HOLDS` を出して 0 を返す。

## 恒久的な検証はどこにあるか

プローブは設計段階の使い捨てで、回帰テストではない。恒久的な保証は
`test/fermion/impurity_blob.cpp`（ctest の `test_fermion_layer`）にあり、
既存 `build_reduced_pair_naive` をオラクルとする elementwise 等価性を
実数・複素の両方、両方向、source 両位置、d=2/4、D=1〜3（d=4 D=3 は
`TENES_RUN_IMPURITY_BLOB_SLOW=1` で有効化）で固定している。

## 残っている作業記録

台帳・レビュー報告・タスク指示書は当時の作業ディレクトリ
`work/fermion/impurity-blob/`（gitignore）に残置。
