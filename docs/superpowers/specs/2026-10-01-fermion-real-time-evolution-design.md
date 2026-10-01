# fermion モードの実時間発展 — 設計書

2026-10-01、ブランチ `fermion-real-time`(develop dbf05f18 から)。

## 1. 目的と範囲

fermion モード(`[parameter.general] fermion = true`)で `mode = "time"`(実時間発展)を使えるようにする。
主な用途は、偶パリティのプロダクト状態、または保存した基底状態から始めるクエンチである。

**やること**

- solver の入力ガードを、実時間発展は通し、有限温度だけ拒否するように狭める(§3.1)。
- 実時間の SU スイープ後にも、基底状態の SU と同じ物理パリティ台帳の不変条件を検査する(§3.2)。
- `tensor_load` で、実数で保存したテンソルを complex の計算に読み込めるようにする(§3.3)。
  ボゾンにも効く既存不具合の修正である。
- simple update(SU)と full update(FU)の両方で実時間発展を通し、厳密解のある参照で検証する(§5)。
- ドキュメント(ja/en)の fermion 制限の記述と NEWS.md を更新する(§3.5)。

**やらないこと**

- 有限温度(`mode = "finite"`)。密度行列(iTPO)のアンシラにパリティの扱いが要り、
  `initialize_tensors_density()` と density CTM の経路も別なので、別件とする。入力時の拒否を続ける。
- 奇パリティのプロダクト初期状態(spinless の CDW `|1010…>` など)。各サイトテンソルを偶に保つ規約に反するため、従来どおり拒否する。
  解除するには、奇サイト同士を奇の仮想ボンドで結ぶダイマー被覆が要る(M1 設計書 §6.1 で M2 送りとしたもの)。別件とする。
- 長距離ハミルトニアンボンドの FU。基底状態と同じく `tenes_std` が拒否する現状を維持する。
- プロダクト状態から大きな D で始める fermion FU のパリティ射影の失敗(§4 P0-2)。基底状態モードでも起きる既存の問題なので、別件とする。
  ドキュメントに回避策(小さい D から始める、または SU で発展させる)を書く。
- 実時間専用の最適化(FU の毎ステップの CTM warm start の調整など)。

## 2. 現状

- `validate_fermion_constraints`(`src/iTPS/load_toml.cpp:656`)が `calcmode != ground_state` を一律に拒否している(`"non-ground-state mode"`)。
  このガードを直接テストしているものはない。
- `time_evolution()`(`src/iTPS/time_evolution.cpp`)には fermion 用の分岐がない。中で呼ぶものは、すべて fermion 経路をすでに持つ。
  - `simple_update(up)`:fermion の graded SU。ただし基底状態の `simple_update()`(`src/iTPS/simple_update.cpp:227`)が
    各スイープ後に行う「`finfo.phys[site]` が `phys_parity[site]` に戻ったか」の検査を、`time_evolution()` は行わない。
    長距離ゲート列はスイープ途中で物理台帳を一時的に書き換えるので、この検査は実時間でも要る。
  - `full_update_in_sweep(up, …)` → `full_update(up)`:fermion の FU(fast / plain)。
    `time_evolution()` は冒頭で cold な `update_CTM()` を呼ばない。FU は平均場環境と併用できない(`PEPS_Parameters::check()`)ので、
    FU のときは冒頭の `measure(t, "TE_")` が必ず `update_CTM()` を呼び、`ctm_valid_ = true` にする。
    fast FU 設計書(2026-09-07 §4.2)は、この経路を見越して `ctm_valid_` を導入済みである。
  - `measure(t, "TE_")`:fermion の CTM / 平均場測定。`time` と prefix の扱いはボゾンと共通である。
  - `fix_local_gauge()` は `Simple_Gauge_Fix = true` のときだけ呼ばれる。fermion はこれを入力時に拒否しているので、通らない。
- `run_timeevolution`(`src/iTPS/main.cpp:378`)は常に complex テンソルで走る。
- `tenes_std` は fermion でも `mode` が `time` で始まれば `tau *= 1.0j` で実時間ゲートを作る。fermion を理由に拒否する箇所はない。
  1 サイトのハミルトニアン項からは 1 サイトゲート(`make_evolution_onesite`)を作る。
- `tenes_simple` の Hubbard は `initial = "cdw"` を持つ。偶パリティ状態だけで構成され、
  **副格子 0(座標和が偶数のサイト)が空、副格子 1 が doublon** である(`tool/tenes_simple.py:399, 1606`)。
- テンソルの読み込み(`load_tensor` と `load_tensors_versioned`、`src/iTPS/saveload_tensors.cpp`)は mptensor の `Tensor::load` を使う。
  mptensor はファイルの値の型(`value_type`)と読み先の型の一致を `assert` でしか確かめない
  (`deps/mptensor/include/mptensor/file_io/load.hpp:118`)。
  このため、実数で保存したテンソルを complex の計算で読むと、Debug ビルドでは abort、Release ビルドでは double の列を complex として読んだ壊れた値になる。
  fermion ではパリティ検査がこれを「テンソルと `fermion.dat` が対応しない」という誤ったメッセージで捕まえる(§4 P0-3)。
  ボゾンでは何の検査もなく、壊れた値のまま計算が進む。
  `is_real` の既定値は false なので、ボゾンの実時間サンプル(`sample/02_time_evolution`)はこの組み合わせを踏まない。
  一方、fermion のサンプル(`sample/08_fermion_hubbard`)は `is_real = true` を使うので、基底状態からのクエンチでは必ず踏む。

## 3. 変更

### 3.1 solver の入力ガード

`src/iTPS/load_toml.cpp` の

```cpp
if (peps_parameters.calcmode != PEPS_Parameters::ground_state) {
  throw_fermion_guard("non-ground-state mode");
}
```

を、`calcmode == finite_temperature` のときだけ `throw_fermion_guard("finite-temperature mode")` を投げるように変える。
他のガード(`Simple_Gauge_Fix`、`Use_RSVD`、マルチサイト、奇パリティのプロダクト初期状態、自己隣接セルなど)は変えない。

### 3.2 実時間 SU の台帳検査

`simple_update()` のスイープ後の検査を `iTPS` のメンバ関数(例:`check_fermion_phys_ledger_restored()`)に切り出し、
`simple_update()` と `time_evolution()` の SU 経路の両方から呼ぶ。`time_evolution()` では、各ステップの全ゲート適用後に呼ぶ。
検査の内容とメッセージ(`std::logic_error`、不変条件違反)は変えない。

### 3.3 実数で保存したテンソルの complex への読み込み

`load_tensor()` と `load_tensors_versioned()` の読み込みで、ファイルの `value_type` を読み先の型と照合する。

- 一致する:従来どおり `Tensor::load` で読む。
- ファイルが実数で、読み先が complex:実数のテンソルとして読み、要素ごとに complex に昇格させる。
- ファイルが complex で、読み先が実数:`tenes::load_error` を投げる。メッセージで、complex で保存したテンソルは `is_real = true` の計算では読めないことを伝える。

`value_type` はベースファイル(`T_0.dat` など)のヘッダ 3 行目(`value_type= 0 (double)` / `value_type= 1 (complex)`)にある。
rank 0 が読み、全ランクに配る(読めない場合の例外も rank 0 だけで投げず、全ランクで揃えて投げる。既存の読み込みの作法に合わせる)。
mptensor 0.2 以前の形式(ヘッダに `value_type` がない)は従来どおり読む。

この変更はボゾンの読み込みも変える。従来は Release で壊れた値を黙って読んでいた組み合わせが、正しく読めるようになる。

### 3.4 ツール

コードの変更は予定しない。`tenes_simple` → `tenes_std` が fermion の `mode = "time"` で `input.toml` を作れることを、テストで固定する(§5 (b2) がこの経路を使う)。

### 3.5 ドキュメント

- `docs/sphinx/{ja,en}/file_specification/parameter_section.rst` の fermion の記述:
  「基底状態計算に対応」を「基底状態計算と実時間発展に対応」に改め、非対応の列挙から実時間発展を外す。
  実時間発展の初期状態は、偶パリティのプロダクト状態、ランダム偶初期化、`tensor_load` の 3 つだと明記する。
  プロダクト状態から大きな D で FU を始めると、パリティ射影で止まることがあると書き、回避策を添える。
- `tensor_load` の説明に、実数で保存したテンソルを complex の計算で読めること、逆はエラーになることを書く。
- NEWS.md:fermion の項目に実時間発展の対応を足し、Limitations から real-time evolution を外す。
  load の修正は、ボゾンにも効く修正として独立の項目にする。

## 4. P0:実現性の事前計測(実施済み)

ガードだけを外した捨てビルド(`work/fermion-real-time/p0/`、Release、`OMP_NUM_THREADS=1`)で測った。
模型はすべて `t = 1`、`U = 0`、`mu = 0` の Hubbard、初期状態は doublon/空の交互配置、`tau = 0.01`。

**P0-1:SU と厳密解**

- 1D 鎖(2 サイトセル `L_sub = [2, 1]`、`skew = 1` に鎖を埋め込み、使わない仮想脚は D = 1)。
  厳密解は n_doublon側 − n_空側 = 2·J0(4T)。
  - 横鎖・縦鎖・右上がり階段・右下がり階段の 4 通りとも、D = 16 で 6 桁まで同じ値になり、T ≤ 0.6 で厳密解との差は 6e-5 以内。
    縦ボンドの graded 経路と、raster 順に対して斜めに走る鎖の符号処理が、実時間・complex で正しく動いている。
  - D = 8(平均場環境)は 1.4 秒で、T ≤ 0.3 で差 1e-5 以内。D = 4 は T = 0.3 で 8e-4。
- 2D 正方格子。厳密解は 2·J0(4T)²(40×40 周期格子の一粒子計算で 1e-15 まで一致を確認)。
  - 平均場環境で測ると、D = 6 と D = 9 の誤差が桁まで一致して頭打ちになる(T = 0.05 で −6.2e-5、T = 0.1 で −1.0e-3)。
    平均場環境は 2D のループの寄与を落とすので、表現できている状態に対しても測定が厳密にならない。
  - CTM 環境では D とともに縮む:T = 0.05 で D = 4(chi = 16)が 1.3e-4、D = 6(chi = 36)が 2.7e-6。
    T = 0.1 では D = 4 が 1.8e-3、D = 6 が 7.2e-5(D = 6 は 10 ステップで 2200 秒かかる)。
    D = 3(chi = 9)は T = 0.1 で 2.3e-3 (2 からのずれ 0.155 の 1.5%)、1.5 秒。
  - D = 2 は T² の係数がちょうど半分になる。up と down が 1 ボンドを独立に 1 回ずつホップするだけで Schmidt ランクが 4 になるため、
    D = 2 では片方のスピンしか表せない。打ち切りの効果であり、不具合ではない。

**P0-2:FU**

- 右上がり階段の鎖で、D = 4 の fast FU と plain FU がともに完走し、SU と同じ値(差 1e-6 程度)を出す。0.7〜1.2 秒。
- D = 8 の fast FU は、プロダクト状態からの 1 ステップ目で
  `fermion full update: parity projection failed`(禁止ブロックの相対量 1.4e-8、上限 1e-8)で止まる。
  同じ入力を基底状態モード(虚時間・実数)で走らせても 2 ステップ目で同じ理由で止まるので、実時間に固有の問題ではない。
  仮想空間のほとんどが空のまま ALS を解くための悪条件と見られる。§1 のとおり別件とする。

**P0-3:load**

- complex で保存した基底状態を complex の実時間計算で読むと、そのまま動く。
- `is_real = true` で保存した基底状態を読むと、`the loaded tensor 0 breaks fermion parity under the loaded parity ledger` で止まる。
  原因は §2 の mptensor の型検査の欠落で、§3.3 で直す。

## 5. 検証

E2E は 1 本 30 秒以内、`OMP_NUM_THREADS=1`、tenes の同時実行は 2 本までとする。
参照値の比較は、golden ファイルではなく、同じビルドで回した値どうしか解析解との比較で行う。
許容誤差は P0 で測った実際の差に余裕を持たせて決め、その根拠をテストのコメントに残す。

### (a) ボソン等価性(実時間の配線)

`test/data/TE_TFI.toml` と同じ模型・ゲート・初期状態・測定で、ボゾン版と fermion 版(`fermion = true`、全サイト `parity = [0, 0]`)を同じビルドで走らせ、`TE_*.dat` が一致することを確かめる。
fermion が拒否するマルチサイト観測量と相関長は、既存の `boson_equivalence_full` と同じく両方の入力から外す。
SU 版と FU 版の 2 本を作る。グレーディングが自明なので、ここで確かめるのは測定 prefix・時刻列・complex テンソル・TE ループと fermion 経路のつなぎ込みである。

### (b1) 鎖を埋め込んだ自由フェルミオン(実時間の graded ゲートの精度)

P0-1 の 2 サイトセルに、横鎖・縦鎖・右上がり階段・右下がり階段の 4 通りの鎖を埋め込み、SU(平均場環境)で n_doublon側 − n_空側 を 2·J0(4T) と比べる。
D = 8、T ≤ 0.3 を目安とする。右上がり階段では FU(D = 4、fast と plain)も 1 本ずつ走らせる。
厳密解に 1e-5 の精度で合うので、ゲートの符号・位相・複素共役の誤りを鋭く検出できる。

### (b2) 2D の自由フェルミオン(ループを含む実時間発展と正規経路)

`tenes_simple` → `tenes_std` → `tenes` の正規経路で、2D 正方格子の Hubbard を `initial = "cdw"` から発展させる。
CTM 環境(D = 3、chi = 9 を目安)で、T = 0.1 での n_doublon側 − n_空側 を 2·J0(4T)² と比べる。
この比較の精度は 1e-3 程度で、主に確かめるのは 2D の両方向のボンドを同時に使う発展と CTM 測定が complex で通ることである。
ループに由来する符号(プラケットを回る交換)の効果は T⁴ の項に入り、この精度では見分けられない。
ループの符号は基底状態のテスト群が、モードに依存しない同じコード(SU と CTM)で押さえている。

### (c) 1 サイトゲートと時間の向き

1 サイトのペア場 H = Δ(|↑↓⟩⟨0| + h.c.) だけを持つ入力で、真空のプロダクト状態から発展させる。
厳密解は |ψ(T)⟩ = cos(ΔT)|0⟩ − i·sin(ΔT)|↑↓⟩ なので、
- 密度 ⟨n⟩ = 2·sin²(ΔT)
- 偶パリティの非対角 1 サイト演算子 |0⟩⟨↑↓| の期待値 = −i·cos(ΔT)·sin(ΔT)(虚部の符号で時間の向きがわかる)

と比べる。fermion の 1 サイトゲート経路(`apply_onesite_gate_fermion`)を complex で通し、
e^{−iHT} と e^{+iHT} の取り違えや複素共役の誤りを検出する。CDW の密度は時間反転で不変なので、(b1)(b2) ではこれを検出できない。
SU と FU の両方で確かめる。

### (d) load からのクエンチ

- fermion の基底状態を `is_real = true` で計算して `tensor_save` し、別のハミルトニアンの実時間計算で `tensor_load` して始める。
  時刻 0 の測定値が、保存元の基底状態計算の測定値と一致し、そのまま完走することを確かめる。complex で保存した場合も同様に確かめる。
- complex で保存したテンソルを `is_real = true` の基底状態計算で読むと、`load_error` で止まることを確かめる。
- ボゾンでも、実数で保存したテンソルを complex の計算で読んだ値が、保存元と一致することを確かめる(§3.3 の回帰テスト)。

### (e) ガードの単体テスト

- fermion で `mode = "time"` の入力を `validate_fermion_constraints` が受け付ける。
- fermion で `mode = "finite"` の入力は、メッセージに `finite-temperature` を含む `input_error` で拒否される。
- 他のガード(例:`Simple_Gauge_Fix = true`)は、`mode = "time"` でも引き続き拒否される。

### 検出できない欠陥クラス

- 2D のループに由来する符号の誤りのうち、実時間でだけ現れるもの。(b2) の精度では見えない。
  SU・CTM・測定のコードは虚時間と共通なので、基底状態のテスト群でカバーされているとみなす。
- §3.2 の台帳検査は、正しい実装では発火しない不変条件の検査なので、E2E では発火を確かめられない。

## 6. リスク

- **R1:FU の悪条件。** P0-2 の D = 8 の失敗は、D = 4 の FU テストでは現れない。
  別件として扱うが、ユーザーが最初に当たりやすいので、ドキュメントに回避策を書く。
- **R2:load の修正がボゾンの挙動を変える。** 従来は壊れた値を黙って読んでいた組み合わせなので、正しく読むようにする変更で既存の正しい結果は変わらない。
  既存の checkpoint 系テスト(`checkpoint_dir.py` など)を全件回して確かめる。
- **R3:平均場環境での測定の誤差。** 2D で平均場環境を使うと、表現できる状態に対しても測定が厳密にならない(P0-1)。
  テストでは 2D の厳密解比較に CTM を使う。ユーザー向けには注意書きを足さない(ボゾンでも同じ性質で、既存のドキュメントの範囲)。

## 7. 進め方

CLAUDE.md の多段手順どおりに進める。

1. この設計書を Codex にレビューさせ、指摘を反映する(初版のレビューは `work/fermion-real-time/review/codex-design-review.md`。
   指摘 4 件のうち 2 件(台帳検査、CDW の向き)を §3.2・§2 に、1 件(1 サイトゲート)を §5 (c) に、1 件(`ctm_valid_` の条件)を §2 に反映した)。
2. 振る舞い契約書(`docs/superpowers/specs/2026-10-01-fermion-real-time-evolution-contract.md`)を書き、テスト作成者が §5 のテストを書く。
   RED を 1 件ずつ確認し、`test/` 全体のスナップショットを取る。
3. Codex が §3 を実装する(タスク分割:T1 ガードと台帳検査、T2 load の型照合、T3 ドキュメント)。テスト不改変を diff で検査し、Claude が ctest 全件を回す。
4. タスクごとのレビューのあと、最後に最上位モデルで全ブランチレビューを 1 回行う。
