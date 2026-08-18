# フェルミオン対応 M1 設計書 — swap gate 方式・正方格子・最近接模型

日付: 2026-08-15
ブランチ: `fermion`
前提調査: `docs/superpowers/notes/2026-08-13-fermion-survey-summary.md`
(swap gate / Z₂グレード / Grassmann の3定式化と既存実装のサーベイ要約)

## 1. 背景と目的

現行 TeNeS はテンソル縮約においてすべての演算子・状態の入れ替え符号を +1 として扱うため、
フェルミオン系(Hubbard、t-J、スピンレスフェルミオン等)を扱えない。本設計は swap gate 方式
(Corboz–Orús–Bauer–Vidal, PRB 81, 165104 (2010))に基づき、密テンソル(mptensor 無改造)の
まま TeNeS リポジトリ内で完結するフェルミオン対応の第一段階(M1)を定める。

計算量への影響: 符号処理は要素ごとのマスク乗算であり、縮約のリーディングオーダーは
ボゾン系と同一(Barthel–Pineda–Eisert, PRA 80, 042333 (2009))。

## 2. スコープ

### M1 でやること

- 正方格子・**最近接**ボンドハミルトニアン・スピンレスフェルミオン
- simple update(虚時間発展)
- CTMRG 環境構築+ **1サイト・最近接2サイト観測量**(エネルギー、占有数、NN ホッピング相関)
- 入力: `input.toml` 拡張(fermion フラグ+物理脚パリティ宣言)+検証用入力生成スクリプト
- 実・複素両テンソル型対応(既存テンプレート構造を維持)
- 合格基準: **自由スピンレスフェルミオン(ホッピング t + 化学ポテンシャル μ)の
  厳密解(運動量積分)に対するエネルギー密度・占有数の D 収束**

### M1 でやらないこと(fermion モード時は load_toml 段階で明示エラー)

- full update / fast full update
- 距離 2 以上の相関関数・multisite 観測量(明示的に指定された場合はエラー)
- 相関長(転送行列): **デフォルト有効のため、エラーではなく強制無効化+警告**(§6.1)
- 有限温度(purification / density 系)・実時間発展
- mean-field 環境での測定
- RSVD(fermion モードではフル SVD を強制)
- simple update のゲージ固定(`Simple_Gauge_Fix`)
- **checkpoint(`tensor_save` / `tensor_load`)**: 既存の保存形式(`saveload_tensors.cpp`)は
  形状とテンソル本体のみでパリティメタデータを持たないため、fermion モードでは両方を
  エラーにする。パリティ直列化は M2(§9)
- `tenes_simple` / `tenes_std` の本格対応(M2)
- U(1) ブロックスパース化(性能最適化として将来課題)

具体的なガード実装は §6.1 のガード表で定める。

## 3. アーキテクチャ: 符号レイヤー

手作業で縮約図ごとに swap 挿入位置を導出する代わりに、**脚ごとのパリティメタデータを持つ
薄いラッパー `ftensor` と、符号を機械的に生成する f-プリミティブ関数群**を導入する
(locally-ordered 簿記。数学的には交差数え上げと同値だが、符号バグの混入面が小さい)。
密テンソルのまま・mptensor 無改造・ブロックストレージなし。

### 3.1 ftensor

```cpp
template <class tensor>
struct ftensor {
  tensor t;                                  // 密テンソル本体(mptensor 型)
  std::vector<std::vector<bool>> leg_parity; // 脚ごと・添字ごとのパリティ
};
```

パリティを「セクター次元 (D_e, D_o)」ではなく**添字ごとの任意 bool ベクトル**で持つ。
reshape(脚の融合)後の非連続パターンもそのまま扱える。

### 3.2 f-プリミティブ(mptensor 自由関数のミラー)

| 関数 | 意味論 |
|---|---|
| `ftranspose(A, axes)` | transpose + 転倒対 (i<j が逆転) ごとの (−1)^{p_i p_j} を要素マスクで適用 |
| `ftensordot(A, B, axesA, axesB)` | 「A の縮約脚を末尾へ、B の縮約脚を先頭へ寄せる置換」の符号+縮約。縮約脚のパリティ一致を実行時検証 |
| `fsvd` / `fqr` | パリティセクターごとの分解。特異値の大きい順にセクター横断で切断し、偶先頭に再ソートして新パリティを返す。縮退時の選択は決定論的(偶優先) |
| `fconj(A)` | 要素共役+脚順序反転(随伴)+元の軸番号への ftranspose。正味は「共役+全脚対マスク (−1)^{Σ_{i<j} p_i p_j}」(最終形は規約固定テストで確定) |
| `apply_swap(A, ax1, ax2)` | 奇⊗奇要素の符号反転(低レベル) |
| `apply_parity(A, ax)` | 指定脚の奇要素に −1(`multiply_vector` で実装可) |
| `check_parity_even(A)` | 禁止ブロックのノルム検査(デバッグ・テスト用) |

要素ごとアクセスは `tensors.cpp:114-125` の `local_size()`/`global_index()`/`set_value()`
ループが前例。

### 3.3 既存コードとの統合(境界ラップ方式)

- `iTPS<tensor>` に `FermionInfo` メンバを追加: 有効フラグ+物理脚パリティ(サイトごと、
  入力から)+仮想ボンドパリティ(動的、`fsvd` のたびに更新)。`Tn`・環境テンソルの格納型は
  既存のまま。
- フェルミオン化するルーチンの呼び出し境界で `ftensor{テンソル, パリティ}` を構築して渡し、
  戻りでメタデータを回収する。テンプレートカーネル自体は `ftensor` でインスタンス化し、
  自由関数のオーバーロード解決で符号付きに切り替える(**カーネルソースはボゾン/フェルミオン
  共通**)。`ftensor` は M1 対象コードが使う API サブセットのみミラーする
  (member: `transpose`/`multiply_vector`/`shape` など)。
- `iTPS<ftensor>` の全面インスタンス化は行わない(未移植機能への波及を避ける)。
- ボゾン経路はコード・golden ともに完全不変。

## 4. 規約

1. **基底順序**: 物理脚のパリティは input.toml で基底状態ごとに宣言(偶先頭を推奨規約)。
   仮想脚は `fsvd`/`fqr` の出力規約として偶先頭ソートで返す(スライスベースの符号マスクを
   高速化する最適化)。ただし**プリミティブは任意の bool ベクトルで正しく動作する**ことを
   要件とし、ソート済みであることに依存してはならない(CTM 環境更新の `extend` による
   χ パディング等でソートが崩れるため。パディング要素のパリティは偶を割り当てる)。
2. **脚の正準順序**: 各 ftensor の格納軸順序がそのテンソルのフェルミオン順序。
   Tn の既存規約 (0=左, 1=上, 2=右, 3=下, 4=物理) を正準順序として採用。
   グローバルなサイト順序は不要。
3. **bra/ket**: `fconj` を §3.2 の式で定義。共役まわりの符号配置には規約の自由度があるため、
   **最終形は規約固定テスト(§7 層2)が決定する**: (i) ⟨ψ|ψ⟩ > 0、(ii) 有限ミニネットワーク
   縮約が JW 順序の厳密参照と一致、(iii) 水平・垂直の整合。
4. **演算子**: 2サイトゲートは既存 rank-4 `op[in1,in2,out1,out2]` 形式のまま、行列要素は
   **順序付き2サイト基底 |n₁n₂⟩ = (c†₁)^{n₁}(c†₂)^{n₂}|0⟩(site 1 = source_site)で定義**。
   NN スピンレスではゲート行列は符号なしで書け、幾何由来の符号はネットワーク側が負担。
   ロード時に宣言パリティに対する偶性を検証し、奇な演算子は M1 では拒否。
   **演算子・ゲートも ftensor として扱い**、脚パリティは対象サイトの物理脚パリティから
   付与する。**source/target 逆向きの測定で必要になる並べ替えは
   `ftranspose(op, {1,0,3,2})` で行う**(現行 `twosite_obs.cpp:190,209` の
   `mptensor::transpose(op.op, {1,0,3,2})` の置き換え。プレーンな transpose では
   奇×奇の符号が落ちる)。
5. **λ**: 実正値対角のまま。`fsvd` が特異値と同時にボンドのパリティベクトルを返すため、
   λ の各要素は対応する添字のパリティと常に整合する(**ソート順には依存しない**。
   `FermionInfo` のボンドパリティと λ は同一の添字順で管理する)。対角・偶なので
   `multiply_vector` は符号を生まない。

## 5. アルゴリズム変更

### 5.1 simple update(`core/simple_update.cpp` `Simple_update_bond`)

| ステップ | 変更 |
|---|---|
| λ 吸収(`multiply_vector`) | 変更なし(対角・偶) |
| 方向別 transpose(4分岐) | `ftranspose`(符号自動。方向別の手動符号解析は不要) |
| QR | `fqr` |
| θ = R1·R2·gate | `ftensordot` ×2 |
| SVD+切断 | `fsvd(θ, ..., dc)`。現行の「svd → slice で打ち切り」パターンは、セクター横断の特異値選択が必要なため **切断次元 dc を引数に取る fsvd に置換**する。ボンドのパリティベクトルが動的に決まり、`FermionInfo` と λ を更新 |
| λ⁻¹ 戻し+逆 transpose | `ftranspose` |

### 5.2 CTMRG(`core/ctm.cpp` + `contract_itps_ctm/ctm.cpp` の1・2サイト測定カーネル)

- 環境テンソル C1..C4 / eT の全脚にパリティベクトルを付与(D 脚はボンドパリティの
  ket/bra コピー、χ 脚は射影子計算時の `fsvd` が返す動的パリティ)。初期化も f-プリミティブで。
- 射影子: `fsvd(ftensordot(...))`。吸収(`Left/Right/Top/Bottom_move`)と格子回転 transpose は
  `ftranspose`/`ftensordot` への機械的移植。
- 収束判定(C の特異値比較)は変更なし。
- bra/ket 符号は `fconj`+呼び出し列から自動発生。**文献規約(Corboz PRB 81, 165104 付録)で
  手導出した1サイト reduced tensor と f-プリミティブ経由の結果を突き合わせる橋渡しテストを
  1本置く**。

### 5.3 測定

- 1サイト(n など偶演算子): 既存カーネルの ftensor 移植。
- NN 2サイト(ホッピング等): `Contract_two_sites_{horizontal,vertical}_op12` の移植。
  演算子は §4-4 の規約で構築されたものを使用。
- `twosite_obs.cpp` の距離 2 以上の SVD 分解経路(`:226`)は fermion モードでガード
  (JW string 挿入が必要なため M2)。

## 6. 入力スキーマ・初期状態

- `[parameter.general] fermion = true`(モードスイッチ)
- `[[tensor.unitcell]] parity = [0, 1]`(長さ physical_dim の 0/1 リスト。fermion 時必須)
- 仮想脚の初期パリティ分割は**脚ごと**にデフォルト D_e = ⌈D_leg/2⌉
  (D_leg = `lattice.virtual_dims[site][leg]`。入力指定不要。以後 `fsvd` が動的に決定)
- 検証: パリティ宣言の整合、演算子の偶性、§6.1 のガード表

### 6.1 ガード表(fermion = true 時)

「沈黙して間違った数値を出さない」を原則に、到達経路ごとに具体的なガードを置く。
デフォルトで無効な機能を明示的に有効化した場合は**エラー**、デフォルトで有効な機能は
**強制無効化+警告**とする。

| 項目 | 現行の到達経路 | 挙動 |
|---|---|---|
| `mode` が ground 以外(実時間・有限温度) | `main.cpp` の mode 分岐 | エラー(load 時) |
| full update(`num_full_step > 0`) | `optimize()` | エラー(load 時) |
| `MeanField_Env = true` | `measure()` の分岐・`twosite_obs.cpp` の MF 分岐 | エラー(CTM 必須) |
| `Simple_Gauge_Fix = true` | `simple_update()` → `local_gauge.cpp` | エラー |
| `Use_RSVD = true` | `core/ctm.cpp:229-242` の射影子 | エラー |
| multisite 演算子の定義 | `measure()` が無条件に `measure_multisite()` を呼ぶ | 演算子定義時点でエラー(未定義なら空で通過) |
| `correlation`(`r_max > 0`) | `measure()` | エラー |
| 相関長(`correlation_length.measure`、**デフォルト true**) | `measure()` → `tmatrix_param.to_calculate` | 強制 false+警告 |
| `tensor_load_dir` 非空 | `load_tensors()`(パリティ復元不能) | エラー |
| `tensor_save_dir` 非空 | `save_tensors()`(パリティ未保存) | エラー |
| 距離 2 以上の 2 サイト演算子 | `twosite_obs.cpp:226` の SVD 分解経路 | 演算子定義時点でエラー |
| パリティ奇の演算子・ゲート | `load_toml` の演算子読込 | エラー(偶性検証) |
| 奇パリティのプロダクト初期状態 | `tensors.cpp` 初期化 | エラー(§9 のダイマー被覆で M2 解除) |
- **初期状態**: ランダム偶初期化(既存ランダム初期化ループに偶ブロックマスク)を推奨経路とする。
  偶パリティのプロダクト状態+偶ノイズも許可。**奇状態のプロダクト初期化はエラー**とし、
  メッセージで理由(奇テンソルは偶制約に反する)と回避策(ランダム偶初期化)を案内する。
- 出力ファイル形式(`output/*.dat`)は変更なし。

## 7. テスト計画(4層)

1. **単体(doctest)**: `apply_swap`/`ftranspose` を素朴な要素ごと参照実装と比較。
   `ftensordot` の縮約順序独立性・パリティ不整合検出。`fsvd`/`fqr` のブロック厳密性・
   セクター横断切断・偶先頭ソート・因子の偶性。
2. **規約固定**: 有限ミニネットワーク(2×2、小 D)の縮約値を JW 順序の numpy/ED 厳密参照と
   比較し、`fconj` の規約を確定する。文献の1サイト reduced tensor との橋渡しテストを含む。
3. **統合(ctest)**:
   - **1D 極限**: t_y=0 でフェルミオン計算とハードコアボソン計算(既存経路)が一致し、
     t_y≠0 で乖離すること(符号機構の狙い撃ち)。
   - **自由フェルミオン D 収束**: t=1、μ 数点、D=2..6。合格基準は
     (i) 最大 D で厳密値(運動量積分)との相対誤差が閾値内、
     (ii) 各 D の値を golden(許容誤差付き)として固定する回帰。
     **単調性は要求しない**(simple update は厳密変分でないため単調収束の保証がない)。
     閾値は実測して決定する。
4. **回帰・ガード**: 既存 ctest 全件が golden 据え置きで不変。fermion+未対応機能の
   組み合わせがエラーになることの検査。

## 8. ファイル構成・見積り・実装順序

### 新規

- `src/fermion/ftensor.hpp` — ftensor 構造体+member API ミラー
- `src/fermion/fops.hpp` / `fops.cpp` — f-プリミティブ群
- `src/fermion/fermion_info.hpp` — FermionInfo
- `test/test_fermion_layer.cpp` — 単体+規約固定テスト
- `test/fermion/` — 自由フェルミオン ctest(入力生成+厳密値スクリプト)

### 変更

- `src/iTPS/iTPS.hpp`(FermionInfo メンバ)、`load_toml.cpp`(パリティ解析+ガード)、
  `tensors.cpp`(偶マスク初期化)、`simple_update.cpp` / `measure.cpp` / `onesite_obs.cpp` /
  `twosite_obs.cpp`(境界ラップ+ガード)、`core/simple_update.cpp` / `core/ctm.cpp` /
  `contract_itps_ctm/ctm.cpp`(ftensor インスタンス化追加)、`PEPS_Parameters`、ドキュメント

### 見積り

約 3,000 行(符号レイヤー+単体テスト 1,500〜2,300、統合+入力+統合テスト 1,000〜1,200)

### 実装順序(M1 内 5 段階)

1. 符号レイヤー TDD(TeNeS 統合なし、doctest のみ)
2. 規約固定テスト(fconj 確定)
3. simple update 移植+1D 極限テスト
4. CTMRG+測定移植+自由フェルミオン D 収束テスト
5. 入力スキーマ・ガード・ドキュメント

## 9. 将来課題(M2 以降)

- **奇パリティ占有配置のプロダクト初期化**(2026-08-15 記録): セルあたり偶数フェルミオンの
  占有配置は「奇仮想ボンド(次元1の奇セクター)で占有サイトをペアにするダイマー被覆」で
  厳密に表現できる(例: 全占有 |11…⟩ は 2×1 セル+奇偶交互ボンド、チェッカーボードは
  2×2 セル+束の素通し)。ダイマー被覆の自動生成を実装し、奇プロダクト初期化の
  エラーを解除する。見積り 150〜300 行。
- **パリティメタデータの直列化**: `tensor_save`/`tensor_load`(checkpoint・restart)の
  fermion 対応。保存ディレクトリに FermionInfo(物理・仮想ボンド・CTM χ 脚のパリティ
  ベクトル)のサイドカーファイルを追加し、フォーマットバージョンを上げる。
- `tenes_std`/`tenes_simple` のフェルミオン対応(フェルミオンモデル定義、
  多ホップ MPO 分解への JW string 挿入)
- 距離 2 以上の相関関数(パリティ演算子挿入付き転送)・相関長のパリティ奇セクター
- full update / fast full update、simple update ゲージ固定、RSVD のパリティ対応
- 汎用 NxM カーネルの fermion 対応: 生成器(`misc/contraction/`)に**平面埋め込みを持つ
  中間表現**を導入(三角格子 iTPS 計画と同じ改修点。幾何一般で設計する)
- 有限温度(purification のアンシラパリティ)、実時間発展
- U(1) ブロックスパース化(実用 D 引き上げ。swap gate 資産は YASTN 型アーキテクチャで温存可)

## 付録A: ftensor API 互換表(M1 呼び出し面の全数)

core 層 3 ファイルの grep 全数調査(2026-08-15、`fermion` ブランチ)に基づく。`ftensor` は
この表の全項目を実装しなければならない。ドライバ層の呼び出し面は実装時に同様の census を
取って本表に追記し、**表外 API の使用が見つかった場合は表を更新してから実装する**。

| API | simple_update.cpp | core/ctm.cpp | contract_itps_ctm/ctm.cpp | パリティ意味論 |
|---|---:|---:|---:|---|
| `tensordot` | 10 | 66 | 139 | `ftensordot`: locally-ordered 符号+縮約脚パリティ検証 |
| `transpose`(member/free) | 14 | 41 | — | `ftranspose`: 転倒対符号 |
| `conj` | — | 23 | 13 | `fconj`: §4-3 の規約 |
| `svd` | 1 | 13 | — | `fsvd`: セクター別分解+横断切断(dc 引数版で svd→slice を置換) |
| `qr` | 2 | — | — | `fqr`: セクター別分解 |
| `slice` | 2 | 42 | — | パリティベクトルを同区間でスライス(符号なし)。simple_update の切断用 slice は fsvd(dc) に吸収 |
| `extend` | — | 24 | — | ゼロパディング。追加要素のパリティは偶を割当(ソート不変条件に依存しない設計が前提、§4-1) |
| `reshape` | — | 12 | — | 脚融合: 合成パリティ(XOR パターン)を計算。脚分割: 分割後パリティを引数で明示指定するオーバーロードを用意 |
| `multiply_vector` | 18 | 8 | — | 対角・偶(λ / 特異値の重み)なので符号なし。転送のみ |
| `trace` | — | — | 11 | 全縮約。locally-ordered 符号 |
| `max_abs` | — | 3 | — | 符号なし(転送のみ) |
| `rsvd` | — | 2 | — | **M1 非実装**(Use_RSVD ガードで到達不能。static_assert かガード済みブランチで除外) |
| `shape` / `Axes` / `Shape` | 多数 | 多数 | 多数 | メタデータ透過 |

注: ドライバ層(`simple_update.cpp` / `measure.cpp` / `onesite_obs.cpp` / `twosite_obs.cpp`)は
原則境界ラップのみだが、**例外として `twosite_obs.cpp` は §4-4 の演算子向き反転
`ftranspose(op, {1,0,3,2})` を直接呼ぶ**。その他の表外 API 使用は実装時の census で
洗い出し、表に追記してから実装する。

## 10. 参考文献

- P. Corboz, R. Orús, B. Bauer, G. Vidal, PRB 81, 165104 (2010), arXiv:0912.0646 — swap gate 方式の原典(fPEPS+CTMRG)
- T. Barthel, C. Pineda, J. Eisert, PRA 80, 042333 (2009), arXiv:0907.3689 — 縮約コスト同一性
- Q. Mortier et al., SciPost Phys. 18, 012 (2025), arXiv:2404.14611 — 統一レビュー(graded 形式との等価性)
- B. Bruognolo et al., SciPost Phys. Lect. Notes 25 (2021), arXiv:2006.08289 — フェルミオン iPEPS+CTMRG 実装手引き
- M. M. Rams et al., SciPost Phys. Codebases 52 (2025), arXiv:2405.12196 — YASTN(swap gate 実装規約)
- 調査レポート: https://claude.ai/code/artifact/653868dd-dff9-4fb7-9333-98a8fe81b69e

## 改訂1(2026-08-15): 脚 arrow(dual フラグ)の導入

**背景**: §4-3 の「arrow を持たず、必要箇所に明示挿入」という簡略化は、閉ループを含む
ネットワーク(CTM 測定の χ リング等)で破綻することが E2E 検証で判明した。パリティ
ベクトルだけでは「ket-ket ボンド縮約」と「ket-bra 物理縮約」を区別できず、奇セクターを
通る閉ループが supertrace の (−1) を拾う(症状: ランダム偶 iPEPS のサイトノルムが負、
密度が [0,1] を逸脱)。開いたネットワーク(Task 9 の JW 検証)や等長性(U†U=1)は
現行実装でも正しい。

**設計**(symmray/YASTN と同型の dual フラグ):

1. `ftensor` に `std::vector<bool> dual;`(脚ごと、false=ket / true=dual)を追加。
2. **縮約規則**: `tensordot`/`trace` は縮約脚対の dual が異なることを実行時検証
   (同じならエラー)。A 側の縮約脚が ket(dual=false)の場合、その脚に
   `apply_parity` を適用してから縮約する(coevaluation の向き補正)。A 側が dual なら
   符号なし。この単一の局所規則が閉ループの supertrace を trace に直す。
3. `conj`: 全脚の dual を反転+要素共役+既存の (−1)^{m(m−1)/2} twist。
4. `transpose`: dual を parity と同時に並べ替え。`slice`/`extend`/`reshape` も同様に転送。
5. `svd`/`qr`: 内部脚は U(最終脚)= ket、VT(先頭脚)= dual として生成
   (U·VT の再構成縮約が規則 2 と整合)。
6. **格子規約**: ket 層 Tn の脚 dual = {左: true, 上: true, 右: false, 下: false,
   物理: false}(右・下が ket、左・上が dual。隣接ボンドは常に ket-dual 対になる)。
   演算子 rank-4 は {in1, in2: dual=false→…} — 具体的には in 脚が ket-Tn の物理
   (ket)と縮約されるため dual=true、out 脚が bra 物理(dual)と縮約されるため
   dual=false。rank-2 も同様。
7. 環境テンソル C/eT の dual は構成(縮約・分解)から自動的に伝播する。

**検証**: (a) supertrace ミニテスト(単一奇ループの向き両方で ± を確認)、
(b) 既存全単体テストのフラグ対応と Task 9 JW 物理値の再固定、
(c) ランダム偶 iPEPS の CTM ノルム正・n∈[0,1] を回帰テスト化(今回の再現ケース)、
(d) 自由フェルミオン E2E(最終判定)。

## 改訂2(2026-08-15): 測定・環境は reduced tensor + 既存 density-CTM 経路へ(改訂1を置換)

**背景**: 改訂1の arrow(dual フラグ)実装でも、ランダム偶 iPEPS の CTM 測定で負ノルムが
解消しなかった(3方式目の規約合わせ失敗)。根本問題は「ボゾン用に書かれた手書き CTM の
全縮約呼び出し列に対して、機械的な符号簿記の規約整合を後付けする」というアプローチ自体の
脆さにある。改訂1(arrow 設計)は**破棄**する。

**確定している事実**: 開いたネットワーク(JW 厳密検証・等長性・simple update)は
pre-arrow の f-プリミティブで正しい。壊れるのはループを閉じる環境・測定のみ。

**新方針**(文献の標準構成):

1. **simple update は現行の ftensor 経路のまま**(変更なし。検証済み)。
2. **測定・環境**: サイトごとに二層 reduced tensor
   a[(l,l̄),(t,t̄),(r,r̄),(b,b̄)] を **f-プリミティブ(開ネットワーク)で明示的に構築**し、
   フェルミオン符号をこの構築関数一箇所に集中させる。粗視化格子上のネットワークは
   交差のない通常のテンソルネットワークになるため、**環境構築・測定は既存の単層
   (density/TPO)CTM 一式(`ctm_single.cpp` + `contract_density_ctm/`)を無改造で流用**する。
3. 観測量: 偶 1 サイト演算子は impurity reduced tensor(演算子挿入版)。NN 2 サイト
   演算子は odd⊗odd 分解の補助脚をボンドに沿わせた impurity 対で扱う。
4. 検証: 有限パッチ(2 サイト水平/垂直、2×2 プラケット=ループ含む)について、
   (a) numpy JW 厳密参照(gen_reference.py の 2D 拡張)、(b) f-プリミティブによる
   開ネットワーク直接縮約、(c) reduced tensor +ボゾン縮約、の 3 者一致を要求する。
   これが reduced 構築関数の規約を一意に固定する。
5. Task 14–15 で入れた ftensor 版 CTM(`ctm.cpp` 等のインスタンス化)は fermion モードでは
   使用しなくなる(ボゾン経路に影響はないため当面残置、M1 仕上げで整理)。

**利点**: ループ系の規約バグのクラスを構造的に根絶。メモリは D⁴ 級に増えるが M1 の
検証規模(D≤4)では問題ない。U(1) ブロック化(将来)とも両立。

## 性能課題(2026-08-18 計測): doubled_cluster の D スケーリング

2サイト reduced pair blob は 6 本の融合脚(各 D²)を持つ D¹² 級オブジェクト
(D=4 で 2.7e8 要素 ≈ 2 GB)。現実装は bra⊗ket の外積(rank-16)を作ってから
要素ごとの get/set ループで joint-swap マスクを適用するため、D=4 の 2 サイト
測定 1 ケースに CPU 6 時間以上を要する(D=2 では数秒)。E2E(FreeFermion ctest)
の実行時間の支配項。改善案: (a) マスク適用を「bra 側要素ループ+ket 側ベクトル
演算」に変えて ~100× 高速化、(b) 融合前に CTM 縮約へストリーム化して D¹² の
実体化自体を回避、(c) ctest では D=4 ケースを nightly ラベルに分離。M1 では
(c) を暫定採用し、(a) を最初のフォローアップとする。

## 改訂3(2026-08-15): 2サイト演算子の graded 表現規約の確定(τ² 残差の根本原因)

**症状**: 1D 極限 λ 軌跡診断(criterion (a))が wrap ボンドで τ² スケール乖離。

**根本原因は 2 つの規約エラーの部分相殺**:

1. **ゲートの素ロード**: 2サイト演算子 op[in1,in2,out1,out2] を素の行列成分
   `<o1 o2|O|i1 i2>` のまま ftensor 化していた。graded tensordot の機械的マスク
   (`(-1)^{|in1||in2|}`、B 側縮約脚の逆順化由来)が **in 脚両方奇のチャネル**
   (例: exp(τh) の `|11><11|` 要素)を反転させる。正しい canonical 表現は
   `fop = (-1)^{|i1||i2|} · <o1o2|O|i1i2>`、実装は **`wrap_twosite_op`**
   (素ロード後に `apply_swap(fop, 0, 1)`)。Task 9 のホッピング演算子はこの
   チャネルに要素を持たないため検出不能だった。**判別テスト**: `n⊗n`(非ゼロ要素が
   両奇チャネルのみ)の 2 サイト適用 vs 検証済み 1 サイト適用合成
   (test_fermion_layer.cpp「two-site operator doubly-odd input channel …」)。
2. **regroup の符号なし転置**: `regroup_theta_for_svd`(ftensor 版)が
   コミット 0d4c82f5 で「素 transpose + parity メタデータ並べ替え」になっていた。
   (out1, aux2) 交差の Koszul 符号が欠落。正しくは graded transpose。

**相殺の構造**: 素ロードの余剰マスク `(-1)^{m1·m2}` と regroup 欠落分
`(-1)^{o1·q2}` は、パリティ偶テンソル上で恒等(対角)チャネルに対しては合成が
「行パリティにのみ依存する対角ゲージ」に退化しスペクトル不変 → update 0 では
一致し、ホッピング(O(τ))チャネルにのみ符号残差 → τ² スケールの λ 乖離として
顕在化していた。両方修正後、criterion (a) は **300 ステップ機械精度(5.8e-14)で
case=A** を達成。

**経路ごとの演算子規約(重要)**:

| 経路 | 規約 | 根拠 |
|---|---|---|
| simple update カーネル(直接 f-プリミティブ縮約) | **wrap_twosite_op(swap ロード)** | nn 判別テスト+λ軌跡 300 步機械精度 |
| reduced pair blob(`build_reduced_pair`)測定 | **素ロード(plain)** | R3 oracle 対査(density_pair/mixed の両奇チャネル込みで 1e-12 一致) |

blob 経路は doubling+joint-swap+gauge の経験的固定に規約が焼き込まれているため
素ロードが正。M2 以降で単一規約への統一を検討(twosite_obs.cpp の NOTE 参照)。
1 サイト演算子(rank-2)は in 脚交差が存在せず素ロードで両経路とも正しい。

**測定経路の 2 バグ(E2E「CDW 崩壊」の真因は測定側)**: 弱2D結合診断
(t_y=0.1、TENES_RUN_WEAK2D_DIAG)で、λ ゲージ平均場推定(f-プリミティブのみ、
CTM 非依存)と本番 CTM 測定を bond ごとに比較して確定:

1. **λ の二重計上**: reduced tensor 構築前に Tn へ full-λ dressing をしていた。
   TeNeS の規約ではカーネルが √Schmidt を bond 両端の Tn に書き込むため
   **状態 = bare Tn の直接縮約**(ボゾン CTM 経路と同じ)。full-λ dressing は
   MeanField 経路の規約であり、CTM に渡すと環境重みが二重になる
   (λ が縮退に近い状態では偶然無害 → 1D 鎖検証をすり抜けた)。
2. **wrap ボンドの blob ペア順序逆転**: `wraps_right`/`wraps_up` 分岐が wrap
   ボンドで build_reduced_pair の (A,B) を逆順にしていた。窓の indices は既に
   幾何役割(col0=左、row0=上)で正しく、無限格子に境界はないため wrap の
   特別扱いは丸ごと不要。逆転は blob と環境テンソルの割り当て不整合を生み
   wrap ボンドの測定値が ~0 に死んでいた(過去の検証は全サイト同一テンソル
   だったため逆転が不可視だった)。

修正後、t_y=0.1・D=2 で全 8 ボンドが MF 推定と一致、E/site = −0.639 vs 厳密
−0.638(0.2%)。なお軌跡(simple update)は開放 2×2 プラケットのカーネル vs
厳密 Trotter 比較(TENES_RUN_PLAQUETTE_TROTTER_DIAG、無切断域で 1e-16)で健全と
証明済み — E2E の見かけの「市松 CDW・負ノルム」はすべて測定アーティファクト
だった。

**n≠0.5 の収束性(バグではなくアルゴリズム的性質)**: 全規約層の検証完了後も
μ=−1 の 2D E2E(固定 τ=0.01、1000 步 = 虚時間 10)は未収束だった。切り分け:
(a) ボゾン対照は同一パラメータで一様収束 → フェルミオン特有、(b) 1D 極限
(横・縦・μ 込み)は case=A 機械精度 → 規約は無罪、(c) τ アニーリング
(0.1→0.05→0.02→0.01、各250步、虚時間45)で D=4 が収束、(d) 固定 τ=0.01 ×
4500 步(虚時間45)の対照も**同一固定点**に収束 → アニーリングは本質でなく
**虚時間不足**。結論: 偶パリティ超選択則下の simple update は n≠0.5 で収束が
大幅に遅い(必要虚時間 ~45 vs μ=0 の ~10)。収束先は周期2の対称性破れ状態で、
エネルギーは厳密値に対し ~0.1%(χ=8 測定)だが site 平均密度は ~20% ずれる
(変分的にエネルギーのみ正確)。E2E は μ=−1 の num_step を 4500 に増やし、
max-D 厳密閾値チェックは energy を両 μ、density は μ=0 のみに適用する。

**ボンド方位の正規化(2D 崩壊の根本原因)**: Tn 脚順 (l,t,r,b,p) の下で、graded
カーネルの canonical ボンド方位は**ラスタ順(左→右、上→下)**: ゲートの第1脚
(source)は JW 順で先のサイト、すなわち横ボンドは左サイト(source_leg==2)、
縦ボンドは**上**サイト(source_leg==3)でなければならない。source_leg∈{0,1} の
更新は driver(`iTPS::simple_update`)で source/target の役割交換+素行列の
(1,0,3,2) 転置により正規化する。判別診断: 縦鎖 λ 軌跡(t_x=0)が source_leg=1
で case=B(τ スケール、1 ステップ目から)、source_leg=3 で case=A(3e-15)。
逆向きに当てると誤った Koszul マスクが掛かり、2D では市松 CDW への崩壊
(hopping ~1e-3、D=4 で負ノルム)として顕在化していた。
`build_reduced_pair` の (top, bottom)/(left, right) 順序もこのラスタ順と整合
(oracle 固定済み)。3 診断(横 leg2・縦 leg1・縦 leg3)すべて case=A を確認。

## 改訂4(2026-08-18): d=4(電子系)検証と blob 演算子規約の d 一般化

**検証**(コミット d365659a, 9fc1f46c): physical_dim=4(parity [0,1,1,0])を
Fock 厳密審判(2〜6モード JW 構成器をテスト内に実装)で層別検証。

- f-プリミティブ層・カーネル(横縦・ループ・無切断域): d=4 で機械精度厳密。
- **JW 双子 λ 判定は d=4 では無効**: 奇セクター多重度 ≥2 では成分の素コピーが
  JW 対応にならない(脚順 (l,r,p)↔(l,p,r) の並べ替えマスク (-1)^{|r||p|} が
  ゲージ吸収不能)。d=4 で case=B が出てもバグの証拠にならない。診断は
  Fock 審判(プラケット3者比較・鎖 θ 成分照合)を使うこと。
- **blob 測定の演算子規約は「in+out 両swap」が正**(`wrap_reduced_pair_op`)。
  d=2 の N 保存演算子では素ロードと縮退(R3 が「素」と固定したのは縮退の
  ため)。分岐チャネル (o,o)→(e,e) は奇多重度 2 で初出現(電子ホッピングの
  線形項)。R5 テスト(r2_convention.cpp、常時実行)が恒久固定。

**演算子規約の最終表(3経路)**:

| 経路 | 規約 | ヘルパ |
|---|---|---|
| simple update カーネル(ゲート) | in-swap | `wrap_twosite_op` |
| reduced pair blob(2サイト測定) | in-swap + out-swap | `wrap_reduced_pair_op` |
| 1サイト(ゲート・測定) | 素 | — |

**電子系の使用法**: physical_dim=4、parity=[0,1,1,0]、局所基底
|0>,|↑>,|↓>,|↑↓>(|↑↓>=c†_↑c†_↓|0>)。演算子・ゲート行列は順序付き Fock 基底で
JW 符号込みで計算する(2サイト4モード構成器: テスト内 electron_bond_hamiltonian /
work/electron-validation/ の numpy 版)。統合検証: U=0 電子鎖 1D 極限
E/site=−1.193(D=4、厳密 −4/π 比 6.3%)、n=0.9998。
