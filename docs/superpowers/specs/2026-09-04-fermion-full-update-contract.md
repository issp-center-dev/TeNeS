# フェルミオン full update 振る舞い契約書

日付: 2026-09-04
ブランチ: `fermion`
対応設計書: `docs/superpowers/specs/2026-09-04-fermion-full-update-design.md`

この文書は **テスト作成者への唯一の入力**である。テストはこの契約書だけから書く。
設計書の内部実装の記述(どの中間テンソルをどの順で作るか)は**テストの根拠にしてはならない** —
実装と同じ手順を写したテストは、実装が間違っていても緑になる。

## 0. テスト作成者への注意

- **実装はまだ存在しない。** テストは最初 RED になるのが正しい。コンパイルエラーで落ちるのではなく、
  「関数が無い」以外の理由(値が違う、例外が出ない)で落ちる形が望ましいが、新規関数を呼ぶテストは
  リンクエラーになるのが自然なので、そこは気にしなくてよい。
- **前提条件は必ずアサートする。** 「この初期状態では奇セクタ振幅が立つ」「ボンド Schmidt rank は 1」
  のような、テストの検出力を支える前提は、テスト自身が明示的に検査すること。
  前提が崩れて検査が空洞化しても緑のままになるのを防ぐ。
- **参照値は導出する。** 「実装を走らせて出た値を焼き込む」ことは禁止。厳密縮約、既存の独立経路、
  解析値のいずれかから作ること。焼き込みが避けられない場合(Fock oracle アンカーなど)は、
  生成スクリプトとその実行方法をコメントに残すこと。
- **契約に誤りを見つけたら、テストを実装に合わせず、誤りとして報告すること。**
- 実験・一時ファイルは `work/fermion/full-update-design/` に置く。リポジトリ直下や build ツリーで
  バイナリを実行しない。

## 1. テストの置き場所と登録

C++ テストは doctest で、`test/test_fermion_layer.cpp` が単一バイナリにまとめている。
新しいテストは:

1. `test/fermion/full_update_env.cpp` と `test/fermion/full_update_bond.cpp` を新規作成する
   (GPL v3 のライセンスヘッダを既存ファイルから複写する)。
2. `test/test_fermion_layer.cpp` の末尾(現在 3932 行の `#include "fermion/fold_geometry.cpp"` の次)に
   `#include "fermion/full_update_env.cpp"` と `#include "fermion/full_update_bond.cpp"` を足す。
   `test/CMakeLists.txt` の変更は不要(`test_fermion_layer` は単一ソースをビルドしている)。
3. これらのファイルは `#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN` を持たない。`TEST_CASE` だけを書く。

Python E2E テストは `test/fermion/<name>.py.in` を作り、`test/CMakeLists.txt` に
`configure_file` + `add_test` + `set_tests_properties(... TIMEOUT ...)` を追加する
(既存の `FreeFermionMF` などの記述を手本にする)。

## 2. テストが呼ぶ公開 API

実装者はこの署名どおりに作る。テストはこの署名に対して書く。
名前空間は `tenes::fermion`(§2.1、2.2)と `tenes::itps::core`(§2.3)。

### 2.1 開放チャネル環境(`src/fermion/full_update_env.hpp`、新規)

```cpp
namespace tenes::fermion {

//! Result of the open-channel fold: the two-site environment of a full-update
//! bond, in both the graded (wrap-form) and the plain (ALS) convention.
template <class tensor>
struct full_update_environment {
  //! Graded coefficient tensor, legs (in_A, in_B, out_A, out_B).
  //! For any parity-even operator O in wrap_twosite_gate form,
  //! <O> equals the plain elementwise sum over N.t and O.t.
  ftensor<tensor> N;
  //! N with the elementwise mask (-1)^{p(in_A) p(in_B)} applied.
  //! This is the metric the bosonic ALS core consumes.
  tensor N_plain;
  //! max |forbidden block| / max |all| before the block was projected out.
  double forbidden_ratio;
};

//! Build the two-site environment for a full-update bond by folding the two
//! environment-side QR factors with the open operator channel.
//!
//! QA and QB are the environment-side factors of the graded QR of the two
//! site tensors (rank 4: three physical-lattice legs plus the internal leg).
//! The environment arguments follow the same window convention as
//! contract_reduced_pair_halves_density_CTM().
//!
//! Throws std::runtime_error if the forbidden parity block of the folded
//! coefficient tensor exceeds 1e-8 relative to its largest element.
template <class tensor>
full_update_environment<tensor> build_full_update_environment(
    const tensor& C1, const tensor& C2, const tensor& C3, const tensor& C4,
    const tensor& eT1, const tensor& eT2, const tensor& eT3, const tensor& eT4,
    const tensor& eT5, const tensor& eT6,
    const ftensor<tensor>& QA, const ftensor<tensor>& QB,
    reduced_pair_direction direction);

}  // namespace tenes::fermion
```

`QA` / `QB` の脚順(direction 別):

| direction | QA | QB |
|---|---|---|
| horizontal | (l, t, b, a) | (t, r, b, β) |
| vertical | (l, t, r, a) | (l, r, b, β) |

環境引数 `C1..C4, eT1..eT6` の選び方は §3.3 に書く。

### 2.2 既存関数の切り出し(挙動不変、`reduced.hpp` / `reduced_measure.hpp`)

```cpp
namespace tenes::fermion {

//! The second half of build_reduced_pair_halves(): fold the two sites with
//! already-factorized gate halves. u carries legs (in_A, out_A, k_A) and vt
//! carries (k_B, in_B, out_B); k_A and k_B are independent dimensions.
template <class tensor>
reduced_pair_halves<tensor> build_reduced_pair_halves_from_factors(
    const ftensor<tensor>& TnA, const ftensor<tensor>& TnB,
    const ftensor<tensor>& u, const ftensor<tensor>& vt,
    reduced_pair_direction direction);

//! Absorb the two halves into their CTM environment, stopping before the
//! join. Both outputs carry legs (boundary, k, boundary).
template <class tensor>
void absorb_reduced_pair_halves(
    const tensor& C1, const tensor& C2, const tensor& C3, const tensor& C4,
    const tensor& eT1, const tensor& eT2, const tensor& eT3, const tensor& eT4,
    const tensor& eT5, const tensor& eT6,
    const reduced_pair_halves<tensor>& halves, tensor& left, tensor& right);

}  // namespace tenes::fermion
```

### 2.3 ボンド更新(`src/iTPS/core/full_update_fermion.cpp`、新規)

```cpp
namespace tenes::itps::core {

//! One full-update bond in fermion mode.
//!
//! Tn1 is the raster-earlier site (left for horizontal, top for vertical) and
//! Tn2 the later one. wrapped_gate must already be in wrap_twosite_gate form
//! with the site roles matching Tn1/Tn2 — the caller does the wrapping and the
//! source swap; this function does neither.
template <class tensor>
void Full_update_bond_fermion(
    const tensor& C1, const tensor& C2, const tensor& C3, const tensor& C4,
    const tensor& eT1, const tensor& eT2, const tensor& eT3, const tensor& eT4,
    const tensor& eT5, const tensor& eT6,
    const tenes::fermion::ftensor<tensor>& Tn1,
    const tenes::fermion::ftensor<tensor>& Tn2,
    const tenes::fermion::ftensor<tensor>& wrapped_gate,
    tenes::fermion::reduced_pair_direction direction,
    const PEPS_Parameters& peps_parameters,
    tenes::fermion::ftensor<tensor>& Tn1_new,
    tenes::fermion::ftensor<tensor>& Tn2_new);

}  // namespace tenes::itps::core
```

`Tn1_new` / `Tn2_new` は rank 5、脚順 (l, t, r, b, s)。台帳は接続ボンド脚だけが更新され、
他の脚は入力と同じ台帳を引き継ぐ。

### 2.4 既存 API(テストが使う)

- `tenes::fermion::ftensor<tensor>` — メンバ `tensor t`、`leg_parities parity`
  (`src/fermion/ftensor.hpp:83`)。`parity_vector = std::vector<bool>`、
  `leg_parities = std::vector<parity_vector>`(`src/fermion/parity.hpp:41,43`)。
- `wrap_twosite_gate(op, p1, p2)` → 脚 (in1, in2, out1, out2)(`fops.hpp:317`)。
- `build_pair_state(TnA, TnB, direction)`、`apply_pair_op(pair, op12)`(`reduced.hpp:379,408`)。
- `build_reduced_pair_halves(TnA, TnB, op12, direction)`(`reduced.hpp:529`)、
  `contract_reduced_pair_halves_density_CTM(C1..C4, eT1..eT6, halves)`
  (`reduced_measure.hpp:396`)。
- `qr(a, rows, cols, q, r)`、`svd_trunc(a, rows, cols, u, s, vt, dc)`、
  `tensordot`、`transpose`、`conj`、`parity_violation`、`max_abs`(`fops.hpp`)。
- `build_reduced_density_tensors(Tn, finfo)`、
  `tenes::itps::core::Calc_CTM_Environment_density(...)`。
- `Create_Environment_two_sites(C1..C4, eT1..eT6, Q1, Q2)`(`full_update.cpp:49`、
  bosonic、`full_update.hpp` に宣言を追加してもらう)。
- `test/fermion/fold_geometry.cpp` の既存治具(有限パッチの厳密縮約、Fock oracle アンカー)。
  同じ doctest バイナリに入るので、そこで定義済みの補助関数は再利用できる。

## 3. 振る舞い契約

### 3.1 用語

- **wrap 形式の二サイト演算子**: `wrap_twosite_gate` が返す rank-4 の graded テンソル、
  脚 (in1, in2, out1, out2)。in が ket 側、out が bra 側。
- **pair state**: `build_pair_state(TnA, TnB, direction)`。接続ボンドを縮約し、
  外側の仮想脚と 2 つの物理脚を開いたままにした rank-8 の graded テンソル。
  接続ボンドのゲージに依存しない。
- **forbidden block**: パリティが偶にならない要素の集まり。偶テンソルでは恒等的にゼロ。
- **相対誤差**: 特記なき限り `max_abs(差) / max(1, max_abs(参照))`。

### 3.2 N の定義的性質(T2)

`build_full_update_environment` が返す `N` は次を満たす。**これが N の定義であり、
実装手順ではない。**

**(T2-i) 開放版と閉鎖版の一致。** `QA`, `QB` を擬似サイトに詰め替えたもの
(内部脚 a / β を物理スロットに置き、接続ボンドスロットに次元 1 の偶ダミーを入れた rank-5)
を `QA′`, `QB′` とする。任意のパリティ偶な wrap 形式演算子 O(脚次元は a, β に合わせる)について

    Σ_x N.t(x) · O.t(x)  ==  contract_reduced_pair_halves_density_CTM(
                                 env, build_reduced_pair_halves(QA′, QB′, O, direction))

が成り立つ。左辺は **plain な要素積和**であり、`fermion::trace` ではない
(`fermion::trace` は奇脚数 4 に対する追加符号を持つので使ってはならない)。

覆う条件: direction は horizontal と vertical の両方、`real_tensor` と `complex_tensor` の両方、
仮想脚台帳は非自明(偶奇混在)なもの、物理次元は d=2 と d=4 の両方、そして
**内部脚次元が左右で異なる形状(nA ≠ nB)を最低 1 ケース**含める。
許容: 相対 1e-12。

**(T2-ii) forbidden block。** 返り値の `forbidden_ratio` が 1e-10 以下であること。
また、環境をわざと壊した(たとえば eT のどれかに偶でない成分を足した)入力では
`build_full_update_environment` が `std::runtime_error` を投げること。

**(T2-iii) N_plain の性質。** `N_plain` はエルミート:
`N_plain(a,β,a′,β′) == conj(N_plain(a′,β′,a,β))`、相対 1e-12。
また半正定値: `eigh` の最小固有値が `-1e-12 * 最大固有値` 以上。

**(T2-iv) 独立 oracle。** 有限パッチの厳密縮約(CTM ではない)を環境として、
`Σ N·O` が独立に計算した真値と一致すること。真値は
`test/fermion/fold_geometry.cpp` にある「厳密に縮約した環境で閉じる」経路、
または Fock oracle から作る。**これが N の符号を検査する唯一のテストである**(§4 参照)ので、
仮想脚台帳が非自明、物理次元 d=4、horizontal と vertical の両方を必ず覆うこと。
許容: 相対 1e-12。

**(T2-v) bosonic 退化。** すべての台帳が偶のとき、`N` は bosonic の
`Create_Environment_two_sites(C1..C4, eT1..eT6, Q1, Q2)` と一致する。
ただし fold 環境の eT は (ket, bra) を融合した D² 脚を持つので、bosonic 関数に渡すには
(ket, bra) の順に reshape して分ける必要がある。
これは **bosonic 回帰**であり、フェルミオン符号の検出力は無いものとして扱う。
許容: 相対 1e-12。

### 3.3 環境窓の選択

`build_full_update_environment` と `Full_update_bond_fermion` に渡す 10 個の環境テンソルは、
測定経路(`src/iTPS/twosite_obs.cpp`)と同じ規約で選ぶ。s1 は raster で先の site
(水平なら左、垂直なら上)、s2 は後の site。

| direction | C1, C2, C3, C4 | eT1, eT2, eT3, eT4, eT5, eT6 |
|---|---|---|
| horizontal | C1[s1], C2[s2], C3[s2], C4[s1] | eTt[s1], eTt[s2], eTr[s2], eTb[s2], eTb[s1], eTl[s1] |
| vertical | C1[s1], C2[s1], C3[s2], C4[s2] | eTt[s1], eTr[s1], eTr[s2], eTb[s2], eTl[s2], eTl[s1] |

### 3.4 ボンド更新の性質(T3)

**(T3-i) 厳密性。** 接続ボンドを横切る Schmidt rank が 1 の初期状態と、
適用後の厳密な pair state の Schmidt rank が D 以下に収まるゲートを与えたとき、
更新後の pair state は厳密解と一致する:

    normalize(build_pair_state(Tn1_new, Tn2_new, direction))
      ==  normalize(apply_pair_op(build_pair_state(Tn1, Tn2, direction), wrapped_gate))

比較は全体ノルムで正規化し、全体位相(実テンソルなら符号)を許す。許容: 相対 1e-8。

このテストが空洞化しないための **必須の前提アサーション**:

1. 初期 pair state の接続ボンドを横切る Schmidt 値が 1 本だけ非ゼロであること。
2. 厳密解(右辺)の接続ボンドを横切る Schmidt rank が D 以下であること。
3. **厳密解が偶・奇どちらのボンドパリティセクタにも非ゼロ振幅を持つこと。**
   これが無いと奇チャネルが立たず、符号の検査にならない。占有数の固有状態どうしの
   単純な積(たとえば両サイトとも空)ではホッピングがゼロ作用になるので不可。
   スピンレスなら「片方が占有、もう片方が空」の重ね合わせ、d=4 なら doublon-holon
   チャネルが立つ配置を選ぶ。
4. `Inverse_Env_cut` と `Full_Inverse_precision` を、必要な方向を落とさない十分小さい値
   (たとえば 1e-14)に設定していること。環境が対象部分空間で正定値であること
   (最小固有値を検査する)。

horizontal と vertical の両方で行う。

**(T3-ii) 恒等ゲート。** 恒等ゲートを与えると pair state は(正規化と位相を除いて)不変。
両方向。許容: 相対 1e-8。

**(T3-iii) パリティ。** `Tn1_new`, `Tn2_new` の `parity_violation` が
`1e-12 * max_abs` 以下。接続ボンドの台帳が両端で一致すること
(`Tn1_new.parity[接続脚] == Tn2_new.parity[対向脚]`)。

**(T3-iv) bosonic 退化。** すべての台帳が偶のとき、bosonic の
`Full_update_bond`(`source_leg` は水平なら 2、垂直なら 3)と pair state が一致する。
ゲージの違いを吸収するため Tn そのものではなく pair state で比較する。許容: 相対 1e-8。

**(T3-v) ゲージ固定なし。** `Full_Gauge_Fix = false` でも (T3-i) が通る。

**(T3-vi) 病的な環境。** 次の 2 つの人工的な環境で、例外を投げずに完走し、
出力の `parity_violation` が (T3-iii) の許容内に収まること:

- 偶セクタと奇セクタで固有値が厳密に縮退している `N_plain`。
- `Inverse_Env_cut` でランク落ちする(固有値のいくつかが厳密にゼロの)`N_plain`。

これらは `build_full_update_environment` を経由せず環境を直接与える必要があるので、
テストからは環境テンソル側を構成して作る。難しければ、この 2 件は
`Full_update_bond_fermion` を通さず、実装が公開するゲージ固定部分
(設計書 §4.4 の `prepare_environment`)に対する単体テストにしてよい。

**(T3-vii) 変異防衛。** 次の変異をそれぞれ 1 つだけ加えた実装のコピーで (T3-i) が
**赤になる**ことを、レビュー時に手で確認する(テストコードには含めない):
`mask_{m s1}` の除去、Θ の graded 合成を plain に置換、Q′ 詰め替えの台帳の破壊、
`Tn_new` の転置の入れ替え。
**N_plain の入力脚マスクの除去はここでは赤にならない**(§4)。

### 3.5 driver とガード(T4, T5, T6)

**(T4) bosonic 同値。** `test/data/AntiferroHeisenberg_real.toml` を元にした入力を、
(a) fermion なし、(b) `fermion = true` かつ全サイトの物理パリティを `[0, 0]`(全偶)にした版、
の 2 通りで走らせ、エネルギー・onesite 観測量・twosite 観測量が一致すること。
許容: 相対 1e-6(fold CTM と素 CTM の収束差を許容する)。

**入力の作り方(2026-09-04 訂正)**:

- `fastfullupdate = false` は **`[parameter.full_update]`** に書く
  (`[parameter.ctm]` ではない。`load_toml.cpp:467` は `[parameter.full_update]` から読む)。
- `validate_fermion_constraints`(`load_toml.cpp:606-654`)が拒否するものを取り除く:
  `[correlation]`(`r_max > 0`)、`[[observable.multisite]]`、`ops = [...]` 形式の
  二サイト観測量(明示 `elements` に書き下す)、`Simple_Gauge_Fix`、`Use_RSVD`。
- **`noise = 0.0` にする(必須)。** 仮想台帳は物理パリティと無関係に
  `even_first_parity(D)`(`tensors.cpp:80-90`)で初期化され、fermion 初期化は総パリティが奇の
  Tn 要素をゼロにする(`tensors.cpp:155-159`)。したがって `noise != 0` では 2 つの run が
  **違う初期テンソルから始まり**、どんな許容値でも一致しない。`noise = 0.0` なら非ゼロ要素が
  仮想インデックス `(0,0,0,0)`(どの台帳でも偶)だけになり初期テンソルがビット一致する。
- この前提(初期テンソル一致、および走行後の台帳が全偶であること)を、テスト自身が
  `fermion.dat` を読んで明示アサートすること。

**(T5) 物理。** 自由フェルミオン(既存 `test/fermion/free_fermion*.py.in` の模型)を
**D=2、小さい CHI、`Max_CTM_Iteration` を明示的に小さく**した設定で、
simple update のあとに full update を数ステップ走らせる。要求:

- エネルギーが増えない: `E_FU <= E_SU + 1e-8`。
- 厳密値に近づく: `|E_FU - E_exact| < |E_SU - E_exact|`。
- 垂直ボンドを必ず通る配置にし、`source_leg` が 0 または 1 になるゲートも
  最低 1 つ含める(driver の raster 正規化を通すため)。
- ctest 登録時に明示 `TIMEOUT` を付ける。既存の `FreeFermionMF`(1800 秒)を上限の目安にする。
  非高速版はボンドごとに CTM を再収束するので、`[[evolution.full]]` の数だけ
  1 ステップあたり CTM 収束が走る。

**(T6) ガードとフォールバック。**

- `fermion = true` かつ full update のステップ数が正の入力が **受理される**
  (従来はここで `input_error` を投げていた)。
- `fermion = true` かつ `fastfullupdate = true`(既定値)の入力は、
  **エラーにならず完走**し、標準エラー出力に `Full_Use_FastFullUpdate` を含む警告が出る。
- `fermion = true` かつ `meanfield_env = true` かつ full update のステップ数が正の入力は
  従来どおり拒否される。
- `source_leg` が 0 または 1 のボンドを含む fermion + full update の実行後、
  仮想パリティ台帳が全ボンドで両端一致していること
  (`tenes::fermion::validate_neighbor_consistency` が投げない)。

**(T7) bosonic 回帰。** 既存テストが無変更で緑であること。新規テストは不要。
対象: doctest `test_full_update`、golden の `AntiferroHeisenberg_real` /
`AntiferroHeisenberg_complex` / `AntiferroHeisenberg_mf` / `Honeycomb` / `J1J2_AFH` /
`RSVD` / `Kitaev`、および既存フェルミオンテスト `FreeFermion` / `FreeFermionMF` /
`FreeFermionSaveLoad` / `FreeFermionSimple` / `test_fermion_layer`。

## 4. 検出力の限界(重要)

**T3-i の厳密性検査は N の符号を検出しない。**
費用関数は ‖ψ(R1,R2) − ψ(Θ)‖²_N である。目標 Θ がボンド次元 D で厳密に表現できるとき、
最小値 0 は N が正定値でありさえすれば計量によらず ψ = Θ で達成される。
したがって N の入力脚マスク・開放 join・crossing mask・transpose 符号を誤っても、
壊れた N がなお正定値である限り T3-i は緑のままになる。

役割分担:

**T4 も符号バグを検出しない(2026-09-04 追記)。** T4 は台帳が全偶なので、
入力脚マスク・crossing mask・graded transpose の符号がすべて +1 に退化する。
T4 が捕まえるのは構造のバグ(driver の配線、窓環境の選択、台帳更新、次元)だけである。

| 壊した箇所 | 検出するテスト |
|---|---|
| N の入力脚マスク、開放 join、crossing mask、N の transpose 符号 | **T2-iv**(独立 oracle)。E2E では **T5 のみ** |
| Q′ 詰め替え、QR 軸 | T2-iv、T3-i |
| Θ の graded 合成、`mask_{m s1}`、初期推定の転置、`Tn_new` の転置 | **T3-i** |
| ゲージ因子のパリティ漏れ | T3-iii、T3-vi |
| driver の配線、窓環境の選択、台帳更新 | T4、T5、T6 |

したがって **T2-iv は省略できない**。ここを共有コード経由の自己整合だけで済ませると、
N の符号バグは T5 まで誰にも捕まらない(そして T5 が落ちても原因の切り分けができない)。

## 5. 契約書チェックリスト(作成者が自分で確認する)

- [ ] 参照値は厳密縮約・独立経路・解析値のいずれかから作った(実装出力の焼き込みでない)
- [ ] 前提条件(Schmidt rank、奇セクタ振幅、環境の正定値性)を明示アサートした
- [ ] 反対称な演算子の符号を、期待値の側でも正しく扱った
- [ ] horizontal と vertical、real と complex、d=2 と d=4 を必要な箇所で覆った
- [ ] `nA != nB` になる形状を T2 に入れた
- [ ] 検査が原理的に不可能なもの(§4 の限界)をテストに要求していない
- [ ] テストは `work/` 以下でのみ一時ファイルを作り、リポジトリ直下を汚さない
