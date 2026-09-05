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
- **テンソルの要素比較は global index 経由で行う(2026-09-04 セッション 2 追加)。**
  `mptensor::Tensor` は `transpose` を遅延評価する(axes_map)ので、shape が同じでも 2 つのテンソルの
  local index `n` が指す global index は一致しない。`a[n]` と `b[n]` を直接比べると、要素の最大値は
  一致するのに並びだけずれて相対差 O(1) に見える。必ず `a.global_index(n)`(または
  `global_index_fast`)で index を取り、`b.get_value(index, v)` で対応する要素を読む。
  本番でのバグ探しがこれで丸 1 セッション空費した。

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

**2026-09-05 追加分(T8、T2-vi complex、T5 改訂)の置き場所**:

- T8-i / T8-ii / T8-iv と T2-vi complex は新規 `test/fermion/ctm_phase.cpp` に書き、
  `test/test_fermion_layer.cpp` の末尾に `#include "fermion/ctm_phase.cpp"` を 1 行足す。
  既存の `full_update_env.cpp` / `full_update_bond.cpp` / `full_update_realctm.cpp` は変更しない
  (治具の再利用は `#include` 経由でなく、同じ翻訳単位にあるので直接呼べる)。
- T8-iii と T8-v は新規 `test/fermion/free_fermion_complex.py.in`(ctest 名 `FreeFermionComplex`)。
- T5 改訂は既存 `test/fermion/free_fermion_full.py.in` を書き換える(このファイルだけは変更を許可する)。
  ctest 名 `FreeFermionFull` と TIMEOUT は `test/CMakeLists.txt` で調整してよい。

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
  //! Phase removed from N so that the window norm sum_x N(x) I_wrap(x) is
  //! real positive (exactly 1 when nothing was changed). Added 2026-09-05,
  //! see §2.5 and §3.6.
  std::complex<double> phase;
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

### 2.5 CTM 環境の位相(2026-09-05 追加、改訂 2。設計書 `2026-09-05-fermion-ctm-phase-design.md`)

```cpp
namespace tenes::itps::core {

// 既存関数に引数を足す(既定値 false。fermion の update_CTM() だけが true を渡す)
template <class tensor>
bool Check_Convergence_CTM_RDM(
    const std::vector<tensor>& C1, const std::vector<tensor>& C2,
    const std::vector<tensor>& C3, const std::vector<tensor>& C4,
    const std::vector<tensor>& eTt, const std::vector<tensor>& eTr,
    const std::vector<tensor>& eTb, const std::vector<tensor>& eTl,
    const std::vector<tensor>& Tn, const SquareLattice lattice,
    std::vector<small_tensor<typename tensor::value_type>>& rdm_old,
    bool& has_rdm_old, const double epsilon, const bool is_density,
    double& rdm_dist, const bool phase_invariant = false);

template <class tensor>
int Calc_CTM_Environment_density(
    std::vector<tensor>& C1, std::vector<tensor>& C2, std::vector<tensor>& C3,
    std::vector<tensor>& C4, std::vector<tensor>& eTt, std::vector<tensor>& eTr,
    std::vector<tensor>& eTb, std::vector<tensor>& eTl,
    const std::vector<tensor>& Tn, const PEPS_Parameters peps_parameters,
    const SquareLattice lattice, bool initialize = true,
    bool phase_invariant = false);   // 戻り値は反復回数

//! Divide a gathered one-site RDM by the phase of its trace so that the
//! trace becomes real positive. |trace| and the Frobenius norm are kept.
//! Returns the phase that was removed (exactly 1 when nothing was changed:
//! |phase - 1| <= 1e-14 leaves the matrix untouched).
//! Throws std::runtime_error if any element is not finite or the trace vanishes.
template <class value_type>
std::complex<double> normalize_rdm_phase(small_tensor<value_type>& rdm);

}  // namespace tenes::itps::core

namespace tenes::fermion {

// §2.1 の構造体にメンバを足す
template <class tensor>
struct full_update_environment {
  ftensor<tensor> N;        // 位相正規化済み: Σ_x N.t(x)·I_wrap.t(x) が実正
  tensor N_plain;           // 正規化済み N から作る
  double forbidden_ratio;
  std::complex<double> phase;  // N から取り除いた位相(何もしなければ exactly 1)
};

}  // namespace tenes::fermion
```

`build_full_update_environment` は、開放 join と forbidden 検査の後、その窓の
`norm = Σ_x N.t(x) · I_wrap.t(x)`(`I_wrap = wrap_twosite_gate(δ_{in_A,out_A} δ_{in_B,out_B}, p_a, p_β)`、
plain 要素積和)の位相 `norm/|norm|` で N を割ってから返す。`N.t` に非有限な要素があるか、`norm` が
非有限または実質ゼロ(`|norm| ≤ 1e-12·max_abs(N)·nA·nB`)なら `std::runtime_error`。
`small_tensor` は `src/tensor.hpp` の非分散テンソル。

改訂 1 にあった `fix_environment_phase` は**存在しない**(環境の位相は窓ごとに異なるので、
環境を一括で回すことはできない。設計書 §1・§2)。

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

**2026-09-05 注記**: 右辺(閉じた測定経路)は環境の位相を含み、左辺の N は §2.5 の位相正規化を
経ているので、実 CTM 環境では両辺は `phase`(返り値)倍だけ異なる。T2-i は有限パッチの環境
(位相 1)で要求する。実 CTM での N の性質は T2-vi と T8-ii で扱う。

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

**(T2-vi) 実 CTM 環境での計量検証(2026-09-04 セッション 2 追加)。** T2-i は共有コードの
自己整合、T2-iv は有限パッチ oracle であり、`Calc_CTM_Environment_density` が作る**実際の CTM 環境**を
通すテストが無かった。本番でのバグ探しを 1 セッション長引かせた穴なので、次を要求する。

環境: `build_reduced_density_tensors(Tn, finfo)` → `Calc_CTM_Environment_density` が作る**実際の CTM
環境**(`iTPS::update_CTM()` と同じ経路)。Tn は次のいずれか:

- 自由フェルミオン(`test/fermion/free_fermion*.py.in` の模型、D=2、chi=8)を simple update で
  数百ステップ収束させた状態(E2E で `tensor_save` したテンソルと `fermion.dat` を読んでもよい)。
- それが難しければ、**パリティ偶なランダム Tn**(2×2 セル、D=2、d=2、偶奇混在台帳)でよい。
  計量検証としての検出力は状態に依らない。

いずれの場合も、CTM が収束していること(`Calc_CTM_Environment_density` の収束、および
`build_full_update_environment` の `forbidden_ratio ≤ 1e-10`)を**前提アサート**する。
complex_tensor は §3.6 の (T2-vi complex) で扱う(環境の位相固定が要る)。
CTM 未収束のパリティ漏れは N_plain の誤差として現れ、比較の許容を破るのは正しい挙動だが、
それをテストの失敗と区別できるようにしておく。
§3.3 の窓選択で 10 個の環境テンソルを取り、`fermion::qr` で QA, RA, QB, RB を作る(bosonic と同じ軸):

| direction | site 1 | site 2 |
|---|---|---|
| horizontal | `qr(Tn1, Axes(0,1,3), Axes(2,4), QA, RA)` → QA(l,t,b,a), RA(a,r,s1) | `qr(Tn2, Axes(1,2,3), Axes(0,4), QB, RB)` → QB(t,r,b,β), RB(β,l,s2) |
| vertical | `qr(Tn1, Axes(0,1,2), Axes(3,4), QA, RA)` → QA(l,t,r,a), RA(a,b,s1) | `qr(Tn2, Axes(0,2,3), Axes(1,4), QB, RB)` → QB(l,r,b,β), RB(β,t,s2) |

パリティ偶なブロック Y(a, β, s1, s2)(台帳 {p_a, p_β, p_s1, p_s2}、p_a / p_β は QR の内部脚台帳)
について:

1. Y を `fermion::svd(Y, Axes(0,2), Axes(1,3), U, s, VT)`(打ち切りなし)で分け、`U` に s を
   軸 2 で掛けたものを R1(a, s1, m)、`transpose(VT, Axes(1,2,0))` を R2(β, s2, m) として、
   次の式(すべて graded 演算)で Tn1(Y), Tn2(Y) を組み立てる:

   | direction | Tn1(Y) | Tn2(Y) |
   |---|---|---|
   | horizontal | `transpose(tensordot(QA, R1, Axes(3), Axes(0)), Axes(0,1,4,2,3))` | `transpose(tensordot(QB, R2, Axes(3), Axes(0)), Axes(4,0,1,2,3))` |
   | vertical | `transpose(tensordot(QA, R1, Axes(3), Axes(0)), Axes(0,1,2,4,3))` | `transpose(tensordot(QB, R2, Axes(3), Axes(0)), Axes(0,4,1,2,3))` |

   **前提アサーション**: Y = X(下記)のとき、この組み立てで得た pair state が
   `build_pair_state(Tn1, Tn2, direction)` と(全体スケールを除き、§0 のとおり global index 経由で)
   一致すること。これが通らなければ組み立て式の写し間違いであり、以降の比較は無意味になる。
2. 測定経路で ⟨ψ(Y)|ψ(Y)⟩ =
   `contract_reduced_pair_halves_density_CTM(env, build_reduced_pair_halves(Tn1(Y), Tn2(Y), I, direction))`
   (I は `wrap_twosite_gate(恒等)`)と、同じ形で ⟨ψ(Y)|G|ψ(Y)⟩(G は wrap 済みの偶ゲート。
   hopping の虚時間ゲートでよい)を計算する。
3. Ỹ = mask_{aβ} ⊙ Y(a と β がともに奇の要素に −1)、Θ_Y = `fermion::tensordot(Y, G, Axes(2,3), Axes(0,1))`、
   Θ̃_Y = mask_{aβ} ⊙ Θ_Y として、bosonic 規約の二次形式
   ⟨Ỹ|N_plain|Z̃⟩ = `trace(tensordot(N_plain, Z̃, Axes(0,1), Axes(0,1)), conj(Ỹ), all, all)`
   で ⟨Ỹ|N_plain|Ỹ⟩ と ⟨Ỹ|N_plain|Θ̃_Y⟩ を計算する。
4. 2 と 3 が一致すること。ノルムは相対 1e-10、⟨G⟩ は比 ⟨Y|G|Y⟩/⟨Y|Y⟩ の形で 1e-10。

Y は次をすべて含める: X = `transpose(tensordot(RA, RB, Axes(1), Axes(1)), Axes(0,2,1,3))`(元の状態)、
Θ = X·G、パリティ偶なランダム 2 種、X + 0.05·ランダム。horizontal と vertical の両方。

検出するもの: N_plain の符号(入力脚マスク、開放 join、crossing mask、transpose)を実 CTM 環境で、
かつ ALS が探索する**一般のブロック**で。X と Θ の組だけでは ⟨X|Θ⟩ = ⟨G⟩ の 1 点しか見ないので
不十分(本番の最初の診断はここで止まっていた)。

**(T2-vii) 擬似サイトの状態同一性(2026-09-04 セッション 2 追加)。** 同じ環境・同じ QR で、
QA′, QB′(T2-i の詰め替え)と X について

    apply_pair_op(build_pair_state(QA′, QB′, direction), X)  ==  build_pair_state(Tn1, Tn2, direction)

要素ごと(§0 のとおり global index 経由)、相対 1e-12、両方向。設計書 §4.2 C1 の production 版。

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

**(T3-viii) 実 CTM 環境での恒等ゲート不変性(2026-09-04 セッション 2 追加)。** T2-vi と同じ
実 CTM 環境で `Full_update_bond_fermion` を呼ぶ。

- 恒等ゲート: 返った `Tn1_new`, `Tn2_new` の pair state が元の pair state と一致。
  **全体スケールは変わるのが正常**(接続ボンドの λ を Σλ²=1 に正規化するため)なので、
  最小二乗で係数 c を合わせて `|pair_new − c·pair_old| / max|pair_old|` を相対 1e-10 で見る。
  位相(実なら符号)も c に吸収させる。両方向。
- hopping ゲート(τ=0.01 程度): 返った pair state の ⟨G⟩(T2-vi の 2 の形で測定経路から計算)が、
  元の pair state の ⟨G⟩ より大きいこと(そのボンドの局所エネルギーが下がる)。全体エネルギーが
  下がることは D=2 では要求しない(§3.5 T5 の注記)。
  **注意(テスト作成者の指摘、2026-09-04)**: これは打ち切りなしなら log-convexity から従うが、
  D への打ち切り後まで含めた定理ではなく状態依存でありうる。SU 収束状態(1.00360 → 1.00368)と
  ランダム偶状態(+1.6e-5 / +2.5e-5)で成立を確認済みだが、別の seed で落ちた場合は
  実装バグと即断せず、状態依存の可能性から切り分けること。
- パラメータは T3-i と同じ(`Inverse_Env_cut` / `Full_Inverse_precision` = 1e-14 など)。

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

**(T5) 物理(2026-09-05 改訂)。** 自由フェルミオン(既存 `test/fermion/free_fermion*.py.in` の模型)を
**D=3、chi=12、`[parameter.ctm] iteration_max = 200`**、SU 1000 step(τ=0.01)のあとに
full update を 1〜4 sweep 走らせる(4 sweep を推奨)。要求:

- エネルギーが下がる: `E_FU < E_SU − 1e-4`(参考実測 chi=12: E_SU = −0.744499、1 sweep −0.744816、
  4 sweep −0.745614)。
- 厳密値に近づく: `|E_FU − E_exact| < |E_SU − E_exact|`。
- 標準出力に「CTM did not converge」が無い。
- 垂直ボンドを必ず通る配置にし、`source_leg` が 0 または 1 になるゲートも最低 1 つ含める。
- ctest 登録時に明示 `TIMEOUT` を付ける(実測: SU 24 秒 + FU 1 sweep 8 秒。4 sweep で 1 分程度。
  TIMEOUT はその 3 倍以上)。

**D=2 の非増加要求は削除する。** D=2 の自由フェルミオンでは SU 収束状態からの full update が
エネルギーを上げることが確定しており(実装は正しく、D=2 の表現力不足で射影虚時間発展の固定点が
SU 固定点より高い。`work/fermion/full-update-design/FINDINGS-task5-energy.md`)、要求として
成立しない。

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

### 3.6 CTM 環境の位相(T8、2026-09-05 追加、改訂 2)

背景: フェルミオン mode の CTM(`build_reduced_density_tensors` → `Calc_CTM_Environment_density`)は、
fold 済みテンソルを一様ベクトルで潰して初期化するため、収束した環境の C, eT は状態依存の位相を持つ
(実数なら ±1、複素なら任意)。位相は**窓ごとに異なり**、複素では反復ごとに回る
(`work/fermion/ctm-phase/report-testauthor.md` §7.1 の計測)。測定値 ⟨O⟩/⟨1⟩ は位相が打ち消されるので
影響しないが、`Check_Convergence_CTM_RDM` の実正値判定、full update の N のエルミート化・正値化、
`onesite_density_matrix.dat` の出力が壊れる。対処は「判定を位相不変にし、環境を消費する側が
自分の窓の ⟨1⟩ で位相を決める」。設計書 `2026-09-05-fermion-ctm-phase-design.md`、API は §2.5。

**(T8-i) 収束判定の位相不変性。** `core::Check_Convergence_CTM_RDM` を、独立な履歴
(`rdm_old`, `has_rdm_old`)を持つ **2 系列**で **2 回ずつ**呼ぶ:

- 基準系列: 環境 E0 → 環境 E1(E1 は E0 と少し違う環境。たとえば E0 の C1..eT に 1% の摂動を足したもの)。
- 変換系列: e^{iφ0}·E0 → e^{iφ1}·E1(位相は全 `C1[i]` に掛ける。φ0 ≠ φ1。実数テンソルでは
  (−1, +1) と (+1, −1) の両方、複素では任意の 2 値、たとえば 0.7 と 2.9)。

要求(`phase_invariant = true`): 各系列の初回は未収束(戻り値 false、`rdm_dist` は NaN)。
2 回目の戻り値が両系列で一致し、`rdm_dist` が相対 1e-12 で一致し、NaN でない。
`is_density = true`(fold 済み Tn)と `false`(bosonic Tn)の両方、real と complex の両方。

追加要求:
- 健全性: 全 `C1[i]` をゼロにした環境、NaN を 1 要素含む環境、Inf を 1 要素含む環境では、
  2 回目も戻り値 false かつ `rdm_dist` が NaN(`phase_invariant` の真偽どちらでも)。
- `phase_invariant = false`(既定)では、実数で −1 を掛けた系列は 2 回目も未収束のまま(現行挙動)。

**(T8-ii) N の位相不変性。** パリティ偶ランダム Tn(全台帳 {偶,奇} など偶奇混在)から
`build_reduced_density_tensors` → `Calc_CTM_Environment_density(..., phase_invariant = true)`
(`count < iteration_max` を前提アサート)で環境を作り、§3.3 の窓で `build_full_update_environment` を
(a) 環境そのまま、(b) 全 `C1[i]` に e^{iφ}(実数は −1)を掛けた環境、(c) 窓に含まれる eT の 1 つ
(たとえば `eTt[s1]`)に e^{iφ} を掛けた環境、の 3 通りで呼ぶ。状態は最低 3 つ:
2×2 実数 D=3、2×2 複素 D=2、**3×2** 複素 D=2。両方向。要求:

1. (a)(b)(c) の `N.t` と `N_plain` が相対 1e-12 で一致する。
2. 返り値 `phase` は (a) に対して (b) では e^{iφ} 倍、(c) でも e^{iφ} 倍(相対 1e-12)。
   (a) の `phase` は |phase| = 1。
3. `Σ_x N.t(x)·I_wrap.t(x)` が実正(`|Im| ≤ 1e-12·|·|`、Re > 0)。
4. `N_plain` がエルミート(相対 1e-12)かつ半正定値(`eigh` の最小固有値 ≥ −1e-12·最大固有値)。
   複素を含む。
5. `forbidden_ratio` は (a)(b)(c) で一致する(値は丸め水準なので **絶対 1e-12** で判定する)。
6. **独立検証**: 正規化前の窓ノルム(`phase` × 正規化後の `Σ_x N.t(x)·I_wrap.t(x)`)が、閉じた測定経路
   `contract_reduced_pair_halves_density_CTM(env, build_reduced_pair_halves(QA′, QB′, I_wrap, direction))`
   (QA′, QB′ は T2-i の擬似サイト)と相対 1e-12 で一致する。**nA ≠ nB になる形状を最低 1 ケース**含める
   (一様な D・d では nA = nB になるので、窓の 2 サイトの物理次元を変える — たとえば市松に d=2 と d=4 —
   ことで作る。I の脚順・wrap マスク・plain 積和が同じ仕方で誤っていると (1)〜(5) だけでは通ってしまう)。
   **さらに(レビュアーの変異実験、2026-09-05)**: wrap マスクを落とした恒等 I_plain(δδ をそのまま)で
   ノルムを計算する誤りは、Σ N·I_plain と Σ N·I_wrap が同じ位相に実数因子(偶ペア重み − 奇奇重み)を
   掛けたものになるため、典型的な状態(奇奇重みが 2〜3 割)では**どのテストも赤にならない**。
   したがって **奇奇チャネルが支配的な窓を最低 1 ケース**含める: 前提として
   `Re(Σ N·I_plain / Σ N·I_wrap) < 0` を明示アサートする(I_plain は wrap なしの恒等)。
   ランダム偶 Tn の奇パリティ成分を大きくスケールする(たとえば奇ボンド index を持つ要素を 3〜5 倍)
   ことで作れる。このケースでは (6) の独立検証が wrap マスクの欠落を捕まえる。
7. 例外: 全 `C1[i]` をゼロにした環境、窓内の C か eT の 1 要素に NaN を注入した環境、Inf を注入した環境の
   それぞれで `std::runtime_error`。密な縮約では非有限が N の全要素と窓ノルムに伝播するので、
   「N の全要素検査」と「norm の検査」をテストで区別することは要求しない(両方が実装にあることは
   設計書 §3.2 の要求であり、レビューで確認する)。

**(T8-iii) 複素 E2E。** 自由フェルミオン D=2、**mu = 0**(半充填)を `is_real = false` で走らせ
(SU のみ、`num_step` は既存 `FreeFermion` の mu = 0 と同じ)、標準出力・標準エラー出力に「CTM did not converge」
「Norm is negative」「Norm is not real」のいずれも無いこと、エネルギーの虚部が 1e-10 以下、密度が既存
`FreeFermion`(mu = 0)と同じ許容内であること。
mu ≠ 0 は複素の初期状態が別の固定点に落ちるので密度比較の対象にしない
(テスト作成者の観察、`report-testauthor.md` §7.2)。

**(T8-iv) full update の位相不変性。** T3-i の入力で、(a) そのまま、(b) `C1` に −1(複素は e^{iφ}、
φ = 2.0)を掛けたもの、で `Full_update_bond_fermion` を呼び、出力の pair state が一致すること
(全体スケールと位相を最小二乗で吸収、相対 1e-10)。両方向。(b) が例外を投げてはならない。

**(T8-v) one-site RDM の位相正規化。** doctest: `core::normalize_rdm_phase` に、
(a) trace が実正の RDM(ランダムなエルミート正定値行列でよい。**trace の虚部が厳密に 0 になるよう
対角要素は実数で作る** — BLAS の丸めで 1e-17 の虚部が残ると |phase − 1| ≤ 1e-14 の分岐が不定になる)を
渡すと要素が完全一致(no-op)で戻り値が exactly 1、(b) 同じ RDM に既知の位相 e^{iφ}(実数は −1)を掛けたものを渡すと (a) の
RDM と全要素が相対 1e-12 で一致し、戻り値が e^{iφ}(相対 1e-12)、|trace| と Frobenius ノルムが
保存される、(c) ゼロ行列、NaN を非対角にだけ含む行列、Inf を含む行列で `std::runtime_error`。
real と complex。加えて E2E: T8-iii の実行が書く `onesite_density_matrix.dat` の各サイトの
trace が実正(`|Im| ≤ 1e-10·|tr|`、Re > 0)。

**(T2-vi complex)** 既存 T2-vi を `complex_tensor` でも行う。環境は
`Calc_CTM_Environment_density(..., phase_invariant = true)`(`count < iteration_max` を前提アサート)で
作る(位相固定は無い)。測定経路のノルム ⟨ψ(Y)|ψ(Y)⟩ は環境の位相を含むので、比較は
**|ノルム|** どうし(相対 1e-10)と、**比** ⟨ψ(Y)|G|ψ(Y)⟩/⟨ψ(Y)|ψ(Y)⟩ どうし(1e-10)で行う。
`N_plain` 側は位相正規化済みなので ⟨Ỹ|N_plain|Ỹ⟩ は実正。

**(T7 追記、検証者が実施)** bosonic golden(`AntiferroHeisenberg_real/complex/mf`、`RSVD`、`J1J2_AFH`、
`Honeycomb`、`Honeycomb_skew`、`TE_TFI`)と finite-T golden(`FT_TFI_square`、`FT_Kitaev`)は、
本変更の前後のバイナリを同一入力・1 thread・MPI off で走らせて `output_*/*.dat`(`time.dat` を除く)の
**sha256 が一致**すること。基準は `work/fermion/ctm-phase/ab/before/sha256.txt`。
これはテストではなく検証者(Claude)が機械的に確認する。

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
| N の入力脚マスク、開放 join、crossing mask、N の graded transpose の Koszul 符号 | **T2-iv**(独立 oracle)、**T2-vi**(実 CTM・一般ブロック)。E2E では **T5 のみ** |
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
- [ ] テンソルの要素比較は global index 経由で行った(§0)
- [ ] T2-vi の Y に X・Θ 以外の一般ブロック(ランダム、X+0.05·ランダム)を含めた
- [ ] T8-i は独立な履歴を持つ 2 系列で 2 回ずつ呼び、φ0 ≠ φ1 にした
- [ ] T8-ii は 3×2 セルと複素を含め、C1 と eT の両方に位相を掛けた場合の N の不変性を検査した
