# フェルミオン full update 実装計画

> **For agentic workers:** この計画は TeNeS の多段エージェント手順(グローバル CLAUDE.md)で実行する。
> テストは**テスト作成者サブエージェントが契約書から先に書き終えている**前提で、実装者(Codex)は
> テストファイルを一切変更しない。ステップはチェックボックス(`- [ ]`)で追跡する。

**Goal:** フェルミオン模式(`fermion = true`)で full update(`num_full_step > 0`)を動かし、
半充填 Hubbard の Mott 崩壊が simple update の局所最適化に起因するかを検証できるようにする。

**Architecture:** 二サイト環境 N を「開放チャネル fold」で構成し(演算子添字を開いたまま
既存の fold CTM 吸収列に通す)、ALS はマスク済み plain 配列 (N_plain, Θ̃, R̃1, R2) の上で
bosonic のコードをそのまま共有する。graded 演算が要るのは Θ の構成・初期推定・後処理の 3 箇所だけ。

**Tech Stack:** C++17、mptensor、doctest、CMake/CTest、Python 3.9+(E2E テスト)

**Spec:** `docs/superpowers/specs/2026-09-04-fermion-full-update-design.md`
**Contract:** `docs/superpowers/specs/2026-09-04-fermion-full-update-contract.md`

## Global Constraints

- 全 C++ は名前空間 `tenes`(ソルバ中核は `tenes::itps`、graded 層は `tenes::fermion`)。
- 各ソースファイル冒頭に GPL v3 のライセンスヘッダを付ける(既存ファイルから複写)。
- テンソル型はテンプレート。実体化は `real_tensor` と `complex_tensor` の 2 つで、
  既存ファイル末尾の明示的インスタンス化の書式に合わせる。
- **テストファイル(`test/` 以下)は変更禁止。** テストが誤っていると思ったら実装を曲げず、
  BLOCKED として報告する。
- **formatter を実行しない。** 新しい行は周囲のスタイルに手で合わせる。
- **リポジトリ直下や build ツリーで `tenes` やテストバイナリを実行しない。**
  実行が必要なら CWD を `work/fermion/full-update-design/` にする。
- ビルドは `cmake --preset gcc && cmake --build --preset gcc`、テストは `ctest --preset gcc`。
- コミットはしない(sandbox では `.git/index.lock` を作れない)。Claude が代理でコミットする。
- 各タスクの最後に報告ファイルを書く。**報告ファイルは成果物であり、無いこと自体が欠陥。**

## File Structure

| ファイル | 責務 |
|---|---|
| `src/iTPS/core/full_update.hpp` | bosonic の宣言。タスク 1 で `prepare_environment` / `als_iterate` / `Create_Environment_two_sites` を公開 |
| `src/iTPS/core/full_update.cpp` | bosonic の実装。タスク 1 で 3 つに分割(演算順序は不変) |
| `src/fermion/reduced.hpp` | fold の前段。タスク 3 で `build_reduced_pair_halves_from_factors` を切り出す |
| `src/fermion/reduced_measure.hpp` | CTM 吸収。タスク 3 で `absorb_reduced_pair_halves` を切り出す |
| `src/fermion/full_update_env.hpp`(新) | 開放チャネル fold による二サイト環境 N の構成だけ |
| `src/iTPS/core/full_update_fermion.cpp`(新) | フェルミオンのボンド更新カーネル |
| `src/iTPS/core/full_update_fermion.hpp`(新) | 同上の宣言 |
| `src/iTPS/full_update.cpp` | driver。タスク 5 でフェルミオン分岐を追加 |
| `src/iTPS/iTPS.cpp`, `src/iTPS/measure.cpp` | タスク 2 で `update_CTM()` をフェルミオン対応に |
| `src/iTPS/load_toml.cpp`, `src/iTPS/main.cpp` | タスク 5 でガード撤去とフォールバック |

---

## Task 1: bosonic core を 3 つに分割する

**Files:**
- Modify: `src/iTPS/core/full_update.cpp:95-391`(`Full_update_bond_horizontal`)
- Modify: `src/iTPS/core/full_update.hpp`(宣言を追加)
- Test: 既存の `test/full_update.cpp`(変更禁止)と golden

**Interfaces:**
- Produces(タスク 4 が使う):

```cpp
namespace tenes::itps::core {

//! Hermitize, positivize and (optionally) gauge-fix the two-site environment.
//! Environment_in has legs (ket_a, ket_b, bra_a, bra_b); Theta_in has legs
//! (a, b, s1, s2). On return Theta_out is the gauged Theta and LR1_inv /
//! LR2_inv are the factors that remove the gauge afterwards (identity-like
//! no-ops when Full_Gauge_Fix is false).
template <class tensor>
void prepare_environment(const tensor &Environment_in, const tensor &Theta_in,
                         const PEPS_Parameters &peps_parameters,
                         tensor &Environment_out, tensor &Theta_out,
                         tensor &LR1_inv, tensor &LR2_inv);

//! Alternating least squares for the two reduced tensors.
//! R1 has legs (envR1, D_connect, m1) and R2 has (envR2, D_connect, m2),
//! both on input (as the initial guess) and on output.
template <class tensor>
void als_iterate(const tensor &Environment, const tensor &Theta,
                 const PEPS_Parameters &peps_parameters, tensor &R1,
                 tensor &R2);

//! Already defined in full_update.cpp; just add the declaration to the header
//! so the tests and the fermion kernel can compare against it.
template <class tensor>
tensor Create_Environment_two_sites(const tensor &C1, const tensor &C2,
                                    const tensor &C3, const tensor &C4,
                                    const tensor &eT1, const tensor &eT2,
                                    const tensor &eT3, const tensor &eT4,
                                    const tensor &eT5, const tensor &eT6,
                                    const tensor &Q1, const tensor &Q2);

}  // namespace tenes::itps::core
```

- [ ] **Step 1: 分割前の基準を取る**

```bash
cmake --preset gcc && cmake --build --preset gcc
ctest --preset gcc -R "test_full_update|AntiferroHeisenberg|Honeycomb|J1J2_AFH|RSVD|Kitaev" --output-on-failure
```

すべて緑であることを確認し、出力を報告ファイルに残す。ここが赤なら BLOCKED として報告する。

- [ ] **Step 2: `prepare_environment` を切り出す**

`full_update.cpp` の現在 122 行(`// Hermite`)から 200 行(`Environment = tensordot(Z, conj(Z), Axes(2), Axes(2));` を含む `else` ブロックの閉じ括弧)までを、そのままの順序で新しい関数本体に移す。

- 入力の `Environment` は `Create_Environment_two_sites(...)` の戻り値、`Theta` は
  `tensordot(tensordot(R1, R2, Axes(1), Axes(1)), op12, Axes(1,3), Axes(0,1))` の戻り値。
- 現在 `envR1`, `envR2` をローカル変数から使っている箇所は、`Environment_in.shape()[0]`,
  `Environment_in.shape()[1]` から取り直す(値は同じ)。
- `Full_Gauge_Fix` が false のときは `Theta_out = Theta_in` とし、`LR1_inv` / `LR2_inv` は
  空のテンソルのままにする(呼び出し側が `Full_Gauge_Fix` を見て使わない)。
- 呼び出し側は `prepare_environment(Environment, Theta, peps_parameters, Environment, Theta, LR1_inv, LR2_inv);`
  のような自己代入を**しない**。別の変数に受けてから代入する。

- [ ] **Step 3: `als_iterate` を切り出す**

現在 236 行(`int count = 0;`)から 358 行(未収束警告と debug 出力の直後、
`if (peps_parameters.Full_Gauge_Fix) { // remove gauge` の直前)までを移す。
`C_phi` / `Old_delta` の初期化、`while` ループ、収束判定、未収束警告、debug 出力を含む。
`D_connect` は `R1.shape()[1]` から取る。

- [ ] **Step 4: 呼び出し側を組み直す**

`Full_update_bond_horizontal` は次の順に読めるようにする(演算は一切変えない):
QR → `Theta` → `Create_Environment_two_sites` → `prepare_environment` →
初期推定(現在 207-235 行の plain `svd` + `slice` + λ)→ `als_iterate` →
ゲージ除去(現在 360-364 行)→ バランシング・組み立て(現在 366-390 行)。

- [ ] **Step 5: 明示的インスタンス化を足す**

`full_update.cpp` 末尾の `template void Full_update_bond_horizontal(...)` に倣い、
`prepare_environment`、`als_iterate`、`Create_Environment_two_sites` を
`real_tensor` と `complex_tensor` の両方で明示的にインスタンス化する。

- [ ] **Step 6: 挙動不変を確認する**

```bash
cmake --build --preset gcc
ctest --preset gcc -R "test_full_update|AntiferroHeisenberg|Honeycomb|J1J2_AFH|RSVD|Kitaev" --output-on-failure
```

Step 1 と同じ結果になること。数値が 1 桁でも動いたら分割が演算順序を変えている。

- [ ] **Step 7: 報告**

`work/fermion/full-update-design/report-task1.md` に、変更ファイル、Step 1 と Step 6 の
ctest 出力、演算順序を変えていないことの根拠(移動した行範囲)を書く。

---

## Task 2: `update_CTM()` をフェルミオン対応にする

**Files:**
- Modify: `src/iTPS/iTPS.cpp:446-451`(`update_CTM`)
- Modify: `src/iTPS/measure.cpp:76-90`
- Test: 既存の `FreeFermion` / `FreeFermionMF` / `FreeFermionSaveLoad` / `FreeFermionSimple` / `test_fermion_layer`(変更禁止)

**Interfaces:**
- Consumes: なし(既存 API のみ)
- Produces: `iTPS<tensor>::update_CTM()` が fermion 模式でも正しい環境を作る。タスク 5 が使う。

- [ ] **Step 1: 基準を取る**

```bash
ctest --preset gcc -R "FreeFermion|test_fermion_layer" --output-on-failure
```

緑であることを確認して報告ファイルに残す。

- [ ] **Step 2: `update_CTM()` に分岐を入れる**

`src/iTPS/iTPS.cpp` の `update_CTM()` を次の形にする。`measure.cpp:76-87` にある
コメント(bare Tn を使う理由)もそのまま移すこと。

```cpp
template <class ptensor>
void iTPS<ptensor>::update_CTM() {
  Timer<> timer;
  if (finfo.enabled) {
    // Bare Tn: the kernel writes sqrt-Schmidt weights into both ends of every
    // bond, so the state is the direct contraction of Tn (same convention the
    // bosonic CTM relies on). Dressing with the full lambda here would
    // double-count the environment weights the CTM itself provides (that is
    // the MeanField-path convention, not the CTM one).
    const std::vector<ptensor> reduced_Tn =
        tenes::fermion::build_reduced_density_tensors(Tn, finfo);
    core::Calc_CTM_Environment_density(C1, C2, C3, C4, eTt, eTr, eTb, eTl,
                                       reduced_Tn, peps_parameters, lattice);
  } else {
    core::Calc_CTM_Environment(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
                               peps_parameters, lattice);
  }
  time_environment += timer.elapsed();
}
```

`iTPS.cpp` に `#include "../fermion/reduced_measure.hpp"` が要る場合は足す。
既存の `update_CTM()` が `time_environment` を計測していなければ、計測は追加せず
既存のままにして、代わりに `measure()` 側で計測している分を壊さないようにする
(どちらでもよいが、二重計上だけはしないこと)。

- [ ] **Step 3: `measure()` を単純化する**

`src/iTPS/measure.cpp:76-90` の `if (!peps_parameters.MeanField_Env && finfo.enabled) { ... } else if (!peps_parameters.MeanField_Env) { update_CTM(); }`
を、次の 3 行に置き換える。

```cpp
  if (!peps_parameters.MeanField_Env) {
    update_CTM();
  }
```

不要になった `#include` は消さない(他の関数が使っている可能性がある)。

- [ ] **Step 4: `update_CTM_density()` は触らない**

`iTPS.cpp:454` の `update_CTM_density()` は有限温度 purification 用で、fermion fold とは
無関係。名前が似ているだけなので混ぜないこと。

- [ ] **Step 5: 確認**

```bash
cmake --build --preset gcc
ctest --preset gcc --output-on-failure
```

**全件**走らせる。共有コードを触っているので範囲を絞らない。

- [ ] **Step 6: 報告**

`work/fermion/full-update-design/report-task2.md` に変更内容と ctest 全件の結果を書く。

---

## Task 3: 開放チャネルで二サイト環境 N を作る

**Files:**
- Modify: `src/fermion/reduced.hpp:529-596`(`build_reduced_pair_halves`)
- Modify: `src/fermion/reduced_measure.hpp:317-406`(吸収列)
- Create: `src/fermion/full_update_env.hpp`
- Test: `test/fermion/full_update_env.cpp`(テスト作成者が用意済み、変更禁止)

**Interfaces:**
- Consumes: `build_pair_state`、`apply_pair_op`、`doubled_pipeline_traced`、
  `reduced_pair_halves`、`fermion::{tensordot, transpose, reshape, conj, max_abs}`
- Produces(タスク 4 が使う): 契約書 §2.1 と §2.2 の署名そのまま。

- [ ] **Step 1: 失敗を確認する**

```bash
cmake --build --preset gcc 2>&1 | tail -40
```

`test/fermion/full_update_env.cpp` が未実装の関数を呼ぶのでリンクまたはコンパイルで落ちる。
それが期待される初期状態。

- [ ] **Step 2: `build_reduced_pair_halves_from_factors` を切り出す**

`reduced.hpp:529-596` の現在の関数から、**gate SVD より後ろ**(`u.multiply_vector(s, 2);` の次の
`const ftensor<tensor> TA6 = ...` から `return {detail::doubled_pipeline_traced(...), ...};` まで)を
新しい関数に移す。既存関数は SVD をして `u.multiply_vector(s, 2)` した u と vt を渡すだけにする。

**一般化が要る点**: 現在のコードは `const std::size_t nk = s.size();` 一つで A 側と B 側の両方を
reshape している。開放経路では k_A と k_B の次元が独立(それぞれ nA², nB²)で、
full update の QR 内部次元は一般に `nA != nB` になる。したがって新しい関数では

```cpp
  const std::size_t nkA = u.parity[2].size();
  const std::size_t nkB = vt.parity[0].size();
```

を取り、TA5 の bundled 軸は `nkA`、TB5 の bundled 軸は `nkB` で作る。
crossing mask は **A 側の bond × k_A にのみ**掛ける(現在の実装も `TA5.multiply_vector(...)` のみで
B 側には掛けていない)。`k_parity` は `u.parity[2]` から取る。
既存の閉じた経路は `nkA == nkB` になるので挙動は変わらない。

- [ ] **Step 3: `absorb_reduced_pair_halves` を切り出す**

`reduced_measure.hpp` の `contract_reduced_pair_halves_horizontal_density_CTM`(317 行)と
`..._vertical_density_CTM`(359 行)から、`joined` を作る直前までを共通の関数に移す。
水平は `left` / `right`、垂直は `top` / `bot` を、それぞれ出力引数 `left` / `right` に返す。
`direction` は `halves.direction` から取る。既存の 2 関数は

```cpp
  tensor left, right;
  absorb_reduced_pair_halves(C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6,
                             halves, left, right);
  tensor joined = mptensor::tensordot(left, right, Axes(1), Axes(1));
  joined = mptensor::transpose(joined, Axes(0, 2, 1, 3));
  return detail::trace_boundary_pairs(joined);
```

の形にする。`ScopedTimer` の名前("measure/twosite/absorb")は既存のまま保つ。

- [ ] **Step 4: `src/fermion/full_update_env.hpp` を書く**

GPL v3 ヘッダとインクルードガード `TENES_SRC_FERMION_FULL_UPDATE_ENV_HPP_` を付ける。
契約書 §2.1 の `full_update_environment` 構造体と `build_full_update_environment` を実装する。
中身は次の順:

1. **Q′ 詰め替え。** `QA`(rank 4)を rank 5 に詰め替える。素の `mptensor::reshape` と
   台帳の手書きのみで、graded 演算は使わない(脚の順序を変えないので Koszul 符号は生じない)。

```cpp
  // horizontal: QA(l,t,b,a) -> QA'(l,t,*,b,a) with a dim-1 even dummy bond
  mptensor::Shape sh = QA.t.shape();
  ftensor<tensor> QA_;
  QA_.t = mptensor::reshape(QA.t, mptensor::Shape(sh[0], sh[1], 1, sh[2], sh[3]));
  QA_.parity = {QA.parity[0], QA.parity[1], parity_vector{false},
                QA.parity[2], QA.parity[3]};
```

   direction 別のダミー位置:

   | direction | QA′ | QB′ |
   |---|---|---|
   | horizontal | (l, t, •, b, a) | (•, t, r, b, β) |
   | vertical | (l, t, r, •, a) | (l, •, r, b, β) |

2. **恒等因子。** `nA = QA.parity[3].size()`、`nB = QB.parity[3].size()` として、
   `I4A(in, out, in′, out′) = δ_{in,in′} δ_{out,out′}`(台帳 `{p_a, p_a, p_a, p_a}`)を作り、
   `u = fermion::reshape(I4A, Shape(nA, nA, nA*nA))` で (in_A, out_A, k_A) にする。
   `I4B(in′, out′, in, out) = δ_{in,in′} δ_{out,out′}`(台帳 `{p_β, p_β, p_β, p_β}`)から
   `vt = fermion::reshape(I4B, Shape(nB*nB, nB, nB))` で (k_B, in_B, out_B) にする。
   `multiply_vector` による特異値の掛け込みは**しない**。
   `fermion::reshape` は隣接脚の融合のみを許すので、この形でなければ落ちる。

3. **fold と吸収。** `build_reduced_pair_halves_from_factors(QA_, QB_, u, vt, direction)` →
   `absorb_reduced_pair_halves(env..., halves, left, right)`。

4. **開放 join。** `tensor M = mptensor::tensordot(left, right, Axes(0, 2), Axes(0, 2));`
   → 脚 (k_A, k_B)。

5. **forbidden block の検査と射影。** `M` を `Shape(nA, nA, nB, nB)` に reshape し、
   台帳 `{p_a, p_a, p_β, p_β}` を付けて `Ntilde` とする。
   `parity_violation(Ntilde) / max(1, max_abs(Ntilde))` を `forbidden_ratio` として記録し、
   `1e-8` を超えたら `std::runtime_error` を投げる(比率と両方の絶対値をメッセージに入れる)。
   閾値以下なら奇成分をゼロにする。**必ず検査してから射影する。**

6. **N と N_plain。** `N = fermion::transpose(Ntilde, Axes(0, 2, 1, 3))` で
   (in_A, in_B, out_A, out_B) にする。`N_plain = N.t` に
   `(-1)^{p(in_A) p(in_B)}` の要素マスク(軸 0 と軸 1)を掛ける。
   マスクは `multiply_vector` では表せない(2 軸の積に依存する)ので、
   `global_index_fast` で局所要素を走査して掛ける。既存の crossing mask の作り方
   (`reduced.hpp:559-568`)に倣い、片方の軸に対する係数ベクトルを作って
   `multiply_vector` を軸ごとに 2 回、という形にはできない点に注意。

- [ ] **Step 5: ビルドとテスト**

```bash
cmake --build --preset gcc
ctest --preset gcc --output-on-failure
```

`test_fermion_layer` の新しい `TEST_CASE`(契約書 T2)が緑になり、
**既存のテストがすべて緑のまま**であること(Step 2, 3 の切り出しが挙動不変であること)。
全件走らせる。

- [ ] **Step 6: 報告**

`work/fermion/full-update-design/report-task3.md` に、切り出しで動かした行範囲、
`nkA` / `nkB` の一般化、ctest 全件の結果を書く。

---

## Task 4: `Full_update_bond_fermion` を書く

**Files:**
- Create: `src/iTPS/core/full_update_fermion.cpp`
- Create: `src/iTPS/core/full_update_fermion.hpp`
- Modify: `src/CMakeLists.txt`(`iTPS_impl` のソース一覧に追加)
- Test: `test/fermion/full_update_bond.cpp`(テスト作成者が用意済み、変更禁止)

**Interfaces:**
- Consumes: タスク 1 の `prepare_environment` / `als_iterate`、
  タスク 3 の `build_full_update_environment`
- Produces(タスク 5 が使う): 契約書 §2.3 の署名そのまま。

- [ ] **Step 1: 失敗を確認する**

```bash
cmake --build --preset gcc 2>&1 | tail -40
```

- [ ] **Step 2: QR と環境**

direction 別の QR 軸:

| direction | site A | site B |
|---|---|---|
| horizontal | `fermion::qr(Tn1, Axes(0,1,3), Axes(2,4), QA, RA)` → QA(l,t,b,a), RA(a,r,s) | `fermion::qr(Tn2, Axes(1,2,3), Axes(0,4), QB, RB)` → QB(t,r,b,β), RB(β,l,s) |
| vertical | `fermion::qr(Tn1, Axes(0,1,2), Axes(3,4), QA, RA)` → QA(l,t,r,a), RA(a,b,s) | `fermion::qr(Tn2, Axes(0,2,3), Axes(1,4), QB, RB)` → QB(l,r,b,β), RB(β,t,s) |

`auto env = fermion::build_full_update_environment(C1..C4, eT1..eT6, QA, QB, direction);`

- [ ] **Step 3: Θ を作る**

```cpp
  // X(a, beta, s1, s2)
  const auto X = fermion::transpose(
      fermion::tensordot(RA, RB, Axes(1), Axes(1)), Axes(0, 2, 1, 3));
  const auto Theta = fermion::tensordot(X, wrapped_gate, Axes(2, 3), Axes(0, 1));
```

`wrapped_gate` は既に wrap 済み・swap 済みで渡ってくる。ここで `wrap_twosite_gate` を
呼んではならない(二重 wrap になる)。

`Theta_tilde = Theta.t` に `(-1)^{p(a) p(β)}`(軸 0 と軸 1)の要素マスクを掛ける。
これは Step 4 のマスクと同じ形なので、共通のヘルパ関数
`apply_pair_mask(tensor&, const parity_vector&, const parity_vector&, int ax1, int ax2)`
をこのファイル内の無名 namespace に置いて使い回す。

- [ ] **Step 4: 環境の準備とパリティ射影**

```cpp
  tensor Env_out, Theta_out, LR1_inv, LR2_inv;
  prepare_environment(env.N_plain, Theta_tilde, peps_parameters, Env_out,
                      Theta_out, LR1_inv, LR2_inv);
```

そのあと **パリティ射影と検査**を行う。`Env_out`(台帳 (p_a, p_β, p_a, p_β))、
`Theta_out` を `mask` で戻したもの(台帳 (p_a, p_β, p_s1, p_s2))、
`Full_Gauge_Fix` が true なら `LR1_inv`(台帳 (p_a, p_a))と `LR2_inv`(台帳 (p_β, p_β))
のそれぞれについて、`parity_violation / max(1, max_abs)` を計算する。

- 1e-8 を超えたら `std::runtime_error`(どのテンソルか、比率をメッセージに入れる)。
- 以下なら奇成分をゼロにする。

これはランク落ちや偶奇縮退で LAPACK がブロック外成分を返す場合への防御であり、
「たぶん大丈夫」で省いてはならない。

- [ ] **Step 5: 初期推定**

ゲージ後の `Theta_out` を `(-1)^{p(a)p(β)}` で戻して `ftensor` にし(台帳 (p_a, p_β, p_s1, p_s2))、
graded `svd_trunc` で分解する。`D_connect` は `Tn1.parity[接続脚].size()`
(水平なら脚 2、垂直なら脚 3)。

```cpp
  fermion::ftensor<tensor> U, VT;
  std::vector<double> s;
  fermion::svd_trunc(Theta_graded, Axes(0, 2), Axes(1, 3), U, s, VT, D_connect);
  // lambda = sqrt(s / ||s||), applied to both sides (same as the bosonic core)
  // svd_trunc returns min(D_connect, rank) values, which can be fewer than
  // D_connect when a parity sector is rank deficient: size the loop off s.
  const std::size_t nkeep = s.size();
  double norm = 0.0;
  for (std::size_t i = 0; i < nkeep; ++i) norm += s[i] * s[i];
  norm = std::sqrt(norm);
  std::vector<double> lambda_c(nkeep);
  for (std::size_t i = 0; i < nkeep; ++i) lambda_c[i] = std::sqrt(s[i] / norm);
  U.multiply_vector(lambda_c, 2);
  VT.multiply_vector(lambda_c, 0);
  auto R1 = fermion::transpose(U, Axes(0, 2, 1));   // (a, D_connect, s1)
  auto R2 = fermion::transpose(VT, Axes(1, 0, 2));  // (beta, D_connect, s2)
```

`svd_trunc` は残す特異値を大きさ順に選んだあと even-first に並べ直すので、
新しいボンド台帳は連続ブロックになる。`plain` の `svd` を使ってはならない
(偶奇の特異値が縮退すると LAPACK が混合ベクトルを返し、台帳が壊れる)。

**セクタ次元の診断**: `U.parity[2]` の偶数個と奇数個を数え、
`simple_update.cpp:33-60` の `TENES_FERMION_SECTOR_LOG` と同じ様式で、
環境変数が設定されているときだけ rank 0 の stderr に出す。
どちらかが 0 になったら `print_level >= PrintLevel::warn` で警告する。

- [ ] **Step 6: ALS**

`R1` を `(-1)^{p(m) p(s1)}`(軸 1 と軸 2)のマスクで `R1_tilde` にし、`R2` はそのまま。
plain 配列にして次を呼ぶ(引数順はタスク 1 の宣言どおり: 環境、Θ、パラメータ、R1、R2)。

```cpp
  tensor R1_plain = R1.t;   // then apply (-1)^{p(m) p(s1)} on axes 1,2
  tensor R2_plain = R2.t;   // unchanged
  als_iterate(Env_out, Theta_out, peps_parameters, R1_plain, R2_plain);
```

- [ ] **Step 7: 後処理**

1. `Full_Gauge_Fix` ならゲージ除去: `R1 = tensordot(LR1_inv, R1, Axes(0), Axes(0));`、
   `R2 = tensordot(LR2_inv, R2, Axes(0), Axes(0));`(plain、bosonic と同じ)。
2. `R1` に `(-1)^{p(m) p(s1)}` を掛け直して graded に戻す。台帳は
   `R1: {p_a, p_m, p_s1}`、`R2: {p_β, p_m, p_s2}`。ここで `p_m` は Step 5 の `U.parity[2]`。
3. `parity_violation(R) / max(1, max_abs(R)) > 1e-8` なら `std::runtime_error`、
   以下なら奇成分をゼロにする。
4. graded バランシング(bosonic 366-386 行と同じ演算を `fermion::` で):

```cpp
  fermion::ftensor<tensor> q1, r1, q2, r2;
  fermion::qr(R1, Axes(0, 2), Axes(1), q1, r1);
  fermion::qr(R2, Axes(0, 2), Axes(1), q2, r2);
  fermion::ftensor<tensor> U2, VT2;
  std::vector<double> s2;
  // fops.hpp has no matrix-shaped svd overload: pass the axes explicitly.
  fermion::svd(fermion::tensordot(r1, r2, Axes(1), Axes(1)), Axes(0), Axes(1),
               U2, s2, VT2);
```

   `s2` はセクタを連結した未ソートの並びだが、`Σ s²` で正規化してから
   `sqrt(s / norm)` を掛けるだけなので順序に依存しない(セクタごとには連続なので
   新しいボンド台帳は連続ブロックのまま)。
   **`s2.size()` は `D_connect` より小さくなり得る**(どちらかのセクタがランク落ちした場合)ので、
   ループは `D_connect` ではなく `s2.size()` で回す。
   `U2.multiply_vector(...)` は軸 1、`VT2.multiply_vector(...)` は軸 0。
   `R1 = fermion::tensordot(q1, U2, Axes(2), Axes(0));` → (a, s1, m′)、
   `R2 = fermion::tensordot(q2, VT2, Axes(2), Axes(1));` → (β, s2, m′)。

5. 組み立て(すべて graded):

| direction | Tn1_new | Tn2_new |
|---|---|---|
| horizontal | `transpose(tensordot(QA, R1, Axes(3), Axes(0)), Axes(0,1,4,2,3))` | `transpose(tensordot(QB, R2, Axes(3), Axes(0)), Axes(4,0,1,2,3))` |
| vertical | `transpose(tensordot(QA, R1, Axes(3), Axes(0)), Axes(0,1,2,4,3))` | `transpose(tensordot(QB, R2, Axes(3), Axes(0)), Axes(0,4,1,2,3))` |

- [ ] **Step 8: 明示的インスタンス化と CMake**

`real_tensor` / `complex_tensor` の両方でインスタンス化する。
`src/CMakeLists.txt` の `iTPS_impl` のソース一覧に `iTPS/core/full_update_fermion.cpp` を足す
(既存の `iTPS/core/full_update.cpp` の隣)。

- [ ] **Step 9: ビルドとテスト**

```bash
cmake --build --preset gcc
ctest --preset gcc --output-on-failure
```

契約書 T3 のテストが緑になり、既存テストが全部緑のままであること。

- [ ] **Step 10: 報告**

`work/fermion/full-update-design/report-task4.md` に、実装の各段の対応、
パリティ射影が実際に発火したかどうか(比率の実測値)、ctest 全件の結果を書く。

---

## Task 5: driver・ガード・フォールバック

**Files:**
- Modify: `src/iTPS/full_update.cpp:28-145`(`iTPS<tensor>::full_update(up)`)
- Modify: `src/iTPS/load_toml.cpp:640-642`(ガード撤去)
- Modify: `src/iTPS/main.cpp:246-249`(フォールバック)
- Test: `test/fermion/boson_equivalence_full.py.in`、`test/fermion/free_fermion_full.py.in`、
  `test/input.cpp` への追記(テスト作成者が用意済み、変更禁止)
- Modify: `test/CMakeLists.txt`(新しい python テストの登録。**テスト本体ではないので変更してよい**)

**Interfaces:**
- Consumes: タスク 2 の `update_CTM()`、タスク 4 の `Full_update_bond_fermion`

- [ ] **Step 1: 失敗を確認する**

`ctest --preset gcc -R "input|BosonEquivalenceFull|FreeFermionFull" --output-on-failure`

- [ ] **Step 2: driver のフェルミオン分岐**

`src/iTPS/full_update.cpp` の `full_update(up)` に、`simple_update.cpp:65-130` と同じ構造の
分岐を入れる。

onesite:

```cpp
  if (up.is_onesite()) {
    const int source = up.source_site;
    if (finfo.enabled) {
      auto fTn = tenes::fermion::wrap_Tn(Tn[source], finfo, source);
      tenes::fermion::ftensor<tensor> fop{
          up.op, {finfo.phys[source], finfo.phys[source]}};
      auto updated = tenes::fermion::tensordot(fTn, fop, mptensor::Axes(4),
                                               mptensor::Axes(0));
      tenes::fermion::unwrap_Tn(updated, Tn[source], finfo, source);
      return;
    }
    Tn[source] =
        tensordot(Tn[source], up.op, mptensor::Axes(4), mptensor::Axes(0));
  } else {
```

これは `simple_update` と同一なので、共通ヘルパ(たとえば
`iTPS<tensor>::apply_onesite_gate_fermion(const EvolutionOperator<tensor>&)`)に括り出して
両方から呼ぶ。

twosite:

```cpp
    if (finfo.enabled) {
      int s1 = source, s2 = target;
      int s1_leg = source_leg, s2_leg = target_leg;
      auto fop = tenes::fermion::wrap_twosite_gate(up.op, finfo.phys[source],
                                                   finfo.phys[target]);
      if (source_leg == 0 || source_leg == 1) {
        std::swap(s1, s2);
        std::swap(s1_leg, s2_leg);
        fop = tenes::fermion::transpose(fop, mptensor::Axes(1, 0, 3, 2));
      }
      auto fTn1 = tenes::fermion::wrap_Tn(Tn[s1], finfo, s1);
      auto fTn2 = tenes::fermion::wrap_Tn(Tn[s2], finfo, s2);
      tenes::fermion::ftensor<tensor> fTn1_work, fTn2_work;
      if (s1_leg == 2) {
        core::Full_update_bond_fermion(
            C1[s1], C2[s2], C3[s2], C4[s1],
            eTt[s1], eTt[s2], eTr[s2], eTb[s2], eTb[s1], eTl[s1],
            fTn1, fTn2, fop,
            tenes::fermion::reduced_pair_direction::horizontal,
            peps_parameters, fTn1_work, fTn2_work);
      } else {
        core::Full_update_bond_fermion(
            C1[s1], C2[s1], C3[s2], C4[s2],
            eTt[s1], eTr[s1], eTr[s2], eTb[s2], eTl[s2], eTl[s1],
            fTn1, fTn2, fop,
            tenes::fermion::reduced_pair_direction::vertical,
            peps_parameters, fTn1_work, fTn2_work);
      }
      finfo.virt[s1][s1_leg] = fTn1_work.parity[s1_leg];
      finfo.virt[s2][s2_leg] = fTn2_work.parity[s2_leg];
      tenes::fermion::unwrap_Tn(fTn1_work, Tn[s1], finfo, s1);
      tenes::fermion::unwrap_Tn(fTn2_work, Tn[s2], finfo, s2);
      tenes::fermion::validate_neighbor_consistency(finfo, lattice);
      update_CTM();
      return;
    }
```

窓環境の選び方(**回転しない**):

| direction | C1, C2, C3, C4 | eT1..eT6 |
|---|---|---|
| horizontal | C1[s1], C2[s2], C3[s2], C4[s1] | eTt[s1], eTt[s2], eTr[s2], eTb[s2], eTb[s1], eTl[s1] |
| vertical | C1[s1], C2[s1], C3[s2], C4[s2] | eTt[s1], eTr[s1], eTr[s2], eTb[s2], eTl[s2], eTl[s1] |

`lambda_tensor` は触らない(bosonic full update も触らない)。
fast full update の `Right_move` 系の分岐には**入らない**(タスク 5 Step 4 で `false` に
落としてあるが、念のためフェルミオン経路は `update_CTM()` に直行する)。

- [ ] **Step 3: ガードを撤去する**

`src/iTPS/load_toml.cpp:640-642`

```cpp
  if (has_positive_steps(peps_parameters.num_full_step)) {
    throw_fermion_guard("full update");
  }
```

を削除する。周囲の他のガード(`meanfield_env` など)は残す。

- [ ] **Step 4: fast full update のフォールバック**

`src/iTPS/main.cpp` の `gen_param()` 直後・`peps_parameters.Bcast(comm)`(249 行)**直前**に
入れる。`load_toml.cpp:557` の `has_positive_steps` は無名 namespace にあり main.cpp からは
見えないので、条件は直接書く。

```cpp
  if (peps_parameters.fermion && peps_parameters.Full_Use_FastFullUpdate &&
      std::any_of(peps_parameters.num_full_step.begin(),
                  peps_parameters.num_full_step.end(),
                  [](int n) { return n > 0; })) {
    if (mpirank == 0) {
      std::cerr << "WARNING: fermion mode disables Full_Use_FastFullUpdate "
                   "because the fast update reuses bare-Tn CTM moves that are "
                   "not fermion-aware in this version"
                << std::endl;
    }
    peps_parameters.Full_Use_FastFullUpdate = false;
  }
```

`<algorithm>` の include が無ければ足す。**`Bcast` の前**に置くこと。後ろに置くと
全 rank に伝わらない。

- [ ] **Step 5: python テストを登録する**

`test/CMakeLists.txt` に、既存の `FreeFermionMF` の記述に倣って追加する。

```cmake
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fermion/boson_equivalence_full.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fermion/boson_equivalence_full.py @ONLY)
add_test(NAME BosonEquivalenceFull
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/boson_equivalence_full.py)
set_tests_properties(BosonEquivalenceFull PROPERTIES TIMEOUT 1800)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fermion/free_fermion_full.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_full.py @ONLY)
add_test(NAME FreeFermionFull
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_full.py)
set_tests_properties(FreeFermionFull PROPERTIES TIMEOUT 1800)
```

テスト作成者が置いたファイル名と一致させること。違っていたら BLOCKED で報告する。

- [ ] **Step 6: ビルドとテスト**

```bash
cmake --build --preset gcc
ctest --preset gcc --output-on-failure
```

全件緑であること。

- [ ] **Step 7: 報告**

`work/fermion/full-update-design/report-task5.md` に、driver の分岐、ガード撤去、
フォールバックの実装地点、ctest 全件の結果を書く。

---

## 実行順序と依存

```
Task 1 ─┐
        ├─→ Task 4 ─→ Task 5
Task 3 ─┘              ↑
Task 2 ────────────────┘
```

Task 1、2、3 は並列に走らせてよい。Task 4 は 1 と 3 の完了後、Task 5 は 2 と 4 の完了後。

## 完了後

1. Claude が全件 ctest を独立に回して検証する(Codex の報告を鵜呑みにしない)。
2. タスクごとに新規サブエージェントでレビュー。変異テスト(契約書 T3-vii、T2-iv の
   N マスク変異)を必ず指示する。
3. 最上位モデルで全ブランチレビューを 1 回。
4. `work/fermion/full-update-mott/` で Hubbard U=4 の Mott 崩壊検証(ctest 外、複数シード)。
