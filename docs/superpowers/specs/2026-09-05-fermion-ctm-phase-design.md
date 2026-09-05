# フェルミオン CTM 環境の位相 — 設計書(改訂 2.1)

日付: 2026-09-05(改訂 2.1。改訂 1 の「C1 を一括で回す」案はテスト作成者の計測で棄却。
改訂 2.1 は Codex の再レビュー `work/fermion/ctm-phase/codex-review-design2.md` の 7 件を反映)
ブランチ: `fermion`
記録: `work/fermion/full-update-design/FINDINGS-task5-energy.md` の (c) 節、
`work/fermion/ctm-phase/codex-review-design.md`(改訂 1 のレビュー)、
`work/fermion/ctm-phase/report-testauthor.md` §7.1(改訂 1 を棄却した計測)
試作: `work/fermion/full-update-design/diagnostics-session2-kept.patch`(環境変数分岐、本設計で撤去)

## 1. 問題

フェルミオン mode の CTM(`iTPS::update_CTM()` → `build_reduced_density_tensors` →
`core::Calc_CTM_Environment_density`)は、fold 済みの単層テンソルを**一様ベクトル**で潰して
C, eT を初期化する(`ctm_single.cpp:789-`)。fold の奇×奇成分には Koszul 符号が入るので、
初期 C, eT の位相(実数なら ±1、複素なら任意)は状態依存である。

計測結果を説明する**局所模型**(projector の位相応答を無視し、各テンソルに 1 つのスカラー位相を
置いたヒューリスティック。実装から一般に導ける不変則ではない): 固定点近傍で C ≈ e^{iθ_c}C*、
eT ≈ e^{iθ_e}eT* とおくと、`Calc_Next_CTM`(C_new ∝ P·C·eT·P)は θ_c → θ_c + θ_e、
`Calc_Next_eT`(eT_new ∝ Tn·eT·P·P)は θ_e を保存する(規格化は `max_abs` で位相を除かない)。
bosonic の iTPS 経路は初期化 Tn⊗conj(Tn) が正定値で θ_e = 0 だが、fold の一様ベクトル初期化では
θ_e ≠ 0 になり得る。実際には方向・サイトごとに異なる eT が窓に入るので位相は 1 つの θ_e では
表せず、それが「窓ごとに位相が違う」計測と整合する。**設計の正当性はこの位相則には置かず、
§2 の「各消費窓が自分の ⟨1⟩ で正規化する」という局所契約に置く。**
計測(テスト作成者、`work/fermion/ctm-phase/testauthor/calib/`、複素 2×2 / 3×2、D=2、chi=8):

1. one-site 窓の位相は全サイトで一致するが、**反復ごとに回る**(0.64 rad/反復)。大きさは 10 桁で固定点。
2. **窓ごとに位相が違う**: two-site 窓は含む eT の数と種類が one-site 窓と異なるので、
   同じ反復でも位相が 0.2〜0.5 rad ずれる窓がある(例: h(0,1) と h(1,0) で −2.565 と −2.896)。
3. 実数では位相は ±1 に限られるので、2×2 / 3×2 の実測(D=2, 3)では全窓が一致していた。
   これは偶然であり、設計の前提にできない。

測定値(⟨O⟩/⟨1⟩)は位相が打ち消されるので**影響なし**。壊れるのは次の 4 箇所:

1. `core::Check_Convergence_CTM_RDM`(`ctm.cpp:1025-`)は one-site RDM の trace に
   「実部 > 0 かつ |虚部| ≤ 1e-6·|trace|」を要求し、満たさないと `rdm_dist = NaN` で
   「未収束」を返し続ける。corner の特異値が 1e-15 で収束していても `iteration_max` まで回る。
2. フェルミオン full update の `prepare_environment`(bosonic 共有)は N_plain をエルミート化
   ((N + N†)/2)してから正値化する。N が e^{iφ}·N_true なら エルミート化で cos φ · N_true になり
   (φ = π で負定値、φ = π/2 でゼロ)、固有値が落ちて Env ≈ 0 → 「zero singular-value norm」。
3. `measure_onesite_rdm` が出力する `onesite_density_matrix.dat`(fermion CTM 経路)が位相付きの
   複素行列になる(trace が実正でない)。実数では符号が負になり得る。
4. one-site / two-site 測定(`onesite_obs.cpp:111-162`、`twosite_obs.cpp:173-231, 373-397`)は各窓の
   生の norm について `min Re ≤ 0` → 「Norm is negative」、`max |Im| > 1e-6` → 「Norm is not real」の
   **偽警告**を出す(複素 fermion の実行で実測。time evolution は `measure()` を反復するので増幅する)。
   出力される正規化観測量の値は正しい。

## 2. 方針

**環境の位相は窓ごとのゲージ自由度**である。環境そのものを canonicalize しようとせず、
判定は位相不変にし、環境を消費する側が**自分の窓の ⟨1⟩ で**位相を決める。

改訂 1(全 `C1[i]` を one-site 窓の位相で一括回転)を採らない理由: 窓ごとに位相が違うので
C1 だけでは全窓を揃えられない(§1 の 2)。eT の位相まで規格化して θ_e = 0 に固定する案は、
「eT の位相」を定義する規約が必要で、それが物理的な位相(⟨1⟩ > 0)と一致する保証がない。
境界 trace で初期化する案の有効性は否定しないが(棄却した計測は C1 一括回転案であって
初期化案ではない)、初期化だけで全反復・全窓の実正性を保証する根拠がないので、それを正しさの
前提にはせず、初期化方式に独立な消費側契約を採る。消費側で直す案は、どの窓でも ⟨1⟩ という
物理量で位相が一意に決まり、CTM の位相構造を仮定しない。

## 3. 変更

### 3.1 `Check_Convergence_CTM_RDM` を位相不変にする(`src/iTPS/core/ctm.cpp`)

現在: 各サイトの one-site RDM ρ_i(plain、`Contract_one_site_RDM_{iTPS,density}_CTM`)について
`trace = Σ_k ρ_i(k,k)`、`valid = |trace| ≥ 1e-12·max|ρ_i| ∧ Re trace > 0 ∧ |Im trace| ≤ 1e-6·|trace|`、
`rdm_dist = max_i max_{kl} |ρ_i^new(k,l) − ρ_i^old(k,l)| / |trace_i^new|`。

変更後: 引数 `const bool phase_invariant = false` を追加する(`Calc_CTM_Environment_density` にも
同名の引数を足し、既定値 `false`。`iTPS::update_CTM()` の fermion 分岐だけが `true` を渡す)。

- 健全性検査(全経路共通): `trace_abs` と `max_abs` と ρ_i の全要素が有限(`std::isfinite`)で、
  `trace_abs > 0` かつ `trace_abs ≥ 1e-12·max_abs`。破れたら `valid = false`。
- `phase_invariant == false`(bosonic、finite-T): 現行どおり `Re trace > 0` と
  `|Im trace| ≤ 1e-6·|trace|` を要求し、ρ_i はそのまま比較する。**コードパスと演算順序を変えない**
  (ビット同一は §5 の A/B で確認)。
- `phase_invariant == true`(fermion): `phase_i = trace_i/|trace_i|` で ρ_i を割ってから比較する
  (無条件に割る)。実正値条件は課さない。位相が反復ごとに回っても(§1 の 1)、正規化後の比較は
  大きさの固定点で収束する。
- `rdm_old` には比較に使った ρ_i(fermion なら正規化済み)を保存する。
- `rdm_dist` の分母は |trace_i|(変更なし)。

### 3.2 `build_full_update_environment` で N の位相を正規化する(`src/fermion/full_update_env.hpp`)

現在: 開放 join → forbidden 検査 → `project_even` → `N = transpose(Ñ, (0,2,1,3))` →
`N_plain = N ⊙ mask_{in_A in_B}`。

変更後(`N` を作った直後、`N_plain` を作る前):

1. `norm = Σ_x N.t(x) · I_wrap.t(x)`(plain 要素積和)。`I_wrap = wrap_twosite_gate(I, p_a, p_β)`、
   I は擬似サイトの物理脚 (a, β) に対する恒等 `δ_{in_A,out_A} δ_{in_B,out_B}`。
   これは擬似サイト対のノルム ⟨pair_Q|pair_Q⟩ をその窓の環境で評価したもので、
   位相が 1 の環境なら実正。分散テンソルなら `allreduce_sum`。
2. 検査: `N.t` の全要素と `max_abs(N)` が有限であること(collective に確認)、`norm` が有限で
   `|norm| > 1e-12·max_abs(N)·(nA·nB)`(実質ゼロでない)こと。破れたら
   `std::runtime_error("build_full_update_environment: ...")`(値を含める)。
3. `phase = norm/|norm|`。`|phase − 1| > 1e-14` なら `N.t *= conj(phase)`(実テンソルでは
   `phase.real()` = ±1)。位相が 1 なら触らない。
4. 返り値 `full_update_environment` に `std::complex<double> phase`(取り除いた位相、
   何もしなければ 1)を追加する。`N_plain` は正規化後の N から作る。

これで `N_plain` は環境の位相に依らず、その窓で ⟨1⟩ > 0 となる規約に固定される。
`prepare_environment` のエルミート化・正値化はそのまま正しく働く。
`Full_update_bond_fermion` は変更しない(改訂 1 の計量ガードは不要)。norm が退化していれば
2 で例外になるが、**位相補正そのものの誤り**(共役の取り忘れ、I の脚順、N の一部だけを回す)は
例外にならないので、T8-ii(独立経路との一致)、T2-vi、T3-i で検出する。

### 3.3 `iTPS::update_CTM()` の fermion 分岐(`src/iTPS/iTPS.cpp`)

`Calc_CTM_Environment_density(..., /*initialize=*/true, /*phase_invariant=*/true)` を渡すだけ。
環境の位相固定は行わない。`update_CTM_density()`(finite-T)と bosonic 分岐は変更しない。

### 3.4 one-site RDM 出力の位相正規化(`src/iTPS/onesite_obs.cpp`、`src/iTPS/core/ctm.hpp`)

位相正規化を `core` の小関数として切り出す(テストから呼べるように):

```cpp
namespace tenes::itps::core {
//! Divide a gathered one-site RDM by the phase of its trace so that the
//! trace becomes real positive. |trace| and the Frobenius norm are kept.
//! Returns the phase that was removed (exactly 1 when nothing was changed:
//! |phase - 1| <= 1e-14 leaves the matrix untouched).
//! Throws std::runtime_error if any element is not finite or the trace vanishes.
template <class value_type>
std::complex<double> normalize_rdm_phase(small_tensor<value_type>& rdm);
}
```

`measure_onesite_rdm` の `is_fermion_ctm` 分岐で、`rdm_local` を作った直後にこれを呼ぶ。
bosonic / finite-T / MF の分岐は触らない。trace の大きさは変えない(bosonic の生出力と同じ規約。
trace 正規化は本設計の外)。

### 3.5 測定経路の norm 診断(`src/iTPS/onesite_obs.cpp`、`src/iTPS/twosite_obs.cpp`)

fermion CTM 経路(`finfo.enabled && !MeanField_Env`)では、各窓の norm を診断する前に
`norm_i / (norm_i/|norm_i|)` = |norm_i| に置き換えてから既存の警告条件(`min Re ≤ 0`、
`max |Im| > 1e-6`)を評価する。ただし norm が非有限またはゼロのときは既存の警告を出す。
**観測量の値の計算(⟨O⟩/⟨1⟩ の除算)は変えない**。bosonic / finite-T / MF の診断経路は触らない。
対象: `onesite_obs.cpp` の measure_onesite(111-162 行付近)、`twosite_obs.cpp` の
measure_twosite(373-397 行付近)。multisite / correlation / transfer-matrix は fermion mode で
入力段階で無効なので対象外。

### 3.6 一時診断の撤去

`diagnostics-session2-kept.patch` の環境変数分岐(`TENES_FU_SIGN_FIX`, `TENES_FU_FORBIDDEN_TOL`,
`TENES_CTM_RDM_ABS`, `TENES_CTM_RDM_LOG`, `TENES_CTM_RDM_PHASE`, `TENES_CTM_PHASE_FIX`,
`TENES_CTM_PHASE_LOG`)はすべて撤去する(実装者に渡す前に Claude が `git checkout src/` で戻す)。
`FORBIDDEN_TOL`(forbidden 閾値)は別課題として残るが、環境変数での上書きはしない。

## 4. 影響と非目標

- 測定値は位相に依らないので既存 golden(fermion / bosonic / finite-T)は変わらない。
  bosonic と finite-T は `phase_invariant = false` で現行コードパスを通るので**ビット同一**
  (§5 の A/B で機械確認)。
- fermion 実数 D=2(既存 golden `FreeFermion`、位相 1): 3.1 で ρ_i を位相 1 で割る演算が入るので
  `rdm_dist` の最下位ビットが変わりうる。収束値は同じで golden は許容内。3.2 と 3.4 は位相 1 なら
  触らない。ビット同一は要求しない。
- 位相 ≠ 1 のケース(複素、D≥3)は CTM の反復回数が減り(判定が効く)、FU が動くようになる。
- `onesite_density_matrix.dat` の fermion 出力は trace が実正になる(bosonic と同じ見た目)。
- 複素 fermion の測定で「Norm is negative / not real」の偽警告が出なくなる。
- CTM の環境テンソル(C, eT)そのものの位相は不定のまま。将来、環境を直接消費する新しい経路
  (multisite の fermion 対応など)を足すときは、その窓の ⟨1⟩ で位相を決めること。
- forbidden 閾値(`build_full_update_environment` の 1e-8)の扱いは本設計の外。
- T5 の D=3 化は本設計の契約に含める(§5)。

## 5. 検証(契約書 `2026-09-04-fermion-full-update-contract.md` の §2.5、§3.6、§3.5 T5)

| # | 検査 | 形式 |
|---|---|---|
| T8-i | `Check_Convergence_CTM_RDM(phase_invariant = true)` の位相不変性。独立な履歴を持つ二系列で 2 回ずつ: 基準 E0 → E1、変換 e^{iφ0}E0 → e^{iφ1}E1(全 C1 に掛ける。φ0 ≠ φ1)。初回は未収束、2 回目の `rdm_dist` と bool が一致(相対 1e-12)。零 RDM・NaN・Inf は `valid = false`。`phase_invariant = false` では位相を掛けた系列が未収束のまま。real/complex、`is_density` 両方 | doctest |
| T8-ii | N の位相不変性: `build_full_update_environment` を、(a) 環境そのまま、(b) 全 `C1[i]` に e^{iφ}(実数は −1)を掛けた環境、(c) 1 つの eT(たとえば eTt[s1])に e^{iφ} を掛けた環境、で呼び、`N` と `N_plain` が相対 1e-12 で一致し、返り値 `phase` が対応して変わる。**独立検証**: 正規化前の窓ノルム(`phase` × 正規化後の Σ N·I_wrap)が閉じた測定経路 `contract_reduced_pair_halves_density_CTM(env, build_reduced_pair_halves(QA′, QB′, I_wrap, dir))` と相対 1e-12 で一致(nA ≠ nB を含める)。`N_plain` はエルミート(相対 1e-12)かつ半正定値。状態: 2×2 実数 D=3、2×2 複素 D=2、3×2 複素 D=2(ランダム偶 Tn)、両方向。CTM は `count < iteration_max` を前提アサート。例外: 全 `C1[i]` をゼロ、I_wrap の台外(in_A ≠ out_A の成分)にだけ NaN、Inf を注入した環境、のそれぞれで `std::runtime_error` | doctest |
| T8-iii | 複素 fermion E2E(`is_real = false`、自由フェルミオン D=2、**mu = 0**、SU のみ)で「CTM did not converge」「Norm is negative」「Norm is not real」のいずれも出ない。エネルギーの虚部 ≤ 1e-10。密度が既存 `FreeFermion`(mu = 0)の許容内 | E2E |
| T8-iv | FU の位相不変性: T3-i の入力で、C1 に −1(複素は e^{iφ})を掛けても `Full_update_bond_fermion` の出力 pair state が不変(相対 1e-10、スケール除く)。両方向 | doctest |
| T8-v | `normalize_rdm_phase`: 既知の位相 e^{iφ} を掛けた RDM に適用すると基準 RDM と全要素一致(相対 1e-12)、|trace| と Frobenius ノルムが保存、戻り値が e^{iφ}(相対 1e-12)。位相 1(実正 trace)では要素完全一致で no-op、戻り値 exactly 1。ゼロ行列・NaN・Inf(非対角のみでも)で `std::runtime_error`。real/complex。加えて E2E(T8-iii)の `onesite_density_matrix.dat` の各サイト trace が実正(|Im| ≤ 1e-10·|tr|、Re > 0) | doctest + E2E |
| T2-vi(complex) | 既存 T2-vi の complex_tensor 版。環境は `Calc_CTM_Environment_density(phase_invariant = true)` で作る(位相固定なし)。測定経路のノルムは位相付きなので、比較は **|ノルム|** と **比 ⟨G⟩** で行う | doctest |
| T5(改訂) | 自由フェルミオン **D=3、chi=12**、`iteration_max = 200`、SU 1000 step 後に FU 1〜4 sweep で `E_FU < E_SU − 1e-4` かつ厳密値に接近、「CTM did not converge」が出ない | E2E |
| T7 | 既存テストが無変更で緑 | 既存 |
| A/B(検証者) | 変更前後のバイナリで、bosonic golden(`AntiferroHeisenberg_real/complex/mf`、`RSVD`、`J1J2_AFH`、`Honeycomb`、`Honeycomb_skew`、`TE_TFI`)と finite-T golden(`FT_TFI_square`、`FT_Kitaev`)を同一入力・1 thread・MPI off で走らせ、`output_*/*.dat`(`time.dat` 除く)の sha256 が一致。基準は `work/fermion/ctm-phase/ab/before/sha256.txt`(315a8622 のバイナリ) | Claude が実施 |

## 6. タスク分割

1 タスク(Codex 1 バッチ): 3.1〜3.6。テストは契約から別担当が先に書く(RED を確認してから実装)。
