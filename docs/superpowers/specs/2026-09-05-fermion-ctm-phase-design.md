# フェルミオン CTM 環境の全体位相の固定 — 設計書

日付: 2026-09-05
ブランチ: `fermion`
記録: `work/fermion/full-update-design/FINDINGS-task5-energy.md` の (c) 節(原因の確認実験、試作の結果)
試作: `work/fermion/full-update-design/diagnostics-session2-kept.patch`(環境変数分岐、本設計で正式化して撤去)

## 1. 問題

フェルミオン mode の CTM(`iTPS::update_CTM()` → `build_reduced_density_tensors` →
`core::Calc_CTM_Environment_density`)は、fold 済みの単層テンソルを**一様ベクトル**で潰して
C, eT を初期化する(`ctm_single.cpp:789-`)。fold の奇×奇成分には Koszul 符号が入っているので、
初期 C, eT の位相(実数なら ±1、複素なら任意)は状態依存である。反復
(`Calc_Next_CTM` / `Calc_Next_eT` は `max_abs` 規格化、projector は PU/PL のペア)は
C, eT について多重線形で規格化が絶対値なので、閉じた窓(C1·eTt·C2·eTr·C3·eTb·C4·eTl·Tn)の
位相は反復の不変量になる。**単一の位相が全サイト・全窓に共通することの証明は与えない。**
実証されているのは次の範囲: 2×2 実数 D=2/3、2×2 複素 D=2/3、3×2 実数 D=3、3×2 複素 D=2 で、
one-site 窓の位相が全サイトで 6 桁一致(`work/fermion/ctm-phase/cells/`、FINDINGS (c))。
fermion mode は入力段階で `L_sub < 2` と `skew ≠ 0` を拒否する(`load_toml.cpp`)ので、
非正方セルのうち検証すべきは Lx, Ly ≥ 2・skew = 0 だけであり、3×2 で覆っている。
したがって収束した環境は「正しい環境 × e^{iφ}」になり、φ は制御できない。
bosonic の iTPS 経路は初期化が Tn⊗conj(Tn) の縮約で正定値なので φ = 0 に固定されるが、
fold の density 経路にはその保証がない。

測定値(⟨O⟩/⟨1⟩)は位相が打ち消されるので**影響なし**。壊れるのは次の 2 箇所だけ:

1. `core::Check_Convergence_CTM_RDM`(`ctm.cpp:1025-`)は one-site RDM の trace に
   「実部 > 0 かつ |虚部| ≤ 1e-6·|trace|」を要求し、満たさないと `rdm_dist = NaN` で
   「未収束」を返し続ける。corner の特異値が 1e-15 で収束していても `iteration_max` まで回る。
   実測: D=3 実数(φ = π)、複素 D=2(φ = −1.1 rad)、ランダム複素 Tn(テスト作成者の観察)。
2. フェルミオン full update の `prepare_environment`(bosonic 共有)は N_plain を正値化する
   (`w/w_max > Inverse_Env_cut` の固有値だけ残す)。φ = π では N_plain が負定値になり
   固有値が全部落ちて Env ≈ 0、Θ′ ≈ 0 → `svd_trunc` で奇セクタ消滅 → 「zero singular-value norm」。

## 2. 方針

**環境の全体位相はゲージ自由度**である。判定は位相不変にし、収束後に物理的に固定する。

- CTM の初期化ベクトルを変える案(境界 trace)は採らない。初期化を変えても位相 1 に落ちる
  **保証**は得られないので、収束後の検査と固定はどのみち必要になる。本案は初期化方式に依存せず
  公開される環境を canonicalize する。
- 判定・正値化の側で位相を「許容」するだけの案(試作 `RDM_ABS` / `SIGN_FIX`)も採らない。
  後段で未正規化の環境を正定値として扱うのは FU の `prepare_environment`
  (`N_plain` の固有値を `w/w_max > Inverse_Env_cut` で切る)であり、そこで直すより
  環境の側で**一度固定する**方が単純で、将来の利用者(複素 fermion の測定・FU)も同じ前提を持てる。

## 3. 変更

### 3.1 `Check_Convergence_CTM_RDM` を位相不変にする(`src/iTPS/core/ctm.cpp`)

現在: 各サイトの one-site RDM ρ_i(plain、`Contract_one_site_RDM_{iTPS,density}_CTM`)について
`trace = Σ_k ρ_i(k,k)`、`valid = |trace| ≥ 1e-12·max|ρ_i| ∧ Re trace > 0 ∧ |Im trace| ≤ 1e-6·|trace|`、
`rdm_dist = max_i max_{kl} |ρ_i^new(k,l) − ρ_i^old(k,l)| / |trace_i^new|`。

変更後: 引数 `bool phase_invariant` を追加する(`Calc_CTM_Environment_density` にも同名の
引数を足し、既定値 `false`。`iTPS::update_CTM()` の fermion 分岐だけが `true` を渡す)。

- 健全性検査(全経路共通、Codex 指摘 1): `trace_abs` と `max_abs` と ρ_i の全要素が有限
  (`std::isfinite`)で、`trace_abs > 0` かつ `trace_abs ≥ 1e-12·max_abs`。破れたら `valid = false`。
  現行は実正値条件が零 RDM を弾いていたが、それを外すとここが唯一の検査になる。
- `phase_invariant == false`(bosonic、finite-T): 現行どおり `Re trace > 0` と
  `|Im trace| ≤ 1e-6·|trace|` を要求し、ρ_i はそのまま比較する。**コードパスを変えない**
  (ビット同一は §5 の A/B で確認)。
- `phase_invariant == true`(fermion): `phase_i = trace_i/|trace_i|` で ρ_i を割ってから比較する
  (無条件に割る。閾値分岐は置かない)。実正値条件は課さない。
- `rdm_old` には比較に使った ρ_i(fermion なら正規化済み)を保存する。
- `rdm_dist` の分母は |trace_i|(変更なし)。

位相不変判定を fermion に限定する理由(Codex 指摘 3): finite-T の density 経路でも trace は
理論上実正だが丸めで虚部が残り、閾値付きの分岐では挙動が変わりうる。フラグで分ければ
bosonic / finite-T は現行コードそのものになる。

### 3.2 `iTPS::update_CTM()` の fermion 分岐で位相を固定する(`src/iTPS/iTPS.cpp`)

位相固定は `core` の公開関数にする(テストから直接呼べるように):

```cpp
namespace tenes::itps::core {
//! Fix the overall phase of a folded (fermion) CTM environment so that the
//! one-site norm trace(Contract_one_site_RDM_density_CTM(...)) is real and
//! positive on every site. Multiplies every C1[i] by the conjugate phase and
//! returns the phase that was removed (exactly 1 when nothing was changed).
//! Throws std::runtime_error if any one-site norm vanishes or is not finite,
//! or if the sites do not share one phase (|phase_i - phase_0| > 1e-8).
template <class tensor>
std::complex<double> fix_environment_phase(
    std::vector<tensor>& C1, const std::vector<tensor>& C2,
    const std::vector<tensor>& C3, const std::vector<tensor>& C4,
    const std::vector<tensor>& eTt, const std::vector<tensor>& eTr,
    const std::vector<tensor>& eTb, const std::vector<tensor>& eTl,
    const std::vector<tensor>& reduced_Tn, const SquareLattice& lattice);
}
```

`iTPS::update_CTM()` の fermion 分岐は `Calc_CTM_Environment_density(..., phase_invariant = true)`
の直後にこれを呼ぶ。中身:

1. 各サイト i について `n_i = trace(Contract_one_site_RDM_density_CTM(C1[i],C2[i],C3[i],C4[i],
   eTt[i],eTr[i],eTb[i],eTl[i], reduced_Tn[i]))`。分散テンソルなので、全サイトの対角和
   (実部・虚部)を 1 本の `std::vector<double>` に詰めて `allreduce_sum`(`src/mpi.hpp:190`)を
   **1 回**だけ呼ぶ。communicator は `rdm.get_comm()`。戻り値が非 0 なら `std::runtime_error`。
2. 検査: すべての i で `n_i` が有限かつ `|n_i| > 0`、かつ `|n_i/|n_i| − n_0/|n_0|| ≤ 1e-8`。
   破れたら `std::runtime_error`(全サイトの n_i を含める)。全サイト一致は §1 の範囲で実証済み。
   一致しない環境をどう回復するかは設計しない(そのような入力は現状存在しないので、
   出たら不変量の理解が誤っていたことになり、設計に戻る)。
3. 固定: `phase = n_0/|n_0|`。`|phase − 1| > 1e-14` のときだけ、すべての i で
   `C1[i] *= conj(phase)`(実テンソルでは `phase.real()`、すなわち ±1)。位相が既に 1 なら
   C1 に触らない(D=2 実数の既存 golden はビット同一)。
   C1 だけを回す根拠: fermion mode で有効な閉じた窓(one-site、two-site 水平/垂直、FU の
   開放チャネル窓)はどれも C1[·] をちょうど 1 つ含む(`contract_density_ctm/ctm.cpp`、
   `reduced_measure.hpp`、`full_update_env.hpp`)。correlation・transfer-matrix・multisite は
   fermion mode では入力段階で無効化されているので対象外(Codex 指摘 5)。
4. 戻り値は取り除いた位相(何もしなければ 1)。

`update_CTM_density()`(finite-T)と bosonic 分岐は変更しない。

### 3.3 FU 核のガード(`src/iTPS/core/full_update_fermion.cpp`)

`build_full_update_environment` の直後、`X = transpose(tensordot(RA, RB, (1),(1)), (0,2,1,3))` から
`X̃ = mask_{aβ} X` を作り、`nrm = ⟨X̃|N_plain|X̃⟩`(bosonic 規約
`trace(tensordot(N_plain, X̃, Axes(0,1), Axes(0,1)), conj(X̃), all, all)`)を計算する。
次のどちらかで `std::runtime_error`(値を含める):
`Re nrm ≤ 1e-10 · max_abs(N_plain) · ‖X̃‖²_F`、または `|Im nrm| > 1e-8 · |nrm|`。
スケール付きにする理由(Codex 指摘 6): 正しい半正定値計量でも丸めで微小負になりうるのは
X が計量の零空間にある場合だけであり、X は現在の状態そのものなのでノルムは O(1) である。
ゼロ空間にある X(状態がゼロ)は異常として扱う。これは §3.2 の防衛であり、符号を直す処理は
**入れない**(試作の `SIGN_FIX` は撤去)。

### 3.4 一時診断の撤去

`diagnostics-session2-kept.patch` の環境変数分岐(`TENES_FU_SIGN_FIX`, `TENES_FU_FORBIDDEN_TOL`,
`TENES_CTM_RDM_ABS`, `TENES_CTM_RDM_LOG`, `TENES_CTM_RDM_PHASE`, `TENES_CTM_PHASE_FIX`,
`TENES_CTM_PHASE_LOG`)はすべて撤去し、上の 3.1〜3.3 を無条件の実装に置き換える。
`FORBIDDEN_TOL`(forbidden 閾値)は別課題として残るが、環境変数での上書きはしない。

## 4. 影響と非目標

- 測定値は位相に依らないので既存 golden(fermion / bosonic / finite-T)は変わらない。
  bosonic と finite-T は `phase_invariant = false` で現行コードパスを通るので**ビット同一**
  (§5 の A/B で機械確認)。
- fermion では位相が 1 のケース(D=2 実数の既存 golden `FreeFermion`)は C1 に触らないが、
  判定側で ρ_i を位相 1 で割る演算が入るので `rdm_dist` の最下位ビットが変わりうる。
  収束値は同じで golden は許容内。ビット同一は要求しない。
  位相 ≠ 1 のケース(複素、D≥3)は CTM の反復回数が減る(判定が効くようになる)だけで
  収束値は同じ。
- forbidden 閾値(`build_full_update_environment` の 1e-8)の扱いは本設計の外。
- T5 の D=3 化は本設計の契約に含める(§5)。

## 5. 検証(契約書 `2026-09-04-fermion-full-update-contract.md` に T8 として追記)

| # | 検査 | 形式 |
|---|---|---|
| T8-i | `Check_Convergence_CTM_RDM(phase_invariant = true)` の位相不変性。独立な履歴を持つ二系列で 2 回ずつ呼ぶ: 基準 E0 → E1、変換 e^{iφ0}E0 → e^{iφ1}E1(全 C1 に掛ける。φ0 ≠ φ1、実数では −1/+1、複素では任意)。初回はどちらも未収束(`rdm_dist` NaN)、2 回目の `rdm_dist` と bool が一致(相対 1e-12)。零 RDM・NaN・Inf を含む環境では `valid = false`(`rdm_dist` NaN、false)。`phase_invariant = false` では位相を掛けた系列が未収束のまま(現行挙動)。real/complex、`is_density` 両方 | doctest |
| T8-ii | `fix_environment_phase` の前後で: (a) one-site norm が全サイトで実正(`|Im| ≤ 1e-10·|n|`、`Re > 0`)、(b) two-site 水平/垂直の ⟨O⟩/⟨1⟩ と FU の N_plain 経由の ⟨G⟩(T2-vi の形)が不変(相対 1e-12)、(c) 未正規化の norm と N_plain は `conj(phase)` 倍、(d) C2/C3/C4/eT は byte 不変、(e) phase = 1 なら C1 も byte 不変、(f) 複素位相の後で N_plain が Hermitian(相対 1e-12)かつ PSD。状態: 2×2 実数 D=3、2×2 複素 D=2、**3×2** 複素 D=2(いずれもランダム偶 Tn)。CTM は `count < iteration_max` を前提アサート | doctest |
| T8-iii | 複素 fermion E2E(`is_real = false`、自由フェルミオン D=2、SU のみ)で「CTM did not converge」が出ない。エネルギーの虚部 ≤ 1e-10。密度が既存 `FreeFermion` の許容内 | E2E |
| T8-iv | FU 核のガード: T3-i の入力で (a) C1 を −1 倍 → `std::runtime_error`、(b) 複素で C1 に e^{iφ}(φ = 2.0)→ `std::runtime_error`(虚部検査)、(c) 変えない入力では投げない | doctest |
| T2-vi(complex) | 既存 T2-vi の complex_tensor 版。環境は `Calc_CTM_Environment_density(phase_invariant = true)` → `fix_environment_phase` で作り、`count < iteration_max` を前提アサート | doctest |
| T5(改訂) | 自由フェルミオン **D=3、chi=12**、`iteration_max = 200`、SU 1000 step 後に FU 1〜4 sweep で `E_FU < E_SU − 1e-4` かつ厳密値に接近、「CTM did not converge」が出ない。D=2 の非増加要求は削除 | E2E |
| T7 | 既存テストが無変更で緑 | 既存 |
| A/B(検証者) | 変更前後のバイナリで、bosonic golden(`AntiferroHeisenberg_real/complex`、`Honeycomb`、`J1J2_AFH`、`Kitaev`)と finite-T golden(`FT_*`)を同一入力・1 thread・MPI off・固定 seed で走らせ、`output_*/*.dat` の sha256 が一致。fermion `FreeFermion`(位相 1)は許容内一致 | Claude が実施 |

## 6. タスク分割

1 タスク(Codex 1 バッチ): 3.1〜3.4 を実装し、`ctest --preset gcc-release` 全件緑。
テストは契約から別担当が先に書く(RED を確認してから実装)。
