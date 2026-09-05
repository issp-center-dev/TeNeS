# フェルミオン CTM 環境の全体位相の固定 — 実装台帳

設計: `docs/superpowers/specs/2026-09-05-fermion-ctm-phase-design.md`
契約: `docs/superpowers/specs/2026-09-04-fermion-full-update-contract.md` §2.5、§3.6、§3.5 T5(改訂)
作業記録: `work/fermion/ctm-phase/`

| 手順 | 担当 | 状態 |
|---|---|---|
| 設計書 | Claude | 改訂 2.1(改訂 1 はテスト作成者の計測で棄却) |
| 設計レビュー | Codex | 済(改訂 1: `codex-review-design.md` 11 件、改訂 2: `codex-review-design2.md` 7 件、いずれも反映) |
| 契約追記 | Claude | 済(改訂 2.1: 49a24ae4, a5e287f6, 2632fd36) |
| テスト作成(T8-i〜v、T2-vi complex、T5 改訂) | テスト作成者 | 済(76249cd4、11 TEST_CASE + 2 E2E、報告 `report-testauthor.md`) |
| RED 確認・スナップショット | Claude | 済(E2E 2 件を現行バイナリで RED 確認、`test-snapshot.sha256`) |
| 変更前バイナリの A/B 基準出力 | Claude | 済(`TeNeS-ab-before` worktree、`work/fermion/ctm-phase/ab/before/sha256.txt`、72 ファイル) |
| 実装(タスク 1) | Codex | 進行中(2026-09-05 09:45〜、`codex-task1-prompt.txt`) |
| 独立検証(ctest 全件、A/B sha256、D=3 FU、複素) | Claude | — |
| レビュー | サブエージェント | — |
| 整形・コミット | Claude | — |

## タスク 1: 位相不変判定・N の位相正規化・RDM 出力(1 バッチ、設計書改訂 2)

対象ファイル: `src/iTPS/core/ctm.cpp`、`src/iTPS/core/ctm.hpp`、`src/iTPS/core/ctm_single.cpp`、
`src/iTPS/iTPS.cpp`、`src/fermion/full_update_env.hpp`、`src/iTPS/onesite_obs.cpp`。

1. `Check_Convergence_CTM_RDM` に `const bool phase_invariant = false` を追加(設計 §3.1)。
   健全性検査(finite、`trace_abs > 0`)は全経路。`phase_invariant == false` の経路は現行と同一の
   演算順序を保つ(bosonic / finite-T のビット同一)。
2. `Calc_CTM_Environment_density` に `bool phase_invariant = false` を追加し、1 に渡す。
   `Calc_CTM_Environment`(bosonic)は変更しない。
3. `build_full_update_environment` で N の位相正規化と窓ノルムのガード(設計 §3.2)。
   `full_update_environment::phase` を追加。
4. `iTPS::update_CTM()` の fermion 分岐: `Calc_CTM_Environment_density(..., true, true)`。
   `update_CTM_density()` と bosonic 分岐は触らない。
5. `measure_onesite_rdm` の fermion CTM 経路で RDM の位相を trace の位相で正規化(設計 §3.4)。
6. `Full_update_bond_fermion` は変更しない。環境変数分岐は入れない。
7. 作業ツリーの一時診断は dispatch 前に Claude が撤去済みであること(`git status --short src/` が空)。

完了条件: `cmake --build --preset gcc-release` が通り、`ctest --preset gcc-release` で
`test_fermion_layer`(T8 を含む)、`FreeFermionComplex`、`FreeFermionFull`(D=3)が緑、
既存テストが緑。`test/` は変更しない。報告 `work/fermion/ctm-phase/report-task1.md`。
