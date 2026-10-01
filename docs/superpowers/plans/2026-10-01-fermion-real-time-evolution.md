# fermion 実時間発展 Implementation Plan

> **For agentic workers:** この計画は CLAUDE.md の多段手順で実行する。テストはテスト作成者が契約書から書き、実装は Codex が行い、Claude が独立に検証する。Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** fermion モードで `mode = "time"`(実時間発展)を SU・FU の両方で使えるようにし、実数で保存した状態を complex の計算で load できるようにする。

**Architecture:** solver の入力ガードを有限温度だけの拒否に狭め、実時間 SU のスイープ後にパリティ台帳の不変条件を検査する。テンソルの読み込みでファイルの値の型を読み先と照合し、実数→complex は変換して読み、complex→実数はエラーにする。実時間発展そのものの計算経路は既存の fermion SU/FU/測定をそのまま使う。

**Tech Stack:** C++17(mptensor v0.5.0、doctest)、Python 3.9+(tenes_simple / tenes_std、E2E テスト)、CMake/ctest。

**Spec:** `docs/superpowers/specs/2026-10-01-fermion-real-time-evolution-design.md`
**Contract(テスト作成者向け):** `docs/superpowers/specs/2026-10-01-fermion-real-time-evolution-contract.md`

## Global Constraints

- 新しい E2E テストは 1 本 30 秒以内(`OMP_NUM_THREADS=1`)。
- tenes の同時実行は 2 本まで、`OMP_NUM_THREADS=1`。実験は `work/fermion-real-time/<テーマ>/` で、CWD をそこにして行う。
- テストファイルは実装担当(Codex)が変更しない。`test/` 全体のスナップショットとの diff で検査する。
- formatter は Codex に実行させない。C++ の整形は Claude が `git add` 後に `git clang-format`(変更行だけ)で行う。Python は black(line-length 88)。
- MPI 有効ビルドで、例外を rank 0 だけで投げない(全ランクで揃える)。
- 有限温度、奇パリティのプロダクト初期状態、長距離ボンドの FU は引き続き拒否する。

## Review Focus

- complex で保存して実数で読む組み合わせは、ボゾンでも明示的な `load_error` になるべきで、黙って読んではならない(Task 2 のテスト §5-3)。
- MPI 2 ランクでの実数→complex 変換。分散配置が保存時と読み込み時で違っても値が一致すること(Task 2 のテスト §5-5)。
- mptensor 0.2 形式(ヘッダに `value_type` がない古い保存)は従来どおり読めること。既存の checkpoint 系テスト全件で確かめる(Task 2 の検証手順)。
- 長距離ボンドのゲート列を含む実時間 SU で、スイープ後の台帳が元に戻ること。正しい実装では発火しない検査なので、長距離ボンドの実時間 SU が完走することで確かめる(Task 1 の検証手順で、P0 ビルドとの比較を含めて Claude が手で確認)。
- 時刻 0 の測定が、load した状態をそのまま測ること(実時間発展の前に状態を変えない)。テスト §5-1 が押さえる。

---

### Task 0: テストの作成(テスト作成者)

**Files:**
- Create: `test/fermion/*.py.in`(E2E、契約書 §2〜§5)、必要なら `test/fermion/*.cpp`(§6)
- Modify: `test/CMakeLists.txt`(登録の追加のみ)
- Report: `work/fermion-real-time/tests/REPORT.md`

- [ ] **Step 1:** テスト作成者(Claude サブエージェント)に契約書だけを渡して dispatch する。計画書・設計書のテスト案のコードは渡さない。依頼文に、作業ツリーの変更禁止(`git stash`/`checkout` 含む)、tenes 同時 2 本、E2E 30 秒、報告ファイル必須を明記する。
- [ ] **Step 2:** 報告ファイルの存在を確認し、RED を 1 件ずつ「正しい理由で落ちているか」確かめる(実装前に期待される理由は契約書の各節に書いてある)。
- [ ] **Step 3:** 検出力の確認(契約書 §7)の結果を読み、赤くならない変異があればテスト作成者に差し戻す。
- [ ] **Step 4:** `test/` 全体のスナップショットを取る:`tar cf work/fermion-real-time/tests/snapshot.tar test && shasum -a 256 work/fermion-real-time/tests/snapshot.tar`、および `git stash` を使わずに `git diff --stat -- test` と未追跡ファイル一覧を記録する。
- [ ] **Step 5:** テストをコミットする(テストは実装前なので ctest の該当分は赤い。コミットメッセージにその旨を書く)。

### Task 1: 入力ガードの縮小と実時間 SU の台帳検査(Codex)

**Files:**
- Modify: `src/iTPS/load_toml.cpp:656-658`
- Modify: `src/iTPS/simple_update.cpp:226-233`、`src/iTPS/iTPS.hpp`(宣言追加)
- Modify: `src/iTPS/time_evolution.cpp:45-58`

**Interfaces:**
- Produces: `void iTPS<tensor>::check_fermion_phys_ledger_restored() const;`(fermion 無効時は何もしない。違反時は既存と同じ文言の `std::logic_error`)

- [ ] **Step 1:** ガードを変える。

```cpp
  if (peps_parameters.calcmode == PEPS_Parameters::finite_temperature) {
    throw_fermion_guard("finite-temperature mode");
  }
```

- [ ] **Step 2:** `simple_update()` のスイープ後の検査を `check_fermion_phys_ledger_restored()` に切り出す。

```cpp
template <class tensor>
void iTPS<tensor>::check_fermion_phys_ledger_restored() const {
  if (!finfo.enabled) {
    return;
  }
  for (int site = 0; site < N_UNIT; ++site) {
    if (finfo.phys[site] != peps_parameters.phys_parity[site]) {
      throw std::logic_error(
          "fermion simple update invariant violated: physical ledger "
          "did not return to its original value after a gate sweep");
    }
  }
}
```

`simple_update()` の該当ブロックをこの呼び出しに置き換える。

- [ ] **Step 3:** `time_evolution()` の SU 経路で、各ステップの全ゲート適用後(`t += dt;` の直前)に `if (su) { check_fermion_phys_ledger_restored(); }` を入れる。
- [ ] **Step 4:** ビルドし、契約書 §2・§3・§4・§6 のテストを回す(`ctest -R <名前>`)。
- [ ] **Step 5:** 報告ファイル `work/fermion-real-time/impl/T1-REPORT.md` を書く(成果物。無いこと自体が欠陥)。
- [ ] **Step 6(Claude):** テスト不改変を検査し、ctest 全件(Debug、`OMP_NUM_THREADS=1`)と MPI ビルドの fermion 系を回す。長距離ボンドを含む fermion の実時間 SU(例:正方格子 t' あり、SU 20 ステップ)を手で走らせ、完走すること、P0 ビルドと同じ出力になることを確かめる。`git clang-format` で整形してコミットする。
- [ ] **Step 7:** タスクレビュー(新規サブエージェント。変異テストを推奨)。

### Task 2: 実数で保存したテンソルの complex への読み込み(Codex)

**Files:**
- Modify: `src/iTPS/saveload_tensors.cpp`(`load_tensor()` と `load_tensors_v0()` 内のラムダ `load`)
- 必要なら Modify: `src/tensor.hpp` / `src/tensor.cpp`(変換関数を置く場合)

**Interfaces:**
- Produces(ファイル内の無名名前空間で可):
  - `int read_saved_value_type(std::string const &path, MPI_Comm comm);` — rank 0 がベースファイルのヘッダを読み、`value_type=` の値(0 = double、1 = complex)を返す。ヘッダが `mptensor` で始まらない古い形式なら −1。結果は全ランクに bcast する。
  - `template <class ptensor> void load_tensor_file(ptensor &A, std::string const &path);` — 型を照合して読む。`load_tensor()` と `load_tensors_v0()` の両方がこれを使う。

- [ ] **Step 1:** `read_saved_value_type` を書く。ヘッダは次の形(`saved/T_0.dat` の実例):

```
mptensor 0.5.0
matrix_type= 1 (LAPACK)
value_type= 0 (double)
comm_size= 1
...
```

- [ ] **Step 2:** `load_tensor_file` を書く。

  - 読み先の型と一致、または −1:従来どおり `A.load(path)`。
  - ファイルが 0(double)で読み先が complex:`real_tensor R(A.get_comm()); R.load(path);` で読み、complex に変換して `A` に入れる。
    変換は分散配置に依存しない方法で行う:全要素数の `std::vector<double>` を 0 で用意し、`R` の各局所要素 n について `R.global_index_fast(n, idx)` から行優先の通し番号を計算して値を置き、`MPI_Allreduce`(和)で全ランクに揃える(各要素はちょうど 1 つのランクが持つので和は厳密)。
    `complex_tensor C(A.get_comm(), R.shape(), R.get_upper_rank());` を作り、`C` の各局所要素に通し番号で値を入れる。`_NO_MPI` ビルドでも通るように、`src/mpi.hpp` の `allreduce_sum(std::vector<double>&, MPI_Comm)` を使う。
  - ファイルが 1(complex)で読み先が実数:全ランクで `tenes::load_error` を投げる。メッセージ例:
    `"ERROR: " + path + " holds a complex tensor, which cannot be loaded into a real-valued calculation (parameter.general.is_real = true). HINT: set is_real = false, or load a checkpoint saved with is_real = true."`

- [ ] **Step 3:** `load_tensor()` の `temp.load(filename.c_str());` と、`load_tensors_v0()` のラムダの `A.load(path.c_str());` を `load_tensor_file` 経由にする。`load_tensor()` の後段(rank の照合、`resize_tensor`)は変えない。
- [ ] **Step 4:** ビルドし、契約書 §5 のテストを回す。既存の checkpoint 系テスト(`ctest -R "checkpoint|SaveLoad|restart|DExpand"`)も回す。
- [ ] **Step 5:** 報告ファイル `work/fermion-real-time/impl/T2-REPORT.md`。
- [ ] **Step 6(Claude):** テスト不改変の検査、ctest 全件(Debug serial と MPI ビルド)、整形、コミット。
- [ ] **Step 7:** タスクレビュー。

### Task 3: ドキュメントと NEWS(Claude)

**Files:**
- Modify: `docs/sphinx/ja/file_specification/parameter_section.rst:76`、`docs/sphinx/en/file_specification/parameter_section.rst:75`
- Modify: `tensor_load` の説明箇所(`parameter_section.rst` の general 節、ja/en)
- Modify: `NEWS.md`

- [ ] **Step 1:** fermion の記述を「基底状態計算と実時間発展に対応」に改め、非対応の列挙から実時間発展を外す。実時間発展の初期状態(偶パリティのプロダクト状態、ランダム偶初期化、`tensor_load`)と、プロダクト状態から大きな D で FU を始めるとパリティ射影で止まることがある点(回避策:小さい D から始める、または SU で発展させる)を書く。
- [ ] **Step 2:** `tensor_load` の説明に、実数で保存したテンソルを complex の計算で読めること、逆はエラーになることを書く。
- [ ] **Step 3:** NEWS.md の fermion の項目に実時間発展の対応を足し、Limitations から real-time evolution を外す。load の修正は独立の項目にする(ボゾンでも、従来 Release で壊れた値を黙って読んでいた組み合わせが正しく読めるようになった)。PR 本文と NEWS には、このブランチ内の経緯ではなく develop に対する差分だけを書く。
- [ ] **Step 4:** Sphinx のビルドで警告が増えないことを確かめ、コミットする。

### Task 4: 全ブランチレビュー

- [ ] **Step 1:** 最上位モデルのサブエージェントで、develop に対するブランチ全体をレビューする(作業ツリーの変更禁止、変異テスト推奨)。
- [ ] **Step 2:** 指摘に対応し、ctest 全件を回す。
- [ ] **Step 3:** PR を作る(本文に「Generated with Claude Code and Codex」)。
