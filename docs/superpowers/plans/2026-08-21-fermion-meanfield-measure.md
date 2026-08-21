# フェルミオン模式の平均場環境対応 実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** フェルミオン模式で `meanfield_env = true` を使えるようにし、2サイト観測量を
D¹²d⁴ の blob ではなく D⁶d² の単層 graded 縮約で評価する。

**Architecture:** ket 層の2サイト状態 `ψ_AB = tensordot_f(A_λ, B_λ)` を組み、
`trace_f(conj_f(ψ_AB), O ∘ ψ_AB)` / `trace_f(conj_f(ψ_AB), ψ_AB)` で期待値を取る。
符号は `fops.hpp` の graded 演算だけが生成し、手書きの swap は置かない。λ の dressing は
既存のボソン MF ブロックをそのまま使う。検証は Fock オラクル(開いた脚の任意ラベル固定に
拡張)を真とし、d=4 の wrap 規約と両端測定、iTPS レベルの配線、E2E で三角測量する。

**Tech Stack:** C++17(mptensor、doctest)、Python 3(numpy、pytest)、CMake + CTest

**Spec:** `docs/superpowers/specs/2026-08-21-fermion-meanfield-measure-design.md`

## Global Constraints

- 役割分担は `~/.claude/CLAUDE.md` の多段エージェント方式。テストはテスト作成者(契約書のみ渡す)、
  実装は Codex(テストファイル変更禁止)、検証・整形・コミットは Claude
- 実装者は `apply_swap` を手で呼ばない。符号は `tensordot` / `trace` / `conj` / `transpose`
  (`src/fermion/fops.hpp`)だけが生成する
- CTM の blob 経路(`build_reduced_pair` の値、`contract_reduced_pair_*_density_CTM`)と
  ボソン経路の数値は変えない。既存の ctest 29 件と `test/data/output_*/` の golden に触らない
- 利用者向け文言に `M1` / `M2` を書かない
- ビルド: `cmake --build out-gcc/build -j 8`、C++ テスト: `out-gcc/build/test/test_fermion_layer`、
  `out-gcc/build/test/test_input`、pytest: `venv/bin/python3 -m pytest test/python`、
  全件: `ctest --test-dir out-gcc/build`
- 整形は Claude がコミット直前に、そのタスクで変更したファイルだけに `clang-format` / `black` を掛ける
- 各タスクの実装者は報告ファイル `work/mf-task<N>-report.md` を必ず作る(存在を Claude が確認する)

---

## ファイル構成

| ファイル | 役割 | 変更 | タスク |
|---|---|---|---|
| `test/fermion/fock_oracle.py` | 開いた脚のラベル固定、`mf_sum`、自己検査、MF 参照値出力 | 変更 | T1 |
| `test/python/test_fock_oracle_mf.py` | オラクル拡張の自己検査を pytest に載せる | 新規 | T1 |
| `src/fermion/reduced.hpp` | `build_pair_state` / `apply_pair_op`(公開)、`build_reduced_pair` がそれを使う | 変更 | T2 |
| `src/fermion/reduced_measure.hpp` | `contract_pair_MF`(norm 版・演算子版) | 変更 | T2 |
| `test/fermion/mf_measure.cpp` | 層1〜3(API vs オラクル、d=4 wrap 規約、線形性) | 新規 | T2 |
| `test/test_fermion_layer.cpp` | 末尾に `#include "fermion/mf_measure.cpp"`、層4(iTPS レベル) | 変更 | T2, T3 |
| `src/iTPS/twosite_obs.cpp` | MF 分岐のフェルミオン経路 | 変更 | T3 |
| `src/iTPS/load_toml.cpp` | `MeanField_Env` ガード削除 | 変更 | T3 |
| `test/input.cpp` | 「rejects mean-field environment」→「accepts」 | 変更 | T3 |
| `test/fermion/free_fermion_mf.py.in` | E2E(MF vs CTM vs 厳密解、golden) | 新規 | T4 |
| `test/CMakeLists.txt` | `FreeFermionMF` 登録 | 変更 | T4 |
| `test/data/output_FreeFermionMF/free_fermion_mf.dat` | golden(初回 GREEN 時に Claude が生成) | 新規 | T4 |
| `docs/sphinx/{ja,en}/file_specification/parameter_section.rst` | 非対応一覧と `meanfield_env` の注記 | 変更 | T4 |
| `NEWS.md` | 項目追加 | 変更 | T4 |

`test/fermion/mf_measure.cpp` を `test_fermion_layer.cpp` から `#include` するのは、
`r2_convention.cpp` と同じ流儀(`test_fermion_layer.cpp` 末尾)。同じ翻訳単位なので
`r2_convention.cpp` の `static` ヘルパ(`make_r4_tensor`、`r4_hop_plain`、`r4_wrap`、
`r2_expect_two`、`r2_norm`)がそのまま使える。`mf_measure.cpp` は `r2_convention.cpp` の
**後**に include する。

---

## Task 1: Fock オラクルの拡張(テスト基盤)

**Files:**
- Modify: `test/fermion/fock_oracle.py`
- Create: `test/python/test_fock_oracle_mf.py`

**Interfaces:**
- Produces: `Oracle(patch, tensors, leg_parities, dangling_labels=None)`、
  `make_case_open(lx, ly, parity, seed=0)`、`mf_sum(patch, tensors, leg_parities, fn)`、
  `Oracle.density_density(i, j)`、`print_mf_case(name, lx, ly, parity, seed, lam)`。
  T2 の C++ テストはこの `print_mf_case` の出力を定数として埋める

担当: テスト作成者(テスト基盤なので実装者には渡さない)。Claude が自己検査を独立に確認する。

### 契約書(テスト作成者に渡す散文)

**背景.** `test/fermion/fock_oracle.py` は、正方格子の小さなパッチ(lx×ly)を Fock 空間で
厳密に構成するオラクルである。各サイトの5脚テンソル `T[l, t, r, b, s]` の添字ラベルは
占有数(パリティ)に対応し、内部ボンドは「ボンド生成子」で両端の補助モードを同時に生成し、
各サイトの射影子がテンソル要素を係数として補助モードを消滅・物理モードを生成する
(`OP_ORDER = (b, r, t, l)` の順)。現在、パッチの外に出る脚(開いた脚)は偶ラベル 0 にしか
固定できない(`apply_physical_projector` で該当モードが無いため奇ラベルの項は落ちる)。

**拡張 1: 開いた脚のラベル固定.** `Oracle.__init__` に省略可能な引数
`dangling_labels` を足す。`{(site, leg_name): label}` の辞書で、パッチの全ての開いた脚
(`leg_name` は `"l"/"t"/"r"/"b"`)を網羅していなければ `ValueError`。与えられた場合:

- 開いた脚ごとに補助モードを 1 つ割り当てる(モード番号は物理・内部ボンドの後ろに、
  `(site, leg_name)` のソート順で)。
- `state()` は真空から始め、**まず** 固定ラベルが奇の開いた脚の補助モードを上記ソート順に
  `apply_create` で生成し、次に内部ボンドの生成子、最後にサイトの射影子を掛ける。
- `apply_physical_projector` は、開いた脚の添字が固定ラベルに等しいテンソル要素だけを
  使う(他は飛ばす)。奇ラベルの開いた脚は補助モードが存在するので通常どおり消滅される。
- `physical_state()` の「補助モードが全て空の成分だけ残す」は変更しない。
- `dangling_labels=None` のときの挙動は現在と完全に同じ(既存テスト `self_check` /
  `task9_check` / R2 参照値の再現が壊れないこと)。

**拡張 2: 密度相関.** `Oracle.density_density(i, j)` を足す:
`⟨ψ| n_i n_j |ψ⟩`(非規格化、`one_body` と同じ流儀で `physical_state()` に作用)。

**拡張 3: 開いた脚が多次元のケース生成.** `make_case_open(lx, ly, parity, seed=0)` を足す。
`make_case` と同じだが、開いた脚にも `parity`(例 `[False, True]`)を使う。テンソルは
`deterministic_tensor(site, leg_parities[site], seed)` のまま(偶パリティ要素だけ非零)。

**拡張 4: MF 和.** モジュール関数
`mf_sum(patch, tensors, leg_parities, fn)`:
開いた脚の全ラベル割り当て `x`(各脚はその脚のパリティベクトルの長さ分)について
`Oracle(patch, tensors, leg_parities, dangling_labels=x)` を作り、`fn(oracle)` を足し合わせて
返す。λ の重みは **取らない**(テンソル側に掛けてから渡す)。MF 期待値は
`mf_sum(..., lambda o: o.one_body(i, j)) / mf_sum(..., lambda o: o.norm())` である。

**拡張 5: 参照値出力.** `print_mf_case(name, lx, ly, parity, seed, lam)`:
`make_case_open` でケースを作り、**開いた脚それぞれに** 長さ `len(parity)` のベクトル
`lam`(例 `[1.0, 0.7]`)を添字ラベルごとに掛けてから(`np.multiply` で該当軸に沿って)、
次を 17 桁で印字する: `norm`(= mf_sum の norm)、各サイトの `n{i}`(規格化済み)、
各内部ボンドの `hop{a}{b}`(= (one_body(a,b)+one_body(b,a)) / norm)、`pair{a}{b}`
(= (pairing(a,b)+pairing(b,a)) / norm)、`nn{a}{b}`(= density_density(a,b) / norm)。
`main()` から、次のケースを印字する:
`mf_horizontal_2site`(2,1,[F,T],seed 0)、`mf_vertical_2site`(1,2,[F,T],seed 0)、
`mf_seed173_horizontal_2site`、`mf_seed173_vertical_2site`、`mf_single_site`(1,1,[F,T],seed 0)、
`mf_seed173_single_site`。いずれも `lam = [1.0, 0.7]`。

**自己検査(`test/python/test_fock_oracle_mf.py` に pytest として書く。
`sys.path` に `test/fermion` を挿して `import fock_oracle`).**

1. `dangling_labels` を省略した `Oracle` は既存 `make_case(2,1,[F,T])` の `observables()` を
   ビット単位で再現する(既存 R2 参照値 `norm = 5.46520500622191588e-04`、
   `n0 = n1 = 5.35874068543760670e-02` と 1e-12 で一致)。
2. 開いた脚が全て 1 次元偶(`make_case`)なら `mf_sum` は割り当てが 1 通りしかなく、
   `observables()` と一致する。
3. 全パリティ偶(仮想 `[False, False]`、物理 `[False, True]` のまま、`make_case_open`)では、
   `mf_sum` の norm が numpy の素の縮約 `Σ_{全添字} |Σ_{内部ボンド} T_A T_B|²` と 1e-12 で一致する
   (ボソン極限)。
4. 奇ラベルの開いた脚の補助モード生成順を逆順にしても(テスト用にモンキーパッチ等で順序を
   差し替える)、`norm`、`n_i`、`one_body(0,1)`、`pairing(0,1)`、`density_density(0,1)` の
   `mf_sum` が 1e-12 で変わらない。
5. `mf_horizontal_2site` と `mf_vertical_2site`(seed 0、`lam=[1.0, 0.7]`)で
   `|hop01| / norm > 1e-3` である(符号を検出できるデータセットであることの担保。
   既存 R3 データセットは真空射影のため `hop = 0`)。
6. `dangling_labels` の網羅漏れ(開いた脚を 1 つ欠く)は `ValueError`。

### 手順

- [ ] **Step 1:** テスト作成者に上の契約書を渡し、`fock_oracle.py` の拡張と
  `test/python/test_fock_oracle_mf.py` を書かせる
- [ ] **Step 2:** Claude が `venv/bin/python3 -m pytest test/python/test_fock_oracle_mf.py -v` と
  `venv/bin/python3 test/fermion/fock_oracle.py` を実行し、自己検査 1〜6 が通ること、
  既存の R2 参照値出力(`horizontal_2site.norm` 等)が 1 桁も変わっていないことを確認する
  (`git stash` は使わず、変更前の出力は `git show HEAD:test/fermion/fock_oracle.py > <scratch>` で
  別ファイルとして実行して取る)
- [ ] **Step 3:** Claude が独立確認: 自己検査 4(生成順不変)を、コードを読んで
  「生成順は x ごとの全体符号にしかならない」ことと整合しているか確認する。
  `mf_single_site` の `n0` を手計算可能な範囲で検算する(1×1、開いた脚 4 本 × D=2 の 16 通りの和)
- [ ] **Step 4:** `black test/fermion/fock_oracle.py test/python/test_fock_oracle_mf.py`、
  `venv/bin/python3 -m pytest test/python -q` 全件緑を確認してコミット
  `"Extend the Fock oracle to mean-field sums over open legs"`

---

## Task 2: 単層縮約 API(`build_pair_state` / `apply_pair_op` / `contract_pair_MF`)

**Files:**
- Modify: `src/fermion/reduced.hpp`(`build_reduced_pair` の前に公開関数 2 つ、`build_reduced_pair` の書き換え)
- Modify: `src/fermion/reduced_measure.hpp`(`contract_pair_MF` 2 オーバーロード)
- Create: `test/fermion/mf_measure.cpp`
- Modify: `test/test_fermion_layer.cpp`(末尾に `#include "fermion/mf_measure.cpp"`)

**Interfaces:**
- Consumes: `fops.hpp` の `tensordot` / `trace` / `conj` / `transpose` / `slice` /
  `wrap_twosite_gate` / `wrap_reduced_pair_op`、`reduced.hpp` の `reduced_pair_direction`
- Produces(T3 が使う):
  ```cpp
  template <class tensor>
  ftensor<tensor> build_pair_state(const ftensor<tensor>& TnA, const ftensor<tensor>& TnB,
                                   reduced_pair_direction direction);           // reduced.hpp
  template <class tensor>
  ftensor<tensor> apply_pair_op(const ftensor<tensor>& pair, const ftensor<tensor>& op12); // reduced.hpp
  template <class tensor>
  typename tensor::value_type contract_pair_MF(const ftensor<tensor>& pair);                // reduced_measure.hpp
  template <class tensor>
  typename tensor::value_type contract_pair_MF(const ftensor<tensor>& pair,
                                               const ftensor<tensor>& op12);                // reduced_measure.hpp
  ```

### 契約書(テスト作成者に渡す散文)

**対象.** `tenes::fermion` 名前空間の 4 関数(上の Interfaces のシグネチャ)。
`build_pair_state` と `apply_pair_op` は `src/fermion/reduced.hpp`、`contract_pair_MF` は
`src/fermion/reduced_measure.hpp` に置かれる。テストは `test/fermion/mf_measure.cpp` に
doctest の `TEST_CASE` として書き、`test/test_fermion_layer.cpp` の末尾
(`#include "fermion/r2_convention.cpp"` の**直後**)に `#include "fermion/mf_measure.cpp"` を足す。
同じ翻訳単位なので `r2_convention.cpp` の `static` ヘルパ(`make_r4_tensor`、`r4_hop_plain`、
`r4_nn_plain`、`r4_wrap`、`r2_expect_two`、`r2_norm`、`r2_hop_op` 等)と
`test_fermion_layer.cpp` の `ft` 型エイリアス・`make_random_ft` が使える。

**記法.** 以下、`tensordot_f` / `trace_f` / `conj_f` / `transpose_f` は `tenes::fermion` 名前空間の
graded 版(`src/fermion/fops.hpp`)を指す。`mptensor::` の素の関数とは区別する。

**振る舞い.**

- `build_pair_state(A, B, horizontal)` は `tensordot_f(A, B, Axes(2), Axes(0))`、脚順
  `(l_A, t_A, b_A, s_A, t_B, r_B, b_B, s_B)`。`vertical` は `tensordot_f(A, B, Axes(3), Axes(1))`、
  脚順 `(l_A, t_A, r_A, s_A, l_B, r_B, b_B, s_B)`。A・B が 5 脚でなければ `std::runtime_error`。
  戻り値のパリティは `tensordot_f` が返すもの(自由脚の順に連結)。
- `apply_pair_op(pair, op12)` は `op12` の脚 `(in_A, in_B, out_A, out_B)` を `pair` の脚 3・7 に
  作用させ、結果を元の脚順に戻す(`tensordot_f(pair, op12, Axes(3,7), Axes(0,1))` の後
  `transpose_f(·, Axes(0,1,2,6,3,4,5,7))`)。`pair` が 8 脚・`op12` が 4 脚でなければ
  `std::runtime_error`。
- `contract_pair_MF(pair)` = `trace_f(conj_f(pair), pair, 全脚, 全脚)`。
- `contract_pair_MF(pair, op12)` = `trace_f(conj_f(pair), apply_pair_op(pair, op12), 全脚, 全脚)`。
- **値の不変条件**: `build_reduced_pair`(blob 経路)の数値は変わらない
  (既存 R3/R5 テストがそのまま通る)。

**テストすべきこと(層1〜3).** すべて `doctest::Approx(...).epsilon(1e-12)`(相対)か、
値が 0 に近いものは絶対 1e-12。

層1(d=2、オラクル参照値): T1 の `print_mf_case` が印字した `mf_horizontal_2site` /
`mf_vertical_2site` / `mf_seed173_*`(`lam = [1.0, 0.7]`)を定数として埋める。C++ 側は
Python の `deterministic_tensor(site, parities, seed)` と同じ式
(`x = (site+2)*(1+seed+Σ_{ax<5}(ax+3+seed%5)*idx[ax])`、`0.19 sin x + 0.13 cos 0.37x`、
偶パリティ要素のみ)で 5 脚 D=2・d=2 の `ft` を作り(`r2_convention.cpp` の `make_r2_tensor` と
同じ式だが、`r2_parities` が境界脚を 1 次元偶にするのに対し全脚を `{false, true}` にした版を
新たに書く)、開いた脚に `lam` を `multiply_vector` で掛け、
`build_pair_state` → `contract_pair_MF` で `norm`、`n0`、`n1`(`r2_density_a_op` /
`r2_density_b_op` を演算子として渡す)、`hop01`(`r2_hop_op`)、`pair01`(`r2_pair_op`)、
`nn01`(要素 `(1,1,1,1)=1` の rank-4)を出し、参照値と一致すること。水平・垂直・seed 2 つ。
演算子は `wrap_twosite_gate(op, p, p)` で包む(d=2 ではこれは `op` に等しいので、
`ftensor{op, {p,p,p,p}}` で作っても同じ。どちらでもよいが、層2 と同じヘルパを使うこと)。

層2(d=4、wrap 規約と両端測定): R5 と同じ `make_r4_tensor(2,1,0,seed)` /
`make_r4_tensor(2,1,1,seed)`(水平)と `(1,2,·)`(垂直)、seed `{0, 173}` で `pair` を組む。
  - `contract_pair_MF(pair, r4_wrap(r4_hop_plain(), true, false)) / contract_pair_MF(pair)`
    が R5 の `hop_ref`(= `r2_expect_two(psi, 3, 7, r4_wrap(r4_hop_plain(), true, false))`、
    `psi` は同じ 2 テンソルの graded tensordot)と一致する。`nn` も同様に `r4_nn_plain()` で一致。
  - 陰性対照 2 件: `r4_wrap(hop, true, true)`(両 swap)と `r4_wrap(hop, false, false)`(wrap なし)
    を渡すと `hop_ref` と **一致しない**(`|差| > 0.1·max(1e-6, |hop_ref|)`)。
  - 両端測定: `o_B = transpose_f(r4_wrap(r4_hop_plain(), true, false), Axes(1,0,3,2))` を渡した
    値が `hop_ref` と一致する(ホッピングはサイト交換対称なので、B 側から同じ行列を与えて
    graded transpose した演算子は A 側のそれに等しい)。陰性対照: graded でない
    `mptensor::transpose(op.t, Axes(1,0,3,2))` でテンソル部だけ入れ替えた `ft`(パリティは同じ)
    を渡すと一致しない。
  - `contract_pair_MF(pair)` が `r2_norm(psi)` と一致する。

層3(線形性・奇ラベルの開いた脚): 層1 のテンソル(d=2、開いた脚 D=2)について、
開いた脚 6 本を `fermion::slice(·, ax, x, x+1)` で各ラベル `x ∈ {0,1}` に固定した
64 通りの `pair(x)` を作り、`Σ_x contract_pair_MF(pair(x))` と `Σ_x contract_pair_MF(pair(x), hop)`
が、それぞれ一括の `contract_pair_MF(pair)` / `contract_pair_MF(pair, hop)` と一致すること
(`slice` は脚のパリティ区間を保つので、奇ラベルに固定した 1 次元脚は奇のまま)。
さらに `build_pair_state(A, B)` を `A` と `B` を個別に slice してから組んだものと、
組んでから slice したものが一致すること。

複素: 層1 のテンソルを `complex_tensor` に詰め替えた(虚部 0)`ftensor<complex_tensor>` でも
`contract_pair_MF` がコンパイルでき、実数版と一致し、虚部が 1e-12 以下であること。

**書いてはいけないこと.** テストの中で `apply_swap` を使って期待値を「手で」作らない
(陰性対照で `r4_wrap` を使うのは可)。`build_reduced_pair` の内部には触れない。

### 実装(Codex に渡す。テストファイル `test/fermion/mf_measure.cpp`、`test/test_fermion_layer.cpp` は変更禁止)

- [ ] **Step 1(Claude):** 契約書からテストを書かせ、RED を確認する。API 未定義のコンパイル
  エラーで赤になるのを避けるため、Claude が先に次のスタブを置く(コミットしない):
  `reduced.hpp` に `build_pair_state` / `apply_pair_op`、`reduced_measure.hpp` に
  `contract_pair_MF` ×2 を、いずれも `throw std::runtime_error("not implemented")` で宣言・定義。
  層1〜3 が「例外で落ちる」こと、R3/R5 が緑のままであることを確認。テストファイルの
  `sha256` を `work/mf-task2-test-snapshot.txt` に記録する

- [ ] **Step 2(Codex):** `src/fermion/reduced.hpp` — `build_reduced_pair` の直前に追加し、
  `detail::apply_reduced_two_site_op` を削除して `apply_pair_op` に置き換える:

```cpp
template <class tensor>
ftensor<tensor> build_pair_state(const ftensor<tensor>& TnA,
                                 const ftensor<tensor>& TnB,
                                 reduced_pair_direction direction) {
  if (TnA.rank() != 5 || TnB.rank() != 5) {
    throw std::runtime_error("build_pair_state expects five-leg Tn ftensors");
  }
  switch (direction) {
    case reduced_pair_direction::horizontal:
      // (l_A, t_A, b_A, s_A, t_B, r_B, b_B, s_B)
      return tensordot(TnA, TnB, mptensor::Axes(2), mptensor::Axes(0));
    case reduced_pair_direction::vertical:
      // (l_A, t_A, r_A, s_A, l_B, r_B, b_B, s_B)
      return tensordot(TnA, TnB, mptensor::Axes(3), mptensor::Axes(1));
  }
  throw std::runtime_error("build_pair_state: invalid direction");
}

// Apply a two-site operator op12 (in_A, in_B, out_A, out_B) to the physical
// legs (3, 7) of a pair state and restore the original leg order.
template <class tensor>
ftensor<tensor> apply_pair_op(const ftensor<tensor>& pair,
                              const ftensor<tensor>& op12) {
  if (pair.rank() != 8) {
    throw std::runtime_error("apply_pair_op expects an eight-leg pair state");
  }
  if (op12.rank() != 4) {
    throw std::runtime_error("apply_pair_op expects a four-leg operator");
  }
  ftensor<tensor> applied =
      tensordot(pair, op12, mptensor::Axes(3, 7), mptensor::Axes(0, 1));
  return transpose(applied, mptensor::Axes(0, 1, 2, 6, 3, 4, 5, 7));
}
```

  `build_reduced_pair` は次のように書き換える(`leg_ids` と以降のゲージ処理は不変):

```cpp
template <class tensor>
tensor build_reduced_pair(const ftensor<tensor>& TnA,
                          const ftensor<tensor>& TnB,
                          const ftensor<tensor>& op12,
                          reduced_pair_direction direction) {
  const ftensor<tensor> ket_ab = build_pair_state(TnA, TnB, direction);
  std::vector<int> leg_ids;
  switch (direction) {
    case reduced_pair_direction::horizontal:
      leg_ids = {0, 1, 3, 1, 2, 3};
      break;
    case reduced_pair_direction::vertical:
      leg_ids = {0, 1, 2, 0, 2, 3};
      break;
    default:
      throw std::runtime_error("doubled_cluster: invalid direction");
  }

  ftensor<tensor> ket_op = apply_pair_op(ket_ab, op12);
  ftensor<tensor> doubled =
      tensordot(conj(ket_ab), ket_op, mptensor::Axes(), mptensor::Axes());
  tensor ret = detail::fuse_doubled_cluster(doubled, leg_ids);
  // Gauge alignment with the single-site doubling convention; measured
  // directly via comparison against build_reduced_op/build_reduced-based
  // direct composition, not derived analytically.
  if (direction == reduced_pair_direction::horizontal) {
    detail::apply_fused_leg_gauge(ret, TnA.parity[3], 2, true);
    detail::apply_fused_leg_gauge(ret, TnB.parity[3], 5, false);
  } else {
    detail::apply_fused_leg_gauge(ret, TnA.parity[0], 0, true);
    detail::apply_fused_leg_gauge(ret, TnB.parity[0], 3, false);
  }
  return ret;
}
```

- [ ] **Step 3(Codex):** `src/fermion/reduced_measure.hpp` — `build_reduced_identity_pair` の
  直後に追加:

```cpp
namespace detail {

inline mptensor::Axes all_axes(int rank) {
  mptensor::Axes axes;
  for (int ax = 0; ax < rank; ++ax) {
    axes.push(ax);
  }
  return axes;
}

}  // namespace detail

// Mean-field norm of a two-site pair state: <pair|pair> as a graded full
// contraction. No environment tensors: the lambda weights on the open legs
// are expected to be multiplied into TnA / TnB beforehand (the same dressing
// the bosonic mean-field path applies in measure_twosite).
template <class tensor>
typename tensor::value_type contract_pair_MF(const ftensor<tensor>& pair) {
  const mptensor::Axes axes = detail::all_axes(pair.rank());
  return trace(conj(pair), pair, axes, axes);
}

// Mean-field expectation value (unnormalized) of op12 on a pair state.
// op12 must be loaded with wrap_twosite_gate (input-leg swap only): this is
// the single-layer convention pinned by the Fock-verified direct path
// (r2_expect_two in the R5 test), not the blob convention
// (wrap_reduced_pair_op).
template <class tensor>
typename tensor::value_type contract_pair_MF(const ftensor<tensor>& pair,
                                             const ftensor<tensor>& op12) {
  const mptensor::Axes axes = detail::all_axes(pair.rank());
  return trace(conj(pair), apply_pair_op(pair, op12), axes, axes);
}
```

  `reduced_measure.hpp` に `detail` 名前空間が既にある(`trace_boundary_pairs`)ので、
  `all_axes` はそこにまとめてよい。

- [ ] **Step 4(Codex):** `cmake --build out-gcc/build -j 8 --target test_fermion_layer` →
  `out-gcc/build/test/test_fermion_layer` 全件緑。報告を `work/mf-task2-report.md` に書く
  (変更ファイル、テスト結果の生出力の要約、疑問点)

- [ ] **Step 5(Claude):** `sha256sum` でテスト不改変を確認、`test_fermion_layer` を自分で実行、
  R3/R5(blob 経路)が緑のままであること、層2 の陰性対照が実際に「落ちる側」で評価されて
  いること(`CHECK(std::abs(...) > ...)` が真)を出力で確認。`clang-format -i` を
  `src/fermion/reduced.hpp src/fermion/reduced_measure.hpp test/fermion/mf_measure.cpp`
  に適用して差分を確認。レビューサブエージェントに変異テスト
  (`contract_pair_MF` の `conj` → 素の複素共役、`wrap_twosite_gate` → `wrap_reduced_pair_op`、
  graded `tensordot` → `mptensor::tensordot`)を依頼。コミット
  `"Add the single-layer mean-field pair contraction for fermions"`

---

## Task 3: `measure_twosite` の配線とガード解除

**Files:**
- Modify: `src/iTPS/twosite_obs.cpp:167-290`
- Modify: `src/iTPS/load_toml.cpp:642-644`
- Modify: `test/input.cpp:512-545`(サブケースの反転)
- Modify: `test/test_fermion_layer.cpp`(層4 の `TEST_CASE` 追加)

**Interfaces:**
- Consumes: T2 の `build_pair_state` / `contract_pair_MF`、`fops.hpp` の `wrap_twosite_gate` /
  `transpose`、`fermion_info.hpp` の `wrap_Tn`
- Produces: `MeanField_Env = true` かつ `finfo.enabled` のときの `measure_twosite()` /
  `measure_onesite()` の値。T4 の E2E が使う

### 契約書(テスト作成者に渡す散文)

**対象.** `tenes::itps::iTPS<tensor>` の `measure_twosite()` と `measure_onesite()` を、
`peps_parameters.MeanField_Env = true` かつ `peps_parameters.fermion = true` で呼んだときの
振る舞い。および `validate_fermion_constraints`(`src/iTPS/load_toml.cpp`)が
`MeanField_Env = true` を受理すること。

**テストの作り方.** `test/test_fermion_layer.cpp` 末尾の
「fermion twosite measurement is translation invariant across wraps」と同じ流儀:
`SquareLattice lattice(2, 2)`、`PEPS_Parameters params` に `fermion = true`、
`phys_parity`、`MeanField_Env = true`、`print_level = none`、`outdir` を設定し、
`iTPS<tensor> state(MPI_COMM_WORLD, params, lattice, updates, {}, onesite_ops, twosite_ops, {}, {}, {})`
を作る。`iTPSTestAccessor::Tn(state)` / `lambda_tensor(state)` / `finfo(state)` で
テンソル・λ・パリティ台帳を直接差し替える。`updates` は空でよい(MF 測定は環境更新を
必要としない)。**`update_reduced_density_environment` は呼ばない**(CTM を使わないことの
確認でもある)。`measure_twosite()` の戻り値は `vector<map<Bond, value>>` で、
`ret[group]` が演算子グループ、`ret.back()` が norm の map(`Bond{source, dx, dy}` がキー)。
`measure_onesite()` は `ret[group][site]`。

**テストすべきこと(層4).** 許容は 1e-12(相対)、0 近傍は絶対 1e-12。

1. **全偶パリティでボソン MF と一致(配線の回帰).** `phys_parity` を全サイト `{false, false}`、
   `finfo.virt` を全サイト・全脚 `{false, false}` にした fermion 状態と、同じ Tn・λ の
   `fermion = false` 状態(`MeanField_Env = true`)で、Tn を決定論的な乱数(`std::mt19937`、
   全要素非零でよい)、λ を脚ごとに異なる値(例 `{1.0, 0.6}`、`{0.9, 0.4}` …)にし、
   - 1サイト演算子(`{{0.2, 0.5}, {0.5, -0.3}}` のような対称行列)の `measure_onesite` が
     全サイトで一致、
   - 2サイト演算子(rank-4、要素を乱数で埋めた一般の実行列)を 8 本のボンド
     (`dx=+1`, `dy=+1` を各サイトから)および `dx=-1`, `dy=-1` を各サイトから、計 16 本
     で測った `measure_twosite` の値と norm が一致すること。
   `real_tensor` と `complex_tensor`(虚部 0)の両方で行い、両者も一致すること。
2. **norm はボソン MF と一致(§3.2).** `phys_parity = {false, true}`、`finfo.virt` は既定の
   `even_first_parity`(偶奇混在)のまま、偶パリティ要素だけ非零の決定論的 Tn で、
   fermion MF の norm(`ret.back()`)が `fermion = false` の同じ Tn・λ の norm と一致する。
3. **両端測定の一致.** 2 と同じ状態で、ホッピング
   `h(0,1,1,0) = h(1,0,0,1) = -1`(添字は `(in1, in2, out1, out2)`)を
   `Bond{site, +1, 0}` と `Bond{lattice.other(site, 1, 0), -1, 0}`(同じボンドを反対側から)で
   測った値が一致する。垂直(`dy=±1`)も同様。さらに d=4(`phys_parity = {false, true, true, false}`、
   `physical_dims = 4`)で R5 の `r4_hop_plain()` を演算子にして同じことを確認する
   (d=2 では graded transpose と素の transpose の差が出ないため)。
4. **並進不変.** 2×2 セルの全サイトに同じ Tn・同じ λ を入れると、水平 4 本の値が互いに
   一致し、垂直 4 本も互いに一致する(MF は CTM 収束を含まないので厳密に一致する)。
5. **1サイト密度がオラクルと一致(§3.4).** 1×1 の窓: T1 の `mf_single_site` /
   `mf_seed173_single_site`(`lam = [1.0, 0.7]`)と同じ `deterministic_tensor` を
   全サイトの Tn に、`lam` を全脚の λ に入れ、`measure_onesite` の密度(演算子 `n = diag(0,1)`)
   が参照値 `n0` と一致する(どのサイトでも同じ値)。
6. **`test/input.cpp`:** `SUBCASE("fermion rejects mean-field environment")` を
   `SUBCASE("fermion accepts mean-field environment")` に置き換え、同じ入力で
   `CHECK_NOTHROW(validate_fermion_constraints(...))` にする。他のガード(full update、
   `Use_RSVD` 等)のサブケースは触らない。
7. **2サイト値がオラクルと一致(配線込み).** 2×2 セルの Tn を T1 の `mf_horizontal_2site`
   の 2 テンソル(サイト 0 に A、サイト 1 に B、サイト 2・3 にも同じく A・B)、λ を
   `{1.0, 0.7}` にしたとき、`Bond{0, +1, 0}` のホッピング(`r2_hop_op` と同じ要素)が
   参照値 `hop01` と一致する。これは λ の dressing が「窓の外側の脚だけ」に掛かっていることの
   確認でもある(内部ボンドに二重に掛かると値が変わる)。
   注: この配線では `lattice.other(0, 1, 0) = 1` で、サイト 1 の左脚がサイト 0 の右脚に
   繋がる。サイト 1 の右脚はサイト 0 の左脚に巻き戻るが、窓の外側なので λ で閉じる。
   垂直も同様に `mf_vertical_2site`(オラクルのサイト 0 が上、サイト 1 が下: `Patch` の
   `internal_bonds` は `(s, "b", site(x, y+1), "t")`)で、A(サイト 0 のテンソル)を
   `lattice.other(s, 0, 1)`(上)に、B をサイト `s`(下)に置き、`Bond{s, 0, +1}`(source が下)と
   `Bond{lattice.other(s, 0, 1), 0, -1}`(source が上)の両方が `hop01` と一致すること。
   水平も `Bond{1, -1, 0}` を加える。`pair01` も同じ配置で確認する。

**書いてはいけないこと.** 期待値を `contract_pair_MF` で作って比べるだけのテスト
(実装の写しになる)。参照はオラクル値・ボソン経路・対称性から取る。

### 実装(Codex に渡す。`test/` 配下は変更禁止)

- [ ] **Step 1(Claude):** テスト作成者にテストを書かせ、RED/GREEN を確認する。
  実装前は fermion + MF がボソン MF 式に落ちるので、**RED になるのは 6(現ガードで例外)と
  7(ホッピングの符号・値が違う)だけ** である。1(全偶)・2(norm)・4(並進)・5(1サイト)は
  ボソン式でも通る回帰・不変性テストであり、実装前から GREEN でよい。3(両端)も d=2 では
  ボソン式で通りうる。7 が RED であることは必ず確認する(値と参照値を記録)。
  `work/mf-task3-test-snapshot.txt` に sha256 を記録

- [ ] **Step 2(Codex):** `src/iTPS/load_toml.cpp` の

```cpp
  if (peps_parameters.MeanField_Env) {
    throw_fermion_guard("MeanField_Env=true");
  }
```

  を削除する。

- [ ] **Step 3(Codex):** `src/iTPS/twosite_obs.cpp` の norm 計算ブロック(`if (norms.count(norm_key) == 0) {` から対応する `}` まで)を次に置き換える:

```cpp
    const auto norm_key = Bond{indices[nrow - 1][0], nrow - 1, ncol - 1};
    if (norms.count(norm_key) == 0) {
      if (finfo.enabled && !is_TPO && nrow * ncol == 2) {
        if (nrow == 2) {
          const int top = indices[0][0];
          const int bottom = indices[1][0];
          const auto fTop = tenes::fermion::wrap_Tn(*(Tn_[0][0]), finfo, top);
          const auto fBottom =
              tenes::fermion::wrap_Tn(*(Tn_[1][0]), finfo, bottom);
          if (is_mf) {
            // Mean field: Tn_ are the lambda-dressed boundary copies, so the
            // single-layer graded contraction already closes the window.
            norms[norm_key] = tenes::fermion::contract_pair_MF(
                tenes::fermion::build_pair_state(
                    fTop, fBottom,
                    tenes::fermion::reduced_pair_direction::vertical));
          } else {
            const ptensor blob = tenes::fermion::build_reduced_identity_pair(
                fTop, fBottom,
                tenes::fermion::reduced_pair_direction::vertical);
            norms[norm_key] =
                tenes::fermion::contract_reduced_pair_vertical_density_CTM(
                    C1[top], C2[top], C3[bottom], C4[bottom], eTt[top],
                    eTr[top], eTr[bottom], eTb[bottom], eTl[bottom], eTl[top],
                    blob);
          }
        } else {
          const int left = indices[0][0];
          const int right = indices[0][1];
          const auto fLeft = tenes::fermion::wrap_Tn(*(Tn_[0][0]), finfo, left);
          const auto fRight =
              tenes::fermion::wrap_Tn(*(Tn_[0][1]), finfo, right);
          if (is_mf) {
            norms[norm_key] = tenes::fermion::contract_pair_MF(
                tenes::fermion::build_pair_state(
                    fLeft, fRight,
                    tenes::fermion::reduced_pair_direction::horizontal));
          } else {
            const ptensor blob = tenes::fermion::build_reduced_identity_pair(
                fLeft, fRight,
                tenes::fermion::reduced_pair_direction::horizontal);
            norms[norm_key] =
                tenes::fermion::contract_reduced_pair_horizontal_density_CTM(
                    C1[left], C2[right], C3[right], C4[left], eTt[left],
                    eTt[right], eTr[right], eTb[right], eTb[left], eTl[left],
                    blob);
          }
        }
      } else if (is_mf) {
        norms[norm_key] = core::Contract_iTPS_MF(Tn_, op_);
      } else if (is_TPO) {
        norms[norm_key] =
            core::Contract_density_CTM(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_);
      } else {
        norms[norm_key] =
            core::Contract_iTPS_CTM(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_);
      }
    }
    auto norm = norms[norm_key];
```

- [ ] **Step 4(Codex):** 同ファイルの値計算、`if (nrow == 2) {` の垂直分岐の
  `if (finfo.enabled && !is_TPO && !is_mf) { ... }` を次に置き換える(`else` のボソン分岐と
  コメントアウトされた旧コードは現状のまま):

```cpp
          if (finfo.enabled && !is_TPO) {
            const int target = top == source ? bottom : top;
            // Two loading conventions, both pinned by the d = 4 R5 oracle:
            //  - blob path: wrap_reduced_pair_op (input and output swaps),
            //  - single-layer mean-field path: wrap_twosite_gate (input swap
            //    only), the convention of the Fock-verified direct path.
            auto o = is_mf ? tenes::fermion::wrap_twosite_gate(
                                 op.op, finfo.phys[source], finfo.phys[target])
                           : tenes::fermion::wrap_reduced_pair_op(
                                 op.op, finfo.phys[source], finfo.phys[target]);
            if (top != source) {
              // Graded transpose: carries the Fock reordering sign
              // |n_B n_A> = (-1)^{n_A n_B} |n_A n_B> on both leg pairs.
              o = tenes::fermion::transpose(o, mptensor::Axes(1, 0, 3, 2));
            }
            const auto fTop = tenes::fermion::wrap_Tn(*(Tn_[0][0]), finfo, top);
            const auto fBottom =
                tenes::fermion::wrap_Tn(*(Tn_[1][0]), finfo, bottom);
            if (is_mf) {
              value = tenes::fermion::contract_pair_MF(
                  tenes::fermion::build_pair_state(
                      fTop, fBottom,
                      tenes::fermion::reduced_pair_direction::vertical),
                  o);
            } else {
              // The window rows already carry geometric roles (row 0 =
              // upper), and the infinite lattice has no boundary, so bonds
              // whose target wraps around the unit cell need no special
              // ordering.
              const ptensor blob = tenes::fermion::build_reduced_pair(
                  fTop, fBottom, o,
                  tenes::fermion::reduced_pair_direction::vertical);
              value =
                  tenes::fermion::contract_reduced_pair_vertical_density_CTM(
                      C1[top], C2[top], C3[bottom], C4[bottom], eTt[top],
                      eTr[top], eTr[bottom], eTb[bottom], eTl[bottom],
                      eTl[top], blob);
            }
          } else {
```

  水平分岐(`} else {  // ncol == 2`)も同型に置き換える(`left` / `right`、
  `left != source` で transpose、`reduced_pair_direction::horizontal`、
  `contract_reduced_pair_horizontal_density_CTM(C1[left], C2[right], C3[right], C4[left], eTt[left], eTt[right], eTr[right], eTb[right], eTb[left], eTl[left], blob)`)。

- [ ] **Step 5(Codex):** `cmake --build out-gcc/build -j 8` → `out-gcc/build/test/test_fermion_layer`、
  `out-gcc/build/test/test_input` 緑。`work/mf-task3-report.md` に報告

- [ ] **Step 6(Claude):** sha256 でテスト不改変確認。`ctest --test-dir out-gcc/build` **全件**
  (ボソン経路・CTM フェルミオン経路の回帰: `AntiferroHeisenberg_mf`、`FreeFermion`、
  `FreeFermionSimple` を含む)。`clang-format -i src/iTPS/twosite_obs.cpp src/iTPS/load_toml.cpp`
  (テストファイルは作成者に整形済みを要求)。レビューサブエージェントに、
  「`is_mf` のボソン分岐・`is_TPO` 分岐・CTM フェルミオン分岐の到達条件が変更前と同じか」
  を分岐表で確認させる。コミット `"Allow the mean-field environment in fermion mode"`

---

## Task 4: E2E テスト、ドキュメント、NEWS

**Files:**
- Create: `test/fermion/free_fermion_mf.py.in`
- Modify: `test/CMakeLists.txt`(`FreeFermion` の登録の直後)
- Create: `test/data/output_FreeFermionMF/free_fermion_mf.dat`(Claude が生成)
- Modify: `docs/sphinx/ja/file_specification/parameter_section.rst:63,164`
- Modify: `docs/sphinx/en/file_specification/parameter_section.rst:62,160`
- Modify: `NEWS.md`

### 契約書(E2E、テスト作成者に渡す散文)

`test/fermion/free_fermion.py.in` を雛形に `test/fermion/free_fermion_mf.py.in` を書く
(`@CMAKE_BINARY_DIR@` 等の置換子、`run_tenes`、`input_toml`、`read_density`、
`two_site_hamiltonian`、`exact_free_fermion` は同じものを使ってよい。共有モジュール化は
しない — `.py.in` は CMake が個別に configure する)。

- 入力は `free_fermion.py.in` の `input_toml` に `[parameter.ctm] meanfield_env = <true|false>`
  を加えたもの。`seed = 11` 固定。`μ = 0`(半充填)のみ。
- 走らせる組: (D=2, CTM)、(D=2, MF)、(D=4, MF)。`num_step = 1000`、`tau = 0.01`。
  CTM 版と MF 版は simple update が同一なので **Tn は同じ**、測定だけが違う。
- 各実行について: 終了コード 0、`output/density.dat` に `Energy` と `n` がある、
  `correlation_length.dat` が作られない、MF 版では `output/time.dat` の
  `time environment` が `1e-6` 未満(CTM を作っていない証拠)。
- 物理の健全性: D=2 で `|E_MF − E_CTM| / |E_CTM| < 0.2`、MF の `n` が `0.5 ± 0.02`
  (D=2, 4 とも)、D=4 MF の `relerr(E_MF, E_exact) < 0.1`、MF のエネルギーが負。
  **これらの閾値は事前見積もりである**。実測が外れたら閾値を黙って緩めず、値を添えて
  Claude に報告する(MF の精度の問題か、符号の問題かを切り分ける)。
- golden: 結果を `output_FreeFermionMF/free_fermion_mf.dat` に
  `# mu D env energy n exact_energy exact_n` の形式で書き、
  `@CMAKE_SOURCE_DIR@/test/data/output_FreeFermionMF/free_fermion_mf.dat` と
  `np.isclose(rtol=1e-3, atol=1e-4)` で比較する。golden が無ければ
  「golden がありません」と印字して失敗する(初回は Claude が結果を確認してから
  golden を置く。E2E はそれまで RED)。
- `test/CMakeLists.txt`: `configure_file(... free_fermion_mf.py.in ... free_fermion_mf.py @ONLY)`、
  `add_test(NAME FreeFermionMF COMMAND ${TENES_PYTHON_EXECUTABLE} .../fermion/free_fermion_mf.py)`、
  `set_tests_properties(FreeFermionMF PROPERTIES TIMEOUT 1800)`。`FreeFermion` の直後に置く。

### ドキュメント・NEWS(Codex に渡す。`test/` 配下は変更禁止)

- [ ] **Step 1(Claude):** E2E をテスト作成者に書かせ、`ctest --test-dir out-gcc/build -R FreeFermionMF`
  が「golden がありません」で RED になることを確認。sha256 を記録
- [ ] **Step 2(Claude):** E2E を 1 回走らせて `output_FreeFermionMF/free_fermion_mf.dat` を得て、
  値(特に D=2 の MF vs CTM、D=4 の厳密解との差)を確認してから
  `test/data/output_FreeFermionMF/free_fermion_mf.dat` にコピーし、再実行で GREEN を確認。
  閾値違反があれば Task 3 に差し戻す(値を添えて)
- [ ] **Step 3(Codex):** `docs/sphinx/ja/file_specification/parameter_section.rst` の fermion の
  非対応一覧(63 行目付近)から「平均場環境、」を外し、`meanfield_env` の行(164 行目付近)の
  説明を
  「CTM ではなく simple update で得られる平均場環境を用いる。フェルミオン模式でも使用でき、
  2サイト観測量は単層縮約で評価されるため CTM 版より大幅に軽いが、精度は simple update 相当」
  にする。英語版(62 行目、160 行目付近)も同内容で:
  "Use mean field environment obtained through simple update instead of CTM. Also available
  in fermion mode, where two-site observables are evaluated by a single-layer contraction
  (much cheaper than the CTM path, with simple-update-level accuracy)"。
  非対応一覧からは "mean-field environment, " を外す
- [ ] **Step 4(Codex):** `NEWS.md` の fermion の項(9 行目付近の Limitations)から
  "mean-field environment, " を外し、その直後に項目を足す:
  "`meanfield_env = true` is now accepted in fermion mode; two-site observables are then
  evaluated by a single-layer graded contraction (cost D^6 d^2 instead of the D^12 d^4
  reduced-pair blob of the CTM path), which makes D > 2 and the Hubbard model (d = 4)
  measurable at simple-update accuracy"。報告を `work/mf-task4-report.md` に
- [ ] **Step 5(Claude):** `black test/fermion/free_fermion_mf.py.in` は `.in` なので掛けない
  (置換子が構文エラーになる)。目視で 88 桁とスタイルを確認。docs の RST を
  `out-docs` プリセットがあればビルドして警告が増えていないことを確認。
  `ctest --test-dir out-gcc/build` 全件緑を確認してコミット
  `"Add the mean-field free-fermion end-to-end test and document meanfield_env for fermions"`

---

## 最終レビュー(全タスク後)

- [ ] 最上位モデルのレビューサブエージェントに、ブランチ全体(T1〜T4 の 4 コミット)を
  1 回レビューさせる。観点: (1) どのテストも通らない経路(例: `is_TPO && finfo.enabled`)、
  (2) `wrap_twosite_gate` / `wrap_reduced_pair_op` の使い分けが値・norm の両方で正しいか、
  (3) λ の dressing が窓の外側だけか、(4) 変異テスト 3 件(Task 2 Step 5)の結果、
  (5) ドキュメント文言に `M1` が無いか
- [ ] `clang-format` / `black` の差分が変更ファイル以外に及んでいないこと
- [ ] `work/mf-task*-report.md` と `work/mf-task*-test-snapshot.txt` は PR に含めない
  (`work/` は実験用、`tenes-work-dir-convention`)
