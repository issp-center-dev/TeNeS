# fermion 長距離ハミルトニアン 実装計画

> **For agentic workers:** この計画は CLAUDE.md の多段エージェント手順で実行する。
> 振る舞い契約書(散文)→ テスト作成者(Claude サブエージェント)→ RED 確認とスナップショット →
> Codex 実装 → テスト不改変の機械検査 → Claude の独立検証 → タスクレビュー → 整形とコミット。
> ステップは `- [ ]` で追跡する。**テストコードはこの計画に書かない**(テスト作成者に計画者のテストを見せないため)。

**Goal:** fermion モードで、最近接より遠いボンドを含むハミルトニアン(正方格子の t'・t''・V'、三角・ハニカム・カゴメ格子)の基底状態を、シンプル更新で計算できるようにする。

**Architecture:** `tenes_std` が長距離ボンドのゲートを、経路順の graded MPO(途中サイトに奇チャネルのパリティ紐、パリティブロックごとの逐次 SVD)として NN ゲート列に分解する。
solver は、物理次元が一時的に変わるゲート列を入力時に模擬実行して各脚のパリティ台帳を推定し、SU の実行時に `finfo.phys` を持ち回す。
SU のボンド処理そのもの(graded 縮約・raster 化)は変えない。

**Tech Stack:** C++17, mptensor, doctest(単体)、Python 3.9+ / pytest / numpy(ツールと E2E)、ctest。

**Spec:** `docs/superpowers/specs/2026-09-30-fermion-longrange-hamiltonian-design.md`(以下「設計書」)。実行者は設計書とこの計画の両方を読むこと。

## Global Constraints

- ブランチ `fermion-longrange-ham`。実験・台帳・ディスパッチ文面・報告は `work/fermion-longrange-ham/<task>/` に置き、CWD をそこにして実行する。リポジトリ直下やビルドツリーで `tenes` やテストバイナリを実行しない。
- C++ はすべて namespace `tenes`(fermion 層は `tenes::fermion`)。新規ファイルの先頭に GPL v3 ヘッダを付ける。
- Python の最低バージョンは 3.9(CI)。3.10 以降の API(`int.bit_count()`、`match` など)を使わない。
- 二サイト演算子の配列規約は `op[in1, in2, out1, out2]`(1 = source、2 = target)。行列の意味は source 先の順序付き Fock 基底 |n_s n_t⟩ = (c†_s)^{n_s}(c†_t)^{n_t}|0⟩(設計書 §3)。
- ゲート列の規約:物理次元が変わるゲートでは、太る脚は target 側の出力 `out2`、消費される脚は source 側の入力 `in1`(設計書 §6.1)。
- 経路順の符号:`T[...] = evo[i_s, i_t, o_s, o_t] · Π_m δ(i_m, o_m) · (−1)^{(Σ_m p(i_m)) · (p(i_t) + p(o_t))}`(設計書 §4.1)。
- 1 ホップのゲートと、ボゾンの分解経路は変えない。ボゾンの `tenes_std` 出力はバイト一致を保つ。
- FU では長距離ゲートを使わない。solver は「`out2` の台帳が target の元の台帳と異なる、または入力脚の次元が元の物理次元と異なるゲート」が `evolution.full` にあり、FU の `num_step` のどれかが正なら拒否する(設計書 §6.1)。
- MPI:推定と検査は全ランクで集団的に行い、拒否は全ランクが同じ判定で投げる。rank 0 だけが投げるとハングする。
- fermion モードは「未対応の入力は読み込み時に拒否する」方針を守る。拒否メッセージは `throw_fermion_guard` の形式に合わせる。
- 整形:C++ は `git add` のあと `git clang-format`(変更行のみ)。ファイル全体に `clang-format -i` を当てない。Python は `black`(line-length 88)。整形は Claude がコミット直前に行う。Codex には実行させない。
- ビルドとテスト:非 MPI は `cmake --preset gcc && cmake --build --preset gcc && ctest --preset gcc`(`OMP_NUM_THREADS=1`)。MPI は `out-gcc-mpi/build`。
- golden ファイルは意図的に再生成する。手で編集しない。

## Review Focus

仕様が暗に要求しているのに、タスク内のテストが自然には踏まない入力。各行のテストは括弧内のタスクに追加してある。

1. **経路が左か上に進む段**:TOML の `source_leg` が 0 か 1 のゲート。raster 化で太る脚が s1 側に移る。台帳の更新を TOML 上の source/target で行う誤りは、右と下にしか進まない経路では見えない(T2)。
2. **経路上で同じ単位胞サイトが再び現れる**:2×2 単位胞の (2, 0) では target が source と同じサイトになる。3 ホップでは途中サイトが source や target と同じサイトになりうる。台帳の持ち回しとテンソルの上書き順が絡む(T2)。
3. **χ の次元が d と等しいゲート**:台帳を「次元が同じなら元の台帳」で決める誤りが通ってしまう。χ = d_m × r なので、`tenes_std` の出力では O_st が積演算子のときだけ起きる(T1 は積演算子で作って確かめ、T2 は同じゲートか手書きのゲートで受理を確かめる)。
4. **identity ゲートの補完**:`complete_ungated_bonds()` が末尾に足す identity ゲートと、ゲート列の共存。経路が単位胞の全ボンドを覆わない入力(T2)。
5. **complex テンソル**:実時間ではなくても、`is_real = false` の入力では complex で走る。conj の抜けは実数テストでは見えない(T2)。

---

## ファイル構成

| ファイル | 役割 | タスク |
|---|---|---|
| `tool/tenes_std.py` | fermion 分解(経路順の符号と、パリティブロックごとの逐次 SVD)、長距離の拒否の撤去、FU の拒否 | T1 |
| `src/operator.hpp` | `EvolutionOperator` に 4 脚のパリティ台帳の欄を追加 | T2 |
| `src/fermion/fops.hpp` | 入出力で別の台帳を取る `wrap_twosite_gate` | T2 |
| `src/iTPS/load_toml.cpp` / `load_toml.hpp` | ゲート列の模擬実行(台帳の推定と検査)、旧ガードの置換、FU の拒否 | T2 |
| `src/iTPS/main.cpp` | 推定の呼び出し(`complete_ungated_bonds()` の後) | T2 |
| `src/iTPS/simple_update.cpp` | 保持した台帳でゲートを包む、raster 化後の台帳で `finfo.phys` を更新、群末の不変条件 | T2 |
| `tool/tenes_simple.py` | `_check_fermion_scope` の緩和(T3:正方の 2・3 近接と三角、T5:ハニカムとカゴメ) | T3, T5 |
| `src/fermion/`・`src/iTPS/` の該当箇所 | D = 1 の脚と空孔で見つかった不具合の修正(あれば) | T4 |
| `test/python/test_tenes_std.py` ほか | ツールのテスト(新規 `test/python/test_fermion_chain.py` を含む) | T1, T3, T5 |
| `test/fermion/longrange_gate.cpp`(新規) | T2 の単体テスト(実行ファイル `test_fermion_longrange_gate`) | T2 |
| `test/fermion/unit_d1.cpp`(新規) | T4 の単体テスト(実行ファイル `test_fermion_d1`) | T4 |
| `test/fermion/*.py.in`(新規)、`test/data/output_*`(新規) | E2E と golden | T3, T5 |
| `test/CMakeLists.txt` | 新しい実行ファイルと E2E の登録 | T2–T5 |
| 既存テスト(`test/python/test_fermion_models.py`、`test/python/test_tenes_std.py`、`test/fermion/fermion_guards.cpp` など) | 許可するものの拒否を確かめている箇所の更新(テスト作成者が行う) | T1–T3, T5 |
| `docs/sphinx/{ja,en}/...`、`NEWS.md` | ドキュメント | T3, T5 |

## 共通手順(全タスク)

各タスクの「手順」欄は、この共通手順を具体化したものである。ディスパッチ文面には以下の定型注意を**最初のタスクから**入れる。

**テスト作成者(Claude サブエージェント、`general-purpose`)への定型注意**
- 渡すのはそのタスクの「振る舞い契約書」節と設計書だけ。計画者のテストコードは存在しない。
- 変更してよいのは、契約書が指定したテストファイル、`test/CMakeLists.txt`、契約書が更新を指示した既存テストだけ。`src/` と `tool/` は変更禁止。
- `git stash` / `git checkout` / `git reset` を含め、作業ツリーを巻き戻す操作は禁止。仮説の検証はファイルを scratchpad にコピーして行う。
- `tenes`・E2E スクリプト・テストバイナリの同時実行は最大 2 本。各 `OMP_NUM_THREADS=1`、ctest は `-j 2` 以下(一斉にバックグラウンドで投げると、ユーザーのマシンが過負荷になる)。
- 契約書の誤り(不可能な自己検査、空洞な参照値、符号の誤りなど)に気づいたら、テストを曲げずに報告する。
- 参照値は、テスト対象の機構を使わずに得ること(Fock 空間での生成消滅演算子による構成、Fock oracle、解析解、検証済みの既存機構)。**`tenes_std` の符号式を参照値の構成に使ってはならない。**
- 報告は `work/fermion-longrange-ham/<task>/test-author-report.md` に書く。存在しないこと自体が欠陥。

**Codex(実装)への定型注意**
- `codex-companion.mjs task --background --fresh --write --model gpt-5.5 "$(cat <prompt file>)"` で起動する。プロンプトは必ずファイル経由で渡す。
- 「確認を求めずに実装完了まで一気に進めること」と明記する。
- テストファイル(スナップショットを取ったもの)は変更禁止。テストが誤りと思ったら実装を曲げず、BLOCKED として報告する。
- formatter は実行しない。新しい行は周囲のスタイルに手で合わせる。
- sandbox はコミットできない。コミットは Claude が代わりに行う。
- `tenes` やテストバイナリをリポジトリ直下で実行しない。同時実行は最大 2 本、各 `OMP_NUM_THREADS=1`、ctest は `-j 2` 以下。
- 報告は `work/fermion-longrange-ham/<task>/codex-report.md` に書く。存在しないこと自体が欠陥。
- 完了後、Claude は `status --json` の running と pid の生存を確認する。`git status` で変更がゼロなら承認待ちを疑い、`--resume` で「はい」を返す。

**RED 確認とスナップショット(Claude)**
1. テストが実行でき(C++ はコンパイルでき)、実装前は**正しい理由で**失敗することを 1 件ずつ確認する。偶然通るテスト(空洞)がないか見る。
2. `test/` 全体のスナップショット(`find test -type f | sort | xargs shasum -a 256`)を `work/fermion-longrange-ham/<task>/test-snapshot.sha256` に取る。列挙型のスナップショットは、既存テストの書き換えをすり抜ける。
3. Codex の完了後、スナップショットと照合し、`git diff --stat -- test/` も見る。不一致なら差し戻す。テストの欠陥はテスト作成者に差し戻す。

**独立検証(Claude)**
- Codex の報告を鵜呑みにしない。該当テストを自分で回す。共有コード(`src/`、`tool/`)に触れたタスクでは、最後に `ctest` 全件を自分で回す。MPI ビルドも回す(n = 1 と、登録されていれば n = 2)。
- 変異テストを自分で 1〜2 件試す(契約書の変異リストから選び、scratchpad のコピーで)。

**タスクレビュー**:タスクごとに新しいレビュアー(`feature-dev:code-reviewer`)。作業ツリーの変更禁止を明記し、変異テストを推奨する。

**整形とコミット**:`git add` → `git clang-format` → `black` → テスト → `git diff --cached --stat` で意図しないファイルがないか確認 → コミット。コミットメッセージの末尾に帰属行を付ける。

**台帳**:各タスクの結果(RED の件数、スナップショット、独立検証の結果、変異の結果、レビュー指摘と処置)を `work/fermion-longrange-ham/LEDGER.md` に記録する。

---

### Task 1: `tenes_std` の fermion 分解

**Files:**
- Modify: `tool/tenes_std.py`(`make_evolution_twosite`、`make_evolution`、`Model.__init__`、`Model._validate_fermion_mode_input`)
- Create: `test/python/test_fermion_chain.py`(テスト作成者)
- Modify: `test/python/test_tenes_std.py`(長距離ボンドの拒否を確かめている既存テストの更新。テスト作成者)

**Interfaces:**
- Consumes(既存):`LatticeGraph.make_path(bond) -> List[Bond]`、`Unitcell.sites[i].parity`(0/1 のリスト、fermion モードでは必須)、`NNOperator`、`_drop_tiny`。
- Produces(T2・T3 が使う):

```python
def make_evolution_twosite(
    hamiltonian: NNOperator,
    graph: LatticeGraph,
    tau: Union[float, complex],
    group: int = 0,
    result_cutoff: float = 1e-15,
    fermion: bool = False,
) -> List[NNOperator]: ...

def make_evolution(
    hamiltonian: Operator,
    graph: LatticeGraph,
    tau: Union[float, complex],
    group: int = 0,
    result_cutoff: float = 1e-15,
    fermion: bool = False,
) -> Union[List[SiteOperator], List[NNOperator]]: ...
```

`fermion=True` かつ経路が 2 ホップ以上のとき、fermion 分解を使う。
返る `NNOperator` の並び、各要素の `bond`(経路の手前側が source)、`elements` の脚の並びは、ボゾンの分解と同じ(設計書 §2)。
`Model` は `parameter.general.fermion` が真のとき `fermion=True` を渡す。

#### 振る舞い契約書(T1)

1. **合成の厳密照合(主条件)**:`make_evolution_twosite(..., fermion=True)` が返すゲート列を、経路順 (s, m1, …, t) の Fock 空間で順に合成する。結果は、O_st ⊗ 1_mid を**生成消滅演算子から直接組んだもの**と、機械精度(絶対誤差 1e-12)で一致する。
   - 参照の O_st ⊗ 1_mid は、設計書 §3 の定義(source 先の順序付き基底での行列要素)と、経路順の Jordan–Wigner 表現から構成する。`tenes_std` の符号式を使わない。
   - ゲート列の合成では、各ゲートを経路順の基底での素の Kronecker 埋め込みとして扱ってよい(偶演算子なので隣接する位置への埋め込みに符号は出ない)。途中サイト m の太った脚は、ゲート列の途中ではサイト m の局所空間として扱う。この扱いの正当性をテストのコメントに一文で書く。
   - 演算子:スピンレス(d = 2)の hopping −t(c†_s c_t + h.c.) と n_s n_t の exp(−τH)。Hubbard(d = 4、`tool/tenes_simple.py` の `HubbardModel` と同じ局所基底 i = n_up + 2 n_dn、パリティ [0, 1, 1, 0])の hopping(両スピン)+ 近接の密度相互作用の exp(−τH)。τ は 0.1 程度で、奇チャネルが無視できない大きさにする。
   - 経路:設計書 §4.3 の一覧をすべて含める。(2, 0)、(0, 2)、(−2, 0)、(0, −2)、(1, 1)、(1, −1)、(−1, 1)、(−1, −1)、(2, 1)、(−2, 1)。単位胞は、経路のサイトが互いに異なる大きさのもの(例:4×4)と、同じ単位胞サイトが経路上に再び現れるもの(2×2 の (2, 0))の両方。
   - 経路は `graph.make_path` が選んだものをそのまま使う(テストが経路を指定するのではない)。各ケースで、実際の経路(ホップ数と向き)をテストの失敗メッセージに出す。
   - 各ケースが左か上に進む段を含むかどうかを記録し、少なくとも 1 ケースが `source_leg` 0 の段、1 ケースが `source_leg` 1 の段を含むことを、テスト自身が確かめる(含まなければテストが失敗する)。
2. **χ のパリティが一意**:各ゲートの `out2`(末尾ゲートを除く)の各添字について、非零要素から求めたパリティ p(in1) + p(in2) + p(out1) mod 2 が一つに決まる。非零要素を持たない添字がない。末尾ゲートは、全脚が既知の台帳で偶である。中間・末尾ゲートの `in1` のパリティ表は、直前のゲートの `out2` と一致する。
3. **χ = d の場合**:契約 1 のケースのうち、少なくとも 1 つのゲートで `out2` の次元が d と等しく、かつそのパリティ表が物理パリティ表と異なるものがあることを確かめる。見つからなければ、そうなる演算子(例えば偶チャネルだけ、あるいは特定の対称性を持つもの)を追加して作る。作れなかった場合は、その旨を報告する(テストを曲げない)。
4. **1 ホップとボゾンは変わらない**:
   - fermion モードの 1 ホップのゲートは、`fermion=False` の出力と要素ごとに一致する。
   - ボゾンの入力(fermion なし)の `Model(...).to_toml()` の出力は、変更前の実装とバイト一致する。変更前の出力は、テスト作成者が現在の `tool/tenes_std.py`(実装前)で生成し、テストデータとして保存する。長距離ボンドを含むボゾン入力を少なくとも 1 つ含める。
5. **入力検査**:
   - fermion モードで長距離ボンドを含み、`parameter.full_update.num_step > 0` なら `RuntimeError`。`num_step` が 0 または未指定なら受理する。
   - fermion モードの長距離ボンドは受理される(既存の拒否テストを、受理を確かめる形に更新する)。
   - 自己隣接セルの拒否、パリティ台帳の必須化、multisite の拒否は変わらない。
6. **変異に対する感度**(テスト作成者は変異を実行しなくてよいが、次の変異で必ず赤くなるケース選定にすること。レビュアーが確認する):
   パリティ紐を外す、紐を (Σ p(i_m))·p(i_t) のように片側だけにする、ブロック SVD を素の SVD に替える、特異値の打ち切りを外して全ゼロの添字を残す、`fermion` フラグを `Model` から渡し忘れる。
7. 実行時間:`pytest test/python/test_fermion_chain.py` が 60 秒以内。

#### 手順(T1)

- [ ] **Step 1: テスト作成者をディスパッチ**。振る舞い契約書(T1)節と設計書を渡す。変更前のボゾン出力の生成(契約 4)を、実装前のこの時点で行わせる。報告 `work/fermion-longrange-ham/t1/test-author-report.md`。
- [ ] **Step 2: RED 確認**。`python -m pytest test/python/test_fermion_chain.py test/python/test_tenes_std.py -q` を `work/fermion-longrange-ham/t1/` から実行する。契約 1〜3 は「fermion 引数がない」か「長距離ボンドを拒否する」理由で落ち、契約 4 は通ることを確認する。契約 1 の参照値の構成に `tenes_std` の符号式が使われていないことを読んで確かめる。スナップショットを取る。
- [ ] **Step 3: Codex に実装させる**。設計書 §4.1・§5.1 と Global Constraints を渡す。報告 `work/fermion-longrange-ham/t1/codex-report.md`。
- [ ] **Step 4: 不改変検査と独立検証**。スナップショット照合 → pytest(該当ファイル、次に `test/python` 全体)→ 変異 2 件(パリティ紐を外す、ブロック SVD を素の SVD に替える)を scratchpad のコピーで試し、赤くなることを確認する → `ctest --preset gcc -R python_unittest`。
- [ ] **Step 5: タスクレビュー → 整形 → コミット**。

```bash
git add tool/tenes_std.py test/python/test_fermion_chain.py test/python/test_tenes_std.py
black tool/tenes_std.py test/python/test_fermion_chain.py test/python/test_tenes_std.py
git commit -m "Decompose long-range fermion gates into a graded chain of nearest-neighbour gates"
```

---

### Task 2: solver のゲート列(台帳の推定と持ち回し)と符号規則のゲート

**Files:**
- Modify: `src/operator.hpp`(`EvolutionOperator`)
- Modify: `src/fermion/fops.hpp`(`wrap_twosite_gate`)
- Modify: `src/iTPS/load_toml.cpp`、`src/iTPS/load_toml.hpp`
- Modify: `src/iTPS/main.cpp`
- Modify: `src/iTPS/simple_update.cpp`
- Create: `test/fermion/longrange_gate.cpp`(テスト作成者)
- Modify: `test/CMakeLists.txt`(実行ファイル `test_fermion_longrange_gate` の登録。MPI ビルドでは n = 2 の登録も。テスト作成者)
- Modify: 既存テストのうち「奇の二サイトゲート」の拒否を確かめているもの(テスト作成者。拒否されることは保ち、メッセージ依存があれば新しい拒否理由に合わせる)

**Interfaces:**
- Consumes:T1 の `tenes_std` の出力(テスト入力の生成に使ってよい)。既存の `iTPSTestAccessor`(`test/test_fermion_common.hpp`:`Tn`、`lambda_tensor`、`finfo`)、`fock_oracle.py`、`build_relay_window`(`src/fermion/relay.hpp`)、`detail::doubled_pipeline`、`core::Contract_density_CTM`。
- Produces(Claude が実装前にスタブとして置く):

```cpp
// src/operator.hpp
template <class tensor>
struct EvolutionOperator {
  // ... existing members ...
  //! Fermion mode only: parity ledgers of the legs (in1, in2, out1, out2)
  //! of a two-site gate, or (in, out) of a one-site gate, in the orientation
  //! of the input (1 = source). Filled by infer_fermion_gate_ledgers();
  //! empty in bosonic runs.
  std::vector<std::vector<bool>> fermion_legs;
};

// src/fermion/fops.hpp
//! Load a two-site gate whose four legs may carry different ledgers.
//! The existing three-argument overload forwards here with (p1, p2, p1, p2).
template <class tensor>
ftensor<tensor> wrap_twosite_gate(const tensor& op, const parity_vector& in1,
                                  const parity_vector& in2,
                                  const parity_vector& out1,
                                  const parity_vector& out2);

// src/iTPS/load_toml.hpp
//! Fermion mode: simulate each group of simple_updates and full_updates in
//! list order, infer the ledger of every gate's out2 leg (design 6.1), store
//! the four ledgers in EvolutionOperator::fermion_legs, and reject inputs
//! that break the chain rules (input_error via throw_fermion_guard, on every
//! rank). No-op unless peps_parameters.fermion. Must be called on the list
//! returned by complete_ungated_bonds().
template <class tensor>
void infer_fermion_gate_ledgers(const PEPS_Parameters& peps_parameters,
                                const SquareLattice& lattice,
                                EvolutionOperators<tensor>& simple_updates,
                                EvolutionOperators<tensor>& full_updates);
```

`main.cpp` は `complete_ungated_bonds()` と `load_full_updates()` の後、`validate_fermion_constraints()` の前に `infer_fermion_gate_ledgers()` を呼ぶ。
`validate_fermion_constraints()` からは、ゲートのパリティ検査(`parity-odd one-site gates`、`parity-odd two-site gates` と full の同等品)を外す。役割は `infer_fermion_gate_ledgers()` が引き継ぐ。

#### 振る舞い契約書(T2)

目的:設計書 §4.3 のゲートを機械的に判定する。このタスクの主条件の合否が、そのまま方式の採否になる。

1. **SU のゲート列の厳密照合(主条件、設計書 §4.3)**:
   - 実際の駆動部 `iTPS<tensor>::simple_update(const EvolutionOperator&)` を、ゲート列の順に呼ぶ(ボンド処理を写したテスト内の関数で代用してはならない)。ゲートは T1 の `tenes_std` の出力(`infer_fermion_gate_ledgers()` を通したもの)。
   - 打ち切りが起きない設定にする。SU はボンド次元を保つので、ゲート列が通るボンドの仮想次元を、演算子を当てた後のランクより大きくとる。初期テンソルは、そのボンドの一部の添字(偶と奇を 1 つずつなど)だけに非零を持たせ、残りをゼロで埋める。λ は 1。
   - 更新後の単位胞テンソルから、経路を含む開放パッチ(外周の脚は固定のラベル)の物理状態 ψ' を作る。これが、同じパッチで更新前の状態 ψ に O_st を直接当てたものと、全体のスカラー倍を除いて一致することを確かめる。
   - 比べる量はスカラー倍に依存しないもの(例:複数の独立な bra 状態 φ_k との重なりの比 ⟨φ_k|ψ'⟩/⟨φ_0|ψ'⟩ を ⟨φ_k|O_st|ψ⟩/⟨φ_0|O_st|ψ⟩ と比べる。k は少なくとも 3)。相対誤差 1e-10 以下。
   - 参照値:d = 2 は `test/fermion/fock_oracle.py` で計算し、定数として埋め込む(生成手順をコメントに書く)。
     **oracle の制約**:物理次元 2(サイトあたり 1 モード)と、内部ボンドあたり偶奇 1 ラベルずつ(次元 2)にしか対応しない。更新前の状態は次元 2 の部分だけに非零を持つように作り、oracle にはその部分を渡す。
     d = 4 は、PR #117 で検証済みの長距離測定の機構(`build_relay_window` など)か、単層の graded 縮約で参照を作ってよい。どちらを使ったか、なぜ独立と言えるかを報告に書く。
   - 経路:設計書 §4.3 の一覧をすべて含める。少なくとも 1 ケースが `source_leg` 0 の段を、1 ケースが `source_leg` 1 の段を含むことを、テスト自身が確かめる。
   - 演算子:d = 2 の hopping と n_s n_t、d = 4 の Hubbard の hopping + 密度相互作用(T1 と同じもの)。real と complex の両方(complex は少なくとも 2 経路)。
2. **経路上で同じ単位胞サイトが再び現れる場合**:2×2 単位胞の (2, 0) を、上と同じ方法で照合する。この場合の開放パッチは、単位胞を展開したものとして作る(同じ単位胞サイトの 2 つの位置は、どちらも更新後のテンソルになる)。展開の扱いを報告に書く。
3. **台帳の持ち回し**:ゲート列の途中で、`iTPSTestAccessor::finfo(state).phys[m]` が推定した χ の台帳になり、末尾ゲートの後で元の台帳に戻る。左か上に進む段では、台帳が raster 化後のサイトに付くことを確かめる(raster 化後のサイトに TOML 上の out1/out2 を付ける、両者を混ぜた誤りで赤くなるケース)。
4. **入力時の推定と拒否**(`infer_fermion_gate_ledgers` を直接呼ぶ単体テスト):
   - 受理:T1 の出力のゲート列。χ の次元が d と等しいゲートを含むもの(T1 契約 3 のケース)。`complete_ungated_bonds()` の identity ゲートが末尾に付いた完全なリスト。
   - 推定値:受理したゲートの `fermion_legs` が、ゲートの非零要素から求めた台帳と一致する。
   - 拒否(それぞれ `tenes::input_error`):
     - `out2` の同じ添字でパリティが混在する
     - 非零要素を一つも持たない `out2` の添字がある
     - `in1` か `in2` の次元が、そのサイトの現在の台帳の長さと一致しない(太った脚に identity ゲートが当たる順序を含む)
     - `out1` の次元が source の元の物理次元と一致しない
     - 群の終わりに、あるサイトの台帳か次元が元に戻っていない(ゲート列が途中で切れている)
     - NN の奇ゲート(c†_s ⊗ 1_t など)。群末の検査で拒否される
     - 奇の一サイトゲート
     - `evolution.full` に物理次元が変わるゲートがあり、`num_full_step` のどれかが正。`num_full_step` がすべて 0 なら受理する
   - ボゾン(`fermion = false`)では何もしない(`fermion_legs` は空のまま、拒否もしない)。
5. **MPI**(n = 2 の登録):契約 4 の拒否が全ランクで起きる(ハングしない)。契約 1 の d = 2 の 2 ケース以上が n = 2 でも通る。テストは rank に依存しない決定的なテンソルを使う。
6. **既存の SU は変わらない**:NN だけのゲート列では、変更前と同じ結果になる(既存の `test_fermion_*`、`test_simple_update`、E2E が通ることで確かめる。新しいテストは不要)。
7. **変異に対する感度**(レビュアーが確認する):
   - `finfo.phys` の更新を止める
   - 台帳の更新で、raster 化後のサイト s1/s2 に TOML 上の `out1`/`out2` を付ける(両者を混ぜる)
   - raster 化の graded transpose で、4 本の台帳の入れ替えを忘れる
   - `out1` を推定に切り替え、`out2` を「次元が同じなら元の台帳」にする
   - 群末の検査を外す
   - 拒否を rank 0 だけで行う
8. 実行時間:Debug ビルド 1 スレッドで 60 秒以内を目安にする。d = 4 は経路の数を絞ってよいが、`source_leg` 0 と 1 の段を含む経路と、3 ホップの経路を 1 つずつ残す。

**ゲート**:契約 1・2・3 がすべて通ること。通らない場合は T3 以降に進まない。原因を特定して設計書を改訂する(統括がユーザーに報告する)。実装を「照合に合わせて」部分的に符号を足す修正は禁止。

#### 手順(T2)

- [ ] **Step 1: スタブを置く**(Claude)。`EvolutionOperator::fermion_legs`、4 台帳版の `wrap_twosite_gate`(本体は `throw std::logic_error("not implemented")`。既存の 3 引数版はまだ転送しない)、`infer_fermion_gate_ledgers`(本体は**何もせずに返す**。fermion モードで例外を投げると既存の fermion E2E がすべて落ちるため)を置き、`main.cpp` から呼ぶ。ビルドが通ること、既存の ctest が全件通ることを確認する。
  このスタブでは `fermion_legs` が空のままなので、契約 3・4 は推定値の欠如で、契約 1・2 は既存の SU が太った脚を固定台帳で包むことによる例外か値の不一致で落ちる。
- [ ] **Step 2: テスト作成者をディスパッチ**。振る舞い契約書(T2)節と設計書を渡す。報告 `work/fermion-longrange-ham/t2/test-author-report.md`。
- [ ] **Step 3: RED 確認**。`cmake --build --preset gcc --target test_fermion_longrange_gate` のあと `work/fermion-longrange-ham/t2/` から実行する。契約 1〜4 が Step 1 に書いた理由で落ちることを 1 件ずつ確認する(拒否テストのうち、現行の `validate_fermion_constraints` がすでに拒否するものは緑でよいが、どれがそうかを台帳に記録する)。oracle の定数が埋め込まれていること、`source_leg` 0・1 の自己検査があることを見る。スナップショットを取る。
- [ ] **Step 4: Codex に実装させる**。設計書 §6 と Global Constraints を渡す。次を明記する。
  - `main.cpp` の呼び出し順(`complete_ungated_bonds()` と `load_full_updates()` の後、`validate_fermion_constraints()` の前)
  - `validate_fermion_constraints` から外す検査(ゲートのパリティ検査)
  - 台帳の更新元:raster 化後の `fTn1_work.parity[4]` と `fTn2_work.parity[4]`(既存の仮想台帳の更新 `finfo.virt[s1][s1_leg] = fTn1_work.parity[s1_leg]` と同じ形)
  - SU の群ループの後に、`finfo.phys` が `peps_parameters.phys_parity` と一致することを確かめる内部不変条件(違反は `std::logic_error`)
  - `log_fermion_sector_dimensions` など、SU の途中で `finfo.phys` を読む処理が太った台帳で壊れないこと、測定・保存・CTM 更新が SU の群ループの内側で走らないことをコードで確かめ、報告に書く

  報告 `work/fermion-longrange-ham/t2/codex-report.md`。
- [ ] **Step 5: 不改変検査と独立検証**。スナップショット照合 → `test_fermion_longrange_gate` を実行 → 変異 2 件(台帳の更新を raster 化前で行う、群末の検査を外す)を試す → `ctest --preset gcc` 全件 → MPI ビルドで `ctest`(n = 1、n = 2)。
- [ ] **Step 6: ゲート判定**。契約 1・2・3 の結果を台帳に記録する。不合格なら停止してユーザーに報告する。
- [ ] **Step 7: タスクレビュー → 整形 → コミット**。

```bash
git add src/operator.hpp src/fermion/fops.hpp src/iTPS/load_toml.cpp src/iTPS/load_toml.hpp \
        src/iTPS/main.cpp src/iTPS/simple_update.cpp test/fermion/longrange_gate.cpp test/CMakeLists.txt
git clang-format
git commit -m "Carry the physical parity ledger through a chain of fermion gates"
```

---

### Task 3: `tenes_simple`(正方の 2・3 近接、三角格子)、E2E、ドキュメント

**Files:**
- Modify: `tool/tenes_simple.py`(`_check_fermion_scope`、`SpinlessFermionModel` と `HubbardModel` の docstring)
- Modify: `test/python/test_fermion_models.py`(三角格子と 2・3 近接の拒否テストを受理に更新。ハニカムとカゴメの拒否は残す。テスト作成者)
- Create: `test/fermion/free_fermion_tprime.py.in`、`test/fermion/free_fermion_triangular.py.in`、`test/fermion/hubbard_triangular.py.in`(テスト作成者)、golden `test/data/output_HubbardTriangular/`(自由フェルミオンの 2 本は厳密値と比べるので golden を持たせるかは作成者が判断)
- Modify: `test/CMakeLists.txt`(E2E の登録)
- Modify: `docs/sphinx/{ja,en}/file_specification/parameter_section.rst`、`simple_format.rst`、`std_format.rst`(該当する記述)、`NEWS.md`

**Interfaces:**
- Consumes:T1・T2 の成果物。既存の E2E の雛形(`test/fermion/free_fermion_longrange.py.in` の `exact_hopping`、`launcher`、golden 比較)。

#### 振る舞い契約書(T3)

1. **`tenes_simple` の受理範囲**:
   - fermion 模型(`spinless`、`hubbard`)で、正方格子の `t'`・`t''`・`v'`・`v''` と、三角格子(1〜3 近接)を受理し、std.toml を生成する。
   - ハニカムとカゴメは、引き続き拒否する(T5 で許可する)。
   - `[correlation_length]` の拒否は変わらない。
2. **パイプライン**:`tenes_simple` → `tenes_std` → `tenes` が、正方格子の t-t'(スピンレス)と三角格子(Hubbard)で最後まで走る。
3. **t-t' 自由フェルミオンの厳密照合**(`FreeFermionTPrime`)。判定点は Step 1 の実測(台帳 2026-09-30)で決めた:
   - スピンレス、正方格子 L = W = 2、t = 1、μ = −1、D = 3、χ = 18、τ = 0.05 × 400 ステップ、`initial = "random"`、シード固定。t' = +0.3 と t' = −0.3 の 2 回走らせる(各 85 秒程度)。
   - **μ = 0 は使わない**:粒子正孔対称性で E(t') = E(−t') となり、符号の誤りが見えない。
   - 参照は k 積分:ε(k) = −2t(cos kx + cos ky) − 4t' cos kx cos ky(`tenes_simple` の t' の符号規約。規約を `tool/tenes_simple.py` で確かめ、コメントに書く)。厳密値:E(+0.3) = −0.53007、E(−0.3) = −0.34272、n(+0.3) = 0.2837、n(−0.3) = 0.4275、⟨c†_0 c_(1,1) + h.c.⟩ = 0.1915(t' = +0.3)/ 0.0086(t' = −0.3)。
   - 判定(符号の誤りを確実に落とし、D = 3 の誤差には耐えるもの):
     - **対角方向の hopping**:std.toml に (1, 1) と (−1, 1) の対角ボンドの hopping 観測量(c†_s c_t + h.c.)を足して測る(`tenes_simple` の hopping 観測量は NN だけなので、E2E が std.toml に追記してよい)。t' = +0.3 で 0.1915 ± 0.08(実測 D = 3 で 0.175 前後。符号が逆の実装では 0.0086 付近になる)。
     - **エネルギー差**:E(+0.3) − E(−0.3) が厳密値 −0.187 の ±25% 以内(実測 D = 3 で −0.193)。
     - 密度:厳密値の ±0.05 以内(実測 D = 3 で 0.262〜0.267、t' = +0.3)。
   - 実装後に統括が、`tenes_std` のパリティ紐を外したコピーで作った入力で、このテストが落ちることを確かめる(変異)。
4. **三角格子の自由フェルミオン**(`FreeFermionTriangular`):
   - スピンレス、三角格子 L = W = 2、t = 1、μ = −1、D = 3、χ = 18、τ = 0.05 × 400、シード固定(70 秒程度)。三角格子の NN の (−1, 1) 方向は 2 ホップのゲート列になる。
   - 厳密値(ε(k) = −2t(cos kx + cos ky + cos(ky − kx))):E = −0.68303、n = 0.304、各 NN 方向の ⟨c†c + h.c.⟩ = 0.329。(−1, 1) の符号を逆にした模型では、その方向が −0.322 になる。
   - 判定:(−1, 1) 方向の hopping が 0.329 ± 0.1(実測 D = 3 でシード 1〜4 が 0.280〜0.290)。(1, 0)・(0, 1) 方向も正で 0.329 ± 0.12。エネルギーは D = 3 でシードによって −0.535〜−0.599 とばらつくので、厳密値との比較は緩く(±0.2)とどめ、符号の判定には使わない。
5. **三角格子 Hubbard**(`HubbardTriangular`):
   - 小さい D(D = 2、3×3 セル、ステップ数は 300 秒以内に収まるよう絞る。200 ステップで 107 秒の実測あり)。golden(`test/fulltest.py.in` の rtol/atol)で比べる。
   - golden だけにしない。符号の健全性として、3 方向すべての NN ボンドで hopping の期待値が正であることを確かめる。**3 方向のボンドエネルギーが等しいことは使わない**(SU の状態が格子の対称性を破るので、実測でも等しくない)。
6. **既存の E2E と golden は変わらない**。
7. **ドキュメント**:設計書 §6.1 のゲート列の規約(`out2` が太り `in1` が消費される)を `std_format.rst` か入力ファイル仕様に明記する。fermion の制限の一覧から「longer-distance Hamiltonian bonds」を外し、FU では使えないこと、ハニカムとカゴメが未対応であること(T5 まで)を書く。ja と en の両方。あわせて次の 2 点を書く。
   - 小さい D(実測では D = 2)でランダムな初期状態から長距離項を入れると、状態が真空に落ちて抜け出せないことがある。NN だけの計算で収束させた状態から始める(`tensor_save` / `tensor_load`)か、D を上げる。
   - 出力脚の添字に非零要素が一つもない二サイトゲート(射影演算子のような非可逆なゲート)は、fermion モードでは NN でも入力時に拒否される。exp(−τh) の形のゲートは影響を受けない。
8. 実行時間:各 E2E が `OMP_NUM_THREADS=1` で 300 秒以内。

#### 手順(T3)

- [ ] **Step 1: 実測**(Claude)。`work/fermion-longrange-ham/t3/probe/` で、T1・T2 の成果物を使って t-t' 自由フェルミオンを数点の D で走らせ、契約 3 の判定点が作れることを確かめる。SU のコスト(Hubbard 2 ホップの 1 ステップの時間)も測り、設計書 §10 のリスクを台帳に記録する。
- [ ] **Step 2: テスト作成者をディスパッチ**。契約書(T3)節、設計書、Step 1 の実測結果を渡す。golden はこの時点では作らない(実装後に統括が生成する)。報告 `work/fermion-longrange-ham/t3/test-author-report.md`。
- [ ] **Step 3: RED 確認**。E2E は `tenes_simple` の拒否で落ち、pytest の受理テストも拒否で落ちることを確認する。スナップショットを取る。
- [ ] **Step 4: Codex に実装させる**(`tool/tenes_simple.py` とドキュメント)。報告 `work/fermion-longrange-ham/t3/codex-report.md`。
- [ ] **Step 5: golden の生成**(Claude)。E2E のスクリプトの生成モードか手順に従って `test/data/output_*` を作る。手で編集しない。スナップショットを更新し、更新理由を台帳に書く。
- [ ] **Step 6: 独立検証**。`ctest --preset gcc` 全件、MPI ビルドで全件。変異:`tenes_std` のパリティ紐を外したコピーで `FreeFermionTPrime` と `FreeFermionTriangular` の入力を作って走らせ、落ちることを確認する。
- [ ] **Step 7: タスクレビュー → 整形 → コミット**。

```bash
git add tool/tenes_simple.py test/python/test_fermion_models.py test/fermion/free_fermion_tprime.py.in \
        test/fermion/hubbard_triangular.py.in test/data/output_FreeFermionTPrime test/data/output_HubbardTriangular \
        test/CMakeLists.txt docs/sphinx NEWS.md
black tool/tenes_simple.py test/python/test_fermion_models.py
git commit -m "Allow beyond-nearest-neighbour terms and the triangular lattice in fermionic simple mode"
```

---

### Task 4: D = 1 の仮想脚と空孔での fermion 経路の検証(フェーズ 2 のゲート)

**Files:**
- Create: `test/fermion/unit_d1.cpp`(テスト作成者)
- Modify: `test/CMakeLists.txt`(実行ファイル `test_fermion_d1`。テスト作成者)
- Modify(不具合があれば):`src/fermion/`・`src/iTPS/` の該当箇所

**Interfaces:**
- Consumes:`fock_oracle.py`(`Patch`、`Oracle`。D = 1 の脚は、パリティ表 `[False]` の脚として扱える)、`build_relay_window`、`detail::doubled_pipeline`、`core::Contract_density_CTM`、`iTPSTestAccessor`、`iTPS::simple_update`。

#### 振る舞い契約書(T4)

D = 1 の仮想脚(台帳 `[偶]`)と、物理次元 1 のサイト(空孔、台帳 `[偶]`、仮想脚はすべて D = 1)を含む配置で、既存の fermion 経路が正しいことを確かめる。ハニカムとカゴメの埋め込み(設計書 §7.1、`tool/tenes_simple.py` の `HoneycombLattice` と `KagomeLattice`)を模した単位胞と開放パッチを使う。

1. **SU のボンド処理**:片方か両方のサイトが D = 1 の脚を持つボンドに、NN の hopping ゲートを当てた結果を、T2 契約 1 と同じ方法(スカラー倍に依らない量、oracle の定数)で照合する。
2. **長距離ゲート列**:ハニカムの NNN(D > 1 のボンドに沿った 2 ホップ)と、カゴメの NN の (−1, 1)(A を経由する 2 ホップ)のゲート列を、T2 契約 1 と同じ方法で照合する。ゲートは `tenes_std` で、ハニカム・カゴメの単位胞の std.toml から生成する(`tenes_simple` はまだ拒否するので、std.toml はテスト作成者が書く)。
3. **測定**:CTM(次元 1 の環境で閉じた窓の厳密縮約でよい)と MF の、NN の二サイト測定。長距離測定のリレーが D = 1 のボンドや空孔を横切る場合(x_first の経路がそうなる配置を選ぶ)。いずれも oracle の ⟨c†_s c_t⟩ と相対誤差 1e-10 以下で一致する。
4. **相関関数**:鎖が D = 1 の脚を通る配置で、`[correlation]` の値が oracle と一致する。
5. **保存と読み込み**:D = 1 の脚と空孔を含む単位胞で、保存して読み込んだ状態の Tn・台帳・一サイト測定値が一致する。
6. 各項目は独立したテストケースにする。一部が赤でも、他の項目の合否が分かるようにする。

**ゲート**:すべて通ること。赤の項目がある場合は、統括が原因を特定する。修正が局所的なら、このタスクで修正(Codex)とテストの再実行を行う。修正が設計に関わる規模なら、ユーザーに報告して、フェーズ 2 を別ブランチに切り出すかを判断してもらう。

#### 手順(T4)

- [ ] **Step 1: テスト作成者をディスパッチ**。契約書(T4)節と設計書 §7 を渡す。報告 `work/fermion-longrange-ham/t4/test-author-report.md`。
- [ ] **Step 2: 実行と分類**(Claude)。このタスクはゲートなので、実装前の RED を前提にしない。テストを実行し、緑の項目と赤の項目を分ける。赤の項目は、テストの欠陥か実装の不具合かを判定する(テストの欠陥はテスト作成者に差し戻す)。スナップショットを取る。
- [ ] **Step 3: 不具合があれば Codex に修正させる**。原因の特定結果と、赤のテストを渡す。報告 `work/fermion-longrange-ham/t4/codex-report.md`。不具合がなければこのステップは飛ばし、台帳にその旨を書く。
- [ ] **Step 4: 不改変検査と独立検証**。スナップショット照合 → `test_fermion_d1` → `ctest --preset gcc` 全件 → MPI ビルドで全件。
- [ ] **Step 5: ゲート判定**。結果を台帳に記録する。
- [ ] **Step 6: タスクレビュー → 整形 → コミット**。

```bash
git add test/fermion/unit_d1.cpp test/CMakeLists.txt   # 修正があれば src/ の該当ファイルも
git clang-format
git commit -m "Verify fermion paths on unit cells with dimension-1 bonds and vacancies"
```

---

### Task 5: `tenes_simple`(ハニカム・カゴメ)、E2E、ドキュメント

**Files:**
- Modify: `tool/tenes_simple.py`(`_check_fermion_scope`)
- Modify: `test/python/test_fermion_models.py`(ハニカムとカゴメの拒否テストを受理に更新。テスト作成者)
- Create: `test/fermion/free_fermion_honeycomb.py.in`、`test/fermion/free_fermion_kagome.py.in`、`test/fermion/hubbard_honeycomb.py.in`(テスト作成者)、対応する golden
- Modify: `test/CMakeLists.txt`、`docs/sphinx/{ja,en}/...`、`NEWS.md`

**Interfaces:**
- Consumes:T3 の E2E の雛形、T4 の結果。

#### 振る舞い契約書(T5)

1. **`tenes_simple` の受理範囲**:fermion 模型で、ハニカムとカゴメ(1〜3 近接)を受理する。空孔には `parity = [0]`、`physical_dim = 1` が出力される。
2. **ハニカムの自由フェルミオン**(`FreeFermionHoneycomb`):スピンレス、NN の t と NNN の t'。参照は 2 バンドの k 積分(グラフェン型の分散に t' の項を足したもの。`tenes_simple` の t' の規約に合わせる)。判定点と許容誤差は T3 契約 3 と同じ考え方で決める(t' の符号反転によるエネルギー差の 1/4 以下)。
3. **カゴメの自由フェルミオン**(`FreeFermionKagome`):スピンレス、NN の t だけ(NN が 2 ホップのゲート列を含むことを使う)。参照は 3 バンド(平坦バンドを含む)の k 積分。μ は平坦バンドを避けた値。許容誤差は、NN の 2 ホップのゲートの符号を反転させたとき(奇チャネルに −1 を掛けた場合)のエネルギー差の 1/4 以下になるよう決める。
4. **ハニカム Hubbard**(`HubbardHoneycomb`):小さい D で走る設定。golden に加えて、厳密に成り立つ関係を 1 つ以上確かめる(3 方向の NN ボンドエネルギーが等しい、など)。
5. **ドキュメント**:格子の制限の記述を更新する(fermion で使える格子に、ハニカムとカゴメを加える)。ja と en の両方。NEWS。
6. 実行時間:各 E2E が `OMP_NUM_THREADS=1` で 300 秒以内。

#### 手順(T5)

- [ ] **Step 1: 実測**(Claude)。`work/fermion-longrange-ham/t5/probe/` で、ハニカムとカゴメの自由フェルミオンを数点の D で走らせ、判定点が作れることを確かめる。
- [ ] **Step 2: テスト作成者をディスパッチ**(契約書(T5)節、設計書、Step 1 の結果)。報告 `work/fermion-longrange-ham/t5/test-author-report.md`。
- [ ] **Step 3: RED 確認**(`tenes_simple` の拒否で落ちること)。スナップショットを取る。
- [ ] **Step 4: Codex に実装させる**。報告 `work/fermion-longrange-ham/t5/codex-report.md`。
- [ ] **Step 5: golden の生成**(Claude)。
- [ ] **Step 6: 独立検証**。`ctest --preset gcc` 全件、MPI ビルドで全件。変異:カゴメの NN の 2 ホップで、`tenes_std` のパリティ紐を外したコピーを使い、`FreeFermionKagome` が落ちることを確認する。
- [ ] **Step 7: タスクレビュー → 整形 → コミット**。

```bash
git add tool/tenes_simple.py test/python/test_fermion_models.py test/fermion/free_fermion_honeycomb.py.in \
        test/fermion/free_fermion_kagome.py.in test/fermion/hubbard_honeycomb.py.in test/data/output_* \
        test/CMakeLists.txt docs/sphinx NEWS.md
black tool/tenes_simple.py test/python/test_fermion_models.py
git commit -m "Allow the honeycomb and kagome lattices in fermionic simple mode"
```

---

## 最終段

- [ ] 最上位モデルで全ブランチのレビューを 1 回行う(タスク単位では見えない欠陥、どのテストも通らない経路を探す)。指摘の処置を台帳に記録する。
- [ ] `ctest` 全件(非 MPI・MPI)を統括が回す。
- [ ] PR 本文を `work/fermion-longrange-ham/pr/PR-BODY.md` に書く。develop に対する差分だけを書き、ブランチ内の経緯は書かない。「Generated with Claude Code and Codex」を明記する。
- [ ] ユーザーの承認を得てから push と PR 作成を行う。
