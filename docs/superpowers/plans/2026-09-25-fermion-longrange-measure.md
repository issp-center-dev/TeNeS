# fermion 長距離二サイト物理量・相関関数 実装計画

> **For agentic workers:** この計画は CLAUDE.md の多段エージェント手順で実行する。
> 振る舞い契約書(散文)→ テスト作成者(Claude サブエージェント)→ RED 確認とスナップショット →
> Codex 実装 → テスト不改変の機械検査 → Claude の独立検証 → タスクレビュー → 整形とコミット。
> ステップは `- [ ]` で追跡する。**テストコードはこの計画に書かない**(テスト作成者に計画者のテストを見せないため)。

**Goal:** fermion モードで、最近接以外の二サイト物理量(窓 4×4 まで)、ops 形式の二サイト物理量、奇の一サイト演算子、相関関数を、CTM 環境と MF 環境の両方で測れるようにする。

**Architecture:** graded SVD の各チャネルの k 脚を次元 1 の補助脚 κ とし、窓の中で仮想ボンドに沿ってリレーする。
各サイトは既存の `doubled_pipeline` で fold され、既存の density カーネル(`Contract_density_CTM`、`Start/Transfer/FinishCorrelation_density_CTM`)に渡る。
MF は χ = 1 の合成環境で同じ経路に乗せる。

**Tech Stack:** C++17, mptensor, doctest(単体)、Python 3.9+ / pytest / numpy(ツールと E2E)、ctest。

**Spec:** `docs/superpowers/specs/2026-09-25-fermion-longrange-measure-design.md`(以下「設計書」)。実行者は設計書とこの計画の両方を読むこと。

## Global Constraints

- ブランチ `fermion-longrange`。実験・台帳・ディスパッチ文面・報告は `work/fermion-longrange/<topic>/` に置き、CWD をそこにして実行する。リポジトリ直下やビルドツリーで `tenes` やテストバイナリを実行しない。
- C++ はすべて namespace `tenes`(fermion 層は `tenes::fermion`)。新規ファイルの先頭に GPL v3 ヘッダを付ける。
- Python の最低バージョンは 3.9(CI)。3.10 以降の API(`int.bit_count()`、`match` など)を使わない。
- fermion モードは「未対応の入力は読み込み時に拒否する」方針を守る(設計書 §5.1)。
- 窓は |dx| ≤ 3 かつ |dy| ≤ 3。(dx, dy) = (0, 0) は拒否。
- 一サイト演算子は偶か奇のどちらかに確定していること。混在は拒否。奇の一サイト演算子の一サイト測定値は厳密に 0。
- A_s B_t の内部配列: `op4[i_s, i_t, o_s, o_t] = (-1)^{p_B · p(i_s)} · A[i_s, o_s] · B[i_t, o_t]`(設計書 §4.2)。
- 交差符号は経路の辺ごとに正確に 1 回、上流側の出口融合脚に `(-1)^{p_bond(b)·p_k}`(設計書 §3.3)。
- 既存の bundled-k(`build_reduced_pair_halves*`)と最近接の測定経路は変更しない。
- 整形: C++ は `git add` のあと `git clang-format`(変更行のみ)。ファイル全体に `clang-format -i` を当てない。Python は `black`(line-length 88)。整形は Claude がコミット直前に行う。Codex には実行させない。
- ビルドとテスト: 非 MPI は `cmake --preset gcc && cmake --build --preset gcc && ctest --preset gcc`。MPI は `out-gcc-mpi/build`。
- golden ファイルは意図的に再生成する。手で編集しない。

## Review Focus

仕様が暗に要求しているのに、タスク内のテストが自然には踏まない入力。各行のテストは括弧内のタスクに追加してある。

1. **単位胞より大きい窓**: 2×2 単位胞で dx = 2 や dx = 3 のように、窓の中で同じサイトが複数回現れる場合。target が source と同じ副格子の場合も含む。環境テンソルの割り当てと紐付きテンソルの取り違えが起きやすい(T2・T3)。
2. **skew のある単位胞**: 長距離の窓は `lattice.other` で組むので、skew で窓の配置が変わる。展開した skew 0 のセルと一致すること(T2)。
3. **complex テンソル**: `complex_tensor` での値。conj の抜けは実数テストでは見えない(T1・T2)。
4. **負の変位**: source が target の右や下にある場合。窓の中で source の位置が反転する(T1 は幾何、T2 は本番の配線)。
5. **偶と奇のチャネルが混在する演算子**: hopping + nn のように、偶チャネルと奇チャネルを両方持つ d = 4 の演算子。チャネルの和の取りこぼしや二重計上(T2)。

---

## ファイル構成

| ファイル | 役割 | タスク |
|---|---|---|
| `src/fermion/relay.hpp`(新規) | 紐リレーの経路、チャネル分解、サイトと窓の builder、A_s B_t の 4 脚行列 | T1, T2 |
| `src/fermion/parity.hpp` | 一サイト演算子のパリティ分類(偶・奇・混在、全ランク一致) | T2 |
| `src/iTPS/iTPS.hpp` | 一サイト演算子のパリティ表(メンバ)、長距離測定の補助メソッド宣言 | T2, T3 |
| `src/iTPS/twosite_obs.cpp` | fermion 長距離の窓測定、ops 形式の変換 | T2, T5 |
| `src/iTPS/onesite_obs.cpp` | 奇の一サイト演算子の 0 出力 | T2 |
| `src/iTPS/correlation_function.cpp` | fermion の相関関数(CTM・MF) | T3, T5 |
| `src/iTPS/load_toml.cpp` / `load_toml.hpp` | 入力ガードの変更 | T2, T3, T5 |
| `src/iTPS/measure.cpp` | `validate_fermion_ctm_measurement` の変更 | T2, T3, T5 |
| `test/fermion/relay_window.cpp`(新規) | T1 の単体テスト | T1 |
| `test/fermion/longrange_measure.cpp`(新規) | T2・T3 の単体テスト | T2, T3 |
| `test/fermion/longrange_mf.cpp`(新規) | T5 の単体テスト | T5 |
| `test/CMakeLists.txt` | 新しい実行ファイル `test_fermion_relay`(T1)、`test_fermion_longrange`(T2・T3・T5)の登録、E2E 登録 | T1–T5 |
| `test/fermion/fermion_guards.cpp` ほか既存テスト | 拒否していた入力を受理に変える既存テストの更新(テスト作成者が行う) | T2, T3 |
| `test/fermion/fold_geometry.cpp` | T14a–e(長距離拒否の回帰テスト)を新仕様に合わせる(テスト作成者が行う) | T2 |
| `tool/tenes_simple.py`, `tool/tenes_std.py` | 演算子の追加、`[correlation]` の受理、ops 形式の受理 | T4 |
| `test/python/test_tenes_simple.py`, `test_tenes_std.py`, `test_fermion_models.py` | ツールのテスト | T4 |
| `test/fermion/free_fermion_longrange.py.in`(新規)、`boson_equivalence_longrange.py.in`(新規)、`hubbard_longrange.py.in`(新規) | E2E | T4 |
| `test/data/output_FreeFermion*/` ほか | golden の再生成 | T4 |
| `docs/sphinx/{ja,en}/...`、`NEWS.md`、`sample/07_*`・`sample/08_*` | ドキュメントとサンプル | T4 |

## 共通手順(全タスク)

各タスクの「手順」欄は、この共通手順を具体化したものである。ディスパッチ文面には以下の定型注意を**最初のタスクから**入れる。

**テスト作成者(Claude サブエージェント、`general-purpose`)への定型注意**
- 渡すのはそのタスクの「振る舞い契約書」節と設計書だけ。計画者のテストコードは存在しない。
- 変更してよいのは契約書が指定したテストファイルと `test/CMakeLists.txt` だけ。`src/` と `tool/` は変更禁止。
- `git stash` / `git checkout` / `git reset` を含め、作業ツリーを巻き戻す操作は禁止。仮説の検証はファイルを scratchpad にコピーして行う。
- 契約書の誤り(不可能な自己検査、空洞な参照値、符号の誤りなど)に気づいたら、テストを曲げずに報告する。
- 参照値は、テスト対象の機構を使わずに得ること(単層 graded 縮約、Fock oracle、解析解)。
- 報告は `work/fermion-longrange/<task>/test-author-report.md` に書く。存在しないこと自体が欠陥。

**Codex(実装)への定型注意**
- `codex-companion.mjs task --background --fresh --write --model gpt-5.5 "$(cat <prompt file>)"` で起動する。プロンプトは必ずファイル経由で渡す。
- 「確認を求めずに実装完了まで一気に進めること」と明記する。
- テストファイル(スナップショットを取ったもの)は変更禁止。テストが誤りと思ったら実装を曲げず、BLOCKED として報告する。
- formatter は実行しない。新しい行は周囲のスタイルに手で合わせる。
- sandbox はコミットできない。コミットは Claude が代わりに行う。
- `tenes` やテストバイナリをリポジトリ直下で実行しない。
- 報告は `work/fermion-longrange/<task>/codex-report.md` に書く。存在しないこと自体が欠陥。
- 完了後、Claude は `status --json` の running と pid の生存を確認する。`git status` で変更がゼロなら承認待ちを疑い、`--resume` で「はい」を返す。

**RED 確認とスナップショット(Claude)**
1. テストがコンパイルでき、実装前は**正しい理由で**失敗することを 1 件ずつ確認する。インターフェースのスタブ(Claude が事前に置く)は「未実装」例外を投げるか空を返すので、失敗理由はそれになるはずである。偶然通るテスト(空洞)がないか見る。
2. `sha256sum` でテストファイルのスナップショットを `work/fermion-longrange/<task>/test-snapshot.sha256` に取る。
3. Codex の完了後、スナップショットと照合して不一致なら差し戻す。テストの欠陥はテスト作成者に差し戻す。

**独立検証(Claude)**
- Codex の報告を鵜呑みにしない。該当テストバイナリを自分で回す。共有コード(`src/iTPS/`、`src/fermion/`、`tool/`)に触れたタスクでは、最後に `ctest` 全件を自分で回す。MPI ビルドも回す。
- 変異テストを自分で 1〜2 件試す(契約書の変異リストから選び、scratchpad のコピーで)。

**タスクレビュー**: タスクごとに新しいレビュアー(`feature-dev:code-reviewer`)。作業ツリーの変更禁止を明記し、変異テストを推奨する。

**整形とコミット**: `git add` → `git clang-format` → `black` → テスト → コミット。コミットメッセージの末尾に帰属行を付ける。

---

### Task 1: 紐リレーの builder と厳密照合のゲート

**Files:**
- Create: `src/fermion/relay.hpp`
- Create: `test/fermion/relay_window.cpp`(テスト作成者)
- Modify: `test/CMakeLists.txt`(実行ファイル `test_fermion_relay` を登録。テスト作成者)

**Interfaces:**
- Consumes(既存): `ftensor<tensor>`、`transpose` / `tensordot` / `svd` / `conj`(`src/fermion/fops.hpp`)、`detail::doubled_pipeline`(`src/fermion/reduced.hpp`)、`wrap_Tn` / `wrap_twosite_gate`。
- Produces(T2・T3・T5 が使う。Claude が実装前にスタブとして置く):

```cpp
namespace tenes::fermion {

//! Cell of a measurement window in Contract_* orientation:
//! row 0 is the top row, col 0 the left column.
struct window_cell {
  int row;
  int col;
};

enum class relay_order { x_first, y_first };
enum class relay_role { source, middle, target };

//! Cells from source to target inclusive, one nearest-neighbour step at a
//! time. x_first walks along the source row to the target column, then
//! along that column; y_first the other way round.
//! @throw std::invalid_argument if source == target.
std::vector<window_cell> relay_path(window_cell source, window_cell target,
                                    relay_order order);

//! Virtual leg (0 = l, 1 = t, 2 = r, 3 = b) of `from` pointing at the
//! adjacent cell `to` (row - 1 is up = t).
//! @throw std::invalid_argument if the cells are not adjacent.
int relay_leg(window_cell from, window_cell to);

//! One channel of a two-site operator: u has legs (in, out, kappa) and vt
//! (kappa, in, out); kappa has dimension 1 and a definite parity; the
//! singular value is folded into u.
template <class tensor>
struct relay_channel {
  ftensor<tensor> u;
  ftensor<tensor> vt;
};

//! Graded SVD of a wrap_twosite_gate()-loaded operator (legs in1, in2,
//! out1, out2), split into dimension-1 channels. Channels with a zero
//! singular value may be omitted.
template <class tensor>
std::vector<relay_channel<tensor>> relay_channels(const ftensor<tensor>& op12);

//! One folded site of the relay network: plain rank-6 tensor
//! ([l lb], [t tb], [r rb], [b bb], s_ket, s_bra).
//! source: entry_leg == -1; target: exit_leg == -1; middle: both in 0..3
//! and distinct. Applies the crossing sign on the exit leg (source, middle).
//! @throw std::invalid_argument on an inconsistent role / leg combination.
template <class tensor>
tensor build_relay_site(const ftensor<tensor>& Tn, relay_role role,
                        int entry_leg, int exit_leg,
                        const relay_channel<tensor>& channel);

//! Folded window for one channel: [row][col] rank-6 tensors. Cells on the
//! path come from build_relay_site(); the others are
//! detail::doubled_pipeline(Tn, Tn).
template <class tensor>
std::vector<std::vector<tensor>> build_relay_window(
    const std::vector<std::vector<ftensor<tensor>>>& Tn, window_cell source,
    window_cell target, const relay_channel<tensor>& channel,
    relay_order order = relay_order::x_first);

}  // namespace tenes::fermion
```

#### 振る舞い契約書(T1)

目的: 設計書 §3.5 のゲートを機械的に判定できるテストを作る。このタスクの合否がそのまま方式の採否になる。

1. **経路** `relay_path`:
   - 両端を含み、長さは |Δrow| + |Δcol| + 1。隣り合うセルは必ず最近接。
   - x_first は source の行を target の列まで進み、そこから target の行まで進む。y_first は逆。
   - source == target は `std::invalid_argument`。
   - `relay_leg` は右隣に 2、左隣に 0、上隣(row − 1)に 1、下隣に 3 を返し、隣接していなければ `std::invalid_argument`。
2. **チャネル分解** `relay_channels`: 全チャネルの u_k と vt_k を κ で graded 縮約して足すと、元の `op12` に戻る(相対誤差 1e-12)。各 κ は次元 1 で、パリティ表の長さは 1。d = 2 の hopping(偶と奇のチャネルを持つ)と、d = 4 の hopping + nn(偶奇混在)で確かめる。
3. **厳密照合(主条件)**: 外周の脚を次元 1(偶)に閉じた開いたパッチの上で、決定的なパリティ偶のサイトテンソル
   (`test/fermion/fold_geometry.cpp` の決定的テンソルと同じ式、実数と複素数)を使う。次の 2 つの量を比べる。
   - 真値 ⟨O⟩ = ⟨ψ|O|ψ⟩ / ⟨ψ|ψ⟩。単層 graded 縮約(`tenes::fermion::tensordot` など、fold を使わない)で計算する。
   - relay 値 = Σ_k [窓 k の厳密縮約] / [恒等窓の厳密縮約]。恒等窓は全セル `doubled_pipeline(Tn, Tn)`。
     窓の厳密縮約は、builder を使わない任意の方法で行ってよい。たとえば、次元 1 の C・eT を 1 で埋めて `core::Contract_density_CTM` に渡し、全セルに恒等の一サイト演算子を与える。あるいは素の mptensor 縮約でもよい。

   両者の相対誤差は 1e-12 以下であること。演算子は d = 2 の ⟨c†_s c_t⟩ と hopping、d = 4 の hopping + nn(Global Constraints の op 規約で与える 4 脚行列)。ケースは次の**全部**を含める:
   - 直進: 1×3 の (2,0)、1×4 の (3,0)、3×1 の (0,2)、1×3 の (−2,0)
   - 曲がり角: 2×2 の (1,1)・(1,−1)・(−1,1)・(−1,−1)、3×2 の (2,1)・(−2,1)
   - 3 本以上の非自明な脚を持つ途中サイト: 3×3 パッチで、中央を**直進**で通る経路(例: 上端中央 → 下端中央)と、中央で**曲がる**経路(例: 左端中央 → 上端中央)
   - 窓の中で source が target の右または下にある配置(負の変位)を上のリストが含んでいることを確かめる
4. **Fock oracle の錨**: d = 2 の少なくとも 3 ケース((2,0)、(1,1)、3×3 の曲がり角)の ⟨c†_s c_t⟩ を `test/fermion/fock_oracle.py`(`Oracle.one_body(i, j)` と `norm()`)で計算し、定数として埋め込む。relay 値との相対誤差は 1e-12 以下。生成手順(コマンドと引数)をテストのコメントに書く。
5. **経路非依存性**: 3 のケースのうち曲がり角のもの全部で、`relay_order::x_first` と `y_first` の relay 値が相対誤差 1e-12 以下で一致する。
6. **最近接の一致**: |dx| + |dy| = 1 の 4 方向(source が左・右・上・下)で、relay 値が既存の bundled-k 構成(`build_reduced_pair_halves` の両半分を縮約したもの)の値と相対誤差 1e-12 以下で一致する。3×3 パッチの中央付近の対も 1 件含める。
7. **入力検査**: `build_relay_site` は、source なのに entry_leg ≠ −1、target なのに exit_leg ≠ −1、middle で entry == exit または範囲外、のそれぞれで `std::invalid_argument` を投げる。
8. **変異に対する感度**(テスト作成者は変異を実行しなくてよいが、次の変異で必ず赤くなるケース選定にすること。レビュアーが確認する):
   κ の graded transpose を素の transpose に替える(全サイト / 途中サイトだけ)、交差符号を落とす、交差符号を両側に掛ける、融合順を (κ, bond) にする、途中サイトで κ_in と κ_out を入れ替える、途中サイトの S を物理脚の前に外積する、奇チャネルを落とす。
   「交差符号を下流側に移す」は値の変わらない等価変異なので対象外。
9. 実行時間: Debug ビルド 1 スレッドで 30 秒以内を目安にする。d = 4 の 3×3 は D = 2 に抑えてよい。

**ゲート**: 3・4・5・6 がすべて通ること。通らない場合は T2 以降に進まない。原因を特定して設計書を改訂する。builder を部分的に「測定で合わせる」修正は禁止(設計書 §3.5)。

#### 手順(T1)

- [ ] **Step 1: スタブを置く**(Claude)。`src/fermion/relay.hpp` に上の宣言を置き、本体はすべて `throw std::logic_error("relay: not implemented")` にする。GPL ヘッダと `@file` の Doxygen を付ける。ビルドが通ることを確認する。
- [ ] **Step 2: テスト作成者をディスパッチ**。振る舞い契約書(T1)節と設計書を渡す。報告 `work/fermion-longrange/t1/test-author-report.md`。
- [ ] **Step 3: RED 確認**。`cmake --build --preset gcc --target test_fermion_relay` のあと `work/fermion-longrange/t1/` から実行し、全ケースが "relay: not implemented" で落ちることを確認する(Fock 錨の定数が埋め込まれていること、ケース一覧が契約の 3 を網羅していることも見る)。スナップショットを取る。
- [ ] **Step 4: Codex に実装させる**。設計書 §3.3 の脚順と交差符号、Global Constraints を渡す。報告 `work/fermion-longrange/t1/codex-report.md`。
- [ ] **Step 5: 不改変検査と独立検証**。スナップショット照合 → `test_fermion_relay` を実行 → 変異 2 件(交差符号を落とす、κ_in と κ_out の入れ替え)を scratchpad のコピーで試し、赤くなることを確認する → `ctest --preset gcc` 全件。
- [ ] **Step 6: ゲート判定**。契約の 3・4・5・6 の結果を台帳 `work/fermion-longrange/LEDGER.md` に記録する。不合格なら停止してユーザーに報告する。
- [ ] **Step 7: タスクレビュー → 整形 → コミット**。

```bash
git add src/fermion/relay.hpp test/fermion/relay_window.cpp test/CMakeLists.txt
git clang-format
git commit -m "Relay the odd channel of a fermionic two-site operator along virtual bonds"
```

---

### Task 2: CTM での長距離二サイト物理量、ops 形式、奇の一サイト演算子、ガード

**Files:**
- Modify: `src/fermion/relay.hpp`(`product_twosite_op` を追加)
- Modify: `src/fermion/parity.hpp`(`operator_parity` を追加)
- Modify: `src/iTPS/iTPS.hpp`, `src/iTPS/iTPS.cpp`(一サイト演算子のパリティ表を構築時に作る)
- Modify: `src/iTPS/twosite_obs.cpp`(fermion の最近接以外と ops 形式)
- Modify: `src/iTPS/onesite_obs.cpp`(奇 → 0)
- Modify: `src/iTPS/load_toml.cpp`, `load_toml.hpp`, `src/iTPS/measure.cpp`(ガード)
- Create: `test/fermion/longrange_measure.cpp`(テスト作成者)
- Modify: `test/CMakeLists.txt`(`test_fermion_longrange` を登録)、既存のガードテスト(`test/fermion/fermion_guards.cpp`、`fold_geometry.cpp` の T14a–e など。テスト作成者)

**Interfaces:**
- Consumes: T1 の全シグネチャ。
- Produces:

```cpp
namespace tenes::fermion {

enum class op_parity { even, odd, mixed };

//! Parity class of a one-site operator op[in, out] under the physical
//! ledger `phys`. The zero operator is even. Collective over the tensor's
//! communicator: every rank returns the same value.
template <class tensor>
op_parity operator_parity(const tensor& op, const parity_vector& phys);

//! Plain two-site operator of the product A_s B_t in the ordered basis
//! |n_s n_t>: op4[i_s, i_t, o_s, o_t] = (-1)^{p_B p(i_s)} A[i_s, o_s]
//! B[i_t, o_t]. B must have a definite parity (odd_B).
template <class tensor>
tensor product_twosite_op(const tensor& A, const tensor& B,
                          const parity_vector& phys_s,
                          const parity_vector& phys_t, bool odd_B);

}  // namespace tenes::fermion
```

- iTPS 側(private、`iTPSTestAccessor` から見える): `std::vector<fermion::op_parity> onesite_parity;`(`onesite_operators` と同じ添字)。

#### 振る舞い契約書(T2)

1. **入力の受理と拒否**(`load_toml` の `validate_fermion_constraints`、`tenes::input_error`):
   - 受理する: 0 < |dx| + |dy| で |dx| ≤ 3 かつ |dy| ≤ 3 の二サイト物理量、ops 形式の二サイト物理量、奇の一サイト演算子。
   - 拒否する: (dx, dy) = (0, 0)、|dx| > 3 または |dy| > 3(メッセージに "4x4" を含む)、偶奇混在の一サイト演算子(メッセージに "mixed parity" を含む)、パリティの異なる 2 つの一サイト演算子を組む ops 形式、奇の一サイト演算子を使う一サイトゲート(従来どおり)、パリティ奇の二サイト物理量(従来どおり)、multisite(従来どおり)。
   - `meanfield_env = true` では、最近接以外の二サイト物理量と、最近接以外を指す ops 形式を拒否する(T5 で解除)。最近接の ops 形式は MF でも受理する。
   - `validate_fermion_ctm_measurement` も同じ条件に揃える。
2. **`operator_parity`**: 偶だけ・奇だけ・混在・零行列(偶)を正しく分類する。d = 2 と d = 4 の両方で確かめる。
3. **`product_twosite_op` の符号**: T1 と同じ開いたパッチの厳密照合の枠組みで、`product_twosite_op(c†, c)` を `wrap_twosite_gate` → `relay_channels` → `build_relay_window` に通す。その値が Fock oracle の ⟨c†_s c_t⟩ に一致すること。さらに `product_twosite_op(c, c†)` の値が −⟨c†_t c_s⟩ に一致すること。(1,0) と (2,0) の両方で確かめる。偶の積(n_s n_t)も 1 件入れる。
4. **一サイト測定**: 奇の一サイト演算子について、`measure_onesite` の値は実部・虚部とも厳密に 0.0。偶の演算子の値は変更前と同じ。
5. **窓測定の配線**: 小さな fermion の `iTPS`(`iTPSTestAccessor` で構築、2×2 単位胞、D = 2、d = 2 と d = 4)で `update_CTM()` を行ったあと、`measure_twosite()` の値を検査する。(dx, dy) ∈ {(2,0), (0,2), (1,1), (−1,2), (3,0), (2,−1)} のそれぞれで、同じ環境テンソルを使って T1 の builder と `core::Contract_density_CTM` を直接呼んだ値(テスト内で窓を組み立てる)と、相対誤差 1e-12 以下で一致すること。環境テンソルの割り当てはボソンの長距離窓(`twosite_obs.cpp` の CTM 分岐)と同じ規約である。
6. **最近接の一致**: (1,0)、(−1,0)、(0,1)、(0,−1) の hopping で、この経路の値と既存の最近接経路の値が相対誤差 1e-10 以下で一致すること。本番はこれまでどおり最近接に既存の経路を使うので、テストは T1 の builder で窓を組んだ値と `measure_twosite()` の値を比べる。
7. **ops 形式**: `ops = [i, j]` で指定した二サイト物理量が、`product_twosite_op` で作った明示形式の同じ物理量と、(1,0) と (2,1) で一致すること。
8. **Review Focus 1(単位胞より大きい窓)**: 2×2 単位胞で (2,0)(target が source と同じ副格子)と (3,1) を測り、5 と同じ直接計算に一致すること。
9. **Review Focus 2(skew)**: skew = 1 の [2,1] セルでの (2,0) と (1,1) の値が、展開した skew 0 のセルでの対応する値と一致すること(`test/fermion/skew_unfold.cpp` と同じ型。CTM は有限 χ の精度で比べ、許容誤差はそのファイルと同じ考え方で決める)。
10. **Review Focus 3(complex)**: 5 の一部を `complex_tensor` でも行う。
11. **Review Focus 5(偶奇混在チャネル)**: d = 4 の hopping + nn の (2,0) で、値が hopping 単独と nn 単独の和に一致すること(相対誤差 1e-12)。
12. **既存テストの更新**: 長距離・ops 形式・奇の一サイト演算子の「拒否」を確かめていた既存のテストは、新仕様(受理、または新しい拒否条件)を確かめるように書き換える。書き換えた箇所を報告に列挙する。

#### 手順(T2)

- [ ] **Step 1: スタブを置く**(Claude)。`operator_parity` と `product_twosite_op` のスタブ(`logic_error`)と、iTPS の `onesite_parity` メンバ(空のまま)。ビルドを確認する。
- [ ] **Step 2: テスト作成者をディスパッチ**。報告 `work/fermion-longrange/t2/test-author-report.md`。
- [ ] **Step 3: RED 確認とスナップショット**。ガードのテストは「まだ拒否される」ことで、測定のテストは未実装例外で落ちることを確認する。
- [ ] **Step 4: Codex に実装させる**。実装方針: fermion で最近接以外なら、`twosite_obs.cpp` の CTM 分岐の窓の組み立て(`indices`、`C_`、`eT*_`)をそのまま使い、`Tn_` だけを `build_relay_window` のチャネルごとの結果に置き換えて `core::Contract` を呼ぶ。`op_` は全セル恒等。ノルムは全セル `build_reduced_density_tensors` 相当の窓で計算する。ops 形式は `product_twosite_op` で明示形式に変換してから、既存の分岐(最近接なら従来の経路)に流す。
- [ ] **Step 5: 不改変検査と独立検証**。`test_fermion_longrange`、`test_fermion_layer`、`test_fermion_fold_full_update` → `ctest --preset gcc` 全件 → MPI ビルドで `ctest` 全件。
- [ ] **Step 6: タスクレビュー → 整形 → コミット**。

---

### Task 3: CTM での相関関数

**Files:**
- Modify: `src/iTPS/correlation_function.cpp`(fermion 分岐)、`src/iTPS/iTPS.hpp`
- Modify: `src/iTPS/load_toml.cpp`, `src/iTPS/measure.cpp`(`r_max > 0` を CTM で受理)
- Modify: `test/fermion/longrange_measure.cpp`(テスト作成者が追記)

**Interfaces:**
- Consumes: `build_relay_site`、`relay_channels`、`product_twosite_op`、`operator_parity`、`onesite_parity`。
- Produces: `measure_correlation()` が fermion の CTM でボソンと同じ `std::vector<Correlation>` を返す。MF は T5 まで拒否。

#### 振る舞い契約書(T3)

1. **入力**: fermion の CTM で `r_max > 0` を受理する。`meanfield_env = true` との組み合わせは拒否する(T5 で解除)。`correlation.operators` の各組は、定義済みの一サイト演算子で、パリティが確定していること。
2. **出力の形**: ボソンと同じく、各左端サイト、各演算子の組、r = 1..r_max、水平(+x)と垂直(+y)について 1 行ずつ。値は ⟨A_s B_t⟩(A が左または下)。
3. **パリティが異なる組**: 値は実部・虚部とも厳密に 0.0。行は出力する。
4. **窓測定との一致**: T2 の 5 と同じ `iTPS`(d = 2 と d = 4)で、r = 1, 2, 3 の相関関数の値が、同じ演算子の組から `product_twosite_op` で作った二サイト物理量を (r, 0) と (0, r) で `measure_twosite()` した値に、相対誤差 1e-10 以下で一致すること。組は (c†, c)、(c, c†)、(n, n) を含める。d = 4 では (c†_↑, c_↑)、(c†_↓, c_↓)、(Sz, Sz) を含める。
5. **r = 1 の一致**: r = 1 の値は既存の最近接経路の値にも一致する(4 から従う。明示的に 1 件確かめる)。
6. **Review Focus 1**: 2×2 単位胞で r_max = 5(単位胞を 2 周以上する)にして、r = 2, 3 が 4 と同じく窓測定に一致すること。
7. **相関長**: fermion では従来どおり無効化され、警告が出る。
8. **MPI**: この実行ファイルの相関関数のケースを MPI n = 2 でも登録し、n = 1 と同じ値(相対誤差 1e-12)になることを確かめる。ランクごとに固定の乱数で状態を作れないフィクスチャなら、決定的テンソルを使う。

#### 手順(T3)

- [ ] **Step 1: テスト作成者をディスパッチ**(スタブは不要。`measure_correlation` は既存で、現状は入力で拒否される)。報告 `work/fermion-longrange/t3/test-author-report.md`。
- [ ] **Step 2: RED 確認とスナップショット**。現状ではガードで落ちることを確認する。
- [ ] **Step 3: Codex に実装させる**。実装方針: 左端の演算子 A と右端の演算子 B の組ごとに `product_twosite_op` → `wrap_twosite_gate` → `relay_channels` を作る。チャネルごとに、Start に source(`build_relay_site`、exit = 右または上)、Transfer に middle(entry = 左、exit = 右。垂直は entry = 下、exit = 上)、Finish に target を渡す。垂直は紐付きテンソルを物理的な向きで作ってから、ボソンと同じ `Axes(3,0,1,2,4,5)` で回す。ノルムの鎖は従来の reduced テンソルで作る。パリティが異なる組は縮約せずに 0 を出力する。
- [ ] **Step 4: 不改変検査と独立検証**(T2 と同じ範囲、MPI を含む)。
- [ ] **Step 5: タスクレビュー → 整形 → コミット**。

---

### Task 4: ツール、E2E、golden、ドキュメント

**Files:**
- Modify: `tool/tenes_simple.py`(`SpinlessFermionModel`、`HubbardModel`、`_check_fermion_scope`、correlation の既定値)、`tool/tenes_std.py`(`_validate_fermion_mode_input`)
- Modify: `test/python/test_tenes_simple.py`, `test_tenes_std.py`, `test_fermion_models.py`(テスト作成者)
- Create: `test/fermion/free_fermion_longrange.py.in`, `test/fermion/boson_equivalence_longrange.py.in`, `test/fermion/hubbard_longrange.py.in` と対応する入力(テスト作成者)
- Modify: `test/CMakeLists.txt`(E2E 登録、1 件は MPI n = 2)
- Regenerate: 一サイト演算子の追加で行が増える既存の golden(`test/data/output_FreeFermion*/` ほか。どれが該当するかは Step 3 で機械的に洗い出す)
- Modify: `docs/sphinx/{ja,en}/file_specification/simple_format.rst`, `parameter_section.rst`, `output_format.rst`, `docs/sphinx/en/algorithm/algorithms.rst`(と ja の対応ファイル)、`NEWS.md`、`sample/07_spinless_fermion/`、`sample/08_fermion_hubbard/`

#### 振る舞い契約書(T4)

1. **`tenes_simple`**:
   - 一サイト演算子の名前と番号: spinless は 0: n, 1: cdag, 2: c。Hubbard は 0: n, 1: n_up, 2: n_dn, 3: Sz, 4: doublon, 5: holon, 6: cdag_up, 7: c_up, 8: cdag_dn, 9: c_dn。
   - cdag, c の行列は、Hubbard のサイト内順序 |up dn⟩ = c†_up c†_dn |0⟩(`local_index_to_occupation`)と整合する(c_dn は c_up を跨ぐ符号を持つ)。
   - fermion モデルで `[correlation]` を受理する。`operators` 省略時の既定値は、spinless `[[0, 0], [1, 2]]`、Hubbard `[[0, 0], [1, 1], [2, 2], [3, 3], [4, 4], [5, 5], [6, 7], [8, 9]]`。
   - `[correlation_length]` は従来どおり拒否。最近接以外の Hamiltonian 項も従来どおり拒否。
2. **`tenes_std`**: fermion モードで ops 形式の二サイト物理量を受理する。multisite と最近接以外の Hamiltonian 項は従来どおり拒否。
3. **E2E: 自由フェルミオンの解析値**(spinless、half filling):
   - (2,0)、(1,1)、(2,1)、(3,0) の ⟨c†_0 c_r⟩ + h.c. と、r_max = 6 の相関関数 ⟨c†_0 c_r⟩ を、運動量積分による解析値 ⟨c†_0 c_r⟩ = ∫ d²k/(2π)² e^{ik·r} n(k) と比べる。
   - 許容誤差は、相殺しない量に錨を置いて決める(`docs/superpowers/specs/2026-09-05-fermion-test-tolerance-contract.md`)。有限 D の誤差を実測してから決め、根拠をテストのコメントに書く。
   - golden との比較(`fulltest.py` と同じ rtol・atol)も行う。
4. **E2E: ボソン等価性**: パリティ表をすべて偶にした fermion の実行と、同じ入力のボソンの実行で、長距離二サイト物理量((2,0), (1,1))と相関関数が一致する(`test/fermion/boson_equivalence_full.py.in` と同じ型・同じ許容誤差の考え方)。
5. **E2E: Hubbard**: 小さい D で、⟨c†_{↑,0} c_{↑,r}⟩ と ⟨c†_{↓,0} c_{↓,r}⟩ が一致する(スピン対称な初期条件と Hamiltonian)。SzSz の相関関数が有限値である。
6. **MPI**: 3 の E2E を MPI n = 2 でも登録し、n = 1 と一致すること。
7. **golden の再生成**: 既存の fermion の E2E で一サイト演算子の行が増えるものは再生成する。**既存の行の値が変わらないこと**を、再生成前後の差分(追加行を除く)が空であることで確かめ、その確認手順を報告に書く。
8. **ドキュメント**: ja と en の両方で、次の内容を反映する。長距離二サイト物理量(窓 4×4 まで)、ops 形式、奇の一サイト演算子(一サイト測定値は 0)、`[correlation]`(値は ⟨A_s B_t⟩、パリティの異なる組は 0)、MF でも使えること(T5 のあとに追記)、相関長と multisite が未対応であること。NEWS は develop に対する差分だけを書く(ブランチ内の経緯は書かない)。

#### 手順(T4)

- [ ] **Step 1: テスト作成者をディスパッチ**(1・2・3・4・5・6 のテスト)。報告 `work/fermion-longrange/t4/test-author-report.md`。
- [ ] **Step 2: RED 確認とスナップショット**。
- [ ] **Step 3: 再生成が必要な golden の洗い出し**(Claude)。fermion の E2E の入力が `tenes_simple` を経由するかを調べ、該当する golden を列挙して台帳に書く。
- [ ] **Step 4: Codex に実装させる**(ツール、golden の再生成、サンプル、ドキュメント、NEWS)。golden は E2E と同じ手順で生成し、7 の差分確認を報告に含めさせる。
- [ ] **Step 5: 不改変検査と独立検証**。`pytest test/python`、`ctest --preset gcc` 全件、MPI の `ctest` 全件、ドキュメントのビルド(`docs/sphinx` の ja・en)。7 の差分確認を Claude が自分でやり直す。
- [ ] **Step 6: タスクレビュー → 整形(black を含む)→ コミット**。

---

### Task 5: MF 環境

**Files:**
- Modify: `src/iTPS/twosite_obs.cpp`, `src/iTPS/correlation_function.cpp`(fermion MF の長距離)
- Modify: `src/fermion/relay.hpp`(χ = 1 の合成環境の生成関数)
- Modify: `src/iTPS/load_toml.cpp`, `src/iTPS/measure.cpp`(MF の拒否を解除)
- Create: `test/fermion/longrange_mf.cpp`(テスト作成者)
- Modify: `test/CMakeLists.txt`、E2E(`free_fermion_longrange.py.in` に MF の場合を追加。テスト作成者)、ドキュメント

**Interfaces:**
- Produces:

```cpp
namespace tenes::fermion {

//! CHI = 1 corner of the delta environment: shape (1, 1), value 1.
template <class tensor>
tensor make_delta_corner(MPI_Comm comm);

//! CHI = 1 edge closing one outer fused leg [x xb] (x fastest, as
//! doubled_pipeline fuses it) with delta_{x, xb}: shape (1, 1, D * D), the
//! density-CTM edge layout (src/iTPS/tensors.cpp). D is the virtual
//! dimension of that leg; legs of one window may differ.
template <class tensor>
tensor make_delta_edge(int D, MPI_Comm comm);

}  // namespace tenes::fermion
```

#### 振る舞い契約書(T5)

1. **入力**: `meanfield_env = true` で、最近接以外の二サイト物理量、ops 形式、`r_max > 0` を受理する。
2. **λ の掛け方**: 設計書 §6 のとおり。窓測定では外周のサイトの外周側の脚に、相関関数(水平)では鎖の全サイトの t・b と、source の l、target の r に掛ける(垂直は対応する脚)。λ 済みのテンソルから ket と bra の両層を作る。
3. **ゲート 1(最近接の一致)**: 最近接の (1,0)、(0,1)、(−1,0)、(0,−1) で、この経路(λ 済みの窓 + δ 環境 + density カーネル)の値が、既存の `contract_pair_MF` の値と相対誤差 1e-12 以下で一致すること。d = 2 と d = 4、実数と複素数。一致しなければ停止して設計を改訂する。
4. **ゲート 2(全偶パリティのボソン一致)**: パリティ表をすべて偶にした状態で、長距離の窓測定 (2,0)・(1,1)・(3,0) と相関関数 r = 1..3 が、ボソンの MF 経路(`Contract_iTPS_MF`、`StartCorrelation_iTPS_MF` 系)と相対誤差 1e-12 以下で一致すること。
5. **窓と相関関数の一致**: MF でも、相関関数 r = 1..3 が同じ組の窓測定に相対誤差 1e-10 以下で一致すること(T3 の 4 の MF 版)。
6. **E2E**: 自由フェルミオンの長距離測定を MF でも走らせ、golden と比べる(MF は近似なので、解析値との比較はしない)。
7. **ドキュメント**: MF でも使えることを ja・en に追記する。

#### 手順(T5)

- [ ] **Step 1: スタブを置く**(`make_delta_corner`、`make_delta_edge`)。
- [ ] **Step 2: テスト作成者をディスパッチ** → **Step 3: RED 確認とスナップショット** → **Step 4: Codex 実装** → **Step 5: ゲート 1・2 の判定を台帳に記録**(不合格なら停止)→ **Step 6: 独立検証**(ctest 全件、MPI)→ **Step 7: タスクレビュー → 整形 → コミット**。

---

## 最終段

- [ ] **全ブランチレビュー**: 最上位モデルのレビュアーで 1 回。どのテストも通らない経路(例: 最近接以外の MF で ops 形式が使われる経路、相関関数で A = B が奇の場合)を探させる。
- [ ] **PR の準備**: PR 本文は develop に対する差分だけを書く。「Generated with Claude Code and Codex」と明記する。`work/fermion-longrange/pr/PR-BODY.md` に下書きを置く。
