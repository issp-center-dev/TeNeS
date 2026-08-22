# フェルミオン二体測定の定数倍高速化 実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
>
> **本リポジトリ固有の運用**: `~/.claude/CLAUDE.md` の多段エージェント方針に従う。
> **テストは専用サブエージェントが「契約」(散文)だけから書く。実装者(Codex)はテストファイルを変更してはならない。**
> 計画者(Claude)はテストコードを書かず、契約のみを渡す。

**Goal:** フェルミオン二体測定の rank-16 符号スイープを 17 パスから 1 パスに畳み、除算を減らし並列化して、数値結果を変えずに D=4 の observable を 1857 s から数十秒に落とす。

**Architecture:** 符号適用を「パリティビットマスク上の二次形式 + 対角因子」として一般化した単一スイープカーネル(`sign_sweep.hpp`)を新設する。索引デコードは `global_index_fast`(軸ごとに libc `div`)ではなく `make_l2g_map` + `global_index_l2g_map`(除算 2 回 + 表引き)を使い、OpenMP で並列化する。あわせて bra-bra swap を外積の前の小テンソルへ移し、rank-16 の冗長コピー 2 本を消す。

**Tech Stack:** C++17、mptensor (fork `yomichi/mptensor` @ `986e21f`)、doctest、OpenMP、CMake/ctest

**Spec:** `docs/superpowers/specs/2026-08-22-fermion-twosite-perf-design.md`

## Global Constraints

以下は全タスクの要件に暗黙に含まれる。

- **数値結果を変えない。** 契約は代数的同値。検証としては、有限値の入力に対して**全出力が bit-identical** であること。差が出たら「許容誤差内だからよい」ではなく実装バグの兆候として扱う(spec §5)。
- **`deps/mptensor` を変更しない。** 公開 API (`prep_local_to_global` / `make_l2g_map` / `global_index_l2g_map` / `local_size` / `operator[]`) の範囲で済ませる(spec §3.3)。
- **OpenMP 領域に入る前に `prep_local_to_global()` と `make_l2g_map()` を呼ぶ。** ScaLAPACK バックエンドの遅延初期化と競合する。`_NO_MPI` ビルドでは lapack 版が no-op なので顕在化せず、**MPI ビルドでのみ壊れる**(spec §3.1-2)。
- **`1 << rank` を評価する前に rank を判定する。** 順序を逆にすると rank ≥ 32 でシフト UB(spec §3.1-5)。`kMaxTableRank = 24`。
- **対角因子を事前に掛け合わせない。** 同一スイープ内で個別に乗算する(spec §3.1-7)。
- **ボゾン経路の挙動を変えない。**
- OpenMP は `#ifndef _NO_OMP` で保護する。`_NO_OMP` ビルドでもコンパイル・ctest が通ること。
- 名前空間は `tenes::fermion`。新規ファイルには GPL v3 のライセンスヘッダを既存ファイルからそのまま複製する。
- 整形は `clang-format`(`.clang-format`)。**実装者は formatter を実行しない**。新しい行は周囲のスタイルに手で合わせる。整形はコミット直前に統括者がそのタスクで変更したファイルだけに適用する。
- 実装者はリポジトリ直下でベンチマークを実行しない。実験は `work/<テーマ>/` 配下で行う。
- **実装者はテストファイル(`test/` 配下)を変更してはならない。** テストが誤っていると考えた場合は実装を曲げず BLOCKED として報告する。

## 参照ベースライン

- タグ `wip/fermion_20260822` が改訂前の実装を指す。層2 の参照実装はここから取る:
  `git show wip/fermion_20260822:src/fermion/reduced.hpp`、`:src/fermion/ftensor.hpp`
- ベンチマークのベースライン出力と時間は `work/fermion-perf/baseline/` にある(`work/` は `.gitignore` 済み)。

## ファイル構成

| ファイル | 責務 |
|---|---|
| `src/fermion/sign_sweep.hpp` (新規) | `SwapForm`、`LegGauge`、`SignEval`、単一スイープカーネル、4 つの入口 |
| `src/fermion/ftensor.hpp` (変更) | `ftensor::transpose` をカーネル利用に置換 |
| `src/fermion/reduced.hpp` (変更) | `doubled_pipeline` / `fuse_doubled_cluster` / `build_reduced_pair` / `apply_fused_leg_gauge` |
| `src/fermion/fops.hpp` (変更) | `apply_swap` はカーネル委譲に(公開 API は据え置き) |
| `test/fermion/sign_sweep.cpp` (新規) | 本計画で追加する doctest。`test/test_fermion_layer.cpp` 末尾から `#include` する(既存の `r2_convention.cpp` 等と同じ方式。CMake 変更は不要) |

`test/fermion/*.cpp` は `test/test_fermion_layer.cpp:3855-3857` で `#include` されている。新規ファイルもそこに 1 行足す。**この 1 行の追加はテスト作成者が行う。**

---

## タスク 1: 符号スイープカーネル

**Files:**
- Create: `src/fermion/sign_sweep.hpp`
- Create: `test/fermion/sign_sweep.cpp`
- Modify: `test/test_fermion_layer.cpp:3857` 付近(`#include "fermion/sign_sweep.cpp"` を追加)
- Modify: `src/fermion/fops.hpp:208-219` (`apply_swap` をカーネル委譲に)

**Interfaces — Produces:**

```cpp
namespace tenes::fermion {

// 符号評価の方式。テストから両経路を強制できるようにするための引数であり、
// automatic 以外を本番コードから渡さない。
enum class SignEval { automatic, table, direct };

// パリティ二次形式:  sign = (-1)^{ Σ_{(x,y) ∈ terms} p_x(idx[x]) · p_y(idx[y]) }
class SwapForm {
 public:
  // 非順序ペア {ax1, ax2} をトグルする。同じペアを偶数回呼ぶと消える。
  // ax1 == ax2 は何もしない。
  void toggle(int ax1, int ax2);
  bool empty() const;
  // 正規形 (first < second)、first,second の昇順にソート済み
  const std::vector<std::pair<int, int>>& terms() const;
};

// 1 本の軸に掛かる対角因子
struct LegGauge {
  int axis;
  std::vector<double> factor;
};

// 1 スイープ: element *= (-1)^{form} × Π_k gauges[k].factor[idx[gauges[k].axis]]
template <class tensor>
void apply_sign_sweep(ftensor<tensor>& a, const SwapForm& form,
                      const std::vector<LegGauge>& gauges,
                      SignEval eval = SignEval::automatic);

template <class tensor>
void apply_swap_form(ftensor<tensor>& a, const SwapForm& form,
                     SignEval eval = SignEval::automatic);

// パリティ台帳を持たない素の tensor 用(fused blob の gauge に使う)
template <class tensor>
void apply_leg_gauges(tensor& a, const std::vector<LegGauge>& gauges);

namespace detail {
// テーブルを使うかどうかの判定。rank を先に見ること。
// ヘッダ内の非テンプレート関数なので inline を付けること。
// テストから直接呼ぶので detail であっても実体を公開ヘッダに置く。
inline bool use_sign_table(std::size_t rank, std::size_t local_size);
constexpr std::size_t kMaxTableRank = 24;
}
}  // namespace tenes::fermion
```

**実装の骨子**(この構造に従うこと。埋めるのは実装者)

```cpp
template <class tensor>
void apply_sign_sweep(ftensor<tensor>& a, const SwapForm& form,
                      const std::vector<LegGauge>& gauges, SignEval eval) {
  const std::size_t rank = a.parity.size();
  const std::size_t n_local = a.t.local_size();
  // 軸ごとのビット寄与 bits[ax][i] = parity[ax][i] ? (1u << ax) : 0
  // eval に応じて table(サイズ 1<<rank の int8_t) を作るか作らないかを決める
  a.t.prep_local_to_global();   // ← OpenMP 領域の外。必須
  a.t.make_l2g_map();           // ← 同上。axes_map の現状を反映した表を作る
#ifndef _NO_OMP
#pragma omp parallel default(shared)
#endif
  {
    std::vector<std::size_t> idx(rank);
#ifndef _NO_OMP
#pragma omp for
#endif
    for (std::size_t n = 0; n < n_local; ++n) {
      a.t.global_index_l2g_map(n, idx.data());
      std::uint32_t m = 0;
      for (std::size_t ax = 0; ax < rank; ++ax) m |= bits[ax][idx[ax]];
      if (sign_of(m) < 0) a.t[n] = -a.t[n];
      for (const LegGauge& g : gauges) a.t[n] *= g.factor[idx[g.axis]];
    }
  }
}
```

`apply_swap(a, ax1, ax2)` は `SwapForm` を 1 ペアだけトグルして `apply_swap_form` に委譲する。既存の呼び出し元(`src/fermion/fops.hpp:235,252`、`test/` 3 ファイル)はシグネチャを変えずに通ること。

### 契約(テスト作成者へ渡す散文。テストコードは渡さない)

新しい符号スイープが満たすべき振る舞い。すべて **`==` による厳密一致**で確かめること(許容誤差を使ってはならない)。テンソルの中身は乱数で埋め、parity 台帳は偶奇の内訳を軸ごとに変えた非対称なものを使うこと。実テンソル・複素テンソルの両方を通すこと。

- **C1-1 (トグルの相殺)**: `SwapForm` に同じ非順序ペアを偶数回トグルすると、そのペアは `terms()` に現れない。奇数回なら 1 度だけ現れ、`first < second` に正規化されている。`(x, y)` と `(y, x)` は同じペアとして扱われる。`toggle(x, x)` は `terms()` を変えない。
- **C1-2 (逐次適用との一致)**: 重複のないペア列を与えたとき、`apply_swap_form` の結果は、同じペア列に対して `apply_swap` を 1 本ずつ順に適用した結果と全要素一致する。ペアの与える順序を入れ替えても結果は同じであること。
- **C1-3 (重複ペアの扱い)**: `(0,2)` と `(2,0)` の両方を含むペア列を与えたとき、`apply_swap_form` の結果は、その 2 本を取り除いた列を逐次 `apply_swap` した結果と全要素一致する(すなわち相殺する)。**これは現行 `apply_joint_swaps` が実際に生成する形なので、必ず 1 ケース以上置くこと。**
- **C1-4 (空の form)**: 空の `SwapForm` を渡してもテンソルの全要素が変化しない。
- **C1-5 (索引デコードの同値)**: rank 2 から 16 までの複数の形状(各軸の次元も揃えないこと)について、`make_l2g_map()` を呼んだ後の `global_index_l2g_map(n, idx)` が返す添字が、同じ `n` に対する `global_index_fast(n, idx)` の返す添字と**全ての局所インデックスで一致**する。これは本タスクの高速化の前提そのものなので、独立したテストケースとして置くこと。
- **C1-6 (評価経路の同値)**: 同一の入力・同一の form に対し、`SignEval::table` と `SignEval::direct` が全要素一致する。加えて `SignEval::automatic` がそのどちらかと一致する。
- **C1-7 (判定関数の境界)**: `detail::use_sign_table(rank, local_size)` は、`rank` が `kMaxTableRank` を超えるとき `local_size` の値によらず false を返す。`rank` が `kMaxTableRank` 以下のときは、`2^rank` が `local_size` を超えるなら false、超えないなら true を返す。**rank = 32、rank = 64 のような大きな値を渡しても未定義動作を起こさず false を返すこと。**
- **C1-8 (対角因子)**: `LegGauge` を複数本渡したとき、結果は同じ因子を `multiply_vector` で 1 本ずつ順に掛けた結果と全要素一致する。因子が ±1 でない値(例えば 2.0 と 0.5)のケースも 1 つ置くこと。
- **C1-9 (既存 API の非回帰)**: `apply_swap(a, ax1, ax2)` は、パリティが奇である軸の組み合わせの要素だけを符号反転する。`ax1` と `ax2` を入れ替えても結果は同じ。

### 手順

- [ ] **Step 1: 契約をテスト作成サブエージェントに渡す**

`test/fermion/sign_sweep.cpp` を新規作成し、`test/test_fermion_layer.cpp` 末尾の `#include` 群に 1 行足させる。渡すのは上の「契約」節の散文のみ。既存の `test/fermion/r2_convention.cpp` を書式の手本として示す。

- [ ] **Step 2: RED を 1 件ずつ確認する(統括者)**

Run: `cmake --build --preset gcc --target test_fermion_layer && ./out-gcc/build/test/test_fermion_layer`
Expected: 新規ケースが**コンパイルエラーではなく**アサーション失敗で落ちる、あるいは未定義シンボルで落ちる。
**各ケースについて「正しい理由で落ちているか」を目視で確認する。** 期待値を書き間違えて落ちているだけのケースはテスト作成者に差し戻す。

- [ ] **Step 3: テストファイルのバイト単位スナップショットを取る(統括者)**

```bash
mkdir -p work/twosite-perf
git hash-object test/fermion/sign_sweep.cpp test/test_fermion_layer.cpp \
  > work/twosite-perf/task1-test-snapshot.txt
```

- [ ] **Step 4: 実装を Codex に投げる**

`src/fermion/sign_sweep.hpp` を新設し、`src/fermion/fops.hpp` の `apply_swap` を委譲に変える。上の Interfaces と実装骨子、Global Constraints を渡す。**テストファイルの変更は禁止**と明記する。報告ファイル `work/twosite-perf/task1-report.md` は成果物であり、無いこと自体が欠陥であると明記する。

- [ ] **Step 5: テスト不改変を機械検査する(統括者)**

```bash
git hash-object test/fermion/sign_sweep.cpp test/test_fermion_layer.cpp \
  | diff - work/twosite-perf/task1-test-snapshot.txt && echo "TESTS UNCHANGED"
```

- [ ] **Step 6: 統括者が独立に検証する**

Run: `cmake --build --preset gcc && ctest --preset gcc`
Expected: 全件 PASS。Codex の報告を鵜呑みにせず自分で回す。

- [ ] **Step 7: `_NO_OMP` 構成でのビルドと ctest**

Run: `cmake -S . -B out-noomp -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-16 -DTENES_PYTHON_EXECUTABLE=$PWD/venv/bin/python3 -DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=ON && cmake --build out-noomp -j && (cd out-noomp && ctest)`
Expected: 全件 PASS。

- [ ] **Step 8: タスクレビュー(新規サブエージェント)**

レビュアーには **作業ツリーの変更を禁止する(`git stash` / `git checkout` を含む)**。仮説の検証はファイルを `$CLAUDE_JOB_DIR/tmp` にコピーして行わせる。変異テストを 1 件以上依頼する: `SwapForm` の XOR 畳み込みを「重複を無視して追加」に変えたコピーで C1-3 が赤くなること。

- [ ] **Step 9: 整形してコミット(統括者)**

```bash
clang-format -i src/fermion/sign_sweep.hpp src/fermion/fops.hpp
git add src/fermion/sign_sweep.hpp src/fermion/fops.hpp \
        test/fermion/sign_sweep.cpp test/test_fermion_layer.cpp
git commit -m "Add the single-sweep parity sign kernel for fermionic tensors"
```

---

## タスク 2: `ftensor::transpose` をカーネルに載せ替える

**Files:**
- Modify: `src/fermion/ftensor.hpp:96-108` (`ftensor::transpose`)
- Modify: `src/fermion/sign_sweep.hpp` (`transpose_with_swap_form` を追加)
- Modify: `src/fermion/fops.hpp:80-91` (`apply_transpose_sign_mask`、恒等置換の短絡)
- Test: `test/fermion/sign_sweep.cpp` (追記)

**Interfaces — Consumes:** タスク 1 の `SwapForm`、`SignEval`、`apply_sign_sweep`
**Interfaces — Produces:**

```cpp
namespace tenes::fermion {
// 1 スイープで form の符号と `axes` に対する graded transpose 符号を掛け、
// そのあと a.t を `axes` で(遅延)置換し、a.parity も並べ替える。
template <class tensor>
void transpose_with_swap_form(ftensor<tensor>& a, const SwapForm& form,
                              const mptensor::Axes& axes,
                              SignEval eval = SignEval::automatic);
}
```

graded transpose 符号は、現行 `src/fermion/ftensor.hpp:19-33` の `transpose_sign`
(奇パリティ脚だけに制限した転倒数)と同じ量である。これをビットマスクの関数として
テーブル化する。`apply_transpose_sign_mask` は置換が恒等のとき常に +1 を返すので、
その場合はスイープ自体を省く。

### 契約(テスト作成者へ渡す散文)

- **C2-1 (旧 transpose との一致)**: 空の `SwapForm` を渡した `transpose_with_swap_form` の結果は、タグ `wip/fermion_20260822` 時点の `ftensor::transpose` の実装(`git show wip/fermion_20260822:src/fermion/ftensor.hpp` から `reference` 名前空間にそのまま写して使う)と、テンソルの全要素および `parity` の並びの両方で一致する。rank 2 から 10 までの複数形状、恒等置換・逆順・巡回を含む複数の `axes` で確かめること。
- **C2-2 (swap との合成)**: 空でない `SwapForm` と `axes` を同時に渡した結果は、「先に `apply_swap_form` を適用してから参照実装の `transpose` を適用」した結果と全要素一致する。
- **C2-3 (評価経路の同値)**: `SignEval::table` と `SignEval::direct` が一致する。
- **C2-4 (恒等置換の短絡)**: 恒等置換を渡したとき、テンソルの全要素と `parity` が変化しない。
- **C2-5 (書き換え後の `ftensor::transpose`)**: 書き換え後の `ftensor::transpose` が C2-1 の参照実装と全要素一致する。**`fermion::transpose` / `tensordot` / `svd` / `qr` は `ftensor::transpose` を経由するので、既存のこれらのテストが緑のままであることも合わせて確認すること。**

### 手順

- [ ] **Step 1: 契約をテスト作成サブエージェントに渡す**(`test/fermion/sign_sweep.cpp` に追記)
- [ ] **Step 2: RED を 1 件ずつ確認する(統括者)**
- [ ] **Step 3: テストファイルのスナップショット**(`work/twosite-perf/task2-test-snapshot.txt`)
- [ ] **Step 4: 実装を Codex に投げる**(報告は `work/twosite-perf/task2-report.md`)
- [ ] **Step 5: テスト不改変の機械検査**
- [ ] **Step 6: 統括者が `ctest --preset gcc` を全件回す。**`ftensor::transpose` は `svd`/`qr`/`tensordot` の共有コードであり simple update も通る。**フェルミオン関連だけに絞ってはならない。**
- [ ] **Step 7: `_NO_OMP` 構成でのビルドと ctest**
- [ ] **Step 8: タスクレビュー(新規サブエージェント)。** 変異テスト: `transpose_sign` のテーブル化で「奇パリティ脚に制限する条件」を落としたコピーで C2-1 が赤くなること。
- [ ] **Step 9: 整形してコミット**

```bash
clang-format -i src/fermion/ftensor.hpp src/fermion/sign_sweep.hpp src/fermion/fops.hpp
git commit -m "Sweep the graded transpose sign in one pass"
```

---

## タスク 3: 二重化パイプラインの再構成

ここが本命。rank-16 テンソルを触るパスを 17 本から 1 本に減らし、冗長コピーを 2 本消す。

**Files:**
- Modify: `src/fermion/reduced.hpp:52-71` (`apply_joint_swaps` → form を返す関数へ)
- Modify: `src/fermion/reduced.hpp:83-108` (`doubled_pipeline`)
- Modify: `src/fermion/reduced.hpp:111-140` (`fuse_doubled_cluster`)
- Modify: `src/fermion/reduced.hpp:196-235` (`build_reduced_pair`)
- Test: `test/fermion/sign_sweep.cpp` (追記)

**Interfaces — Consumes:** タスク 1・2 の全て
**Interfaces — Produces:**

```cpp
namespace tenes::fermion::detail {
// apply_joint_swaps が適用していた swap 群を、適用先で二分して返す。
//   cross : ket 側と bra 側にまたがるペア。外積後の doubled テンソルに適用する。
//   bra   : bra 側だけで閉じるペア。外積の前に第一因子へ適用する。
// どちらも SwapForm なので重複ペアは自動的に相殺される。
struct JointSwapForms {
  SwapForm cross;
  SwapForm bra;
};
JointSwapForms joint_swap_forms(const std::vector<int>& bra_axes,
                                const std::vector<int>& ket_axes,
                                const std::vector<int>& leg_ids);

// 変更後のシグネチャ。外積を関数の内側で行う。
template <class tensor>
tensor fuse_doubled_cluster(const ftensor<tensor>& bra_pair,
                            const ftensor<tensor>& ket_pair,
                            const std::vector<int>& leg_ids);
}
```

**実装の骨子**

`joint_swap_forms` は現行 `apply_joint_swaps`(`src/fermion/reduced.hpp:52-71`)の
二重ループをそのまま辿り、`apply_swap(a, ket_axes[ix], bra_axes[iy])` に対応するペアを
`cross` へ、`apply_swap(a, bra_axes[ix], bra_axes[iy])` に対応するペアを `bra` へ
トグルする。**判定条件(`kDoubledJointMask`、`joint_bit`)は一切変えない。**

`fuse_doubled_cluster` は:

```cpp
const std::vector<int> cluster_axes = {0, 1, 2, 4, 5, 6};
// bra_axes = cluster_axes, ket_axes = cluster_axes の各要素 + 8
const auto forms = joint_swap_forms(bra_axes, ket_axes, leg_ids);
ftensor<tensor> bra = bra_pair;
apply_swap_form(bra, forms.bra);                        // rank 8、16k 要素
ftensor<tensor> doubled =
    tensordot(bra, ket_pair, mptensor::Axes(), mptensor::Axes());  // rank 16
transpose_with_swap_form(doubled, forms.cross, interleaved);       // ← 唯一の rank-16 スイープ
tensor fused = mptensor::reshape(doubled.t, sh);
return mptensor::contract(fused, mptensor::Axes(6, 7), mptensor::Axes(8, 9));
```

`build_reduced_pair` の呼び出しは
`detail::fuse_doubled_cluster(conj(ket_ab), ket_op, leg_ids)` になる。
**`conj()` は今までどおり先に適用する**(bra-bra form は conj の後の因子に掛かる)。
`doubled_pipeline` も同じ形に直す(`bra_axes = {0,1,2,3}`、`ket_axes = {5,6,7,8}`、
`leg_ids = {0,1,2,3}`、bra form は rank-5 の `conj(bra_Tn)` に適用)。

**注意(spec §3.2)**: bra-bra swap を外積の前に移せるのは、その `tensordot` の縮約軸が
**空**だからである。graded `tensordot` は `apply_transpose_sign_mask` で被演算子に符号を
掛けるが、`Axes(), Axes()` では左右どちらの置換も恒等なのでマスクが全て +1 になる。
**縮約軸が空でない `tensordot` に同じ移動を一般化してはならない。**

`ftensor<tensor> prepared = doubled;` と `transpose(prepared, interleaved)` の
値返しコピーはこの書き換えで消える。**新たなコピーを増やさないこと**(D=4 で 1 本 2.15 GB)。

### 契約(テスト作成者へ渡す散文)

参照実装は、タグ `wip/fermion_20260822` 時点の `src/fermion/reduced.hpp`
(`git show wip/fermion_20260822:src/fermion/reduced.hpp`)から
`apply_joint_swaps` / `doubled_pipeline` / `fuse_doubled_cluster` /
`build_reduced_op` / `build_reduced_pair` を `reference` 名前空間へそのまま写して使う。

- **C3-1 (`build_reduced_op` の同値)**: 乱数 `Tn` について、新実装の `build_reduced_op` の出力が参照実装の出力と**全要素 `==` 一致**する。形状(rank と各軸の次元)も一致する。物理次元 d=2 と d=4、仮想次元 D=2 と D=3、実テンソルと複素テンソルの組み合わせを網羅すること。
- **C3-2 (`build_reduced_pair` の同値、horizontal)**: 乱数 `TnA`, `TnB` と乱数の二体演算子について、`reduced_pair_direction::horizontal` で新実装と参照実装の出力が全要素一致する。演算子は少なくとも (a) 恒等、(b) 偶パリティの一般の演算子、(c) 奇×奇の項を含むパリティ保存演算子(ホッピング `c†_A c_B + h.c.` に相当するもの)の 3 種を通すこと。
- **C3-3 (`build_reduced_pair` の同値、vertical)**: C3-2 と同じことを `vertical` で行う。**`leg_ids` が方向ごとに違うので、片方だけのテストでは取りこぼす。**
- **C3-4 (`joint_swap_forms` の網羅性)**: `joint_swap_forms` が返す `cross` と `bra` を使って参照実装の `apply_swap` を「cross は doubled に、bra も doubled 上の bra 側軸に」適用した結果が、参照実装の `apply_joint_swaps` を doubled に適用した結果と全要素一致する。すなわち**二分と相殺で符号を落としていない**こと。`doubled_pipeline` 用の引数(rank 10)と `fuse_doubled_cluster` 用の引数(rank 16)の両方、horizontal と vertical の両方で確かめること。
- **C3-5 (`build_reduced` / `build_reduced_identity_pair` の非回帰)**: これらは `build_reduced_op` の上に載っているので、参照実装との全要素一致を 1 ケースずつ置くこと。

### 手順

- [ ] **Step 1: 契約をテスト作成サブエージェントに渡す**
- [ ] **Step 2: RED を 1 件ずつ確認する(統括者)**
- [ ] **Step 3: テストファイルのスナップショット**(`work/twosite-perf/task3-test-snapshot.txt`)
- [ ] **Step 4: 実装を Codex に投げる**(報告は `work/twosite-perf/task3-report.md`)
- [ ] **Step 5: テスト不改変の機械検査**
- [ ] **Step 6: 統括者が `ctest --preset gcc` を全件回す**
- [ ] **Step 7: `_NO_OMP` 構成でのビルドと ctest**
- [ ] **Step 8: bit-identical の実地確認(統括者)**

```bash
cd work/fermion-perf
../../out-gcc-release/build/src/tenes input_D2.toml
cmp out_D2/twosite_obs.dat baseline/out_D2/twosite_obs.dat
cmp out_D2/onesite_obs.dat baseline/out_D2/onesite_obs.dat
cmp out_D2/density.dat     baseline/out_D2/density.dat
```

Expected: 3 本とも差分なし。**差が出たらここで止めて原因を特定する。**「小さいから許容」としてはならない。

- [ ] **Step 9: タスクレビュー(新規サブエージェント)。** 変異テスト: `joint_swap_forms` で `bra` 側のトグルを 1 本落としたコピーで C3-4 と C3-2 が赤くなること。
- [ ] **Step 10: 整形してコミット**

```bash
clang-format -i src/fermion/reduced.hpp
git commit -m "Fuse the doubled-cluster sign passes into a single sweep"
```

---

## タスク 4: fused leg gauge の 1 スイープ化

優先度は低い(spec §2.4 のとおり `multiply_vector` は既に OpenMP 並列で、D=4 で 4 秒程度)。タスク 3 までで目標を達していれば省いてよい。**省く場合は spec §7 のタスク 4 に「見送り」と理由を追記すること。**

**Files:**
- Modify: `src/fermion/reduced.hpp:68-81` (`apply_fused_leg_gauge`)
- Modify: `src/fermion/reduced.hpp:220-233` (`build_reduced_pair` の呼び出し 2 箇所)
- Test: `test/fermion/sign_sweep.cpp` (追記)

**Interfaces — Consumes:** タスク 1 の `LegGauge` / `apply_leg_gauges`

現行の `apply_fused_leg_gauge` が組み立てる符号ベクトルの計算はそのまま残し、
`multiply_vector` を 2 回呼ぶ代わりに 2 本の `LegGauge` を作って
`apply_leg_gauges` に一度で渡す形にする。

### 契約(テスト作成者へ渡す散文)

- **C4-1**: `build_reduced_pair` の出力が、タスク 3 完了時点の実装(すなわち gauge を `multiply_vector` で 2 回掛ける版)の出力と全要素 `==` 一致する。horizontal と vertical の両方。**±1 しか掛からないので厳密一致でなければならない。**

### 手順

- [ ] **Step 1: 契約をテスト作成サブエージェントに渡す**(`test/fermion/sign_sweep.cpp` に追記)
- [ ] **Step 2: RED を 1 件ずつ確認する(統括者)**

Run: `cmake --build --preset gcc --target test_fermion_layer && ./out-gcc/build/test/test_fermion_layer`
Expected: 新規ケースがアサーション失敗で落ちる。**正しい理由で落ちているかを目視で確認する。**

- [ ] **Step 3: テストファイルのバイト単位スナップショットを取る(統括者)**

```bash
git hash-object test/fermion/sign_sweep.cpp test/test_fermion_layer.cpp \
  > work/twosite-perf/task4-test-snapshot.txt
```

- [ ] **Step 4: 実装を Codex に投げる**(報告は `work/twosite-perf/task4-report.md`。テストファイル変更禁止を明記)
- [ ] **Step 5: テスト不改変を機械検査する(統括者)**

```bash
git hash-object test/fermion/sign_sweep.cpp test/test_fermion_layer.cpp \
  | diff - work/twosite-perf/task4-test-snapshot.txt && echo "TESTS UNCHANGED"
```

- [ ] **Step 6: 統括者が `ctest --preset gcc` を全件回す**
- [ ] **Step 7: `_NO_OMP` 構成でのビルドと ctest**

Run: `cmake --build out-noomp -j && (cd out-noomp && ctest)`

- [ ] **Step 8: bit-identical の実地確認(統括者)**

```bash
cd work/fermion-perf
../../out-gcc-release/build/src/tenes input_D2.toml
cmp out_D2/twosite_obs.dat baseline/out_D2/twosite_obs.dat
cmp out_D2/onesite_obs.dat baseline/out_D2/onesite_obs.dat
cmp out_D2/density.dat     baseline/out_D2/density.dat
```

Expected: 3 本とも差分なし。

- [ ] **Step 9: タスクレビュー(新規サブエージェント)。** レビュアーには**作業ツリーの変更を禁止する**。変異テスト: 2 本の `LegGauge` のうち片方を落としたコピーで C4-1 が赤くなること。
- [ ] **Step 10: 整形してコミット**

```bash
clang-format -i src/fermion/reduced.hpp
git add src/fermion/reduced.hpp test/fermion/sign_sweep.cpp
git commit -m "Apply both fused-leg gauges in one sweep"
```

---

## タスク 5: 再測定と設計書の実測値差し替え

**Files:**
- Modify: `docs/superpowers/specs/2026-08-22-fermion-twosite-perf-design.md` (§4 と §9)

**このタスクは統括者が自分で行う。** Codex には投げない。

- [ ] **Step 1: D=2/3/4 を再実行して時間と RSS を記録する**

```bash
cd work/fermion-perf
for D in 2 3 4; do
  /usr/bin/time -l ../../out-gcc-release/build/src/tenes input_D$D.toml 2> rss_D$D.txt
  cat out_D$D/time.dat
done
```

- [ ] **Step 2: 全 D で出力のバイト一致を確認する**

```bash
for D in 2 3 4; do
  for f in twosite_obs.dat onesite_obs.dat density.dat; do
    cmp out_D$D/$f baseline/out_D$D/$f || echo "MISMATCH D=$D $f"
  done
done
```

Expected: `MISMATCH` が 1 件も出ない。

- [ ] **Step 3: D=4 の observable フェーズを再プロファイルする**

```bash
cd work/fermion-perf && ../../out-gcc-release/build/src/tenes input_D4.toml > run_D4.log 2>&1 &
until grep -q "twosite operators" run_D4.log; do sleep 5; done
sample $(pgrep -n -f "tenes input_D4.toml") 20 -f /tmp/samp_after.txt
```

新しい支配項を特定する。spec §4 の予想は「`contract_reduced_pair_*_density_CTM` の BLAS 縮約(D¹²·χ² ≈ 1.7e10 flops/回 × 24 回)」である。**予想と違ったら、違ったことを書く。**

- [ ] **Step 4: spec §4 を実測値で置き換え、§9(Phase 2 申し送り)を再プロファイル結果で更新する**
- [ ] **Step 5: コミット**

```bash
git add docs/superpowers/specs/2026-08-22-fermion-twosite-perf-design.md
git commit -m "Record the measured speedup and the new dominant term"
```

---

## タスク 6: 全ブランチの最終レビュー

- [ ] **Step 1: 最上位モデルのサブエージェントに全差分を 1 回でレビューさせる**

タスク単位では構造的に見えない欠陥(どのテストも通らない経路、タスク間で食い違う規約など)は全体レビューだけが捕まえる。レビュアーには**作業ツリーの変更を禁止する**。

観点として最低限これを渡す:
- OpenMP 領域に入る前の `prep_local_to_global()` / `make_l2g_map()` が全経路で呼ばれているか。**`_NO_MPI` ビルドでは検出できない**ことを踏まえて目視で追うこと
- `1 << rank` を評価する前に rank を判定しているか
- rank-16 テンソルのコピーが増えていないか(D=4 で 1 本 2.15 GB)
- bra-bra swap の前倒しが、縮約軸が空でない `tensordot` に波及していないか
- テストが空洞でないか(変異テストで確認)

- [ ] **Step 2: 指摘を反映し、`ctest --preset gcc` を全件回してからコミット**

---

## 完了条件

1. `ctest --preset gcc` が全件 PASS(通常構成)
2. `_NO_OMP` 構成でもビルドと ctest が全件 PASS
3. `work/fermion-perf/` の D=2/3/4 で `twosite_obs.dat` / `onesite_obs.dat` / `density.dat` が `baseline/` と**バイト一致**
4. テストファイルがタスクごとのスナップショットから変更されていないこと(実装者による改変がないこと)
5. spec §4 が実測値に置き換わっていること
6. 全ブランチレビューの指摘が解消していること
