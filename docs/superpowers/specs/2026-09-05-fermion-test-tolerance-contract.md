# フェルミオンテストの許容誤差 — 振る舞い契約(2026-09-05)

原因究明: `work/fermion/ci-tolerance/FINDINGS.md`
発端: push `4b466d10` に対する GitHub Actions で `test_fermion_layer` が 7 ジョブ中 7 で赤。
落ちたのは 2 箇所だけで、どちらも**テスト側の許容誤差の較正**の問題であり、実装の誤りではない。

この契約は既存のテストの**判定基準だけ**を変える。テストが何を主張しているか(covered condition、
反空洞化の前提、ケース集合)は変えない。実装(`src/`)は変えない。

---

## 規約 1: 縮約の順序違いを比べる判定は、相殺しうる結果ではなく計算のスケールに錨を打つ

### 根拠(実測)

`test/fermion/fold_geometry.cpp` の T15、ケース
`d=2 vp=eo dir=v op=hopping chi=3 real op_seed=31 env_seed=140 site_seeds=(61,62)` で:

- C-1(演算子側)の参照値 `-7.81267e-04`、C-2(identity 側)の参照値 `-1.64733e-07`。
- **両者の `|got - ref|` は完全に同一の 1.0842e-19。** 丸め誤差の絶対値は縮約の最大中間量で決まり、
  結果の大きさとは無関係である。C-2 だけ結果が 4700 倍小さいので、結果に比例する許容誤差では
  余裕が 1.5 倍しか残らず、コンパイラを変えると落ちる(CI で実際に 5 ジョブが 3.25e-19、
  macOS g++-13 が 2.71e-19 を記録)。
- 実スケール 7.8e-4 に対して測れば `4e-16` = 2 ulp である。

同じファイルの `fg_check_scalar`(`fold_geometry.cpp:541`)はすでに正しい設計で、
`1e-10 * max(|ref_norm| * op_max, |ref|)` とノルムに錨を打っており実測余裕は 7141 倍ある。
**この設計に揃える。**

### 要求

以下の 3 つのヘルパの許容誤差を、そのケースで計算された**同種の閉じ値のうち最大のもの**に錨を打つ形へ改める。
相対係数(1e-12 等)は変えない。

| # | 場所 | ヘルパ / テスト | 現状の最悪余裕 |
|---|---|---|---|
| 1-a | `test/fermion/fold_geometry.cpp:2111` | `fg15_check_rel`(T15 C-1 / C-2) | **1.52**(CI 失敗) |
| 1-b | `test/fermion/full_update_env.cpp:155` | `fue_check_rel`(T2-i / T2-iv) | 43.2 |
| 1-c | `test/fermion/ctm_phase.cpp:1017` | T8-ii (6)「T2-i note」 | 123 |

具体的には、`tol = rel * max(|got|, |ref|)` を `tol = rel * max(|got|, |ref|, scale)` に改める。
`scale` はそのケースで比較される閉じ値の絶対値の最大(T15 なら `max(|ref_op|, |ref_id|)`)とする。
`scale` は判定される側(`got`)からではなく**参照側から**取ること。実装が壊れて `got` が巨大になったときに
許容誤差まで一緒に膨らむ経路を作らないため。

### 反空洞化(必須)

- 上の変更で緩むのは、同じケース内で相対的に小さい方の閉じ値だけである。大きい方(T15 なら C-1)の
  判定は現状と厳密に同じ値になること。**これをテストで示す**: `scale` を導入する前後で、
  C-1 側の `tol` が変わらないことをケース単位で確認した旨を報告に書く。
- 変異検査: `fg15_check_rel` の `scale` を `1e300` のような無意味に大きい値へ置き換えたコピーで
  T15 が**緑のまま**になることを確認し(= 検出力が消えることを示し)、正しい `scale` では
  既知の変異(たとえば `contract_reduced_pair_halves_density_CTM` の入力脚を 1 本入れ替える)で
  **赤になる**ことを確認すること。

---

## 規約 2: ALS を経由する 2 経路の比較は、ALS が保証している精度で判定する

### 根拠(実測、lldb で `-O0 -g -ffp-contract=off` ビルドの内部を直接読んだ)

`als_iterate`(`src/iTPS/core/full_update.cpp:287`)の停止条件は

```cpp
if (std::abs(Old_delta - delta) / std::abs(C_phi) < peps_parameters.Full_Convergence_Epsilon)
```

で、**コストの変化**を見ており**反復点の変化**を見ていない。T3-iv の 5 ケースの実測:

| ケース | フェルミオン経路 | ボゾン経路 | 状態のずれ |
|---|---|---|---|
| 4100 exact h | conv=1 count=2 | conv=1 count=2 | 5.03e-16 |
| 4200 exact v | conv=1 count=2 | conv=1 count=2 | 3.33e-16 |
| 4300 generic h | conv=1 count=13 | conv=1 count=13 | 1.22e-15 |
| **4400 generic v** | conv=1 **count=37**, delta=1.3044117405847579e-12 | conv=1 **count=38**, delta=1.3044117405847622e-12 | **7.67e-15**(ローカル)/ **4.92e-08**(CI・`-ffp-contract=off`) |

落ちるケースだけ反復数が 1 回ずれ、**コストは相対 3.3e-15 で一致している**。
37 反復かかるのは最適点が平坦だからで、平坦な谷でコストを基準に止める限り、解は `√(ε_cost)` の
精度でしか決まらない(`√1e-14 = 1e-7`、観測 4.92e-8 と整合)。

**したがって「2 経路の解が一致する」は ALS が保証していない主張である。** 保証しているのはコストの一致。

### 要求

| # | 場所 | テスト | 現状 |
|---|---|---|---|
| 2-a | `test/fermion/full_update_bond.cpp:539` | T3-iv(フェルミオン vs ボゾン) | tol=1e-8、余裕 1.05e6(ローカル)/ 赤(CI) |
| 2-b | `test/fermion/boson_equivalence_full.py.in` | BosonEquivalenceFull(同じ 2 経路の E2E) | TOL=1e-6、実測余裕 45 |
| 2-c | `test/fermion/full_update_realctm.cpp:883` | `CHECK(r_new > r_old)` | 許容誤差なしの裸の不等号 |

**2-a**: 状態の比較の許容誤差を
`std::max(1.0e-8, 10.0 * std::sqrt(params.Full_Convergence_Epsilon))` とする
(この設定では `√1e-14 = 1e-7` なので 1e-6、観測 5.78e-8 に対して余裕 17 倍)。
この式は `fub_check_proportional` の呼び出し側で作り、**なぜその値かをコメントに書く**こと
(「ALS はコストで止まるので解は √ε_cost しか決まらない」)。

同時に**検出力を落とさないための追加要求**: T3-iv のケース集合のうち、
**exact 行(厳密表現可能なゲートで残差ゼロ、実測 5e-16 以下)については 1e-12 の厳しい許容誤差で判定する**こと。
平坦性は truncation を伴う generic 行だけの性質であり、exact 行を緩める理由は無い。
すなわち T3-iv は「exact 行は 1e-12、generic 行は `max(1e-8, 10√ε)`」の 2 段構えにする。

**2-b**: 許容誤差 `TOL` は**動かさない**。代わりに、比較する 2 つの run の入力 toml に
`[parameter.full_update]` の収束設定を明示的に書き、既定値(`Full_Convergence_Epsilon = 1e-6`、
`Full_Inverse_precision = 1e-12`)への依存を断つこと。既定のままだと `ε_cost = TOL` かつ
`√(Full_Inverse_precision) = TOL` で、比較の不定性が許容誤差と同じ大きさになる。
`Full_Convergence_Epsilon` を 1e-12 以下に絞ったうえで観測量の相対誤差を実測し、
**報告に「絞る前 / 絞った後の最悪相対誤差」を数値で載せる**こと。絞った後の余裕が 1000 倍に届かない場合は、
TOL を動かす前に必ず報告して指示を仰ぐこと(勝手に TOL を緩めない)。

**2-c**: 裸の `>` に、比較している量の水準に見合った下限を入れる。実測差は h: 1.6e-5、v: 2.5e-5 なので、
`CHECK(r_new > r_old + 1.0e-7)` のように、②で観測された経路乖離(5e-8)より大きく、
実測差より 2 桁以上小さい下限とすること。下限の値の根拠(実測差)をコメントに書く。

### 反空洞化(必須)

- 2-a の 2 段構えで、exact 行が 1e-12 で通ることを実測で示すこと。
- 変異検査: フェルミオン経路の符号を 1 箇所反転させた(たとえば `wrap_twosite_gate` の
  swap mask を落とす)コピーで T3-iv が**赤になる**ことを確認すること。緩めた許容誤差でも
  符号の誤りは O(1) のずれを生むので検出できるはずである。これが確認できなければ緩めすぎである。

---

## 検証(このタスクの完了条件)

1. `cmake --build --preset gcc-release` と `ctest --preset gcc-release` が **37/37 緑**。
2. **`-ffp-contract=off` ビルドでも緑**であること。手順:
   ```
   cmake -S . -B <tmp>/build-nofma -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-16 \
     -DCMAKE_CXX_FLAGS="-ffp-contract=off" -DTesting=ON -DGIT_SUBMODULE=OFF \
     -DTENES_PYTHON_EXECUTABLE=<repo>/venv/bin/python3
   cmake --build <tmp>/build-nofma --target test_fermion_layer -j 8
   <tmp>/build-nofma/test/test_fermion_layer
   ```
   これは CI(macOS g++-13)の②の失敗をビット一致で再現する構成である。
   **これが緑にならないうちは完了ではない。**
3. `-O0 -g -ffp-contract=off` ビルドでも `test_fermion_layer` が緑。
4. ①は手元のどのビルドでも再現しない(2 ulp の話なので)。したがって①については、
   **変更後の `tol` が CI で観測された `|diff|`(3.25261e-19)を上回ることを数値で示す**こと。
5. `src/` を変更していないこと(`git status --short src/` が空)。
6. 報告 `work/fermion/ci-tolerance/report-testauthor.md`。上の実測値をすべて載せること。
