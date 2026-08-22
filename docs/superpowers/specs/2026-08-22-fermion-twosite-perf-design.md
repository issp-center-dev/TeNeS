# フェルミオン二体測定の定数倍高速化 設計書

日付: 2026-08-22
ブランチ: `fermion`(タグ `wip/fermion_20260822` が改訂前の現状を指す)
前提: `docs/superpowers/specs/2026-08-15-fermion-m1-design.md`(graded 規約)、
`docs/superpowers/notes/2026-08-19-fermion-implementation-guide.md`(実装解説)。

改訂1 (2026-08-22): Codex レビューと追加のソース確認を反映。主な変更は
(a) OpenMP 領域の前に `prep_local_to_global()` / `make_l2g_map()` を呼ぶ必要がある
(§3.1 — mptensor の遅延初期化と競合する)、(b) **`global_index_fast` ではなく
`global_index_l2g_map` を使う**。要素あたりの除算が 16 回から 2 回になり、当初設計より
大きく効く(§3.1、§3.3)、(c) `1 << rank` のシフト UB を避けるガード順序(§3.1)、
(d) bit-identical の契約を「代数的同値、有限値では bit-identical」に精密化(§5)、
(e) Phase 1 後の支配項は `contract_reduced_pair_*_density_CTM` の BLAS になる見込みで、
§2.4 の gauge は既に OpenMP 並列なので優先度が下がる(§2.4、§4)、
(f) §3.4 のメモリ内訳の帰属を訂正(遅延 transpose ではなく値コピーと符号スイープが原因)。

改訂2 (2026-08-23、実装完了後): §4 を予測から**実測**に差し替え、§7 のタスク4 を
見送りとし、§9 に実測で判明した最優先項を追加した。実装は
`docs/superpowers/plans/2026-08-22-fermion-twosite-perf.md` に従って完了している
(コミット `010fd560`, `0d59c0bc`, `c987f451`, `71d25c3e`, `7297bc4c`, `bf6bce98`)。
予測と実測がずれた点を 3 つ残しておく:
(a) D=4 の observable は 30〜90 秒と見込んだが **138.8 秒**だった。支配項が
`contract_reduced_pair_*_density_CTM` の BLAS ではなく **mptensor 自身の
`reshape`(52%)と `contract`(27%)** に移っていたためである(§4、§9.0)。
(b) §3.4 のメモリ見積もり 4.3 GB は、mptensor が確保する索引スクラッチ
(12 B/要素 = D=4 で 3.0 GiB)を数え落としていた。**削減幅**(rank-16 の実体 2 本ぶん)は
予測どおりで、ピークヒープは 11.13 → 7.13 GiB(§4)。
(c) §2.4 で「4 秒程度」と見た gauge は、observable が 138.8 秒まで落ちた結果
**2.5% 未満**になり、タスク4 を見送った(§7)。

## 1. 目的とスコープ

フェルミオン模式の二体演算子測定が D の増加に対して実用不能な速度になっている。
本設計は **数値結果を変えずに**、その定数倍を落とす。

### 実測(本設計の出発点)

`work/fermion-perf/` に spinless free fermion (2×2 セル、χ = 2D²、simple update 200 step)
のベンチ入力を置いた。Release ビルド、g++-16、macOS、`_NO_MPI`、シングルプロセス。
`work/` は `.gitignore` 済みでコミットされないため、`work/fermion-perf/gen.py`
(ビルド済み `test/fermion/free_fermion.py` の `input_toml()` を再利用)で再生成する。
ベースライン出力と時間は `work/fermion-perf/baseline/` に保存してある。

| | D=2 (χ=8) | D=3 (χ=18) | D=4 (χ=32) | D=3→4 比 |
|---|---|---|---|---|
| simple update (200 step) | 4.7 s | 4.9 s | 5.3 s | 1.1× |
| environment (CTM) | 0.97 s | 2.48 s | 74.7 s | 30.1× |
| **observable (二体測定)** | **0.52 s** | **58.0 s** | **1857 s** | **32.0×** |
| 最大 RSS | — | — | 7.96 GB | |

- observable の D=3→4 比 32.0 は (4/3)¹² = 31.6 に一致する。**D¹² スケーリングである**。
- environment の 30.1 倍は射影子 SVD の行列サイズ χ·D² に対する (χD²)³ = 31.6 倍に一致する。
  すなわち **CTM は SVD 律速であり、フェルミオン固有のペナルティを負っていない**
  (ボゾンの density CTM と同じコスト構造)。D=4 で全体の 4% にすぎない。

### やること

`build_reduced_pair` 経路の要素走査を削る。数学は変えない(§5 の契約)。

### やらないこと

- **blob(rank-16 中間テンソル)の撤廃**。サイト毎の `build_reduced_op` を合成して
  ボゾンの `Contract_two_sites_*_op12_density_CTM` に流せば D¹² は計算からもメモリからも
  消えるが、二体演算子に対する符号規約の再導出が必要になる。**Phase 2 として分離する**
  (ユーザー決定、2026-08-22)。§9 に申し送る。
- CTM の高速化。上記のとおり伸びしろがない。次の一手は RSVD 解禁
  (現在 `src/iTPS/load_toml.cpp:646` で禁止)だが精度の議論を伴うため **Phase 3** とする。
- `deps/mptensor` の変更。公開 API の範囲内で済ませる(§3.3)。
- ボゾン経路の挙動変更。

## 2. ボトルネックの同定

### 2.1 プロファイル(macOS `sample`、D=4 の各フェーズ)

CTM フェーズ: 1451/1513 サンプル (96%) が
`Calc_projector_left_block_single` → `mptensor::svd` → LAPACK。

二体測定フェーズ: 14779/14779 サンプル (100%) が単一の経路にある。

```
measure_twosite
└ fermion::build_reduced_pair
  └ fermion::detail::fuse_doubled_cluster
    └ fermion::detail::apply_joint_swaps          ← 100%
      └ fermion::apply_swap
        └ mptensor::Tensor::global_index_fast
          └ div  (libsystem_c.dylib / DYLD-STUB$$div)
```

BLAS は一度も現れない。**この測定は浮動小数点演算律速ではなく、整数除算律速である。**

### 2.2 なぜそうなるか

`build_reduced_pair` (`src/fermion/reduced.hpp:196`) は

```
ket_ab  = tensordot(TnA, TnB)            rank 8   (l,t,b,s | t,r,b,s)
ket_op  = apply_pair_op(ket_ab, op12)    rank 8
doubled = tensordot(conj(ket_ab), ket_op, Axes(), Axes())   ← rank 16
```

で **rank-16 テンソル**を実体化する。D=4, d=2 で 2.68e8 要素 = 2.15 GB。
続く `fuse_doubled_cluster` (`src/fermion/reduced.hpp:111`) が
`apply_joint_swaps` を呼び、その中で `apply_swap` が **16 回**呼ばれる。

`apply_swap` (`src/fermion/fops.hpp:209`) は毎回全要素を走査し、要素ごとに
`global_index_fast` を呼ぶ。`global_index_fast`
(`deps/mptensor/include/mptensor/tensor_impl.hpp:539`) は軸ごとに `std::div` を呼ぶ。

> 16 パス × 2.68e8 要素 × 16 軸 ≈ **6.9e10 回の `div` ライブラリ呼び出し**、シングルスレッド。

その後 `ftensor::transpose` (`src/fermion/ftensor.hpp:97`) がもう 1 パス走り、
こちらは要素ごとに `transpose_sign` が O(rank²) = 256 回の `std::vector<bool>` アクセスを行う。

### 2.3 swap ペアの数え上げ(コードから機械的に再現、Codex 確認済み)

`kDoubledJointMask` が立てる方向ペアは `(0,3), (1,0), (2,3), (3,0)` の 4 本。
`apply_joint_swaps` はそれぞれについて `leg_ids` の一致位置を全て走査し、
**cross** (`ket_axes[ix]`, `bra_axes[iy]`) と **bra-bra** (`bra_axes[ix]`, `bra_axes[iy]`)
の 2 本ずつを適用する。

| 呼び出し元 | rank | cross | bra-bra | 合計パス |
|---|---|---|---|---|
| `doubled_pipeline`(サイト毎) | 10 | 4 | 4 | 8 |
| `fuse_doubled_cluster` horizontal | 16 | 8 | 8 | **16** |
| `fuse_doubled_cluster` vertical | 16 | 8 | 8 | **16** |

具体的なペア(horizontal、bra = 軸 0..7、ket = 軸 8..15):

```
cross  : (8,2) (8,6) (9,0) (12,0) (13,2) (13,6) (10,0) (14,0)
bra-bra: (0,2) (0,6) (1,0) (4,0)  (5,2)  (5,6)  (2,0)  (6,0)
```

**`apply_swap(a, x, y)` は `parity[ax1] && parity[ax2]` しか見ないので x, y について
対称であり、bra-bra には同じ swap が 2 回現れて相殺する組がある。** horizontal では
`{0,2}` と `{2,0}`、`{0,6}` と `{6,0}` がそのまま打ち消し合い、実効は
`{0,1} {0,4} {2,5} {5,6}` の 4 本。vertical も同様に 8 → 4 本。
`doubled_pipeline` は 4 → 2 本。**現状コードは無駄なパスを走らせている。**

### 2.4 副次的なホットスポット(優先度は低い)

`build_reduced_pair` 末尾の `apply_fused_leg_gauge`
(`src/fermion/reduced.hpp:68`)は rank-6 の blob (D¹² = 1.67e7 要素) に対して
`multiply_vector` を **2 回**呼ぶ。

ただし mptensor の `multiply_vector` (`tensor_impl.hpp:712-730`) は
**既に OpenMP 並列化されており**、`prep_local_to_global()` を先に呼ぶ正しい形になっている。
D=4 では 1.67e7 × 6 軸 × 2 回 × 16 blob ≈ 3.2e9 除算を 8 スレッドで割って **4 秒程度**と
見積もる(改訂前は逐次と誤認して 30 秒としていた)。§3.1 のカーネルで 2 スイープを
1 スイープに畳み、除算も減らせば 1 秒未満になるが、**支配項ではない**。

## 3. 設計

### 3.1 中核: 1 スイープの符号カーネル

新規ヘッダ `src/fermion/sign_sweep.hpp` に、局所要素を **1 回だけ**走査して
次を一括適用するカーネルを置く。

```
element *= (-1)^{ Σ_{(x,y) ∈ pairs} p_x(idx[x]) · p_y(idx[y]) }     … swap 群
         × (-1)^{ transpose_sign(parity, idx, axes) }               … graded transpose(任意)
         × Π_k  vec_k[ idx[ax_k] ]                                   … 脚ごとの対角因子(任意)
```

実装の要点:

1. **索引デコードを `global_index_l2g_map` に替える**(本改訂の最大の変更)。
   `Tensor::make_l2g_map()` は行・列それぞれについて軸分解表を先に作り
   (`tensor_impl.hpp:265-319`)、`global_index_l2g_map(lindex, gindex)`
   (同 330-350)は行・列インデックスを取り出す **2 回の除算**と表引きだけで済ませる。
   `global_index_fast` の **軸ごと 16 回**の `std::div` ライブラリ呼び出しが消える。
   - 表のサイズは要素数ではなく行数・列数に比例する。D=4 の rank-16
     (2.68e8 要素、l_row ≈ l_col ≈ 16384)で約 2 MB、構築コストは 2.6e5 除算で無視できる。
   - どちらも `mptensor::Tensor` の public メンバである(`tensor.hpp:144-145`)。
   - `make_l2g_map()` はキャッシュを持たず毎回作り直すので、**遅延 transpose で
     `axes_map` が変わった後に呼べば正しい表になる**。スイープの直前に呼ぶ。
   - 素性の確認: 本 submodule は `fe0540b`(`yomichi/fix/l2g-map-empty-range` のマージ)を
     指しており、`make_l2g_map` の空レンジ UB は修正済みである。
     `global_index_l2g_map` は mptensor 内部での使用例が少ないため、
     §6 層1 で `global_index_fast` と**全要素一致**を直接検証する。
2. **並列化の前処理**(Codex 指摘、必須)。`global_index()` は ScaLAPACK バックエンドで
   `prep_local_to_global()` を**遅延的に**呼び `global_row`/`global_col` を書き込む
   (`matrix_scalapack_impl.hpp:164-171, 238-`)。OpenMP 領域の中でこれが起きると競合する。
   mptensor 自身の並列カーネルは例外なく領域の外で先に呼んでいる。
   本カーネルも **`prep_local_to_global()` と `make_l2g_map()` を `#pragma omp parallel`
   の前に呼ぶ**。現行の `apply_swap` は逐次なのでこの前処理を持たない。そのまま並列化すると
   バグる。
3. **ビットマスク化**: 軸 `ax`、添字 `i` に対し `bit[ax][i] = parity[ax][i] ? (1u<<ax) : 0`
   を `std::vector<std::vector<std::uint32_t>>` で前計算する。`std::vector<bool>` の
   ビット取り出しを内側ループから追い出す。
4. **符号テーブル**: 要素ごとにマスク `m` (rank ビット) を組み立て、
   前計算した `std::vector<std::int8_t> table(std::size_t(1) << rank)` を引く。
   `transpose_sign` の O(rank²) = 256 演算が **1 回のテーブル参照**になる。
5. **テーブルを作らない分岐**: 判定は **必ず rank を先に見る**
   (`rank > kMaxTableRank || (std::size_t(1) << rank) > local_size()`)。
   `1 << rank` を先に評価すると rank ≥ 32 でシフト UB になる(Codex 指摘)。
   `kMaxTableRank = 24` とする。分岐した側はマスクから直接評価する。
   **両経路は同一の値を返さねばならない**(§6 層1 の検証対象)。
6. **OpenMP 並列化**: 要素は互いに独立なので `#pragma omp parallel for` を掛ける
   (`#ifndef _NO_OMP` で保護、索引バッファはスレッドごとに確保)。
   測定経路に外側の OpenMP 領域はない(確認済み)。
7. **対角因子は掛け合わせない**。複数の `vec_k` を事前に積んでから 1 回掛けると丸めが
   変わりうる。同一スイープ内で**個別に乗算**する(いまの用途では全て ±1 なので厳密だが、
   一般の因子に備えて規約として決めておく)。
8. **小テンソルは直列経路に落とす**(実装時に追加、2026-08-23)。項目 1・6 の仕掛けには
   要素数によらない固定費がある(`make_l2g_map` と OpenMP 領域の生成)。`tensordot` は
   `apply_transpose_sign_mask` を 1 回につき 2 回呼ぶので、simple update のように
   数百要素のテンソルを大量に扱う経路ではこの固定費が支配的になる。**実測で 8 スレッド、
   local 16 要素なら 500 倍遅、256 要素で 53 倍遅。実際に simple update が 20% 退行した。**

   対策として `detail::kSerialSweepThreshold`(現在 **8192**)を置き、`local_size` が
   これ未満なら **`make_l2g_map` を呼ばず、OpenMP 領域も作らず、`global_index_fast` で
   索引をデコードする直列経路**に落とす。`#pragma omp parallel if(...)` では不十分である
   (`make_l2g_map` は mptensor 側の固定費で、その中にも OpenMP 領域がある)。
   直列経路でも結果は変わらない — `global_index_fast` と `global_index_l2g_map` が
   同じ添字を返すことは §6 層1 が全要素で固定している。

   **8192 は機械依存の較正値である。** 8 スレッドの実機で分岐点の両側を測って決めた
   (4096 では直列 0.142ms vs 並列 0.212ms で直列が勝ち、8192 では直列 0.298ms vs
   並列 0.237ms で並列が勝つ)。4 スレッドでは真の分岐点が ~3000 なので、
   3000〜8192 の帯では 1 回あたり ~0.05ms を取り逃すが、多スレッド側の +0.24ms/回 の
   利得が上回るため現状維持としている。**本設計のベンチで実際に現れる要素数は
   256 / 2916 / 16384 で、この帯は実質空である。**

公開する入口は 4 つ:

| 関数 | 用途 |
|---|---|
| `apply_swap_form(a, form)` | swap 群のみ(外積前の小テンソル) |
| `transpose_with_swap_form(a, form, axes)` | swap 群 + graded transpose 符号 + 置換(rank-16 の本命) |
| `apply_leg_gauges(a, gauges)` | 脚ごとの対角因子をまとめて(§2.4 の blob) |
| `ftensor::transpose`(書き換え) | svd/qr/tensordot 経由の既存経路すべてが恩恵を受ける |

**未接続の入口について(2026-08-23 追記)**: §7 のタスク4 を見送った結果、
`apply_leg_gauges` と `LegGauge`、および `apply_sign_sweep` の `gauges` 引数は
**本番コードからの参照が 0 件**である。テスト(§6 の C1-8)は網羅しているので動作は
保証されているが、実際に使われるのは §7 タスク4 を再開したときになる。
`SignEval` も本番は `automatic` しか渡さず、`table` / `direct` を明示するのは
テスト(C1-6、C2-3)だけである — これは両経路が一致することを検査するために
意図的に公開している。

**`ftensor::transpose` の定義位置**: 宣言は `src/fermion/ftensor.hpp` にあるが、
定義は `src/fermion/sign_sweep.hpp` の末尾にある。`ftensor.hpp` が末尾で
`sign_sweep.hpp` を include し `sign_sweep.hpp` が先頭で `ftensor.hpp` を include する
循環を避けるためで、実装中に一度この循環を作って `sign_sweep.hpp` の単独 include が
壊れた(§8 のリスク表を参照)。`ftensor.hpp` だけを include する翻訳単位は
リンクエラーになるので、増やす場合は `sign_sweep.hpp` を include すること。

`SwapForm` は `{x,y}` を正規化(`x<y`)して **XOR 畳み込み**で保持する。
これにより §2.3 の相殺が自動的に効く。

### 3.2 bra-bra swap を外積の前に移す

`apply_swap` は要素ごとの ±1 倍であり、外積(`tensordot(·,·,Axes(),Axes())`)と可換である。
bra-bra swap は第一因子 `conj(ket_ab)` の軸だけに作用するので、
**外積の前に rank-8 テンソル(16384 要素)に適用できる。**

**この可換性が成り立つ根拠は「軸が空である」ことに依存する**(Codex 指摘)。
graded `tensordot` は `apply_transpose_sign_mask` で被演算子に符号マスクを掛ける
(`src/fermion/fops.hpp:288-293`)が、`Axes(), Axes()` では左右どちらの置換も恒等なので
マスクは全て +1 になる。**縮約軸が空でない tensordot に同じ移動を一般化してはならない。**

そのため `fuse_doubled_cluster` のシグネチャを変える:

```cpp
// 変更前
tensor fuse_doubled_cluster(const ftensor<tensor>& doubled,
                            const std::vector<int>& leg_ids);
// 変更後
tensor fuse_doubled_cluster(const ftensor<tensor>& bra_pair,
                            const ftensor<tensor>& ket_pair,
                            const std::vector<int>& leg_ids);
```

外積を関数の内側で行い、`bra_pair` に bra-bra form を先に適用する。
呼び出し元は `build_reduced_pair` の 1 箇所のみ
(`fuse_doubled_cluster(conj(ket_ab), ket_op, leg_ids)`)。
bra-bra form は `conj()` の**後**に適用する(現状も conj 後の因子から作られた doubled に
適用されているため)。`doubled_pipeline` は既に `(bra_Tn, ket_Tn)` を受け取っているので
同じ処理を内側で行う。

**結果**: rank-16 テンソルを触るのは cross form 8 本 + transpose 符号 = **1 スイープ**だけになる。

### 3.3 自前の索引デコーダを書かない理由

`Mat.global_index(lindex, g_row, g_col)` の後の軸分解は `axes_map` / `upper_rank` に
依存し、どちらも `mptensor::Tensor` の private メンバである。
非分散(`_NO_MPI`)なら局所添字が大域添字の単純な列優先分解なのでオドメータで
除算を完全に消せるが、ScaLAPACK バックエンドではブロックサイクリック分割のため成り立たない。
バックエンドで分岐する実装は MPI ビルドだけ検証が手薄になる典型的な事故源なので採らない。

代わりに §3.1 の項目 1 のとおり **mptensor の公開 API (`make_l2g_map` /
`global_index_l2g_map`) を使う**。これはバックエンドに依存せず、除算を軸ごと 16 回から
行・列の 2 回に減らす。`deps/mptensor` 自体は変更しない。

### 3.4 冗長コピーの削除

rank-16 テンソルは現状 3 回複製される。内訳は次のとおりで、
**`mptensor::Tensor::transpose(const Axes&)` メンバが遅延評価であることとは無関係**である
(Codex 指摘)。実体化の原因は値コピーと符号スイープ、および `reshape` の再分配である。

| 出所 | 根拠 |
|---|---|
| `ftensor<tensor> prepared = doubled;` | `src/fermion/reduced.hpp:114` |
| `fermion::transpose(prepared, axes)` の値返し | `src/fermion/fops.hpp:280-285`(`ftensor ret = a;` してから符号スイープ) |
| `mptensor::reshape(ordered.t, sh)` | `tensor_impl.hpp:1163`(再分配を伴い必ず実体化する) |

§3.2 の書き換えで外積結果を関数内のローカル変数として **その場で書き換える**形になり、
前二者が消える。

D=4 でのピークメモリ見積もり:

| | 現状 | 変更後 |
|---|---|---|
| `doubled` (rank 16) | 2.15 GB | 2.15 GB |
| `prepared` コピー | 2.15 GB | — |
| `fermion::transpose` の値返し | 2.15 GB | — |
| `reshape` の結果 | 2.15 GB | 2.15 GB |
| 合計(概算) | ~8.6 GB(実測 7.96 GB) | **~4.3 GB** |

## 4. 効果(実測、2026-08-23)

実装後に測り直した結果である。改訂前は「D=4 observable が 30〜90 s」と見積もっていたが、
実測は 138.8 s だった。**外れた理由も含めて下に記す。**

測定条件は §1 と同じ(`work/fermion-perf/`、Release、g++-16、macOS、`_NO_MPI`、
シングルプロセス、8 performance core)。他のジョブが走っていないことを確認して測定した。

| | D=2 (χ=8) | | D=3 (χ=18) | | D=4 (χ=32) | |
|---|---|---|---|---|---|---|
| | 変更前 | 変更後 | 変更前 | 変更後 | 変更前 | 変更後 |
| simple update | 4.718 | 4.67–4.74 | 4.868 | 4.995 | 5.276 | 5.059 |
| environment (CTM) | 0.973 | 0.963–0.985 | 2.478 | 2.456 | 74.705 | 73.915 |
| **observable** | **0.517** | **0.0767–0.0783** | **57.977** | **4.016** | **1857.22** | **138.85** |
| 全体 | 6.213 | ~5.7 | 65.329 | 11.471 | 1937.2 | 217.8 |
| **observable の倍率** | | **6.7×** | | **14.4×** | | **13.4×** |

- **D=2 / D=3 / D=4 の出力 9 ファイル(`twosite_obs.dat` / `onesite_obs.dat` / `density.dat`)は
  すべて変更前とバイト一致した。** §5 の契約は満たされている。
- environment は全 D で変化なし。§1 のとおり CTM は SVD 律速で、本変更の対象外である。
- simple update も変化なし。ただし途中で 20% の退行を作り込み、レビューで検出して直した
  (小テンソルに OpenMP の固定費が乗る問題。§3.1 の項目 8 に置いた
  `kSerialSweepThreshold` がその対策)。

### メモリ

| | 変更前 | 変更後 |
|---|---|---|
| D=4 ピークヒープ(`operator new` 差し替えで測定) | 11.13 GiB | **7.13 GiB** |
| rank-16 実体 1 本(2.0 GiB)に対する倍率 | 5.56× | 3.56× |

**削減幅はちょうど rank-16 の実体 2 本ぶんで、§3.4 の予測どおりである。**
一方 §3.4 の「合計 4.3 GB」という絶対値は誤っていた。mptensor の `reshape` / `contract` が
確保する索引スクラッチ(`vector<size_t> local_position` 8 B/要素 +
`vector<int> dest_mpirank` 4 B/要素 = D=4 で 3.0 GiB)を数え落としていたためである。

なお `/usr/bin/time -l` の最大 RSS は変更後で 10.18 GB だが、§1 に記した変更前の 7.96 GB は
実行中の `ps` サンプリングで取った値であり、**両者は測り方が違うので比較してはならない。**
比較可能なのは上のピークヒープ測定である。

### なぜ observable が予測の 30〜90 s ではなく 138.8 s だったか

再プロファイル(D=4、`sample` 25 秒、observable フェーズ)の結果、**支配項が移動していた**。

| | サンプル | 割合 |
|---|---|---|
| `mptensor::reshape` | 7943 | **52%** |
| `mptensor::contract` | 4129 | **27%** |
| 自前の `apply_sign_sweep` | 1228 | 8% |

改訂前の §4 は「Phase 1 後の支配項は
`contract_reduced_pair_*_density_CTM` の BLAS 縮約(D¹²·χ²)」と予想していたが、
実際にはそこに到達する前の **mptensor 自身の `reshape` と `contract`** が
79% を占めていた。どちらも内部で `global_index_fast`(軸ごとに libc の `std::div`)を
要素ごとに呼んでいる。§3.3 で「`deps/mptensor` は変更しない」と決めた領域である。

要素ごとの全走査の回数はレビューで実測されており、**24 回から 6 回**に減っている
(内訳は tensordot の書き込み 1 + 符号スイープ 1 + `reshape` 2 + `contract` 2)。
自前のスイープに限れば **17 回 → 1 回**で、これは設計どおりである。

## 5. 正しさの契約

本変更は符号(厳密な ±1)の **適用順序と走査方法**だけを変える。適用する符号の集合は
不変である。

> **契約: 変更前後で代数的に同値であること。
> 検証としては、有限値の入力に対して全出力が bit-identical であることを要求する。**

bit-identical を要求できる根拠: IEEE-754 では符号反転は厳密であり、乗算は符号対称
(`(-x)*y` と `-(x*y)` は指数・仮数が同一)なので、±1 因子を積の前後どちらに置いても
丸めが変わらない。複素数の積や FMA を経ても同様である。OpenMP 化も要素独立なので
加算順序を変えない。§3.1 項目 7 のとおり対角因子を事前に積まないのも同じ理由による。

例外は NaN のペイロード符号ビットのみで、これは NaN が出た時点で計算が破綻しているため
実務上問題にならない(Codex 指摘に対する判断)。

`test/data/output_*/` のゴールデンファイルとも厳密一致するはずである。
差が出た場合は「許容誤差内だからよい」ではなく**実装バグの兆候**として扱う。

## 6. 検証戦略

### 層1: 参照実装との突き合わせ(doctest、`test/test_fermion_layer.cpp`)

新カーネルの検証は「最適化前の素朴な実装をテストファイル内に参照として残し、
乱数テンソルで一致を見る」形にする。これは最適化リファクタに対する唯一の直接的な網である。

- `global_index_l2g_map` と `global_index_fast` が**全要素で同じ添字を返す**こと
  (§3.1 項目 1 の前提そのもの。rank と形状を変えて複数ケース)
- `apply_swap` の逐次適用 vs `apply_swap_form` の 1 スイープ
- 旧 `transpose_sign` (O(rank²) 版) vs テーブル版、および直接評価版
- テーブル経路と直接評価経路の一致(`local_size` を跨ぐ形状を両方用意する)
- rank 2〜16、偶奇の内訳を非対称にした parity ベクトル、実/複素の両方
- **一致は `==` で見る(許容誤差ではない)**

### 層2: 合成経路の同値性(doctest)

`build_reduced_op` / `build_reduced_pair` について、変更前の実装を参照として
テストファイルに写し取り、乱数 `Tn` で **bit-identical** を確認する。
horizontal / vertical の両方、d=2 と d=4 の両方。

### 層2b: 入力検証ガードの検査(2026-08-23 追加、最終レビューの指摘による)

**この層は当初の設計に無く、最終レビューで抜けが見つかって追加した。**

§3.1 の項目 5・8 や実装中の修正ラウンドで、「軸番号や台帳長が範囲外なら黙って壊れるのではなく
`std::runtime_error` を投げる」という規約を決めた。ところがこの規約は 3 つのタスクの
修正ラウンドにまたがって決まったため、**どの契約書にも検査項目として書かれず、
ガードが 1 つもテストされないまま残っていた。**

実測(テスト作成者が `throw` を 1 文ずつ削った 17 ビルドで確認):

- ガードは `sign_sweep.hpp` に 7 文、`reduced.hpp` に 9 文の**計 16 文**ある
- **16 文すべてを同時に無効化しても、層1・層2 の数値テストは 41 ケース
  4767 アサーションが緑のまま**だった。数値テストは正しい引数しか渡さないので、
  原理的にガードを見ない
- 16 文のうち **7 文は削除するとプロセスが SIGABRT / SIGTRAP で落ちる**。
  すなわちガードは実際にメモリ安全性を守っていた
- 既存のテストは 16 文のうち 14 文を見逃していた

**型だけを見る `CHECK_THROWS_AS` では不十分である。** `build_reduced_op` の rank 検査は、
削除しても下流の `validate_axes` が同じ型を投げるため緑のままだった。
**全アサーションを `CHECK_THROWS_WITH_AS` にしてメッセージまで厳密一致させること。**

教訓: **規約を修正ラウンドで決めたら、その場で契約書に検査項目を書き足すこと。**
そうしないと、規約は実装されるがテストされない。

### 層3: 既存テストの全件(`ctest`)

`ftensor::transpose` は `fermion::tensordot` / `svd` / `qr` が使う共有コードであり、
simple update も通る。したがって **フェルミオン関連だけでなく全 ctest を回す**
(共有コードを触って絞った範囲だけで破損を見逃した事例があるため)。
特に `FreeFermion`, `FreeFermionMF`, `FreeFermionSaveLoad`, `FreeFermionSimple`,
`test_fermion_layer`, `r2_convention`, `mf_measure`, `saveload`。

`_NO_OMP` ビルドでもビルドと ctest が通ることを確認する。

### 層4: ベンチマークと再プロファイル

`work/fermion-perf/` の D=2/3/4 を変更後に実行し、`baseline/` と比較する。

- `out_D4/twosite_obs.dat`, `onesite_obs.dat`, `density.dat` が `baseline/` と**バイト一致**
- `time.dat` と最大 RSS を記録し、§4 を実測値で置き換える
- D=4 の observable フェーズを `sample` で再プロファイルし、新しい支配項を特定して
  §9(Phase 2 申し送り)を更新する

### 変異テスト(レビュアー向け)

- `SwapForm` の XOR 畳み込みを「重複を無視して追加」に変えると層1・層2 が赤くなること
- テーブル構築の分岐条件を反転させても層1 が緑のままなら、その分岐はテストされていない
- `prep_local_to_global()` / `make_l2g_map()` の呼び出しを OpenMP 領域の中に移しても
  層1・層3 が緑のままなら、MPI ビルドでの競合を検出できていない
- OpenMP の `#pragma` を消しても結果が変わらないこと(変わるなら競合がある)

## 7. タスク分割(実装計画の骨子)

1. `src/fermion/sign_sweep.hpp` の新設 — カーネルと `SwapForm`、4 つの入口。
   単体で層1 のテストが通る状態にする。
2. `ftensor::transpose` をカーネル利用に書き換え。層3 の全 ctest を回す。
3. `doubled_pipeline` / `fuse_doubled_cluster` の書き換え(§3.2、§3.4)。層2 を通す。
4. ~~`apply_fused_leg_gauge` を `apply_leg_gauges` の 1 スイープに統合(§2.4、優先度低)。~~
   **見送った(2026-08-23)。** タスク3 完了後に D=4 を再プロファイルしたところ、
   `multiply_vector` は 400 サンプル以上の項に一度も現れず、**2.5% 未満**だった。
   §2.4 の「D=4 で 4 秒程度」という見積もりは observable が 1857 秒だった頃の相対値であり、
   138.8 秒まで落ちた現在は誤差の範囲である。§4 のとおり、いま手を入れるべきは
   `reshape` / `contract` であって gauge ではない。
5. ベンチマーク再測定・再プロファイルと設計書 §4/§9 の実測値差し替え。

各タスクの完了条件に **全 ctest 緑** と **bit-identical** を含める。

## 8. リスクと判断

| リスク | 対処 |
|---|---|
| OpenMP 領域内で `prep_local_to_global` が競合(MPI ビルドのみ顕在化) | 領域の前で呼ぶ。変異テストで検出可能性を確認 |
| `global_index_l2g_map` は mptensor 内での使用例が少ない | 層1 で `global_index_fast` と全要素一致を直接検証 |
| `make_l2g_map` が遅延 transpose 後に陳腐化 | キャッシュを持たない実装であることを確認済み。スイープ直前に呼ぶ |
| `SwapForm` の XOR 畳み込みで符号を落とす | 層1・層2 の参照実装比較(bit-identical) |
| OpenMP 導入で競合 | 要素独立。層1 を並列有無の両方で回す |
| テーブル経路の分岐が片方しかテストされない | 層1 で `local_size` を跨ぐ形状を両方用意する |
| `1 << rank` のシフト UB | rank を先に判定する順序を守る。境界を層1 で突く |
| `_NO_OMP` ビルドの破損 | `#ifndef _NO_OMP` で保護し、両構成でビルド・ctest |
| MPI ビルドでの破綻 | バックエンド依存の最適化を入れない(§3.3) |

## 9. Phase 2 への申し送り(本設計の範囲外)

### 9.0 実測で判明した最優先項: `mptensor::reshape` + `contract`(2026-08-23 追記)

§4 の再プロファイルのとおり、Phase 1 後の D=4 observable 138.8 秒のうち **79%** が
`mptensor::reshape`(52%)と `mptensor::contract`(27%)である。改訂前の §9 は
blob 撤廃を最優先と見込んでいたが、**実測はそれより手前にもっと大きな項があることを示した。**

どちらも `fuse_doubled_cluster` の末尾で rank-16 の実体(D=4 で 2.68e8 要素)に対して
呼ばれている:

```cpp
tensor fused = mptensor::reshape(doubled.t, sh);                       // rank 16 -> rank 10
return mptensor::contract(fused, mptensor::Axes(6, 7), mptensor::Axes(8, 9));
```

**有望な一手: 縮約を reshape の前に回す。** `contract` が潰すのは 4 本の物理脚であり、
仮想脚の融合(reshape)とは独立である。先に物理脚を縮約すれば、reshape が扱う要素数は
2.68e8 から **1.67e7 へ 16 分の 1** になる。時間だけでなく、§4 に記した索引スクラッチ
(12 B/要素)も同じ比で減るのでピークメモリにも効く。

**ただし bit-identical が崩れる可能性がある。** `contract` は物理脚について 4 項の和を取るが、
その加算順序は要素の記憶配置に依存し、reshape はその配置を変える。§5 の契約を
維持できるかどうかを最初に確かめること。崩れるなら、この案は採らないか、
契約を「代数的同値 + 許容誤差」に緩める判断が別途要る。

### 9.0b 本設計の範囲外として残した既存コードの課題(2026-08-23)

いずれも `wip/fermion_20260822` 時点から存在するもので、本変更が持ち込んだのではない。
性能改善の PR で無関係な既存コードに手を入れるのはスコープの逸脱なので、記録に留める。

- **`src/fermion/reduced.hpp` の `build_reduced_pair` にある `switch` の `default:` は
  到達不能**である。直前で呼ぶ `build_pair_state` が同じ判定で先に投げるため、
  どんなテストでもこの行を殺せない。`default:` の削除か `switch` の順序入れ替えで
  直せるが、`enum class` に対する `switch` から `default:` を外すと別の議論を呼ぶ。
- **`src/fermion/fops.hpp` にも同種の検証ガードが複数あり、同じく無テストである
  可能性が高い**(層2b と同じ穴)。本設計では `sign_sweep.hpp` と `reduced.hpp` の
  16 文だけを塞いだ。

その先は `deps/mptensor` 側の課題である。`reshape` / `contract` / `multiply_vector` は
いずれも要素ごとに `global_index_fast`(軸ごとに libc の `std::div` を**関数呼び出し**)を
使っている。本設計 §3.1 で TeNeS 側のスイープに対して行ったのと同じ置き換え
(`make_l2g_map` + `global_index_l2g_map`)を mptensor 内部に適用できれば、
TeNeS 全体が速くなる。fork(`yomichi/mptensor`)があるので上流課題として扱える。

### 9.1 blob 撤廃

blob 撤廃の見通しをここに残す。

- cross form は GF(2) 上で階数 2 に潰れる。horizontal では
  `(k₈⊕k₁₃)(b₂⊕b₆) + (k₉⊕k₁₀⊕k₁₂⊕k₁₄)·b₀`。したがって rank-16 への符号ループは
  原理的に完全に消せる(小テンソル 4 組の外積の和になる)。
- さらに進めて、サイト毎の `build_reduced_op` を合成し、符号ドレスした `op12` を
  ボゾンの `Contract_two_sites_horizontal_op12_density_CTM`
  (`src/iTPS/core/contract_density_ctm/ctm.cpp:199`)に流せば、**D¹²·χ² の環境縮約も
  D¹² のメモリも消える**。既存の blob 実装が参照解になるので検証は容易である。
- ノルム経路は既に合成方式(`build_reduced_identity_pair`)を使っており、
  `build_reduced_pair` の `apply_fused_leg_gauge` は「blob 方式を合成方式の規約に
  合わせる」ためのものである。つまり恒等演算子については両方式の一致が既に取れている。
  残る課題は演算子挿入時の符号ドレスの導出だけである。
- 演算子が偶パリティ同士の積に分解できる項では追加符号は出ず、奇×奇の項でのみ
  ボンド脚を跨ぐ対角符号が出ると見込まれる。`op12` を定パリティ項に分解して扱う。
