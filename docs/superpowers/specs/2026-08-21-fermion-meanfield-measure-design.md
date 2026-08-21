# フェルミオン模式の平均場環境(MeanField_Env)対応 設計書

日付: 2026-08-21
ブランチ: `fermion`
前提: `docs/superpowers/specs/2026-08-15-fermion-m1-design.md`(C++ 側の graded 規約)、
`docs/superpowers/notes/2026-08-19-fermion-implementation-guide.md`(実装解説)、
`docs/superpowers/specs/2026-08-20-fermion-tools-design.md`(ツール側)。

改訂1 (2026-08-21): Codex レビューを反映。(a) `build_pair_state` / `apply_pair_op` を
`reduced.hpp` に置く(`reduced_measure.hpp` → `reduced.hpp` の include 方向を保つ、§4.1)、
(b) 層2 に「wrap なし」の陰性対照を追加(§5)、(c) 相関長は入力ガードではなく実行時無効化で
あることを明記(§1)、(d) `save_density` の既存の脆さはスコープ外として記録(§8)。
数理(§3)には反論なし。

## 1. 目的とスコープ

フェルミオン模式の2サイト観測量は、現行では CTM 環境に載せる「reduced pair blob」
(`src/fermion/reduced.hpp` の `build_reduced_pair`)を経由する。この blob は bra⊗ket の
外積(rank 16、要素数 D¹²d⁴)なので、D>2 では測定が計算全体を支配し、d=4(Hubbard)の
D=4 はメモリ的に実行できない。本設計は **フェルミオン模式で `meanfield_env = true` を
使えるようにし**、2サイト観測量を単層縮約(要素数 D⁶d²)で評価できるようにする。

### 実測(spinless d=2、2×2 セル、測定ボンド1本、Release ビルド、macOS Apple Silicon)

| D | χ | simple update 20 step | CTM 環境 | twosite 1ボンド(norm+値) | ピーク RSS |
|---|---|---|---|---|---|
| 3 | 9 | 0.5 s | 0.7 s | 7.4 s | — |
| 4 | 16 | 0.6 s | 9.5 s | 237 s | 7.6 GB |

`sample` によるプロファイルではサンプル期間の 100% が
`build_reduced_pair → fuse_doubled_cluster → apply_joint_swaps → apply_swap` で、
rank-16 テンソルに `global_index_fast` の要素ループを 16 パスかけている部分である。
Hubbard(d=4)では blob の要素数が D=3 で 1.4×10⁸、D=4 で 4.3×10⁹(>30 GB)になる。

### やること

- `validate_fermion_constraints` の `MeanField_Env=true` 拒否を外す
- `measure_twosite` の MF 分岐にフェルミオン経路を追加する(norm と値の両方)
- そのための単層縮約 API を `src/fermion/reduced.hpp` / `reduced_measure.hpp` に追加する
- 1サイト観測量は既存のボソン MF 式をそのまま使う(§3.4 で理由を述べる)。
  ただしテストで担保する
- Fock オラクル(`test/fermion/fock_oracle.py`)を「開いた脚を任意ラベルに固定できる」
  ように拡張し、MF の値を Fock 空間の厳密計算と突き合わせる
- ドキュメント(ja/en の非対応一覧)と `NEWS.md` の更新

### やらないこと

- CTM blob 経路の高速化(`apply_joint_swaps` の 1 パス化、不純物分解による blob 廃止)。
  別タスク。本設計は CTM 経路のコードに触れない
- MF での相関関数・相関長・マルチサイト・有限温度・実時間。M1 のガードを維持する。
  相関長だけは入力ガードではなく `measure()` 冒頭の実行時無効化(`measure.cpp:89-97`)で
  止まっており、これは MF 分岐より前に評価されるので MF でも変わらない
- `tenes_simple` / `tenes_std` の変更。`meanfield_env` は `[parameter.ctm]` の通し項目で、
  ツールは fermion と MF の組を現在も拒否していない(変更不要)
- MF での full update(`PEPS_Parameters` が MF + full update を既に拒否、fermion も
  full update を拒否)

## 2. 全体方針

MF 測定は「2サイト窓の外側を λ² の対角重みで閉じる」近似であり、窓の中に閉路がない。
graded テンソル代数では **縮約の順序によらず値が一意** なので、blob(bra⊗ket を先に
組んで環境で閉じる)の代わりに、ket 層だけを組んで bra と直接内積を取ればよい:

```
ψ_AB  = tensordot_f(A_λ, B_λ)                 # rank 8、graded tensordot
value = trace_f(conj_f(ψ_AB), O ∘ ψ_AB)      # 全脚の graded trace
norm  = trace_f(conj_f(ψ_AB), ψ_AB)
```

この式は **既に Fock オラクルで検証済みの「direct path」そのもの** である
(`test/fermion/r2_convention.cpp` の `r2_expect_two` / `r2_norm`。R2 テストで Python の
Fock オラクルと、R5 テストで d=4 の電子ホッピングと一致を確認している)。ただし既存の
検証は開いた脚が 1 次元(真空)の場合に限られる。MF では開いた脚に奇ラベルが入るので、
そこを新たに検証する(§5)。

### 却下した代案

- **blob を組んで「恒等環境」で閉じる**: D¹²d⁴ を作る時点で目的に反する。却下。
- **ボソン MF 式(`Contract_two_sites_*_op12_iTPS_MF`)に符号補正を後付けする**:
  符号は `ψ_AB` の構成と演算子の作用(脚の交差)で生じ、bra と ket の接続では生じない
  (§3.2)。補正表方式は M1 で根絶したバグクラスの再導入である。却下。
- **CTM 経路の高速化だけで済ませる**: 16 パス→1 パス化で数倍は縮むが D¹²d⁴ のメモリは
  残り、d=4 の D=4 は救えない。MF とは目的が違う(近似精度 vs 速度)ので別タスクとする。

## 3. 数理と規約

### 3.1 λ の扱い

既存のボソン MF 経路(`src/iTPS/twosite_obs.cpp:117-148`)が窓の外側の脚に
`lambda_tensor` を掛けたコピー(`boundaries`)を作る。フェルミオン経路は **このコピーを
そのまま消費する**。λ は対角の実重みでパリティ偶なので、graded 構造に影響しない。
simple update の λ 規約(`core/simple_update.cpp:234-236`、`λ = sqrt(s/Σs)` を両端に吸収)は
ボソンとフェルミオンで共通のカーネルなので、dressing も共通でよい。

### 3.2 符号の所在

`trace_f(conj_f(ψ), φ)` で全脚を縮約すると、`conj_f` の符号 `(-1)^{m(m-1)/2}`(m は奇脚数)と
`tensordot_right_perm` の全反転符号が打ち消し合い、**素の内積 Σ conj(ψ)·φ** になる
(`test_fermion_layer.cpp` の「trace matches manual swap for reversed contracted ordering」
「norm of even vector tensor is positive via conj」が根拠)。したがって

- **norm**: `ψ_AB` 構成時の符号マスクは 2 乗で消えるので、フェルミオン MF の norm は
  ボソン MF の norm と**厳密に一致**する(正定値)。テストで固定する。
- **値**: 符号は `ψ_AB` の graded tensordot(A の脚 2 を末尾へ動かす交差)と、演算子を
  `s_A, s_B` に作用させる graded tensordot(`s_A` が B の仮想脚を跨ぐ)からのみ生じる。
  これが「ホッピングの JW 弦が隣のサイトの仮想脚を跨ぐ」符号に対応する。

### 3.3 演算子の wrap 規約(重要)

blob 経路は `wrap_reduced_pair_op`(入力脚と出力脚の両方に swap)を使うが、direct path は
**`wrap_twosite_gate`(入力脚のみ swap)** を使う。R5 テストの参照値
`hop_ref = r2_expect_two(psi, 3, 7, r4_wrap(hop, /*swap_in=*/true, /*swap_out=*/false))`
がこれを固定している。d=2 の粒子数保存演算子では両者は区別がつかない
((奇,奇)→(偶,偶) チャネルが無い)ので、**d=4 のホッピングでテストしないと取り違えが
検出できない**(§5 層2)。

この in-swap 規約は、単層での graded ゲート作用が d=4 の電子系でも厳密な Fock 空間の
JW 時間発展と一致することを確認した「plaquette kernel vs exact trotter」テスト
(`test_fermion_layer.cpp:2066`、`kernel_vs_fock` / `fprim_vs_fock`)でも裏付けられている。

source が窓の 2 番目のサイトである場合(`dx<0` / `dy>0`)は、blob 経路と同じく
wrap 後に **graded** transpose `(1,0,3,2)` を掛ける。これは単なる添字の入れ替えではない:
順序付き Fock 基底 `|n_B n_A⟩ = (-1)^{n_A n_B} |n_A n_B⟩` の並べ替え符号を
入力脚・出力脚の両方に付けたものに等しく、`wrap_twosite_gate` の in-swap と合わせると
「(A,B) 順序で書き直した行列に in-swap を掛けたもの」と厳密に一致する。
d=2 の粒子数保存演算子では素の transpose と区別がつかず、d=4 のホッピング
((奇,奇)→(偶,偶) チャネル)で初めて差が出るので、**d=4 で両端から測って一致する**
テストが必要(§5 層2)。

### 3.4 1サイト観測量

1サイト演算子はパリティ偶に限られ(ガード済み)、物理脚は末尾で交差がないので、
`tensordot_f(A, op, Axes(4), Axes(0))` に符号マスクは付かない。よって
`Contract_one_site_iTPS_MF`(ボソン式)がそのまま正しい。CTM のフェルミオン経路でも
1サイト演算子は wrap なしで載せている(`onesite_obs.cpp:87-102`)のと整合する。
コードは変えず、テスト(§5 層4)で Fock オラクルの密度と一致することを確認する。

## 4. 実装設計

### 4.1 単層縮約 API(新規関数)

`build_pair_state` / `apply_pair_op` は `build_reduced_pair` からも使うので `src/fermion/reduced.hpp` に、
`contract_pair_MF` は `src/fermion/reduced_measure.hpp` に置く。

```cpp
namespace tenes::fermion {

// ket 層の2サイト状態。blob 経路の build_reduced_pair と同じ脚順:
//   horizontal: tensordot(TnA, TnB, Axes(2), Axes(0))
//               → (l_A, t_A, b_A, s_A, t_B, r_B, b_B, s_B)
//   vertical:   tensordot(TnA, TnB, Axes(3), Axes(1))
//               → (l_A, t_A, r_A, s_A, l_B, r_B, b_B, s_B)
template <class tensor>
ftensor<tensor> build_pair_state(const ftensor<tensor>& TnA,
                                 const ftensor<tensor>& TnB,
                                 reduced_pair_direction direction);

// 2サイト演算子 op12 (in_A, in_B, out_A, out_B) を s_A, s_B に作用させる。
// 既存 detail::apply_reduced_two_site_op と同一(tensordot + transpose(0,1,2,6,3,4,5,7))。
template <class tensor>
ftensor<tensor> apply_pair_op(const ftensor<tensor>& pair,
                              const ftensor<tensor>& op12);

// MF norm:  trace_f(conj_f(pair), pair)
template <class tensor>
typename tensor::value_type contract_pair_MF(const ftensor<tensor>& pair);

// MF 期待値(非規格化): trace_f(conj_f(pair), apply_pair_op(pair, op12))
template <class tensor>
typename tensor::value_type contract_pair_MF(const ftensor<tensor>& pair,
                                             const ftensor<tensor>& op12);

}  // namespace tenes::fermion
```

- `build_reduced_pair` の冒頭(`ket_ab` の構成)は `build_pair_state` に置き換える。
  脚順の定義を 1 箇所にするための最小限のリファクタで、blob 経路の値は変わらない
  (既存 R3/R5 テストが回帰を捕まえる)。
- `detail::apply_reduced_two_site_op` は `apply_pair_op` に改名して公開する(呼び出し元は
  `build_reduced_pair` のみ)。
- 実装は `fops.hpp` の `tensordot` / `trace` / `conj` / `transpose` だけで書く。
  **`apply_swap` を手で呼ばない**(符号を手書きしない)。

### 4.2 `src/iTPS/twosite_obs.cpp` — 分岐構造

現在の分岐は `finfo.enabled && !is_TPO && !is_mf`(blob)と、それ以外(ボソン)の 2 択で、
fermion + MF は**ボソン式に落ちて黙って間違う**構造になっている(ガードで到達不能に
しているだけ)。これを次のように組み替える:

```
norm:
  if (finfo.enabled && !is_TPO && nrow*ncol == 2)
     if (is_mf)  norms[key] = contract_pair_MF(build_pair_state(wrap_Tn(*Tn_[0][0]), wrap_Tn(*Tn_[0][1]), dir))
     else        norms[key] = (既存 blob + CTM)
  else if (is_mf) (既存ボソン MF)
  else            (既存 CTM)

値(nrow*ncol == 2、op.ops_indices.empty()):
  if (finfo.enabled && !is_TPO)
     o = is_mf ? wrap_twosite_gate(op.op, p_src, p_tgt)
               : wrap_reduced_pair_op(op.op, p_src, p_tgt)
     if (source が窓の 2 番目) o = transpose(o, Axes(1,0,3,2))
     if (is_mf)  value = contract_pair_MF(build_pair_state(...), o)
     else        (既存 blob + CTM)
  else (既存ボソン)
```

- `wrap_Tn(*(Tn_[0][0]), finfo, left)` は λ を掛けたコピー(`boundaries`)を指す。
  既存の blob 分岐が同じポインタを使っているので、扱いは揃う。
- `pair` は norm と値で同じものなので、ボンドごとに 1 回だけ組んでよい(最適化は任意。
  D=4 でも rank-8 は 6.5×10⁴ 要素なので組み直しても問題ない)。
- `is_TPO`(有限温度)は fermion ガードで到達不能だが、条件には残す。

### 4.3 `src/iTPS/load_toml.cpp` — ガード解除

`validate_fermion_constraints` から `if (peps_parameters.MeanField_Env) throw_fermion_guard("MeanField_Env=true");`
を削除する。他のガード(full update、RSVD、Gauge_Fix、相関、多サイト、skew、1 幅セル、
`ops` 形式、距離 2 以上、パリティ奇)はそのまま。

### 4.4 `src/iTPS/measure.cpp` / `onesite_obs.cpp`

変更なし。`measure()` は `MeanField_Env` のとき環境更新を飛ばす既存の分岐で足りる。
1サイトは §3.4 のとおりボソン MF 式を使う。

### 4.5 ファイル一覧

| ファイル | 変更 |
|---|---|
| `src/fermion/reduced.hpp` | `build_pair_state` / `apply_pair_op` 追加(`detail::apply_reduced_two_site_op` の改名・公開)。`build_reduced_pair` がそれを使う(値は不変) |
| `src/fermion/reduced_measure.hpp` | `contract_pair_MF`(norm 版・演算子版)追加 |
| `src/iTPS/twosite_obs.cpp` | MF 分岐のフェルミオン経路 |
| `src/iTPS/load_toml.cpp` | `MeanField_Env` ガード削除 |
| `test/fermion/fock_oracle.py` | 開いた脚のラベル固定、MF 重み付き和、自己検査 |
| `test/fermion/r2_convention.cpp` または新規 `test/fermion/mf_measure.cpp` | 層1〜3 |
| `test/test_fermion_layer.cpp` | 層4(iTPS レベル) |
| `test/input.cpp` | 「fermion rejects mean-field environment」を「accepts」に置換 |
| `test/fermion/free_fermion.py.in` / `test/CMakeLists.txt` | 層5(E2E、MF 版) |
| `docs/sphinx/{ja,en}/file_specification/parameter_section.rst` | 非対応一覧から平均場環境を外し、精度の注意を追記 |
| `NEWS.md` | 項目追加 |

## 5. 検証戦略

符号規約の取り違えは「黙って間違う」ので、テストは Fock 空間の厳密計算を真とする。
各層は独立した根拠を持ち、合わせて三角測量になる。

### 層0: Fock オラクルの拡張(`test/fermion/fock_oracle.py`)

- `Oracle` に開いた脚ごとの固定ラベル `dangling_labels[(site, leg)] = x` を受け取らせる。
  奇ラベルの開いた脚には補助モードを割り当て、状態構成の最初に生成しておく
  (どの順で生成しても x ごとの全体符号にしかならず、bra と ket で打ち消す)。
  `apply_physical_projector` はその脚の添字を x に固定する。
- `mf_expectation(patch, tensors, parities, lambdas, op)` =
  `Σ_x Π_legs λ²_leg(x_leg) ⟨ψ(x)|O|ψ(x)⟩ / Σ_x Π λ² ⟨ψ(x)|ψ(x)⟩`。
- **自己検査**: (a) 開いた脚を全て 1 次元偶にすると既存 `observables()` と一致、
  (b) 全パリティ偶(ボソン扱い)にすると `plain_boson_norm` 系の値と一致、
  (c) 開いた脚の奇ラベルの生成順を入れ替えても期待値が変わらない。
- 参照値は `gen_reference.py` / `print_case` と同じ流儀で C++ テストに定数として埋める。
  **ホッピングの期待値が非零になるデータセットを選ぶ**(既存 R3 データセットでは
  真空射影のためホッピングが 0 で、符号を検出できない)。

### 層1: 単層縮約 API vs オラクル(d=2)

`build_pair_state` + `contract_pair_MF` を、λ を掛けた決定論的テンソル(水平・垂直、
seed 複数)で評価し、層0 の参照値(norm、密度、ホッピング、ペアリング、混合項)と
1e-12 で一致すること。source を 2 番目に置いた(transpose 付き)場合も含む。

### 層2: wrap 規約の固定(d=4)

R5 の `make_r4_tensor` / `r4_hop_plain` を使い、`contract_pair_MF(pair, wrap_twosite_gate(hop))`
が R5 の `hop_ref`(Fock 検証済み)と一致し、`wrap_reduced_pair_op(hop)`(両 swap)を渡しても、
wrap なしの `ftensor{hop, {p,p,p,p}}` を渡しても **一致しない**こと(取り違え検出の
陰性対照 2 件)。さらに同じ d=4 ホッピング(サイト交換対称)を
source を 2 番目に置いて(wrap 後に graded transpose)測った値が A 端からの値と一致し、
素の `mptensor::transpose` に置き換えると一致しないこと(§3.3)。

### 層3: 線形性・順序不変(開いた脚が奇の領域を C++ 側で閉じる)

開いた脚を `fermion::slice` で各ラベル x に固定した `ψ(x)` の direct path 値を λ² 重みで
足し合わせたものが、`contract_pair_MF` の一括 trace と一致すること。層0〜1 が Python、
層2 が 1 次元脚なので、この層で「奇ラベルの開いた脚を graded trace で閉じる」ことの
C++ 側の整合を確認する。

### 層4: iTPS レベル(`test_fermion_layer.cpp`、`iTPSTestAccessor` 経由)

- 全パリティ偶の `fermion = true` 状態で `MeanField_Env = true` の `measure_twosite` /
  `measure_onesite` が、同じ Tn・λ の `fermion = false` の結果と一致する(配線の回帰)。
  仮想脚の初期パリティは `even_first_parity`(偶奇混在、`tensors.cpp:81`)なので、
  `iTPSTestAccessor::finfo` 経由で `finfo.virt` も全偶に差し替える。
- フェルミオン状態で MF の norm がボソン MF の norm と一致する(§3.2)。
- 同一ボンドを両端から測った値(dx=+1 を A から、dx=−1 を B から、対称な演算子)が一致する。
- 2×2 セルで並進不変(既存の「translation invariant across wraps」と同型の MF 版)。
- MF の 1サイト密度が、窓を 1 サイトにした Fock オラクルの密度と一致する(§3.4 の担保)。
- `test/input.cpp`: `fermion = true` + `meanfield_env = true` が `validate_fermion_constraints`
  を通る(既存サブケースの反転)。他のガードが生きていることは既存サブケースが担保。

### 層5: E2E(`ctest`)

`test/fermion/free_fermion.py.in` に MF 版を追加する。spinless 自由フェルミオン半充填、
D=2、`meanfield_env = true` で、エネルギーが厳密解に対して **緩い許容(10% 程度、MF の
精度を反映)** に入り、かつ同じ Tn の CTM 版と符号が一致すること。実行時間は数秒
(CTM 環境を作らないので既存の E2E より短い)。

### 変異テスト(レビュアー向け)

- `contract_pair_MF` の `conj_f` を素の複素共役に差し替える → 層1 のホッピングが落ちること
- `wrap_twosite_gate` を `wrap_reduced_pair_op` に差し替える → 層2 が落ちること
- `build_pair_state` の graded tensordot を mptensor の素の tensordot に差し替える →
  層1 のホッピングが落ちること(norm は落ちない。§3.2 のとおり)

## 6. ドキュメント・NEWS

- `parameter_section.rst`(ja/en): fermion の非対応一覧から「平均場環境」を外す。
  `meanfield_env` の説明に「フェルミオン模式でも使える。2サイト観測量は単層縮約で
  評価され CTM 版より大幅に軽いが、精度は simple update 相当」と注記する。
- `NEWS.md`: 「フェルミオン模式で `meanfield_env = true` が使えるようになった。D>2 の
  2サイト観測量のコストが D¹²d⁴ から D⁶d² になる」旨。

## 7. タスク分割(実装計画の骨子)

| # | 内容 | 担当 |
|---|---|---|
| T1 | 層0: Fock オラクル拡張 + 自己検査 + 参照値生成 | テスト作成者(Claude が自己検査を独立確認) |
| T2 | 層1〜3 のテスト → `reduced_measure.hpp` の API 実装、`reduced.hpp` のリファクタ | テスト作成者 → Codex |
| T3 | 層4 のテスト + `input.cpp` 反転 → `twosite_obs.cpp` 配線、`load_toml.cpp` ガード削除 | テスト作成者 → Codex |
| T4 | 層5 E2E、docs(ja/en)、NEWS | テスト作成者 → Codex |

T1 はテスト基盤なので実装者には渡さない。T2 の RED 確認では「API 未定義によるコンパイル
エラー」ではなく、スタブ(例外を投げる)を置いて値の不一致で落ちることを確認する。

## 8. リスクと判断

- **オラクルと C++ が一致しない場合**: オラクルを真とし、まず層0 の自己検査と層3 の
  線形性で「どちら側の問題か」を切り分ける。graded 代数の順序不変性が壊れているなら
  `fops.hpp` の問題で本設計の範囲外(その場合は中断して報告)。
- **精度**: MF はボソンと同じく simple update 相当。利用者への注意はドキュメントに書く。
  CTM 版との使い分け(スキャンは MF、最終値は CTM)は従来どおり。
- **`save_density` の既存の脆さ**: 1サイト観測量が 1 つも無い入力では
  `onesite_operator_names[0]`(`density.cpp:127`)で落ちる。ボソン・CTM でも同じで本設計と
  無関係なのでスコープ外(別途修正候補として記録)。
- **CTM 経路のコスト**は本設計で解決しない。`apply_joint_swaps` の 1 パス化(数倍)と、
  不純物分解による blob 廃止(根本解)は別設計とする。
