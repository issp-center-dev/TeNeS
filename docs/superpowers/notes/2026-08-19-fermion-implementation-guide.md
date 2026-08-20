# フェルミオン対応 実装解説(レビュー用ガイド)

対象: `fermion` ブランチ(develop からの差分 33 コミット)。
コードレビューの下敷きとして、**アルゴリズムの考え方 → 実装の構造 → ファイル・コードの対応 → 各所の検証根拠**の順に説明する。
規約の確定経緯は設計書 `docs/superpowers/specs/2026-08-15-fermion-m1-design.md`(改訂1〜4)が正。

---

## 0. 何をやっているのか(3行)

フェルミオンの反交換符号は、テンソルネットワークでは「脚の交差」に局在する。
本実装は**脚ごとの Z₂ パリティを持つ薄いラッパー `ftensor` を導入し、縮約・転置・分解のたびに交差符号を機械生成する**。
既存のテンプレートカーネル(simple update)はそのまま `ftensor` でインスタンス化するだけで済み、ボゾン経路のコードと golden は一切変わらない。

---

## 1. アルゴリズムの原理

### 1.1 なぜ符号が出るか

第二量子化の状態は生成演算子の**順序付き積**で定義される:

```
|n₁ n₂ … n_N⟩ = (c†₁)^{n₁} (c†₂)^{n₂} … (c†_N)^{n_N} |0⟩
```

テンソルネットワークの縮約は、この積の中身を並べ替える操作に対応する。
奇パリティ(= フェルミオン数が奇数)の対象を2つ入れ替えるたびに −1 が付く。
これがすべての符号の源で、実装上は次の1本の規則に集約される:

> **軸の置換で順序が反転する脚対のうち、両方が奇である対の個数を k とすると符号は (−1)^k**

要素ごとに見れば、`idx` における脚 a の parity と脚 b の parity の積で判定できる(`ftensor.hpp` の `transpose_sign`)。

### 1.2 Z₂ 超選択則

物理的な状態・演算子は**全脚パリティの合計が偶**でなければならない(奇成分は非物理)。
実装ではこれを不変条件として扱い、`parity_violation()` で検査する。
初期化・分解・更新のすべてで偶が保たれることをテストで確認している。

### 1.3 なぜ「1サイト測定」と「2サイト測定」で扱いが違うか

CTM 環境は**閉ループ**を含む。開いたネットワーク(状態の縮約)では上記の局所規則で完結するが、
ループを閉じると supertrace(グレード付きトレース)の (−1) が混入し、素朴な機械簿記では規約が破綻する
(設計書 改訂1・2 に、この破綻を実測して arrow 方式を破棄した経緯がある)。

そこで M1 の最終方針は:

- **simple update**(開ネットワーク)→ `ftensor` で直接計算
- **測定・環境**(ループあり)→ サイトごとに bra/ket 2層を融合した **reduced tensor** をまず構築し、
  そこから先は**交差の無い通常のテンソルネットワーク**として既存のボゾン用 density-CTM をそのまま使う

つまり**フェルミオン符号を reduced tensor の構築関数1箇所に閉じ込める**。これが本実装の中心的な設計判断。

---

## 2. レイヤー構造とファイル対応

```
  入力層     load_toml.cpp        fermion フラグ・parity 読込・ガード
                 ↓
  簿記層     fermion_info.hpp     FermionInfo(物理/仮想脚パリティ台帳)+ wrap/unwrap
                 ↓
  符号層     parity.hpp           パリティ純関数(fuse / count_odd / parity_sort_perm)
             ftensor.hpp          ftensor 本体 + transpose_sign(符号の定義そのもの)
             fops.hpp             f-プリミティブ(tensordot/conj/qr/svd/svd_trunc/…)
                 ↓
  発展       core/simple_update.cpp   ftensor でインスタンス化(カーネル共通)
                 ↓
  測定       reduced.hpp          reduced tensor 構築(符号はここに集中)
             reduced_measure.hpp  reduced blob + density CTM の縮約関数
```

| ファイル | 行数 | 役割 |
|---|---|---|
| `src/fermion/parity.hpp` | 73 | `parity_vector` / `leg_parities` 型と純関数3つ |
| `src/fermion/ftensor.hpp` | 117 | `ftensor<tensor>` 構造体、member `transpose`(符号適用の本体)、`transpose_sign` |
| `src/fermion/fops.hpp` | 706 | f-プリミティブ全部。mptensor と**同名**の自由関数でオーバーロード解決 |
| `src/fermion/fermion_info.hpp` | 95 | `FermionInfo`(phys/virt パリティ)、`wrap_Tn`/`unwrap_Tn`、隣接整合検証 |
| `src/fermion/reduced.hpp` | 212 | reduced tensor(1サイト・2サイトblob)の構築 |
| `src/fermion/reduced_measure.hpp` | 154 | λ dressing、reduced 群の生成、blob と CTM 環境の縮約 |

---

## 3. 符号層の詳細(レビュー要点つき)

### 3.1 `ftensor`(`ftensor.hpp`)

```cpp
template <class tensor>
struct ftensor {
  tensor t;              // 密テンソル本体(mptensor 型、無改造)
  leg_parities parity;   // parity[axis][index] = その添字のパリティ
};
```

パリティを「セクター次元 (D_e, D_o)」ではなく**添字ごとの bool ベクトル**で持つのが要点。
reshape(脚の融合)後の非連続パターンや、CTM の `extend` によるパディングでも表現できる。

member `transpose` は「**符号マスクを掛けてから** `mptensor::transpose`」の順で、parity も同時に並べ替える。

> **レビュー観点**: member と自由関数 `transpose` が同じ規約であること(自由関数は member に委譲)。

### 3.2 graded tensordot(`fops.hpp:289`)

意味論は「A の縮約脚を末尾へ、B の縮約脚を**逆順で**先頭へ寄せる」置換の符号 + 素の tensordot。
実装は置換を実行せず**符号マスクだけ**を掛けて `mptensor::tensordot` に渡す(等価かつ高速):

```cpp
a_masked = apply_transpose_sign_mask(a, tensordot_left_perm(...));   // (自由脚…, 縮約脚…)
b_masked = apply_transpose_sign_mask(b, tensordot_right_perm(...));  // (縮約脚逆順…, 自由脚…)
ret.t = mptensor::tensordot(a_masked.t, b_masked.t, axes_a, axes_b);
```

縮約脚対のパリティ一致は `validate_contracted_parity` が実行時検証する(規約ミスの早期検出)。

> **レビュー観点**: `tensordot_right_perm` の逆順が B 側にのみ適用されていること。ここが規約の核。

### 3.3 `conj`

**グレード付きのエルミート共役**を、脚順序は元のまま(ラベルを戻した形)で表したもの。
定義は「要素共役 ∘ 全脚の順序反転 ∘ 元の軸番号への転置」で、この合成の正味が要素ごとの
(−1)^{m(m−1)/2}(m = その要素で奇な脚の数)になる。符号は**導出されたもの**で、
m 個の奇脚を反転するのに要する互換の数 C(m,2) がそのまま指数。

**JW 厳密参照テストが確定するのは「規約の組み合わせ」**である点に注意。
twist を bra 側に置くか ket 側に置くか、`tensordot` の逆順を A 側にするか B 側にするかは
個別には任意で、物理的に意味を持つのは組み合わせだけ。Task 9 のテストは組み立てた全体が
厳密解を再現することを保証している。

### 3.4 graded QR / SVD(`fops.hpp:410` / `:503` / `:600`)

手順(QR も SVD も共通):

1. `transpose(a, rows + cols)` で行・列に整列(graded、符号込み)
2. 行・列を融合し `fuse` でパリティを計算
3. `parity_sort_perm` の**置換行列との tensordot** で偶先頭に整列(分散安全。直接インデックス書き込みをしない)
4. 偶ブロック・奇ブロックを `slice` で切り出し、各々 `mptensor::qr` / `svd`
5. 直和で組み直し、置換を戻す
6. 内部脚のパリティ(偶 k_e 個 + 奇 k_o 個)を返す

`svd_trunc`(`:600`)はさらに**セクター横断の切断**を行う: 全特異値を降順(同値なら偶優先)で dc 本選び、
選択後に偶先頭へ再ソート。選択は**選択行列との tensordot** で実装(同じく分散安全)。

**行と列が交互になる分割**(例 rows=(0,2), cols=(1,3))は、規約上まったく問題ない。
`svd(a, rows, cols)` は 1 行目で `transpose(a, rows + cols)` を行うため、interleaved 指定は
「先に graded transpose してから連続分割する」のと**同一の計算**である
(テスト `graded svd is invariant under regrouping to a contiguous split` が固定)。
かつて非交差化が効いたように見えたのは、同時に導入した regroup が**符号を落とす転置**で、
その欠落が当時のゲート符号の欠落と相殺していたため。両方修正した現在、カーネルは
interleaved の軸指定を `svd_trunc` にそのまま渡している。

> **レビュー観点**: 分散(ScaLAPACK)対応のため、行・列の並べ替えと切断列の選択は
> 要素の直接書き込みではなく**置換行列・選択行列との tensordot** で実装している。
> `make_perm_matrix` の規約は `P[new, old] = δ`。
> また `qr`/`svd` は偶偶・奇奇の対角ブロックしか見ないため、パリティ偶でない入力を渡すと
> 非対角成分が黙って捨てられる。Debug ビルドでは `validate_block_diagonal` が検出する。

---

## 4. simple update(発展)

### 4.1 カーネル(`src/iTPS/core/simple_update.cpp`)

**ボゾンと完全に同一のソース**を `ftensor` でインスタンス化する。カーネル側の変更は1点のみ:

- 「svd → slice で切断」を `svd_trunc(...)` に置換(ボゾン用オーバーロードを `src/tensor.hpp` に追加して共通化)

θ の脚順序は (aux1, aux2, out1, out2) で、二分割はサイト1 = (aux1, out1) 対 サイト2 = (aux2, out2)。
この軸指定を `svd_trunc` にそのまま渡す(§3.4 のとおり、必要な graded transpose は分解側が行う)。
分解後に `enforce_even_parity` で超選択則を検査(閾値超過は例外、微小な数値ゴミは 0 にクリップ)。

### 4.2 ドライバ境界(`src/iTPS/simple_update.cpp:60-`)

`finfo.enabled` の分岐で wrap → カーネル → unwrap を行う。ここで**2つの規約**が入る:

```cpp
// (a) ボンド方位の正規化: canonical はラスタ順(左→右、上→下)
if (source_leg == 0 || source_leg == 1) {   // 右→左 / 下→上 で来た場合
  std::swap(s1, s2); std::swap(s1_leg, s2_leg);
  op12 = mptensor::transpose(up.op, Axes(1, 0, 3, 2));
}
// (b) ゲートのロード: graded tensordot の縮約脚逆順規約を打ち消す
auto fop = wrap_twosite_gate(op12, finfo.phys[s1], finfo.phys[s2]);
```

- (a) がないと、JW 順で後のサイトからゲートを当てることになり誤ったマスクが掛かる
  (症状: 2D で市松 CDW への崩壊、D=4 で負ノルム)
- (b) がないと `|11⟩⟨11|` チャネル(全 Trotter ゲートが持つ)が反転する
  (症状: τ² スケールの λ 乖離。ホッピングのみの参照テストでは検出できなかった)

**(b) は追加の物理ではない。** ユーザーが与える `op12` の行列要素には既に反交換関係が
入っている(順序付き Fock 基底の定義)。ここで足す swap は、graded tensordot が
**B 側(演算子)の縮約脚を逆順に寄せる**ために生じるマスクを打ち消すためのもの:

- `tensordot(Theta_before, op12, Axes(1,3), Axes(0,1))` に対し
  `tensordot_right_perm(4, (0,1)) = (1,0,2,3)` → マスクは「in1, in2 が共に奇の要素で −1」
- `wrap_twosite_gate` の `apply_swap(fop, 0, 1)` は**同一のマスク**

両者は厳密に相殺し、正味では**書いたとおりの行列**が θ に作用する。つまり
「行列要素という意味論」と「graded テンソルという意味論」の間のアダプタである。
一方、同じ tensordot が A 側に掛けるマスク (−1)^{p_{aux2} p_{phys1}}(サイト2の補助脚が
サイト1の物理脚を横切る)は**本物の幾何的交差**であり、これは残る。

> **レビュー観点**: `wrap_twosite_gate` は入力脚のみ、測定側の `wrap_reduced_pair_op` は
> 入力・出力の両方。後者では演算子の出力脚も bra 層に対して閉じられるため。§5.3 参照。

---

## 5. 測定(reduced tensor + density CTM)

### 5.1 考え方

サイトごとに bra 層と ket 層を融合した4本脚テンソル a[(l l̄), (t t̄), (r r̄), (b b̄)] を作る。
この粗視化格子の上では**交差が存在しない**ので、以降は既存のボゾン用 density-CTM 一式
(`Calc_CTM_Environment_density` + `contract_density_ctm/`)を**無改造で**使える。

### 5.2 構築パイプライン(`reduced.hpp`)

#### 1サイト: `doubled_pipeline`(`:92`)

共通の心臓部。ket 側テンソル `T[l,t,r,b,s]` と bra 側 `T̄` から4段階で作る:

| 段階 | 操作 | 脚(rank) |
|---|---|---|
| 0 | 入力 | bra `(l̄,t̄,r̄,b̄,s̄)`、ket `(l,t,r,b,s)` |
| 1 | `tensordot(conj(bra), ket, Axes(), Axes())`(外積) | `(l̄,t̄,r̄,b̄,s̄, l,t,r,b,s)` = axes 0–4 が bra、5–9 が ket(rank 10) |
| 2 | `apply_joint_swaps(doubled, {0,1,2,3}, {5,6,7,8}, {0,1,2,3})` | 形は不変。bra/ket 層をまたぐ交差の符号(= swap gate)を要素マスクで挿入 |
| 3 | `transpose` で (ket, bra) インターリーブ | `(l,l̄, t,t̄, r,r̄, b,b̄, s, s̄)`。**graded transpose なので追加の Koszul 符号もここで発生** |
| 4 | `mptensor::reshape` で隣接ペアを融合 | `([l l̄], [t t̄], [r r̄], [b b̄], s_ket, s_bra)`(rank 6) |

段階 2–3 で全てのフェルミオン符号が数値に焼き込まれるため、**出力は普通の(符号を持たない)テンソル**
になる。だからこれ以降の CTM はボゾン用実装をそのまま使える。

`apply_joint_swaps` は `kDoubledJointMask` が指定する脚ペア (x,y) について
`apply_swap(a, ket_axes[ix], bra_axes[iy])` と `apply_swap(a, bra_axes[ix], bra_axes[iy])` を掛ける。
このマスク(現在は (0,3), (1,0), (2,3), (3,0))は**解析的に導出したものではなく、
Fock oracle を審判にした探索で同定した**もの(コメントに明記)。
YASTN の `fuse_layers()` 規約とゲージ同値であることを確認している。**レビュー最重要点**。

用途別の入口:

| 関数 | 実装 | 用途 |
|---|---|---|
| `build_reduced_op`(`:160`) | `doubled_pipeline(Tn, Tn)` | 物理脚 (s_ket, s_bra) を開いたまま。1サイト演算子の挿入用 |
| `build_reduced`(`:168`) | 上を `contract(Axes(4), Axes(5))` | 物理脚をトレース。ノルム・環境構築用 |
| `build_reduced_pair`(`:174`) | 下記 `doubled_cluster` | 2サイト演算子用の 6 本脚 blob |

#### 2サイト: `doubled_cluster`

2サイトをまとめて doubling する。ペアにまたがる演算子は**ket 層に挿入してから**二重化する
(奇演算子は仮想脚を通過する際に反交換するため、二重化後の挿入では表現できない):

| 段階 | 操作 | 脚 |
|---|---|---|
| 1 | `tensordot(TnA, TnB, …)` で ket 側クラスタ | 横なら A の右脚と B の左脚を縮約 → 外部 6 脚 + 物理 2 脚(rank 8)|
| 2 | `apply_reduced_two_site_op`(`:69`) | `tensordot(psi, op12, Axes(3,7), Axes(0,1))` で両物理脚に演算子を作用 → out 脚を元位置へ `transpose` |
| 3 | `tensordot(conj(ket_ab), ket_op, Axes(), Axes())` | 外積(rank 16)|
| 4 | `fuse_doubled_cluster`(`:120`) | joint swap → インターリーブ → 6 本の外部脚を融合 → 物理ペアを `contract` で閉じる |

**演算子は段階 2 で in 脚、段階 4 で out 脚がそれぞれ閉じられる**。これが測定側で
入力・出力の**両方**に補正が要る理由(§5.3)。最後にサイト間ボンドの融合ゲージを
`apply_fused_leg_gauge` で 1 サイト版の規約に合わせる(これも oracle 固定)。

### 5.3 演算子ロードの3規約(**最重要の注意点**)

| 経路 | 規約 | ヘルパ | 根拠 |
|---|---|---|---|
| simple update ゲート | 入力脚 swap | `wrap_twosite_gate` | d=2 判別テスト + λ軌跡 + Fock 審判 |
| 2サイト測定 blob | 入力 + 出力 swap | `wrap_reduced_pair_op` | R3(d=2)+ **R5(d=4)** oracle |
| 1サイト(ゲート・測定) | 素 | — | 入力脚の交差が存在しない |

d=2 では blob の「素ロード」と「両swap」が**N保存演算子上で縮退**するため、当初 R3 は「素」と誤って固定していた。
d=4(電子系)の (odd,odd)→(even,even) チャネル(スピンありホッピングの線形項)で初めて分岐し、R5 で修正した。

> **レビュー観点**: この3規約が混ざると符号が黙って狂う。新しい経路を足すときは必ず oracle テストを先に書く。

### 5.4 測定ドライバ(`src/iTPS/measure.cpp`, `onesite_obs.cpp`, `twosite_obs.cpp`)

- `measure.cpp:45-` : fermion 時、**bare Tn**(λ dressing なし)から reduced 群を作って density CTM を回す
  - TeNeS 規約ではカーネルが √Schmidt を bond 両端の Tn に書き込むため、状態は bare Tn の直接縮約
  - full-λ dressing は MeanField 経路の規約。CTM に渡すと環境重みの二重計上になる(実際にバグだった)
- `onesite_obs.cpp:87-` : reduced impurity + `Contract_one_site_density_CTM`
- `twosite_obs.cpp:174-`(ノルム)、`:217-`(縦ボンド)、`:259-`(横ボンド):
  blob + 専用の窓関数(`contract_reduced_pair_{horizontal,vertical}_density_CTM`)

---

## 6. 入力層とガード(`load_toml.cpp`)

- `[parameter.general] fermion`、`[[tensor.unitcell]] parity` の読込
- `validate_fermion_constraints`(`:602-`)が M1 の非対応機能を**入力時に**まとめて拒否
  (非 ground state / full update / MeanField / Simple_Gauge_Fix / RSVD / 相関関数 / multisite /
  距離2以上の2サイト / tensor_save・load / パリティ奇の演算子 / 奇プロダクト初期状態)
- 相関長のみ「デフォルト有効」なので、エラーではなく**強制無効化 + 警告**

> **レビュー観点**: 「沈黙して間違った数値を出さない」が原則。ガードの網羅性は仕様書 §6.1 の表と対照。

---

## 7. テストの地図(何がどの層を保証しているか)

### 層1: 単体(`test/test_fermion_layer.cpp` 前半、行番号は現在の値)

| 対象 | テスト |
|---|---|
| parity 純関数 | `fuse …column-major with XOR`(:95)ほか |
| マスク演算 | `apply_swap negates odd-odd elements only`(:115) |
| graded transpose | `applies Koszul signs`(:154)、`transpose twice returns original`(:170) |
| tensordot | `rejects parity mismatch`(:197)、`contraction-order independent`(:207)、`matches manual swap`(:228) |
| conj | `norm … positive via conj`(:260)、`involution`(:279) |
| QR/SVD | `graded QR reconstructs`(:418)、`graded SVD reconstructs`(:445)、`svd_trunc selects across sectors`(:468)、`rejects interleaved row-column splits`(:736) |

### 層2: 規約固定(厳密参照との対査)

| テスト | 何を固定するか |
|---|---|
| `JW four-site reference matches f-primitive contractions`(:588) | `conj` と tensordot の規約(JW 厳密値、rtol 1e-12) |
| `two-site operator doubly-odd input channel …`(:636) | `wrap_twosite_op`(ゲート規約) |
| `manual four-swap reduced tensor …`(:672) | 文献規約(Corboz)との橋渡し |
| `r2_convention.cpp` R2/R3 | Fock oracle vs f-プリミティブ vs reduced(d=2、3者一致) |
| `r2_convention.cpp` **R5**(:619) | **d=4 の blob 演算子規約**(Fock 検証済み直接経路が審判) |

### 層3: 統合(常時実行)

- `reduced density CTM measurement gives positive fermion norms`(:764)
- `fermion twosite measurement is translation invariant across wraps`(:3027)
- ctest `FreeFermion`: E2E(既定は D=2 のみの quick モード、48秒)

### 層4: 診断(環境変数ゲート、常時は skip)

| 環境変数 | 用途 |
|---|---|
| `TENES_RUN_LAMBDA_TRAJECTORY_DIAG` | 横鎖 λ 軌跡を JW 双子ボゾンと比較(case=A/B 判定)|
| `TENES_RUN_VERTICAL_LAMBDA_DIAG` | 同・縦鎖(leg 1/3 で向き規約を判別)|
| `TENES_RUN_PLAQUETTE_TROTTER_DIAG` | 開放 2×2 パッチでカーネル vs 厳密 Trotter vs Fock(3者)|
| `TENES_RUN_WEAK2D_DIAG` | 弱2D結合で CTM 測定 vs CTM 非依存の平均場推定(測定バグの切り分け)|
| `TENES_RUN_ELECTRON_CHAIN_DIAG` | d=4 電子鎖(θ 成分ダンプ、Fock 審判用)|

> **注意**: JW 双子 λ 判定は **d=2 でのみ有効**。d=4 では奇セクター多重度が 2 になり成分の素コピーが
> JW 対応でなくなるため、case=B が出てもバグの証拠にならない(設計書 改訂4)。

---

## 8. レビューで特に見てほしい箇所

1. **`reduced.hpp` の `kDoubledJointMask` と `apply_fused_leg_gauge`** — 唯一「解析導出でなく実験的に固定した」部分。
   oracle テストで縛ってはいるが、導出を与えられるならそれが理想
2. **演算子ロードの3規約**(§5.3)— 経路ごとに違うのは設計として気持ち悪い。統一できないか
3. **`svd_trunc` の切断規約** — セクター横断の降順選択が変分的に正しいか(縮退時の偶優先は決定論性のため)
4. **`validate_noninterleaved_split` の責任分界** — ガードで拒否 + 呼び出し側が regroup、で妥当か
5. **性能**: `doubled_cluster` が D¹² スケーリング(D=4 の2サイト測定に6時間超)。
   改善案は設計書の性能課題節に記載(bra 要素ループ + ket ベクトル演算で ~100× 見込み)
6. **`iTPS.hpp` の `finfo` メンバと分岐の入り方** — ボゾン経路への影響がゼロであることの確認

---

## 9. 既知の制限(M1 スコープ外)

- expert モードのみ(`tenes_simple`/`tenes_std` は未対応)
- 距離2以上の観測量・相関関数・相関長・有限温度・実時間発展・full update・MeanField 環境
- `tensor_save`/`tensor_load`(FermionInfo の直列化が未実装)
- n ≠ 0.5(μ≠0)では simple update の収束が大幅に遅い(必要虚時間 ~45 vs 半充填 ~10)
- d=4(電子系)は規約固定済み・1D 極限検証済みだが、**2D の量的検証は未実施**
