# tenes_simple / tenes_std のフェルミオン対応 設計書

日付: 2026-08-20
ブランチ: `fermion`
前提: `docs/superpowers/specs/2026-08-15-fermion-m1-design.md`(改訂1〜4)で確定した
C++ 側の規約。特に**改訂4の「演算子規約の最終表」**と、
`docs/superpowers/notes/2026-08-19-fermion-implementation-guide.md`(実装解説)。

改訂1 (2026-08-20): Codex レビューと Claude の独立確認を反映。主な変更は
(a) `tenes_std` の長距離 MPO 分解を fermion で拒否(§6.2 — 黙って間違う唯一の経路)、
(b) 「やらないこと」を実際に止めるガード節を新設(§6.3)、
(c) 明示 rank-4 2サイト観測量のチャネルを追加(§4.3 — 既存機構では `hopping` を表現できない)、
(d) `parity` はパース時点で落ちている(§6.2)、
(e) d = 4 専用のテスト層を追加(§8 層1b)、
(f) 利用者向けメッセージに `M1` を出さない文言規約(§6.3)。

## 1. 目的とスコープ

M1 でフェルミオンの C++ 経路(simple update + reduced-tensor 密度 CTM 測定)が動くように
なったが、入力は expert モード(`input.toml` 直書き)に限られ、Trotter ゲート行列を
利用者が自分で構成する必要がある。本設計は `tenes_simple` / `tenes_std` を対応させ、
**3段パイプラインでフェルミオン模型を解けるようにする**。

### やること

- `tenes_simple` に2つのフェルミオン模型を追加(**正方格子・最近接のみ**)
  - `type = "spinless fermion"`: physical_dim 2
  - `type = "hubbard"`: physical_dim 4(スピンあり電子)
- ボンドハミルトニアン・2サイト観測量を**順序付き Fock 基底で JW 符号込みに構成**
- `[parameter.general] fermion = true` と `[[tensor.unitcell]] parity` の出力
- `tenes_std` の `parity` 保持・通し(現在は**パース時点で落ちている**)と fermion 時の入力検証
- 「やらないこと」を**実際に止めるガード**の実装(方針を書くだけでは足りない。§6.3)
- C++ 側の穴埋め: fermion + `ops = [i,j]` 形式の2サイト観測量をロード時に拒否

### やらないこと(M1 のガードどおりエラー)

- 正方格子以外(三角・蜂の巣・カゴメ)、距離2以上のボンド
- 相関関数・相関長・有限温度・実時間発展・full update・平均場環境
- パリティ奇のプロダクト初期状態(`|1⟩`、Néel `|↑⟩/|↓⟩` 等)。仮想ボンドの奇セクターを
  使うダイマー被覆構成が必要で M1/M2 の課題(前設計書 §9)
- t-J 模型(拘束空間の扱いが別課題)

## 2. 全体方針

フェルミオン符号の責任分界は M1 で決めた線をそのまま守る:

| 層 | 責任 |
|---|---|
| `tenes_simple` | **順序付き2サイト Fock 基底での行列要素**を作る。すなわちサイト内(スピン間)と、そのボンドを構成する2サイト間の反交換符号 |
| `tenes_std` | 行列に対する `expm(-τh)`。**符号の知識は持たない**(通し役) |
| C++ (`tenes`) | そのボンドが2次元ネットワークに埋め込まれる際に生じる符号(他の脚との交差)を graded 機構で生成 |

すなわちツール側の新規性は「Fock 構成器」1点に集約され、それ以外はスキーマの受け渡しに
過ぎない。この構成なら `tenes_std` の数値ロジックは無変更で済む。

### 却下した代案

- **積形式 + 符号補正表**: 既存の `ops = [i,j]`(1サイト演算子の直積)機構に符号表を
  後付けする案。積形式が誤りであることこそ M1 で根絶したバグクラスであり、補正表は
  同じ誤りを別の場所に再導入する。却下。
- **tenes_std でフェルミオン化**: `std.toml` のハミルトニアンは既に構造を失った行列で、
  どの部分がフェルミオン的かを `tenes_std` は知り得ない。「書いた行列がそのまま作用する」
  という M1 の規約も壊れる。却下。

## 3. Fock 構成器(共有モジュール関数)

`tool/tenes_simple.py` にモジュール関数として置く。**符号を生む唯一の場所**であり、
テストもここに集中させる。

```python
def fermion_modes(nspin: int) -> int:
    """2サイト系のモード数。順序は (site1, spin...), (site2, spin...)。"""
    return 2 * nspin

def fock_cop(dagger: bool, mode: int, nmodes: int) -> np.ndarray:
    """JW string 込みの生成/消滅演算子(2^nmodes 次元)。

    基底は占有数ビット列 g で、|n_0 ... n_{M-1}> = (c†_0)^{n_0} ... (c†_{M-1})^{n_{M-1}}|0>。
    符号は先行モードの占有数の偶奇 (-1)^{sum_{k<mode} n_k}。

    軸慣習は NumPy 標準の **`mat[out, in]`**(すなわち `mat @ vec` が作用)。
    これは `tenes_simple` の1サイト演算子(`Sz` 等)および `model_sitehamiltonian` が
    `ham[input, output]` 慣習なのと**逆**なので、`bond_matrix` / `onesite_matrix` で
    明示的に転置して吸収する(下記)。取り違えは符号ではなく転置として現れ、
    エルミートな観測量では見えないため、テストで固定する(§8 層1)。
    """

def local_index_to_occupation(i: int, nspin: int) -> List[int]:
    """局所基底 i を (n_up[, n_dn]) に。spinless: i = n。
    spinful: i = n_up + 2 n_dn、すなわち |0>, |up>, |dn>, |up dn>。
    サイト内順序は |up dn> = c†_up c†_dn |0> に固定(file_specification と同一規約)。
    """

def bond_matrix(fock_op: np.ndarray, nspin: int) -> np.ndarray:
    """2サイト Fock 空間の行列を rank-4 の op[in1, in2, out1, out2] に並べ替える。

    op[i1,i2,o1,o2] = <o1 o2| O |i1 i2> = fock_op[out_global, in_global]、
    ここで in_global / out_global は §3 の占有数ビット列で、
    局所添字 (i1,i2) / (o1,o2) と `local_index_to_occupation` で相互変換する。

    脚順 `[in1, in2, out1, out2]` は `tenes_simple` の既存2サイト出力
    (`np.einsum("ij,kl -> ikjl", oi, oj)`、`oi[in,out]`)と同一であり、
    C++ の `wrap_twosite_gate` / `wrap_reduced_pair_op` が受ける
    「素の行列要素 <out1 out2|O|in1 in2>」とも一致する(swap は C++ 側が足す)。
    """

def onesite_matrix(fock_op_1site: np.ndarray, nspin: int) -> np.ndarray:
    """1サイト演算子(必ずパリティ偶)を rank-2 `op[in, out]` に。"""
```

`bond_matrix` に渡す `fock_op` は `fock_cop` の積和で組む。例(spinful ホッピング):

```python
h = 0
for s in range(nspin):
    m1, m2 = s, nspin + s          # site1 spin s, site2 spin s
    h += -t * (cd(m1) @ c(m2) + cd(m2) @ c(m1))
```

**符号を手で書くコードは一切書かない。**

## 4. 模型クラス

`SpinModel` / `BoseHubbardModel` と同じく `Model` を継承する薄い2クラス。共有の危険部分は
§3 のモジュール関数にあるため、クラス側は分けたほうが読みやすい(パラメータ集合・観測量・
初期状態・エラーメッセージが模型ごとに異なるため、1クラスに畳むと `if nspin == 2:` が
散在する)。

### 4.1 `SpinlessFermionModel`

- `physical_dim = 2`、局所基底 `|0>, |1>`、`parity = [0, 1]`
- パラメータ: `t`(ホッピング)、`v`(最近接斥力)、`mu`(化学ポテンシャル)
  - 近傍レベル・ボンドタイプ記法(`t1`, `t'` 等)は `BoseHubbardModel.read_params` を踏襲。
    ただし M1 では最近接のみ有効で、2近傍以上が非ゼロなら明示エラー
- ボンドハミルトニアン(z = 配位数、onsite 項は z 等分):

  H_bond = −t (c†₁c₂ + h.c.) + V n₁n₂ − (μ/z)(n₁ + n₂)

- 1サイト観測量: `n`
- 2サイト観測量: `hopping` = ⟨c†₁c₂ + h.c.⟩(明示 rank-4)、`nn` = ⟨n₁n₂⟩(直積)
  - `bond_hamiltonian` は `tenes_simple()` が group 0 として**自動で出力する**ので
    模型側では持たない(現行実装 `tool/tenes_simple.py` の `for ham in hams:` ループ。
    既に `elements` 形式なのでフェルミオンでもそのまま正しい)

### 4.2 `HubbardModel`

- `physical_dim = 4`、局所基底 `|0>, |up>, |dn>, |up dn>`、`parity = [0, 1, 1, 0]`
- パラメータ: `t`、`u`(オンサイト斥力)、`v`(最近接斥力)、`mu`、`h`(Zeeman、任意)
- ボンドハミルトニアン:

  H_bond = −t Σ_σ (c†₁σc₂σ + h.c.) + V n₁n₂
           + (1/z)[ U(n↑n↓)₁ + U(n↑n↓)₂ − μ(n₁+n₂) − h(S^z₁+S^z₂) ]

- 1サイト観測量: `n`, `n_up`, `n_dn`, `Sz`, `doublon` = n↑n↓, `holon` = (1−n↑)(1−n↓)
- 2サイト観測量: `hopping`(明示 rank-4)、`nn`、`SzSz`(いずれも後2つは直積)。
  `bond_hamiltonian` は §4.1 同様に自動出力

### 4.3 明示 rank-4 2サイト観測量のチャネル(既存機構の拡張)

現行の `Model.twosite_ops: List[Tuple[int, int]]` は**1サイト演算子の添字ペア**しか持てず、
出力側(`tool/tenes_simple.py` の `for i, j in model.twosite_ops:`)は
`np.einsum("ij,kl -> ikjl", oi, oj)` すなわち**直積しか作れない**。
`hopping` = ⟨c†₁c₂ + h.c.⟩ は直積に分解できないため、この機構では表現できない。

そこで `Model` に並列のチャネルを追加する:

- `twosite_ops_explicit: List[Tuple[str, np.ndarray]]` — 名前と rank-4 配列
  `op[in1, in2, out1, out2]` の組。既定は空リスト

`tenes_simple()` の2サイト観測量出力ループの直後に、このリストを走査して
`elements` 形式で書き出す分岐を足す(group 番号は既存カウンタ `k` を継続)。
既存模型は空リストのままなので**挙動不変**。

`nn` / `SzSz` のような真の直積は既存の `twosite_ops` に載せてよいが、フェルミオンでは
§4.4 の理由で結局 `elements` に落とす。

### 4.4 共通の制約

- **1サイト観測量はパリティ偶のもののみ出力する**。`c†`, `c` 単体は奇で、C++ 側の
  ロードガードに弾かれるため出力しない(相関関数も M1 では非対応)
- **2サイト観測量は必ず明示 `elements` で出力する**。`ops = [i, j]` 形式(1サイト演算子の
  直積)は fermion モードで未対応(§6.2、§7)。したがってフェルミオン模型では
  `is_complex` 分岐で `ops` 形式を選ぶ経路に入らないよう、模型が fermion のときは
  常に `elements` 側へ倒す

## 5. 初期状態

既存の `[lattice] initial` キーを使い、フェルミオン模型では値を差し替える。
`"ferro"` / `"antiferro"` は**受け付けない** — 実際に作れるのは電荷秩序であり、
特に `"antiferro"` は SDW を想起させて誤解を招くため。

| `initial` | spinless fermion | hubbard |
|---|---|---|
| `"random"`(既定) | 乱数偶初期化(現行どおり `initial_state = [0.0]` を出力) | 同左 |
| `"vacuum"` | 全 `\|0>` | 全 `\|0>` |
| `"full"` | **エラー** | 全 `\|up dn>`(n = 2) |
| `"cdw"` | **エラー** | `\|0>` / `\|up dn>` の市松(n = 1) |
| `"ferro"`, `"antiferro"` | エラー | エラー |

`"vacuum"` は `noise` が全要素(偶マスク付き)に乗るため発展の不動点にはならず、
低密度側から出発する実用的な初期値になる。

### 5.1 モードをモデルへ渡す仕組み

現在の `tenes_simple()` は `lattice.initial_states` の文字列を見て
「`"random"` ならゼロ配列、`"ferro"` なら `st[0,:]`、それ以外は `st[i,:]`」と分岐しており、
モデル側はモードを知らない(`initial_states(num_sublattice)` がパターンを返すだけ)。
フェルミオン模型はモードごとに異なるパターンと固有のエラーが必要なので、`Model` に

- `initial_state_vectors(mode: str, num_sublattice: int) -> Optional[np.ndarray]`

を追加する。戻り値 `None` は「乱数初期化(ゼロ配列を出力)」を意味する。既定実装は
現行の挙動(`"random"` → `None`、`"ferro"` → 全副格子に `st[0]`、それ以外 → `st[i]`)を
そのまま再現し、`tenes_simple()` の分岐をこの1呼び出しに置き換える。
`SpinModel` / `BoseHubbardModel` は既定実装のままで**挙動不変**。フェルミオン2模型は
これを上書きし、§5 の表に無い値には理由付きの `RuntimeError` を投げる。

**添字の規約(現行コードの暗黙前提を明示化する)**: 戻り値の第0軸は
**非 vacancy 副格子を出現順に詰めた添字**とする。現行の `tenes_simple()` は
`num_sublattice` を非 vacancy のみで数えておきながら、書き出しループでは
全副格子を走る `i` で `st[i, :]` を引いており、vacancy が末尾でない格子では
添字がずれる(現状の格子実装ではたまたま vacancy が末尾なので露見していない)。
フックの導入に合わせて、書き出し側は**非 vacancy 用の別カウンタ**を進めて参照する。
この修正は既存格子では出力が変わらないので**挙動不変**である。

spinless の `"full"` / `"cdw"` と Hubbard の Néel 型が作れないのは、`|1>` や `|up>` が
パリティ奇で、TeNeS のプロダクト初期化が仮想脚の添字 0(偶)に状態ベクトルを置くため
サイトテンソルの全脚パリティ合計が奇になるから。エラーメッセージにこの理由と回避策
(`"random"`)を書く。

## 6. スキーマと出力

### 6.1 `tenes_simple` → `std.toml`

- `[parameter.general] fermion = true` を模型種別から**自動注入**する(利用者は書かない)。
  利用者が明示的に `fermion = false` と書いていたら、矛盾として明示エラー
- `[[tensor.unitcell]]` に `parity = [...]` を追加出力。vacancy 副格子(`physical_dim = 1`)は `parity = [0]`(偶)
- 2サイト観測量は `ops = [i,j]` を使わず `elements` を出力(`is_complex` 分岐に関わらず)
- `[correlation]` / `[correlation_length]` が書かれていたら明示エラー。前者は C++ が
  ロード時に弾くので二重ガードだが、後者は C++ では**警告して黙って無効化される**だけなので、
  ツール層で止めないと「指定したのに出ない」になる

### 6.2 `tenes_std` → `input.toml`

`parity` は**書き出しだけでなくパース時点で落ちている**。`LocalTensor.__init__` は
`physical_dim` と `virtual_dim` しか読まないので、修正は2箇所になる。

- `LocalTensor` に `parity` フィールドを追加してパースし、`check()` で
  「長さ == `physical_dim`」「各要素が 0 か 1」を検証する
- `[[tensor.unitcell]]` の書き出しに `parity` を追加
- `[parameter]` は既にテーブルごと素通しなので `fermion` は自動的に伝わる
- fermion 時の検証を、**時間発展演算子を生成する前**に実施する(生成後では手遅れ):
  - 全 unitcell に `parity` があること(欠けていたら明示エラー)
  - `[[observable.twosite]]` が `ops` 形式なら明示エラー(理由と `elements` 形式への
    案内を含める)
  - **`graph.make_path(bond)` の長さが 1 でないボンド項があれば明示エラー**(下記)
- **`expm(-τh)` の計算ロジックは変更しない**

#### 長距離ボンドの MPO 分解を止める(最重要)

`make_evolution_twosite()` は `nhops > 1` のとき、ゲートに中間サイトの恒等演算子を
テンソル積してから SVD で最近接ゲート列へ分解する経路を持つ。この分解は
**中間サイトに JW string を置かない**ため、フェルミオンでは誤りである。しかも分解後の
各ゲートは最近接なので、C++ 側の「距離2以上」ガードには到達しない。つまり
**黙って間違った数値が出る唯一の経路**であり、`tenes_std` で止めるほかない。

`tenes_simple` は正方格子・最近接のみを出すのでこの経路には入らないが、
`tenes_std` は `std.toml` を直書きする利用者にとって独立の入口であるため、
ここでのガードが必須になる。

### 6.3 「やらないこと」を実際に止めるガード

§1 の「やらないこと」は方針の記述であって、現行コードには対応する拒否が無い。
`make_lattice()` は honeycomb / triangular / kagome をそのまま構築し、
`hamiltonians()` は `params_twosite` に入った**全 neighbor level** のボンドを出力する。
フェルミオン模型を足すだけでは、非正方格子や `t'` / `v'` が黙って `std.toml` へ流れる。

`Model` に `is_fermion: bool = False` を持たせ(フェルミオン2模型で `True`)、
`tenes_simple()` の**格子構築直後**に:

- `SquareLattice` 以外なら明示エラー(`lattice.type` ではなくクラスで判定する。
  `"square lattice"` 以外の名前で `SquareLattice` に落ちる経路を取りこぼさないため)
- `params_twosite[n][t]` の `n >= 1`(2近傍以上)に非ゼロ値があれば明示エラー
- `[correlation]` / `[correlation_length]` があれば明示エラー(§6.1)

エラーメッセージには「現行版のフェルミオン対応は正方格子・最近接のみ」という理由を書く。

#### 利用者向けメッセージの文言規約

`M1` / `M2` は開発上のマイルストーン区分であり、**利用者向けのエラー・警告文には出さない**。
「現行版では未対応」「this version does not support ...」のように、利用者が読んで
意味の取れる表現にする。既存の C++ 側にも漏れていた2箇所
(`load_toml.cpp` の `throw_fermion_guard`、`measure.cpp` の correlation_length 警告)は
本設計の実装に先立って修正済み。設計書・実装解説など開発者向け文書では M1 / M2 を
そのまま使ってよい。

## 7. C++ 側の修正(必須)

`load_toml.cpp` の `validate_fermion_constraints` に、
**fermion モードで `ops_indices` が非空の2サイト観測量**を拒否するガードを追加する。

現状 `validate_fermion_constraints` の2サイト観測量ループは
`if (op.ops_indices.empty())` のときだけパリティ検査をしており、`ops` 形式は素通しで
測定時に `twosite_obs.cpp` のボゾン経路(`core::Contract(...)`)へ落ちる。

これはツール層の拒否(§6.2)だけでは塞げない。`input.toml` を直書きする expert
利用者が同じ経路に入れるからである。したがってツール側と独立の**必須修正**として、

```cpp
if (!op.ops_indices.empty()) {
  throw_fermion_guard("two-site observables in ops form");
}
```

を同ループに追加する。fermion モードの環境テンソルは reduced tensor から作られており
脚の形が違うため実行時エラーになる見込みだが、意図の伝わらない失敗であり、
「沈黙して間違った数値を出さない」原則からも入力時に理由付きで止めるべき。

## 8. テスト計画

| 層 | 内容 |
|---|---|
| 1. Fock 構成器(d = 2) | `h_bond` のエルミート性、パリティ保存(奇⇄偶をつなぐ要素がゼロ)、`t = V = μ = 0` で零行列、脚順 `[in1,in2,out1,out2]` の固定(非エルミートな `c†₁c₂` 単体で転置の取り違えを検出) |
| 1b. Fock 構成器(d = 4) | 既知要素 ⟨↑↓, 0\|c†₁↑c₂↑\|↓, ↑⟩ = −1、サイト内順序 `\|up dn> = c†_up c†_dn \|0>`、`parity = [0,1,1,0]` と行列のブロック構造が整合すること |
| 2. 同値性(強い回帰) | spinless の生成物が手書きサンプル `sample/07_spinless_fermion/input.toml` のゲート・ハミルトニアン行列と数値一致 |
| 2b. τ 展開 | `tenes_std` の出力ゲートが τ → 0 で `(expm(−τh) − I)/(−τ) → h` に戻ること(d = 2 と d = 4 の両方)。`tenes_simple` の `h` と `tenes_std` のゲートを結ぶ唯一の検査 |
| 3. パイプライン | `simple.toml` → `tenes_simple` → `tenes_std` → 生成された `input.toml` に `fermion = true` と `parity` が入っていること、`ops =` が現れないこと |
| 3b. ガード | §6.3 の各拒否が実際に発火すること: 非正方格子、2近傍非ゼロ、`fermion = false` 明示、`[correlation_length]`、`tenes_std` の長距離ボンド、`ops` 形式の2サイト観測量 |
| 4. E2E(ctest 1件) | 自由フェルミオン(t = 1, μ = 0, D = 2, χ = 8)を3段パイプラインで解き、厳密値と比較。既存 `FreeFermion` と同じ閾値・実行時間(1分以内) |
| 5. 回帰 | 既存 ctest 28件が不変(ボゾン経路・golden に触れない) |

層 1b を独立に立てるのは、**d = 2 では規約の違いが見えないため**である。実装解説の
三経路の表のとおり、reduced-pair blob の「素ロード」と「入力+出力 swap」は
数保存演算子の場合 d = 2 で縮退し、d = 4 の R5 オラクルで初めて差が出た。
spinless サンプル一致(層2)は Hubbard の rank-4 出力の規約を**捕まえられない**。

Hubbard の 2D 定量検証は D = 4 の2サイト測定コスト(D¹² スケーリング)のため本設計の
スコープ外とし、性能改善後に行う(前設計書の性能課題節)。1D 極限の検証は
`work/electron-validation/` の手順を流用できる。

## 9. ドキュメント

- `docs/sphinx/{ja,en}/file_specification/simple_format.rst`: 模型一覧に2種を追加。
  `"hubbard"` はフェルミオン Hubbard、Bose-Hubbard は `"boson"` である旨を明記。
  パラメータ表、`initial` の取りうる値、観測量一覧
- `sample/` に3段パイプラインのサンプルを追加(既存 `07_spinless_fermion` は
  expert モードの例として残す)
- `NEWS.md` に項目追加

## 10. 実装順序

1. Fock 構成器 + 単体テスト(層1)
2. `Model` の2つの拡張点 — `twosite_ops_explicit`(§4.3)と
   `initial_state_vectors`(§5.1)。既存模型の出力が1バイトも変わらないことを
   既存 ctest で確認してから次へ進む
3. `SpinlessFermionModel` + スキーマ出力 + 同値性テスト(層2)
4. §6.3 のガード群 + ガードテスト(層3b)
5. `tenes_std` の `parity` パース/検証/出力 + 長距離ボンド拒否 + パイプラインテスト(層3、3b)
6. τ 展開テスト(層2b)
7. E2E ctest(層4)
8. `HubbardModel` + 層1b
9. C++ ガード追加(§7)
10. ドキュメント・サンプル・NEWS

ステップ2を独立させるのは、フェルミオン模型と既存模型の**回帰リスクを分離する**ため。
ここで既存 ctest 28件が緑なら、以降の失敗はすべてフェルミオン側に切り分けられる。
