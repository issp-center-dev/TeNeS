# tenes_simple / tenes_std のフェルミオン対応 設計書

日付: 2026-08-20
ブランチ: `fermion`
前提: `docs/superpowers/specs/2026-08-15-fermion-m1-design.md`(改訂1〜4)で確定した
C++ 側の規約。特に**改訂4の「演算子規約の最終表」**と、
`docs/superpowers/notes/2026-08-19-fermion-implementation-guide.md`(実装解説)。

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
- `tenes_std` の `parity` 通し(現在は落ちている)と fermion 時の入力検証
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
| `tenes_simple` | **順序付き2サイト Fock 基底での行列要素**(サイト内・ボンド内の反交換符号込み)を作る |
| `tenes_std` | 行列に対する `expm(-τh)`。**符号の知識は持たない**(通し役) |
| C++ (`tenes`) | サイト**間**・2次元幾何由来の交換符号を graded 機構で生成 |

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
    """

def local_index_to_occupation(i: int, nspin: int) -> List[int]:
    """局所基底 i を (n_up[, n_dn]) に。spinless: i = n。
    spinful: i = n_up + 2 n_dn、すなわち |0>, |up>, |dn>, |up dn>。
    サイト内順序は |up dn> = c†_up c†_dn |0> に固定(file_specification と同一規約)。
    """

def bond_matrix(fock_op: np.ndarray, nspin: int) -> np.ndarray:
    """2サイト Fock 空間の行列を rank-4 の op[in1, in2, out1, out2] に並べ替える。
    op[i1,i2,o1,o2] = <o1 o2| O |i1 i2>(順序付き基底)。
    """

def onesite_matrix(fock_op_1site: np.ndarray, nspin: int) -> np.ndarray:
    """1サイト演算子(必ずパリティ偶)を rank-2 に。"""
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
- 2サイト観測量: `bond_hamiltonian`、`hopping` = ⟨c†₁c₂ + h.c.⟩、`nn` = ⟨n₁n₂⟩

### 4.2 `HubbardModel`

- `physical_dim = 4`、局所基底 `|0>, |up>, |dn>, |up dn>`、`parity = [0, 1, 1, 0]`
- パラメータ: `t`、`u`(オンサイト斥力)、`v`(最近接斥力)、`mu`、`h`(Zeeman、任意)
- ボンドハミルトニアン:

  H_bond = −t Σ_σ (c†₁σc₂σ + h.c.) + V n₁n₂
           + (1/z)[ U(n↑n↓)₁ + U(n↑n↓)₂ − μ(n₁+n₂) − h(S^z₁+S^z₂) ]

- 1サイト観測量: `n`, `n_up`, `n_dn`, `Sz`, `doublon` = n↑n↓, `holon` = (1−n↑)(1−n↓)
- 2サイト観測量: `bond_hamiltonian`、`hopping`、`nn`、`SzSz`

### 4.3 共通の制約

- **1サイト観測量はパリティ偶のもののみ出力する**。`c†`, `c` 単体は奇で、C++ 側の
  ロードガードに弾かれるため出力しない(相関関数も M1 では非対応)
- **2サイト観測量は必ず明示 `elements` で出力する**。`ops = [i, j]` 形式(1サイト演算子の
  直積)は fermion モードで未対応(§6)

## 5. 初期状態

既存の `[lattice] initial` キーを使い、フェルミオン模型では値を差し替える。
`"ferro"` / `"antiferro"` は**受け付けない** — 実際に作れるのは電荷秩序であり、
特に `"antiferro"` は SDW を想起させて誤解を招くため。

| `initial` | spinless fermion | hubbard |
|---|---|---|
| `"random"`(既定) | 乱数偶初期化(`initial_state = [0.0, 0.0]`) | 乱数偶初期化 |
| `"vacuum"` | 全 `\|0>` | 全 `\|0>` |
| `"full"` | **エラー** | 全 `\|up dn>`(n = 2) |
| `"cdw"` | **エラー** | `\|0>` / `\|up dn>` の市松(n = 1) |
| `"ferro"`, `"antiferro"` | エラー | エラー |

`"vacuum"` は `noise` が全要素(偶マスク付き)に乗るため発展の不動点にはならず、
低密度側から出発する実用的な初期値になる。

spinless の `"full"` / `"cdw"` と Hubbard の Néel 型が作れないのは、`|1>` や `|up>` が
パリティ奇で、TeNeS のプロダクト初期化が仮想脚の添字 0(偶)に状態ベクトルを置くため
サイトテンソルの全脚パリティ合計が奇になるから。エラーメッセージにこの理由と回避策
(`"random"`)を書く。

## 6. スキーマと出力

### 6.1 `tenes_simple` → `std.toml`

- `[parameter.general] fermion = true` を模型種別から**自動注入**する(利用者は書かない)。
  利用者が明示的に `fermion = false` と書いていたら、矛盾として明示エラー
- `[[tensor.unitcell]]` に `parity = [...]` を追加出力
- 2サイト観測量は `ops = [i,j]` を使わず `elements` を出力(`is_complex` 分岐に関わらず)

### 6.2 `tenes_std` → `input.toml`

- `[[tensor.unitcell]]` の書き出しに `parity` を追加(**現状は落ちている** — 唯一の必須修正)
- `[parameter]` は既にテーブルごと素通しなので `fermion` は自動的に伝わる
- fermion 時の検証を追加:
  - `[[observable.twosite]]` が `ops` 形式なら明示エラー(理由と `elements` 形式への
    案内を含める)
  - `parity` の長さが `physical_dim` と一致すること
- **`expm(-τh)` の計算ロジックは変更しない**

## 7. C++ 側の修正(1点)

`load_toml.cpp` の `validate_fermion_constraints` に、
**fermion モードで `ops_indices` が非空の2サイト観測量**を拒否するガードを追加する。

現状これは未ガードで、測定時に `twosite_obs.cpp` のボゾン経路
(`core::Contract(...)`)へ落ちる。fermion モードの環境テンソルは reduced tensor から
作られており脚の形が違うため実行時エラーになる見込みだが、意図の伝わらない失敗であり、
「沈黙して間違った数値を出さない」原則からも入力時に理由付きで止めるべき。

## 8. テスト計画

| 層 | 内容 |
|---|---|
| 1. Fock 構成器 | `h_bond` のエルミート性、パリティ保存(奇⇄偶をつなぐ要素がゼロ)、既知要素 ⟨↑↓,0\|c†₁↑c₂↑\|↓,↑⟩ = −1、`t = U = V = μ = 0` で零行列 |
| 2. 同値性(強い回帰) | spinless の生成物が手書きサンプル `sample/07_spinless_fermion/input.toml` のゲート・ハミルトニアン行列と数値一致 |
| 3. パイプライン | `simple.toml` → `tenes_simple` → `tenes_std` → 生成された `input.toml` に `fermion = true` と `parity` が入っていること、`ops =` が現れないこと |
| 4. E2E(ctest 1件) | 自由フェルミオン(t = 1, μ = 0, D = 2, χ = 8)を3段パイプラインで解き、厳密値と比較。既存 `FreeFermion` と同じ閾値・実行時間(1分以内) |
| 5. 回帰 | 既存 ctest 28件が不変(ボゾン経路・golden に触れない) |

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
2. `SpinlessFermionModel` + スキーマ出力 + 同値性テスト(層2)
3. `tenes_std` の `parity` 通し + 検証 + パイプラインテスト(層3)
4. E2E ctest(層4)
5. `HubbardModel`(層1のテストを d = 4 に拡張)
6. C++ ガード追加
7. ドキュメント・サンプル・NEWS
