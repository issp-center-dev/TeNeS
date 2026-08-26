# フェルミオン2サイト blob の不純物分解 — 設計書

- 日付: 2026-08-25
- ブランチ: fermion（起点タグ: `wip/fermion_20260825`）
- 置き場所: `misc/fermion_twosite_blob/`（当時の作業ディレクトリは
  `work/fermion/impurity-blob/`。台帳・レビュー報告・タスク指示書はそちらに残置）
- 状態: **改訂2（アーキテクチャのピボット）**。実装完了・マージ済み（`37a66c82`）。
  改訂1 の「graded QR 分解 +
  不純物二重化 + join」は、プローブ実験により join の符号 dressing が
  開脚パリティの 2 次形式（サイト間 joint マスク項）+ wrap 規約補正を含む
  ことが判明し（`probe_join.cpp`）、同じ実験系からより単純で等価な解 —
  **direct fuse（外積の直接縮約化）** — が見つかったため置き換える。
  direct fuse の恒等式は `probe_direct.cpp` で 22 ケース機械精度検証済み。
- 経緯: 改訂1 = Codex レビュー指摘 5 件反映。T1 前提（graded QR 分割）は
  `probe_split.cpp` で検証済み（QR 分解自体は正しいが、本改訂で不要になった。
  blob を作らないストリーミング化という将来課題の部品として記録に残す）。


> **改名の注記**: この文書中の `build_reduced_pair` / `fuse_doubled_cluster` は、
> その後 `..._naive`（rank-16 を作るリファレンス）に、`..._factorized` は
> `..._direct`（外積を経由せず物理脚を直接縮約する）に改名された。本文は当時の記録としてそのまま残してある。

## 1. 目的と背景

フェルミオン模式の 2 サイト観測量（CTM 経路、`meanfield_env = false`）は
`build_reduced_pair`（`src/fermion/reduced.hpp:246`）が **rank-16 の全外積**
（`fuse_doubled_cluster` 内 `tensordot(bra, ket_pair, Axes(), Axes())`、
`reduced.hpp:168`）を実体化する。要素数は (D⁶d²)² = D¹²d⁴。

実測（2026-08-25、Hubbard t=1 μ=0 U=4、2×2 セル、χ=40、`is_real=true`）:

| | 1コピー | peak RSS | 2サイト測定時間 (serial) |
|---|---|---|---|
| D=3, d=4 | 1.09 GB | 4.85 GB（≈3.8 コピー同時） | 521 s |
| D=4, d=4 | 34.4 GB | 実機 176 GB (np=4, PBS 会計) | 外挿 ≈4.6 h |

D=4 は 1 ノード 240 GB で cgroup OOM kill（実機で確認、2 ノード分散で回避済み）。
測定 1 回あたり blob を 32 回作る（4 演算子 × 8 ボンド）。
本設計はこの外積を排し、**中間テンソルの最大サイズを rank-6 blob 本体
（(D²)⁶ = D=4 で 134 MB）まで落とす**。norm 側
（`build_reduced_identity_pair`）は既にサイト分解型で安価なので触らない。

## 2. 現状の構造（読解結果・規約表）

### 2.1 呼び出し経路

`src/iTPS/twosite_obs.cpp` の fermion CTM 分岐のみが対象:

```
o = wrap_reduced_pair_op(op.op, p_src, p_tgt)        // in+out 両 swap
if (source が窓の2番目) o = transpose(o, (1,0,3,2))   // graded transpose
blob = build_reduced_pair(fA, fB, o, direction)       // ← 本設計の置換対象
value = contract_reduced_pair_{horizontal,vertical}_density_CTM(..., blob)
```

MF 経路（`contract_pair_MF`）、ボソン経路、density(TPO) 経路は対象外。

### 2.2 脚規約（現物から確定）

- wrap 済み `Tn`: rank 5、`(0=l, 1=t, 2=r, 3=b, 4=s)`。
- `build_pair_state`（`reduced.hpp:212`）:
  - horizontal: `tensordot(A,B,Axes(2),Axes(0))` → `(l_A,t_A,b_A,s_A, t_B,r_B,b_B,s_B)`
  - vertical:   `tensordot(A,B,Axes(3),Axes(1))` → `(l_A,t_A,r_A,s_A, l_B,r_B,b_B,s_B)`
- `apply_pair_op`（`reduced.hpp:231`）: op12 の脚は `(in_A, in_B, out_A, out_B)`。
  pair の物理脚 (3,7) に縮約後、元の脚順に戻す。
- `doubled_pipeline`（`reduced.hpp:117`）: 単サイト二重化。出力
  `([l lb],[t tb],[r rb],[b bb], s_ket, s_bra)`、interleave は (ket, bra) 順、
  凍結 joint-swap マスク `kDoubledJointMask`（YASTN fuse_layers 相当、oracle で pin）。
- `build_reduced_op` = 物理脚 open の rank-6、`build_reduced` = その trace の rank-4。
- `fuse_doubled_cluster`（`reduced.hpp:149`）: pair の 6 外部脚
  （horizontal は leg_ids `{0,1,3,1,2,3}` = l_A,t_A,b_A,t_B,r_B,b_B）で
  joint-swap を計算し、**bra 側 swap を外積前に適用**（rank-8 のうちに済ませる
  最適化、2026-08-23）、外積 → interleave 転置 → fuse → 物理 trace。
- blob（rank 6）の脚順:
  - horizontal: `(L_A, T_A, B_A, T_B, R_B, B_B)`（fused D² each）
  - vertical:   `(L_A, T_A, R_A, L_B, R_B, B_B)`
- `apply_fused_leg_gauge`（op blob のみ、identity pair には無い）:
  - horizontal: leg 2 (=B_A) に `TnA.parity[3]`, ket_odd_bra_even=true /
    leg 5 (=B_B) に `TnB.parity[3]`, false
  - vertical: leg 0 (=L_A) に `TnA.parity[0]`, true / leg 3 (=L_B) に
    `TnB.parity[0]`, false
  - コード注記のとおり「解析導出ではなく比較で測定」された規約。
- `build_reduced_identity_pair`（`reduced_measure.hpp:62`）: rank-4 reduced 同士の
  1 縮約。gauge なし。**blob と同じ脚順**なので CTM 縮約関数を共有している。

### 2.3 既存の検証資産（今回の契約の土台）

- `test/fermion/sign_sweep.cpp:1741-` に現行 `build_reduced_pair` の**逐語コピー**
  （ss_reference3）があり、本体との一致を pin している。
- `test/fermion/r2_convention.cpp`: R3（d=2 oracle、blob とサイト不純物の整合）、
  R5（d=4、Fock 検証済み direct path との一致）。**新経路はこの網を素通しで
  通らなければならない**（blob の値自体が同一になるため、テスト変更は不要のはず）。
- 入力検証（`load_toml.cpp:604` `validate_fermion_constraints`）が
  **パリティ奇の 2 サイト演算子・ops 形式・距離 2 以上を既に拒否**。
  よって本番の op12 は常に総パリティ偶で、graded QR のブロック対角前提は
  solver 経路では入力検証が保証する。

## 3. 新設計（改訂2）: direct fuse — 外積 + trace を直接縮約に置き換える

### 3.1 恒等式

現行 `fuse_doubled_cluster`（`reduced.hpp:149`）は

1. bra へ `forms.bra` を適用（bra 内部の swap 項）
2. **rank-16 外積**（graded tensordot、縮約軸なし = マスクなしの純外積）
3. `transpose_with_swap_form(doubled, forms.cross, interleaved)` —
   `forms.cross` + 転置の Koszul 項を**対角符号マスク**として掛けてから
   plain 転置
4. fuse（reshape）して物理脚 4 本を **plain trace**（δ 縮約）

の順で blob を作る。ここで手順 2〜4 は数学的に次と同値である:

- 全符号項（`forms.cross` + `transpose_sign_form(interleaved)`）は
  脚ペア (i,j) 上の対角マスク (-1)^{p_i p_j} の集合。これを所属で振り分ける:
  - **bra 内部ペア** → bra テンソルへのマスク
  - **ket 内部ペア** → ket テンソルへのマスク
  - **trace 対脚を含むペア** → δ trace は添字を同一視するので、
    **もう一方の脚の所属に合わせて** bra→ket または ket→bra のどちらにも
    twin 書き換えしてよい（(bra s × ket開脚) ≡ (ket s' × ket開脚) は ket
    内部へ、(ket s' × bra開脚) ≡ (bra s × bra開脚) は bra 内部へ）。
    twin 同士のペア (bra s × ket s') は p(s)² = p(s) の**線形パリティ
    マスク**に落ちる（probe_direct.cpp の分類ロジックがこの正）
  - **open×open のサイト間ペア** → 縮約後の rank-12 結果へのマスク
- マスク振り分け後、物理脚 2 対を **plain tensordot で直接縮約**
  （`tensordot(bra', ket', Axes(3,7), Axes(3,7))`）すれば rank-12
  （= D¹² = blob と同要素数）しか実体化しない。
- 残りは plain の interleave 転置 + fuse reshape（従来と同じ）。

**この恒等式は `probe_direct.cpp` で検証済み**（2026-08-25）: 方向 {h,v} ×
op {nn, hopping, pairing, random 偶, source 第2サイト} × d {2,4} × D {1,2,3}
の 22 ケース全てで既存 `build_reduced_pair` と elementwise 一致
（maxdiff ≤ 4.4e-16）。**新しい符号規約は一切導入していない** — 既存の
凍結マスクの項集合を F₂ で振り分けただけである。

### 3.2 メモリ・時間

D=4, d=4, blob 1 個あたりの最大単一テンソル:

| 段階 | 要素数 | メモリ |
|---|---|---|
| 現行: rank-16 外積 ×(4〜5 コピー) | 4.29e9 | 140〜180 GB |
| 新: rank-12 直接縮約結果（+転置コピー） | 1.68e7 | 134 MB ×2〜3 |
| 新: blob（従来どおり） | 1.68e7 | 134 MB |

rank-12 は blob と同要素数なので、**blob を出力とする限りこれが下限**。
縮約は GEMM（D⁶d² × D⁶d²、d² 縮約 = D=4 で 2.7e8 flops）で、支配項は
マスク掛けと転置の O(D¹²) スイープ数回。D=3 実測 521 s の 2 サイト測定は
数秒台、D=4 でも measurement 全体で CTM 環境縮約が支配項になる見込み
（T4 で実測）。ピーク RSS は測定パス全体で数百 MB〜1 GB 程度を見込む。

### 3.3 API（改訂2）

- 現行 `build_reduced_pair` は**名前も実装も無変更**（オラクル。
  ss_reference3 のビット厳密 pin も無傷）。
- 新関数 `build_reduced_pair_factorized`（同一シグネチャ、スタブ導入済み）を
  direct fuse で実装する。名前の由来: 外積を経由せず物理脚を直接縮約する、の意
  （当初は `_lean`（メモリが lean）としたが、Lean 定理証明器と紛らわしく
  分野の用語でもないため `_direct` に改めた）。**QR 分解・build_reduced_impurity・
  split_pair_op は改訂2 で不要になった**（スタブと T1 実装は撤去する）。
- 実装完了時に `twosite_obs.cpp` の呼び出し 2 箇所だけを factorized に
  切り替える。
- 実装の参照コード: `probe_direct.cpp` の
  `direct_fuse` / `direct_pair`（設計検証済みの試作。製品化では elementwise
  ループを ftensor の `apply_swap_form` / `multiply_vector` 系に置き換え、
  項の振り分けは同一のロジックを用いる）。

## 4. 符号リスク（改訂2 で消滅）

改訂1 の主リスク「join の符号 dressing」は、direct fuse では**存在しない**。
新規約はゼロで、既存 `joint_swap_forms` の出力と転置 Koszul 項の
機械的な振り分けのみ。振り分けの正しさは
(1) probe_direct の 22 ケース機械精度一致（設計段階、済）、
(2) T3 等価性テストの全掃引（実装後、作成者のテスト）、
(3) 既存 C3（ビット厳密）/ R3 / R5 / E2E の全回帰
の三重で検証される。複素テンソル経路（conj の取り扱い・型依存分岐）は
probe が実数のみのため、契約追補2 の小型 complex ケースで pin する
（r2 レビュー important-1）。

参考（記録）: 改訂1 ルートを probe_join.cpp で調べた結果 —
naive join は nn D=1 で全体 −1（wrap の in/out-swap マスクが単脚縮約では
相殺しない）、hopping では対角 gauge で埋まらない構造的不一致
（サイト間 joint 項は開脚パリティの 2 次形式）。この知見が direct fuse への
ピボットの直接の動機である。
## 5. 検証契約の要点（契約書は別ファイルで作成）

契約書は `contract.md`（同ディレクトリ。v1 発行済み・テスト作成済み）。改訂2 に伴う
修正は**契約追補**として作成者へ渡す: **T1（split_pair_op）と
T2（build_reduced_impurity）のテストケースを削除**（API 撤去のため）。
**T3-1 / T3-2 は一字も変えない**（`build_reduced_pair_factorized` 対
既存オラクルの等価性テストは実装方式に依存しないため、direct fuse でも
そのまま受け入れ基準になる）。T3-3（dressing pin）は dressing 自体が
消滅したため**不要になった**（追補で正式に取り下げる）。v1 の柱（記録）:

- **T1 分解**: `split_pair_op` の再結合一致（graded tensordot で
  opA·opB == 元 op、elementwise、d=2/d=4、ランダム偶演算子＋物理演算子
  （hopping/nn/SzSz/pairing）、k パリティ序列 `[F..,T..]`、奇入力は
  Debug で検査例外）。
- **T2 不純物サイト**: `split_pair_op` の半身は QR の基底選択に依存する
  gauge 自由度を持つため、**半身単独の数値は契約対象にしない**
  （レビュー指摘 minor-2）。T2 の契約は (i) shape・k パリティ序列・
  ゼロセクタ構造の検査、(ii) opA/opB を**必ずペアで**使った検査
  （identity の分解ペアを doubled_impurity 2 つに通して join した結果が
  `build_reduced_identity_pair` に一致する等）に限定し、
  数値等価性の本体は T3 に寄せる。
- **T3 等価性（本丸）**: `build_reduced_pair_factorized`（新）==
  `build_reduced_pair`（既存、無変更）elementwise（許容 1e-12 相対）。
  掃引は**直交させる**（レビュー指摘 important-1）:
  方向 {h,v} × **source 位置 {窓の第1サイト, 第2サイト}** × d {2,4} ×
  D × パリティパターン（全偶 / 混合 / bond 全奇）× op 種
  {identity, hopping, nn, pairing, ランダム偶}。
  source 第2サイトのケースは本番と同順で
  `wrap_reduced_pair_op(op, p_src, p_tgt)` → graded transpose (1,0,3,2)
  を適用してから両経路に渡す（transpose は op 種の一つではなく独立軸）。
  D の掃引: d=2 は D {1,2,3}、d=4 は D {1,2} を全組合せ、
  **加えて d=4, D=3 を最小セット**（direction {h,v} × op {hopping,
  transpose 済みランダム偶} × 混合パリティ）で入れる
  （レビュー指摘 important-2。d=4 の奇セクタ多重度 × 偶奇非対称 bond の
  交差は D=2 以下では踏めない。旧経路 1 コピー ≈1 GB なので、この
  ケースだけ slow 指定を許す）。乱数 Tn は既存テストの make_r2_tensor 流儀。
  （改訂2 注: この段落末尾の dressing pin テストは v1 の記録であり、
  追補1 で取り下げ済み。追補1・追補2 が正。）
- **回帰**: 既存 ctest 全件（fermion R2-R5、E2E、ボソン系含む）が無変更で緑。
  golden data の再生成は**不要のはず**（blob 同一のため）。必要になったら
  それは等価性の破れ = バグ。
- **変異耐性**: レビュアー向け推奨 — join の gauge を 1 箇所落とすと
  T3 のどのケースが赤くなるか確認（混合パリティ×奇 op ケースが検出器）。

## 6. タスク分割（改訂2）

| # | 内容 | 主ファイル | 受け入れ |
|---|---|---|---|
| U1 | 撤去: split_pair_op（T1 実装済み分）・build_reduced_impurity スタブ。作成者による T1/T2 テスト削除後に統括が実施 | fops.hpp, reduced.hpp | ビルド緑、残テスト構成が想定どおり |
| T3' | `build_reduced_pair_factorized` を direct fuse で実装、呼び出し切替（twosite_obs.cpp 2箇所） | reduced.hpp, twosite_obs.cpp | 契約 T3-1/T3-2 緑 + 全 ctest 緑 |
| T4 | 性能実測（D=3/D=4、ピーク RSS・時間）、記録 | work/ | D=4 χ=40 が mac 1 プロセスで完走、実測記録 |

手順: 契約追補 → 作成者が T1/T2 削除（それ以外バイト不変）→ 統括が
RED 再確認（T3 の 2 ケースが logic_error で赤のまま）+ スナップショット v2
→ U1 → Codex 実装（テスト変更禁止）→ バイト検査 → 統括の独立検証
（全 ctest + 変異検査: direct fuse のマスク項を 1 つ落として T3 が赤くなる
ことを確認）→ タスクレビュー → 最終全体レビュー。

## 7. スコープ外

- MF 経路・density(TPO) 経路・ボソン経路（無変更）。
- norm 側 `build_reduced_identity_pair`（既に安価）。
- multisite / 相関関数のフェルミオン対応（別 M2 項目）。
- CTM 環境縮約関数群（blob 形式が同一なので無変更）。
- blob を作らず環境へストリームする更なる最適化（(D²)⁶ の 134 MB すら
  避ける案）。今回の修正で律速が環境縮約側に移ってから判断。

## 8. 検討済み代替案

- **B案: fused パリティ ftensor で doubled_pipeline 全体を graded 化**:
  新旧の一致点が増え、既存の凍結マスク規約と二重帳簿になるため不採用。
  join のみ graded にする本案が変更面積最小。
- **SVD による分解**: QR と等価だが s の再配分と truncation ノブが増えるだけ。
  QR 採用。
- **op12 を redB へ直接縮約（分解なし）**: (D²)⁶d⁴ の中間が再出現し
  本末転倒。不採用。
