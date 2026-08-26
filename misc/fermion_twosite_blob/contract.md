# 振る舞い契約書: フェルミオン2サイト blob の省メモリ化

（`test/fermion/impurity_blob.cpp` の仕様。追補1 で当初案の T1/T2 を取り下げ、
追補2 で複素ケース、追補3 でガード検査を追加した。最終的に実装されたのは
direct fuse 方式 — `design.md` 改訂2 を参照）

- 日付: 2026-08-25
- 読者: テスト作成者（この文書**のみ**を仕様として受け取る。設計書や
  計画者のテストコードは渡されない）
- テストの置き場所（当時の指示）: `test/fermion/impurity_blob.cpp` を新規作成し、
  `test/test_fermion_layer.cpp` 末尾の `#include "fermion/..."` 群に
  1 行追加して組み込む（既存の r2_convention.cpp などと同じ方式。
  doctest、実行バイナリは test_fermion_layer）。
- 疑わしい点はソース（`src/fermion/*.hpp`）と照合してよい。契約と
  ソースが矛盾したら、**テストを曲げず統括に差し戻す**こと。


> **改名の注記**: この文書中の `build_reduced_pair` / `fuse_doubled_cluster` は、
> その後 `..._naive`（rank-16 を作るリファレンス）に、`..._factorized` は
> `..._direct`（外積を経由せず物理脚を直接縮約する）に改名された。本文は当時の記録としてそのまま残してある。

## 0. 背景（1 段落だけ）

フェルミオン模式の 2 サイト観測量 blob を作る既存関数
`tenes::fermion::build_reduced_pair`（`src/fermion/reduced.hpp`）は正しいが
巨大な中間テンソルを作る。これと**同一のテンソル**を省メモリで返す新実装
`build_reduced_pair_factorized` と、その部品 2 つが今回の対象。既存関数は
オラクル（正解生成器）としてテストから自由に呼んでよい。

## 1. テスト対象 API（スタブは統括が用意する。RED はコンパイルエラー
ではなく実行時の `std::logic_error` で出る）

すべて `namespace tenes::fermion`、`src/fermion/fops.hpp` /
`src/fermion/reduced.hpp` 宣言。テンソル型は既存テストと同じ
`tenes::real_tensor`（`_NO_MPI` ビルド）でよい。

```cpp
// (1) 演算子分解
template <class tensor>
std::pair<ftensor<tensor>, ftensor<tensor>> split_pair_op(
    const ftensor<tensor>& op12);

// (2) 不純物サイト二重化
template <class tensor>
tensor build_reduced_impurity(const ftensor<tensor>& Tn,
                              const ftensor<tensor>& op_half);

// (3) 新 blob（既存 build_reduced_pair と同一シグネチャ・同一戻り値）
template <class tensor>
tensor build_reduced_pair_factorized(const ftensor<tensor>& TnA,
                                     const ftensor<tensor>& TnB,
                                     const ftensor<tensor>& op12,
                                     reduced_pair_direction direction);
```

## 2. 規約（テストデータを作るのに必要な事実）

- 物理パリティ: d=2 スピンレスは `[偶, 奇]`、d=4 電子は `[偶, 奇, 奇, 偶]`。
- サイトテンソル `Tn` は rank 5、脚順 `(l, t, r, b, s)`（0..4）。仮想脚の
  パリティは任意に与えてよい（例: D=2 で `[偶,奇]`、全偶、全奇）。
  ランダム ftensor の作り方は既存テスト（`test/fermion/r2_convention.cpp` の
  `make_r2_tensor` まわり）を参考にしてよいが、コードの共有はしない
  （このファイルは同一 TU に include されるため、ヘルパは自前の
  無名 namespace に書く）。**パリティ保存則は課さなくてよい**（既存の
  等価性テストもランダム密テンソルで pin している。ただし §5 の
  物理演算子ケースは別）。
- 2 サイト演算子 `op12` は rank 4、脚順 `(in_A, in_B, out_A, out_B)`、
  パリティ `(p1, p2, p1, p2)`。本番と同じロード規約は
  `wrap_reduced_pair_op(op_plain, p1, p2)`（`fops.hpp`）。
  総パリティ**偶**の演算子のみが契約対象（総パリティ奇の挙動は未定義。
  テスト不要）。
- source（演算子の第 1 脚が乗るサイト）が窓の第 2 サイトになる本番ケースは、
  wrap 後に graded transpose `transpose(o, Axes(1,0,3,2))` を掛けてから
  blob 関数に渡す（`src/iTPS/twosite_obs.cpp` の fermion CTM 分岐と同順）。
- `reduced_pair_direction` は `horizontal`（A が左、接続は r_A×l_B）と
  `vertical`（A が上、接続は b_A×t_B）。
- blob（rank 6、fused D² 脚）の脚順は既存 `build_reduced_pair` の出力と
  同じ。**等価性テストでは脚順を自分で解釈する必要はない**（新旧を
  同一引数で呼び elementwise 比較するだけでよい）。
- graded QR `tenes::fermion::qr(a, rows, cols, q, r)` の内部脚パリティは
  `[偶×n_e, 奇×n_o]` の序列（`fops.hpp` 参照）。

## 3. 許容誤差

- elementwise: `|new − old| ≤ 1e-12 × max(1, max|old|)`。
- 全要素を走査すること（ノルムや和だけの比較は不可 — 符号の相殺を
  見逃す）。`local_size()` ループ + `global_index` / `get_value` で
  全要素比較するのが確実。

## 4. T1: `split_pair_op` の契約

対象演算子（すべて wrap 済みを渡す。d=2 と d=4 の両方）:

1. identity（d×d 単位行列の 2 サイト版: `op[a,b,c,e] = δ_ac δ_be`）
2. hopping: d=2 では `c†_A c_B + c†_B c_A`
   （非零要素 `op[0,1,1,0] = op[1,0,0,1] = 1`）
3. 密度 n_A n_B（d=2: `op[1,1,1,1] = 1`）
4. ペアリング `c_B c_A + c†_A c†_B`（d=2: `op[1,1,0,0] = op[0,0,1,1] = 1`。
   反対称演算子の符号規約に注意 — 要素表の通りに作ること。
   自分で第二量子化から導出し直して符号を「修正」しない）
5. ランダム総パリティ偶（パリティが偶になる要素だけ非零、複数シード）

検査項目:

- **再結合**: `tensordot(opA, opB, Axes(2), Axes(2))`（graded 版、
  `tenes::fermion::tensordot`）→ 脚順 `(in_A, out_A, in_B, out_B)` →
  graded `transpose(Axes(0,2,1,3))` → 元の wrap 済み op12 と
  elementwise 一致（§3 の許容誤差）。
- **脚形状**: opA・opB とも rank 3、脚順 `(in, out, k)`。k 次元は
  両者equal かつ ≤ d²。
- **k パリティ序列**: opA と opB の k 脚パリティベクトルが同一で、
  `[偶...偶, 奇...奇]` の序列（偶が前）。
- **セクタ整合**: hopping では k の奇セクタに、n_A n_B では偶セクタに
  重みが乗る（反対セクタの opA 列は全要素ゼロ）。「重みが乗る」の
  判定は該当セクタの max|要素| > 0.1 でよい（QR は正規化するため
  厳密値は基底依存 — **厳密値をテストしない**こと）。

## 5. T2: `build_reduced_impurity` の契約

QR の基底自由度があるため、**半身単独の数値は仕様ではない**。
検査はペア使用と構造のみ:

- **形状**: 戻り値は rank 5 plain tensor、脚 `(L,T,R,B,k)`、
  L..B の次元は各仮想脚次元の 2 乗、k は op_half の k 次元。
- **identity 合成 = 既存 norm blob**: identity op12 を `split_pair_op` で
  割り、両半身をそれぞれ `build_reduced_impurity` に通し、
  `build_reduced_pair_factorized` と同じ join …は T2 では使えないので、
  代わりに**identity ケースの全経路検査**を T3 に置く（下記 T3-1）。
  T2 単体では次を検査する:
  - k trace 検査: op_half として「identity の分解半身」ではなく
    **自明半身**（`op_half[i,o,k] = δ_io δ_k0`、k 次元 1、k パリティ偶）を
    手で作って渡すと、結果の k=0 スライスが既存
    `build_reduced(Tn)`（`reduced.hpp`、rank 4）に elementwise 一致する。
    これは基底自由度に依らない厳密な等式である。
  - 物理脚 open 版との整合: 自明半身で k 次元 1 の代わりに
    `op_half[i,o,k] = δ_io δ_k,k0`（k 次元 2、k0=0、k パリティ全偶）でも
    k=0 スライスは同じで k=1 スライスは全零。
- 乱数 Tn は D=1 と D=2（混合パリティ）、d=2 と d=4。

## 6. T3: 等価性（本丸）

すべて `build_reduced_pair_factorized(...)` と `build_reduced_pair(...)` を
**同一引数**で呼んで elementwise 比較（§3）。

**T3-1 基本掃引**（全組合せ）:

- 方向: {horizontal, vertical}
- source 位置: {第 1 サイト（wrap のみ）, 第 2 サイト（wrap + graded
  transpose (1,0,3,2)）} — **これは op 種とは独立の直交軸**
- d: {2, 4}
- D: d=2 は {1, 2, 3}、d=4 は {1, 2}
- 仮想脚パリティ: {全偶, 混合, 全奇}（D=1 は自明なので全偶のみ）
- op 種: {identity, hopping, nn, ランダム偶}（d=2 は pairing も）

組合せ爆発は seed を固定して各軸から代表を選ぶ形で間引いてよいが、
**「方向 × source 位置 × (d,D) × 混合パリティ × hopping」の全組合せは
省略不可**（奇セクタ k と bond 交差符号の検出器のため）。

**T3-2 d=4, D=3 ケース**（環境変数 `TENES_RUN_IMPURITY_BLOB_SLOW` が
設定されたときだけ実行、未設定なら doctest で skip/早期 return。
既存の env-gated 診断と同じ流儀）:

- direction {h,v} × op {hopping, transpose 済みランダム偶} × 混合パリティ。
  旧経路が 1 ケースあたり約 1 GB / 数十秒かかるための隔離である。

**T3-3 dressing の pin**（実装が確定させる符号 dressing の発生位置を
恒久固定する。実装完了後に統括から dressing の最終仕様が契約追補として
渡されるので、このテストは**追補待ちのプレースホルダを作らず**、追補が
来てから追加してよい。T3-1/T3-2 が先）。

## 7. 回帰と変異耐性（テスト作成者への注意）

- 既存テスト・既存 golden data は一切変更しない。新ファイル追加と
  `test_fermion_layer.cpp` への include 1 行だけが許される変更。
- 各テストは「新実装が正しく実装されたら緑、スタブ（logic_error 送出）の
  間は赤」であること。スタブ段階で緑になるテストは検出力ゼロなので
  書かない（例: 例外が出たら成功、のような書き方は禁止）。
- 乱数はシード固定で再現可能に。doctest の `INFO()` に seed・軸の値を
  出し、落ちたケースが一意に特定できること。

---

## 追補1（2026-08-25、改訂2）: T1/T2 の取り下げ

設計変更により、`split_pair_op` と `build_reduced_impurity` は API から
撤去されることになった。契約の変更は以下のみ:

1. **削除**: テストケース
   `impurity blob T1: split_pair_op recombines and orders k sectors` と
   `impurity blob T2: build_reduced_impurity with a trivial half matches
   build_reduced` を丸ごと削除する。これらのケースだけが使っていた
   ヘルパ（ib_trivial_half、ib_check_impurity_k_slice 系など）も併せて
   削除してよい（未使用警告を残さない）。
2. **不変**: T3-1（2 ケース）と T3-2 は**一字も変えない**。
   `build_reduced_pair_factorized` のシグネチャ・意味論・許容誤差は
   §1/§3/§6 のまま有効である（実装方式が変わるだけで、既存
   `build_reduced_pair` との elementwise 等価という契約は不変）。
3. **取り下げ**: §6 の T3-3（dressing pin、追補待ちだったもの）は
   不要になった。プレースホルダは元々無いので作業は発生しない。
4. 期待される結果: フルスイートは 136 cases / 134 passed / 2 failed
   （T3-1 の 2 ケースのみ logic_error で赤、T3-2 は env 未設定で緑）。

## 追補2（2026-08-25、改訂2 レビュー反映）: 複素テンソルの小型ケース追加

T3-1 に **`tenes::complex_tensor` 版の等価性ケース**を追加する（常用、
env ゲート不要の小型のみ）:

- 型: `tenes::complex_tensor`（`tenes::fermion::ftensor<tenes::complex_tensor>`）。
  既存ヘルパの real 固定部分は complex 対応の薄い複製で構わない
  （サイトテンソルの要素は実部・虚部とも決定論的擬似乱数で埋める。
  演算子は (a) hopping — 行列要素は実のまま型だけ complex、
  (b) ランダム総パリティ偶 — 複素要素、の 2 種）。
- ケース: d=2, D=2, 混合パリティで
  {hopping, ランダム複素偶} × {horizontal, vertical} × source 第1サイトの
  4 ケース + source 第2サイト × ランダム複素偶 × 1 方向の 1 ケース、計 5。
- 比較は同一引数で新旧を呼び、複素差の絶対値
  `|new − old| ≤ 1e-12 × max(1, max|old|)` を全要素走査。
- 期待 RED: 追補1 適用後のスイートに complex ケースが加わり、スタブ段階では
  これも logic_error で赤（doctest ケースは T3-1 の複素版として独立に
  1 本立てるとよい。136 → 137 cases / 134 passed / 3 failed になる）。

## 追補3（2026-08-25、最終レビュー反映）: parity 不一致ガードの検査

`src/fermion/reduced.hpp` の `detail::fuse_doubled_cluster_factorized` には
trace 対脚のパリティが bra/ket で食い違う入力を拒否するガードがある
（`std::runtime_error("fuse_doubled_cluster_factorized: traced-leg parity
mismatch")`）。本番入力からは構造的に到達しないが、ガード自体の検査を
1 件追加する:

- `tenes::fermion::detail::fuse_doubled_cluster_factorized` を直接呼ぶ。
  bra 側 rank-8 ftensor と ket 側 rank-8 ftensor を小さく作り
  （D=1, d=2 で十分。値はゼロでもよい）、ket 側の trace 脚
  （軸 3 または 7）のパリティベクトルだけ bra 側と食い違わせる。
  leg_ids は horizontal の `{0,1,3,1,2,3}` でよい。
- `CHECK_THROWS_WITH_AS`（または同等）で、上記メッセージの
  `std::runtime_error` が投げられることを検査する。
- 正常パリティなら投げないことも同じケース内で確認する（同じ形状で
  パリティを一致させて呼び、例外が出ないこと。戻り値の中身は検査不要）。
- 追加は独立の小さな TEST_CASE 1 本。既存ケースは一字も変えない。
- 期待: 追加後スイートは 138 cases / 138 passed（実装済みガードなので
  最初から緑。検出力はガード除去変異で統括が確認する）。
