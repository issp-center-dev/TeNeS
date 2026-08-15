# フェルミオン対応 M1 設計書 — swap gate 方式・正方格子・最近接模型

日付: 2026-08-15
ブランチ: `fermion`
前提調査: フェルミオン iPEPS アルゴリズム調査(2026-08-13、swap gate / Z₂グレード / Grassmann の3定式化と既存実装のサーベイ)

## 1. 背景と目的

現行 TeNeS はテンソル縮約においてすべての演算子・状態の入れ替え符号を +1 として扱うため、
フェルミオン系(Hubbard、t-J、スピンレスフェルミオン等)を扱えない。本設計は swap gate 方式
(Corboz–Orús–Bauer–Vidal, PRB 81, 165104 (2010))に基づき、密テンソル(mptensor 無改造)の
まま TeNeS リポジトリ内で完結するフェルミオン対応の第一段階(M1)を定める。

計算量への影響: 符号処理は要素ごとのマスク乗算であり、縮約のリーディングオーダーは
ボゾン系と同一(Barthel–Pineda–Eisert, PRA 80, 042333 (2009))。

## 2. スコープ

### M1 でやること

- 正方格子・**最近接**ボンドハミルトニアン・スピンレスフェルミオン
- simple update(虚時間発展)
- CTMRG 環境構築+ **1サイト・最近接2サイト観測量**(エネルギー、占有数、NN ホッピング相関)
- 入力: `input.toml` 拡張(fermion フラグ+物理脚パリティ宣言)+検証用入力生成スクリプト
- 実・複素両テンソル型対応(既存テンプレート構造を維持)
- 合格基準: **自由スピンレスフェルミオン(ホッピング t + 化学ポテンシャル μ)の
  厳密解(運動量積分)に対するエネルギー密度・占有数の D 収束**

### M1 でやらないこと(fermion モード時は load_toml 段階で明示エラー)

- full update / fast full update
- 距離 2 以上の相関関数・multisite 観測量・相関長(転送行列)
- 有限温度(purification / density 系)・実時間発展
- mean-field 環境での測定
- RSVD(fermion モードではフル SVD を強制)
- simple update のゲージ固定(`Simple_Gauge_Fix`)
- `tenes_simple` / `tenes_std` の本格対応(M2)
- U(1) ブロックスパース化(性能最適化として将来課題)

## 3. アーキテクチャ: 符号レイヤー

手作業で縮約図ごとに swap 挿入位置を導出する代わりに、**脚ごとのパリティメタデータを持つ
薄いラッパー `ftensor` と、符号を機械的に生成する f-プリミティブ関数群**を導入する
(locally-ordered 簿記。数学的には交差数え上げと同値だが、符号バグの混入面が小さい)。
密テンソルのまま・mptensor 無改造・ブロックストレージなし。

### 3.1 ftensor

```cpp
template <class tensor>
struct ftensor {
  tensor t;                                  // 密テンソル本体(mptensor 型)
  std::vector<std::vector<bool>> leg_parity; // 脚ごと・添字ごとのパリティ
};
```

パリティを「セクター次元 (D_e, D_o)」ではなく**添字ごとの任意 bool ベクトル**で持つ。
reshape(脚の融合)後の非連続パターンもそのまま扱える。

### 3.2 f-プリミティブ(mptensor 自由関数のミラー)

| 関数 | 意味論 |
|---|---|
| `ftranspose(A, axes)` | transpose + 転倒対 (i<j が逆転) ごとの (−1)^{p_i p_j} を要素マスクで適用 |
| `ftensordot(A, B, axesA, axesB)` | 「A の縮約脚を末尾へ、B の縮約脚を先頭へ寄せる置換」の符号+縮約。縮約脚のパリティ一致を実行時検証 |
| `fsvd` / `fqr` | パリティセクターごとの分解。特異値の大きい順にセクター横断で切断し、偶先頭に再ソートして新パリティを返す。縮退時の選択は決定論的(偶優先) |
| `fconj(A)` | 要素共役+脚順序反転(随伴)+元の軸番号への ftranspose。正味は「共役+全脚対マスク (−1)^{Σ_{i<j} p_i p_j}」(最終形は規約固定テストで確定) |
| `apply_swap(A, ax1, ax2)` | 奇⊗奇要素の符号反転(低レベル) |
| `apply_parity(A, ax)` | 指定脚の奇要素に −1(`multiply_vector` で実装可) |
| `check_parity_even(A)` | 禁止ブロックのノルム検査(デバッグ・テスト用) |

要素ごとアクセスは `tensors.cpp:114-125` の `local_size()`/`global_index()`/`set_value()`
ループが前例。

### 3.3 既存コードとの統合(境界ラップ方式)

- `iTPS<tensor>` に `FermionInfo` メンバを追加: 有効フラグ+物理脚パリティ(サイトごと、
  入力から)+仮想ボンドパリティ(動的、`fsvd` のたびに更新)。`Tn`・環境テンソルの格納型は
  既存のまま。
- フェルミオン化するルーチンの呼び出し境界で `ftensor{テンソル, パリティ}` を構築して渡し、
  戻りでメタデータを回収する。テンプレートカーネル自体は `ftensor` でインスタンス化し、
  自由関数のオーバーロード解決で符号付きに切り替える(**カーネルソースはボゾン/フェルミオン
  共通**)。`ftensor` は M1 対象コードが使う API サブセットのみミラーする
  (member: `transpose`/`multiply_vector`/`shape` など)。
- `iTPS<ftensor>` の全面インスタンス化は行わない(未移植機能への波及を避ける)。
- ボゾン経路はコード・golden ともに完全不変。

## 4. 規約

1. **基底順序**: 物理脚のパリティは input.toml で基底状態ごとに宣言(偶先頭を推奨規約)。
   **仮想脚は常に偶先頭ソート済み**を不変条件とする(`fsvd`/`fqr` が構成的に保証)。
2. **脚の正準順序**: 各 ftensor の格納軸順序がそのテンソルのフェルミオン順序。
   Tn の既存規約 (0=左, 1=上, 2=右, 3=下, 4=物理) を正準順序として採用。
   グローバルなサイト順序は不要。
3. **bra/ket**: `fconj` を §3.2 の式で定義。共役まわりの符号配置には規約の自由度があるため、
   **最終形は規約固定テスト(§7 層2)が決定する**: (i) ⟨ψ|ψ⟩ > 0、(ii) 有限ミニネットワーク
   縮約が JW 順序の厳密参照と一致、(iii) 水平・垂直の整合。
4. **演算子**: 2サイトゲートは既存 rank-4 `op[in1,in2,out1,out2]` 形式のまま、行列要素は
   **順序付き2サイト基底 |n₁n₂⟩ = (c†₁)^{n₁}(c†₂)^{n₂}|0⟩(site 1 = source_site)で定義**。
   NN スピンレスではゲート行列は符号なしで書け、幾何由来の符号はネットワーク側が負担。
   ロード時に宣言パリティに対する偶性を検証し、奇な演算子は M1 では拒否。
5. **λ**: 実正値対角のまま。仮想脚が偶先頭ソート済みなので自動的にパリティ整合し、
   `multiply_vector` は符号を生まない。

## 5. アルゴリズム変更

### 5.1 simple update(`core/simple_update.cpp` `Simple_update_bond`)

| ステップ | 変更 |
|---|---|
| λ 吸収(`multiply_vector`) | 変更なし(対角・偶) |
| 方向別 transpose(4分岐) | `ftranspose`(符号自動。方向別の手動符号解析は不要) |
| QR | `fqr` |
| θ = R1·R2·gate | `ftensordot` ×2 |
| SVD+切断 | `fsvd`。ボンドの (D_e, D_o) が動的に決まり、`FermionInfo` と λ を更新 |
| λ⁻¹ 戻し+逆 transpose | `ftranspose` |

### 5.2 CTMRG(`core/ctm.cpp` + `contract_itps_ctm/ctm.cpp` の1・2サイト測定カーネル)

- 環境テンソル C1..C4 / eT の全脚にパリティベクトルを付与(D 脚はボンドパリティの
  ket/bra コピー、χ 脚は射影子計算時の `fsvd` が返す動的パリティ)。初期化も f-プリミティブで。
- 射影子: `fsvd(ftensordot(...))`。吸収(`Left/Right/Top/Bottom_move`)と格子回転 transpose は
  `ftranspose`/`ftensordot` への機械的移植。
- 収束判定(C の特異値比較)は変更なし。
- bra/ket 符号は `fconj`+呼び出し列から自動発生。**文献規約(Corboz PRB 81, 165104 付録)で
  手導出した1サイト reduced tensor と f-プリミティブ経由の結果を突き合わせる橋渡しテストを
  1本置く**。

### 5.3 測定

- 1サイト(n など偶演算子): 既存カーネルの ftensor 移植。
- NN 2サイト(ホッピング等): `Contract_two_sites_{horizontal,vertical}_op12` の移植。
  演算子は §4-4 の規約で構築されたものを使用。
- `twosite_obs.cpp` の距離 2 以上の SVD 分解経路(`:226`)は fermion モードでガード
  (JW string 挿入が必要なため M2)。

## 6. 入力スキーマ・初期状態

- `[parameter.general] fermion = true`(モードスイッチ)
- `[[tensor.unitcell]] parity = [0, 1]`(長さ physical_dim の 0/1 リスト。fermion 時必須)
- 仮想脚の初期パリティ分割はデフォルト D_e = ⌈D/2⌉(入力指定不要。以後 `fsvd` が動的に決定)
- 検証: パリティ宣言の整合、演算子の偶性、§2 の未対応機能ガード
- **初期状態**: ランダム偶初期化(既存ランダム初期化ループに偶ブロックマスク)を推奨経路とする。
  偶パリティのプロダクト状態+偶ノイズも許可。**奇状態のプロダクト初期化はエラー**とし、
  メッセージで理由(奇テンソルは偶制約に反する)と回避策(ランダム偶初期化)を案内する。
- 出力ファイル形式(`output/*.dat`)は変更なし。

## 7. テスト計画(4層)

1. **単体(doctest)**: `apply_swap`/`ftranspose` を素朴な要素ごと参照実装と比較。
   `ftensordot` の縮約順序独立性・パリティ不整合検出。`fsvd`/`fqr` のブロック厳密性・
   セクター横断切断・偶先頭ソート・因子の偶性。
2. **規約固定**: 有限ミニネットワーク(2×2、小 D)の縮約値を JW 順序の numpy/ED 厳密参照と
   比較し、`fconj` の規約を確定する。文献の1サイト reduced tensor との橋渡しテストを含む。
3. **統合(ctest)**:
   - **1D 極限**: t_y=0 でフェルミオン計算とハードコアボソン計算(既存経路)が一致し、
     t_y≠0 で乖離すること(符号機構の狙い撃ち)。
   - **自由フェルミオン D 収束**: t=1、μ 数点、D=2..6 でエネルギー密度・占有数が厳密値に
     単調接近。閾値は実測して golden 方式で固定。
4. **回帰・ガード**: 既存 ctest 全件が golden 据え置きで不変。fermion+未対応機能の
   組み合わせがエラーになることの検査。

## 8. ファイル構成・見積り・実装順序

### 新規

- `src/fermion/ftensor.hpp` — ftensor 構造体+member API ミラー
- `src/fermion/fops.hpp` / `fops.cpp` — f-プリミティブ群
- `src/fermion/fermion_info.hpp` — FermionInfo
- `test/test_fermion_layer.cpp` — 単体+規約固定テスト
- `test/fermion/` — 自由フェルミオン ctest(入力生成+厳密値スクリプト)

### 変更

- `src/iTPS/iTPS.hpp`(FermionInfo メンバ)、`load_toml.cpp`(パリティ解析+ガード)、
  `tensors.cpp`(偶マスク初期化)、`simple_update.cpp` / `measure.cpp` / `onesite_obs.cpp` /
  `twosite_obs.cpp`(境界ラップ+ガード)、`core/simple_update.cpp` / `core/ctm.cpp` /
  `contract_itps_ctm/ctm.cpp`(ftensor インスタンス化追加)、`PEPS_Parameters`、ドキュメント

### 見積り

約 3,000 行(符号レイヤー+単体テスト 1,500〜2,300、統合+入力+統合テスト 1,000〜1,200)

### 実装順序(M1 内 5 段階)

1. 符号レイヤー TDD(TeNeS 統合なし、doctest のみ)
2. 規約固定テスト(fconj 確定)
3. simple update 移植+1D 極限テスト
4. CTMRG+測定移植+自由フェルミオン D 収束テスト
5. 入力スキーマ・ガード・ドキュメント

## 9. 将来課題(M2 以降)

- **奇パリティ占有配置のプロダクト初期化**(2026-08-15 記録): セルあたり偶数フェルミオンの
  占有配置は「奇仮想ボンド(次元1の奇セクター)で占有サイトをペアにするダイマー被覆」で
  厳密に表現できる(例: 全占有 |11…⟩ は 2×1 セル+奇偶交互ボンド、チェッカーボードは
  2×2 セル+束の素通し)。ダイマー被覆の自動生成を実装し、奇プロダクト初期化の
  エラーを解除する。見積り 150〜300 行。
- `tenes_std`/`tenes_simple` のフェルミオン対応(フェルミオンモデル定義、
  多ホップ MPO 分解への JW string 挿入)
- 距離 2 以上の相関関数(パリティ演算子挿入付き転送)・相関長のパリティ奇セクター
- full update / fast full update、simple update ゲージ固定、RSVD のパリティ対応
- 汎用 NxM カーネルの fermion 対応: 生成器(`misc/contraction/`)に**平面埋め込みを持つ
  中間表現**を導入(三角格子 iTPS 計画と同じ改修点。幾何一般で設計する)
- 有限温度(purification のアンシラパリティ)、実時間発展
- U(1) ブロックスパース化(実用 D 引き上げ。swap gate 資産は YASTN 型アーキテクチャで温存可)

## 10. 参考文献

- P. Corboz, R. Orús, B. Bauer, G. Vidal, PRB 81, 165104 (2010), arXiv:0912.0646 — swap gate 方式の原典(fPEPS+CTMRG)
- T. Barthel, C. Pineda, J. Eisert, PRA 80, 042333 (2009), arXiv:0907.3689 — 縮約コスト同一性
- Q. Mortier et al., SciPost Phys. 18, 012 (2025), arXiv:2404.14611 — 統一レビュー(graded 形式との等価性)
- B. Bruognolo et al., SciPost Phys. Lect. Notes 25 (2021), arXiv:2006.08289 — フェルミオン iPEPS+CTMRG 実装手引き
- M. M. Rams et al., SciPost Phys. Codebases 52 (2025), arXiv:2405.12196 — YASTN(swap gate 実装規約)
- 調査レポート: https://claude.ai/code/artifact/653868dd-dff9-4fb7-9333-98a8fe81b69e
