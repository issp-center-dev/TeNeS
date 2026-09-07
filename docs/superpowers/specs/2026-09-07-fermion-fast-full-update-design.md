# フェルミオン fast full update — 設計書

日付: 2026-09-07
ブランチ: `fermion`、HEAD = `e6ab497f`
実測: `work/fermion/fastfu/FINDINGS.md`(試作は `work/fermion/fastfu/spike.patch`、計測後に撤去済み)
関連: `docs/superpowers/specs/2026-09-05-fermion-ctm-phase-design.md`(CTM 環境の位相)、
`work/fermion/full-update-design/HANDOFF.md` の残課題 4(D=4 chi=16 の CTM 停滞)

## 1. 問題

`src/iTPS/main.cpp:250-261` はフェルミオン模式で `Full_Use_FastFullUpdate` を強制的に false にする。

```
WARNING: fermion mode disables Full_Use_FastFullUpdate because the fast update
reuses bare-Tn CTM moves that are not fermion-aware in this version
```

その結果、フェルミオンの full update はボンド 1 本ごとに `update_CTM()` を呼び、
`Calc_CTM_Environment_density(..., initialize = true, ...)` が毎回**一様ベクトル初期化から
CTMRG を収束させ直す**。

これが 2 つの問題を生んでいる。

1. **遅い。** 最終レビューの実測で D=4 d=4 chi=16 の 1 sweep が 61 秒、うち `environment` が 88%。
2. **simple update が収束していない状態から D >= 3 に入ると完走しない。**
   simple update 50 step の自由フェルミオン D=3 では chi=12 でも chi=24 でも、
   `build_full_update_environment` の forbidden parity ガードに引っかかって異常終了する。
   chi=12 では `rdm_dist` が 8.23657e-05 に固着し、`iteration_max` を 50 から 200 に増やしても
   **同じ値**から動かない。反復不足ではなく、毎ボンドのゼロスタートで CTMRG が別の固定点に入る。

   **ただし D >= 3 が常に落ちるわけではない**(初版の「完走した実績は D=2 chi=8 だけ」という
   記述を訂正)。既存 ctest `FreeFermionFull` は D=3 chi=12 を `fastfullupdate = false` で
   回して緑であり、違いは simple update の step 数(1000)である。同じ入力の SU を 1000 step に
   すると非 fast も完走する(177.65 秒、ただし 80 ボンド中 3 回は CTM 未収束の警告が出る)。
   fast はこの固着を踏まないが、「fast でなければ D >= 3 が動かない」ではない。

「bare-Tn の move を使うので fermion-aware でない」という警告の理由は、実は成り立たない。
フェルミオン経路の `update_CTM()` は

```cpp
const std::vector<ptensor> reduced_Tn =
    tenes::fermion::build_reduced_density_tensors(Tn, finfo);
core::Calc_CTM_Environment_density(..., reduced_Tn, ..., true, true);
```

と、**符号を折り込んだ reduced density tensor** を作って密度行列経路に渡している。
`Calc_CTM_Environment_density` は内部で `Make_single_tensor_density` して `core::*_move_single` を
回すだけで、move 自体はボゾンの有限温度経路と同一コードである。フェルミオン性は reduced tensor の
構築に閉じており、**move に渡すテンソルさえ reduced にすれば fast full update は成立する**。

位相についても追加の考慮は要らない。`phase_invariant` は `Check_Convergence_CTM_RDM` の中でしか
使われず、move 本体には位相処理がない。位相は消費側(FU の `build_full_update_environment`、
測定の `normalize_rdm_phase`)が自分の窓の ⟨1⟩ で決めるという既存の契約
(位相設計書 §2)がそのまま効く。

## 2. 実測の要約

試作(`spike.patch`)で確認済み。詳細は `work/fermion/fastfu/FINDINGS.md`。

| ケース | 現状 | fast | 速度比 | Energy 差(相対) |
|---|---|---|---|---|
| D=2 chi=8, FU 10 sweep | 78.2 s | 5.6 s | 14x | 5.8e-4 |
| D=3 chi=12, FU 10 sweep(SU 50 step) | throw | 7.4 s | — | 非 fast が落ちる |
| D=3 chi=24, FU 10 sweep(SU 50 step) | throw (283.9 s) | 14.2 s | — | 非 fast が落ちる |
| D=3 chi=12, FU 10 sweep(SU 1000 step) | 177.7 s | 未測定 | — | **A/B は実装後 (T5) に取る** |
| D=4 chi=16, FU 1 sweep | 91.6 s | 21.5 s | 4.3x | **2.4e-6** |

- **D=4 の 1 sweep で相対 2.4e-6 の一致**が fast 経路の正しさの直接の裏付け。
- forbidden parity ガードの誤発火という事前の懸念は外れた。fast の forbidden ratio は base より
  1〜2 桁**小さい**(D=2 で max 1.70e-10 対 3.60e-09、閾値は 1e-6)。**閾値の緩和は不要。**
- FU 中の move と reduced tensor 再構築の合計は、D=4 chi=16 の 1 sweep で 1.5 秒程度、
  base の full update 72.5 秒に対して 2%。**差分更新のキャッシュは不要。**

## 3. 方針

1. フェルミオン経路に、reduced density tensor を `core::*_move_single` に渡す fast full update を
   実装し、`Full_Use_FastFullUpdate` の強制 OFF を撤去する。既定値はボゾンと同じ true。
2. 非 fast 経路(`fastfullupdate = false`)の `update_CTM()` に warm start を入れ、
   2 回目以降は前回の環境から続きを収束させる。fast が使えない run の受け皿であり、
   D >= 3 の固着に対する保険でもある。

**前提にしてはならないこと**: フェルミオン模式が現在拒否している機能(skew セル、実時間発展、
1 幅セル、`Use_RSVD`、相関関数・相関長など)の多くは、以前の試みが失敗して先送りされたもので
あって恒久的な制約ではない。**それらの不在を前提とした実装・最適化・簡略化を入れない。**
本タスクで具体的に効くのは skew(§4.1)と実時間発展(§4.2)。

## 4. 設計

### 4.1 フェルミオン fast full update

`iTPS` に 1 メソッドを足す。

```cpp
//! Absorb the row/column touched by one full-update bond into the CTM
//! environment, the fermionic counterpart of the bosonic fast full update.
void update_CTM_fast_fermion(int source, int target, int source_leg);
```

実装(`src/iTPS/full_update.cpp`):

```cpp
template <class tensor>
void iTPS<tensor>::update_CTM_fast_fermion(int source, int target,
                                           int source_leg) {
  Timer<> timer;
  const std::vector<tensor> Tn_single = core::Make_single_tensor_density(
      tenes::fermion::build_reduced_density_tensors(Tn, finfo));
  // 引数の並びは 4 つの move すべてで
  //   (C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn_single, i, peps_parameters, lattice)
  // で共通(const 修飾だけが move ごとに違う)。以下の "..." はこの並びの省略。
  // 方向の対応はボゾン分岐(full_update.cpp:160-186)と同一。
  if (source_leg == 0) {
    core::Right_move_single(..., source % LX, ...);
    core::Left_move_single (..., target % LX, ...);
  } else if (source_leg == 1) {
    core::Bottom_move_single(..., source / LX, ...);
    core::Top_move_single   (..., target / LX, ...);
  } else if (source_leg == 2) {
    core::Left_move_single (..., source % LX, ...);
    core::Right_move_single(..., target % LX, ...);
  } else {
    core::Top_move_single   (..., source / LX, ...);
    core::Bottom_move_single(..., target / LX, ...);
  }
  time_environment += timer.elapsed();
}
```

呼び出し側は `full_update(up)` のフェルミオン分岐末尾:

```cpp
tenes::fermion::validate_neighbor_consistency(finfo, lattice);
if (peps_parameters.Full_Use_FastFullUpdate) {
  update_CTM_fast_fermion(source, target, source_leg);
} else {
  update_CTM(/*warm_start=*/true);
}
return;
```

注意点:

- **move には `source` / `target` / `source_leg` を渡す。** フェルミオン分岐は内部で
  `s1` / `s2` に入れ替える(`source_leg` が 0 か 1 のとき swap して gate を transpose する)が、
  これは Full_update_bond_fermion の引数規約であって格子上の位置ではない。
  行・列の決定には入れ替え前の値を使う。
- ボゾン分岐との共通化はしない。ボゾンは rank-5 の `Tn` を `core::*_move` に、フェルミオンは
  rank-4 の単層テンソルを `core::*_move_single` に渡すので、共有できるのは
  `source_leg` から move の組を選ぶ対応表だけであり、そのために型を捻じ曲げる価値はない。
  **ただし対応表は 2 箇所に重複するので、片方を直したらもう片方も直す** という注意をコメントに残す。
- reduced tensor は毎ボンド全サイト作り直す(§2 の実測)。
- move が **吸収する行・列と書き換える行・列は 1 セルずれる**(§5 の契約 4 の表)。渡す引数は
  「吸収させたい行・列」であり、ボゾン分岐と同じ `source % LX` / `source / LX` でよい。
  この非対称性は move 側の既存規約で、本タスクでは触らない。
- **skew を前提にした実装を入れないこと。** 添字はボゾン fast 分岐と同一の
  `source % LX` / `source / LX` を使い、skew の折り返しは move 側が `lattice.index()` /
  `top()` / `bottom()` を通じて処理する。`LY_noskew == LY` を仮定した最適化や、
  `LX_noskew` / `LY_noskew` を `LX` / `LY` に読み替える簡略化を書いてはならない
  (`skew != 0` では `LY_noskew = LY * lcm(LX, skew) / skew > LY`。`SquareLattice.cpp:63-70`)。

  現在フェルミオン模式は skew セルを入力読み込みで拒否するが、その理由は
  「skewed unit cells (measured to give wrong fermionic numbers)」(`load_toml.cpp:621-623`)、
  すなわち**測定がフェルミオン数を誤る**という CTM move とは無関係の問題であり、
  いずれ解禁される見込みの先送りである。解禁時に fast full update 側で追加作業が生じないよう、
  skew に依存しない形で書く。

  構成要素はどちらも skew 対応の実績を持つ。ボゾンの fast full update + skew は
  `test/data/Honeycomb_skew.toml`(skew = -1、full update 10 step、`fastfullupdate` 未指定 =
  既定 true)が、`*_move_single` + skew は `test/data/FT_Kitaev.toml`(skew = 1)が、
  それぞれ既存の回帰テストとして通っている。両者の組み合わせだけが未検証。
- 1 サイトの full ゲート(`up.is_onesite()`)は現状 CTM を更新せずに return する。
  この挙動は変えない(ボゾンも同じ)。

### 4.2 warm start つき `update_CTM()`

```cpp
//! @param warm_start Reuse the current environment as the CTMRG initial guess
//!        instead of rebuilding it from uniform vectors.  Ignored (treated as
//!        false) until the environment has been built at least once.
void update_CTM(bool warm_start = false);
```

`iTPS` に `bool ctm_valid_ = false;` を持たせ、

- `initialize_tensors()` / `initialize_tensors_density()` の末尾で false
- `load_tensors()`(`src/iTPS/saveload_tensors.cpp:172`、宣言は `iTPS.hpp:315`)の末尾で false。
  読み込んだ環境は現状どのモードでも捨てられるので、warm start でうっかり使わないため
- `update_CTM()` の末尾で true

とする。`initialize` に渡す値は `!(warm_start && ctm_valid_)`。

**状態フラグを持つ理由**: `full_update(up)` は `full_update()` からだけでなく
`time_evolution.cpp:54` からも呼ばれ、後者は冒頭に cold な `update_CTM()` を持たない。
「full_update() が先に環境を作っているはず」という前提を置くと、実時間発展 + full update の経路で
未初期化の環境から CTMRG を始めることになる。フェルミオンは実時間発展を入力読み込みで拒否するので
現時点では踏まないが、前提に依存しない形にしておく。

warm start を有効にするのは **フェルミオンの非 fast 経路だけ**とする。

- `full_update()` 冒頭の `update_CTM()`: 現状どおり cold。simple update 直後で Tn が大きく動いており、
  古い環境が良い初期値である保証がない。
- `measure()` の `update_CTM()`: 現状どおり cold。「独立に収束させた環境で測る」という
  現在の契約を変えない。
- ボゾンの非 fast 経路: 現状どおり cold。同じ恩恵を受けられるはずだが、
  `test/data/output_*/` の golden 値が収束経路の変更でずれるリスクを本タスクでは取らない。
  **ボゾンへの展開は別タスク**とし、NEWS に「フェルミオンのみ」と明記する。

### 4.3 フラグ・既定値・ドキュメント

- `src/iTPS/main.cpp:250-261` の強制 OFF と警告を撤去する。
- `Full_Use_FastFullUpdate` の既定値は true のまま。**フェルミオン模式の既定の挙動が変わる**
  (非 fast → fast)。非 fast が simple update 未収束の D >= 3 で固着することと、
  速度差(§2)を踏まえれば妥当な既定。
- `docs/sphinx/{ja,en}/file_specification/parameter_section.rst` の
  「``fastfullupdate = true``(既定値)は非対応です」の箇条書きを差し替える。
  `meanfield_env = true` との組み合わせがエラーである点は変えない。
- `NEWS.md` に、フェルミオン fast full update の追加と、非 fast 経路の warm start を書く。

## 5. 振る舞い契約(要点)

正式な契約書は別途 `2026-09-07-fermion-fast-full-update-contract.md` に書き、テスト作成者へ渡す。
要点のみ:

1. **完走**: フェルミオン模式で `fastfullupdate = true` を指定しても警告が出ず、false に落とされない。
2. **一致**: D=4 chi=16 の自由フェルミオンで FU 1 sweep を回したとき、
   `fastfullupdate` の true と false でエネルギーが相対 1e-5 以内で一致する。
3. **ボゾン等価**: parity を全て 0 にした fermion の fast full update が、同じ入力の boson の
   fast full update と一致する(`test/fermion/boson_equivalence_full.py.in` の枠組みを
   `fastfullupdate = true` に拡張)。これは非 fast 経路が完走しない D >= 3 でも成立する検証手段。
4. **局所性**: 各 `source_leg` について、呼ばれる move の種類と引数が正しいこと。
   4 つの move は書き換えるテンソルの種別が互いに素で、しかも **吸収する行・列と書き換える行・列が
   1 セルずれる**(`ctm_single.cpp:555-566, 618-633, 686-701, 754-771` で確認済み)。

   | move | 吸収する行・列 | 書き換えるテンソル | 書き換わる行・列 |
   |---|---|---|---|
   | `Left_move_single(ix)` | 列 `ix` | `C1` `C4` `eTl` | `right(ix)` |
   | `Right_move_single(ix)` | 列 `ix` | `C2` `C3` `eTr` | `left(ix)` |
   | `Top_move_single(iy)` | 行 `iy` | `C1` `C2` `eTt` | `bottom(iy)` |
   | `Bottom_move_single(iy)` | 行 `iy` | `C3` `C4` `eTb` | `top(iy)` |

   種別の集合 {C1,C4,eTl} / {C2,C3,eTr} / {C1,C2,eTt} / {C3,C4,eTb} は 4 つとも異なるので、
   1 ボンドの fast 更新の前後で環境テンソルのスナップショットを取れば、**変化した種別の集合から
   move の種類が一意に決まり**、変化した添字から行・列引数が逆算できる。

   **テストは「吸収した行・列が変わる」と書いてはならない。** 実際に変わるのは隣のセルであり、
   そう書くと正しい実装を赤くする。期待値は上の表の右 2 列から作ること。

   `source_leg` と move の対応表を 1 つずらす変異、および `source` と `target` を入れ替える変異が
   赤くなること。後者は §4.1 の「入れ替え前の `source` / `target` を使う」という規約の防衛にあたる。
5. **パリティ**: fast 経路で `build_full_update_environment` の forbidden ratio が閾値を超えない。
6. **warm start**: 非 fast 経路で warm start あり・なしの最終結果が CTM の収束許容誤差内で一致する。
   未初期化の環境に対して warm start を要求しても cold start にフォールバックする。
7. **skew 非依存**: 実装に `skew == 0` / `LY_noskew == LY` を仮定した分岐や簡略化が無いこと。
   フェルミオン模式が skew を拒否している以上フェルミオンでの直接のテストは書けないので、
   これはレビューで担保する項目とし、`skew` および `*_noskew` を含む行を実装差分から
   拾って確認する。ボゾンの `Honeycomb_skew` と有限温度の `FT_Kitaev` が引き続き緑であることも
   併せて確認する(共有コードに触るため)。

## 6. タスク分割

1. **T1** `update_CTM_fast_fermion` の追加と `full_update(up)` からの呼び出し。
2. **T2** `main.cpp` の強制 OFF 撤去。
3. **T3** `update_CTM(bool warm_start)` と `ctm_valid_` の導入、非 fast 経路からの呼び出し。
4. **T4** ドキュメント(ja/en)と `NEWS.md`。
5. **T5** 検証: D=2 / D=3 / D=4 の A/B と ctest 全件。

T1 と T3 は独立。T5 は Claude が独立に実施する(Codex の報告を鵜呑みにしない)。
テストは契約書のみを渡したテスト作成者が書き、RED を 1 件ずつ確認してスナップショットを取ってから
実装に入る。タスクごとにサブエージェントでレビューし、最後に最上位モデルで全ブランチレビューを 1 回。
共有コード(`update_CTM` のシグネチャ変更)に触るので、検証はテストスイート全件を Claude が回す。

## 7. 非目標

- reduced density tensor の差分更新キャッシュ(§2 の実測で不要と判断)。
- ボゾン経路の warm start(§4.2)。
- 1 サイト full ゲートでの CTM 更新(現状の挙動を維持)。
- D=3 で非 fast 経路の CTMRG が `rdm_dist` の別固定点に入る現象そのものの解決
  (HANDOFF 残課題 4)。fast はこれを回避するが、原因の究明は本タスクの範囲外。
- MPI 環境での性能・正当性の確認(HPC で別途。HANDOFF minor-3)。
- `core::*_move_single` 自体の添字規約(吸収する行・列と書き換える行・列が 1 セルずれること)の変更。
  本タスクはボゾン fast 分岐と同じ move API を同じ規約で呼ぶだけで、規約を直す作業ではない。
- フェルミオン模式での skew セルの解禁(`load_toml.cpp:621-623`)。これは測定がフェルミオン数を
  誤る問題であって本タスクの範囲外。ただし **§4.1 のとおり、fast full update 側に skew を前提と
  した実装を残さない**こと。解禁作業が fast full update の書き直しを伴ってはならない。

## 8. リスクと未解決

- **fast は CTM の収束判定を通らない。** 環境は「1 move 動かしただけ」であり、収束の保証がない。
  これはボゾンの fast full update と同じ性質で、fast full update の定義そのものである。
  D=4 の A/B で相対 2.4e-6 の一致が取れているが、より大きな D / chi での劣化は未確認。
  §5 の契約 3(ボゾン等価)が、非 fast 経路の完走に依存しない検証手段として効く。
- **D=2 の A/B 差が相対 5.8e-4 と大きい。** 10 sweep 分の蓄積であり、かつ D=2 は
  「FU がエネルギーを上げる」病的領域(`parameter_section.rst:68`)。1 sweep 当たりで見れば
  D=4 と整合するが、契約 2 の許容誤差を D=2 に適用してはならない。
- **警告文の削除で、非 fast にフォールバックしていた既存ユーザーの結果が変わる。**
  NEWS に明記する。
