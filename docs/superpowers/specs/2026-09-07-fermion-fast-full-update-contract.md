# フェルミオン fast full update — 振る舞い契約(2026-09-07)

設計書: `docs/superpowers/specs/2026-09-07-fermion-fast-full-update-design.md`
事前実測: `work/fermion/fastfu/FINDINGS.md`(試作での A/B、再現手順つき)

## 0. テスト作成者への注意

- **実装はまだ無い。** この契約だけを読んでテストを書く。実装を先に読んで、その構造に
  合わせたテストを書かないこと。
- 契約に書かれた性質のうち、**根拠を自分で確かめられるもの(既存コードの現在の振る舞い、
  既存テストの有無、ファイル名や行番号)は確かめてから使うこと。** 契約が間違っていたら
  その旨を報告する。過去に契約書の誤りをテスト作成者が見つけた例が複数ある。
- **参照値は「別経路の計算結果」か「解析的に決まる値」から取る。** この契約の中心は
  「2 つの独立な経路が一致すること」であって、特定の数値を焼き込むことではない。
- テストは `test/fermion/` に置く。C++ 単体テストは doctest、E2E は `.py.in` を
  `test/CMakeLists.txt` から configure して登録する既存の作法に従う。
- **実装がまだ無いので、書いたテストは失敗するかコンパイルが通らないかのどちらかになる。
  それが正しい状態である。** どのテストがどちらになるかを報告書に書くこと。
  通ってしまうテストがあれば、それは何も検証していない疑いがあるので報告すること。

## 1. 何が変わるか

フェルミオン模式では現在 `src/iTPS/main.cpp` が `Full_Use_FastFullUpdate` を強制的に false に
落とし、警告を出す。そのため full update はボンド 1 本ごとに `iTPS::update_CTM()` を呼び、
CTM を一様ベクトル初期化から収束させ直している。

変更後は、フェルミオン模式でも `fastfullupdate`(既定 true)が効くようになる。fast のとき、
1 ボンドの更新後に行うのは **CTM 全体の再収束ではなく、更新された 2 サイトを含む行または列を
1 回だけ吸収する部分 move 2 回**である。

あわせて、非 fast 経路(`fastfullupdate = false`)の `update_CTM()` に warm start が入り、
2 回目以降は現在の環境を CTMRG の初期値に使う。

## 2. 前提として確認済みの事実

以下は実装前に確認済み。テスト作成者も裏を取ってよい。

- フェルミオン経路の `iTPS::update_CTM()` は `tenes::fermion::build_reduced_density_tensors(Tn, finfo)`
  で符号を折り込んだテンソルを作り、`core::Calc_CTM_Environment_density` に渡している。
  `Calc_CTM_Environment_density` は内部で `core::Make_single_tensor_density` を通してから
  `core::*_move_single` を回すだけで、move 自体はボゾンの有限温度経路と同一コードである。
  **フェルミオン性は reduced tensor の構築に閉じており、move には入っていない。**
  fast 経路が move に渡すテンソルも同じ
  `Make_single_tensor_density(build_reduced_density_tensors(Tn, finfo))` である。
- `phase_invariant` は `Check_Convergence_CTM_RDM` の中でしか使われない。move 本体に位相処理はない。
- 4 つの move は、**吸収する行・列と書き換える行・列が 1 セルずれる**。

  | move | 吸収する行・列 | 書き換えるテンソル | 書き換わる行・列 |
  |---|---|---|---|
  | `Left_move_single(ix)` | 列 `ix` | `C1` `C4` `eTl` | `right(ix)` |
  | `Right_move_single(ix)` | 列 `ix` | `C2` `C3` `eTr` | `left(ix)` |
  | `Top_move_single(iy)` | 行 `iy` | `C1` `C2` `eTt` | `bottom(iy)` |
  | `Bottom_move_single(iy)` | 行 `iy` | `C3` `C4` `eTb` | `top(iy)` |

  根拠は `src/iTPS/core/ctm_single.cpp` の各 move の後半(更新ループ)。**「吸収した行・列が
  変わる」と書いたテストは正しい実装を赤にする。**
- 4 つの「書き換えるテンソル種別の集合」{C1,C4,eTl} / {C2,C3,eTr} / {C1,C2,eTt} / {C3,C4,eTb}
  は互いに異なる。したがって環境テンソルのスナップショット差分から **move の種類が一意に決まる**。
- 試作での A/B(`work/fermion/fastfu/FINDINGS.md`):
  D=4 chi=16 の自由フェルミオンで full update 1 sweep を回すと、fast と非 fast の
  エネルギーが**相対 2.4e-6** で一致した(base 92 秒 / fast 22 秒)。
  同じ系で D=2 chi=8 の 10 sweep では相対 5.8e-4 まで開く。D=2 は full update が
  エネルギーを上げる病的領域であり、sweep 数を積むほど差が開く。
- **simple update が収束していない状態から D >= 3 の full update に入ると**、非 fast 経路は
  毎ボンドの cold start で CTMRG が別の固定点に固着し、forbidden parity ガードで異常終了する
  (SU 50 step の D=3 では chi=12 / chi=24 のいずれでも落ちる)。
  **ただし D >= 3 が常に落ちるわけではない**(初版の記述を 2026-09-08 に訂正)。既存 ctest
  `FreeFermionFull` は D=3 chi=12 を `fastfullupdate = false` で回して緑であり、違いは
  simple update の step 数(1000)。同じ入力の SU を 1000 step にすると非 fast も完走する
  (177.65 秒、80 ボンド中 3 回は CTM 未収束の警告が出る)。

## 3. 要求

### R1: `fastfullupdate = true` がフェルミオンで有効になる

フェルミオン模式かつ `full_update` の step 数が 1 以上のとき、`fastfullupdate = true` を
指定しても

- `Full_Use_FastFullUpdate` を false に落とさない
- `Full_Use_FastFullUpdate` に関する警告を出さない
- run が正常終了する

`meanfield_env = true` との組み合わせが入力読み込みでエラーになる現在の振る舞いは変えない。

**既存テストとの矛盾**: `test/fermion/fermion_guards.cpp` の T6 は現在「fallback 警告が
**出ること**」を要求している。R1 はこれを禁じるので、実装が入るとこの既存テストが赤になる。
実装者はテストファイルを変更できないので、**テスト作成者が反転させた**(2026-09-08 対応済み)。
実装後もここが赤のままなら、警告が別の場所からも出ている。

### R2: fast のときに呼ばれる move

`source_leg` を、更新するボンドの source 側の脚番号(0 = 左, 1 = 上, 2 = 右, 3 = 下)、
`source` / `target` をボンドの 2 サイトとする。fast のとき、1 ボンドの更新後に呼ばれるのは
次の 2 つの move だけである。

| `source_leg` | 1 番目 | 2 番目 |
|---|---|---|
| 0(左) | `Right_move_single(source の列)` | `Left_move_single(target の列)` |
| 1(上) | `Bottom_move_single(source の行)` | `Top_move_single(target の行)` |
| 2(右) | `Left_move_single(source の列)` | `Right_move_single(target の列)` |
| 3(下) | `Top_move_single(source の行)` | `Bottom_move_single(target の行)` |

これはボゾンの fast full update が呼ぶ move と 1 対 1 で同じである(ボゾンは
`core::*_move`、フェルミオンは `core::*_move_single` を呼ぶ点だけが違う)。

**注意**: フェルミオンのボンド更新は内部で 2 サイトを canonical な向きに入れ替えることがある
(`source_leg` が 0 か 1 のとき)。move に渡す行・列は**入れ替え前の** `source` / `target` から
決めなければならない。入れ替え後の値を使うと、`source_leg` 0 と 2、1 と 3 で move が逆になる。

CTM の再収束(`update_CTM()` 相当)は fast のときには行われない。

### R3: 非 fast 経路の warm start

`fastfullupdate = false` のフェルミオン run で、ボンド更新後の CTM 再収束は
**現在の環境を初期値として**行う。ただし環境がまだ一度も構築されていない場合は
一様ベクトル初期化(cold start)にフォールバックする。

warm start を使うのはフェルミオンのボンド更新後の再収束だけである。以下は現状どおり cold:

- full update 開始時の 1 回目の環境構築
- `measure()` の環境構築
- ボゾン経路のすべての環境構築

環境が有効かどうかの状態は、環境テンソルの形が変わりうる操作(テンソルの初期化、
チェックポイントの読み込み `iTPS::load_tensors()`)で無効に戻ること。

**シグネチャは `void update_CTM(bool warm_start)` とすること**(既定値 `= false` を付けてよい)。
enum や `update_CTM_warm()` のような別の綴りにしない。テストがこの形を直接呼ぶ。

### R4: skew を前提にしない

実装に `skew == 0` や `LY_noskew == LY` を仮定した分岐・簡略化を入れないこと。
フェルミオン模式は現在 skew セルを拒否するが、その理由は測定がフェルミオン数を誤ることであり、
CTM move とは無関係の先送りである。

`skew != 0` では `LY_noskew` は `LY` より大きくなりうる。move に渡す行・列は
ボゾンの fast full update と同じ計算(unit cell 内の x 座標・y 座標)でよく、
skew の折り返しは `SquareLattice` 側が処理する。

### R5: ドキュメント

`docs/sphinx/ja/file_specification/parameter_section.rst` と英語版の対応箇所にある
「`fastfullupdate = true`(既定値)は非対応です」の記述を、実態に合わせて書き換える。
`NEWS.md` に、フェルミオン fast full update の追加と、非 fast 経路の warm start が
フェルミオンのみであることを書く。

## 4. テストへの要求

### S-1: fast と非 fast が一致する(E2E)

自由フェルミオンで、`fastfullupdate` の true と false 以外はまったく同じ入力を 2 つ作り、
エネルギーが一致することを確認する。

- **条件の基準は D=4 chi=16、full update 1 sweep**。試作での実測は相対 2.4e-6、
  所要は base 92 秒 / fast 22 秒。既定の ctest に入れるには重いので、
  `TENES_FERMION_FULL_TEST` / `TENES_FREE_FERMION_FULL` と同じ切り替えの下に置いてよい。
- **より軽い条件を選ぶ場合は、その条件で実際に A/B を回して差を実測し、許容誤差を
  実測から決めること。** 契約の 1e-5 をそのまま軽い条件に流用しない。
- **D=2 を使うなら full update は 1 sweep に絞ること。** 10 sweep 積むと相対 5.8e-4 になる。
- **非 fast 側の CTM が収束していることを要求しないこと。** 非 fast 経路はボンドごとの
  cold start のうち 1 本が別の固定点に固着することがあり(D=2 chi=8 / SU 200 step では
  `iteration_max` を 400、`convergence_epsilon` を 1e-10 にしても `rdm_dist` が 7.06e-4 で
  固着する)、これは **R3 が緩和しようとしている病理そのもの**であって S-1 の主題ではない。
  fast 側の収束は要求してよい。非 fast 側は回数を印字するに留めること。
- **D >= 3 でも SU を十分回した状態なら非 fast は完走する**ので A/B は取れるが、
  条件によっては上の固着を踏む。条件は実測で選ぶこと。

### S-2: grading が自明なとき boson と一致する(E2E)

`test/fermion/boson_equivalence_full.py.in` は「parity をすべて 0 にした fermion の full update は
boson の full update に一致する」ことを `fastfullupdate = false` で検証している。
**同じ検証を `fastfullupdate = true` でも行うこと。**

これは S-1 と違い、非 fast 経路の完走に依存しない。フェルミオン fast 経路
(reduced tensor + `*_move_single`)とボゾン fast 経路(bare Tn + `*_move`)という
**別実装どうしの比較**である点に価値がある。許容誤差は既存テストが
`fastfullupdate = false` で使っているものを出発点にし、必要なら実測で調整する。

### S-3: 呼ばれる move が正しい(単体)

R2 の表のとおりに move が呼ばれることを、**環境テンソルの変化から**確認する。
実装に「どの move を呼んだか」を問い合わせるのではなく、外から観測すること。

§2 の表のとおり、4 つの move は書き換えるテンソル種別の集合が互いに異なり、
書き換わる行・列は吸収した行・列の隣である。1 ボンドの fast 更新の前後で
`C1` `C2` `C3` `C4` `eTt` `eTr` `eTb` `eTl` のスナップショットを取れば、
どの move がどの引数で呼ばれたかを復元できる。

**アクセスの手当て**: `iTPS::full_update(up)` と `iTPS::update_CTM()` は public なので直接呼べるが、
環境テンソル `C1` `C2` `C3` `C4` `eTt` `eTr` `eTb` `eTl` は private である
(`src/iTPS/iTPS.hpp` の 318 行目以降)。テストからは friend の `iTPSTestAccessor` 経由で読む。
この struct は `test/test_fermion_layer.cpp` と `test/input.cpp` にそれぞれ定義があり、
現在は `Tn` と `finfo` のアクセサしか持たないので、**環境テンソル用のアクセサを足すこと**。
`iTPS` 本体に public なアクセサを追加してはならない(テストのためだけに公開範囲を広げない)。

4 つの `source_leg` すべてについて確認すること。次の変異が赤くなること:

- **(a)** `source_leg` と move の対応表を 1 つずらす(例: leg 2 で Left と Right を入れ替える)
- **(b)** move に渡す行・列に、入れ替え**後**の 2 サイトを使う
- **(c)** 2 回の move のうち 1 回を落とす

### S-4: 警告が出ない(E2E)

R1 の 3 条件を確認する。とくに **`Full_Use_FastFullUpdate` に関する警告が標準エラーに
出ないこと**。既存の `boson_equivalence_full.py.in` に、`fastfullupdate = false` を
指定したのに警告が出たら失敗にする検査があるので、その逆向きの検査になる。

### S-5: warm start(単体)

R3 により `iTPS::update_CTM()` は warm start を要求する引数を取るようになる。public なので
テストから cold / warm の両方を呼び分けられる。**これが A/B の取り方である**
(入力パラメータからは制御できないので E2E では書けない)。

- 同じ状態から cold と warm で環境を収束させ、結果が CTM の収束許容誤差内で一致すること。
  simple update 後に一度環境を作り、`Tn` をわずかに変えてから両方を比べる形でよい。
  どちらも収束する限り同じ固定点に行くはずである。
- 環境が一度も構築されていない状態で warm start を要求しても cold start に落ちて
  正常に動くこと(未初期化の環境から CTMRG を始めて壊れないこと)。
  これは warm start の**フォールバックが効いているか**を見るもので、
  「たまたま動いた」で通らないよう、結果が cold start と一致することまで確認する。
- ボゾン run の結果が warm start の導入で変わらないこと。**既存のボゾン回帰テスト
  (`AntiferroHeisenberg_real` など)がそのまま緑であること**で担保してよい。

### S-6: forbidden parity ガードに掛からない(S-1 に含めてよい)

fast 経路で `build_full_update_environment` の forbidden parity 比が閾値
(`max(1e-8, 100 * CTM_Convergence_Epsilon)`)を超えないこと。
`TENES_FERMION_FULL_UPDATE_LOG` を設定すると `N_forbidden_ratio=` の行が標準エラーに出る
(`src/iTPS/core/full_update_fermion.cpp` の既存の診断)。試作では fast の値は
非 fast より 1〜2 桁小さかった(D=2 で max 1.7e-10 対 3.6e-9)。

## 5. 検出力の限界(重要)

- **fast full update は近似である。** fast と非 fast が厳密に一致することはない。
  S-1 が確認できるのは「近似が想定の範囲に収まっていること」だけで、
  「fast 経路が正しい」ことの証明ではない。S-2 の boson 等価が、別実装との比較という
  意味でより強い。
- **skew はフェルミオンでテストできない**(入力読み込みで拒否される)。R4 は
  レビューで担保する。実装差分から `skew` と `noskew` を含む行を拾って確認すること。
  ボゾンの `Honeycomb_skew` と有限温度の `FT_Kitaev` が緑のままであることも確認する
  (どちらも skew を使い、後者は `*_move_single` を通る)。
- **fast が「非 fast より良い」ことは直接は示せない。** 速度差は測れるが、
  「毎ボンドの cold start が踏む固着を fast は踏まない」ことは、固着する条件を選んで
  はじめて見える性質であり、安定した回帰テストにしにくい。実装後の検証 (T5) で
  Claude が A/B を取る項目とする。
- **1 サイトの full ゲートは CTM を更新しない**(現状の振る舞い、ボゾンも同じ)。
  この経路は本契約の対象外で、既存のテストも無い。

## 6. 契約書チェックリスト(作成者が自分で確認する)

- 参照値が「別経路の計算結果」か「解析的に決まる値」から来ているか。実装自身の出力を
  参照値にしていないか。
- 期待値が §2 の表の**右 2 列**(書き換わる行・列)から作られているか。吸収した行・列を
  期待値にしていないか。
- S-3 の変異 (a)(b)(c) を実際に当てて、赤くなることを確認したか。
  ガード節を消したコピーでテストが赤くなるかを見る形でよい。
- 許容誤差を、その条件での実測から決めたか。別条件の数値を流用していないか。
- この契約に書かれていない前提(フェルミオンが現在拒否している機能の不在など)を
  テストが仮定していないか。
