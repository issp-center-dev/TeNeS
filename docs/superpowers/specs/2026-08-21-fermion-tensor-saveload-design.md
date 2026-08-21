# フェルミオン模式の tensor_save / tensor_load 対応 設計書

日付: 2026-08-21
ブランチ: `fermion`
前提: `docs/superpowers/specs/2026-08-15-fermion-m1-design.md`(C++ 側の graded 規約)、
`docs/superpowers/specs/2026-08-21-fermion-meanfield-measure-design.md`(平均場環境対応)、
`docs/superpowers/notes/2026-08-19-fermion-implementation-guide.md`(実装解説)。

改訂1 (2026-08-21): Codex レビューと scratchpad のスパイクを反映。主な変更は
(a) 格子の同一性(`L_sub`/`skew`)も保存・照合する(§3.1、§4 の V6b — `validate_neighbor_consistency`
だけでは同じ `N_UNIT` の別トポロジを通してしまう)、(b) 検証を**テンソル読み込みの前後に分割**
(`resize_tensor` が先に走ると次元不一致を検出する前にデータが壊れる、§4/§5.3)、
(c) V7 のパリティ違反は `allreduce_max` で全プロセスにまたがって取る(§4)、
(d) 分割が既定のままでも壊れる経路の実測を追加(§2.4)。

## 1. 目的とスコープ

フェルミオン模式では `[parameter.general] tensor_save` / `tensor_load` が入力読み込み時に
拒否される(`src/iTPS/load_toml.cpp:651-653`)。本設計は **モデルパラメータを変えながらの
連続計算(スキャン)** のために、フェルミオン模式でチェックポイントの保存・再開を
できるようにする。想定する運用は「(t, U, μ₁) で収束させて保存 → (t, U, μ₂) の初期状態として
読み込む」であり、**仮想ボンド次元 D は保存時と同じ**とする(ユーザー決定、2026-08-21)。

### やること

- `finfo.virt`(仮想脚の偶奇台帳)を保存ディレクトリに永続化し、読み込み時に復元する
- 読み込み時の整合性検証(下記 §4)。すべて**エラーで止める**(黙って間違わない)
- `validate_fermion_constraints` の `tensor load/save` 拒否を外す
- ドキュメント(ja/en の非対応一覧、`tensor_save` / `tensor_load` の説明)と `NEWS.md`

### やらないこと

- **D(`virtual_dim`)の変更を伴う再開**。保存時と異なる `virtual_dim` はエラーにする(§3.3)
- CTM 環境テンソルの再利用(§2.3 のとおり、そもそも保存されていないし再計算される)
- 有限温度・実時間発展での再開(フェルミオン模式自体が非対応のまま)
- ボソン経路の挙動変更。ボソンの「次元が違えば警告して resize」は現状のまま

## 2. なぜガードを外すだけでは足りないか

### 2.1 `finfo.virt` は計算中に変化する状態である

`FermionInfo`(`src/fermion/fermion_info.hpp:32-36`)は次を持つ:

| メンバ | 由来 | 保存の要否 |
|---|---|---|
| `enabled` | `[parameter.general] fermion` | 不要(入力から再構成) |
| `phys` | `[[tensor.unitcell]] parity` | 不要だが**照合に使う**(§4) |
| `virt` | 初期値 `even_first_parity(D)` → **simple update で更新** | **必要** |

`virt` は `src/iTPS/simple_update.cpp:111-112` で毎ボンド更新される。更新後の値は
`svd_trunc`(`src/fermion/fops.hpp:599-700`)が返す内部脚のパリティで、**偶奇セクターを
またいで特異値上位 D 本を選ぶ**ため、偶奇の内訳は原理的にデータ依存である。

**実測(`TENES_FERMION_SECTOR_LOG=1`、Release)**: 実際に試した 5 構成
(spinless μ=−1 の D=2/D=3、spinless μ=0 の D=4、Hubbard の D=3/D=4)では、
分割は最後まで `even_first_parity(D)` のままだった。初期値が `even_first_parity` で、
更新がその内訳を保つ(自己保存的である)ためと見られる。したがって
**今日のコードに限れば** 台帳を復元しない実装でも偶然動く。

しかし分割は「保たれる」のであって「既定に戻る」のではない。初期分割を 3/1 に変えた
ビルドで D=4 を 300 ステップ回すと分割は 3/1 のまま保たれ、その保存を既定(2/2)の
ビルドで読み込むと **エラーも警告も出さずにエネルギーが −0.6355 → −0.5895(7.2% のずれ)**
になった(§2.4)。初期化を変える将来の変更(U(1) 疎化、別の初期状態、非既定分割からの
再開)や、`svd_trunc` が実際に内訳を動かすパラメータ領域が一つでもあれば、同じ壊れ方をする。
台帳の保存は数百バイトのコストで、この経路を構造的に塞ぐ。

### 2.2 復元しないと黙って壊れる

保存されるのは `Tn` と `lambda` だけ(`src/iTPS/saveload_tensors.cpp:60-81`)。読み込み後に
`finfo.virt` が既定の `even_first_parity(D)` のままだと、テンソルの中身(どの添字が奇か)と
ラベルがずれる。結果:

- graded 演算(`tensordot` / `conj` / `transpose`)が誤った Koszul 符号を掛ける
- `enforce_even_parity`(`src/iTPS/core/simple_update.cpp`)が、ラベル上は奇総和に
  なった**実振幅をゼロに落とす**

いずれも例外にならず、値だけが間違う。M1 で根絶したバグクラスと同型である。

### 2.4 スパイク(scratchpad、2026-08-21 実測)

ガードだけを外した Release ビルドで確認した:

| 実験 | 結果 |
|---|---|
| spinless D=2 μ=0 を 200 step 保存 → `num_step = 0` で読み込み | `density.dat` が**ビット単位で一致**。テンソル・λ 以外に失われる状態は無い |
| 保存ディレクトリの中身 | `T_*` / `El,Et,Er,Eb_*` / `C1..C4_*` / `lambda_*` / `params.dat`。`eT` は rank 4 のままで読み込みも通る(§2.3 の前提を確認) |
| 初期分割 3/1 の D=4 を 300 step 保存 → 既定(2/2)ビルドで読み込み | エネルギー −6.35467e-01 → −5.89509e-01(7.2%)。1サイト密度は 5.053754e-01 → 5.053754e-01 で 7 桁一致 |

最後の行が本設計の存在理由である。**密度は合ったまま結合エネルギーだけがずれる**ので、
利用者が気づける見た目の異常が出ない。なお §4 の V7(パリティ違反の実測)は、
このケースを実際に捕まえる: 3/1 で偶だった要素が 2/2 のラベルでは奇総和になるため。

### 2.3 CTM 環境は保存対象ではない(現状の帰結)

`run_groundstate`(`src/iTPS/main.cpp:150-158`)は `optimize()` → `save_tensors()` →
`measure()` の順に呼ぶ。フェルミオン模式では full update が禁止されているので
`optimize()` は simple update のみで CTM に触れず、CTM 環境は `measure()` の中で
`Calc_CTM_Environment_density(..., initialize = true)`(`src/iTPS/measure.cpp:45-56`)により
**毎回ゼロから作り直される**。したがって

- 保存される `C*` / `eT*` は `initialize_tensors()` が確保しただけの未使用テンソルで、
  形は iTPS 形(`eTt` は rank 4: `(CHI, CHI, D, D)`)のまま
- 読み込んでも `measure()` が上書きするので、値は結果に影響しない

この結論は**基底状態モードのフェルミオン実行**(現在フェルミオンで許されている唯一のモード。
有限温度・実時間・full update は `load_toml.cpp:636,639` のガードで禁止)に限った話である。

**設計上の判断**: 環境テンソルの保存・読み込みコードは**変更しない**(フェルミオン用の
分岐を足さない)。理由は (a) 差分を最小に保つ、(b) 形が rank 4 のまま一貫しており
`load_tensor` の rank 検査を通る(§2.4 で実測)、(c) 将来フェルミオンで full update や
CTM のウォームスタートを入れる場合に配管が残っている。
なお本設計で外すのは `tensor load/save directories` のガード **だけ** であり、
モード・full update・その他のフェルミオンガードはすべて残す。

**この判断は §2.3 の呼び出し順序に依存する**(もし将来 `measure()` の後に保存するように
変わると、`eT` は density 形の rank 3 になって読み込み時に rank 不一致で落ちる)。
この依存を `save_tensors()` にコメントで明記し、E2E テスト(§5 層4)で往復を固定する。

## 3. ファイル形式

### 3.1 `<save_dir>/fermion.dat`(新規、フェルミオン模式のときだけ書く)

`params.dat` と同じ「値 + `# コメント`」形式。`params.dat` の形式バージョンは **1 のまま**
(ボソンの保存ディレクトリは一切変わらない)。

```
1 # Fermion_Format_Version
4 # N_UNIT
2 2 # L_sub
0 # skew
0 1 1 0 # parity of the physical leg of Tn[0]
0 1 # parity of the virtual leg 0 of Tn[0]
0 1 # parity of the virtual leg 1 of Tn[0]
0 1 # parity of the virtual leg 2 of Tn[0]
0 1 # parity of the virtual leg 3 of Tn[0]
0 1 1 0 # parity of the physical leg of Tn[1]
...
```

- 各パリティ行は `0`(偶)/`1`(奇)を空白区切りで、その脚の次元と同じ個数だけ並べる
- サイトは 0 から `N_UNIT-1` の順、各サイトは物理脚 → 仮想脚 0,1,2,3 の順(計 5 行)
- 書き出しは `mpirank == 0` のみ(`lambda_*.dat` と同じ流儀)
- `L_sub` と `skew` を持つのは、パリティ台帳が「どのサイトのどの脚か」でしか意味を持たず、
  同じ `N_UNIT` でも格子のトポロジが変われば `Tn[i]` と `lambda[i][leg]` の解釈が変わるため。
  `validate_neighbor_consistency` は現在の隣接関係で台帳の一致を見るだけなので、
  台帳がたまたま一致する別トポロジを通してしまう(Codex レビュー指摘)。
  `params.dat`(形式バージョン 1)は `N_UNIT` / `CHI` / 各テンソルの形しか持たない

パリティベクトルを丸ごと書くのは、`svd_trunc` が返す台帳が現状「偶が先、奇が後」に
整列している(`fops.hpp` の 2 段目の `stable_sort`)からといって、偶の個数だけを保存すると
その不変条件に依存してしまうため。全要素を書けば形式は自己記述的で、将来
整列規約が変わっても壊れない。

### 3.2 `params.dat` は変更しない

形式バージョン 1 を維持する。フェルミオンかどうかは `fermion.dat` の有無で判定する。

### 3.3 D 変更を許さない

`fermion.dat` の各パリティ行の長さが入力の `virtual_dim` / `physical_dim` と違えばエラー。
メッセージには「フェルミオン模式の再開は保存時と同じ `virtual_dim` を要求する。
D を変える場合はテンソルを読み込まずに新規に始めること」を含める。

理由: `load_tensor` の `resize_tensor`(`saveload_tensors.cpp:133`)は各脚の**末尾に
ゼロ詰め/切り詰め**するだけで、偶ブロック・奇ブロックの位置を保たない。保存された
台帳 `[偶×a, 奇×b]` を D_new に伸ばすと、新しいスロットの偶奇をどう定義しても
既存データの位置は動かないため、偶奇を保った埋め込みには**並べ替えを伴う散布**が要る
(実装可能だが本設計の範囲外。ユーザー判断で v1 は同一 D に限定)。

## 4. 読み込み時の検証(すべてエラー)

検証は **テンソル読み込みの前と後の 2 段**に分ける。前段が必要なのは、`load_tensor` が
`A = resize_tensor(temp, A.shape())`(`saveload_tensors.cpp:133`)で即座に代入し、
`load_tensors_v1` が `lambda_tensor[i][j].resize(vdim[j])`(同 270 行付近)で λ も詰め直すため。
`resize_tensor`(`src/tensor.cpp:139`)は各脚の末尾で extend / slice するので、
**次元不一致を検出する前に読み込んだデータが壊れる**(Codex レビュー指摘)。

### 前段(テンソルを読む前)

| # | 条件 | メッセージの要点 |
|---|---|---|
| V1 | フェルミオン模式で `fermion.dat` が無い | この保存はフェルミオン模式のものではない |
| V2 | `Fermion_Format_Version != 1` | 未知の形式バージョン |
| V3 | `N_UNIT` 不一致 | 入力と保存でユニットセルのサイト数が違う |
| V4 | 物理脚パリティが入力の `parity` と不一致(長さ・内容) | 物理基底が違う |
| V5 | 仮想脚パリティの長さが入力の `virtual_dim` と不一致 | §3.3 の D 固定。「保存時と同じ `virtual_dim` が必要」と書く |
| V6a | 復元後 `validate_neighbor_consistency(finfo, lattice)` が失敗 | 隣接ボンドで台帳が食い違う |
| V6b | `L_sub` / `skew` が入力と不一致 | 格子のトポロジが違うので同じ台帳でも意味が変わる |
| V8 | **ボソン実行**(`fermion = false`)で `fermion.dat` が存在 | フェルミオン模式の保存なので `fermion = true` が必要 |

前段が通ったら `finfo.virt` を復元してからテンソルを読む。

### 後段(テンソルを読んだ後)

| # | 条件 | メッセージの要点 |
|---|---|---|
| V7 | 読み込んだ `Tn[i]` にパリティ違反要素がある | 台帳とテンソルが不整合 |

- V7 は `tenes::fermion::parity_violation(wrap_Tn(Tn[i], finfo, i))`(`fops.hpp:266-279`)で
  実測する。ただしこの関数は `local_size()` しか走査せず**プロセスローカルの最大値**を返すので、
  `allreduce_max` で全プロセスにまたがって取る(既存の同型パターンが
  `src/iTPS/load_toml.cpp:583` にある)。
- 閾値は `Tn[i]` の最大絶対値に対する相対値で **1e-12 × max_abs**(絶対ゼロ比較は丸めで落ちる)。
  `max_abs` も `allreduce_max` で取る。
- V7 は §2.4 のスパイク(3/1 の保存を 2/2 のラベルで読む)を実際に捕まえる網であり、
  V1〜V6 をすり抜けた不整合に対する最後の防壁でもある。

MPI: `fermion.dat` の読み出しは `mpirank == 0` で行い、`std::vector<int>` に平坦化して
`bcast`(`src/mpi.hpp:101`)する。`std::vector<bool>` はビット詰めなので直接 bcast しない。

## 5. 実装設計

### 5.1 変更ファイル

| ファイル | 変更 |
|---|---|
| `src/iTPS/saveload_tensors.cpp` | `save_fermion_parity()` / `load_fermion_parity()` を追加し、`save_tensors()` / `load_tensors()` から呼ぶ |
| `src/iTPS/iTPS.hpp` | 上記 2 メソッドの宣言(private) |
| `src/iTPS/load_toml.cpp` | `tensor load/save directories` ガードの削除 |
| `docs/sphinx/{ja,en}/file_specification/parameter_section.rst` | 非対応一覧から除外、`tensor_save`/`tensor_load` に D 固定の注記 |
| `NEWS.md` | 項目追加 |

### 5.2 `save_tensors()` への追加

末尾で、`finfo.enabled` なら `save_fermion_parity(save_dir)` を呼ぶ(中で `mpirank == 0` を判定)。
`finfo.phys[i]`、`finfo.virt[i][leg]`、`lattice.LX` / `lattice.LY` / `lattice.skew` を §3.1 の
形式で書く。あわせて「環境テンソルは fermion 模式では意味のある値を持たない(§2.3)」旨と、
`save_tensors()` が `measure()` の前に呼ばれる前提に依存していることをコメントで残す。

### 5.3 `load_tensors()` の分割

```cpp
template <class ptensor>
void iTPS<ptensor>::load_tensors() {
  std::string const &load_dir = peps_parameters.tensor_load_dir;
  // ... 既存の isdir 検査と format version の読み出し ...

  load_fermion_ledger(load_dir);   // 前段: V1-V6, V8。finfo.virt を復元する

  if (tensor_format_version == 0) {
    load_tensors_v0();
  } else if (tensor_format_version == 1) {
    load_tensors_v1();
  } else { ... }

  validate_loaded_fermion_tensors();  // 後段: V7(finfo.enabled のときだけ)
}
```

- `load_fermion_ledger` は `finfo.enabled` が false でも呼ぶ(V8 の判定のため)。
- 前段が `finfo.virt` を書き換えるので、`load_tensors_v1` 内のテンソル読み込みは
  復元後の台帳のもとで行われる。`load_tensors_v1` 自体は変更しない。
- `validate_loaded_fermion_tensors` は `finfo.enabled` のときだけ中身を実行する。

### 5.4 ボソン経路への影響

`fermion.dat` を書かない・読まない(V8 の存在確認のみ)。既存の ctest(`restart` を含む)の
数値は不変。

## 6. ドキュメント

- `parameter_section.rst`(ja/en): フェルミオンの非対応一覧から `tensor_save` / `tensor_load` を
  外し、代わりに「フェルミオン模式では保存時と同じ `virtual_dim` でのみ読み込める。
  仮想脚の偶奇台帳が `fermion.dat` に保存され、読み込み時に検証される」旨を書く。
- `NEWS.md`: 「フェルミオン模式で `tensor_save` / `tensor_load` が使えるようになった
  (パラメータスキャンの再開用)。仮想ボンドのパリティ台帳が一緒に保存され、
  読み込み時に物理脚パリティ・次元・隣接整合・テンソルのパリティ違反が検証される。
  `virtual_dim` を変えての再開は非対応」。

## 7. 検証戦略

Fock オラクル(`test/fermion/fock_oracle.py`)は不要。ここで守るべき性質は
「保存した状態を読み込むと同じ状態になる」という往復の同一性と、
「壊れた組合せは必ずエラーになる」というガードなので、参照は**同一実行内の値**と
**例外**である。

### 層1: 台帳の往復(`test/test_fermion_layer.cpp` に doctest)

`iTPSTestAccessor` 経由で `finfo.virt` を**既定と異なる分割**(例: D=2 で `{false, false}`、
D=4 で `{f,f,f,t}`)に差し替え、`save_tensors()` → 別インスタンスで `load_tensors()` →
`finfo.virt` が要素まで一致すること。既定と同じ分割だけで試すと、台帳を読まない実装でも
通ってしまう(空洞)ので、**必ず既定と違う分割**を使う。

### 層2: 測定値の往復

層1 と同じ状態で `measure_twosite()` / `measure_onesite()` の値が保存前後で一致すること。
台帳を復元しない実装では、既定と違う分割にしてあるので値が変わる(＝これが RED になる)。

### 層3: ガード(V1〜V8)

各条件について `CHECK_THROWS_AS(..., tenes::load_error)`(または適切な例外型)。
特に V7 は「ボソン実行で作った保存ディレクトリに手で `fermion.dat` を足す」ではなく、
**パリティ違反を持つ `Tn` を保存して読み込む**形で作る(V1〜V6 をすり抜ける経路を通す)。
V5(D 不一致)は「エラーになること」に加えて、**エラー後に `Tn` / `lambda` が
resize されていないこと**も見る(§4 前段の存在理由)。

**注意**: §2.4 のとおり、既定の初期化では分割が `even_first_parity` のままなので、
**E2E の往復(層4)は台帳を復元しない実装でも通ってしまう**。台帳そのものを守るのは
層1・層2 であり、そこでは必ず既定と異なる分割を作ること。

### 層4: E2E(`ctest`)

`test/fermion/free_fermion_saveload.py.in`(新規)または既存 `restart.py.in` に倣う:

1. spinless free fermion(D=2、`num_step` 少なめ)を `tensor_save` 付きで実行
2. 同じ入力から `num_step = 0`、`tensor_load` で再実行
3. `output/density.dat` が 1 と 2 で `np.isclose(rtol=1e-8)` で一致すること
   (`num_step = 0` なので状態は変わらず、測定だけが走る。スパイクではビット一致した)
4. さらに μ を変えた 3 回目を `tensor_load` から流し、正常終了して
   エネルギーが有限であること(スキャン運用の煙テスト)
5. 保存ディレクトリに `fermion.dat` が存在し、行数・各行の長さが `N_UNIT`・
   `virtual_dim`・`physical_dim` と整合すること(パースできること)

`fermion.dat` の中身が `even_first_parity` と違うことを E2E で要求してはいけない。
§2.4 のとおり既定の初期化では一致するのが正常だからである。

E2E が RED になるのは 5(`fermion.dat` の存在)と、実装前のガードによる異常終了である。
3 の一致自体は §2.4 のとおり台帳なしでも成り立つので、E2E を台帳の網と見なさないこと。

### 変異テスト(レビュアー向け)

- `load_fermion_parity` の `finfo.virt` 代入を消す → 層1・層2・E2E が赤
- V7 のパリティ違反検査を消す → 層3 の該当ケースが赤
- 保存を偶の個数だけにして `even_first_parity` で再構成する → 層1(非偶先頭の分割)が赤

## 8. タスク分割(実装計画の骨子)

| # | 内容 | 担当 |
|---|---|---|
| T1 | 層1〜3 のテスト → `save_fermion_parity` / `load_fermion_ledger` / `validate_loaded_fermion_tensors` 実装、ガード削除 | テスト作成者 → Codex |
| T2 | 層4 E2E、docs(ja/en)、NEWS | テスト作成者 → Codex |

T1 の RED 確認では、実装前に `tensor_save`/`tensor_load` がガードで拒否されるため、
テストは `validate_fermion_constraints` を経由しない `iTPS` 直接構築で書く(層1〜3)。
E2E(層4)はガードのエラーで RED になる。

## 9. リスクと判断

- **順序依存**(§2.3): `save_tensors()` が `measure()` の前に呼ばれる前提。E2E の往復が
  この前提を固定する。将来 measure 後に保存するようになれば E2E が落ちる。
- **格子の同一性検査は `L_sub` / `skew` まで**。同じ `L_sub` で `[[tensor.unitcell]]` の
  `index` 割り当てだけを組み替えた入力は検出できない(ボソン経路も同じ制約)。
- **D 固定の制約**: ユーザーの用途(パラメータスキャン)には十分。D を上げたい場合は
  偶奇を保った埋め込みが別途必要で、本設計の範囲外(§3.3)。
- **CTM のウォームスタートが効かない**: 読み込み後の最初の測定でも CTM は最初から
  作り直される。これは既存の挙動であり本設計では変えない。フェルミオンの CTM 構築は
  D=4 で 10 秒程度(実測)なのでスキャン運用の障害にはならない。
- **平均場環境との併用**: `meanfield_env = true` なら CTM 自体が不要なので、
  スキャンは MF で回し、確認したい点だけ CTM で測る運用ができる。
