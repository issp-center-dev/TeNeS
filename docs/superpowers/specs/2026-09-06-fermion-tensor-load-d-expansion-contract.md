# フェルミオン系で tensor_load による D 拡張を可能にする — 振る舞い契約(2026-09-06)

前提の実測: `work/fermion/dexpand/FINDINGS.md`(使い捨てパッチ、再現入力つき)

## 足すこと

fermion モードでは、保存したテンソルを **`virtual_dim` を変えて** 読み込むことができない。
小さい `D` で得た状態を大きい `D` の計算の初期値にする、という標準的な使い方ができない。

**拡張(`D_new > D_old`)だけを可能にする。** 縮小は従来どおり拒否する。

**対象は fermion の基底状態計算経路に限る。** `validate_fermion_constraints` は
`calcmode != ground_state` を拒否する(`src/iTPS/load_toml.cpp:634`)ので、
有限温度 / 実時間発展 は fermion モード自体が正規入口で弾かれる。
`initialize_tensors_density()`(`src/iTPS/tensors.cpp:185`)はそもそも `finfo` を構築せず、
`load_fermion_ledger` は `!finfo.enabled` の枝に入る。density 経路の D 拡張は今回の対象外。

## 前提として確認済みの事実

### F1: 障壁は台帳の長さチェック1箇所だけ

拒否しているのは `src/iTPS/saveload_tensors.cpp` の `load_fermion_ledger` で、保存された
偶奇台帳の長さが入力の `virtual_dim` と一致しない場合に投げる箇所(2026-09-06 時点で 332-342 行)。
テンソル本体・`lambda` は**既に**新しい形へ合わせる仕組みが動いている:

- `load_tensor` が `resize_tensor` を通す(`saveload_tensors.cpp:405`)。
  `resize_tensor` は伸びる軸を `mptensor::extend` でゼロ詰めする(`src/tensor.cpp:115`)。
  **新しいインデックスは各脚の末尾に付く。**
- `lambda` は保存時の長さ(`params.dat` の shape)だけ読んでから新しい `virtual_dim` へ
  `resize` される(`saveload_tensors.cpp:533-538`)。伸びた分は 0。
- 0 になった `lambda` の逆数は `Inverse_lambda_cut` が 0 に潰す
  (`src/iTPS/core/simple_update.cpp:155,165`)。
- `params.dat` の shape 不一致は WARNING のみでエラーにしない(`saveload_tensors.cpp:455`)。
- fermion モードの CTM 環境はもともとプレースホルダとして保存され毎回作り直される
  (`saveload_tensors.cpp:66-79`)。D が変わっても失われる情報はない。

### F2: 読み込み時の台帳は数スイープだけの初期値である

simple update はボンドを更新するたびに graded `svd_trunc` の結果で台帳を上書きする
(`src/iTPS/simple_update.cpp:128-129`、`src/fermion/fops.hpp` の `svd_trunc` が
内部脚の `internal` を生成)。**したがって読み込み時のパディング規則が結果を支配することはない。**

実測(`FINDINGS.md`「パディング規則の感度」)でも、バランス型と全偶パディングの最終エネルギーは
-0.791099251 と -0.791099382 で 1.3e-7 しか違わなかった。

### F3: ゼロ詰めのままで拡張は効く。ノイズ注入は不要

自由フェルミオン(スピンレス、正方格子、t=1、mu=0、2x2、chi=8、SU 1000 ステップ、CTM 環境、
seed 11)で、D=2 の -0.723448 が D=4 引き継ぎで **-0.791099** になる(厳密値 -0.810564)。
拡張分の振幅は厳密にゼロだが、ボンドのランクは自分自身ではなく側面 3 脚が張る空間から来るため、
最初のスイープで埋まる。**ノイズ注入は不要。**

**引き継ぎがスクラッチより良いとは言えない**(初版の契約書はそう書いていたが誤り)。
D=4 スクラッチはシードで大きく散り、7 シード中 4 回は引き継ぎと同じ ≈ -0.7911 に達し、
3 回は -0.7415 / -0.7550 / -0.7578 の悪い固定点に落ちる。seed 11 の -0.741466 は外れ側で、
ステップ数を 8000 まで伸ばしても -0.743592 で頭打ちなので収束不足ではない。
一方 **引き継ぎは 5 シードすべて -0.791098 〜 -0.791104(ばらつき 5.5e-6)** に収まる
(D=2 自体がシードにほぼ依らないため)。数値は `FINDINGS.md` の「結果 3」。

したがって機能の価値は「エネルギーが良い」ことではなく **「良い固定点を再現的に与える」** こと。
**この性質はシードの統計に依るのでテストの判定条件にはできない**(S-2 参照)。

### F4: 拡張分のパリティは何を割り当ててもパリティ保存を破らない

拡張されたインデックスの振幅は厳密に 0 なので、`validate_loaded_fermion_tensors` の
パリティ違反検査は割り当てによらず通る。実際に効く制約は2つだけ:

- 既存インデックスの偶奇を変えないこと(変えると既存の非ゼロ要素がパリティ違反になる)
- ボンドの両端で台帳が一致すること(`validate_neighbor_consistency`)

新しいインデックスは末尾に付くので、**台帳も末尾に追加すれば前者は自動的に満たされる**。
後者は、追加規則が決定的で、両端の既存台帳が(保存時点で)一致している以上、
入力の `virtual_dim` がボンドの両端で同じであれば自動的に満たされる。両端で食い違う入力は
`validate_neighbor_consistency` が従来どおり捕まえる(その検査は今回変更しない)。

## 要求

### R1: `extend_parity` を追加する

`src/fermion/fermion_info.hpp` に、`even_first_parity` と並べて

```cpp
parity_vector extend_parity(const parity_vector& p, std::size_t new_dim);
```

を追加する。`parity_vector` の要素は **false = 偶、true = 奇**(以下 0 = 偶、1 = 奇と書く)。

- `new_dim < p.size()` のとき `std::runtime_error` を投げる。
- `new_dim == p.size()` のとき `p` をそのまま返す。
- `new_dim > p.size()` のとき、**`p` の全要素を順序どおり保ったまま**末尾に
  `new_dim - p.size()` 個を追加して返す。

追加規則は、1 個足すごとに「**足した後の**長さ `n` に対する偶の目標数 `ceil(n/2)` に、
足す前の偶の数が足りていなければ偶、足りていれば奇」を足す。
`n` は追加後の長さである(追加前の長さと読むと下の表と矛盾する)。
**期待値(テストはこの値を直接ピンする)**:

| 入力 | `new_dim` | 出力 |
|---|---|---|
| `{0,1}` | 4 | `{0,1,0,1}` |
| `{0}` | 4 | `{0,1,0,1}` |
| `{0,0,0}` | 4 | `{0,0,0,1}` |
| `{1,1}` | 4 | `{1,1,0,0}` |
| `{0,1}` | 2 | `{0,1}` |
| `{0,1}` | 5 | `{0,1,0,1,0}` |

この性質は `extend_parity` の Doxygen コメントにも 1 行書くこと
(`even_first_parity` の「偶が先、奇が後」という並びは保たれない。例:
`extend_parity({0,0,1,1}, 6)` = `{0,0,1,1,0,1}`。現状これに依存する箇所は無いが、
将来 even-first を仮定するコードが書かれると壊れる。2026-09-07 追記)。

**注意**: `extend_parity(p, n)` は一般に `even_first_parity(n)` **とは一致しない**。
並び順が違う(`extend_parity({0,1}, 4)` = `{0,1,0,1}` に対し `even_first_parity(4)` = `{0,0,1,1}`)
だけでなく、既存部分が既に片方に偏っていれば **偶奇の個数も一致しない**
(`extend_parity({0,0,0}, 4)` = `{0,0,0,1}` は偶 3・奇 1、`even_first_parity(4)` は偶 2・奇 2)。
既存要素は動かせないので、目標配分には**足せる範囲で**近づくだけである。
F2 のとおり、どちらの差も結果を左右しない。

### R2: `load_fermion_ledger` が拡張を受理する

`load_fermion_ledger` は、`fermion.dat` から読んだ脚の台帳 `p` について

- `p.size() < lattice.virtual_dims[i][leg]` のとき、`extend_parity` で伸ばして採用する
- `p.size() > lattice.virtual_dims[i][leg]` のとき、従来どおり `tenes::load_error` を投げる。
  メッセージの HINT を、**縮小のみが不可**であることが分かる文面に改める
  (現行は「`virtual_dim` を変えられない」と読める)
- 一致するときの挙動は**1 バイトも変えない**

`fermion.dat` の**保存**形式は変えない(バージョンも据え置き)。読む側だけが緩む。

**台帳の長さは「保存されたテンソルの次元」と突き合わせること(2026-09-06 追記)。**

初版の R2 は「台帳の長さを入力の `virtual_dim` と比べて、短ければ拡張する」としか書いておらず、
**台帳が保存テンソルと食い違っている壊れたチェックポイントを検出する道を塞いでいた**。
既存テスト `SL V5 a refused load leaves the tensors untouched`(`test/fermion/saveload.cpp:572`)は、
D=4 で保存したチェックポイントの `fermion.dat` を長さ 2 に壊し、`virtual_dim = 4` で読ませて
**テンソルが読まれる前に**拒否されることを検査している。初版どおりに実装すると 2 < 4 が
「拡張」として受理され、後段の `validate_loaded_fermion_tensors` まで落ちない。

正しい不変条件は次の 2 つである。

1. **`fermion.dat` の各仮想脚の台帳の長さは、`params.dat` の `Shape of Tn[i]` 行が言う
   保存時の脚次元と一致しなければならない。** 一致しなければ壊れたチェックポイントであり
   `tenes::load_error` を投げる。
2. 入力の `virtual_dim` は保存時の脚次元**以上**でなければならない。小さければ従来どおり
   `tenes::load_error`(凍結された ERROR 行の文面)。

拡張は「保存時の脚次元 → 入力の `virtual_dim`」で行う。

**`tensor_format_version = 0` の経路では検査 1 を行わないこと(2026-09-07 追記)。**
v0 には `params.dat` が無いので保存時の脚次元が分からない。現状は lattice の現在値を
`saved_shape` として渡しているため、別の `D` で読むと
`has 2 entries but the saved tensor has dimension 4`(4 は入力値であって保存テンソルの次元ではない)
と `fermion.dat and params.dat do not describe the same saved tensors`(v0 に `params.dat` は無い)
という **二重に誤ったメッセージ**が出る。v0 では変更前の挙動(台帳長を入力の `virtual_dim` と
比べ、違えば凍結された ERROR 行で拒否。拡張もしない)に戻すこと。
契約の「拡張は正規の v1 チェックポイントに限る」と一致する。

**この 2 つの検査は、テンソルファイルを 1 つも読む前に完了していなければならない。**
`load_tensors_v1` は `load_tensor` のループに入る前に `params.dat` の shape を読んで
bcast 済みなので、そこに検査と拡張を置けば追加のパース無しで条件を満たせる
(構造の選択は実装者に委ねるが、この位置なら成立することは確認済み)。

**`load_tensors()` に巻き戻し(状態のフルコピーと try/catch による復元)を入れてはならない。**
既存の契約は「**前段**の検証が失敗したときテンソルは無傷」であって、後段の失敗
(パリティ違反など)では従来から上書きされる。巻き戻しはその保証を勝手に広げるうえ、
`load_tensors()` は boson と共通なので、`Tn` と CTM 環境 8 本と `lambda_tensor` の
フルコピーが全ユーザーの読み込み時ピークメモリを倍にする
(CHI=200・D=8・2x2 セルで環境だけ 328MB)。上の 1. を守れば巻き戻しは不要になる。

**ERROR 行の文言は変えないこと。** 縮小を拒否するときの
`ERROR: the virtual dimension of the leg <n> of the tensor <i> is <A> but the saved tensors have <B>.`
という文面に、既存の `test/fermion/free_fermion_saveload.py.in`(run4)と新しいユニットテスト
(DE layer3)が依存している。改めてよいのはその後ろの HINT だけである
(HINT 中の `virtual_dim` という語も、run4 が以前それに依存していた経緯があるので、
文言を変える場合はテストが ERROR 行側に張り替わっていることを前提にする)。

### R3: 防御的な範囲修正

`src/iTPS/core/simple_update.cpp` の正規化ループ(2026-09-06 時点で 238-245 行、`norm` の 2 つのループ)は
`lambda_c` を `for (int i = 0; i < dc; ++i)` で舐めるが、`lambda_c` の長さは
`svd_trunc` が返す `nkeep = min(dc, full_s.size())`(`src/fermion/fops.hpp`)であって
`dc` とは限らない。ループ上限を `lambda_c.size()` に改める。

**挙動は変えない**(`nkeep == dc` のときは同じ計算)。実測では全偶パディングという偏った
条件でも `nkeep < dc` は発現しなかったので、これは未定義動作を残さないための処置であり、
このバグを踏むテストは要求しない。

**注意**: `core::Simple_update_bond` は boson と fermion で共通のテンプレートである。
`nkeep < dc` を返し得ると分かっているのは fermion の `svd_trunc` だけだが、
修正は共通カーネルのループ上限を `lambda_c.size()` に変えるものであり、boson の
インスタンス化にも掛かる。`nkeep == dc` である限りどちらも計算は同一。

### R5: `params.dat` の検証を全ランクで行う(2026-09-07 追加)

`load_tensors_v1` の `params.dat` 解析は `if (mpirank == 0)` の中で行われ、
`format_version` / `N_UNIT` / 物理脚次元の不一致を **ランク 0 だけが** 投げる。他のランクは
直後の `bcast(loaded_shape[i], 0, comm)` で待つのでハングする。

この構図自体は既存(boson は以前からそう)だが、**今回 `load_fermion_ledger` の呼び出しを
`load_tensors_v1` の内部に移したことで、fermion モードでも露出するようになった**。
移動前は `load_fermion_ledger` が先に走り、`fermion.dat` の内容を全ランクに bcast してから
`N_UNIT` / `L_sub` / `skew` / 物理パリティを全ランクで検査していたので、食い違うチェックポイントは
クリーンなエラーで止まっていた。**クリーンなエラーからハングへの後退であり、直すこと。**

`params.dat` から読んだ値を bcast したうえで、検証を全ランクで行うようにする。
`loaded_shape` は既に bcast されているので、`format_version` と `N_UNIT` を同様に配り、
throw する箇所を `if (mpirank == 0)` の外へ出せばよい。

ctest には多ランク実行が 1 件も無い(すべて `-n 1`)ので、**この後退はテストでは捕まらない**。
テストの追加は求めない。

### R4: ドキュメント

`docs/sphinx/ja/file_specification/parameter_section.rst:69` の
「``virtual_dim`` を変えての読み込みは非対応です(偶奇ブロックの構造が保てないため)」と、
`docs/sphinx/en/file_specification/parameter_section.rst:68` の
"Restarting with a different ``virtual_dim`` is not supported, because resizing a leg does not
preserve its even/odd blocks" を、**拡張は可、縮小は不可**に改める。

`NEWS.md` の 11 行目にも同じ主張
("restarting with a different ``virtual_dim`` is not supported")があり、
R4 の初版はこれを見落としていた(2026-09-07 追記)。あわせて直すこと。
拡張分がゼロ詰めであること、したがって効果を出すには読み込み後に simple update を
回す必要があること(測定だけでは広がらない)を1文添える。

## テストへの要求(テスト作成者が書く)

### U: ユニット(C++ doctest)

`test/fermion/saveload.cpp` の層構成に合わせて足すか、`test/test_fermion_layer.cpp` 側に足すかは
テスト作成者の判断でよい。

- **(U-1)** `extend_parity` が R1 の表の 6 例をそのまま返すこと。**表の値を直接書くこと**
  (実装を呼んで比べるのではなく、期待値をリテラルで置く)。
- **(U-2)** `extend_parity` が `new_dim < p.size()` で投げること。
- **(U-3)** `extend_parity` が既存要素を破壊しないこと。長さ 1..6、全パターンの `p` を
  総当たりし、`new_dim` を `p.size()`..8 まで振って、
  「先頭 `p.size()` 要素が `p` と一致」「長さが `new_dim`」を検査する。
- **(U-4)** D_old = 2 の `fermion.dat` を D_new = 4 の入力で読み込むと `load_fermion_ledger`
  が成功し、`finfo.virt` の各脚が長さ 4 になり、`validate_neighbor_consistency` を通ること。
  **反空洞化**: 保存する台帳は `even_first_parity(2)` = `{0,1}` 以外にもう1通り
  (例: `{1,0}`)を試し、既存部分がその通り保たれることを見ること。
- **(U-5)** D_old = 4 の `fermion.dat` を D_new = 2 の入力で読み込むと `tenes::load_error`
  で落ちること(縮小の拒否が残っていること)。

### S: E2E(Python、`test/fermion/` に新規 `.py.in` + `test/CMakeLists.txt` に登録)

`sample/07_spinless_fermion/input.toml` と同じ模型(スピンレス自由フェルミオン、正方格子、
t = 1、mu = 0、`L_sub = [2,2]`、`seed = 11`、`[parameter.simple_update] tau = 0.01`、
`num_step = 1000`)で、**両 run とも `[parameter.ctm] dimension = 8`**、
`meanfield_env` は既定(false、CTM 環境)。

- **run A**: `virtual_dim = [2,2,2,2]`、`tensor_save`
- **run B**: `virtual_dim = [4,4,4,4]`、run A の保存先を `tensor_load`
  (run B の入力は run A と `virtual_dim` / `output` / `tensor_save` / `tensor_load` 以外
  1 行も違わないこと。tau・num_step・chi・seed を取り違えたまま比較する事故を塞ぐ)
**(S-1)** run B が正常終了すること(修正前は exit 1 で
`ERROR: the virtual dimension of the leg 0 of the tensor 0 is 4 but the saved tensors have 2.`
で落ちる。これが RED)。

**(S-2)** `E_B < E_A - 0.03` かつ `|E_B - (-0.791099)| < 1.0e-3` かつ `|n_B - 0.5| < 1.0e-6`。

実測は `E_A = -0.723448`、`E_B = -0.791099` で差は 0.0677、しきい値はその半分未満。
`E_B` のピン留めは、引き継ぎが 5 シードで -0.791098 〜 -0.791104 に収まる(ばらつき 5.5e-6)
という実測に基づく。`n_B` の実測は 0.500000000000000。

**D=4 のスクラッチ run を対照として比較してはならない。** D=4 スクラッチはシードで
-0.7415 〜 -0.7913 に散る(F3)ので、`E_B < E_scratch` のような判定はシード次第で
偽陽性にも偽陰性にもなる。E2E が検査するのは「拡張が実際に効いて D=2 より良くなること」だけで、
「引き継ぎがスクラッチより優れること」ではない。後者は成り立たない。

**(S-3)** **反空洞化**: run A が本当に D=2 で走ったことを、run A の保存先の `params.dat` の
shape 行が `2 2 2 2 2` であることで確かめること。これを見ないと、run A が何かの拍子に D=4 で
走った場合に S-2 が意味を失う。
run B 側も `tensor_save` を付け、その `params.dat` の shape 行が `4 4 4 4 2`、
`fermion.dat` の各仮想脚の台帳が長さ 4 であることを確かめること。

**(S-4)** `chi = 8 < D*D = 16` なので run B は
`WARNING: CTM may be too small (chi < D*D) for iTPS` を出す。**これは想定どおり**であり、
テストはこの警告で落ちてはならない。

**このテストが検出できないこと(テスト作成者への注意)**: F2 のとおり台帳は最初のボンド更新で
書き直されるので、**E2E は `extend_parity` の規則が間違っていても緑になる**。
パディング規則の正しさは U-1/U-3/U-4 だけが担保する。E2E をそのための検査に使わないこと。

### 既存テスト

**(S-5)** boson の save/load、および fermion で D を変えない save/load の挙動が
変わらないこと。

**2026-09-06 訂正**: 初版は「既存テストが緑のままであることで足りる」と書いていたが、**これは誤り**。
`test/fermion/free_fermion_saveload.py.in` の run4 は、`D = 2` で保存したチェックポイントを
`virtual_dim = 3` で読み込んで **拒否されること** を検査していた。3 > 2 なので R2 の後はこれが
受理され、既存テストが赤くなる。テスト作成者がこれを検出し、run4 を `D - 1`(= 1、縮小)に
差し替えたうえで、判定の手掛かりを HINT 中の `virtual_dim` から ERROR 行の
`virtual dimension` に張り替えた。**この差し替えは意図的であり、実装者が元に戻してはならない。**
差し替え後の `FreeFermionSaveLoad` が実装前の時点で緑であることは確認済み。

### レビュー(タスク 8)で見つかった追加要求(2026-09-07)

**(S-6)** `test/fermion/free_fermion_saveload.py.in` の run4 の判定文字列を張り替えること。
現在の

```python
    if "virtual dimension" not in message4:
```

は **空洞** である。`virtual dimension` という語は凍結された ERROR 行だけでなく、
D を変えると必ず出る WARNING(`src/iTPS/saveload_tensors.cpp:478`
`WARNING: virtual dimension of the leg ... is ... but loaded tensor has ...`)にも含まれる。
縮小の実測では WARNING が 16 行出たうえで ERROR が 1 行出るので、**ERROR が出なくても
この判定は通る**。実際、変異テスト(「要求 virtual_dim >= 保存 shape」の検査を削除)で
`tenes` が別の例外で死んでも `FreeFermionSaveLoad` は緑のままだった。

ERROR 行の接頭辞まで含めて突き合わせること:

```python
    if "ERROR: the virtual dimension of the leg" not in message4:
```

**(U-6)** 「台帳が保存 shape より**長い**」壊れたチェックポイントを拒否することのテストが無い。
R2 の追記節の検査 1 は今回新設したガードなので、テストを 1 件足すこと。
`fermion.dat` の台帳だけを長くした checkpoint が `tenes::load_error` になり、
そのメッセージが検査 1 のもの(`the virtual parity ledger of the leg`)であること、
かつ **テンソルが読まれる前に**落ちること(V5b と同じマーカー方式)を見ること。

**(U-7)** DE layer2 が `lambda_tensor` を見ていない。拡張後に各脚の `lambda` の長さが
新しい `virtual_dim` になり、保存されていた先頭部分が保たれ、増えた分が 0 であることを
足すこと。

## 完了条件

1. 修正前に S-1 と U-1 が **RED** であることを確認していること(S-1 は上記の
   `ERROR: the virtual dimension ...` で、U-1 は `extend_parity` 未定義で)。
2. `ctest --preset gcc-release` が全件緑。
3. `work/fermion/dexpand/` の 3 入力(`in_A_d2.toml` / `in_C_d4.toml` /
   `in_B_d2to4.toml`)の手動再現が FINDINGS.md の表と一致すること。
4. 報告 `work/fermion/dexpand/report-*.md`。

## 範囲外(残課題)

- **縮小(`D_new < D_old`)**: `resize_tensor` は先頭 `D_new` 個を取るが、graded の並びは
  「偶を大きい順、次に奇を大きい順」であって特異値の大域順ではない。素朴に切ると奇セクタを
  丸ごと落としかねない。物理的に正しい切り方の設計が別途要る。
- **`CHI` 変更**: 今回は触らない(fermion の CTM 環境はどのみち作り直される)。
- **WARNING の冗長さ**: D を変えると `N_UNIT × 4` 行の WARNING が出る(2x2 で 16 行)。
  情報としては正しいので今回は放置する。
- **`tensor_format_version = 0` のチェックポイント**: 正規に保存された v0 は `fermion.dat` を
  持たないので、fermion モードでは既に明確なエラーで停止する。ただし `load_tensors_v0` は
  `load_tensor` / `resize_tensor` を通さずテンソルを保存時の shape のまま `A.load()` する
  (`saveload_tensors.cpp:563` 以降)。`params.dat` を欠いたまま `fermion.dat` だけ足した
  手製データを食わせると、台帳だけが拡張され Tn と lambda が古い shape のまま残り得る。
  **これは非対応**とし、今回の拡張は正規の v1 チェックポイントに限る。
