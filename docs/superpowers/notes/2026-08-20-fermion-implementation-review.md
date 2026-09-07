# 2026-08-20: Fermion implementation review

2026-08-19-fermion-implementation-review.md に基づきレビューを行った。

## 0

## 1

## 2

## 3

parity はphysical bond に紐づくfermion operator のparity から来るので、virtual bond のparity は初期値がevenということでOK?
even / odd を表現する値については明示的に示しておきたい。下記でよろしいか。

| even  | odd  |
|-------|------|
| 1     | -1   |
| 0     | 1    |
| false | true |

### 3.1

`transpose` では軸順序の入れ替えに伴い要素に符号マスクをかける（要素の符号を反転する）。その後、 `transpose` する。

`transpose_sign` はこのsignを計算する関数。odd 要素が何個の odd要素を追い越したかを数え、偶奇に応じて1 or -1 を返す。

- 符号マスクについて
  - `transpose_sign` の結果でif 分岐しているが、1 or -1 を返すので、直接この値を要素にかければ良い
  - `get_value` と `set_value` を使っているが、これはglobal index から local index への変換を伴う。 `operator[]` を使えばlocal indexから直接要素にアクセスできるはず。
  - まとめると、 `ret.t[n] *= transpose_sign(ret.parity, ret.t.global_index(n), axes);` と書ける。
  - `transpose` 以外の符号マスクでも同様。

### 3.2

`tensordot_left_perm(rank, axes)` は、 `axes` を最後に持ってくるように並べ替えた `Axes` を返す。
例： `tensordot_left_perm(5, {0,2}) -> {1,3,4,0,2}`
`tensordot_right_perm(rank, axes)` は逆に、 `axes` を逆順にして最初に持ってくる。
例： `tensordot_right_perm(5, {0,2}) -> {2,0,1,3,4}`

ftensor をこのpermにしたがってtransposeしたあとにtensordotをすると、tensordotではparityを考慮しなくても良い。これが意味論。
`transpose` と同様に符号マスクを計算・適用することで、 `transpose` を行わずに `tensordot` を実行できる。これが実装。
`apply_transpose_sign_mask` は符号マスクを計算・適用したftensorを返す関数。

### 3.3

複素共役+軸順序反転。つまりエルミート共役？
「JW 厳密参照テストが唯一の確定根拠」というのは、 `conj` がエルミート共役であることがどこにも明示されていないということ？

符号の計算は、 odd 要素を2つ選ぶ組み合わせの数になるため、この式で問題ない。
local/global indexについて、符号マスクと同様の注意が成り立つ。

### 3.4

`qr` などのテンソル分解では、テンソルの脚を並べ替え・束ねて行列を作る。
まずは row, col を用いてtranspose。
次にreshapeで行列へと変換するが、parity情報と要素 `t` とを別個に管理する。
`t` はbosonic reshape で単純に変換 (`mat`)。

row, col それぞれのparity情報は `fuse_axes` で融合して、even / odd の値を計算する。

`fuse_axes(parities, axes)` は `axes` に対応するlegのparityを `fuse` して新しい `parity_vector` を作る。

`fuse(a,b)` は2つの `parity_vector` を融合して（束ねて）新しい `parity_vector r` を作る。
r の要素数は a と b の要素数の積であり、各要素は a と b の要素の排他的論理和である(parityとしての和)。

ここで偶奇のセクターに分けるため、even 要素を前に、 odd 要素を後に並べ替える。

`parity_sort_perm(p)` は、 even 要素を前に、 odd 要素を後に並べ替えるためのindex list `perm` を返す。
`make_perm_matrix(perm)` は、 `perm` に従ってindexを並べ替えるための行列を返す。

`row`, `col` それぞれで変換行列 `prow, pcol` を作り、 `mat` と tensordot することで、 `sorted` をつくる。
`row`, `col` 中の even要素数を用いて `slice` することで even 要素のみの行列と odd 要素のみの行列に分け、それぞれで `qr` を行う。
分解結果は全体のQやR行列 `q_sorted, r_sorted` の even / odd ブロックとして保存する。
`q_sorted, r_sorted` について、外側の軸を `prow` と `pcol` を用いて逆変換することで、全体のQやR行列 `q_mat, r_mat` を得る。
QとRをつなぐ脚は、 evenが前、oddが後ろにあるように並んだまま。
それぞれのparity情報を適切に並び替える。

`svd` など、他のテンソル分解でも同様。

questions:

- `fuse_axes` でaxes が空の場合には `{false}` を返しているが、 `{}` ではない理由はあるか？
- even / odd セクターが混じることはない？
- interleaved_split だと失敗する理由がよくわからない。実際にsimple update中では、 `interleaved_split` を避けるようにあらかじめ軸を並び替えたらうまくいくようだが、何が変わるのか？

## 4

### 4.1

`log_theta_blocks` はデバッグ用のログ出力？

### 4.2

(a) の正規化について、 `source` のサイト番号が `target` のサイト番号よりも小さくしたいということと思われる。
ボンド方位のみを見る場合、unitcell をまたぐ時にサイト番号が逆転することがあるのでは？素直にサイト番号を比べれば良さそう？
(b) がよくわからない。 `op12` で保持している演算子（の行列表現）は、すでに交換関係が考慮されているのでは？

## 5

### 5.2

`doubled_pipeline` の `braTn` と `ketT` で、 `T` と `Tn` の違いはないはず？どちらかに命名を統一。

ここではswap gate を導入することで reduced tensor (=one-site density matrix) に fermion signを反映させている？

## 8

1. reduced tensor の作り方が完全には理解できていない。 5.2 へのコメントも参照。
2. simple updateでは入力のみで、2サイト測定blobの作成時には入出力でswapが必要というのは、前者は出力legが開いたままだけれど、後者は閉じているからだろうか？ `wrap_twosite_op` という名前よりも `wrap_twosite_evolution` という名前のほうが適切かも知れない。
3. 特異値が縮退した時に even / odd のどちらを優先するかということだろうか。パラメータで指定できるようにする？
4. よくわかっていない。
5. 性能についてはひとまずそのままでよい。
6. 特に問題は見つかっていない。
