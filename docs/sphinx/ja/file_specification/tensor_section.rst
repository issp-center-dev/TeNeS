.. highlight:: none

「ユニットセル」の情報を記述します。
ユニットセルは ``Lx`` かける ``Ly`` の大きさをもつ長方形の形をしています。
また、サブセクション ``unitcell`` を持ちます。

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 15, 30, 20, 10

   ``L_sub``, "ユニットセルの大きさ", 整数または整数のリスト, "--"
   ``skew``, "skew 境界条件におけるシフト値", 整数, 0


``L_sub`` として2つの整数からなるリストを渡した場合、はじめの要素が ``Lx`` に、もう片方が ``Ly`` になります。
3つ以上の要素からなるリストを渡した場合にはエラー終了します。
``L_sub`` として整数を渡した場合、 ``Lx`` と ``Ly`` とが等しくなります。

ユニットセル内のサイトは0から順番に番号付けされます。 x 方向から順に並びます。

``L_sub = [2,3]`` としたときの例::

 y
 ^     4 5
 |     2 3
 .->x  0 1


``skew`` は y 方向にユニットセル1つ分動いたときのx 方向のズレです。
``L_sub = [3,2], skew = 1`` としたときの例 (罫線はユニットセルの区切り)::

 y    ---------
 ^     5|3 4 5|
 |     2|0 1 2|
 |    ---------
 |    |3 4 5|
 .->x |0 1 2| 
      -------


ボンドの情報は ``hamiltonian`` (``tenes_std``) や ``evolution`` (``tenes``) で与えられます。


``tensor.unitcell`` サブセクション
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

各サイトのテンソルがあらわす初期状態を指定します。
全体の初期状態は、各サイトの初期状態の直積で与えられます。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``index``, "サイト番号", 整数 or 整数のリスト
   ``physical_dim``, "サイトにあるテンソルの持つ physical bond の次元", 整数
   ``virtual_dim``, "サイトにあるテンソルの持つ virtual bond の次元",  整数 or 整数のリスト
   ``initial_state``, "初期状態の係数", 実数のリスト
   ``noise``, "テンソル要素のゆらぎの大きさ", 実数


``index`` にリストを渡すことによって、複数のサイトを同時に指定できます。
空のサイトは全サイトを意味します。

``virtual_dim`` にリストを渡すことで、4方向のボンド次元を個別に指定できます。
順番は、左(-x)、上(+y)、右(+x)、下(-y) の順番です。

``initial_state`` では
:math:`|\psi\rangle_i = \sum_\alpha A_\alpha |\alpha\rangle_i` 
の :math:`A_\alpha` の値を指定します。
ゼロのみからなる配列を渡した場合、乱数初期化します。
実際には、すべてのvirtual ボンドインデックスが0 である要素が、 :math:`T_{0,0,0,0}^\alpha = A_\alpha` のように初期化されます。
他の要素には ``[-noise, noise)`` の一様乱数が入力されます。

たとえば、 :math:`S=1/2` のとき、 :math:`S^z` 方向に向いた状態を初期値にしたい場合には `initial_state = [1.0, 0.0]` に、
:math:`S^x` 方向に向いた状態を初期値にしたい場合には `initial_state = [1.0, 1.0]` とします。
