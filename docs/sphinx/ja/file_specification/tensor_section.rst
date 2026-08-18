.. highlight:: none

「ユニットセル」の情報を記述します(ボンドの情報は ``hamiltonian`` (``tenes_std``) や ``evolution`` (``tenes``) で与えます)。
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

.. figure:: ../../img/tensor_sec_fig1.*
   :width: 150px

   ``L_sub = [2,3]`` としたときの例

``skew`` は y 方向にユニットセル1つ分動いたときのx 方向のズレです。


.. figure:: ../../img/tensor_sec_fig2.*
   :width: 400px

   ``L_sub = [3,2], skew = 1`` としたときの例 (罫線はユニットセルの区切り)



``tensor.unitcell`` サブセクション
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

サイトテンソル :math:`T_{ijkl\alpha}^{(n)}` の情報を指定します。
ここで :math:`i,j,k,l` は virtual bond のインデックス、 :math:`\alpha` は physical bond のインデックス、 :math:`n` はサイト番号を意味します。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``index``,         "サイト番号",                                     整数 or 整数のリスト
   ``physical_dim``,  "サイトテンソルの physical bond の次元",          整数
   ``virtual_dim``,   "サイトテンソルの virtual bond の次元 :math:`D`", 整数 or 整数のリスト
   ``initial_state``, "初期状態",                                       実数のリスト
   ``noise``,         "初期テンソルのゆらぎの大きさ",                   実数
   ``parity``,        "物理基底のフェルミオンパリティ(fermion モード時のみ)", 整数のリスト


``index`` にリストを渡すことによって、複数のサイトを同時に指定できます。
空のサイト ``[]`` は全サイトを意味します。

``parity`` は ``parameter.general.fermion = true`` のとき必須で、
``physical_dim`` 個の 0(偶)または 1(奇)を並べます。
``parity[i]`` は物理基底の第 ``i`` 状態に含まれるフェルミオン数の偶奇です。
たとえばスピンレスフェルミオンの基底 :math:`\{|0\rangle, |1\rangle\}` は
``parity = [0, 1]`` です。

1サイトに複数のフェルミオンモードがある場合(例: スピンあり電子、
``physical_dim = 4``)は、**サイト内の生成演算子の順序を1つ固定し、
以後すべてで同じ順序を使ってください**。各基底状態は、固定した内部順序で
生成演算子を真空に作用させたものとして定義します。たとえば内部順序を
:math:`(\uparrow, \downarrow)` とすると

.. math::
   |0\rangle,\quad
   |{\uparrow}\rangle = c^\dagger_\uparrow |0\rangle,\quad
   |{\downarrow}\rangle = c^\dagger_\downarrow |0\rangle,\quad
   |{\uparrow\downarrow}\rangle = c^\dagger_\uparrow c^\dagger_\downarrow |0\rangle

であり、``parity = [0, 1, 1, 0]`` となります。
選んだ生成演算子順序は ``parity`` のリスト自体には現れませんが、
演算子・ゲートの行列要素の符号を決めます
(``parameter.general`` セクションの ``fermion`` の項を参照)。
すべての演算子・ゲート・初期状態を同一の規約で書く必要があり、
規約が混在すると符号が黙って間違います。

初期プロダクト状態は各サイトでパリティ偶である必要があり、
パリティ奇の ``initial_state`` はエラーになります。

``virtual_dim`` にリストを渡すことで、4方向のボンド次元を個別に指定できます。
順番は、左(-x)、上(+y)、右(+x)、下(-y) の順番です。

系全体の初期状態 :math:`|\Psi\rangle` は、各サイト :math:`i` の初期状態 :math:`|\Psi_i\rangle` の直積で与えられます。

.. math::
   |\Psi\rangle = \otimes_i |\Psi_i\rangle

サイトテンソルはこの直積状態を表現するように初期化されます。
``initial_state`` では各サイト :math:`i` の初期状態
:math:`|\Psi_i\rangle = \sum_\alpha A_\alpha |\alpha\rangle_i` における実展開係数 :math:`A_\alpha` の値を指定します。
係数は自動的に規格化されます。
テンソル自体は、 すべてのvirtual ボンドインデックスが0 である要素が、 :math:`T_{0000\alpha} = A_\alpha` のように初期化されます。
他の要素には ``[-noise, noise)`` の一様乱数が互いに独立に入力されます。
たとえば、 :math:`S=1/2` のとき、 :math:`S^z` 方向に向いた状態 :math:`|\Psi_i\rangle = |\uparrow\rangle = |0\rangle` を初期値にしたい場合には ``initial_state = [1.0, 0.0]`` に、
:math:`S^x` 方向に向いた状態 :math:`|\Psi_i\rangle = \left(|\uparrow\rangle + |\downarrow\rangle\right)/\sqrt{2}` を初期値にしたい場合には ``initial_state = [1.0, 1.0]`` とします。

``initial_state`` にゼロのみからなる配列を渡した場合、テンソルのすべての要素が独立に ``[-noise, noise)`` で乱数初期化されます。

