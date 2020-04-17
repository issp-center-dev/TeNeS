.. highlight:: none

.. _sec-simple-format:

``tenes_simple`` の入力ファイル
---------------------------------

-  ファイルフォーマットは
   `TOML <https://github.com/toml-lang/toml/blob/master/versions/ja/toml-v0.5.0.md>`__ 形式

-  ``model``, ``lattice``, ``parameter``, ``correlation``
   の4つのセクションを持ちます。
   
   - ``parameter`` セクションはそのままスタンダードモードの入力へとコピーされます。


``model`` セクション
==========================

計算する模型を指定します。
スピン系 (spin) が定義済みです。

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``type``, 模型の種類, 文字列, --


模型の種類によって相互作用などのパラメータ名が変わります。

スピン系 spin
~~~~~~~~~~~~~~~~~~~~~

スピン系

.. math ::

 \mathcal{H} = \sum_{\langle ij \rangle}\left[\sum_\alpha^{x,y,z} J^\alpha_{ij} S^\alpha_i S^\alpha_j + B \left(\vec{S}_i\cdot\vec{S}_j\right)^2 \right] - \sum_i \left[ h S^z_i + \Gamma S^x_i - D \left(S^z_i\right)^2 \right]


一体項のパラメータは次の通り。

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``S``, "局所スピンの大きさ",       実数（整数もしくは半整数）, 0.5
   ``h``, "縦磁場 h",                 実数, 0.0
   ``G``, "横磁場 :math:`\Gamma` ",   実数, 0.0
   ``D``, "オンサイトスピン異方性 D", 実数, 0.0


交換相互作用 :math:`J` にはボンド依存性をもたせることができます。

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``J0``, "最近接・第0方向ボンドの交換相互作用", 実数, 0.0
   ``J1``, "最近接・第1方向ボンドの交換相互作用", 実数, 0.0
   ``J2``, "最近接・第2方向ボンドの交換相互作用", 実数, 0.0
   ``J0'``, "次近接・第0方向ボンドの交換相互作用", 実数, 0.0
   ``J1'``, "次近接・第1方向ボンドの交換相互作用", 実数, 0.0
   ``J2'``, "次近接・第2方向ボンドの交換相互作用", 実数, 0.0
   ``J0''``, "三次近接・第0方向ボンドの交換相互作用", 実数, 0.0
   ``J1''``, "三次近接・第1方向ボンドの交換相互作用", 実数, 0.0
   ``J2''``, "三次近接・第2方向ボンドの交換相互作用", 実数, 0.0


ボンドの方向は lattice セクションで定義される格子に依存します。
例えば正方格子では x方向 (0) と y方向 (1) の2種類のボンド方向ごとに定義できます。
方向を示す番号を省略することで、すべての方向について一度に指定することもできます。
また、最後に `xyz` のうち一文字を追加するとイジング的な相互作用を指定できます。
同一ボンド・同一成分を2回以上指定するとエラー終了します。

.. image:: ../../img/J.*
   :width: 300px
   :align: center

双二次相互作用 :math:`B` も :math:`J` と同様にボンド依存性をもたせらることもできます。

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``B0``, "最近接・第0方向ボンドの双二次相互作用", 実数, 0.0
   ``B1``, "最近接・第1方向ボンドの双二次相互作用", 実数, 0.0
   ``B2``, "最近接・第2方向ボンドの双二次相互作用", 実数, 0.0
   ``B0'``, "次近接・第0方向ボンドの双二次相互作用", 実数, 0.0
   ``B1'``, "次近接・第1方向ボンドの双二次相互作用", 実数, 0.0
   ``B2'``, "次近接・第2方向ボンドの双二次相互作用", 実数, 0.0
   ``B0''``, "三次近接・第0方向ボンドの双二次相互作用", 実数, 0.0
   ``B1''``, "三次近接・第1方向ボンドの双二次相互作用", 実数, 0.0
   ``B2''``, "三次近接・第2方向ボンドの双二次相互作用", 実数, 0.0


物理量測定に使われる1サイト物理量として、 :math:`S^z` と :math:`S^x` 、(``parameter.general.is_real = false`` ならば) :math:`S^y` を自動的に定義されます。
また、2サイト物理量として、ボンドハミルトニアン

.. math ::

 \mathcal{H}_{ij} = \left[\sum_\alpha^{x,y,z} J^\alpha_{ij} S^\alpha_i S^\alpha_j + B \left(\vec{S}_i\cdot\vec{S}_j\right)^2 \right] 
 - \frac{1}{z} \left[ H \left(S^z_i + S^z_j \right) + \Gamma \left(S^x_i + S^x_j\right) - D \left(\left(S^z_i\right)^2 + \left(S^z_j\right)^2 \right) \right],

および最近接ボンド上の相関 :math:`S_i^\alpha S_j^\alpha` ( :math:`\alpha=x, y, z` )が自動的に定義されます。

``lattice`` セクション
==========================

計算する格子を指定します。
正方格子 (square) と 三角格子 (triangular), 蜂の巣格子 (honeycomb) , かごめ格子 (kagome) が定義されています。

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``type``, "格子名 (square, triangular, honeycomb, もしくは kagome)", 文字列, --
   ``L``, ユニットセルのx 方向の大きさ, 整数, --
   ``W``, ユニットセルのy 方向の大きさ, 整数, ``L``
   ``virtural_dim``, ボンド次元, 整数, --
   ``initial``, 初期状態, 文字列, "random"
   ``noise``, 初期テンソルの揺らぎ, 実数, 1e-2


``initial`` と ``noise`` は波動関数の初期状態を決めるパラメータです。
なお、 ``parameter.general`` で ``tensor_load`` が設定されている場合には、そちらが優先され、テンソルをファイルから読み込みます。

- ``initial``

   - ``"ferro"``

     - 強磁性状態。 各サイトで :math:`S^z=S`` となる状態。

   - ``"antiferro"``

     - 反強磁性状態。
       正方格子、蜂の巣格子では :math:`S^z = S` と :math:`S^z = -S` が互いに並んだNeel 秩序。
       三角格子、かごめ格子では スピンが :math:`(\theta, \phi) = (0, 0), (2\pi/3, 0), (2\pi/3, \pi)` 方向に向いた120 度秩序。

   - ``"random"``

     - 各サイトバラバラなランダム状態。

- ``noise``

  - テンソルの要素に付与されるゆらぎの大きさ。

正方格子 square lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

正方格子 ``type = "square lattice"`` では、 サイトが :math:`(1,0)` 方向に ``L`` 個、 :math:`(0,1)` 方向に ``W`` 個並びます。
具体例として、 ``L=3, W=3`` のときのサイトの並びを :numref:`fig_square_lattice` (a) に示します。
また、最近接、次近接、三次近接のボンドタイプの定義を
:numref:`fig_square_lattice` (b), (c), (d) にそれぞれ示します。
青線は ``bondtype = 0`` のボンドを, 赤線は ``bondtype = 1`` のボンドを表します。


.. figure:: ../../img/SquareLattice.*
   :width: 550px
   :align: center
   :name: fig_square_lattice

   正方格子 (square) のサイト・ボンド。
   (a) ``L=3, W=3`` としたときのサイトの並び。
   (b) 最近接ボンド。 ``bondtype=0`` (青) は 0 度方向に、 ``bondtype=1`` (赤) は 90度方向に伸びる。
   (c) 次近接ボンド。 ``bondtype=0`` (青) は 45 度方向に、 ``bondtype=1`` (赤) は -45度方向に伸びる。
   (d) 三次近接ボンド。 ``bondtype=0`` (青) は 0 度方向に、 ``bondtype=1`` (赤) は 90度方向に伸びる。


三角格子 triangular lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

三角格子 ``type = "triangular lattice"`` では、 サイトが :math:`(1,0)` 方向に ``L`` 個、 :math:`(1/2, \sqrt{3}/2)` 方向に ``W`` 個並びます。
具体例として、 ``L=3, W=3`` のときのサイトの並びを :numref:`fig_triangular_lattice` (a) に示します。
また、最近接、次近接、三次近接のボンドタイプの定義を
:numref:`fig_triangular_lattice` (b), (c), (d) にそれぞれ示します。
青線は ``bondtype = 0`` のボンドを、 赤線は ``bondtype = 1`` のボンドを、 緑線は ``bondtype = 2`` のボンドを、それぞれ表します。

.. figure:: ../../img/TriangularLattice.*
   :width: 550px
   :align: center
   :name: fig_triangular_lattice

   三角格子 (triangular) のサイト・ボンド。
   (a) ``L=3, W=3`` としたときのサイトの並び。
   (b) 最近接ボンド。 ``bondtype=0`` (青) は 0 度方向に、 ``bondtype=1`` (赤) は 60度方向に、 ``bondtype=2`` (緑)は 120度方向にそれぞれ伸びる。
   (c) 次近接ボンド。 ``bondtype=0`` (青) は 90 度方向に、 ``bondtype=1`` (赤) は -30度方向に、 ``bondtype=2`` (緑)は 30度方向にそれぞれ伸びる。
   (d) 三次近接ボンド。 ``bondtype=0`` (青) は 0 度方向に、 ``bondtype=1`` (赤) は 60度方向に、 ``bondtype=2`` (緑)は 120度方向にそれぞれ伸びる。


蜂の巣格子 honeycomb lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

蜂の巣格子 ``type = "honeycomb lattice"`` では、 座標 :math:`(0, 0)` と :math:`(\sqrt{3}/2, 1/2)` の2つのサイトからなるユニットが、 :math:`(\sqrt{3},0)` 方向に ``L`` 個、 :math:`(1/2, 3/2)` 方向に ``W`` 個並びます。
具体例として、``L=3, W=3`` のときのサイトの並びを :numref:`fig_honeycomb_lattice` (a) に示します。
破線はユニットを表します。
また、最近接、次近接、三次近接のボンドタイプの定義を
:numref:`fig_honeycomb_lattice` (b), (c), (d) にそれぞれ示します。
青線は ``bondtype = 0`` のボンドを、 赤線は ``bondtype = 1`` のボンドを、 緑線は ``bondtype = 2`` のボンドを、それぞれ表します。

.. figure:: ../../img/HoneycombLattice.*
   :width: 550px
   :align: center
   :name: fig_honeycomb_lattice

   蜂の巣格子 (honeycomb) のサイト・ボンド。
   (a) ``L=3, W=3`` としたときのサイトの並び。 破線で表されるユニットが ``L`` かける ``W`` 個並ぶ。
   (b) 最近接ボンド。 ``bondtype=0`` (青) は 30 度方向に、 ``bondtype=1`` (赤) は 150度方向に、 ``bondtype=2`` (緑)は -90度方向にそれぞれ伸びる。
   (c) 次近接ボンド。 ``bondtype=0`` (青) は 120 度方向に、 ``bondtype=1`` (赤) は 60度方向に、 ``bondtype=2`` (緑)は 0度方向にそれぞれ伸びる。
   (d) 三次近接ボンド。 ``bondtype=0`` (青) は -30 度方向に、 ``bondtype=1`` (赤) は -150度方向に、 ``bondtype=2`` (緑)は 90度方向にそれぞれ伸びる。


かごめ格子 kagome lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

かごめ格子 ``type = "kagome lattice"`` では、 座標 :math:`(0, 0)`, :math:`(1, 0)`, :math:`(1/2, \sqrt{3}/2)` の3つのサイトからなるユニット（上向き三角）が、 :math:`(2,0)` 方向に ``L`` 個、 :math:`(1, \sqrt{3})` 方向に ``W`` 個並びます。
具体例として、``L=3, W=3`` のときのサイトの並びを :numref:`fig_kagome_lattice` (a) に示します。
破線はユニットは破線を表します。
また、最近接、次近接、三次近接のボンドタイプの定義を
:numref:`fig_kagome_lattice` (b), (c), (d) にそれぞれ示します。
青線は ``bondtype = 0`` のボンドを、 赤線は ``bondtype = 1`` のボンドを表します。


.. figure:: ../../img/KagomeLattice.*
   :width: 550px
   :align: center
   :name: fig_kagome_lattice

   かごめ格子 (kagome) のサイト・ボンド。
   (a) ``L=3, W=3`` としたときのサイトの並び。 破線で表されるユニットが ``L`` かける ``W`` 個並ぶ。
   (b) 最近接ボンド。 ``bondtype=0`` (青) は 上向き三角形を、 ``bondtype=1`` (赤) は 下向き三角形を作る。
   (c) 次近接ボンド。
   (d) 三次近接ボンド。 ``bondtype=0`` (青) はサイトを横切り、 ``bondtype=1`` (赤) は横切らない。

``parameter`` セクション
=========================

``tenes_simple`` では使われず、 ``tenes_std`` の入力ファイルとしてそのままコピーされます。

.. include:: ./parameter_section.rst

``correlation`` セクション
==========================

``tenes_simple`` では相関関数 ``C = <A(0)B(r)>`` はデフォルトでは計算されません。
相関関数を計算したい場合は、 ``tenes`` の入力ファイルで指定する ``correlation`` セクションと共通のフォーマットで指定することができます。詳細は、:doc:`expert_format` の ``correlation`` セクションをご覧ください。
