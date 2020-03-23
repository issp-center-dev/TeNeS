.. highlight:: none

.. _sec-simple-format:

``tense_simple`` の入力ファイル
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
方向を示す番号を省略することで、すべての方向について一度に指定できます。
また、最後に `xyz` のうち一文字を追加するとイジング的な相互作用を指定できます。
同一ボンド・同一成分を2回以上指定するとエラー終了します。

.. image:: ../../img/J.*
   :width: 400px
   :align: center

双二次相互作用 :math:`B` も :math:`J` と同様にボンド依存性をもたせられます。

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
また、2サイト物理量として、ボンドハミルトニアンおよび最近接ボンド上の :math:`S^z` 相関、 :math:`S^x` 相関、 :math:`S^y` 相関が自動的に定義されます。

``lattice`` セクション
==========================

計算する格子を指定します。
正方格子 (square lattice) と 蜂の巣格子 (honeycomb lattice) , 三角格子 (triangular lattice) が定義されています。

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``type``, "格子名 (square lattice もしくは honeycomb lattice)", 文字列, --
   ``L``, ユニットセルのx 方向の大きさ, 整数, --
   ``W``, ユニットセルのy 方向の大きさ, 整数, ``L``
   ``initial``, 初期テンソル, 文字列, "random"
   ``noise``, 初期テンソル, 実数, 1e-2


ユニットセルは ``L`` かける ``W`` の大きさをもつ長方形の形をしています。
ユニットセル内のサイトは0から順番に番号付けされます。 x 方向から順に並びます。

``L = 2, W = 3`` としたときの例::

 y
 ^     4 5
 |     2 3
 .->x  0 1


``initial`` と ``noise`` は波動関数の初期状態を決めるパラメータです。
``initial`` としては ``"ferro", "antiferro", "random"`` が指定可能です。
``noise`` はテンソルの要素に付与されるゆらぎの大きさです。なお、 ``parameter.general`` で ``tensor_load`` が設定されている場合には、そちらが優先されます。

正方格子 square lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

正方格子 ``type = "square lattice"`` で ``L=2`` のときの
テンソルの並びと最近接、次近接、三次近接のボンドタイプの定義を
:numref:`fig_square_1st`, :numref:`fig_square_2nd`, :numref:`fig_square_3rd` に示します。
四角はテンソルを表し、細線が最近接ボンドを表します。
太矢印と数字は2種類のボンドを表します。
破線はひとつのユニットセルを表します。

.. figure:: ../../img/Square_1st.*
   :width: 200px
   :align: center
   :name: fig_square_1st

   正方格子の最近接ボンド。

.. figure:: ../../img/Square_2nd.*
   :width: 200px
   :align: center
   :name: fig_square_2nd

   正方格子の次近接ボンド。

.. figure:: ../../img/Square_3rd.*
   :width: 200px
   :align: center
   :name: fig_square_3rd

   正方格子の三次近接ボンド。


蜂の巣格子 honeycomb lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

蜂の巣格子 ``type = "honeycomb lattice"`` で ``L=2`` のときの
テンソルの並びと最近接、次近接、三次近接のボンドタイプの定義を
:numref:`fig_honeycomb_1st`, :numref:`fig_honeycomb_2nd`, :numref:`fig_honeycomb_3rd` に示します。
蜂の巣格子ではx 方向のユニットセルの大きさ ``L`` は偶数でなければいけません。
四角はテンソルを表し、細線が最近接ボンドを表します。
太矢印と数字は3種類のボンドを表します。
破線はひとつのユニットセルを表します。

.. figure:: ../../img/Honeycomb_1st.*
   :width: 200px
   :align: center
   :name: fig_honeycomb_1st

   蜂の巣格子の最近接ボンド。

.. figure:: ../../img/Honeycomb_2nd.*
   :width: 200px
   :align: center
   :name: fig_honeycomb_2nd

   蜂の巣格子の次近接ボンド。

.. figure:: ../../img/Honeycomb_3rd.*
   :width: 200px
   :align: center
   :name: fig_honeycomb_3rd

   蜂の巣格子の三次近接ボンド。

三角格子 triangular lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

三角格子 ``type = "triangular lattice"`` で ``L=2`` のときの
テンソルの並びと最近接、次近接、三次近接のボンドタイプの定義を
:numref:`fig_triangular_1st`, :numref:`fig_triangular_2nd`, :numref:`fig_triangular_3rd` に示します。
四角はテンソルを表し、細線が最近接ボンドを表します。
太矢印と数字は3種類のボンドを表します。
破線はひとつのユニットセルを表します。

.. figure:: ../../img/Triangular_1st.*
   :width: 200px
   :align: center
   :name: fig_triangular_1st

   三角格子の最近接ボンド。

.. figure:: ../../img/Triangular_2nd.*
   :width: 200px
   :align: center
   :name: fig_triangular_2nd

   三角格子の次近接ボンド。

.. figure:: ../../img/Triangular_3rd.*
   :width: 200px
   :align: center
   :name: fig_triangular_3rd

   三角格子の三次近接ボンド。


``parameter`` セクション
=========================

``tenes_simple`` では使われず、 ``tenes_std`` の入力ファイルとしてそのままコピーされます。

.. include:: ./parameter_section.rst

``correlation`` セクション
==========================

``tenes_simple`` では相関関数 ``C = <A(0)B(r)>`` はデフォルトでは計算されません。
相関関数を計算したい場合は、 ``tenes`` の入力ファイルで指定する ``correlation`` セクションと共通のフォーマットで指定することができます。詳細は、:doc:`expert_format` の ``correlation`` セクションをご覧ください。
