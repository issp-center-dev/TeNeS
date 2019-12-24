.. highlight:: none

``tense_simple`` の入力ファイル
---------------------------------

-  ファイルフォーマットは
   `TOML <https://github.com/toml-lang/toml/blob/master/versions/ja/toml-v0.5.0.md>`__
   形式
-  ``model``, ``parameter``, ``lattice``, ``observable``, ``correlation``
   の5つのセクションを持ちます。 なお、 ``observable``, ``correlation`` セクションは、
   エキスパートモードの入力ファイルと共通です。

   -  将来的にはファイル名を指定することで分割可能にする予定です。


``parameter`` セクション
==========================

``parameter`` セクションの内容は、そのまま出力の ``parameter`` セクションにコピーされます。

また、 simple update と full update における虚時間発展演算子の虚時間刻み幅を
サブセクション ``simple_update``, ``full_update`` で指定できます。

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``tau``, "虚時間発展演算子における虚時間の刻み幅", 実数, 0.01

以下は ``tenes`` の入力ファイルと共通のパラメータになります。

``parameter.tensor``
~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``D``,        "中心テンソルがもつ virtual ボンドの次元",  整数,   2
   ``CHI``,      "角転送行列の virtual ボンドの次元",        整数,   4
   ``save_dir``, "最適化後のテンソルを書き込むディレクトリ", 文字列, \"\"
   ``load_dir``, "初期テンソルを読み込むディレクトリ",       文字列, \"\"


- ``save_dir``

  - 最適化後のテンソルをこのディレクトリ以下に保存します
  - 空文字列の場合は保存しません

- ``load_dir``

  - 各種テンソルをこのディレクトリ以下から読み込みます
  - 保存したときと同じ並列度である必要があります
  - 空文字列の場合は読み込みません

``parameter.simple_update``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``num_step``,              "simple update の回数",                       整数, 0
   ``lambda_cutoff``, "simple update でゼロとみなす平均場のcutoff", 実数, 1e-12

``parameter.full_update``
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``num_step``,            "full update の回数",                                                 整数,   0
   ``env_cutoff``,          "full update で環境テンソルを計算する際にゼロとみなす特異値のcutoff", 実数,   1e-12
   ``inverse_precision``,   "full update で擬似逆行列を計算する際にゼロとみなす特異値のcutoff",   実数,   1e-12
   ``convergence_epsilon``, "full update でtruncationの最適化を行う際の収束判定値",               実数,   1e-12
   ``iteration_max``,       "full update でtruncationの最適化を行う際のiterationの最大回数",      整数,   1000
   ``gauge_fix``,           "テンソルのゲージを固定するかどうか",                                 真偽値, true
   ``fastfullupdate``,      "Fast full update にするかどうか",                                    真偽値, true

``parameter.ctm``
~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``projector_cutoff``,         "CTMのprojectorを計算する際にゼロとみなす特異値のcutoff",         実数,   1e-12
   ``convergence_epsilon``,      "CTMの収束判定値",                                                実数,   1e-10
   ``iteration_max``,            "CTMの収束iterationの最大回数",                                   整数,   100
   ``projector_corner``,         "CTMのprojector計算で1/4角のテンソルのみを使う",                  真偽値, true
   ``use_rsvd``,                 "SVD を 乱択SVD で置き換えるかどうか",                            真偽値, false
   ``rsvd_oversampling_factor``, "乱択SVD 中に計算する特異値の数の、最終的に用いる数に対する比率", 実数,   2.0


``parameter.random``
~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``seed``, "テンソルの初期化や乱択SVD に用いる疑似乱数生成器のシード", 整数, 11

MPI 並列において、各プロセスは ``seed`` にプロセス番号を足した数を実際のシードとして持ちます。

例
~~

::

    [parameter]
    [parameter.tensor]
    D  = 4     # tensor_dim
    CHI  = 16  # env_dim

    [parameter.simple_update]
    num_step = 1000

    [parameter.full_update]
    num_step = 1

    [parameter.ctm]
    iteration_max = 5


``lattice`` セクション
==========================

計算する格子を指定します。
正方格子 (square lattice) と 蜂の巣格子 (honeycomb lattice) が定義されています。

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``type``, "格子名 (square lattice もしくは honeycomb lattice)", 文字列, --
   ``L_sub``, ユニットセルの大きさ, 整数もしくは2つの整数からなるリスト, --
   ``initial``, 初期テンソル, 文字列, "random"
   ``noise``, 初期テンソル, 実数, 1e-2


ユニットセルは ``Lx`` かける ``Ly`` の大きさをもつ長方形の形をしています。
``L_sub`` として2つの整数からなるリストを渡した場合、はじめの要素が ``Lx`` に、もう片方が ``Ly`` になります。
3つ以上の要素からなるリストを渡した場合にはエラー終了します。
``L_sub`` として整数を渡した場合、 ``Lx`` と ``Ly`` とが等しくなります。

ユニットセル内のサイトは0から順番に番号付けされます。 x 方向から順に並びます。

``L_sub = [2,3]`` としたときの例::

 y
 ^     4 5
 |     2 3
 .->x  0 1


``initial`` と ``noise`` は波動関数の初期状態を決めるパラメータです。
``initial`` としては ``"ferro", "antiferro", "random"`` が指定可能です。
``noise`` はテンソルの要素に付与されるゆらぎの大きさです。

正方格子 square lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ボンドは水平方向 (0) と垂直方向 (1) の2種類あります（下図の ``-`` と ``|`` に対応）。

``L_sub = 2`` のときのユニットセルは次の通り::

 0   1
 |   |
 2 - 3 - 2
 |   | 
 0 - 1 - 0


蜂の巣格子 honeycomb lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ユニットセルの大きさ (``L_sub`` の各要素)は偶数でなければなりません。

ボンドはx (0), y (1), z (2) の3種類あります（それぞれ、下図の ``-``, ``~``, ``|`` に対応）。
偶数番のサイトは右(x)、左(y)、上(z) に伸びるボンドを持ち、
奇数番のサイトは左(x)、右(y)、下(z) に伸びるボンドを持ちます。

``L_sub = 2`` のときのユニットセルは次の通り::

 0   1
     |
 2 ~ 3 - 2
 |   
 0 - 1 ~ 0


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

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``S``,  "局所スピンの大きさ",       実数,                     0.5
   ``Jx``, "交換相互作用J のx 成分",   実数もしくは実数のリスト, 1.0
   ``Jy``, "交換相互作用J のy 成分",   実数もしくは実数のリスト, 1.0
   ``Jz``, "交換相互作用J のz 成分",   実数もしくは実数のリスト, 1.0
   ``BQ``, "双二次相互作用B",          実数もしくは実数のリスト, 0.0
   ``h``,  "縦磁場 h",                 実数,                     0.0
   ``G``,  "横磁場 :math:`\Gamma` ",   実数,                     0.0
   ``D``,  "オンサイトスピン異方性 D", 実数,                     0.0


交換相互作用および双二次相互作用としてリストを与えることで、格子ボンドの種類ごとに相互作用の大きさを変えることができます。
リストの要素数が格子ボンドの種類より少ない場合、足りない分は指定された最後の要素で埋められます。


``observable`` セクション
==========================

``tenes_simple`` ではデフォルトでは、物理量測定に使われる局所物理量として、 :math:`S^z` と :math:`S^x` が定義されます。 より詳細な物理量測定は、 ``tenes`` の入力ファイルで指定する ``observable`` セクションと共通のフォーマットで指定することができます。詳細は、:doc:`expert_format` の ``observable`` セクションをご覧ください。

``correlation`` セクション
==========================

``tenes_simple`` では相関関数 ``C = <A(0)B(r)>`` はデフォルトでは計算されません。
相関関数を計算したい場合は、 ``tenes`` の入力ファイルで指定する ``correlation`` セクションと共通のフォーマットで指定することができます。詳細は、:doc:`expert_format` の ``correlation`` セクションをご覧ください。
