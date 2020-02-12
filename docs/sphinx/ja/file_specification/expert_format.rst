.. highlight:: none

``tenes`` の入力ファイル
---------------------------------

-  ファイルフォーマットは
   `TOML <https://qiita.com/minoritea/items/c0de47b8beb813c655d4>`__
   形式
-  ``parameter``, ``lattice``, ``evolution``, ``observable``, ``correlation``
   の5つのセクションを持ちます。

   -  将来的にはファイル名を指定することで分割可能にする予定です。

``parameter`` セクション
========================

更新回数など、種々のパラメータを記述します。このセクションでは、各々のパラメータにデフォルト値が存在します。
サブセクションとして ``tensor``, ``simple_update``, ``full_update``,
``ctm``, ``random`` を持ちます。

``parameter.tensor``
~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

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

   ``CHI``,                      "角転送行列の virtual ボンドの次元",                              整数,   4
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

    [parameter.simple_update]
    num_step = 1000

    [parameter.full_update]
    num_step = 1

    [parameter.ctm]
    CHI  = 16  # env_dim
    iteration_max = 5


``tensor`` セクション
========================

「ユニットセル」の情報を記述します。
ユニットセルは ``Lx`` かける ``Ly`` の大きさをもつ長方形の形をしています。
また、サブセクション ``unitcell`` を持ちます。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``L_sub``, "ユニットセルの大きさ", 整数または整数のリスト


``L_sub`` として2つの整数からなるリストを渡した場合、はじめの要素が ``Lx`` に、もう片方が ``Ly`` になります。
3つ以上の要素からなるリストを渡した場合にはエラー終了します。
``L_sub`` として整数を渡した場合、 ``Lx`` と ``Ly`` とが等しくなります。

ユニットセル内のサイトは0から順番に番号付けされます。 x 方向から順に並びます。

``L_sub = [2,3]`` としたときの例::

 y
 ^     4 5
 |     2 3
 .->x  0 1


ボンドの情報は ``evolution`` や ``observable`` で与えられます。


``tensor.unitcell`` サブセクション
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


``observable`` セクション
==========================

物理量測定に関する諸々を記述します。
``onesite`` と ``twosite`` の2種類のサブセクションを持ちます。


``observable.onesite``
~~~~~~~~~~~~~~~~~~~~~~~~~

ひとつのサイト上で定義される物理量を示す一体演算子を定義します。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``group``,    "演算子の識別番号",   整数
   ``sites``,    "サイト番号",         整数 or 整数のリスト
   ``dim``,      "演算子の次元",       整数
   ``elements``, "演算子の非ゼロ要素", 文字列

``group`` はonesite 演算子の識別番号です。同じ番号を持つ演算子を複数定義することもできます。

``sites`` は演算子が作用するサイト番号です。リストを渡すことで複数同時に定義できます。
空リストは全サイトを意味します。

``dim`` は演算子の次元です。

``elements`` は演算子の非ゼロ要素を指定する文字列です。
1つの要素は2つの整数と2つの浮動小数点数を空白区切りからなる1つの行からなります。
最初の2つはそれぞれ演算子が作用する前と後の状態番号を、
あとの2つはそれぞれ演算子の要素の実部と虚部を示します。

例えば S=1/2 のSz 演算子は次のように定義されます。::

  [[observable.onesite]]
  group = 0
  sites = []
  dim = 2
  elements = """
  0 0  0.5 0.0
  1 1  -0.5 0.0
  """


``observable.twosite``
~~~~~~~~~~~~~~~~~~~~~~~~~

ふたつのサイト上で定義される物理量を示す演算子を定義します。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``group``,    "演算子の識別番号",   整数
   ``bonds``,    "ボンド",             文字列
   ``dim``,      "演算子の次元",       整数のリスト
   ``elements``, "演算子の非ゼロ要素", 文字列

``group`` は twosites 演算子の識別番号です。同じ番号を持つ演算子を複数定義することもできます。

``bonds`` は演算子が作用するサイト対の集合を表す文字列です。
4つの整数からなる1行が1つのサイト対を意味します。
最初の2つの整数は それぞれのサイト (source, target) のサイト番号です。
あとの2つの整数は target site が属するユニットセルの x, y 方向のオフセットです。

``dim`` は演算子の次元、すなわち作用するサイトの取りうる状態数です。

``elements`` は演算子の非ゼロ要素を指定する文字列です。
1つの要素は4つの整数と2つの浮動小数点数を空白区切りからなる1つの行からなります。
最初の2つは演算子が作用する前の source site, target site の状態番号を、
つぎの2つは演算子が作用した後の source site, target site の状態番号を、
最後の2つはそれぞれ演算子の要素の実部と虚部を示します。

例えば S=1/2 のハイゼンベルグハミルトニアンは次のように定義されます。::

  [[observable.twosite]]
  group = 0
  dim = [2, 2]
  bonds = """
  0 1 0 0
  0 2 0 0
  1 0 1 0
  1 3 0 0
  2 3 0 0
  2 0 0 1
  3 2 1 0
  3 1 0 1
  """
  elements = """
  0 0 0 0  0.25 0.0
  1 0 1 0  -0.25 0.0
  0 1 1 0  0.5 0.0
  1 0 0 1  0.5 0.0
  0 1 0 1  -0.25 0.0
  1 1 1 1  0.25 0.0
  """


``evolution`` セクション
========================

simple update, full update で使う虚時間発展演算子を記述します。
次のようなフィールドを持つ ``simple``, ``full`` の2つのサブセクションを持ちます。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``source_site``, "source site の番号",                      整数
   ``source_leg``,  "source site から見た target site の方向", 整数
   ``dimensions``,  "虚時間発展演算子テンソルの次元",          整数のリスト
   ``elements``,    "虚時間発展演算子テンソルの非ゼロ要素",    文字列


例 ::

  [evolution]
  [[evolution.simple]]
  source_site = 0
  source_leg = 2
  dimensions = [2, 2, 2, 2]
  elements = """
  0 0 0 0  0.9975031223974601 0.0
  1 0 1 0  1.0025156589209967 0.0
  0 1 1 0  -0.005012536523536871 0.0
  1 0 0 1  -0.005012536523536871 0.0
  0 1 0 1  1.0025156589209967 0.0
  1 1 1 1  0.9975031223974601 0.0
  """


``correlation`` セクション
==========================

相関関数 ``C = <A(0)B(r)>`` を指定するセクション

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``r_max``,     "相関関数の距離 :math:`r` の最大値", 整数
   ``operators``, "相関関数を測る演算子 A, B の番号", 整数のリストのリスト

演算子は ``observable`` セクションで指定したものが用いられます。

例
~~

::

    [correlation]
    r_max = 5
    operators = [[0,0], [0,1], [1,1]]
