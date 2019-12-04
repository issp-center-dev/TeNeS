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

更新回数など、種々のパラメータを記述します。このセクションのみ、各々のパラメータにデフォルト値が存在します。
サブセクションとして ``tensor``, ``simple_update``, ``full_update``,
``ctm``, ``random`` を持ちます。

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
========================

「ユニットセル」の情報を記述します。
ユニットセルは ``Lx`` かける ``Ly`` の大きさをもつ長方形の形をしています。

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


``evolution`` セクション
========================

simple update, full update で使う虚時間発展演算子を記述します。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``matrix``,        "虚時間発展演算子の行列表現",                                                   文字列のリスト
   ``simple_update``, "simple update における、虚時間発展演算子のインデックスと作用するボンドの順番", 文字列
   ``full_update``,   "full update における、虚時間発展演算子のインデックスと作用するボンドの順番",   文字列

``matrix``
~~~~~~~~~~

-  ひとつの文字列がひとつの行列を意味します。
-  列は１つ以上の空白で区切られ、行は１つ以上の改行で区切られます。
-  定義した順番がそのまま行列の番号に対応し、 ``*_update``
   での指定で使われます (0-origin)。

``simple_update``, ``full_update``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  1行が1回の演算子作用を表します。
-  各行は ``int int char int`` の４つのフィールドからなります。

   1. ボンドがつながるサイト
   2. ボンドがつながるサイト
   3. 横方向 (h) か縦方向 (v) か
   4. 演算子番号 (0-origin)

例
~~

.. code:: 

    [evolution]
    simple_update = """
    0 1 h 0
    3 2 h 0
    2 3 h 0
    1 0 h 0
    0 2 v 0
    3 1 v 0
    2 0 v 0
    1 3 v 0
    """

    full_update = """
    0 1 h 0
    3 2 h 0
    2 3 h 0
    1 0 h 0
    0 2 v 0
    3 1 v 0
    2 0 v 0
    1 3 v 0
    """

    matrix = [
    """
    0.9975031223974601 0.0 0.0 0.0
    0.0 1.0025156589209967 -0.005012536523536887 0.0
    0.0 -0.005012536523536888 1.0025156589209967 0.0
    0.0 0.0 0.0 0.9975031223974601
    """
    ]

``observable`` セクション
==========================

物理量測定に関する諸々を記述します。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``local_operator``,    "サイト演算子 (ex. Sz)",                          文字列のリスト
   ``hamiltonian``,       "ボンドハミルトニアン",                           文字列のリスト
   ``hamiltonian_bonds``, "ボンドハミルトニアンの種類と作用するボンドの組", 文字列

``local_operator``, ``hamiltonian``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``evolution.matrix`` と同様。
定義した順番がそのまま演算子・ハミルトニアンのインデックスに対応します。

``hamiltonian_bonds``
~~~~~~~~~~~~~~~~~~~~~

``evolution.simple_update`` と同様。

例
~~

::

    [observable]
    local_operator = [
    """
      0.5  0.0
      0.0 -0.5
    """,
    """
      0.0 0.5
      0.5 0.0
    """,
    ]

    hamiltonian_bonds = """
    0 1 h 0
    3 2 h 0
    2 3 h 0
    1 0 h 0
    0 2 v 0
    3 1 v 0
    2 0 v 0
    1 3 v 0
    """

    hamiltonian = [
    """
      0.25   0.0    0.0     0.0
      0.0   -0.25   0.5     0.0  
      0.0    0.5   -0.25    0.0  
      0.0    0.0    0.0     0.25
    """,
    ]

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
