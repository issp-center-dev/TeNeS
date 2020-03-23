.. highlight:: none


更新回数など、 計算にあらわれる種々のパラメータを記述します。
サブセクションとして ``general``, ``simple_update``, ``full_update``,
``ctm``, ``random`` を持ちます。

simple update およびfull updateの虚時間刻み ``parameter.simple_update.tau`` と ``parameter.full_update.tau`` のみ、 ``tenes`` 本体ではなくスタンダードモード ``tenes_std`` で使われるパラメータです。

``parameter.general``
~~~~~~~~~~~~~~~~~~~~~~~~~~

全般的な設定をします。

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 20, 30, 10, 10

   ``is_real``, "すべてのテンソルを実数に制限するかどうか", 真偽値, false
   ``iszero_tol``, "テンソルの読み込みにおいてゼロとみなす絶対値カットオフ", 実数, 0.0
   ``output``, "物理量などを書き込むディレクトリ", 文字列, \"output\"
   ``tensor_save``, "最適化後のテンソルを書き込むディレクトリ", 文字列, \"\"
   ``tensor_load``, "初期テンソルを読み込むディレクトリ",       文字列, \"\"


- ``is_real``

  - ``true`` にするとテンソルの要素を実数に制限して計算を行います
  - 一つでも複素演算子があると計算が始まりません

- ``iszero_tol``
  - テンソル要素の実部・虚部の読み込みにおいて、絶対値が ``iszero_tol`` 以下はゼロとみなします

- ``output``

  - 物理量などの計算結果をこのディレクトリ以下に保存します
  - 空文字列の場合はカレントディレクトリ以下に保存します

- ``tensor_save``

  - 最適化後のテンソルをこのディレクトリ以下に保存します
  - 空文字列の場合は保存しません

- ``tensor_load``

  - 各種テンソルをこのディレクトリ以下から読み込みます
  - 保存したときと同じ並列度である必要があります
  - 空文字列の場合は読み込みません
  - 

``parameter.simple_update``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``tau``,           "虚時間発展演算子における虚時間刻み",         実数, 0.01
   ``num_step``,      "simple update の回数",                       整数, 0
   ``lambda_cutoff``, "simple update でゼロとみなす平均場のcutoff", 実数, 1e-12



``parameter.full_update``
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``tau``,                 "虚時間発展演算子における虚時間刻み",         実数, 0.01
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

   ``dimension``,                "角転送行列の virtual ボンドの次元",                              整数,   4
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
