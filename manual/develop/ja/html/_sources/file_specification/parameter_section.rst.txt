.. highlight:: none


更新回数など、 計算にあらわれる種々のパラメータを記述します。
サブセクションとして ``general``, ``simple_update``, ``full_update``,
``ctm``, ``random`` を持ちます。

``parameter.general``
~~~~~~~~~~~~~~~~~~~~~~~~~~

``tenes`` の全般的な設定パラメータ

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 20, 30, 10, 10

   ``mode``,        "計算モード",                                                   文字列, \"ground state\"
   ``is_real``,     "すべてのテンソルを実数に制限するかどうか",                     真偽値, false
   ``iszero_tol``,  "演算子テンソルの読み込みにおいてゼロとみなす絶対値カットオフ", 実数,   0.0
   ``measure``,     "基底状態計算において、物理量測定をするかどうか",               真偽値, true
   ``measure_interval``, "実時間発展・有限温度計算において物理量を測定する頻度",    整数 or 整数のリスト, 10 
   ``output``,      "物理量などを書き込むディレクトリ",                             文字列, \"output\"
   ``tensor_save``, "最適化後のテンソルを書き込むディレクトリ",                     文字列, \"\"
   ``tensor_load``, "初期テンソルを読み込むディレクトリ",                           文字列, \"\"


- ``mode``

  - 計算モードを指定します
  - ``"ground state"``

    - 基底状態計算
    - ``tenes_std`` は虚時間発展演算子 :math:`U(\tau) = \exp(-\tau \mathcal{H})` を計算します

  - ``"time evolution"``

    - 実時間発展計算
    - ``tenes_std`` は実時間発展演算子 :math:`U(t) = \exp(-it \mathcal{H})` を計算します

  - ``"finite temperature"``

    - 有限温度計算
    - ``tenes_std`` は虚時間発展演算子 :math:`U(\tau) = \exp(-\tau \mathcal{H})` を計算します

- ``is_real``

  - ``true`` にするとテンソルの要素を実数に制限して計算を行います
  - 一つでも複素演算子があるとエラー終了します

- ``iszero_tol``

  - 各種演算子テンソル要素の実部・虚部の読み込みにおいて、絶対値が ``iszero_tol`` 以下はゼロとみなします

- ``measure``

  - 基底状態計算において、 ``false`` にすると物理量計算・保存をスキップします
  - 実行時間 ``time.dat`` は常に保存されます

- ``measure_interval``

  - 実時間発展計算および有限温度計算において、 物理量を測定する頻度を指定します
  - ``measure_interval`` ステップ計算した後に物理量を測定します

- ``output``

  - 物理量などの計算結果をこのディレクトリ以下に保存します
  - 空文字列の場合はカレントディレクトリに保存します

- ``tensor_save``

  - 最適化後のテンソルをこのディレクトリ以下に保存します
  - 空文字列の場合は保存しません

- ``tensor_load``

  - 初期テンソルをこのディレクトリ以下から読み込みます
  - 空文字列の場合は読み込みません


``parameter.simple_update``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

simple update に関するパラメータ

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``tau``,                       "(虚)時間発展演算子における(虚)時間刻み :math:`\tau`",         実数 or 実数のリスト,   0.01
   ``num_step``,                  "simple update の回数",                                        整数 or 整数のリスト,   0
   ``lambda_cutoff``,             "simple update において平均場 :math:`\lambda` の切り捨て閾値", 実数,   1e-12
   ``gauge_fix``,                 "テンソルのゲージを固定するかどうか",                          真偽値, false
   ``gauge_maxiter``,             "ゲージ固定操作のループ最大数",                                整数,   100
   ``gauge_convergence_epsilon``, "ゲージ固定操作の収束判定値",                                  実数,   1e-2


- ``tau``

  - (虚)時間発展演算子における(虚)時間刻み :math:`\tau` を指定します

    - ``tenes_std`` では時間発展演算子を計算するために用いられます
    - ``tenes`` では各ステップでの経過時間・逆温度を求めるために用いられます

      - For finite temperature calculation, note that the inverse temperature increase :math:`2\tau` at a step because :math:`\rho(\beta + 2\tau) = U(\tau)\rho(\beta)\bar{U}(\tau)`
      - 有限温度計算の場合、 :math:`\rho(\beta + 2\tau) = U(\tau)\rho(\beta)\bar{U}(\tau)` なので、 ステップごとに逆温度は :math:`2\tau` だけ増加することに注意してください。

  - リストを指定すると、時間発展演算子のグループごとに刻み幅を変えることができます

- ``num_step``

  - simple update の回数を指定します
  - リストを指定すると、時間発展演算子のグループごとに回数を変えることができます


``parameter.full_update``
~~~~~~~~~~~~~~~~~~~~~~~~~

full update に関するパラメータ

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``tau``,                 "(虚)時間発展演算子における(虚)時間刻み :math:`\tau`",                実数 or 実数のリスト,   0.01
   ``num_step``,            "full update の回数",                                                 整数 or 整数のリスト,   0
   ``env_cutoff``,          "full update で環境テンソルを計算する際にゼロとみなす特異値のcutoff", 実数,   1e-12
   ``inverse_precision``,   "full update で擬似逆行列を計算する際にゼロとみなす特異値のcutoff",   実数,   1e-12
   ``convergence_epsilon``, "full update でtruncationの最適化を行う際の収束判定値",               実数,   1e-6
   ``iteration_max``,       "full update でtruncationの最適化を行う際のiterationの最大回数",      整数,   100
   ``gauge_fix``,           "テンソルのゲージを固定するかどうか",                                 真偽値, true
   ``fastfullupdate``,      "Fast full update にするかどうか",                                    真偽値, true

``parameter.ctm``
~~~~~~~~~~~~~~~~~

角転送行列 (CTM) に関するパラメータ

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``dimension``,                "CTM のボンド次元 :math:`\chi`",                                  整数,   4
   ``projector_cutoff``,         "CTMのprojectorを計算する際にゼロとみなす特異値のcutoff",         実数,   1e-12
   ``convergence_epsilon``,      "CTMの収束判定値",                                                実数,   1e-6
   ``iteration_max``,            "CTMの収束iterationの最大回数",                                   整数,   100
   ``projector_corner``,         "CTMのprojector計算で1/4角のテンソルのみを使う",                  真偽値, true
   ``use_rsvd``,                 "SVD を 乱択SVD で置き換えるかどうか",                            真偽値, false
   ``rsvd_oversampling_factor``, "乱択SVD 中に計算する特異値の数の、最終的に用いる数に対する比率", 実数,   2.0
   ``meanfield_env``,            "CTM ではなく simple update で得られる平均場環境を用いる",        真偽値, false

乱拓SVDを用いたテンソル繰り込み群の手法については、 S. Morita, R. Igarashi, H.-H. Zhao, and N. Kawashima, `Phys. Rev. E 97, 033310 (2018) <https://journals.aps.org/pre/abstract/10.1103/PhysRevE.97.033310>`_ を参照してください。


``parameter.random``
~~~~~~~~~~~~~~~~~~~~~

疑似乱数生成器に関するパラメータ

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 30, 30, 10, 10

   ``seed``, "テンソルの初期化や乱択SVD に用いる疑似乱数生成器のシード", 整数, 11

MPI 並列において、各プロセスは ``seed`` にプロセス番号を足した数を実際のシードとして持ちます。

例
~~

::

  [parameter]
  [parameter.general]
  is_real = true
  [parameter.simple_update]
  num_step = 100
  tau = 0.01
  [parameter.full_update]
  num_step = 0  # No full update
  tau = 0.01
  [parameter.ctm]
  iteration_max = 10
  dimension = 9 # CHI
