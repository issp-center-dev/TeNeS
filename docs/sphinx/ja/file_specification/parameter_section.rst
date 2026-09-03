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
   ``fermion``,     "フェルミオン系として扱うかどうか(実験的機能)",               真偽値, false
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

- ``fermion``

  - **実験的機能です。** ``true`` にすると、サイトテンソルをフェルミオン(Z2 グレーディング付き)テンソルとして扱い、ボンド更新・測定時のフェルミオン交換符号を自動的に生成します
  - 各サイトの物理基底のパリティを ``tensor.unitcell`` セクションの ``parity`` キーで指定する必要があります。パリティはサイト内の生成演算子順序を固定した上で定義します(``tensor.unitcell`` セクション参照)
  - **演算子・ゲートの行列要素**: 2サイト演算子の ``elements``(添字 ``in1 in2 out1 out2``)は、順序付き2サイト Fock 基底

    .. math::
       |i_1 i_2\rangle = A^\dagger_1(i_1)\, A^\dagger_2(i_2)\, |0\rangle

    (:math:`A^\dagger(i)` は固定した内部順序による局所基底 :math:`i` の生成演算子列、サイト1 = ``source_site``)における行列要素 :math:`\langle o_1 o_2 | O | i_1 i_2 \rangle` です。この行列要素は**フェルミオンの反交換関係込みで**評価してください。スピンレス(``physical_dim = 2``)の最近接ホッピングでは余分な符号は現れませんが、1サイトに複数モードがある場合は現れます。例: スピンあり電子では :math:`c_{2\uparrow}` が :math:`c^\dagger_{1\downarrow}` を跨ぐため :math:`\langle \uparrow\downarrow, 0 |\, c^\dagger_{1\uparrow} c_{2\uparrow} \,| \downarrow, \uparrow \rangle = -1` となります。符号を手で書くのではなく、関与するモードの2サイト Fock 空間を数値的に構成(Jordan-Wigner)して行列要素を読み出すことを推奨します
  - サイト**間**の交換符号(2次元幾何由来のものを含む)はすべて TeNeS が生成します。ユーザーが与えるのは上記の符号込み局所行列のみです
  - 演算子・時間発展ゲートはフェルミオンパリティを保存する必要があり、パリティ奇の演算子は入力読み込み時にエラーになります
  - 現バージョンでは simple update による基底状態計算のみに対応しています(環境は CTM または平均場)。full update、``Use_RSVD``、``Simple_Gauge_Fix``、有限温度計算、実時間発展、マルチサイト演算子、``ops`` 形式2サイト観測量、距離2以上の2サイト演算子、skew セル、1幅セル(``LX < 2`` または ``LY < 2``)、相関関数、相関長は非対応で、入力読み込み時にエラーになります(相関長は強制的に無効化されます)
  - ``tensor_save`` / ``tensor_load`` はフェルミオン模式でも使えます。仮想ボンドの偶奇台帳が保存先の ``fermion.dat`` に書き出され、読み込み時に物理脚のパリティ・``virtual_dim``・``L_sub``・``skew``・テンソル自身のパリティが検証されます。``virtual_dim`` を変えての読み込みは非対応です(偶奇ブロックの構造が保てないため)

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
   ``use_onesite_rdm_convergence``, "CTMの収束判定に1サイト縮約密度行列の反復間距離も用いる",     真偽値, true
   ``iteration_max``,            "CTMの収束iterationの最大回数",                                   整数,   100
   ``projector_corner``,         "CTMのprojector計算で1/4角のテンソルのみを使う",                  真偽値, true
   ``use_rsvd``,                 "SVD を 乱択SVD で置き換えるかどうか",                            真偽値, false
   ``rsvd_oversampling_factor``, "乱択SVD 中に計算する特異値の数の、最終的に用いる数に対する比率", 実数,   2.0
   ``meanfield_env``,            "CTM ではなく simple update で得られる平均場環境を用いる。フェルミオン模式でも使用でき、2サイト観測量は単層縮約で評価されるため CTM 版より大幅に軽いが、精度は simple update 相当",        真偽値, false

``use_onesite_rdm_convergence`` が ``true`` の場合、角転送行列の特異値スペクトルに加えて、1サイト縮約密度行列の反復間距離を用いて CTM の収束を判定します。
距離は全サイト・全行列要素に対する max 要素ノルムを、その密度行列のトレースで割ったものです（密度行列の形の変化とノルムの相対変化の両方が含まれます）。
次元1の仮想ボンドによってテンソルネットワークが独立な部分に分解する構成、例えば横方向の仮想ボンド次元を1とした擬一次元的な計算では、角スペクトルだけでは収束を誤検知することがあるため、その対策として使われます（アルゴリズムの説明も参照してください）。
なお、そのような擬一次元的な計算を推奨しているわけではありません。TeNeS が対象としているのは2次元格子です。
``false`` にすると、従来の角スペクトルのみの収束判定に戻ります。
収束時の観測量の残差はおおむね ``convergence_epsilon`` 程度のスケールになるため、より高精度が必要な場合は ``convergence_epsilon`` を小さくしてください。

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
