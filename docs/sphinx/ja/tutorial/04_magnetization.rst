.. highlight:: none

三角格子・正方格子ハイゼンベルク模型の磁化過程
-----------------------------------------------

次に三角格子上で定義されたスピン\ :math:`S=1/2`\ の量子ハイゼンベルク模型の磁化過程の計算を紹介します。
ハミルトニアンは以下のようになります:

.. math::

   \begin{aligned}
   H = J \sum_{\langle i,j \rangle}\sum_{\alpha}^{x,y,z} {S}_i^{\alpha} {S}_j^{\alpha} - h \sum_i S_i^z\end{aligned}

ここで\ :math:`\langle i, j\rangle`\ は隣接サイトの組を表し、\ :math:`h`\ は\ :math:`z`\ 方向にかけられた外部磁場の大きさを表します。
この模型の基底状態を計算し、ユニットセルの平均磁化 \ :math:`\langle S_z \rangle\equiv \frac{1}{N_u}\sum_i^{N_u} \langle S_i^z \rangle`\ を磁場\ :math:`h`\ の関数として求めてみましょう ( :math:`N_u` はユニットセル内のサイト数)。

この計算を行うには、 ``sample/04_magnetization`` のディレクトリ内にある toml ファイル ``basic.toml`` と、pythonスクリプト ``tutorial_magnetization.py`` を利用します。 ``basic.toml`` ファイルには、模型の設定やパラメータなどが記述されています。

.. literalinclude:: ../../../../sample/04_magnetization/basic.toml

``lattice`` セクションで三角格子を指定しており、ユニットセルの大きさは\ :math:`3\times 3`\ を指定しています。
ここでは計算を軽くするために、 ``simple update`` だけを行っており、虚時間の刻み幅\ :math:`\tau`\ は\ :math:`\tau = 0.01`\ としています。また簡単のため、\ :math:`J=1`\ としています。この基本設定ファイルを用いて、 ``tutorial_magnetization.py`` では磁場を掃引したときの磁化を計算します。

.. literalinclude:: ../../../../sample/04_magnetization/tutorial_magnetization.py

このスクリプトでは、磁場\ :math:`h`\ を0から5まで0.25刻みで変化させ、基底状態のエネルギーと\ :math:`\langle S_z \rangle`\ を計算して、 ``energy.dat`` および ``magnetization.dat`` に出力します。 ``simple update`` の時間ステップ数を\ :math:`100`, :math:`200`, :math:`500`, :math:`1000`, :math:`2000`\ と変化させたときの様子を見るために、各磁場でステップ数を変えた計算も行っています。
計算量を減らすために、少ないステップ数で得られた波動関数の情報を ``tensor_save`` に保存し、それをより多いステップ数の計算の初期状態としてとっています。例えば、最初に時間ステップ数を100とした計算を行って結果を出力したあと、ステップ数100の計算終了時の波動関数からさらにステップ数100の計算を行って、結果的にステップ数200の計算結果を得ます。

実際に実行してみましょう。あらかじめ ``tenes`` などにパスを通した上で

::

    python tutorial_magnetization.py

により計算を実行します。ノートPC(シングルプロセッサ)では数時間程度の計算量となります。
計算が終了したら、gnuplotを起動し、

::

    load 'plot.gp'

とすれば、:numref:`fig_tutorial5_triangular` の右図のような磁化カーブが得られます。同様に

::

    load 'plot_ene.gp'

とすれば、:numref:`fig_tutorial5_triangular` の左図のような基底エネルギーのグラフが得られます。

十分なステップ数(例えばステップ数2000)の計算結果からわかるように、磁化過程には飽和磁化\ :math:`\langle S_z \rangle=0.5`\ の\ :math:`1/3`\ の磁化のところで、プラトー構造が生じます。
このプラトー上では、3つの格子上のスピンが\ :math:`\uparrow`,
:math:`\uparrow`,
:math:`\downarrow`\ と磁化した周期構造を形成し、スピンギャップが生じています。このプラトー構造は三角格子特有のものです。
実際に計算精度がでているかどうかをみるには、エネルギーのステップ依存性が参考になります。
理想的にはステップ数を増やすほど基底エネルギーが下がるはずですが、一部の磁場領域では逆に基底エネルギーが増加します。これは計算精度があまりでていない兆候です。
ボンド次元を増やすなどして、より計算精度を高める必要があると推測されます。

.. figure:: ../../img/tutorial_5_triangular.*
     :name: fig_tutorial5_triangular
     :width: 600px

     三角格子量子ハイゼンベルク模型のエネルギー(左図)と磁化過程(右図)

では正方格子でも同じことをやってみましょう。 ``sample/04_magnetization`` のディレクトリ中にある toml ファイル ``basic_square.toml`` と、pythonスクリプト ``tutorial_magnetization_square.py`` を利用します。
``basic_square.toml`` は、 ``lattice`` セクションが以下のように変更されているほかは ``basic.toml`` と同じ内容です。
::

    [lattice]
    type = "square lattice"
    L = 2
    W = 2

実際に計算を行うには、 ::

   python tutorial_magnetization_square.py

とします。計算が終了したら、gnuplotを起動し、
::

    load 'plot_square.gp'

とすれば、 :numref:`fig_tutorial5_square` の右図のような磁化カーブが得られます。同様に
::

    load 'plot_square_ene.gp'

とすれば、:numref:`fig_tutorial5_square` の左図のような基底エネルギーのグラフが得られます。

.. figure:: ../../img/tutorial_5_square.*
     :name: fig_tutorial5_square
     :width: 600px

     正方格子量子ハイゼンベルク模型のエネルギー(左図)と磁化過程(右図)

ステップ数2000でほぼ収束しており、三角格子ハイゼンベルク模型と異なり、プラトー構造は現れないことがわかります。
エネルギーは概ね、ステップ数を増加させると減少するため、ある程度計算精度がでていると推測されます。
