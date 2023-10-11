.. highlight:: none

横磁場イジング模型の有限温度計算
-----------------------------------------------

ここでは正方格子のイジングモデルに対して、横磁場 ``hx`` をかけた場合の有限温度の計算例について紹介します。
このチュートリアルで使用する入力ファイルおよびスクリプトファイルは ``sample/08_finitetemperature`` に格納されています。
以下が入力ファイルの例です (``simple_ft_strong.toml`` ファイル)。ここでは、

.. literalinclude:: ../../../../sample/08_finitetemperature/simple_ft_strong.toml

とします( ```Jz = -1.0`` なので、 ``hx`` 基底状態は強磁性状態になります)。
実時間発展は ``mode`` を ``finite`` にすることで行うことができます。
ここでは、横磁場を ``hx = 2.0`` 、 ``tau = 0.01`` (逆温度の刻み幅は 2 ``tau`` )として有限温度計算をしています。
時間発展の様子をいくつかの横磁場で見るために、 ``simple_ft_middle.toml`` 、 ``simple_te_weak.toml`` 、``simple_ft_zero.toml`` という入力ファイルもサンプルとして用意しています。
また、これらを一括して計算するためのスクリプト ``run.sh`` も準備しています。
あらかじめ ``tenes`` などにパスを通した上で

::

    sh run.sh

により計算を実行します。1分程度で計算が終わります。
計算結果を可視化するため、エネルギー、比熱、磁化 ( :math:`S_x` , :math:`S_z` )を表示するためのスクリプト ``plot_e.plt`` 、 ``plot_c.plt``、 ``plot_mx.plt`` 、 ``plot_mz.plt`` を用意しています。

::

    gnuplot -persist plot_e.plt
    gnuplot -persist plot_c.plt
    gnuplot -persist plot_mx.plt
    gnuplot -persist plot_mz.plt

とすれば、 エネルギー、比熱、磁化 ( :math:`S_x` , :math:`S_z` ) が、それぞれプロットされます。
その結果を :numref:`fig_tutorial8_finitetemperature` に示します。
比較のため、量子モンテカルロ法を用いて計算した結果も一緒に表示しています ( ``ALPS/looper`` を使用)。

.. figure:: ../../img/tutorial_08_finitetemperature.*
	:name: fig_tutorial8_finitetemperature
	:width: 600px

	イジング模型の有限温度計算の図: (a) エネルギー、(b) 比熱、(c) :math:`S_x` 、(d) :math:`S_z`. 縦軸は物理量、横軸は温度を表す。
