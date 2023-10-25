.. highlight:: none

横磁場イジング模型の実時間発展
-----------------------------------------------

ここでは正方格子のイジングモデルに対して、横磁場 ``hx`` をかけた場合の実時間発展の計算例について紹介します。
このチュートリアルで使用する入力ファイルおよびスクリプトファイルは ``sample/07_timeevolution`` に格納されています。

最初に初期状態として、基底状態の計算を行います (``simple.toml`` ファイル)。ここでは、

.. literalinclude:: ../../../../sample/07_timeevolution/simple.toml

とします(``Jz = -1.0`` なので、 基底状態は強磁性状態になります)。
初期状態として基底状態を使用するため、 ``tensor_save = "save_tensor"`` として、基底状態のテンソルを保存しておきます。

次に、実時間発展を行うための入力ファイルを用意します。
実時間発展は ``mode`` を ``time`` にすることで行うことができます。
以下、入力ファイル例です( ``simple_te_strong.toml`` ファイル)。

.. literalinclude:: ../../../../sample/07_timeevolution/simple_te_strong.toml

ここでは、横磁場を ``hx = 2.0`` 、実時間発展の刻み幅を ``tau = 0.01`` として時間発展させています。
また、初期状態として先ほどの基底状態を用いるため、 ``tensor_load = "save_tensor"`` として、基底状態のテンソルを読み込みます。
時間発展の様子をいくつかの横磁場で見るために、 ``simple_te_middle.toml`` 、 ``simple_te_weak.toml`` という入力ファイルもサンプルとして用意しています。
また、これらを一括して計算するためのスクリプト ``run.sh`` も準備しています。
あらかじめ ``tenes`` などにパスを通した上で

::

    sh run.sh

により計算を実行します。数十秒で計算が終わります。
計算が終了したら、gnuplotを起動し、

::

    load 'plot.plt'

とすれば、 磁化 :math:`S_z` の時間発展の様子がプロットされます。
その結果を :numref:`fig_tutorial7_timeevolution` に示します。

.. figure:: ../../img/tutorial_07_timeevolution.*
	:name: fig_tutorial7_timeevolution
	:width: 600px

	イジング模型の実時間発展の図. 縦軸は磁化、横軸は時間を表す。


なお、時間発展が進むにつれてエンタングルメントが大きくなり、ある時点でテンソルネットワークの容量が波動関数を表現するのに足りなくなります。
今の場合、 ``t=4.25`` の ``hx=2.0`` の場合のジャンプがこの問題を示しています。
実際に使用する場合には、このような不連続性が存在しないかを確認の上、ジャンプが存在する場合にはテンソルネットワークの容量が足りなくならないように、
``virtual_dimension`` を大きくするなどの対策を行う必要があります。
例えば、 ``virtual_dimension = 10`` に変更して上述の計算を行うと、
:numref:`fig_tutorial7_te_D10` に示すように不連続性が消えることがわかります。

.. figure:: ../../img/tutorial_07_timeevolution_D10.*
	:name: fig_tutorial7_te_D10
	:width: 600px

	イジング模型の実時間発展の図. 縦軸は磁化、横軸は時間を表す。
	``virtual_dimension = 10`` とした場合の結果。
