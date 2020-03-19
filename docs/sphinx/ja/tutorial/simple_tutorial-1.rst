.. highlight:: none

横磁場イジングモデル
----------------------------

ここでは横磁場イジングモデルに対して、横磁場を変化させた場合の計算例について紹介します。
入力ファイルの変数 ``G`` を用いることで横磁場の大きさを調整することが可能です。例えば、横磁場が0の場合には、

.. code::

   [parameter]
   [parameter.general]
   is_real = true

   [parameter.simple_update]
   num_step = 10000
   tau = 0.01

   [parameter.full_update]
   num_step = 0
   tau = 0.01

   [parameter.ctm]
   iteration_max = 10
   dimension = 10

   [lattice]
   type = "square lattice"
   L = 2
   W = 2
   virtual_dim = 2
   initial = "ferro"

   [model]
   type = "spin"
   Jz = -1.0
   Jx = 0.0
   Jy = 0.0
   G  = 0.0


とします(``Jz = -1.0`` なので、 ``G=0`` では強磁性状態になります)。入力ファイルを ``simple.toml`` とした場合、
   
.. code:: bash

   $ tenes_simple simple.toml
   $ tenes_std std.toml
   $ tenes input.toml

を実行することで計算が開始されます。
計算を実行すると、

.. code:: bash
   Number of Processes: 1
   Number of Threads / Process: 1
   Tensor type: real
   Start simple update
     10% [100/1000] done
     20% [200/1000] done
     30% [300/1000] done
     40% [400/1000] done
     50% [500/1000] done
     60% [600/1000] done
     70% [700/1000] done
     80% [800/1000] done
     90% [900/1000] done
     100% [1000/1000] done
   Start calculating observables
     Start updating environment
     Start calculating onesite operators
       Save onesite observables to output/onesite_obs.dat
     Start calculating twosite operators
       Save twosite observables to output/twosite_obs.dat
       Save observable densities to output/density.dat
       Save elapsed times to output/time.dat

   Onesite observables per site:
     Sz          = 0.297866964052 0
     Sx          = 0.386024172907 0
   Twosite observables per site:
     hamiltonian = -0.75730305866 0
     SzSz        = 0.21686921659 0
     SxSx        = 0.319350111777 0
     SySy        = -0.0477650003168 0
   Wall times [sec.]:
     simple update = 0.691005566
     full update   = 0
     environmnent  = 0.269068351
     observable    = 0.030814304

   Done.

のように計算が実行されます。
最初に並列化の情報およびテンソルの実虚が表示されます。
次に計算プロセスの実行状況が表示されます。
計算終了後、 1サイト演算子 ``Sz``,   ``Sx`` およびハミルトニアン ``hamiltonian`` , 最近接相関 ``SzSz``, ``SxSx``, ``SySy`` のサイトあたりの期待値が出力されます。
最後にフェーズごとの計算時間が出力されます(単位は秒)。
計算終了後は ``output`` フォルダに
``density.dat, parameters.dat, time.dat, onesite_obs.dat, twosite_obs.dat``
がそれぞれ出力されます。各出力ファイルの詳細は、ファイルフォーマットをご覧ください。
例えば ``<Sz>`` の値は、 ``onesite_obs.dat`` の2行目から読み取ることが可能です。
``G`` をパラメータとして0.1刻みで0-3.0まで振ったときの結果を下図に表示します。

なお、サンプルスクリプトの例として、 ``sample/01_transverse_field_ising`` の ``tutorial_example.py`` , ``tutorial_read.py`` があります。
あらかじめ ``tenes`` などにパスを通した上で

.. code::

   $ python tutorial_example.py

として実行できます(MacBook2017, 1.4 GHz Intel Core i7で数分程度で計算が全て終了します)。
得られた結果は

.. code::

   $ python tutorial_read.py

とすることで集計でき、 ``G``, エネルギー、 ``<Sz>`` 、 ``<Sx>`` が出力されます。


.. figure:: ../../img/tutorial_1_Sz_vs_G.*
   :width: 400px
   :align: center

図から ``G`` が大きくなるにつれ、 ``<Sz>`` が ``0.5`` から徐々に小さくなり最終的には0になることがわかります。
