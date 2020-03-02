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
   10% done
   20% done
   30% done
   40% done
   50% done
   60% done
   70% done
   80% done
   90% done
   100% done
   Start calculating observables
     Start updating environment
     Start calculating local operators
       Save onesite observables to output/onesite_obs.dat
     Start calculating NN correlation
       Save twosite observables to output/twosite_obs.dat
       Save energy density and onesite observable densities to output/energy.dat
       Save elapsed times to output/time.dat

   Energy density = -0.5
   Onesite operator[0] density = 0.5
   Onesite operator[1] density = -5.19665527844e-92

   time simple update = 30.756653582
   time full update   = 0
   time environmnent  = 0.10026456
   time observable    = 0.030680172


のように計算が実行されます。
最初に並列化の情報およびテンソルの実虚が表示されます。
次に計算プロセスの実行状況が表示されます。
計算終了後、 ``Energy`` と局在演算子 ``Onesite operator [0]`` (``<Sz>``),   ``Onesite operator [1]`` (``<Sx>``)がそれぞれ出力されます。最後に ``time`` でどの程度計算時間がかかったか出力されます(単位は秒)。
計算終了後は ``output`` フォルダに
``energy.dat, parameters.dat, time.dat, onesite_obs.dat, twosite_obs.dat``
がそれぞれ出力されます。各出力ファイルの詳細は、ファイルフォーマットをご覧ください。
例えば ``<Sz>`` の値は、 ``onesite_obs.dat`` の2行目から読み取ることが可能です。
``G`` をパラメータとして0.1刻みで0-3.0まで振ったときの結果を下図に表示します。

なお、サンプルスクリプトの例として、 ``sample/01_transverse_field_ising`` の ``tutorial_example.py`` , ``tutorial_read.py`` があります。実行は、

.. code::

   $ python tutorial_example.py

でできます(MacBook2017, 1.4 GHz Intel Core i7で数分程度で計算が全て終了します)。
得られた結果は

.. code::

   $ python tutorial_read.py

とすることで, 標準出力に、G, エネルギー、 ``<Sz>`` 、 ``<Sx>`` が出力されます。


.. figure:: ../../img/tutorial_1_Sz_vs_G.*
   :width: 400px
   :align: center

図から ``G`` が大きくなるにつれ、 ``<Sz>`` が ``0.5`` から徐々に小さくなり最終的には0になることがわかります。
