.. highlight:: none

横磁場イジングモデル
----------------------------

ここでは横磁場イジングモデルに対して、横磁場を変化させた場合の計算例について紹介します。
入力ファイルの変数 ``G`` を用いることで横磁場の大きさを調整することが可能です。例えば、横磁場が0の場合には、

.. code::

   [parameter]
   [parameter.tensor]
   D  = 2     # tensor_dim
   CHI  = 10  # env_dim

   [parameter.simple_update]
   num_step = 1000
   tau = 0.01

   [parameter.ctm]
   iteration_max = 10

   [lattice]
   type = "square lattice"
   L_sub = [2,2]

   [model]
   type = "spin"
   Jz = -1.0
   Jx = 0.0
   Jy = 0.0
   G = 0.0

とします(``Jz = -1.0`` なので、 ``G=0`` では強磁性状態になります)。入力ファイルを ``simple.toml`` とした場合、
   
.. code:: bash

   $ tenes_simple simple.toml
   $ tenes input.toml

を実行することで計算が開始されます。
計算を実行すると、

.. code:: bash

	  Start simple update
	  Start calculating observables
	  Start updating environment
	  Start calculating local operators
	  Save site observables to output/site_obs.dat
	  Start calculating energy
	  Save energy to output/energy.dat
	  Start calculating NN correlation
	  Save NN correlation to output/neighbor_obs.dat
	  Save elapsed times to output/time.dat

	  Energy = -0.5
	  Local operator 0 = 0.5
	  Local operator 1 = 1.90794709356e-11

	  time simple update = 3.21127
	  time full update   = 0
	  time environmnent  = 0.875561
	  time observable    = 0.132412
	  
のように計算が実行されます。
最初に各プロセスの実行状況が表示されます。
計算終了後、 ``Energy`` と局在演算子 ``Local operator 0`` (``<Sz>``),   ``Local operator 1`` (``<Sx>``)がそれぞれ出力されます。最後に ``time`` でどの程度計算時間がかかったか出力されます(単位は秒)。
計算終了後は ``output`` フォルダに
``energy.dat, parameters.dat, time.dat, neighbor_obs.dat, site_obs.dat``
がそれぞれ出力されます。各出力ファイルの詳細は、ファイルフォーマットをご覧ください。
``<Sz>`` の値は、 ``site_obs.dat`` の ``0 0`` 成分もしくは標準出力の ``Local operator 0`` に続く値から抽出することが可能で、 
``G`` をパラメータとして0.1刻みで0-2.0まで振り、得られた結果を下図に表示します。

なお、サンプルスクリプトの例として、 ``sample/01_transverse_field_ising`` の ``tutorial_example.py`` , ``tutorial_read.py`` があります。実行は、

.. code::

   $ python tutorial_example.py

でできます(MacBook2017, 1.4 GHz Intel Core i7で数分程度で計算が全て終了します)。
得られた結果は

.. code::

   $ python tutorial_read.py

とすることで, 標準出力に、G, エネルギー、 ``<Sz>`` 、 ``<Sx>`` が出力されます。


.. image:: ../../img/tutorial_1_Sz_vs_G.pdf
   :width: 400px
   :align: center

図から ``G`` が大きくなるにつれ、 ``<Sz>`` が ``0.5`` から徐々に小さくなり最終的には0になることがわかります。
