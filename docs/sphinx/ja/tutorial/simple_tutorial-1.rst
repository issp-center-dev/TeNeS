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
   num_step = 1000
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
（あらかじめ TeNeS をインストールしたのち、環境変数 PATH を適切に設定してください。）
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
	  Save onesite observables to output_0/onesite_obs.dat
	  Start calculating twosite operators
	  Save twosite observables to output_0/twosite_obs.dat
	  Save observable densities to output_0/density.dat
	  Save elapsed times to output_0/time.dat

	  Onesite observables per site:
	  Sz          = 0.5 0
	  Sx          = -1.28526262482e-13 0
	  Twosite observables per site:
	  hamiltonian = -0.5 0
	  SzSz        = 0.5 0
	  SxSx        = -1.7374919982e-18 0
	  SySy        = 1.73749202733e-18 0
	  Wall times [sec.]:
	  simple update = 3.545813509
	  full update   = 0
	  environmnent  = 0.123170523
	  observable    = 0.048149856

	  Done.
	  
のように計算が実行されます。
最初に並列化の情報およびテンソルの実虚が表示されます。
次に計算プロセスの実行状況が表示されます。
計算終了後、 1サイト演算子 ``Sz``,   ``Sx`` およびハミルトニアン ``hamiltonian`` , 最近接相関 ``SzSz``, ``SxSx``, ``SySy`` のサイトあたりの期待値が出力されます。
最後にフェーズごとの計算時間が出力されます(単位は秒)。
計算終了後は ``output`` ディレクトリに
``density.dat, parameters.dat, time.dat, onesite_obs.dat, twosite_obs.dat``
がそれぞれ出力されます。各出力ファイルの詳細は、 :ref:`sec-output-format` をご覧ください。
例えば ``<Sz>`` の値は、 ``onesite_obs.dat`` から読み取ることが可能です。
``G`` をパラメータとして0.2刻みで0-3.0まで振ったときの結果を下図に表示します。

なお、サンプルスクリプトの例として、 ``sample/01_transverse_field_ising`` の ``tutorial_example.py`` , ``tutorial_read.py`` があります。
あらかじめ ``tenes`` などにパスを通した上で

.. code::

   $ python tutorial_example.py

として実行できます(MacBook2017, 1.4 GHz Intel Core i7で数分程度で計算が全て終了します)。
得られた結果は

.. code::

   $ python tutorial_read.py

とすることで集計でき、 ``G``, エネルギー、 ``<Sz>`` 、 ``<Sx>`` が出力されます。


.. figure:: ../../img/tutorial_1_Sz_vs_G.pdf
     :name: fig_transverse
     :width: 400px
     :align: center

     ``<Sz>`` , ``<Sx>`` の ``G`` 依存性


:numref:`fig_transverse` から ``G`` が大きくなるにつれ、 ``<Sz>`` が ``0.5`` から徐々に小さくなり最終的には0になる一方、 ``<Sx>`` は ``0`` から大きくなり最終的には ``0.5`` になることが分かります。

三角格子・正方格子ハイゼンベルク模型の磁化過程
-----------------------------------------------

次に三角格子上で定義されたスピン\ :math:`S=1/2`\ の量子ハイゼンベルク模型の磁化過程の計算を紹介します。
ハミルトニアンは以下のようになります:

.. math::

   \begin{aligned}
   H = J \sum_{\langle i,j \rangle}\sum_{\alpha}^{x,y,z} {S}_i^{\alpha} {S}_j^{\alpha} - \sum_i h S_i^z\end{aligned}

ここで\ :math:`\langle i, j\rangle`\ は隣接格子の組を表し、\ :math:`h`\ は\ :math:`z`\ 方向にかけられた外部磁場の大きさを表します。
この模型の基底状態を計算し、ユニットセルの平均磁化 \ :math:`\langle S_z \rangle\equiv \frac{1}{N_u}\sum_i^{N_u} \langle S_i^z \rangle`\ を磁場\ :math:`h`\ の関数として求めてみましょう ( :math:`N_u` はユニットセル内のサイト数)。

この計算を行うには、 ``sample/05_magnetization`` のディレクトリ中にある toml ファイル ``basic.toml`` と、pythonスクリプト ``tutorial_magnetization.py`` を利用します。 ``basic.toml`` ファイルには、模型の設定やパラメータなどが記述されています。

::

    [parameter]
    [parameter.general]
    is_real = true

    [parameter.simple_update]
    num_step = 200
    tau = 0.01

    [parameter.full_update]
    num_step = 0
    tau = 0.01

    [parameter.ctm]
    iteration_max = 10
    dimension = 10

    [lattice]
    type = "triangular lattice"
    L = 3
    W = 3
    virtual_dim = 4
    initial = "random"

    [model]
    type = "spin"
    J = 1.0

``lattice`` セクションで三角格子を指定しており、ユニットセルの大きさは\ :math:`3\times 3`\ を指定しています。
ここでは計算を軽くするために、 ``simple update`` だけを行っており、虚時間の刻み幅\ :math:`\tau`\ は\ :math:`\tau = 0.01`\ としています。また簡単のため、\ :math:`J=1`\ としています。この基本設定ファイルを用いて、 ``tutorial_magnetization.py`` では磁場を掃引したときの磁化を計算します。

::

    import subprocess
    from os.path import join
    import numpy as np
    import toml

    num_h = 21
    min_h = 0.0
    max_h = 5.0
    num_step_table = [100, 200, 500, 1000, 2000]

    fout = open("magnetization.dat","w")
    for idx, h in enumerate(np.linspace(min_h, max_h, num=num_h)):
        print("Caclulation Process: {}/{}".format(idx+1, num_h))
        inum = 0
        num_pre = 0
        fout.write("{} ".format(h))
        for num_step in num_step_table:
            ns = num_step - num_pre
            print("Step numter: {}".format(num_step))
            with open("basic.toml") as f:
                dict_toml = toml.load(f)
            dict_toml["parameter"]["general"]["output"] = "output_{}_{}".format(idx,num_step)
            dict_toml["parameter"]["general"]["tensor_save"] = "tensor_save_{}_{}".format(idx,num_step)
            dict_toml["model"]["H"] = float(h)
            dict_toml["parameter"]["simple_update"]["num_step"] = ns
            if inum > 0:
                dict_toml["parameter"]["general"]["tensor_load"] = "tensor_save_{}_{}".format(idx,num_pre)
            with open("simple_{}_{}.toml".format(idx,num_step), 'w') as f:
                toml.dump(dict_toml, f)
            cmd = "tenes_simple simple_{}_{}.toml -o std_{}_{}.toml".format(idx,num_step,idx,num_step)
            subprocess.call(cmd.split())
            cmd = "tenes_std std_{}_{}.toml -o input_{}_{}.toml".format(idx,num_step,idx,num_step)
            subprocess.call(cmd.split())
            cmd = "tenes input_{}_{}.toml".format(idx,num_step)
            subprocess.call(cmd.split())
            with open(join("output_{}_{}".format(idx,num_step), "density.dat")) as f:
                lines = f.readlines()
                mag_sz = lines[0].split('=')[1].strip()
            fout.write("{} ".format(mag_sz))
            inum = inum + 1
            num_pre = num_step
        fout.write("\n")
    fout.close()

ここスクリプトでは、磁場\ :math:`h`\ を0から5まで0.25刻みで変化させ、基底状態のエネルギーと\ :math:`\langle S_z \rangle`\ を計算して、 ``energy.dat`` および ``magnetization.dat`` に出力します。 ``simple update`` の時間ステップ数を\ :math:`100`, :math:`200`, :math:`500`, :math:`1000`, :math:`2000`\ と変化させたときの様子を見るために、各磁場でステップ数を変えた計算も行っています。
計算量を減らすために、少ないステップ数で得られた波動関数の情報を ``tensor_save`` に保存し、それをより多いステップ数の計算の初期状態としてとっています。例えば、最初に時間ステップ数を100とした計算を行って結果を出力したあと、ステップ数100の計算終了時の波動関数からさらにステップ数100の計算を行って、結果的にステップ数200の計算結果を得ます。

実際に実行してみましょう。あらかじめ ``tenes`` などにパスを通した上で

::

    python tutorial_magnetization.py

により計算を実行します。ノートPC(シングルプロセッサ)では数時間程度の計算量となります。
計算が終了したら、gnuplotを起動し、

::

    load 'plot.gp'

とすれば、:numref:`fig_triangular` の右図のような磁化カーブが得られます。同様に

::

    load 'plot_ene.gp'

とすれば、:numref:`fig_triangular` の左図のような基底エネルギーのグラフが得られます。

十分なステップ数(例えばステップ数2000)の計算結果からわかるように、磁化過程には飽和磁化\ :math:`\langle S_z \rangle=0.5`\ の\ :math:`1/3`\ の磁化のところで、プラトー構造が生じます。
このプラトー上では、3つの格子上のスピンが\ :math:`\uparrow`,
:math:`\uparrow`,
:math:`\downarrow`\ と磁化した周期構造を形成し、スピンギャップが生じています。このプラトー構造は三角格子特有のものです。
実際に計算精度がでているかどうかをみるには、エネルギーのステップ依存性が参考になります。
理想的にはステップ数を増やすほど基底エネルギーが下がるはずですが、一部の磁場領域では逆に基底エネルギーが増加します。これは計算精度があまりでていない兆候です。
ボンド次元を増やすなどして、より計算精度を高める必要があると推測されます。

.. figure:: ../../img/Fig_Triangular.pdf
     :name: fig_triangular
     :width: 800px

     三角格子量子ハイゼンベルク模型のエネルギー(左図)と磁化過程(右図)

では正方格子でも同じことをやってみましょう。 ``sample/05_magnetization`` のディレクトリ中にある toml ファイル ``basic_square.toml`` と、pythonスクリプト ``tutorial_magnetization_square.py`` を利用します。
``basic_square.toml`` は、 ``lattice`` セクションが以下のように変更されているほかは ``basic.toml`` と同じ内容です。

::

    [lattice]
    type = "square lattice"
    L = 2
    W = 2
    \begin{lstlisting}

    実際に計算を行うには、
    \begin{lstlisting}
    python tutorial_magnetization.py

とします。計算が終了したら、gnuplotを起動し、

::

    load 'plot_square.gp'

とすれば、 :numref:`fig_square` の右図のような磁化カーブが得られます。同様に

::

    load 'plot_ene_square.gp'

とすれば、:numref:`fig_square` の左図のような基底エネルギーのグラフが得られます。

.. figure:: ../../img/Fig_square.pdf
     :name: fig_square
     :width: 800px

     正方格子量子ハイゼンベルク模型のエネルギー(左図)と磁化過程(右図)

ステップ数2000でほぼ収束しており、三角格子ハイゼンベルク模型と異なり、プラトー構造は現れないことがわかります。
エネルギーは概ね、ステップ数を増加させると減少するため、ある程度計算精度がでていると推測されます。


