.. highlight:: none

``tenes_std`` の使用方法
----------------------------

``tenes_std`` は与えられたハミルトニアンから虚時間発展演算子を導出し、 ``tenes`` の入力ファイルを生成するツールです。

.. code:: bash

   $ tenes_std std.toml


- 引数としてファイルを取ります
- ``tenes`` の入力ファイルを出力します
- コマンドラインオプションは以下の通りです
   - ``--help``
      - ヘルプメッセージの表示
   - ``--version``
      - バージョン番号の表示
   - ``--output=filename``
      - 出力するファイルの名前 ``filename`` を指定します
      - デフォルトは ``input.toml``
      - 入力ファイル名と同じファイル名にすることはできません

入力ファイルを生成・編集することで、定義されていない模型・格子での計算が行なえます。
入力ファイルの詳細は :ref:`sec-std-format` を参照してください。
以下、正方格子上で定義されたスピン1/2のハイゼンベルグ模型の入力ファイル例です。

::
   [parameter]
   [parameter.general]
   is_real = true   # limit tensors as real-valued ones
   [parameter.simple_update]
   num_step = 1000  # number of steps
   tau = 0.01       # imaginary time step
   [parameter.full_update]
   num_step = 0     # number of steps
   tau = 0.01       # imaginary time step
   [parameter.ctm]
   dimension = 9    # bond dimension

   [tensor]
   type = "square lattice"
   L_sub = [2, 2]   # unitcell size
   skew = 0         # boundary condition

   # tensors in unitcell
   [[tensor.unitcell]]
   index = [0, 3]   # index of tensors
   physical_dim = 2 # physical bond dimension
   virtual_dim = [3, 3, 3, 3]
                    # virtual bond dimension
   noise = 0.01     # noise in initial tensor
   initial_state = [1.0, 0.0]
                    # initial state

   [[tensor.unitcell]]
   index = [1, 2]
   physical_dim = 2
   virtual_dim = [3, 3, 3, 3]
   noise = 0.01
   initial_state = [0.0, 1.0]


   # (bond) hamiltonian
   [[hamiltonian]]
   dim = [2, 2]    # physical bond dimensions
   bonds = """     # bond information
   0 1 0           # first: index of one site
   1 1 0           # second: x coord of the other
   2 1 0           # third:  y coord of the other
   3 1 0
   0 0 1
   1 0 1
   2 0 1
   3 0 1
   """
   elements = """     # nonzero elements of tensor
   0 0 0 0 0.25 0.0   # first:  initial state of one site
   1 0 1 0 -0.25 0.0  # second: initial state of the other
   0 1 1 0 0.5 0.0    # third:  final state of one site
   1 0 0 1 0.5 0.0    # fourth: final state of the other 
   0 1 0 1 -0.25 0.0  # fifth:  real part
   1 1 1 1 0.25 0.0   # sixth:  imag part
   """

   # observables
   [observable]
   [[observable.onesite]]
   name = "Sz"        # name
   group = 0          # index
   sites = []         # sites to be acted
   dim = 2            # dimension
   elements = """     # nonzero elements
   0 0 0.5 0.0
   1 1 -0.5 0.0
   """

   [[observable.twosite]]
   name = "hamiltonian"
   group = 0
   dim = [2, 2]
   bonds = """
   0 1 0
   1 1 0
   2 1 0
   3 1 0
   0 0 1
   1 0 1
   2 0 1
   3 0 1
   """
   elements = """
   0 0 0 0 0.25 0.0
   1 0 1 0 -0.25 0.0
   0 1 1 0 0.5 0.0
   1 0 0 1 0.5 0.0
   0 1 0 1 -0.25 0.0
   1 1 1 1 0.25 0.0
   """

   [[observable.twosite]]
   name = "SzSz"
   group = 1
   dim = [2, 2]
   bonds = """
   0 1 0
   1 1 0
   2 1 0
   3 1 0
   0 0 1
   1 0 1
   2 0 1
   3 0 1
   """
   ops = [0, 0]  # index of onesite operators

