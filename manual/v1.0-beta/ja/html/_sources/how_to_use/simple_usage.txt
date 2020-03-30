.. highlight:: none

``tenes_simple`` の使用方法
----------------------------

``tenes_simple`` は定義済みの模型、格子に対する ``tenes_std`` の入力ファイルを生成するツールです。

.. code:: bash

   $ tenes_simple simple.toml

- 引数としてファイルを取ります
- ``tenes_std`` の入力ファイルを出力します
- コマンドラインオプションは以下の通りです
   - ``--help``
      - ヘルプメッセージの表示
   - ``--version``
      - バージョン番号の表示
   - ``--output=filename``
      - 出力するファイルの名前 ``filename`` を指定します
      - デフォルトは ``std.toml``
      - 入力ファイル名と同じファイル名にすることはできません

現在定義されている模型・格子は次の通り。

- 模型
   - スピン系
- 格子
   - 正方格子
   - 三角格子
   - 蜂の巣格子

模型・格子や入力ファイルの詳細は :ref:`sec-simple-format` を参照してください。
::
   [lattice]
   type = "square lattice" # type of lattice
   L = 2                   # size of unitcell
   W = 2                   # size of unitcell
   virtual_dim = 3         # bond dimension
   initial = "antiferro"   # initial state

   [model]
   type = "spin" # type of model
   J = 1.0       # Heisenberg interaction

   [parameter]
   [parameter.general]
   is_real = true # use real tensor

   [parameter.simple_update]
   num_step = 1000  # number of steps
   tau = 0.01       # imaginary time step

   [parameter.full_update]
   num_step = 0    # number of steps
   tau = 0.01      # imaginary time step

   [parameter.ctm]
   dimension = 9       # bond dimension

