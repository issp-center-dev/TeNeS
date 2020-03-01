.. highlight:: none

``tenes_std`` の使用方法
----------------------------

``tenes_std`` は ``tenes`` の入力ファイルを生成するツールです。

.. code:: bash

   $ ./tenes_std --help
   usage: tenes_std [-h] [-o OUTPUT] input

   Simple input generator for TeNeS

   positional arguments:
     input                 Input TOML file

   optional arguments:
     -h, --help            show this help message and exit
     -o OUTPUT, --output OUTPUT
                           Output TOML file


- 引数として入力ファイル名を取ります。
- コマンドラインオプションは、以下の通りです。
   - ``help``
      - ヘルプメッセージの表示
   - ``output``
      - 出力ファイル名
      - デフォルトは ``input.toml``
      - 入力ファイル名と同じファイル名にすることはできません
