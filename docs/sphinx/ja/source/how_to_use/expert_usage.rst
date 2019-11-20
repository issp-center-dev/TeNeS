.. highlight:: none

``tenes`` の使用方法
------------------------------

``tenes`` の実行は以下のように行うことができます。

.. code:: bash

 $ tenes --help
 TeNeS: TEnsor NEtwork Solver for 2D quantum lattice system
 
   Usage:
     tenes [--quiet] <input_toml>
     tenes --help
     tenes --version
 
   Options:
     -h --help       Show this help message.
     -v --version    Show the version.
     -q --quiet      Do not print any messages.


-  引数として入力ファイル名を取ります。
-  コマンドラインオプションは、以下の通りです。

   -  ``help``
     - ヘルプメッセージの表示
   -  ``version``
     - バージョン情報の表示
   -  ``quiet``
     - 標準出力に何も書き出さないようにします

.. code:: bash

   tenes input.toml


