.. highlight:: none

エキスパートモードの使用方法
------------------------------

エキスパートモードの実行は以下のように行うことができます。

.. code:: bash

    $ tenes --help
    TeNeS: PEPS+CTM method for solving 2D quantum lattice system
    Usage: tenes [OPTIONS] input_toml

    Positionals:
      input_toml TEXT REQUIRED    Input TOML file

    Options:
      -h,--help                   Print this help message and exit
      -v,--version                Show version information

-  引数として入力ファイル名を取ります。
-  コマンドラインオプションは、以下の通りです。

   -  ``help``
   -  ``version``

.. code:: bash

    tnsolve input.toml

