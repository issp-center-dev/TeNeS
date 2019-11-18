.. highlight:: none

シンプルモードの入力ファイル
---------------------------------

-  ファイルフォーマットは
   `TOML <https://qiita.com/minoritea/items/c0de47b8beb813c655d4>`__
   形式
-  ``model``, ``parameter``, ``lattice``, ``observable``,  ``correlation``
   の5つのセクションを持ちます。 なお、 ``parameter``, ``lattice``, ``observable`` , ``correlation`` セクションは、
   エキスパートモードの入力ファイルと共通ですので、そちらを参照してください。ここでは, ``model`` セクションについて解説します。

   -  将来的にはファイル名を指定することで分割可能にする予定です。


``model`` セクション
==========================

.. csv-table::
   :header: "name", "desc", "default"
   :widths: 15, 30, 10

   ``type``, "(spin, bose)", "spin"
   ``S`` , "スピンの大きさ", 0.5
   ``Jx`` , "Jxの大きさ", 1
   ``Jxy`` , "Jxyの大きさ", 1
   ``Jy`` , "Jyの大きさ", 0
   ``Jz`` , "Jz", 1
   ``G`` , "G", 0
