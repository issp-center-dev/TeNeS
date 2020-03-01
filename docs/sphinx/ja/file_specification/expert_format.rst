.. highlight:: none

``tenes`` の入力ファイル
---------------------------------

-  ファイルフォーマットは
   `TOML <https://qiita.com/minoritea/items/c0de47b8beb813c655d4>`__
   形式
-  ``parameter``, ``lattice``, ``evolution``, ``observable``, ``correlation``
   の5つのセクションを持ちます。

``parameter`` セクション
========================

.. include:: ./parameter_section.rst


``tensor`` セクション
========================

.. include:: ./tensor_section.rst


``observable`` セクション
==========================

.. include:: ./observable_section.rst

``evolution`` セクション
========================

simple update, full update で使う虚時間発展演算子を記述します。
次のようなフィールドを持つ ``simple``, ``full`` の2つのサブセクションを持ちます。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``source_site``, "source site の番号",                      整数
   ``source_leg``,  "source site から見た target site の方向", 整数
   ``dimensions``,  "虚時間発展演算子テンソルの次元",          整数のリスト
   ``elements``,    "虚時間発展演算子テンソルの非ゼロ要素",    文字列


例 ::

  [evolution]
  [[evolution.simple]]
  source_site = 0
  source_leg = 2
  dimensions = [2, 2, 2, 2]
  elements = """
  0 0 0 0  0.9975031223974601 0.0
  1 0 1 0  1.0025156589209967 0.0
  0 1 1 0  -0.005012536523536871 0.0
  1 0 0 1  -0.005012536523536871 0.0
  0 1 0 1  1.0025156589209967 0.0
  1 1 1 1  0.9975031223974601 0.0
  """


``correlation`` セクション
==========================

.. include:: ./correlation_section.rst
