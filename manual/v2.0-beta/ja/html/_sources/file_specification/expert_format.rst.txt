.. highlight:: none

.. _sec-expert-format:

``tenes`` の入力ファイル
---------------------------------

-  ファイルフォーマットは
   `TOML <https://qiita.com/minoritea/items/c0de47b8beb813c655d4>`__
   形式
-  ``parameter``, ``tensor``, ``evolution``, ``observable``, ``correlation``
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

simple update, full update で使う(虚)時間発展演算子を記述します。
1サイトおよび隣接2サイトに関する(虚)時間発展を定義できます。
次のようなフィールドを持つ ``simple``, ``full`` の2つのサブセクションを持ちます。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``group``,       "演算子のグループ",                        整数 (0-)
   ``site``,        "site の番号",                             整数 (0-)
   ``source_site``, "source site の番号",                      整数 (0-)
   ``source_leg``,  "source site から見た target site の方向", 整数 (0-3)
   ``dimensions``,  "虚時間発展演算子テンソルの次元",          整数のリスト
   ``elements``,    "虚時間発展演算子テンソルの非ゼロ要素",    文字列

``group`` は時間発展演算子のグループを指定します（省略した場合 0とみなされます）。
``parameter.simple_update`` および ``parameter.full_update`` における ``tau`` や ``num_steps`` で刻み幅やステップ数をリストを用いて複数指定したときに、そのインデックスに対応します。

``site`` は1サイト演算子に、 ``source_site`` と ``source_leg`` は2サイト演算子に使用します。

``source_leg`` は 0 から3までの整数で指定します。
-x 方向から順番に時計回りに、 ``0:-x, 1:+y, 2:+x, 3:-y`` として定義されています。

``dimensions`` は ``observable`` の ``dim`` と異なり、すべての足の次元を指定する必要があります。
足の順番は ``elements`` と同様に、 ``source_initial, target_initial, source_final, target_final`` の順番です。

例 ::

    [evolution]
  
    # One site
    [[evolution.simple]]
    site = 0
    dimensions = [2, 2]
    elements = """
    0 0  1.0012507815756226 0.0
    1 1  0.9987507809245809 0.0
    """
  
    # Two site
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


``correlation_length`` セクション
======================================

.. include:: ./correlation_length_section.rst
