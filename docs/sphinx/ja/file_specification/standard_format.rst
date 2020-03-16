.. highlight:: none

``tenes_std`` の入力ファイル
---------------------------------

-  ファイルフォーマットは
   `TOML <https://qiita.com/minoritea/items/c0de47b8beb813c655d4>`__ 形式
-  ``parameter``, ``tensor``, ``hamiltonian``, ``observable``, ``correlation``
   の5つのセクションを持ちます。

   - ``hamiltonian`` 以外の4つは以下に挙げる例外を除き、 ``tenes`` の入力ファイルフォーマットと同一であり、そのままtenes の入力ファイルとしてコピーされます。
  
       - ``parameter.simple_update.tau`` および ``parameter.full_update.tau`` に実数を渡すことで、虚時間発展演算子における虚時間刻みを指定できます。

       - ``observable.twosite`` について、 ``group=0`` の演算子はボンドハミルトニアンで上書きされます。
        
           - ``tenes_std`` に ``--no-hamiltonian-observe`` オプションを渡すことでこの動作を無効化できます。

``hamiltonian`` セクション
==========================

twosite ハミルトニアンを定義します。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``bonds``,    "ボンド",             文字列
   ``dim``,      "演算子の次元",       整数のリスト
   ``elements``, "演算子の非ゼロ要素", 文字列

``bonds`` は演算子が作用するサイト対の集合を表す文字列です。
3つの整数からなる1行が1つのサイト対を意味します。
最初の整数は 始点サイト (source) の番号です。
あとの2つの整数は source site からみた終点サイト (target) の座標 (dx, dy) です。

``dim`` は演算子の次元、すなわち作用するサイトの取りうる状態数です。

``elements`` は演算子の非ゼロ要素を指定する文字列です。
1つの要素は4つの整数と2つの浮動小数点数を空白区切りからなる1つの行からなります。
最初の2つは演算子が作用する前の source site, target site の状態番号を、
つぎの2つは演算子が作用した後の source site, target site の状態番号を、
最後の2つはそれぞれ演算子の要素の実部と虚部を示します。

例えば S=1/2 のハイゼンベルグハミルトニアンは次のように定義されます。::

  [[hamiltonian]]
  dim = [2, 2]
  bonds = """
  0 0 1
  0 1 0
  1 0 1
  1 1 0
  2 0 1
  2 1 0
  3 0 1
  3 1 0
  """
  elements = """
  0 0 0 0  0.25 0.0
  1 0 1 0  -0.25 0.0
  0 1 1 0  0.5 0.0
  1 0 0 1  0.5 0.0
  0 1 0 1  -0.25 0.0
  1 1 1 1  0.25 0.0
  """


``observable`` セクション
=============================

.. include:: ./observable_section.rst


- ``observable.twosite`` について、 ``group=0`` の演算子はボンドハミルトニアンで上書きされます。

   - ``tenes_std`` に ``--no-hamiltonian-observe`` オプションを渡すことでこの動作を無効化できます。



``tensor`` セクション
===========================

.. include:: ./tensor_section.rst


``correlation`` セクション

.. include:: ./correlation_section.rst
