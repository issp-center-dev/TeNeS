.. highlight:: none

.. _sec-std-format:

``tenes_std`` の入力ファイル
---------------------------------

-  ファイルフォーマットは
   `TOML <https://qiita.com/minoritea/items/c0de47b8beb813c655d4>`__ 形式
-  ``parameter``, ``tensor``, ``hamiltonian``, ``observable``, ``correlation``
   の5つのセクションを持ちます。

   - ``hamiltonian`` 以外の4つは以下に挙げる例外を除き、 ``tenes`` の入力ファイルフォーマットと同一であり、そのままtenes の入力ファイルとしてコピーされます。
  
   - ``parameter.simple_update.tau`` および ``parameter.full_update.tau`` に実数を渡すことで、虚時間発展演算子における虚時間刻みを指定できます。


``parameter`` セクション
===========================

.. include:: ./parameter_section.rst


``tensor`` セクション
===========================

.. include:: ./tensor_section.rst


``observable`` セクション
=============================

.. include:: ./observable_section.rst


``hamiltonian`` セクション
==============================

ハミルトニアン全体をボンドハミルトニアン(2サイトハミルトニアン) の和

.. math::
  \mathcal{H} = \sum_{i,j} \mathcal{H}_{ij}

であると捉えて、個々のボンドハミルトニアンを定義します。
定義のやりかたは ``observable.twosite`` でなされる2サイト演算子と同様です。

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
例として、2つの :math:`S=1/2` スピンの相互作用の場合は、 ``dim = [2,2]`` です。

``elements`` は演算子の非ゼロ要素を指定する文字列です。
1つの要素は4つの整数と2つの浮動小数点数を空白区切りからなる1つの行からなります。

- 最初の2つは演算子が作用する **前** の source site, target site の状態番号を示します。
- つぎの2つは演算子が作用した **後** の source site, target site の状態番号を示します。
- 最後の2つはそれぞれ演算子の要素の実部と虚部を示します。


``correlation`` セクション
===========================

.. include:: ./correlation_section.rst


``correlation_length`` セクション
====================================

.. include:: ./correlation_length_section.rst
