.. highlight:: none

.. _sec-input-summary:

TeNeS の入力ファイルの簡易まとめ
-----------------------------------

TeNeS の入力ファイルは `TOML <https://github.com/toml-lang/toml/blob/master/versions/ja/toml-v0.5.0.md>`__ 形式 で書かれており、
入力ファイルはいくつかのセクションに分かれています。
``tenes_simple`` と ``tenes_std`` は自分が必要とするセクションの情報を入力として読み取り、
それぞれ ``tenes_std`` と ``tenes`` の入力ファイルを生成します。
``tenes`` は入力ファイルの各セクションに書かれた情報を元に実際の計算を行います。

例えば ``tenes_simple`` は ``model`` と ``lattice`` の情報から ``tensor``, ``observable``, ``hamiltonian`` の情報を生成し、
さらに ``parameter``, ``correlation``, ``correlation_length`` はそのままコピーして、 ``tenes_std`` の入力ファイルとして出力します。

次表は各セクションの簡単な説明および各ツールがどう扱うかを示しています。

.. csv-table::
  :header: "セクション名", "説明", ``tenes_simple``, ``tenes_std``, ``tenes``
  :widths: 15, 15, 15, 15, 10

  ``parameter``,          "計算パラメータ",   "copy", "in / copy", "in"
  ``model``,              "模型パラメータ",   "in",   "",          ""
  ``lattice``,            "格子パラメータ",   "in",   "",          ""
  ``tensor``,             "テンソル",         "out",  "in / copy", "in"
  ``observable``,         "測定する演算子",   "out",  "copy",      "in"
  ``correlation``,        "相関関数",         "copy", "copy",      "in"
  ``correlation_length``, "相関長",           "copy", "copy",      "in"
  ``hamiltonian``,        "ハミルトニアン",   "out",  "in",        ""
  ``evolution``,          "虚時間発展演算子", "",     "out",       "in"

- "in"

  - ツールはこのセクションの情報を利用します

- "out"

  - ツールはこのセクションを新たに生成し、出力します

- "copy"

  - ツールはこのセクションを変更せずにそのまま出力します

