.. highlight:: none

相関関数 ``C = <A(0)B(r)>`` を指定するセクション

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``r_max``,     "相関関数の距離 :math:`r` の最大値", 整数
   ``operators``, "相関関数を測る1体演算子 A, B の番号", 整数のリストのリスト

演算子は ``observable.onesite`` セクションで指定したものが用いられます。

例
~~

::

    [correlation]
    r_max = 5
    operators = [[0,0], [0,1], [1,1]]
