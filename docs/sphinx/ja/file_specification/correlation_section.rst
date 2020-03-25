.. highlight:: none

相関関数 ``C = <A(0)B(r)>`` に関する情報を指定するセクションです。
TeNeS は x ないし y 軸に平行な方向にのみ相関関数を計算します。
つまり、 :math:`\langle A(x,y) B(x+r, y) \rangle` という形式のみ計算し、
斜め方向 :math:`\langle A(x,y) B(x+r_x, y+r_y)` は計算しません。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``r_max``,     "相関関数の距離 :math:`r` の最大値", 整数
   ``operators``, "相関関数を測る1体演算子 A, B を表す番号", 整数のリストのリスト

演算子は ``observable.onesite`` セクションで指定したものが用いられます。

例
~~

例えば :math:`S^z` が0 番で、 :math:`S^x` が1 番として定義されている場合、

::

    [correlation]
    r_max = 5
    operators = [[0,0], [0,1], [1,1]]

では :math:`S^z(0)S^z(r), S^z(0)S^x(r), S^x(0)S^x(r)` が測定されます。
