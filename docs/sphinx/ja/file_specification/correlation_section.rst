.. highlight:: none

サイト演算子の相関関数 :math:`C = \left\langle A(\boldsymbol{r}_0)B(\boldsymbol{r}_0+\boldsymbol{r})\right\rangle` に関する情報を指定するセクションです。
本セクションを省略した場合、相関関数は計算されません。

座標は正方格子TNS の座標系で測られます。すなわち、右隣のテンソルは :math:`\boldsymbol{r} = (1,0)` で、真上は :math:`\boldsymbol{r} = (0,1)` です。 
中心座標 :math:`\boldsymbol{r}_0` として、ユニットセル内のすべてのサイトが用いられます。
また、:math:`\boldsymbol{r}` は :math:`x` ないし :math:`y` 軸に平行な方向に、正の向きにのみ動きます。すなわち、

.. math::
   \boldsymbol{r} = (0,0), (1,0), (2,0), \dots, (r_\text{max}, 0), (0,1), (0,2), \dots, (0, r_\text{max})

です。

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

では相関関数 :math:`S^z(0)S^z(r), S^z(0)S^x(r), S^x(0)S^x(r)` が、 :math:`0 \le r \le 5` の範囲で測定されます。
