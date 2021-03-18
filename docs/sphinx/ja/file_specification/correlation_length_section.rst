.. highlight:: none

相関長 :math:`\xi` の計算に関する情報を指定するセクションです。

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 15, 30, 20, 20

   ``measure``,                  "相関長を測るかどうか",                               真偽値, false
   ``maxdim_dense_eigensolver``, "密行列の対角化手法を用いる最大行列サイズ",           整数,   200
   ``arnoldi_maxdim``,           "Arnoldi 法で生成する Hessenberg 行列の次元",         整数,   100
   ``arnoldi_restartdim``,       "Arnoldi 法のリスタートで生成する初期ベクトルの本数", 整数,   100
   ``arnoldi_maxiterations``,    "Arnoldi 法の最大イテレーション回数",                 整数,   1
   ``arnoldi_rtol``,             "Arnoldi 法で目指す相対残差",                         実数,   1e-10

相関長は転送行列の固有値から計算されます。
行列サイズが ``maxdim_dense_eigensolver`` 以下のときには密行列対角化(``?geev`` ルーチン)による対角化を、
そうでない場合は Implicit Restart Arnoldi (IRA)法による対角化を用いて固有値を計算します。

IRA 法では、 Arnoldi 過程によって大きさ ``arnoldi_maxdim`` のHessenberg 行列を生成し、その固有値を計算します。
収束していない場合は、新たに ``arnoldi_restartdim`` 本の初期ベクトルを作成し、 Arnodi 過程をやり直します (restart)。
転送行列の場合は、多くの場合で restart をする必要はありません (``arnoldi_maxiterations = 1``)。
