.. highlight:: none

相関長 :math:`\xi` の計算に関する情報を指定するセクションです。

.. csv-table::
   :header: "名前", "説明", "型", "デフォルト"
   :widths: 15, 30, 20, 20

   ``measure``,                  "相関長を測るかどうか",                               真偽値, true
   ``num_eigvals``,              "計算する転送行列の固有値の数",                       整数,   4
   ``maxdim_dense_eigensolver``, "密行列の対角化手法を用いる最大行列サイズ",           整数,   200
   ``arnoldi_maxdim``,           "Arnoldi 法で生成する Hessenberg 行列の次元",         整数,   0 (自動)
   ``arnoldi_restartdim``,       "Arnoldi 法のリスタートで生成する初期ベクトルの本数", 整数,   0 (自動)
   ``arnoldi_maxiterations``,    "Arnoldi 法の最大イテレーション回数",                 整数,   0 (自動)
   ``arnoldi_rtol``,             "Arnoldi 法で目指す相対残差",                         実数,   1e-10
   ``eigensolver``,              "転送行列の固有値ソルバ(auto / arpack / builtin)",   文字列, auto


相関長は転送行列の固有値から計算されます。
行列サイズが ``maxdim_dense_eigensolver`` 以下のときには密行列対角化(``?geev`` ルーチン)による対角化を、
そうでない場合は Implicit Restart Arnoldi (IRA)法による対角化を用いて固有値を計算します。

IRA 法では、 Arnoldi 過程によって大きさ ``arnoldi_maxdim`` のHessenberg 行列を生成し、その固有値を計算します。
収束していない場合は、新たに ``arnoldi_restartdim`` 本の初期ベクトルを作成し、 Arnodi 過程をやり直します (restart)。

``arnoldi_maxdim`` ・ ``arnoldi_restartdim`` ・ ``arnoldi_maxiterations`` に
0 以下の値を指定すると(デフォルト)、使用するソルバに応じた値が自動的に使われます。

- ARPACK-NG: ``arnoldi_maxdim`` は :math:`\max(2 \times \text{num\_eigvals} + 1, 25)` 、
  ``arnoldi_maxiterations`` は 10。ARPACK-NG は restart により小さい部分空間でも
  確実に収束できるため、小さめの部分空間と多めの restart 回数が効率的です。
- 組み込みソルバ: ``arnoldi_maxdim`` は :math:`\max(2 \times \text{num\_eigvals} + 1, 50)` 、
  ``arnoldi_maxiterations`` は 1。restart に頼らず 1 回の Arnoldi 過程で
  収束する大きさの部分空間を使います。
- ``arnoldi_restartdim`` は :math:`\max(\text{num\_eigvals} + 1, \text{arnoldi\_maxdim} / 2)` です。

収束に必要な部分空間の次元は固有値の分布(ギャップ)で決まり、行列サイズには
ほとんど依存しません。相関長が大きい(臨界点に近い)系では転送行列のスペクトルが
密になり収束が遅くなるため、 ``arnoldi_maxdim`` を増やすか、restart 回数
``arnoldi_maxiterations`` を増やしてください。
ARPACK-NG 使用時に部分空間が不足すると、 ``arnoldi_maxdim`` が自動(0 以下)の
場合は部分空間を 2 倍(上限は行列サイズ)にしながら収束するまで自動的に再試行します。
``arnoldi_maxdim`` を明示的に指定した場合は再試行せず、未収束の固有値は NaN として
出力され、相関長 :math:`\xi` の列も NaN になります
(組み込みソルバは未収束でも近似値をそのまま返します)。

行列サイズが ``maxdim_dense_eigensolver`` より大きい場合に使う反復固有値ソルバは
``eigensolver`` で選べます。

- ``"auto"`` (デフォルト): ARPACK-NG 付きでビルドされていれば ARPACK-NG を、
  そうでなければ組み込みの IRA 法を使います。
- ``"arpack"``: ARPACK-NG を使います。ARPACK-NG なしでビルドされたバイナリでは
  入力エラーになります(CMake オプション ``-DENABLE_ARPACK=ON`` でビルドしてください)。
- ``"builtin"``: 組み込みの IRA 法を使います。

ARPACK-NG を使う場合、 ``arnoldi_maxdim`` は Krylov 部分空間の次元 (ncv)、
``arnoldi_maxiterations`` は implicit restart の最大回数、 ``arnoldi_rtol`` は
Ritz 値の相対残差の許容値として引き継がれます。
``arnoldi_restartdim`` は組み込みソルバでのみ意味を持ちます。
