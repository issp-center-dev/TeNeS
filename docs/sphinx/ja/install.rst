
インストール方法
-------------------


ダウンロード
===================
TeNeS のソースコードは `GitHub page <https://github.com/issp-center-dev/TeNeS>`_ からダウンロードできます。gitがインストールされている環境で、以下のコマンドを打つとダウンロードが開始されます。

``$ git clone https://github.com/issp-center-dev/TeNeS``


必要なライブラリ・環境
======================
TeNeSをコンパイルするには以下のライブラリ・環境が必要です。

1. C++11 compiler
2. CMake (>=3.6.0)

TeNeSは以下のライブラリに依存していますが、自動でダウンロードおよびビルドがされます。

1. `mptensor <https://github.com/smorita/mptensor>`_ 
2. `cpptoml <https://github.com/skystrife/cpptoml>`_
3. `sanitizers-cmake <https://github.com/arsenm/sanitizers-cmake>`_

MPI および ScaLAPACK を利用することでテンソル演算を並列化できます。
MPI, ScaLAPACKについては自身でインストールする必要があります。
たとえば Debian GNU/Linux （やその派生ディストリビューション）をお使いで、root 権限をお持ちの場合は、

.. code::

   sudo apt install openmpi-bin libopenmpi-dev libscalapack-mpi-dev

でインストールすることが可能です。それ以外の方は、Open MPI などのMPI実装ならびに ScaLAPACKのホームページを参照の上、インストールをしてください。

入力ファイル作成ツールの使用には Python (>= 3.0.0) および
以下の Python パッケージが必要です。

1. numpy
2. scipy
3. toml

   
インストール
======================

1. TeNeS のディレクトリに移動したのち、以下の手順に従ってビルドを行います。

::

  $ mkdir build
  $ cd build
  $ cmake -DCMAKE_INSTALL_PREFIX=<path to install to> ../
  $ make

上記のコマンドで ``build/src`` ディレクトリに実行ファイル ``tests`` が作成されます。
``<path to install to>`` のデフォルト値は ``/usr/local`` です。

（ CentOS など、環境によっては ``cmake3`` とする必要があります。）

2. 次にインストールします。

::
  $ make install
 
実行ファイル ``tenes``, ``tenes_std``, ``tenes_simple`` が ``<path to install to>/bin`` にインストールされます。 


.. admonition:: MPI/ScaLAPACK 並列化の無効化
  
  MPI および ScaLAPACK を利用しない場合には、 ``-DENABLE_MPI=OFF`` オプションを ``cmake`` コマンドに追加してください。
  macOS では ScaLAPACK の一部関数とシステムのBLAS, LAPACK とで相性が悪く、エラー終了するのを確認しています。
  MPI 並列の無効化を推奨しています。

.. admonition:: コンパイラの指定

   CMake では自動でコンパイラを検出してビルドを行います。コンパイラを指定したい場合には, 以下のようにオプションを追加してください。
   ::

      $ cmake -DCMAKE_CXX_COMPILER=<path to your compiler> ../


.. admonition:: Python インタープリタの指定

   CMake では自動でPython インタープリタを検出してビルドを行います。指定したい場合には, 以下のようにオプションを追加してください。
   ::

      $ cmake -DTENES_PYTHON_EXECUTABLE=<path to your interpreter> ../


.. admonition:: mptensor の指定

   TeNeS は並列テンソル演算ライブラリ mptensor を利用しています。
   CMake は自動で mptensor をダウンロード・ビルドしますが、
   ユーザーが事前にインストールした mptensor を使用したい場合には、以下のようにオプションを追加してください。
   ::

      $ cmake -DMPTENSOR_ROOT=<path to mptensor> ../
