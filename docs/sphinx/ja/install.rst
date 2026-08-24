
インストール方法
-------------------


ダウンロード
===================
TeNeS のソースコードは `GitHub page <https://github.com/issp-center-dev/TeNeS>`_ からダウンロードできます。gitがインストールされている環境で、以下のコマンドを打つとダウンロードが開始されます。

``$ git clone https://github.com/issp-center-dev/TeNeS``


必要なライブラリ・環境
======================
TeNeSをコンパイルするには以下のライブラリ・環境が必要です。

1. C++17 compiler

   - Intel コンパイラの場合、古い ``icpc`` では mptensor をビルドできないため、新しい ``icpx`` を利用してください。
     なお、 ``icpx`` はデフォルトで ``-fp-model=fast`` が有効になっており Inf/NaN の扱いが壊れるため、
     TeNeS のビルドでは自動的に ``-fp-model=precise`` を付加します。

2. CMake (>=3.8.0)
3. BLAS, LAPACK

TeNeSは以下のライブラリに依存していますが、自動でダウンロードおよびビルドがされます。

1. `mptensor <https://github.com/smorita/mptensor>`_ (>= v0.5.0)
2. `toml11 <https://github.com/ToruNiina/toml11>`_ (>= v4.0.0)

MPI および ScaLAPACK を利用することでテンソル演算を並列化できます。
MPI, ScaLAPACKについては自身でインストールする必要があります。
たとえば Debian GNU/Linux （やその派生ディストリビューション）をお使いで、root 権限をお持ちの場合は、

.. code::

   sudo apt install openmpi-bin libopenmpi-dev libscalapack-mpi-dev

でインストールすることが可能です。それ以外の方は、Open MPI などのMPI実装ならびに ScaLAPACKのホームページを参照の上、インストールをしてください。

また、TeNeS は任意の依存ライブラリとして `ARPACK-NG <https://github.com/opencollab/arpack-ng>`_ を利用できます。
インストールされていれば CMake が自動的に検出・リンクし、相関長計算における転送行列の
固有値ソルバとして使われます(入力ファイルの ``correlation_length`` セクションの
``eigensolver`` パラメータを参照)。なくてもビルド・実行に支障はありません
(組み込みのソルバが使われます)。
Debian GNU/Linux (やその派生ディストリビューション)では

.. code::

   sudo apt install libarpack2-dev

macOS (Homebrew) では

.. code::

   brew install arpack

でインストールすることが可能です。

入力ファイル作成ツールの使用には Python (>= 3.0.0) および
以下の Python パッケージが必要です。

1. numpy
2. scipy
3. toml

   
インストール
======================

1. TeNeS のディレクトリに移動したのち、以下の手順に従ってビルドを行います
（ CentOS など、環境によっては ``cmake3`` とする必要があります）。

::

  $ mkdir build
  $ cd build
  $ cmake -DCMAKE_INSTALL_PREFIX=<path to install to> ../
  $ make

``<path to install to>`` のデフォルト値は ``/usr/local`` です。

.. admonition:: 並列ビルド
  
  ``make`` コマンドに ``-j <num>`` オプションを追加し、 ``<num>`` 個のプロセスを用いた並列ビルドを行うと、 TeNeS を高速にビルド可能です。

なお、上記のコマンドで ``build/src`` ディレクトリに実行ファイル ``tenes`` が作成されます。

::

  $ make test

と打つとテストを実行することができます。


2. 次にインストールします。

::

  $ make install
 
実行ファイル ``tenes``, ``tenes_std``, ``tenes_simple`` が ``<path to install to>/bin`` にインストールされます。 


.. admonition:: MPI/ScaLAPACK 並列化の無効化

  MPI および ScaLAPACK を利用しない場合には、 ``-DENABLE_MPI=OFF`` オプションを ``cmake`` コマンドに追加してください。
  macOS では ScaLAPACK の一部関数とシステムのBLAS, LAPACK とで相性が悪く、エラー終了する場合があります。
  そのため macOS では ``ENABLE_MPI`` の既定値が ``OFF`` になっています。
  Homebrew で ``open-mpi`` と ``scalapack`` を導入した環境では
  ``-DENABLE_MPI=ON`` で動作することを確認しています
  (Homebrew の ScaLAPACK は Accelerate ではなく OpenBLAS にリンクされているため)。

.. admonition:: macOS で OpenMP を使う

  Apple clang には OpenMP ランタイムが同梱されていないため、別途 libomp が必要です。

  ::

    $ brew install libomp

  CMake は ``brew --prefix libomp`` の結果と Homebrew の既定の場所を順に探します。
  自動検出に失敗する場合は、libomp の場所を明示してください。

  ::

    $ cmake -DOpenMP_ROOT=$(brew --prefix libomp) ../

  なお Homebrew の ``mpicxx`` は Apple clang のラッパーなので、
  ``-DCMAKE_CXX_COMPILER=mpicxx`` を指定した場合もこの経路になります。
  Apple clang を使う場合は CMake 3.12 以降が必要です
  (``OpenMP_ROOT`` によるライブラリ検索が CMake 3.12 で追加されたため)。

.. admonition:: OpenMP スレッド数の指定

  macOS では OpenMP のバリア同期がシステムコールになるため、
  1 ノード内の細粒度な並列がかえって大きなオーバーヘッドになります。
  実行時には ``OMP_NUM_THREADS=1`` を設定し、並列化は MPI 側で行うことを推奨します。

  ::

    $ OMP_NUM_THREADS=1 mpiexec -np 4 tenes input.toml

.. admonition:: コンパイラの指定

   CMake では自動でコンパイラを検出してビルドを行います。コンパイラを指定したい場合には, 以下のようにオプションを追加してください。
   ::

      $ cmake -DCMAKE_CXX_COMPILER=<path to your compiler> ../


.. admonition:: ScaLAPACK ライブラリの指定

    CMake では自動で ScaLAPACK ライブラリを検出しますが、見つけられなかった場合などに、
    ``<path>/lib/libscalapack.so`` を利用したい場合には以下のようにオプションを追加してください。
  ::

    $ cmake -DSCALAPACK_ROOT=<path> ../


.. admonition:: ARPACK-NG の検出の制御

   CMake は自動で ARPACK-NG を検出します(見つからない場合は組み込みソルバのみでビルドされます)。
   検出の挙動は ``-DENABLE_ARPACK=AUTO/ON/OFF`` で制御できます
   (デフォルトは ``AUTO`` 。 ``ON`` では見つからない場合にエラーになり、
   ``OFF`` では検出を行いません)。
   特定の場所にインストールされた ARPACK-NG (``<path>/lib/libarpack.so``)を
   利用したい場合には以下のようにオプションを追加してください。
   ::

      $ cmake -DARPACK_ROOT=<path> ../


.. admonition:: mptensor の指定

   TeNeS は並列テンソル演算ライブラリ mptensor を利用しています。
   CMake は自動で mptensor をダウンロード・ビルドしますが、
   ユーザーが事前にインストールした mptensor (``<path>/lib/mptensor.a``)を使用したい場合には、以下のようにオプションを追加してください。
   ::

      $ cmake -DMPTENSOR_ROOT=<path> ../


.. admonition:: Python インタープリタの指定

   TeNeS ツールは Python3 で書かれており、 パスの通った ``python3`` コマンドを自動的に起動します。
   ツールの実行がうまく行かない場合には、 ``type python3`` などを利用して、 ``python3`` コマンドにパスが通っているかどうか確認してください。

   使うインタープリタを固定したい場合（あるいは ``/usr/bin/env`` コマンドが実行できずにエラーが出る場合）には、 以下のようにCMake オプションを追加してください。
   ::

      $ cmake -DTENES_PYTHON_EXECUTABLE=<path to your interpreter> ../
