
インストール方法
-------------------


ダウンロード
===================
TeNeS のソースコードは `GitHub page <https://github.com/issp-center-dev/TeNeS>`_ からダウンロードできます。gitがインストールされている環境で、以下のコマンドを打つとダウンロードが開始されます。

``$ git clone https://github.com/issp-center-dev/TeNeS``


必要なライブラリ・環境
======================
TeNeSのコンパイルするには以下のライブラリ・環境が必要です。

1. C++11 compiler
2. CMake (>=2.8.14)
3. MPI と ScaLAPACK

TeNeSは以下のライブラリに依存していますが、自動でダウンロードおよびビルドがされます。

1. `mptensor <https://github.com/smorita/mptensor>`_ 
2. `cpptoml <https://github.com/skystrife/cpptoml>`_
3. `sanitizers-cmake <https://github.com/arsenm/sanitizers-cmake>`_

インストール
======================

1. 以下の手順に従ってビルドを行います。

::

  $ mkdir build
  $ cd build
  $ cmake ..
  $ make

上記のコマンドで ``build/src`` ディレクトリに実行ファイル ``tests`` が作成されます。
  
2. 次にインストールを実行します。

::

  $ cmake -DCMAKE_INSTALL_PREFIX=<path to install to> ../
  $ make
  $ make install
 
上の例では、実行ファイル ``tenes`` が ``<path to install to>/bin`` にインストールされます。 ``<path to install to>`` のデフォルト値は ``/usr/local`` です。


.. admonition:: コンパイラの指定

   ``CMake`` では自動でコンパイラを検出してビルドを行います。コンパイラを指定したい場合には, 以下のコマンドを打ってください。
   ::

      $ cmake -DCMAKE_CXX_COMPILER=<path to your compiler> ../


.. admonition:: プリインストールされた ``mptensor`` の利用

   ``TeNeS`` でプリインストールされた ``mptensor`` を使用したい場合には、以下のコマンドを打ってください。
   ::

      $ cmake -DMPTENSOR_ROOT=<path to mptensor> ../
