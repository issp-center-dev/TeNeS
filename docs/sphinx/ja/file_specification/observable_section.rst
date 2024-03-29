.. highlight:: none


物理量測定に関する諸々を記述します。
``onesite``, ``twosite`` と ``multisite`` の3種類のサブセクションを持ちます。


``observable.onesite``
~~~~~~~~~~~~~~~~~~~~~~~~~

ひとつのサイト上で定義される物理量を示す一体演算子を定義します。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``name``,     "演算子の名前",               文字列
   ``group``,    "演算子の識別番号",           整数
   ``sites``,    "サイト番号",                 整数 or 整数のリスト
   ``dim``,      "演算子の次元",               整数
   ``elements``, "演算子の非ゼロ要素",         文字列
   ``coeff``,    "演算子にかかる係数（実部）", 実数
   ``coeff_im``, "演算子にかかる係数（虚部）", 実数

``name`` は演算子の名前です。

``group`` はonesite 演算子の識別番号です。

``sites`` は演算子が作用するサイト番号です。
リストを渡すことで複数同時に定義できます。
空リスト ``[]`` は全サイトを意味します。

``dim`` は演算子の次元です。

``elements`` は演算子の非ゼロ要素を指定する文字列です。
1つの要素は、空白で区切られた2つの整数と2つの浮動小数点数からなる1つの行で表されます。

- 最初の2つはそれぞれ演算子が作用する前と後の状態番号を示します。
- あとの2つはそれぞれ演算子の要素の実部と虚部を示します。

``coeff``, ``coeff_im`` は演算子にかかる係数の実部と虚部です。
省略した場合はそれぞれ 1.0, 0.0 になります。

例
....

S=1/2 のSz 演算子

.. math::
  S^z = \left(\begin{array}{cc} 0.5 & 0.0 \\ 0.0 & -0.5 \end{array}\right)

を具体例として説明します。

まず、名前は ``name = "Sz"`` として、識別番号は ``group = 0`` としておきます。

次に、演算子の作用するサイトですが、すべてのサイトで同一の演算子を用いる場合には
``sites = []`` とします。
そうではない場合、例えばスピンの大きさが異なるサイトがある場合には、
``sites = [0,1]`` などと具体的なサイト番号を指定します。

演算子の次元は、上に示した行列表示のサイズなので、 ``dim = 2`` です。

最後に演算子の要素です。
非ゼロ要素について、そのインデックス（ゼロ始まり）と要素を順番に並べれば良いので、
::

  elements = """
  0 0   0.5 0.0
  1 1  -0.5 0.0
  """

となります。

結果として、 S=1/2 の Sz 演算子は次のように定義されます。
::

  [[observable.onesite]]
  name = "Sz"
  group = 0
  sites = []
  dim = 2
  elements = """
  0 0  0.5 0.0
  1 1  -0.5 0.0
  """


``observable.twosite``
~~~~~~~~~~~~~~~~~~~~~~~~~

ふたつのサイト上で定義される物理量を示す演算子を定義します。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``name``,     "演算子の名前",             文字列
   ``group``,    "演算子の識別番号",         整数
   ``bonds``,    "ボンド",                   文字列
   ``dim``,      "演算子の次元",             整数のリスト
   ``elements``, "演算子の非ゼロ要素",       文字列
   ``ops``,      "onesite 演算子の識別番号", 整数のリスト
   ``coeff``,    "演算子にかかる係数（実部）", 実数
   ``coeff_im``, "演算子にかかる係数（虚部）", 実数

``name`` は演算子の名前です。

``group`` は twosites 演算子の識別番号です。

``bonds`` は演算子が作用するサイト対の集合を表す文字列です。
3つの整数からなる1行が1つのサイト対を意味します。

- 最初の整数は 始点サイト (source) の番号です。
- あとの2つの整数は source site からみた終点サイト (target) の座標 (dx, dy) です。

  - dx, dy ともに :math:`-3 \le dx \le 3` の範囲に収まる必要があります。 

``dim`` は演算子の次元、すなわち作用するサイトの取りうる状態数です。
例として、2つの :math:`S=1/2` スピンの相互作用の場合は、 ``dim = [2,2]`` です。

``elements`` は演算子の非ゼロ要素を指定する文字列です。
1つの要素は4つの整数と2つの浮動小数点数を空白区切りからなる1つの行からなります。

- 最初の2つは演算子が作用する **前** の source site, target site の状態番号を示します。
- つぎの2つは演算子が作用した **後** の source site, target site の状態番号を示します。
- 最後の2つはそれぞれ演算子の要素の実部と虚部を示します。

``ops`` を使うと ``observable.onesite`` で定義した1体演算子の直積として2体演算子を定義できます。
例えば ``observable.onesite`` の ``group=0`` として :math:`S^z` を定義していた場合には、
``ops = [0,0]`` として :math:`S^z_iS^z_j` を表現できます。

``elements`` と ``ops`` を同時に定義した場合にはエラー終了します。

``coeff``, ``coeff_im`` は演算子にかかる係数の実部と虚部です。
省略した場合はそれぞれ 1.0, 0.0 になります。

例
....
ここでは具体例として、``Lsub=[2,2]`` の正方格子 S=1/2 ハイゼンベルグ模型のボンドハミルトニアンのエネルギーを求めるため、
ハミルトニアン

.. math::
  \mathcal{H}_{ij} = S_i^z S_j^z + \frac{1}{2} \left[S_i^+ S_j^- + S_i^- S_j^+ \right]

を2体演算子として設定する例を説明します。

まず、名前と識別番号はそれぞれ ``name = "hamiltonian"`` と ``group = 0`` としておきます。
それぞれのサイトの状態は :math:`|\uparrow\rangle` と :math:`|\downarrow\rangle` の2状態の重ね合わせとなるため、次元は 2 となり、
``dim = [2,2]`` となります。

次にボンドです。サイトは :numref:`bond_22` のように並んでいます。
0 番と 1 番をつなぐボンドは、 1番は 0 番から見て (1,0) の位置にあるので ``0 1 0`` と表現されます。
同様に 1 番と 3 番をつなぐボンドは、 3 番が 1 番から見て (0,1) の位置にあるので ``1 0 1`` と表現されます。

.. figure:: ../../img/obs_sec_fig1.*
   :name: bond_22
   :width: 150px

   ``Lsub=[2,2]`` の正方格子 S=1/2 ハイゼンベルグ模型のサイトの並び順

最後に演算子の要素です。
まずはサイトの基底を番号付ける必要がありますが、ここでは :math:`|\uparrow\rangle` を0, :math:`|\downarrow\rangle` を 1 とします。
この基底と番号を用いると、
例えば対角項の1つ :math:`\left\langle \uparrow_i \uparrow_j | \mathcal{H}_{ij} | \uparrow_i \uparrow_j \right\rangle = 1/4` は
``0 0 0 0 0.25 0.0`` と表現されます。
他に、非対角項の1つ :math:`\left\langle \uparrow_i \downarrow_j | \mathcal{H}_{ij} | \downarrow_i \uparrow_j \right\rangle = 1/2` は
``1 0 0 1 0.5 0.0`` と表現されます。

結果として、 S=1/2 のハイゼンベルグハミルトニアンは次のように定義されます。
::

  [[observable.twosite]]
  name = "hamiltonian"
  group = 0
  dim = [2, 2]
  bonds = """
  0 0 1
  0 1 0
  1 0 1
  1 1 0
  2 0 1
  2 1 0
  3 0 1
  3 1 0
  """
  elements = """
  0 0 0 0  0.25 0.0
  1 0 1 0  -0.25 0.0
  0 1 1 0  0.5 0.0
  1 0 0 1  0.5 0.0
  0 1 0 1  -0.25 0.0
  1 1 1 1  0.25 0.0
  """


``observable.multisite``
~~~~~~~~~~~~~~~~~~~~~~~~~

みっつ以上のサイト上で定義される物理量を示す演算子を定義します。
サイトごとの1体演算子の直積として定義されます。

.. csv-table::
   :header: "名前", "説明", "型"
   :widths: 15, 30, 20

   ``name``,       "演算子の名前",             文字列
   ``group``,      "演算子の識別番号",         整数
   ``multisites``, "サイトの組み合わせ",       文字列
   ``ops``,        "onesite 演算子の識別番号", 整数のリスト
   ``coeff``,    "演算子にかかる係数（実部）", 実数
   ``coeff_im``, "演算子にかかる係数（虚部）", 実数

``name`` は演算子の名前です。

``group`` は multisites 演算子の識別番号です。

``multisites`` は演算子が作用するサイト群の集合を表す文字列です。
整数からなる1行が1つのサイト群を意味します。

- 最初の整数は 始点サイト (source) の番号です。
- のこりの整数は source site からみた他サイト (target) の座標の組 (dx, dy) を並べたものです。

  - Nサイトの場合、 ``source_site dx2 dy2 dx3 dy3 ... dxN dyN`` という形式です。
  - :math:`4 \times 4` の正方形内に収まる必要があります。

``ops`` を用いて ``observable.onesite`` で定義した1体演算子の直積として演算子を定義します。
例えば ``observable.onesite`` の ``group=0`` として :math:`S^z` を定義していた場合には、
``ops = [0,0,0]`` として :math:`S^z_iS^z_jS^z_k` を表現できます。

``coeff``, ``coeff_im`` は演算子にかかる係数の実部と虚部です。
省略した場合はそれぞれ 1.0, 0.0 になります。
