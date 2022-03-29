.. highlight:: none

スタンダードモードを用いた格子・模型・演算子の定義
----------------------------

シンプルモードは定義済み模型・格子のパラメータからハミルトニアンやユニットセル情報を生成するツールであるが、スタンダードモードでは格子・模型・演算子を定義することができる。ここではスタンダードモードの使い方を説明する。

.. figure:: ../../img/ja_tutorial_1_Network.*
     :name: fig_transverse
     :width: 400px
     :align: center

     4サイトユニットセルのテンソルネットワーク

ユニットセル(格子)の定義
----------------------------

.. figure:: ../../img/ja_tutorial_1_Tensor.*
     :name: fig_transverse
     :width: 400px
     :align: center

     [tensor]と[[tensor.unitcell]]

ユニットセルの定義は[tensor]と[[tensor.unitcell]]を用いて行う:

.. code::

   [tensor]/*どういう格子を定義したか？ 
   type = "square lattice" # 無視される (simple.toml の名残)
   L_sub = [2, 2] # \ :math:`2\times 2`\ unitcell
   skew = 0 # \ :math:`y`\ 方向の境界を越えた時の\ :math:`x`\方向のずれ

   [[tensor.unitcell]]
   virtual_dim = [4, 4, 4, 4] # ボンド次元 (\ :math:`\leftarrow~\uparrow~rightarrow~downarrow`\の順)
   index = [0, 3] #　ユニットセル中のどのテンソルかを示す番号
   physical_dim = 2 # 物理ボンドの次元
   initial_state = [1.0, 0.0] # 初期状態の係数
   noise = 0.01 # 初期テンソルのゆらぎ


全体の初期状態\ :math:`\psi`\はサイトごとの初期状態\ :math:`\psi_i`\の直積状態\ :math:`| \Psi \rangle = \otimes_i |\Psi_i\rangle`\
\ :math:`\psi_i`\は、initial_state配列の要素を前から\ :math:`a_0,a_1,\cdots,a_{d-1}`\とすると、

.. math::

   \begin{aligned}
   |\Psi_i\rangle \propto \sum_k^{d-1}a_k|k\rangle\end{aligned}



模型(ボンドハミルトニアン)の定義
----------------------------

.. figure:: ../../img/ja_tutorial_1_Hamiltonian.*
     :name: fig_transverse
     :width: 400px
     :align: center

     [[hamiltonian]]


TeNeSが扱うハミルトニアンはボンドハミルトニアン (2サイトハミルトニアン) の和
(磁場などのサイトハミルトニアンも近くのボンドハミルトニアンを取り入れる)

.. math::

   \begin{aligned}
   mathcal{H} = \sum_{i,j}\mathcal{H}_{i,j}\end{aligned}

ボンドはsourceサイトとtargetサイトの組であると考える

ボンドハミルトニアンは、その行列要素と作用するボンドで規定する
行列要素を定義すれば模型を定義できる
ボンドを定義すれば格子を定義できる
source, targetが同じ番号のテンソルになるのは禁止する


std.tomlでのボンドハミルトニアンの定義

ボンドハミルトニアンの作用するボンドの定義
.. code::

   [[hamiltonian]]
   dim = [2, 2] # 作用するボンド [source, target] の取りうる状態数の対 
   bonds = """ # 作用するボンドの集合　(1行1ボンド)
   0 1 0 # 1列目: ユニットセル内のsourceの番号
   1 1 0 # 2列目: sourceからみたtargetの\ :math:`x`\座標(変位)
   2 1 0 # 3列目: sourceからみたtargetの\ :math:`y`\座標(変位)
   3 1 0
   0 0 1
   1 0 1
   2 0 1
   3 0 1
   """

0 1 0 は0番と右隣(1) (\ :math:`x+=1, y+=0`\)
1 0 1 は1番と上隣(3) (\ :math:`x+=0, y+=1`\)
1 1 0 は1番と右隣(1) 


ボンドハミルトニアン演算子の行列要素の定義
.. code::
   elements = """ # ハミルトニアンの(非ゼロな)行列要素(1行1要素)
   0 0 0 0 0.25 0.0 # 1列目: 作用前のsourceの状態
   1 0 1 0 -0.25 0.0 # 2列目: 作用前のtargetの状態
   0 1 1 0 0.5 0.0 # 3列目: 作用後のsourceの状態
   1 0 0 1 0.5 0.0 # 4列目: 作用後のtargetの状態
   0 1 0 1 -0.25 0.0 # 5列目: 要素の実部
   1 1 1 1 0.25 0.0 # 6列目: 要素の虚部
   """

0 0 0 0 0.25 0.0は\ :math:`\langle 00|\mathcal{H}_b|00\rangle=0.25`\
0 1 1 0 0.25 0.0は\ :math:`\langle 10|\mathcal{H}_b|01\rangle=0.5`\



演算子の定義
----------------------------

.. figure:: ../../img/ja_tutorial_1_Observable.*
     :name: fig_transverse
     :width: 400px
     :align: center

     [[observable.onesite]]


最終的に期待値を計算する演算子の定義
現在は1サイト演算子と2サイト演算子を計算可能

エネルギー演算子　= ボンドハミルトニアンも改めて指定する必要がある
(tenes_stdが0番の2サイト演算子として自動でコピーしてくれる)

1サイト演算子の数式は

.. math::

   \begin{aligned}
   S^z = \begin{pmatrix}
   0.5 & 0.0 \\ 0.0 & -0.5
   \end{pmatrix}\end{aligned}
 
である。

.. code::

   [observable]
   [[observable.onesite]] # 1サイト演算子
   name = "Sz" # 名前
   group = 0 # 1サイト演算子の識別番号
   sites = [] # 1サイト演算子が作用するテンソルの番号 ([]はすべてを意味する)
   dim = 2 # 1サイト演算子の次元
   elements = """ # 1サイト演算子行列の非ゼロ要素 (1行1要素)
   0 0 0.5 0.0 # 1,2列目: 作用前後の状態
   1 1 -0.5 0.0 # 3,4列目: 要素の実部・虚部
   """
   


最終的に期待値を計算する演算子の定義
現在は1サイト演算子と2サイト演算子を計算可能

エネルギー演算子　= ボンドハミルトニアンも改めて指定する必要がある
(tenes_stdが0番の2サイト演算子として自動でコピーしてくれる)

2サイト演算子の数式は

.. math::

   \begin{aligned}
   S^z_i S^z_j
   \end{aligned}
 
である。

.. code::

   [[observable].twosite]] # 2サイト演算子
   name = "SzSz" # 名前
   group = 1 # 2サイト演算子の識別番号 (1サイトとは独立)
   dim = [2, 2] # 次元
   bonds = """ # 作用するボンド (サイト対)
   0 1 0
   1 1 0
   2 1 0
   3 1 0
   0 0 1
   1 0 1
   2 0 1
   3 0 1
   """
   ops = [0, 0] # 1サイト演算子の直積で書ける場合、その識別番号
                # 今回は"Sz"が0番の1サイト演算子
                # elementsとして行列要素を陽に書くことも可能
                # (ボンドハミルトニアンと同様の書式)
   


   
   
反強磁性体の2次元ハイゼンベルグ模型のハミルトニアン
----------------------------

.. figure:: ../../img/ja_tutorial_1_2DHeisenberg.*
     :name: fig_transverse
     :width: 400px
     :align: center

     反強磁性体の2次元ハイゼンベルグ模型

.. code::

   [tensor]/*どういう格子を定義したか？ 
   type = "square lattice" # 無視される (simple.toml の名残)
   L_sub = [2, 2] # \ :math:`2\times 2`\ unitcell
   skew = 0 # \ :math:`y`\ 方向の境界を越えた時の\ :math:`x`\方向のずれ

   [[tensor.unitcell]]
   virtual_dim = [4, 4, 4, 4] # ボンド次元 (\ :math:`\leftarrow~\uparrow~rightarrow~downarrow`\の順)
   index = [0, 3] #　ユニットセル中のどのテンソルかを示す番号
   physical_dim = 2 # 物理ボンドの次元
   initial_state = [1.0, 0.0] # 初期状態の係数
   noise = 0.01 # 初期テンソルのゆらぎ
   
   [[tensor.unitcell]]
   virtual_dim = [4, 4, 4, 4] # ボンド次元 (\ :math:`\leftarrow~\uparrow~rightarrow~downarrow`\の順)
   index = [1, 2] #　ユニットセル中のどのテンソルかを示す番号
   physical_dim = 2 # 物理ボンドの次元
   initial_state = [0.0, 1.0] # 初期状態の係数
   noise = 0.01 # 初期テンソルのゆらぎ
   
   [[hamiltonian]]
   dim = [2, 2] # 作用するボンド [source, target] の取りうる状態数の対 
   bonds = """ # 作用するボンドの集合　(1行1ボンド)
   0 1 0 # 1列目: ユニットセル内のsourceの番号
   1 1 0 # 2列目: sourceからみたtargetの\ :math:`x`\座標(変位)
   2 1 0 # 3列目: sourceからみたtargetの\ :math:`y`\座標(変位)
   3 1 0
   0 0 1
   1 0 1
   2 0 1
   3 0 1
   """
   elements = """ # ハミルトニアンの(非ゼロな)行列要素(1行1要素), J=-1, h=1とする
   0 0 0 0 0.0 0.0 # 1列目: 作用前のsourceの状態
   1 0 1 0 0.25 0.0 # 2列目: 作用前のtargetの状態
   0 1 1 0 0.5 0.0 # 3列目: 作用後のsourceの状態
   1 0 0 1 0.5 0.0 # 4列目: 作用後のtargetの状態
   0 1 0 1 -0.75 0.0 # 5列目: 要素の実部
   1 1 1 1 0.0 0.0 # 6列目: 要素の虚部
   """
   
