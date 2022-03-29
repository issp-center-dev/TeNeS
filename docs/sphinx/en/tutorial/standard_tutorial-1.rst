.. highlight:: none

Definition of lattices, models, and operators using the standard mode
----------------------------

While the simple mode is a tool to generate Hamiltonian and unit cell information from parameters of predefined models and lattices, the standard mode enables the definition of lattices, models, and operators. In this section, we explain how to use the standard mode.

.. figure:: ../../img/en_tutorial_1_Network.*
     :name: fig_transverse
     :width: 400px
     :align: center

     Tensor network of 4-site unit cell

Definition of unit cell (lattice)
----------------------------

.. figure:: ../../img/en_tutorial_1_Tensor.*
     :name: fig_transverse
     :width: 400px
     :align: center

     [tensor] and [[tensor.unitcell]]

Unit cells are defined using [tensor] and [[tensor.unitcell]]:

.. code::

   [tensor]/*What kind of lattice did you define? 
   type = "square lattice" # Ignored (remnant of simple.toml)
   L_sub = [2, 2] # \ :math:`2\times 2`\ unitcell
   skew = 0 # Displacement in \ :math:`x`\-direction 

   [[tensor.unitcell]]
   virtual_dim = [4, 4, 4, 4] # Bond dimensions (\ :math:`\leftarrow~\uparrow~rightarrow~downarrow`\ order)
   index = [0, 3] #　Number indicating which tensor is in the unit cell
   physical_dim = 2 # Physical bond dimensions
   initial_state = [1.0, 0.0] # Initial state coefficients
   noise = 0.01 # Initial tensor fluctuation


全体の初期状態\ :math:`\psi`\はサイトごとの初期状態\ :math:`\psi_i`\の直積状態で書ける。\ :math:`| \Psi \rangle = \otimes_i |\Psi_i\rangle`\
\ :math:`\psi_i`\は、initial_state配列の要素を前から\ :math:`a_0,a_1,\cdots,a_{d-1}`\とすると、

.. math::

   \begin{aligned}
   |\Psi_i\rangle \propto \sum_k^{d-1}a_k|k\rangle\end{aligned}



Definition of model (bond Hamiltonian)
----------------------------

.. figure:: ../../img/en_tutorial_1_Hamiltonian.*
     :name: fig_transverse
     :width: 400px
     :align: center

     [[hamiltonian]]


Hamiltonians handled by TeNeS are the sum of bonded Hamiltonians (2-site Hamiltonians)
(Site Hamiltonians such as magnetic fields are also incorporated as bond Hamiltonians)

.. math::

   \begin{aligned}
   mathcal{H} = \sum_{i,j}\mathcal{H}_{i,j}\end{aligned}

Consider a bond to be a pair of source and target sites

The Bond Hamiltonian is defined by its matrix elements and the bond it acts on.
Defining a matrix element enables us to define a model.
Defining a bond enables us to define a lattice.
It is prohibited for source and target to be tensors of the same number.


Definition of Bond Hamiltonian in std.toml

Definition of Bond Hamiltonian acting bond
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

.. figure:: ../../img/en_tutorial_1_Observable.*
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

.. figure:: ../../img/en_tutorial_1_2DHeisenberg.*
     :name: fig_transverse
     :width: 400px
     :align: center

     反強磁性体の2次元ハイゼンベルグ模型

.. code::

   [tensor]/*What kind of lattice did you define? 
   type = "square lattice" # Ignored (remnant of simple.toml)
   L_sub = [2, 2] # \ :math:`2\times 2`\ unitcell
   skew = 0 # Displacement in \ :math:`x`\-direction 

   [[tensor.unitcell]]
   virtual_dim = [4, 4, 4, 4] # Bond dimensions (\ :math:`\leftarrow~\uparrow~rightarrow~downarrow`\　order)
   index = [0, 3] #　Number indicating which tensor is in the unit cell
   physical_dim = 2 # Physical bond dimensions
   initial_state = [1.0, 0.0] # Initial state coefficients
   noise = 0.01 # Initial tensor fluctuation
   
   [[tensor.unitcell]]
   virtual_dim = [4, 4, 4, 4] 
   index = [1, 2] 
   physical_dim = 2 
   initial_state = [0.0, 1.0] 
   noise = 0.01 
   
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
