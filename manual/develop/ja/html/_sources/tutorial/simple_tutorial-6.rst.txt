.. highlight:: none

三角格子上のハードコアボゾン模型の相図
-----------------------------------------------

最後に三角格子上で定義されたハードコアボゾン模型の絶対零度のおける相図の計算を紹介します。
模型のハミルトニアンは以下のようになります:

.. math::

   \begin{aligned}
   H = \sum_{\langle i,j \rangle} \Bigl[ -t (b_i^\dagger b_j + b_j^\dagger b_i) + V n_i n_j \Bigr] -\mu \sum_i n_i
   \end{aligned}

ここで\ :math:`\langle i, j\rangle`\ は隣接サイトの組を表し、\ :math:`\mu`\ は化学ポテンシャル、\ :math:`t`\ はホッピングエネルギー、\ :math:`V`\ は隣接サイト間の相互作用の大きさを表します。
ハードコアボゾン模型では、各サイトのボソン数は0か1に制限されます。
この模型では、2つのタイプの長距離秩序によって特徴づけられるいくつかの秩序相が現れることが知られています :ref:`[Wessel] <Ref-Wessel>` 。
1/3フィリングでは、 :numref:`fig_tutorial6_hardcore_boson` の挿入図に示すような\ :math:`\sqrt{3}\times\sqrt{3}`\ 超格子構造をもつ固体相が出現します。これは波数\ :math:`\bm{Q}=(4\pi/3,0)`\の構造因子\ :math:`S(\bm{Q})`\を計算することで特徴づけることができます。
もう一つの秩序は消滅演算子の期待値の大きさ\ :math:`|\langle b \rangle|`\で特徴づけられる超流動秩序です。

この模型の計算を行うには、 ``sample/06_hardcore_boson_triangular`` ディレクトリ中にある toml ファイル ``basic.toml`` , ``nn_obs.toml`` と、Pythonスクリプト ``run.py`` を利用します。 ``basic.toml`` ファイルには模型の設定やパラメータなどが記述されています。このファイルの記述は、前節の三角格子ハイゼンベルク模型とほぼ同じであるため、内容は割愛します。唯一の変更点は最後の ``model`` セクションだけです。ここで、 ``type = "boson"`` によってハードコアボゾン模型を指定しており、 ``t = 0.1`` , ``V = 1`` によってホッピングおよび隣接サイト間相互作用の大きさを指定しています。

``nn_obs.toml`` ファイルには計算する構造因子の情報が記述されています。

.. literalinclude:: ../../../../sample/06_hardcore_boson_triangular/nn_obs.toml

このファイルによって、構造因子\ :math:`S(\bm{Q})`\を計算することができます。

実際にスクリプト ``run.py`` を使って計算を実行してみましょう。
あらかじめ ``tenes`` などにパスを通した上で

::

    python run.py

により計算を実行します。数分〜十数分程度で計算が終わります。
計算が終了したら、gnuplotを起動し、

::

    load 'plot.gp'

とすれば、 :numref:`fig_tutorial6_hardcore_boson` (a)のような構造因子\ :math:`S(\bm{Q})`\と超流動秩序パラメータ\ :math:`|\langle b \rangle|`\のグラフが得られます。
なお、ここで計算に用いたボンド次元は小さい(ボンド次元2)ため、計算精度はそれほど高くありません。
スクリプト ``run.py`` の冒頭でコメントになっている4行をコメントアウトすることで、時間がかかりますが、より高精度の計算を行うこともできます。その結果を :numref:`fig_tutorial6_hardcore_boson` (b)に示します。
この図から、系の基底状態は3種類の相、つまり(a) 超流動相(\ :math:`-0.5 \lesssim \mu/V \lesssim -0.2`\), (b) 固体相(\ :math:`-0.2 \lesssim \mu/V \lesssim 2.4`\), および(c) 超流動固体相(\ :math:`2.4 \lesssim \mu/V`\)が現れることがわかります。これは先行研究の結果と整合します :ref:`[Wessel] <Ref-Wessel>` 。

.. figure:: ../../img/tutorial_6_hardcore_boson.*
	:name: fig_tutorial6_hardcore_boson
	:width: 600px

	三角格子上のハードコアボゾン模型の基底状態相図. (a) ボンド次元が2のとき. (b) ボンド次元が5のとき. 挿入図: 固体相のときの粒子配置. 

.. rubric:: 参考文献

.. _Ref-Wessel:

[Wessel] 
S. Wessel, M. Troyer, *Supersolid hard-core bosonson the triangular lattice*, Phys. Rev. Lett. **95**, 127205 (2005). `link <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.95.127205>`__.
