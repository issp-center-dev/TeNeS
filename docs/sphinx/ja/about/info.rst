概要
=================
TeNeS (**Te** nsor **Ne** twork **S** olver) はテンソルネットワーク法に基づく多体量子状態計算のためのオープンソースのプログラムパッケージです。
二次元格子上で定義された量子スピン模型などの多体ハミルトニアンについて、磁化や相関関数などの物理量を計算します。
基底状態計算の他に、有限温度計算や時間発展計算も可能です。
あらかじめ定義された模型・格子に対しては、ユーザーが簡単に入力ファイルを作成するためのツールがあり、気軽に体験できます。
OpenMP/MPI ハイブリッド並列に対応しており、大規模計算機による大規模計算が可能です。

開発者
==================
TeNeS は以下のメンバーで開発しています。

- 大久保 毅 (東京大学大学院 理学系研究科)
- 森田 悟史 (慶応大学大学院 理工学研究科)
- 本山 裕一 (東京大学 物性研究所)
- 吉見 一慶 (東京大学 物性研究所)
- 青山 龍美 (東京大学 物性研究所)
- 加藤 岳生 (東京大学 物性研究所)
- 川島 直輝 (東京大学 物性研究所)

バージョン履歴
==================

- ver. 2.1.2: 2025-06-30 にリリース。
- ver. 2.1.0: 2024-02-28 にリリース。
- ver. 2.0.0: 2023-11-17 にリリース。
- ver. 2.0-beta: 2023-10-25 にリリース。
- ver. 1.3.4: 2023-09-13 にリリース。
- ver. 1.3.3: 2023-07-14 にリリース。
- ver. 1.3.2: 2023-06-08 にリリース。
- ver. 1.3.1: 2022-10-21 にリリース。
- ver. 1.3.0: 2022-10-20 にリリース。
- ver. 1.2.0: 2021-12-13 にリリース。
- ver. 1.1.1: 2020-11-09 にリリース。
- ver. 1.1.0: 2020-07-09 にリリース。
- ver. 1.0.0: 2020-04-17 にリリース。
- ver. 1.0-beta: 2020-03-30 にリリース。
- ver. 0.1: 2019-12-04 にリリース。

ライセンス
==================

本ソフトウェアのプログラムパッケージおよびソースコード一式はGNU General Public License version 3（GPL v3）に準じて配布されています。

論文
========

TeNeS を用いた結果を出版する場合には、次の論文を引用していただけると幸いです。

- `Y. Motoyama, T. Okubo, K. Yoshimi, S. Morita, T. Kato, and N. Kawashima, "TeNeS: Tensor Network Solver for Quantum Lattice Systems", Comput. Phys. Commun. 279, 108437 (2022) <https://www.sciencedirect.com/science/article/pii/S0010465522001564>`_
- `Y. Motoyama, T. Okubo, K. Yoshimi, S. Morita, T. Aoyama, T. Kato, and N. Kawashima, "TeNeS-v2: Enhancement for real-time and finite temperature simulations of quantum many-body systems", Comput. Phys. Commun. **315**, 109692 (2025) <https://www.sciencedirect.com/science/article/pii/S0010465525001948>`_

コピーライト
==================

© *2019- The University of Tokyo. All rights reserved.*

本ソフトウェアの一部は2019年度 東京大学物性研究所 ソフトウェア高度化プロジェクトの支援を受け開発されており、その著作権は東京大学が所持しています。
