# フェルミオン iPEPS 対応調査 要約(2026-08-13)

4 系統の並列調査(定式化理論/既存 OSS 実装/応用・実用詳細/TeNeS コードベース影響範囲)
の要約。詳細レポート(artifact): https://claude.ai/code/artifact/653868dd-dff9-4fb7-9333-98a8fe81b69e

## 主要結論

1. **フェルミオン iPEPS は確立技術**(2009–2010 年確立)。符号処理を加えても縮約計算量の
   リーディングオーダーはボゾン系と同一(Barthel–Pineda–Eisert, PRA 80, 042333 (2009),
   arXiv:0907.3689)。
2. **3 つの等価な定式化**: swap gate(Corboz–Orús–Bauer–Vidal, PRB 81, 165104 (2010),
   arXiv:0912.0646)/ Z₂グレード(統一レビュー: Mortier et al., SciPost Phys. 18, 012
   (2025), arXiv:2404.14611)/ Grassmann(格子場理論の TRG 系で使用、iTPS には利点なし)。
   共通核は (i) パリティ偶テンソル制約、(ii) 脚の基準順序、(iii) Koszul 符号則
   ((−1)^{|i||j|})の 3 要素。Jordan–Wigner の 2D 直接適用は非局所項を生むため
   generic PEPS に不適。
3. **CTMRG の骨格は不変**。変更は reduced tensor 構築時の符号・環境テンソルの偶性維持・
   射影子 SVD のセクター整合・観測量の JW string に局在。
4. **TeNeS の実装戦略は 2 択**:
   - (a) 密テンソル+swap gate: mptensor/ScaLAPACK 資産温存。TeNeS リポジトリ内で完結。
   - (b) Z₂グレード・ブロックスパース型: 生成カーネル無変更で通るが、TeNeS が使う
     mptensor API(自由関数 ~16 種+メンバ ~14 種 ≒ 30 API)の新実装+分散再設計が
     独立プロジェクト級。「ScaLAPACK 分散 × graded」の既製 C++ ライブラリは世界的に不在。
   - **採用: (a)**。U(1) ブロックスパース化は性能最適化として将来課題(正しさには不要。
     Z₂ 最小実装で D≲8–10 の検証は可能、文献水準の Hubbard 計算 D=12–18 には U(1) が
     デファクト前提)。
5. **コードベースの事実**: `tensordot` 呼び出し 1,562 箇所、うち 87% は自動生成カーネル
   (`src/iTPS/core/contract_*/` 54 ファイル)。生成器(`misc/contraction/` の
   cont*.jl + tdt.py/netcon.py)の中間表現はボンド名のみで**平面埋め込み(交差)情報を
   持たない**。M1 の主戦場(simple update・CTM 環境更新・1/2 サイト測定)は手書きコードで
   生成器の対象外。汎用 NxM カーネルの fermion 化(M2)で生成器改修が必要になり、
   これは三角格子 iTPS 計画と同じ改修点(幾何一般の中間表現で一度に設計すべき)。
6. **C++/MPI でフェルミオン iPEPS+CTMRG を備えた公開コードは不在**(YASTN・PEPSKit.jl は
   Python/Julia・単ノード)。TeNeS が対応すれば希少ポジション。
7. **検証の標準経路**: 自由スピンレスフェルミオン(厳密解)→ t-V → 小 D の Hubbard/t-J。
   golden data は YASTN / PEPSKit.jl と突き合わせ可能。

## 実装の設計文書として読む文献(優先順)

1. Mortier et al., arXiv:2404.14611 — 統一レビュー(swap gate と graded の等価性)
2. Corboz–Orús–Bauer–Vidal, arXiv:0912.0646 — fPEPS+CTMRG 原典
3. Bruognolo et al., SciPost Phys. Lect. Notes 25 (2021), arXiv:2006.08289 — 実装手引き
4. YASTN: arXiv:2405.12196 — swap gate 実装規約(平面射影+正準順序+交差 swap)
5. symmray(github.com/jcmgray/symmray)+ arXiv:2410.02215 — 符号簿記の最小実装見本

関連 spec: `docs/superpowers/specs/2026-08-15-fermion-m1-design.md`
