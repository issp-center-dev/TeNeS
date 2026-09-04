# フェルミオン full update 設計書

日付: 2026-09-04
ブランチ: `fermion`(HEAD `d6ab1902`、develop マージ済み)
前提: `docs/superpowers/specs/2026-08-15-fermion-m1-design.md`(graded 規約)、
`docs/superpowers/notes/2026-08-28-fermion-folding-methodology.md`(規約の導出・検証方法論)、
`src/fermion/fops.hpp`、`src/fermion/reduced.hpp`、`src/fermion/reduced_measure.hpp`。
設計時 rig: `work/fermion/full-update-design/`(`rig_block_ip.cpp`、`FINDINGS.md`)。

改訂ノート:

- 2026-09-04 初版。
- 2026-09-04 改訂1(Codex レビュー `work/fermion/full-update-design/codex-review-spec.md` 反映):
  §3.2 に k_A/k_B を独立次元とする helper 契約、§3.3 に forbidden block の合否判定、
  §4.3 にゲージ因子のパリティ射影と検査、§4.3 にセクタ次元 D_e/D_o の診断、
  §5.3 の fallback 地点を `main.cpp` の `Bcast` 前に修正、
  §7 に独立 oracle 脚・方向別検査・T3(i) の具体条件・T5 のコスト契約を追加。

## 1. 目的とスコープ

フェルミオン模式(`[parameter.general] fermion = true`)で full update
(`[parameter.full_update] num_step > 0`)を動かす。主目的は半充填 Hubbard の
**Mott 崩壊が simple update の局所最適化に起因するのか**を full update で検証すること
(`tenes-fermion-mott-collapse` メモ)。

スコープ内:

- 非高速版 full update(毎ボンド `update_CTM()` で環境を再収束)+ `Full_Gauge_Fix`。
- 水平ボンド(source_leg 2)と垂直ボンド(source_leg 3)。leg 0/1 は site 役割の入れ替えで正規化。
- ground-state 模式のみ(有限温度・実時間は既存ガードのまま)。

スコープ外:

- fast full update(`Full_Use_FastFullUpdate`)。素の Tn を使う `Left_move` 系は fold 環境と
  互換でない。fermion 模式では警告して非高速版に落とす(§5.3)。
- CTM のウォームスタート(`initialize = false`)。bosonic 非高速版と同じく毎回ゼロから収束させる。
  実行コストの支配項になるが、本設計では触らない。
- `lambda_tensor` の更新。bosonic full update も触らない。

## 2. 記法と既存規約

サイトテンソル Tn(l, t, r, b, s)。パリティ台帳 `finfo.phys[site]`、`finfo.virt[site][leg]`。
`mask_{xy}` は要素マスク (−1)^{p(x)p(y)}(x, y は脚)。`plain[...]` は mptensor の素の演算、
`fermion::` は fops の graded 演算。

依拠する既存規約(すべて検証済み):

- graded `tensordot` / `trace` / `transpose` の Koszul 符号、`conj` の (−1)^{m(m−1)/2}
  (m = 奇脚数)、`qr` / `svd` のパリティブロック分解と even-first 内部台帳、
  `svd_trunc` の大きさ順打ち切り後 even-first 並べ直し(`fops.hpp`)。
- `wrap_twosite_gate(op, p1, p2)`: 脚 (in1, in2, out1, out2)、入力脚スワップマスク事前適用。
- MF 内積: `contract_pair_MF(pair, op12) = trace(conj(pair), apply_pair_op(pair, op12), all, all)`。
- `build_pair_state(TnA, TnB, dir)` / `apply_pair_op(pair, op12)`(`reduced.hpp`)。
- `build_reduced_pair_halves(TnA, TnB, op12, dir)` と
  `contract_reduced_pair_halves_{horizontal,vertical}_density_CTM`(`reduced_measure.hpp`)。
- fold CTM 環境: `build_reduced_density_tensors(Tn, finfo)` →
  `Calc_CTM_Environment_density`(`measure.cpp:76-87`)。eT の bond 脚は (ket, bra) を融合した D²。

## 3. 環境 N の構成(開放チャネル fold)

### 3.1 QR と擬似サイト Q′

graded QR で環境側 Q と更新側 R に分ける(bosonic と同じ軸)。

| 方向 | site A(s1) | site B(s2) |
|---|---|---|
| 水平 | `qr(fTnA, Axes(0,1,3), Axes(2,4))` → QA(l,t,b,a), RA(a,r,s) | `qr(fTnB, Axes(1,2,3), Axes(0,4))` → QB(t,r,b,β), RB(β,l,s) |
| 垂直 | `qr(fTnA, Axes(0,1,2), Axes(3,4))` → QA(l,t,r,a), RA(a,b,s) | `qr(fTnB, Axes(0,2,3), Axes(1,4))` → QB(l,r,b,β), RB(β,t,s) |

内部脚 a, β の台帳 p_a, p_β は even-first。

Q を「内部脚を物理スロットに置き、ボンドスロットに次元 1 の偶ダミーを置いた」rank-5 の
擬似サイト Q′ に詰め替える。**素の reshape と台帳の手書きのみ、符号なし**(脚の順序を
変えないので Koszul 符号は生じない):

| 方向 | QA′ | QB′ |
|---|---|---|
| 水平 | (l, t, •, b, a) = `reshape(QA.t, Shape(l,t,1,b,nA))`、台帳 {p_l, p_t, {0}, p_b, p_a} | (•, t, r, b, β) |
| 垂直 | (l, t, r, •, a) | (l, •, r, b, β) |

### 3.2 開放チャネル fold

`build_reduced_pair_halves(QA′, QB′, op12, dir)`(`reduced.hpp:529-595`)は op12 を
graded SVD(`svd(op12, Axes(0,2), Axes(1,3))` → u(in_A,out_A,k), vt(k,in_B,out_B))して
両 half に振り分ける。ここでは SVD 因子の代わりに恒等因子を使い、演算子添字 k を
**開放したまま**環境に吸収する。

- 恒等因子: `uA = fermion::reshape(I4A, Shape(nA, nA, nA*nA))`
  (I4A(in,out,in′,out′) = δ_{in,in′}δ_{out,out′}、台帳 {p_a,p_a,p_a,p_a}、偶テンソル)
  → (in_A, out_A, k_A)、k_A の台帳は fuse(p_a, p_a)。
  `vB = fermion::reshape(I4B, Shape(nB*nB, nB, nB))`(I4B(in′,out′,in,out) = δδ)→ (k_B, in_B, out_B)。
  s は掛けない(`u.multiply_vector(s, 2)` を省く)。
- **k は奇セクタを含む**ので、TA6/TA5/TB6/TB5 の構成(graded tensordot、transpose、reshape、
  bond×k の crossing mask)は既存関数と**同一処理**でなければならない。重複実装を避けるため、
  `build_reduced_pair_halves` の SVD 後半を
  `build_reduced_pair_halves_from_factors(TnA, TnB, u, vt, dir)` に切り出し、既存関数と
  本設計の両方から呼ぶ。
- helper の契約(既存コードからの一般化): 既存経路は `nk = s.size()` 一つで A 側と B 側の
  両方を reshape しているが、開放経路では **k_A(次元 nA²)と k_B(次元 nB²)が独立**であり、
  full update の QR 内部次元は一般に nA ≠ nB(`envR1 ≠ envR2`)である。したがって helper は
  `nkA = u.parity[2].size()`、`nkB = vt.parity[0].size()` を**別々に**取り、
  TA5 の bundled 軸を nkA、TB5 の bundled 軸を nkB で作る。crossing mask は
  A 側 bond × k_A にのみ掛かる(既存実装も `TA5.multiply_vector(crossing_mask, ...)` のみで
  B 側には掛けていない)。閉じた経路だけが `nkA == nkB` を前提にできる。
- `PA = doubled_pipeline_traced(QA′, TA5)`(脚 2(水平)/脚 3(垂直)の次元が 1·nA²)、
  `PB = doubled_pipeline_traced(QB′, TB5)`(脚 0(水平)/脚 1(垂直)の次元が 1·nB²)。
- 環境吸収は `contract_reduced_pair_halves_{horizontal,vertical}_density_CTM`
  (`reduced_measure.hpp:317-385`)の **join 手前まで同一**。これも
  `absorb_reduced_pair_halves_{horizontal,vertical}(env×10, halves) → (left, right)` に切り出し、
  既存の閉じた版はその戻り値を join する。水平: left → (χ, k_A, χ)、right → (χ, k_B, χ)。
  垂直: top → (χ, k_A, χ)、bot → (χ, k_B, χ)。
- 開放 join: `M = plain::tensordot(left, right, Axes(0,2), Axes(0,2))` → (k_A, k_B)。
  閉じた版は `tensordot(left, right, Axes(1), Axes(1))` → `transpose(0,2,1,3)` →
  `trace_boundary_pairs`(idx[0]==idx[1] かつ idx[2]==idx[3] の和)であり、M はそれの k を
  開放した同じ縮約である。
- 線形性の根拠: 閉じた版は偶 O = Σ_k u_k ⊗ vt_k について Σ_k left(k)·right(k) を返す。
  pipeline は u と vt それぞれに線形で、crossing mask は k のパリティにしか依らないので、
  M(k_A, k_B) = 閉じた版の値 at O = u_{k_A} ⊗ vt_{k_B}(p(k_A) = p(k_B) のとき)。

### 3.3 N と N_plain

- `Ñ = ftensor{reshape(M, Shape(nA, nA, nB, nB)), {p_a, p_a, p_β, p_β}}` (in_A, out_A, in_B, out_B)。
- **forbidden block は合否判定に使う**。厳密な偶ネットワークでは p(k_A) ≠ p(k_B) のブロックは
  恒等的にゼロなので、非ゼロは開放 join・crossing mask・台帳のいずれかの異常信号である。
  `max_abs(M_forbidden) / max_abs(M) > 1e-8` なら `std::runtime_error`(値を報告)。
  **検査してからゼロ射影する**(丸め対策)。黙ってゼロ化すると符号バグを N 構成時に消してしまう。
- `N = fermion::transpose(Ñ, Axes(0,2,1,3))` → (in_A, in_B, out_A, out_B)、符号 (−1)^{p(out_A)p(in_B)}。
  根拠: graded SVD は `tensordot(u, vt, Axes(2), Axes(0))` を (in_A,out_A,in_B,out_B) の順で
  与え、wrap 形式 (in_A,in_B,out_A,out_B) に戻す graded transpose が同じ符号を持つ。
- 性質: 任意の偶 O(wrap 形式)で ⟨O⟩ = Σ_x N(x)·O(x)(**plain 要素積和**)。
  `fermion::trace(N, O, all, all)` は m = 2 の (−1)^{m(m−1)/2} = −1 が入るので使わない。
- **N_plain = N ⊙ mask_{in_A in_B}**(入力脚マスク)。これが ALS の環境になる(§4)。
- §3.2〜3.3 は rig では検証していない(rig は MF レベルで N を基底評価した)。
  CTM 経路での正しさは T2(i) が担う。

サイズ: D = 4, d = 4 で nA = nB = 16、k 次元 256、M は 65536 要素。

### 3.4 窓環境の選択(回転しない)

bosonic `Full_update_bond` は connect1 で Tn を回転して水平に揃えるが、fold 環境の回転は
符号規約が方向依存なので採らない。測定経路(`twosite_obs.cpp:215-300`)と同じ選択を使う:

| 方向 | C1, C2, C3, C4 | eT1..eT6 |
|---|---|---|
| 水平(s1 = 左, s2 = 右) | C1[s1], C2[s2], C3[s2], C4[s1] | eTt[s1], eTt[s2], eTr[s2], eTb[s2], eTb[s1], eTl[s1] |
| 垂直(s1 = 上, s2 = 下) | C1[s1], C2[s1], C3[s2], C4[s2] | eTt[s1], eTr[s1], eTr[s2], eTb[s2], eTl[s2], eTl[s1] |

## 4. ALS の graded 転写

### 4.1 なぜ literal 転写ができないか

graded `conj` は積に対して準同型でない:
conj(tensordot(A, B over m)) = tensordot(conj A, P_m·conj B)(奇ボンド成分に −1)。
bosonic ALS は conj(R1), conj(R2) を分離して使うので、fops 演算で literal に写すと誤る。
そこで費用関数を全体テンソルの内積から導出した。

### 4.2 定理(rig で検証済み、`work/fermion/full-update-design/FINDINGS.md`)

ブロック X(a, β, s1, s2) = `fermion::transpose(fermion::tensordot(R1, R2, Axes(1), Axes(1)), Axes(0,2,1,3))`、
R1(a, m, s1)、R2(β, m, s2)。pair_Q = `build_pair_state(QA′, QB′, dir)`。
N は MF レベル(CTM なし)で基底評価した bare-bra 汎関数
O ↦ `trace(conj(pair_Q), apply_pair_op(pair_Q, O))` の係数テンソル。
水平方向のみ、real/complex × 3 パリティ構成({D2,D3} × {d2,d4} の 3 組)× 2 シード、
264 チェック、unexpected 0。垂直方向は C2〜C7 が方向非依存(汎関数と R の代数だけ)なので
rig では走らせず、方向依存部(§3.1 の QR 軸・Q′ 詰め替え・N 構成)は T2/T3 で押さえる。

| # | 主張 | 誤差 |
|---|---|---|
| C1 | `build_pair_state(TnA, TnB)` == `apply_pair_op(pair_Q, X)` | 1e-15 |
| C2 | ⟨X\|Y⟩_E = Σ N·O_{X,Y}、O_{X,Y}(a,β,a′,β′) = mask_{a′β′} Σ_{s1,s2} Y(a,β,s1,s2)·conj X(a′,β′,s1,s2)。添字順は N(a,β,a′,β′) と一致 | 6e-14 |
| C3 | X = mask_{aβ} ⊙ plain[R̃1·R2]、R̃1 = R1 ⊙ mask_{m s1}、R2 は無変更 | 0 |
| C4 | N_plain はエルミート・半正定値、⟨X\|Y⟩_E = ⟨X̃\|N_plain\|Ỹ⟩(X̃ = mask_{aβ}X) | 4e-15 |
| C5 | Θ = `fermion::tensordot(X, G, Axes(2,3), Axes(0,1))` で apply(apply(P,X),G) = apply(P,Θ)。bosonic 順 `tensordot(tensordot(R1,R2,(1),(1)), G, Axes(1,3), Axes(0,1))` と同値 | 1e-15 |
| C6 | graded `svd_trunc(Θ, Axes(0,2), Axes(1,3))`(打ち切りなし)→ R1 = `transpose(U·s, Axes(0,2,1))`, R2 = `transpose(VT, Axes(1,0,2))` が Θ を再構成し、C3 のマスクで plain 再構成 | 4e-15 |
| C7 | N_plain を plain でエルミート化・eigh・正値化し bosonic ゲージ固定(`qr(Z, Axes(2,1), Axes(0))` 等)→ LR1/LR2 はブロック対角、Θ̃′ = LR1·LR2·Θ̃ を mask_{aβ} で戻した Θ′ と Env′ = Z·conj(Z) はパリティ偶を保つ | 5e-15 |
| M | C2 の mask_{a′β′}、C3 のマスクを落とすと O(1) で不一致(変異防衛) | — |

### 4.3 帰結: bosonic ALS コアを plain 配列でそのまま実行する

費用関数 ‖ψ(R1R2) − ψ(Θ)‖²_E は、マスク済み変数 (N_plain, Θ̃ = mask_{aβ}Θ, R̃1, R2) で
**bosonic と同一の二次形式**になる(C2〜C4)。したがって:

- エルミート化、eigh、正値化(`Inverse_Env_cut`)、ゲージ固定(`Full_Gauge_Fix`)、擬似逆、
  ALS 反復、収束判定、ゲージ除去は **bosonic コードを無改変で共有**する。`fermion::eigh` は不要。
- パリティ保存の根拠: LR1†LR1 = Σ_β Env⁺ は基底非依存でブロック対角なので、**フルランクなら**
  Cholesky/QR の一意性により LR1 もブロック対角(C7)。擬似逆と Env⁺ はスペクトル関数なので
  固有ベクトルの混合に影響されない。ALS の線形解は P 対称問題の一意解。
- ただしこの一意性は**フルランクでのみ**成り立つ。`Inverse_Env_cut` で正値化がランク落ちした場合、
  および偶/奇セクタで固有値が縮退して `eigh` が混合固有ベクトルを返した場合、上三角 R は
  一意でなくなり、LAPACK がブロック外成分を持つ因子を返し得る。したがって
  **`prepare_environment` の直後にパリティ射影と検査を置く**: LR1, LR2, LR1_inv, LR2_inv,
  Θ′, Env′ の forbidden block について `max_abs(forbidden) / max_abs(全体) > 1e-8` なら
  `std::runtime_error`、以下ならゼロ射影してから先に進む。射影誤差は違反量で抑えられる。
- graded が必要なのは次の 3 箇所だけ:
  1. **Θ の構成**: `Θ = fermion::tensordot(X, wrapped_gate, Axes(2,3), Axes(0,1))`、Θ̃ = mask_{aβ}Θ。
     `wrapped_gate` は driver が作る(§5.1): `wrap_twosite_gate(op, p1, p2)` と source swap 用の
     `fermion::transpose(·, Axes(1,0,3,2))` は **driver 側の責務**で、core は
     「wrap 済み・swap 済みのゲートだけを受ける」。simple update と同じ分担にして二重 wrap を防ぐ。
  2. **初期推定**: `prepare_environment` が返す Θ̃′ = LR1·LR2·Θ̃(plain)を mask_{aβ} で戻し
     ftensor 化する(Θ′、台帳 (p_a, p_β, p_s, p_s)。LR がブロック対角なので偶、C7)。
     graded `svd_trunc(Θ′, Axes(0,2), Axes(1,3), U, s, VT, D)` → bosonic と同じく
     λ = sqrt(s / ‖s‖) を U(軸 2)と VT(軸 0)に掛ける → R1 = `transpose(U, Axes(0,2,1))`、
     R2 = `transpose(VT, Axes(1,0,2))`(C6)→ R̃1_init = R1 ⊙ mask_{m s1}、R2_init = R2。
     plain SVD を使わない理由: 偶/奇セクタ間で特異値が偶然縮退すると LAPACK が混合ベクトルを返し、
     ボンド台帳が壊れる。**ボンド台帳のセクタ分割 D = D_e + D_o はここで決まり**、
     ALS とバランシングはそれを保存する(simple update と同じ大きさ順選択)。
     診断: 各ボンドの (D_e, D_o) を `print_level >= debug` でログし、いずれかが 0 になったら
     `print_level >= warn` で警告する。奇セクタが消えるとフェルミオン的な絡み合いが失われ、
     Mott 崩壊の判定に直結し得るため、Hubbard の検証では警告の有無を必ず確認する。
     「各セクタ最低 1 本を保持する」ガードは**入れない** — 変分原理が選んだ結果を人為的に
     曲げることになり、simple update とも規約が食い違うため。セクタ枯渇が起きたら
     それ自体を物理的な観測結果として扱う。
  3. **後処理**: ゲージ除去(`LR1_inv·R̃1`, `LR2_inv·R2`、plain)→ R1 = R̃1 ⊙ mask_{m s1} で戻す →
     ftensor 化(台帳 (p_a, p_m, p_s), (p_β, p_m, p_s))→ `parity_violation` を監視 →
     graded バランシング(`fermion::qr(R1, Axes(0,2), Axes(1))`、同 R2、`fermion::svd(r1·r2)`、
     s は未ソートだが Σs² 正規化なので順序無関係)→ graded で Tn_new を組み立て(§4.5)。

パリティ監視: `parity_violation(R) / max_abs(R) > 1e-8` なら `std::runtime_error`
(想定外の混合)、以下なら奇成分をゼロ化して先に進む(下流の `validate_block_diagonal` を守る)。

代替案(全 graded 化 + P_m 規則で conj の非準同型性を手で補正)は符号敏感な新規コードが
ALS 全体に広がるため採らない。

### 4.4 bosonic core の分割

`Full_update_bond_horizontal`(`src/iTPS/core/full_update.cpp:95-391`)を演算順序を変えずに
3 つに分ける:

- `prepare_environment(Environment, Theta, params) → (Environment′, Theta′, LR1_inv, LR2_inv)`:
  エルミート化(124 行)〜ゲージ固定(200 行)。`Full_Gauge_Fix = false` のときは
  LR_inv を恒等にし、`Environment = tensordot(Z, conj(Z), Axes(2), Axes(2))` の分岐を保つ。
- `als_iterate(Environment′, Theta′, R1_init, R2_init, params) → (R1′, R2′)`:
  C_phi・Old_delta の初期化から ALS ループと未収束警告まで(236〜358 行)。
  費用関数は plain `conj(R1)`, `conj(R2)`, `conj(Theta)` を使うが、これは C4 の
  ⟨X̃\|N_plain\|Ỹ⟩ そのものなので fermion 経路でもそのまま正しい。
- 呼び出し側固有: 初期推定(bosonic は 207〜235 行の plain `svd` + `slice`)、ゲージ除去(360〜364 行)、
  バランシング・組み立て(366〜390 行)。

bosonic 経路は既存 doctest `test/full_update.cpp` と golden(fast 経路も同じ core を通る)で
挙動不変を担保する。

### 4.5 Tn_new の組み立て(方向別)

バランシング後の R1 は (a, s, m′)、R2 は (β, s, m′)(bosonic と同じ順)。すべて graded 演算。

| 方向 | Tn1_new | Tn2_new |
|---|---|---|
| 水平 | `transpose(tensordot(QA, R1, Axes(3), Axes(0)), Axes(0,1,4,2,3))` | `transpose(tensordot(QB, R2, Axes(3), Axes(0)), Axes(4,0,1,2,3))` |
| 垂直 | `transpose(tensordot(QA, R1, Axes(3), Axes(0)), Axes(0,1,2,4,3))` | `transpose(tensordot(QB, R2, Axes(3), Axes(0)), Axes(0,4,1,2,3))` |

新しいボンド台帳は `fermion::svd` の内部台帳(even-first)で、Tn1_new/Tn2_new の該当脚に載る。

## 5. driver・ガード・台帳更新

### 5.1 `iTPS::full_update(up)` のフェルミオン分岐

`simple_update.cpp:63-130` を手本にする。

- onesite: `wrap_Tn` → `fermion::tensordot(fTn, fop, Axes(4), Axes(0))` → `unwrap_Tn`。
  simple update と同じコードなので共通ヘルパに括り出す。
- twosite: driver が `fop = wrap_twosite_gate(up.op, phys[s1], phys[s2])` を作る。
  raster 正規化 — `source_leg ∈ {0, 1}` なら site 役割を swap し、ゲートを
  `fermion::transpose(fop, Axes(1,0,3,2))`。以後 s1 = 左/上、s2 = 右/下、
  direction = horizontal(leg 2)/vertical(leg 3)。core に渡すのは **wrap 済み・swap 済み**の
  ゲートのみ(§4.3)。
- §3.4 の窓環境と `wrap_Tn(Tn[s1])`, `wrap_Tn(Tn[s2])`, fop, direction, `peps_parameters` を
  `core::Full_update_bond_fermion` に渡し、戻り値で
  `finfo.virt[s1][s1_leg] = fTn1_work.parity[s1_leg]`、`finfo.virt[s2][s2_leg] = fTn2_work.parity[s2_leg]`、
  `unwrap_Tn`。`unwrap_Tn` は台帳を書き戻さないのでこの明示更新は必須。
- 台帳更新の直後に `fermion::validate_neighbor_consistency(finfo, lattice)`
  (`fermion_info.hpp:75`)を呼ぶ。source/target leg の取り違えは次の `wrap_Tn` まで潜伏するので、
  ボンド単位で捕まえる。コストは O(N_UNIT·4) の比較なので常時有効でよい。
- そのあと `update_CTM()`(§5.2)。`lambda_tensor` は触らない。

### 5.2 `update_CTM()` のフェルミオン対応

`measure.cpp:76-87` のロジック(`build_reduced_density_tensors(Tn, finfo)` →
`Calc_CTM_Environment_density`)を `iTPS::update_CTM()`(`iTPS.cpp:446`)に移し、
`measure()` は CTM 時に無条件で `update_CTM()` を呼ぶ。full update の初回と毎ボンドの
環境再収束が測定と同じ一本の経路になる。

契約:

- `update_CTM()` は `finfo.enabled` が真のときだけ reduced density tensor を作って
  `Calc_CTM_Environment_density` を呼び、偽なら従来どおり `Calc_CTM_Environment` を呼ぶ。
- `initialize = true` は bosonic と同じ(ウォームスタートはスコープ外、§1)。
- 有限温度用の `update_CTM_density()`(`iTPS.cpp:454`)とは**別物のまま**にする。
  こちらは purification の密度行列 iTPS 用で、fermion fold とは無関係。
- 既存のフェルミオン測定テスト(FreeFermion 系、`mf_measure`、`fold_geometry` の CTM ケース)が
  この移動の挙動不変を担保する。

### 5.3 ガード

- `load_toml.cpp:640-642` の `throw_fermion_guard("full update")` を撤去する。
- `MeanField_Env && num_full_step > 0` は既存 `PEPS_Parameters::check` がそのまま効く。
- `Full_Use_FastFullUpdate`(既定値 true)は fermion 模式では**警告して非高速版に落とす**。
  **地点は `src/iTPS/main.cpp` の `gen_param()` 直後・`peps_parameters.Bcast(comm)`(249 行)直前**。
  `validate_fermion_constraints` は `const PEPS_Parameters &` を取り(`load_toml.hpp:103`)、
  しかも `Bcast` の後(`main.cpp:310`)に呼ばれるので、そこで書き換えても全 rank に伝わらない。
  `peps_parameters.fermion` と `num_full_step` はどちらも `gen_param` の時点で確定しており
  (`PEPS_Parameters.hpp:98,116,136`)、`mpirank` も 212 行で取得済みなので、この地点なら
  条件判定・警告・書き換えがすべて揃う:

  ```cpp
  if (peps_parameters.fermion && peps_parameters.Full_Use_FastFullUpdate &&
      std::any_of(peps_parameters.num_full_step.begin(),
                  peps_parameters.num_full_step.end(),
                  [](int n) { return n > 0; })) {
    if (mpirank == 0) {
      std::cerr << "WARNING: fermion mode disables Full_Use_FastFullUpdate "
                   "because the fast update reuses bare-Tn CTM moves that are "
                   "not fermion-aware in this version" << std::endl;
    }
    peps_parameters.Full_Use_FastFullUpdate = false;
  }
  ```

  書き換えを `Bcast` の前に置くので全 rank が同じ値を持つ。既定値のままで動くことを優先する。
  (`load_toml.cpp:557` の `has_positive_steps` は無名 namespace 内で main.cpp からは見えないので、
  条件は上のように直接書く。)
- `tenes_std` は fermion 入力でも `[[evolution.full]]` を既に出力する
  (`_validate_fermion_mode_input` は full update を見ない)。tool 側の変更は不要。

### 5.4 台帳と保存

`finfo.virt` は §5.1 で両端に書くので `tensor_save` の `fermion.dat` にそのまま反映される。
`svd_trunc` は残す特異値を even-first に並べ直すので、ボンド台帳は常に連続ブロックである。

## 6. コード配置

| 新規/変更 | 内容 |
|---|---|
| `src/fermion/reduced.hpp` | `build_reduced_pair_halves` の SVD 後半を `build_reduced_pair_halves_from_factors(TnA, TnB, u, vt, dir)` に切り出す(挙動不変) |
| `src/fermion/reduced_measure.hpp` | 吸収列を `absorb_reduced_pair_halves_{horizontal,vertical}(env×10, halves) → (left, right)` に切り出す(挙動不変) |
| `src/fermion/full_update_env.hpp`(新) | §3: Q′ 詰め替え、恒等因子、開放 join、`N`、`N_plain`。`build_full_update_environment(env×10, QA, QB, dir) → N_plain` と、テスト用に `N` も返す入口 |
| `src/iTPS/core/full_update.{cpp,hpp}` | §4.4 の分割。`prepare_environment` / `als_iterate` を hpp で公開 |
| `src/iTPS/core/full_update_fermion.cpp`(新) | `Full_update_bond_fermion(env×10, fTn1, fTn2, fop, direction, params, fTn1_new, fTn2_new)`。bosonic core と同じ明示的インスタンス化 |
| `src/iTPS/full_update.cpp` | §5.1 の分岐 |
| `src/iTPS/iTPS.cpp`, `src/iTPS/measure.cpp` | §5.2 |
| `src/iTPS/load_toml.cpp` | §5.3 のガード撤去 |
| `src/iTPS/main.cpp` | §5.3 の fast full update フォールバック(`Bcast` 前) |

## 7. 検証行列

検出力の役割分担を先に決める。閉じた fold 経路(`build_reduced_pair_halves` →
`contract_reduced_pair_halves_*_density_CTM`)は `test/fermion/fold_geometry.cpp` の
Fock oracle アンカー(T5 群)と有限パッチ厳密縮約で**独立に固定済み**なので、
T2(i) の「開放版 == 閉鎖版」は自己参照ではない。ただし T2(i) は Q′ 詰め替え(§3.1)を
両辺で共有するため、そこだけは T2(iv) の独立 oracle が受け持つ。

**重要(検出力の限界)**: T3(i) の厳密性検査は **N の符号を検出しない**。
費用関数 ‖ψ(R1R2) − ψ(Θ)‖²_N の目標 Θ がボンド次元 D で厳密に表現できる場合、
最小値 0 は N が正定値でありさえすれば計量によらず ψ = Θ で達成される。
したがって N の符号を誤っても(正定値である限り)T3(i) は緑のままになる。
役割分担は次のとおりで、契約書にもこの限界を明記する:

| 壊した箇所 | 検出するテスト |
|---|---|
| N の入力脚マスク、開放 join、crossing mask、N の transpose 符号 | **T2(iv) の独立 oracle**(および T4/T5 の E2E) |
| Q′ 詰め替え、QR 軸 | T2(iv)、T3(i) |
| Θ の graded 合成、`mask_{m s1}`、初期推定の転置、Tn_new 転置 | **T3(i)** |
| ゲージ因子のパリティ漏れ | T3(iii)、T3(vii) |

| # | 対象 | 検査 | 合格基準 | 形式 |
|---|---|---|---|---|
| T2 | N 構成(§3) | (i) **開放 == 閉鎖**: ランダム偶 O(wrap 形式)で Σ N·O == `contract_reduced_pair_halves_density_CTM(env, build_reduced_pair_halves(QA′,QB′,O,dir))`。水平・垂直 × real/complex × 非自明台帳(D=2 と D=3、d=2 と d=4)× `nA ≠ nB` になる形状。(ii) **forbidden block**: `max_abs(M_forbidden)/max_abs(M)` を報告し閾値内。(iii) N_plain のエルミート残差・最小固有値 ≥ −tol。(iv) **独立 oracle**: `fold_geometry.cpp` の治具(3×2 パッチの厳密縮約 / Fock oracle)を環境に使い、Σ N·O が厳密真値と一致。CTM ではなく有限窓の厳密縮約を真値にすることで、Q′ 詰め替えと開放 join を共有経路の外から押さえる。(v) **bosonic 退化**: 全偶台帳で N == `Create_Environment_two_sites`(fold 環境の eT を (ket,bra) 順に reshape)。これは bosonic 回帰であり、fermion 符号の検出力は無いものとして扱う | (i)(iv)(v) 1e-12 相対、(ii) 1e-8 相対 | doctest `test/fermion/full_update_env.cpp` |
| T3 | `Full_update_bond_fermion`(§4) | (i) **厳密性**(主変異防衛)。条件を契約で固定する: 初期状態はボンド Schmidt rank 1 で、**ゲート適用後に奇チャネル振幅が実際に立つ**基底の重ね合わせ(単なる占有数固有状態は不可 — ホッピングがゼロ作用になり検出力を失う)。ゲートは適用後の厳密 pair state の Schmidt rank が D 以下、かつ各パリティセクタ内でも保持次元に収まるもの。`Inverse_Env_cut` と `Full_Inverse_precision` は必要方向を落とさない十分小さい値に設定する。合格条件は出力の pair state == `apply_pair_op(pair_old, G)`(正規化後)。**水平・垂直の両方**で行う。(ii) 恒等ゲート → pair state 不変(両方向)。(iii) 出力の `parity_violation` ≤ 1e-12、および §4.3 のゲージ因子 forbidden block が閾値内。(iv) 全偶台帳で bosonic `Full_update_bond` と pair state 一致(ゲージ差を吸収するため Tn ではなく pair state で比較)。(v) 一般ゲートで bra≠ket fold oracle による費用 ‖ψ_new − ψ_Θ‖²_E が初期推定の費用以下。(vi) `Full_Gauge_Fix = false` でも (i) が通る。(vii) 人工的な N_plain(偶奇で固有値が縮退するもの、`Inverse_Env_cut` でランク落ちするもの)で §4.3 の射影・検査が働く。(viii) **変異テストを契約に含める**: `mask_{m s1}` / Θ の graded 合成 / Q′ 詰め替え / Tn_new 転置 のいずれか 1 つを外したコピーで (i) が赤になること(N_plain の入力マスクは上表のとおり T3(i) では捕まらないので、その変異は T2(iv) で確認する) | (i)(ii)(iv) 1e-8、(iii) 1e-12 | doctest `test/fermion/full_update_bond.cpp` |
| T4 | driver、全偶同値 | AFH 2×2(`test/data/AntiferroHeisenberg_real.toml` から `[correlation]` を除き `fastfullupdate = false`)を fermion=false と fermion=true + parity=[0,0] で走らせ、エネルギー・onesite・twosite を比較 | 1e-6 相対(fold CTM と素 CTM の収束差を許容) | `test/fermion/boson_equivalence_full.py.in` |
| T5 | 物理 | 自由フェルミオン **D=2、小さい CHI、`Max_CTM_Iteration` を明示的に小さく**した入力で simple update 後に full update 数ステップ → エネルギー非増加かつ厳密値へ接近。source_leg 0/1 の swap を含む垂直ボンドを必ず通す配置にする。ctest には**明示 TIMEOUT** を付ける | E_FU ≤ E_SU + 1e-8、\|E_FU − E_exact\| < \|E_SU − E_exact\| | `test/fermion/free_fermion_full.py.in` |
| T6 | ガード(§5.3) | fermion + `num_full_step > 0` が受理される。fermion + `fastfullupdate = true` は警告して非高速版で完走(出力に警告文字列)。`meanfield_env = true` + full は従来どおり拒否。source_leg 0/1 swap 後に `validate_neighbor_consistency` が通ること | — | `test/input.cpp` 追記 + python |
| T7 | bosonic 回帰(§4.4) | 既存 doctest `test/full_update.cpp`(`full_update.dat`, 1e-8)と golden(`AntiferroHeisenberg_*`, `Honeycomb`, `J1J2_AFH`)、および既存 fermion 測定テスト(FreeFermion 系、`mf_measure`、`fold_geometry`)が無変更で緑 | 既存許容 | 既存テスト |

### 7.1 コストの見積もりとテスト時間

非高速版はボンドごとに CTM を再収束させる。2×2 の AFH 入力は `[[evolution.full]]` が
**8 個**(`test/data/AntiferroHeisenberg_real.toml`)なので、**1 full step あたり
`update_CTM()` が 8 回**(加えて `full_update()` 冒頭で 1 回)。これがコストの支配項である。

N 構成自体は軽い。D=4・d=4 では QR 内部次元 nA = nB = min(D³, D·d) = 16、開放 k は 256、
M は 65,536 要素。Q′ half(PA)は (l², t², nk_A, b²) = (16, 16, 256, 16) ≈ 1.05×10⁶ 要素
(実数で約 8 MB)なので、メモリは問題にならない。

したがって:

- **T5 は D=2・小 CHI に固定**し、明示 TIMEOUT を付ける(既存 fermion ctest の TIMEOUT 指定に倣う)。
- **D=4・d=4 の Hubbard は ctest に入れない**。`work/fermion/full-update-mott/` の性能・物理検証に分離する(§9)。

## 8. タスク分割(1 タスク = 1 Codex バッチ)

| # | タスク | 触る場所 | テスト |
|---|---|---|---|
| 1 | bosonic core を `prepare_environment` / `als_iterate` / 後処理に分割(挙動不変) | `src/iTPS/core/full_update.{cpp,hpp}` | T7 |
| 2 | `update_CTM()` をフェルミオン対応にし `measure.cpp:76-90` を委譲 | `src/iTPS/iTPS.cpp`, `measure.cpp` | 既存 FreeFermion 系 |
| 3 | `reduced.hpp` / `reduced_measure.hpp` の切り出し(k_A/k_B 独立化を含む)+ 開放チャネル N 構成 + N_plain + forbidden block 検査 | `src/fermion/reduced.hpp`, `reduced_measure.hpp`, `full_update_env.hpp` | T2、既存 fermion 系(切り出しの挙動不変) |
| 4 | `Full_update_bond_fermion`(QR・Θ・初期推定・ALS 呼出・ゲージ因子のパリティ射影と検査・セクタ次元診断・バランシング・組み立て、水平/垂直) | `src/iTPS/core/full_update_fermion.cpp` | T3 |
| 5 | driver 分岐・隣接整合性検査・ガード撤去・fast フォールバック警告 | `src/iTPS/full_update.cpp`, `load_toml.cpp`, `main.cpp` | T4, T5, T6 |

依存: 1 → 4、3 → 4、2 → 5、4 → 5。1〜3 は並列可。

## 9. 進行手順

設計書(本書)をコミット → Codex レビュー(**完了**、指摘は改訂1 で反映。
`work/fermion/full-update-design/codex-review-spec.md`)→ 契約書(散文、テストコード無し)→
テスト作成者サブエージェント → Claude が RED を 1 件ずつ確認しバイト単位スナップショット →
Codex 実装(テスト変更禁止・formatter 不実行・報告ファイル必須・確認を求めず完了まで)→
Claude 独立検証(ctest 全件)→ タスクレビュー → 最上位モデルで全ブランチレビュー 1 回。
そのあと本来の目的である **Hubbard U=4 の Mott 崩壊検証**を `work/fermion/full-update-mott/` で
実施する(ctest 外、複数シード)。

## 10. 残課題(本設計の外)

- CTM ウォームスタート(`initialize = false`)による毎ボンド再収束コストの削減。
- fast full update のフェルミオン対応(fold 環境の `Left_move` 系)。
- full update 後の `lambda_tensor` の陳腐化(MF 測定との整合)。bosonic と同じ状態。
- ボンド台帳のセクタ枯渇(D_e = 0 または D_o = 0)。§4.3 の診断で観測するだけにとどめ、
  対処(セクタ保持ガード、初期推定の与え方)は物理的な必要性が実測で示されてから検討する。
