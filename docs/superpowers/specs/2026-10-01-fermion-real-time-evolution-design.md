# fermion モードの実時間発展 — 設計書

2026-10-01、ブランチ `fermion-real-time`(develop dbf05f18 から)。

## 1. 目的と範囲

fermion モード(`[parameter.general] fermion = true`)で `mode = "time"`(実時間発展)を使えるようにする。
主な用途は、偶パリティのプロダクト状態、または保存した基底状態から始めるクエンチである。

**やること**

- solver の入力ガードを、実時間発展は通し、有限温度だけ拒否するように狭める。
- simple update(SU)と full update(FU)の両方で実時間発展を通す。
- 正しさを、ボゾン等価性・自由フェルミオンの厳密解・load からのクエンチの 3 種類の参照で確かめる(§5)。
- ドキュメント(ja/en)の fermion 制限の記述と NEWS.md を更新する。

**やらないこと**

- 有限温度(`mode = "finite"`)。密度行列(iTPO)のアンシラにパリティの扱いが要り、別件とする。入力時の拒否を続ける。
- 奇パリティのプロダクト初期状態(spinless の CDW `|1010…>` など)。各サイトテンソルを偶に保つ規約に反するため、従来どおり拒否する。
  解除するには、奇サイト同士を奇の仮想ボンドで結ぶダイマー被覆が要る(M1 設計書 §6.1 で M2 送りとしたもの)。別件とする。
- 長距離ハミルトニアンボンドの FU。基底状態と同じく、`tenes_std` が拒否する現状を維持する。
- 実時間専用の最適化(FU の毎ステップの CTM warm start など)。必要だという根拠が出るまで作らない。

## 2. 現状

- `validate_fermion_constraints`(`src/iTPS/load_toml.cpp:656`)が `calcmode != ground_state` を一律に拒否している(`"non-ground-state mode"`)。
  このガードを直接テストしているものはない。
- `time_evolution()`(`src/iTPS/time_evolution.cpp`)には fermion 用の分岐がない。中で呼ぶものは、すべて fermion 経路をすでに持つ。
  - `simple_update(up)`:fermion の graded SU。
  - `full_update_in_sweep(up, …)` → `full_update(up)`:fermion の FU(fast / plain)。
    `time_evolution()` は冒頭で cold な `update_CTM()` を呼ばないが、冒頭の `measure(t, "TE_")` が `update_CTM()` を呼んで `ctm_valid_ = true` にする。
    fast FU 設計書(2026-09-07 §4.2)は、この経路を見越して `ctm_valid_` を導入済みである。
  - `measure(t, "TE_")`:fermion の CTM / 平均場測定。`time` と prefix の扱いはボゾンと共通である。
  - `fix_local_gauge()` は `Simple_Gauge_Fix = true` のときだけ呼ばれる。fermion はこれを入力時に拒否しているので、通らない。
- `run_timeevolution`(`src/iTPS/main.cpp:378`)は常に complex テンソルで走る。fermion の complex 経路は、
  `free_fermion_complex`(基底状態、複素ホッピング)と長距離ゲートの単体テスト(`is_complex`)が通している。
  ただし **complex の fermion FU を通したテストはない**(§6 のリスク R1)。
- `tenes_std` は fermion でも `mode` が `time` で始まれば `tau *= 1.0j` で実時間ゲートを作る。fermion を理由に拒否する箇所はない。
- `tenes_simple` の Hubbard は `initial = "cdw"`(doublon と holon の市松、各サイト偶パリティ)を持つ。§5 の厳密解テストにそのまま使える。
- ボゾンの実時間サンプル(`sample/02_time_evolution`)は、実数で保存した基底状態を complex の実時間計算で `tensor_load` している。
  fermion の load は `fermion.dat` のパリティ台帳を検証するが、実数で保存して complex で読む組み合わせはテストされていない(リスク R2)。

## 3. 変更

### 3.1 solver の入力ガード

`src/iTPS/load_toml.cpp` の

```cpp
if (peps_parameters.calcmode != PEPS_Parameters::ground_state) {
  throw_fermion_guard("non-ground-state mode");
}
```

を、`calcmode == finite_temperature` のときだけ `throw_fermion_guard("finite-temperature mode")` を投げるように変える。
他のガード(`Simple_Gauge_Fix`、`Use_RSVD`、マルチサイト、奇パリティのプロダクト初期状態、自己隣接セルなど)は変えない。

### 3.2 ツール

コードの変更は予定しない。`tenes_simple` → `tenes_std` が fermion の `mode = "time"` で `input.toml` を作れることを、テストで固定する(§5 (b) がこの経路を使う)。
P0(§4)や後続のテストで、ツールに実時間を妨げる箇所が見つかった場合だけ直す。

### 3.3 ドキュメント

- `docs/sphinx/{ja,en}/file_specification/parameter_section.rst` の fermion の記述:
  「基底状態計算に対応」を「基底状態計算と実時間発展に対応」に改め、非対応の列挙から実時間発展を外す。
  実時間発展の初期状態は、偶パリティのプロダクト状態、ランダム偶初期化、`tensor_load` の 3 つだと明記する。
- NEWS.md:fermion の項目に実時間発展の対応を足し、Limitations から real-time evolution を外す。

## 4. P0:実現性の事前計測(Claude が行う)

契約書を書く前に、ガードだけを外した捨てビルドで次を測り、テストの D・ステップ数・許容誤差を決める。
作業場所は `work/fermion-real-time/p0/`。成果物は計測記録だけで、コードは残さない。

1. §5 (b) の Hubbard CDW クエンチを SU で走らせ、D = 2, 3, 4 で厳密解との差を時刻ごとに記録する。
   差が D とともに縮むこと(打ち切り誤差であること)を確かめる。縮まないなら符号・位相の不具合を疑い、設計を見直す。
2. 同じ問題を FU で短く走らせ、完走すること(R1)と、SU と同程度に厳密解に沿うことを確かめる。
3. 実数の基底状態を `tensor_save` し、complex の実時間計算で `tensor_load` できることを確かめる(R2)。
4. 各テスト候補の実行時間を `OMP_NUM_THREADS=1` で測り、1 本 30 秒以内に収まるパラメータを選ぶ。

tenes の同時実行は 2 本までとする。

## 5. 検証

E2E は 1 本 30 秒以内、`OMP_NUM_THREADS=1` で測る。参照値の比較は、golden ファイルではなく同じビルドで回した値どうしか、解析解との比較で行う。

### (a) ボソン等価性(実時間の配線)

`test/data/TE_TFI.toml` と同じ模型・ゲート・初期状態・測定で、ボゾン版と fermion 版(`fermion = true`、全サイト `parity = [0, 0]`)を同じビルドで走らせ、`TE_*.dat` が一致することを確かめる。
fermion が拒否するマルチサイト観測量と相関長は、既存の `boson_equivalence_full` と同じく両方の入力から外す。
SU 版と FU 版の 2 本を作る。グレーディングが自明なので、ここで確かめるのは測定 prefix・時刻列・complex テンソル・TE ループと fermion 経路のつなぎ込みである。
許容誤差の決め方は `boson_equivalence_full` に倣い、実測の差に対して十分な余裕を持たせ、その根拠をテストのコメントに残す。

### (b) 自由フェルミオンの厳密解(非自明なグレーディング + complex ゲート)

2D 正方格子の Hubbard(`t = 1`、`U = 0`、`mu = 0`)を、`initial = "cdw"`(A 副格子が doublon、B 副格子が空)から実時間発展させる。
各スピン成分が独立な自由フェルミオンなので、サイトあたりの全密度の差は厳密に

  n_A(T) − n_B(T) = 2 · J0(4 t T)²

である(J0 は第一種ベッセル関数)。40×40 周期格子の一粒子計算で、T ≤ 0.8 の範囲で 1e-15 まで一致することを確認済み。
ボンド・時刻・D の選び方と許容誤差は P0 の計測で決め、契約書に書く。SU を主とし、FU も短いステップ数で 1 本作る。
入力は `tenes_simple` → `tenes_std` → `tenes` の正規経路で作る(§3.2 の固定を兼ねる)。

この比較は、2D の閉ループを回るホッピングでフェルミオン符号が効く点と、実時間ゲートが complex である点を同時に試す。
符号や位相を落とすと、短時間でも厳密解からずれる。

### (c) load からのクエンチ

fermion の基底状態計算(実数)で `tensor_save` し、別のハミルトニアンの実時間計算(complex)で `tensor_load` して始める。
確かめることは次の 2 点。
- load して時刻 0 で測った値が、保存元の基底状態計算の測定値と一致する。
- そのまま実時間発展が完走する。

### (d) ガードの単体テスト

- fermion で `mode = "time"` の入力を `validate_fermion_constraints` が受け付ける。
- fermion で `mode = "finite"` の入力は、メッセージに `finite-temperature` を含む `input_error` で拒否される。
- 他のガード(例:`Simple_Gauge_Fix = true`)は、`mode = "time"` でも引き続き拒否される。

## 6. リスク

- **R1:complex の fermion FU が未検証。** FU の ALS・パリティ違反検査・位相固定は、実数で調整されてきた。
  P0-2 で完走しない、または厳密解から大きくずれる場合は、原因を調べてから範囲を判断する(SU のみに縮める選択肢も含めて、ユーザーに相談する)。
- **R2:実数で保存して complex で読む fermion load。** `fermion.dat` の検証やテンソルのパリティ検査が複素数で通るかを P0-3 で確かめる。
- **R3:測定のコスト。** 実時間発展は毎回の測定で CTM を作り直す。d = 4 の Hubbard では、(b) の SU は平均場環境で測り、FU 版と (a) は小さい D と少ないステップ数で 30 秒に収める。
- **R4:FU の環境。** 実時間の FU は、毎ステップの CTM を前ステップの環境から warm start する(fermion の plain 経路)か、fast FU で移動させる。
  虚時間ほど状態がゆっくり動かないので、warm start の前提が弱まる可能性がある。P0-2 で収束状況を見る。

## 7. 進め方

CLAUDE.md の多段手順どおりに進める。

1. この設計書を Codex にレビューさせ、指摘を反映する。
2. P0(§4)を行い、結果で §5 のパラメータを確定する。
3. 振る舞い契約書(`docs/superpowers/specs/2026-10-01-fermion-real-time-evolution-contract.md`)を書き、テスト作成者が (a)〜(d) を書く。RED を 1 件ずつ確認し、`test/` 全体のスナップショットを取る。
4. Codex が §3 を実装する。テスト不改変を diff で検査し、Claude が ctest 全件を回す。
5. タスクごとのレビューのあと、最後に最上位モデルで全ブランチレビューを 1 回行う。
