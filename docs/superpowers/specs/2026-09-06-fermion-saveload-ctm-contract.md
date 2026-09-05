# フェルミオン CTM のチェックポイントを読み戻せるようにする — 振る舞い契約(2026-09-06)

原因究明: `work/fermion/saveload-ctm/FINDINGS.md`(最小再現つき)

## 直すこと

fermion モードで**フルアップデートを回した run の `tensor_save` は読み戻せない**。
`update_CTM()` が fermion 分岐で CTM 辺テンソルを畳んだ単層形(rank 3、`(CHI, CHI, D^2)`)に
するのに対し、`initialize_tensors()` は fermion でも rank 4 `(CHI, CHI, D, D)` で確保するため、
`load_tensor` が rank 不一致で投げる。SU のみの run では `optimize()` が `update_CTM()` を
呼ばないので露見しない。

## 前提として確認済みの事実

`Calc_CTM_Environment` / `Calc_CTM_Environment_density` の `initialize` は既定 `true` で、
`iTPS::update_CTM()` はどちらの分岐でもそのまま呼ぶ。**したがって読み込まれた CTM 環境は
boson でも fermion でも必ず捨てられ、ゼロから作り直される。** チェックポイント中の
`C1..C4` / `eTt,eTr,eTb,eTl` は現状どのモードでも再利用されない。
よって「保存する環境の中身を捨てる」変更で失われる情報は無い。

(**残課題**: そもそも環境テンソルを保存する必要があるのか — boson を含めた設計の見直し。
この契約の範囲外。)

## 要求

### R1: fermion モードの保存は環境を初期化形で書く

`iTPS::save_tensors()` は、`finfo.enabled` のとき `C1..C4` および `eTt/eTr/eTb/eTl` を、
**その run が `tensor_load` を使わずに開始した直後に持っていたのと同じ形・同じ内容**で
書くこと。すなわち `initialize_tensors()` が確保する形

- `C1..C4`: `Shape(CHI, CHI)`
- `eTt`: `Shape(CHI, CHI, vdim[1], vdim[1])`、`eTr`: `Shape(CHI, CHI, vdim[2], vdim[2])`、
  `eTb`: `Shape(CHI, CHI, vdim[3], vdim[3])`、`eTl`: `Shape(CHI, CHI, vdim[0], vdim[0])`

で、全要素ゼロ。ライブの(畳まれた)環境は書かない。**ゼロにするのは決定性のため**であり、
新規 run の確保直後の中身に依存しないようにする。

boson モード・有限温度モードの保存内容は**1 バイトも変えないこと**。

### R2: 読み戻せること

fermion + **CTM 環境**(`meanfield_env = false`)+ **フルアップデートあり**の run が
`tensor_save` したディレクトリを、同じ設定の run が `tensor_load` で読み込めること。
読み込んだ run の観測量が、保存した run の観測量と一致すること。

### R3: メッセージの空白落ち

`src/iTPS/saveload_tensors.cpp` の rank 不一致メッセージが `butloaded` と出力される
(`"legs, but"` と `"loaded one has"` の連結で空白が落ちている)。`but loaded` に直す。

## テストへの要求(テスト作成者が書く)

`test/fermion/free_fermion_saveload.py.in` は `meanfield_env = true` のため、
**fermion + CTM の save/load 往復をどのテストも通していない**。これがこのバグを通した穴である。

**(S-1)** fermion + CTM 環境 + フルアップデートありの save → load 往復の E2E を 1 件足すこと。
- run1: `tensor_save`、`meanfield_env = false`(既定)、`full_update.num_step >= 1`、
  `fastfullupdate = false`(fermion では必須)。
- run2: run1 の保存先を `tensor_load` して同じ設定で走らせる。
- 判定: run2 が正常終了し、**run1 と run2 の観測量が一致する**こと。許容誤差は、
  比較する 2 つの run の収束条件と同じ大きさにならないよう選ぶこと
  (`docs/superpowers/specs/2026-09-05-fermion-test-tolerance-contract.md` の規約 2 と同じ理由。
  入力に収束設定を明記し、選んだ許容誤差と実測の最悪相対誤差を報告に書く)。
- **反空洞化**: run1 がフルアップデートを実際に 1 回以上実行したこと(= `update_CTM()` を
  通ったこと)を前提として確認すること。SU のみに落ちると、このテストは修正前でも緑になり
  何も検査しなくなる。保存された `El_0.dat` の `ndim` を読むなどして機械的に確認すること。

**(S-2)** 保存された環境テンソルが R1 の形であること(`El_0.dat` の `ndim` が 4)を、
S-1 の中で直接検査すること。

**(S-3)** boson の保存内容が変わらないこと。既存の save/load テストが緑のままであることで足りる。

## 完了条件

1. 修正前に S-1 が **RED**(rank 不一致で落ちる)であることを確認していること。
2. `ctest --preset gcc-release` が全件緑。
3. `work/fermion/saveload-ctm/in_save_fu.toml` → `in_load_fu.toml` の手動再現が通ること。
4. 報告 `work/fermion/saveload-ctm/report-*.md`。
