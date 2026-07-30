# ARPACK-NG による転送行列固有値ソルバの設計

- 日付: 2026-07-30
- 目的: 相関長計算で使う転送行列の少数固有値計算を、自前 implicit restart Arnoldi
  (`src/arnoldi.cpp`)から ARPACK-NG に置き換え可能にし、信頼性と速度を向上させる。
- 背景: 外部ライブラリ調査(`docs/superpowers/notes/2026-07-20-external-library-survey.md`
  の項目 3)。自前実装は再直交化(DGKS 補正)がなく、収束判定が nev 番目の残差のみ、
  restart 時の Krylov 基底更新が要素単位の `get_value`/`set_value` ループで
  MPI 分散テンソルでは通信律速になる、という弱点を持つ。

## 要件(確定事項)

1. **依存はオプショナル+フォールバック**。ARPACK-NG が見つかればリンクし、
   なければ従来の自前 Arnoldi のみでビルドできる。ビルド障壁を上げない。
2. **実行時切替**。入力 toml でソルバを選択できる(既定は ARPACK があれば ARPACK)。
   同一バイナリで A/B 比較・トラブル時の切り分けを可能にする。
3. **検証は ctest ゴールデン一致+ベンチスイート追加**。相関長計算が支配的な
   ベンチケースを `benchmark/suites/` に追加し、ハーネスで再現可能な A/B 比較を残す。

## アプローチ(採用: 案 A)

**全ランク冗長シリアル ARPACK + 分散 matvec**。全ランクが同一の reverse
communication ループを冗長実行する。ベクトル長は N = χ²(CTM)または D²(MF)で
高々 10³–10⁴ 程度なので、ベクトルのフルコピーと ARPACK 内部の密演算の冗長実行は
コスト無視可能。重い matvec(環境テンソルとの縮約)は既存の分散実装を無変更で使う。

不採用案: rank 0 のみ実行+ワーカーループ(非対称制御フローと終了プロトコルが
必要でコード量に見合わない)、PARPACK 完全分散(N が小さく過剰。mptensor の
ブロックサイクリック分割との対応付けが難しく環境依存も増える)。

バインディングは ICB ヘッダ(`arpack.hpp`、ビルドフラグ依存)を使わず、
必要な 4 ルーチン `dnaupd_` / `dneupd_` / `znaupd_` / `zneupd_` を自前で
`extern "C"` 宣言して直接呼ぶ。どのビルドの libarpack でも動く。

## 全体構成

```
TransferMatrix::eigenvalues()          (src/iTPS/transfer_matrix.cpp)
  ├─ N ≤ maxdim_dense_eigensolver → 密ソルバ(現状のまま)
  └─ それ以外
       ├─ eigensolver = "auto"(既定): ARPACK ビルドなら ARPACK、なければ builtin
       ├─ eigensolver = "arpack":  ARPACK(未ビルドなら入力読込時に input_error)
       └─ eigensolver = "builtin": 現行 Arnoldi(src/arnoldi.cpp、無変更で維持)
```

## コンポーネント

### 1. `src/arpack_solver.{hpp,cpp}`(新規)

```cpp
template <class ptensor>
std::vector<std::complex<double>> arpack_eigenvalues(
    std::function<void(ptensor&, ptensor const&)> matvec,
    ptensor const& initial_vec, size_t nev,
    int ncv, int maxiter, double tol);
```

- 全ランクが同一の naupd ループを冗長実行。実/複素は `ptensor::value_type` で
  ディスパッチ(実: `dnaupd_`/`dneupd_`、複素: `znaupd_`/`zneupd_`)。
- 分散 ptensor ⇄ ローカル配列の変換ヘルパ:
  - 分散→ローカル: 所有要素だけ書き込んだゼロ埋め配列を `MPI_Allreduce(SUM)`。
    各要素の所有ランクは一意なので x+0 の加算のみとなりビットレベルで決定的。
  - ローカル→分散: 所有要素のみ `set_value`。
  - `_NO_MPI` ビルドでは単なるコピーに退化する。
- 固有ベクトルは不要(相関長は固有値の絶対値比のみ使う)ため `rvec = false`。
  実版は dr/di 配列から複素数を組み立てる。|λ| 降順にソートして返す
  (既存 `eigenvalues()` と同じ規約)。
- ARPACK の事前制約 `nev < ncv - 1`、`ncv ≤ N` は呼び出し前にクランプする。
  N が小さいケースは既存の `maxdim_dense_eigensolver` 分岐が密ソルバに逃がすため、
  クランプ不能な状況は原理的に到達しない(防御的にチェックはする)。
- 本体は `TENES_USE_ARPACK` 定義時のみコンパイル。
  `constexpr bool arpack_available()` を提供し、呼び出し側の分岐に使う。

### 2. ARPACK パラメータの対応

| ARPACK | 意味 | 対応 |
|---|---|---|
| `nev` | 求める固有値数 | `num_eigvals`(既存) |
| `ncv` | Krylov 部分空間次元 | `arnoldi_maxdim`(既存) |
| `tol` | Ritz 値の相対残差判定 | `arnoldi_rtol`(既存) |
| `iparam(3)` | 最大 restart 回数 | `arnoldi_maxiter`(既存) |
| `which` | 狙う固有値 | `LM` 固定(転送行列は常に絶対値最大) |
| `iparam(7)`/`bmat` | 問題のモード | mode 1(標準固有値問題)固定 |
| `iparam(1)` | シフト戦略 | exact shifts 固定 |
| `info`(入力) | 初期ベクトル | ユーザー指定固定(既存 `initial_vector` を渡す) |

- 新しい入力パラメータは追加しない(`eigensolver` を除く)。
- `arnoldi_restartdim` は ARPACK に対応物がなく builtin 使用時のみ有効
  (ドキュメントに明記)。
- 挙動差: ARPACK の tol は全固有対の Ritz 残差を見るため、builtin の
  「nev 番目のみ」より保守的。同じ rtol 値で反復が数回増え得るが許容する。

### 3. 入力仕様(`load_toml.cpp` ほか)

- `correlation_length` セクションに `eigensolver = "auto" | "arpack" | "builtin"`
  (既定 `"auto"`)を追加。
- enum(int)化して `TransferMatrix_Parameters` に保持し、`Bcast` の int 配列に載せる。
- `"arpack"` 指定かつ未ビルドは**入力読込時**に `input_error`(測定時まで遅延させない)。
- 未知の文字列も `input_error`。
- ドキュメント `docs/sphinx/{ja,en}/file_specification/correlation_length_section.rst`
  を更新(`eigensolver` の説明、`arnoldi_*` の適用範囲の注記)。

### 4. CMake

- 3 値オプション `ENABLE_ARPACK=AUTO|ON|OFF`(既定 AUTO)。
  - AUTO: 見つかれば使う。見つからなければ警告なしで builtin のみ。
  - ON: 必須。見つからなければ FATAL_ERROR。
  - OFF: 探さない。
- 検出: `find_package(arpack-ng CONFIG QUIET)` → 失敗時
  `find_library(arpack)` フォールバック。`ARPACK_ROOT` ヒントを受け付ける
  (`SCALAPACK_ROOT` と同じ流儀)。
- 見つかったら `TENES_USE_ARPACK` をコンパイル定義に追加し、`iTPS_impl` に
  リンク。共有ライブラリの libarpack なら Fortran ランタイムは自動で解決される
  (`enable_language(Fortran)` は不要)。
- CMake 設定サマリ出力に ARPACK の有無・パスを表示。

### 5. エラー処理

- `naupd` の `info < 0`: info 値と意味(ARPACK ドキュメントの対応表)を添えて
  `tenes::runtime_error`。
- `info == 1`(maxiter 到達): 現行 Arnoldi と同様、警告を出して得られた
  Ritz 値を返す(計算は継続)。
- `info == 3`(シフト適用不能): `arnoldi_maxdim`(ncv)を増やす提案付きエラー。
- `neupd` のエラーも同様に info 値付きで例外。

### 6. タイマーとベンチマーク

- `iTPS::measure_transfer_matrix_eigenvalues`(`src/iTPS/correlation_length.cpp`)
  に `ScopedTimer("measure/correlation_length")` を追加(既存の階層キー流儀)。
- 新スイート `benchmark/suites/correlation_length.toml`: χ を大きめにして
  N = χ² が `maxdim_dense_eigensolver` を超える(= Arnoldi/ARPACK パスに入る)
  設定の同一物理ケースを、`eigensolver = "builtin"` / `"arpack"` の 2 ケースとして
  並べる。builtin/arpack の比較は同一 run 内の 2 ケースの
  `measure/correlation_length` タイマー値の突き合わせで行う
  (`bench.py compare` はコミット間比較用であり、ケース間比較には使わない)。
  なお `"arpack"` ケースは ARPACK ありビルドが前提(なしビルドでは入力エラーで
  落ちる。スイートの README コメントに明記する)。

### 7. テスト

- **単体テスト**(新規、doctest): `TENES_USE_ARPACK` 時のみ登録。
  既知スペクトルの小行列(実・複素)で builtin / ARPACK / 厳密値の三者一致を確認。
- **回帰テスト**: 既存 ctest のゴールデン比較(`correlation_length.dat`、
  rtol 1e-3 / atol 1e-4)がそのまま検証になる。ARPACK 既定(auto)ビルドで
  全 ctest が通ること。`eigensolver = "builtin"` でも通ること。
- **MPI**: ENABLE_MPI ビルドでの ctest(既存の MPI 実行構成)で確認。

## スコープ外

- `src/arnoldi.cpp` の改善(再直交化追加、restart の要素アクセス最適化など)は
  行わない。フォールバックとして無変更で維持する。
- PARPACK による完全分散化。
- Arnoldi/ARPACK 以外の固有値ソルバ(LOBPCG 等)の追加。
- shift-invert 等のスペクトル変換(転送行列では不要)。

## 受け入れ基準

1. ARPACK なし環境: 従来どおりビルド・全 ctest パス(挙動不変)。
2. ARPACK あり環境: 既定(auto)で ARPACK が使われ、全 ctest パス。
3. `eigensolver = "arpack"` 指定 + ARPACK なしビルド: 入力読込時にエラー。
4. 単体テストで builtin / ARPACK / 厳密値が一致。
5. `benchmark/bench.py run` で新スイートが動き、builtin / arpack 両ケースの
   `timers.json` に `measure/correlation_length` が記録される。

## 関連資料

- 調査ノート: `docs/superpowers/notes/2026-07-20-external-library-survey.md`
- ベンチマークハーネス: `benchmark/README.md`
- ARPACK-NG: https://github.com/opencollab/arpack-ng
