# フェルミオン模式の tensor_save / tensor_load 対応 実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** フェルミオン模式で `[parameter.general] tensor_save` / `tensor_load` を使えるようにし、
モデルパラメータを変えながらのスキャンでチェックポイント再開ができるようにする。

**Architecture:** 仮想脚の偶奇台帳 `finfo.virt` を `<save_dir>/fermion.dat` に永続化する。
読み込みは 2 段: テンソルを読む**前**に台帳と格子同一性を検証して `finfo.virt` を復元し、
読んだ**後**にテンソル自身のパリティ違反を実測する。壊れた組合せはすべてエラーで止める。
仮想ボンド次元 D の変更は非対応(エラー)。

**Tech Stack:** C++17(mptensor、doctest)、Python 3(pytest / ctest)、CMake

**Spec:** `docs/superpowers/specs/2026-08-21-fermion-tensor-saveload-design.md`

## Global Constraints

- 役割分担は `~/.claude/CLAUDE.md` の多段エージェント方式。テストはテスト作成者(契約書のみ渡す)、
  実装は Codex(テストファイル変更禁止)、検証・整形・コミットは Claude
- **ボソン経路の数値挙動を変えない。** 既存 ctest 30 件(`restart` を含む)は緑を維持する
- フェルミオンのガードで外すのは `tensor load/save directories` **だけ**。モード・full update・
  RSVD・Gauge_Fix・相関・多サイト・skew・1 幅セル・距離 2 以上・`ops` 形式はすべて残す
- 利用者向け文言に `M1` / `M2` を書かない
- MPI: `fermion.dat` は rank 0 が読んで**生文字列のまま bcast** し、全ランクが同一に構文解析する
  (どれか 1 ランクだけが例外を投げると MPI がデッドロックするため)。パリティ違反の判定は
  `allreduce_max` で全ランクにまたがって取る
- ビルド: `cmake --build out-gcc/build -j 8`、C++ テスト: `out-gcc/build/test/test_fermion_layer`、
  全件: `ctest --test-dir out-gcc/build`
- 整形は Claude がコミット直前に、そのタスクで変更したファイルだけに `clang-format` / `black` を適用
- 各タスクの実装者は報告ファイル `work/sl-task<N>-report.md` を必ず作る(存在を Claude が確認する)

---

## ファイル構成

| ファイル | 役割 | 変更 | タスク |
|---|---|---|---|
| `src/iTPS/iTPS.hpp` | 3 メソッドの宣言 | 変更 | T1 |
| `src/iTPS/saveload_tensors.cpp` | `save_fermion_parity` / `load_fermion_ledger` / `validate_loaded_fermion_tensors` | 変更 | T1 |
| `src/iTPS/load_toml.cpp` | `tensor load/save` ガードの削除 | 変更 | T1 |
| `test/fermion/saveload.cpp` | 層1〜3(台帳の往復、測定値の往復、ガード) | 新規 | T1 |
| `test/test_fermion_layer.cpp` | 末尾に `#include "fermion/saveload.cpp"` を 1 行 | 変更 | T1 |
| `test/fermion/free_fermion_saveload.py.in` | 層4 E2E | 新規 | T2 |
| `test/CMakeLists.txt` | `FreeFermionSaveLoad` 登録 | 変更 | T2 |
| `docs/sphinx/{ja,en}/file_specification/parameter_section.rst` | 非対応一覧と注記 | 変更 | T2 |
| `NEWS.md` | 項目追加 | 変更 | T2 |

`test/fermion/saveload.cpp` を `test_fermion_layer.cpp` から `#include` するのは
`r2_convention.cpp` / `mf_measure.cpp` と同じ流儀(末尾に追加、`mf_measure.cpp` の後)。
同じ翻訳単位なので `iTPSTestAccessor` と既存ヘルパが使える。

---

## Task 1: 台帳の永続化と検証

**Files:**
- Modify: `src/iTPS/iTPS.hpp`(`load_tensors_v0()` の宣言の近く)
- Modify: `src/iTPS/saveload_tensors.cpp`
- Modify: `src/iTPS/load_toml.cpp:651-653`
- Create: `test/fermion/saveload.cpp`
- Modify: `test/test_fermion_layer.cpp`(末尾に `#include` を 1 行)

**Interfaces:**
- Consumes: `finfo`(`tenes::fermion::FermionInfo`)、`lattice`(`LX` / `LY` / `skew` / `N_UNIT`)、
  `Tn`、`util::path_exists` / `util::drop_comment` / `util::split`、`bcast(std::string&, 0, comm)`、
  `tenes::allreduce_max(std::vector<double>&, comm)`、`tenes::fermion::wrap_Tn` /
  `parity_violation` / `max_abs`、`tenes::load_error`
- Produces(T2 が使う): `<save_dir>/fermion.dat`(§3.1 の形式)、
  および fermion 模式で `tensor_save` / `tensor_load` が受理されること

### 契約書(テスト作成者に渡す散文)

**対象.** `tenes::itps::iTPS<tensor>` の `save_tensors()` / `load_tensors()` を
`peps_parameters.fermion = true` で使ったときの振る舞いと、
`validate_fermion_constraints`(`src/iTPS/load_toml.cpp`)が `tensor_save` / `tensor_load` を
受理すること。

**背景(読むもと).**
- `src/iTPS/saveload_tensors.cpp`: 既存の `save_tensors` / `load_tensors` / `load_tensors_v1`。
  保存されるのは `T_*`(5 脚のサイトテンソル)、`El/Et/Er/Eb_*`、`C1..C4_*`、`lambda_*`、`params.dat`。
- `src/fermion/fermion_info.hpp`: `FermionInfo{enabled, phys, virt}`、`even_first_parity(dim)`、
  `wrap_Tn`、`validate_neighbor_consistency`。
- `src/iTPS/tensors.cpp:64-90`: `finfo` の初期化。`virt` の初期値は `even_first_parity(vdim[leg])`。
- `test/test_fermion_layer.cpp`: `iTPSTestAccessor`(`Tn` / `lambda_tensor` / `finfo` にアクセスできる)、
  末尾付近の `MF layer4 ...` 群(2×2 セルで `iTPS` を組み立てて測る手本)。
- `test/fermion/mf_measure.cpp`: 決定論テンソルの作り方(`make_mf2_tensor` 等)。

**テストの置き場所.** `test/fermion/saveload.cpp`(新規)。`test/test_fermion_layer.cpp` の末尾、
`#include "fermion/mf_measure.cpp"` の**直後**に `#include "fermion/saveload.cpp"` を 1 行足す。
それ以外のファイルは変更しない。保存ディレクトリは
`/private/tmp/claude-501/-Users-yomichi-source-github-com-issp-center-dev-TeNeS/57ccf338-dfff-45c3-812f-37392110a70f/scratchpad/`
配下ではなく、テストの実行時 CWD 直下に `output_test_fermion_saveload_<ケース名>/` を作って使い、
**ケースの冒頭で消してから作る**(既存テストが `params.outdir` に使っている流儀と同じ)。

**保存側の振る舞い.**
- `finfo.enabled` のとき、`save_tensors()` は `<save_dir>/fermion.dat` を追加で書く。
  `mpirank == 0` のみが書く。
- 形式(`params.dat` と同じ「値 + `# コメント`」):

```
1 # Fermion_Format_Version
4 # N_UNIT
2 2 # L_sub
0 # skew
0 1 # parity of the physical leg of Tn[0]
0 1 # parity of the virtual leg 0 of Tn[0]
0 1 # parity of the virtual leg 1 of Tn[0]
0 1 # parity of the virtual leg 2 of Tn[0]
0 1 # parity of the virtual leg 3 of Tn[0]
0 1 # parity of the physical leg of Tn[1]
...
```

  パリティは `0`(偶)/`1`(奇)を空白区切りで脚の次元の個数だけ。サイトは 0 から順、
  各サイトは物理脚 → 仮想脚 0,1,2,3 の 5 行。
- ボソン模式(`fermion = false`)では `fermion.dat` を**書かない**。

**読み込み側の振る舞い.** `load_tensors()` は次の順で動く:
1. 前段の検証(下表 V1〜V6b, V8)を行い、通れば `finfo.virt` を `fermion.dat` の内容で置き換える
2. 既存のテンソル読み込み(`load_tensors_v0` / `v1`)
3. 後段の検証(V7)

| # | 条件 | 例外 |
|---|---|---|
| V1 | フェルミオン模式で `fermion.dat` が無い | `tenes::load_error` |
| V2 | 先頭の形式バージョンが 1 以外 | 同 |
| V3 | `N_UNIT` が入力と不一致 | 同 |
| V4 | 物理脚パリティが入力の `parity` と不一致(長さまたは内容) | 同 |
| V5 | 仮想脚パリティの長さが入力の `virtual_dim` と不一致 | 同 |
| V6a | 復元後の `validate_neighbor_consistency` が失敗 | 同(型は問わない。`std::runtime_error` 由来でもよい) |
| V6b | `L_sub` または `skew` が入力と不一致 | `tenes::load_error` |
| V7 | 読み込んだ `Tn[i]` が復元後の台帳のもとでパリティを破っている | 同 |
| V8 | ボソン模式(`fermion = false`)なのに `fermion.dat` がある | 同 |

**テストすべきこと.**

層1(台帳の往復、これが本命):
1. 2×2 セル・`physical_dim = 2`・`virtual_dim = 4`(全脚)・`parity = {false, true}` の
   フェルミオン `iTPS` を作り、`iTPSTestAccessor::finfo` で `virt` を**既定と異なる分割**
   `{false, false, false, true}`(偶3・奇1。既定 `even_first_parity(4)` は偶2・奇2)に
   全サイト・全脚で置き換える。
2. `Tn` を、その分割のもとで**パリティ偶の要素だけ非零**な決定論テンソルで埋める
   (`tenes::fermion::count_odd(parity, idx) % 2 == 0` の要素にだけ値を入れる。値は
   `mf_measure.cpp` の決定論式でも `std::mt19937` でもよい)。λ も脚ごとに違う値にする。
3. `params.tensor_save_dir` を設定して `save_tensors()`。
4. 別インスタンス(同じ入力パラメータ、`tensor_load_dir` を設定)を構築する。
   構築時に `initialize_tensors()` → `load_tensors()` が走る。
5. 読み込み側の `finfo.virt` が 1 で設定した分割と**要素まで**一致すること。
   `finfo.phys` も一致すること。
   **既定と同じ分割だけで試すと、台帳を読まない実装でも通ってしまう**(空洞)。

層2(測定値の往復):
6. 層1 と同じ状態で、保存側と読み込み側の `measure_onesite()` / `measure_twosite()` の
   値と norm が 1e-12 で一致すること。演算子は 1 サイト(対称な実行列)と
   2 サイトのホッピング(`(0,1,1,0) = (1,0,0,1) = -1`)を各サイトの `dx=+1` / `dy=+1` に置く。
   **`MeanField_Env = true` を使うこと**(CTM 環境テンソルは保存されても意味を持たないので、
   CTM で測ると保存・読み込みと無関係な要因が混ざる)。
7. `lambda_tensor` が保存前後で**厳密に**(`operator==` で)一致すること。
   フィクスチャの λ には有効 6 桁では往復しない値(例 `1.0/3.0`、`std::sqrt(2.0)/2.0`)を使う。
   現状の `save_tensors()` は既定精度(有効 6 桁)で λ を書くので、これは Step 5b の
   精度修正を要求するテストになる。

層3(ガード): V1〜V8 のそれぞれについて、その条件だけを満たす保存ディレクトリを作り、
読み込みが例外を投げること。作り方の指針:
- V1: ボソン模式で保存した(= `fermion.dat` が無い)ディレクトリをフェルミオン模式で読む
- V2: 正常に保存してから `fermion.dat` の 1 行目を `2` に書き換える
- V3: `N_UNIT` の行を書き換える(または `L_sub` の違う入力で読む)
- V4: 物理脚パリティの行を `1 0` に書き換える
- V5: `virtual_dim` を変えた入力で読む(**エラーになるだけでなく、エラー後に
  `Tn` / `lambda_tensor` が resize されていないことも確認する**。前段検証の存在理由)
- V6b: `L_sub` の行を `3 3` に書き換える
- V7: 台帳と矛盾する `Tn`(例: 保存後に `T_0.dat` を、パリティを破る値で上書きするか、
  保存側で `finfo.virt` だけを別の分割に差し替えてから `save_tensors()` する)を読む
- V8: フェルミオンで保存したディレクトリをボソン模式(`fermion = false`)で読む

層3b(ガード解除): `test/input.cpp` は変更しない。代わりに `saveload.cpp` の中で
`validate_fermion_constraints` を `tensor_save` / `tensor_load` 付きの入力に対して呼び、
`CHECK_NOTHROW` になること。他のガード(`Use_RSVD` など)が生きていることは既存テストが担保する。

**書いてはいけないこと.** 期待値を実装の関数(`save_fermion_parity` など)で作って比べるテスト。
参照は「保存前の値」「入力パラメータ」「例外」から取る。

### 実装(Codex に渡す。`test/` 配下は変更禁止)

- [ ] **Step 1(Claude):** テスト作成者にテストを書かせ、RED を確認する。実装前は
  層1・層2 が「`finfo.virt` が既定に戻る」ことで、層3 の大半が「例外が飛ばない」ことで、
  層3b が「ガードで例外」で赤になる。テストファイルの sha256 を
  `work/sl-task1-test-snapshot.txt` に記録する

- [ ] **Step 2(Codex):** `src/iTPS/iTPS.hpp` — `void load_tensors_v0();` の直後に追加:

```cpp
  //! save the fermionic parity ledger of the virtual bonds
  void save_fermion_parity(std::string const &save_dir) const;
  //! load and validate the fermionic parity ledger (before reading tensors)
  void load_fermion_ledger(std::string const &load_dir);
  //! check the loaded tensors against the restored ledger (after reading them)
  void validate_loaded_fermion_tensors() const;
```

- [ ] **Step 3(Codex):** `src/iTPS/saveload_tensors.cpp` — `save_tensors()` の末尾、
  `if (peps_parameters.print_level >= PrintLevel::info)` の**直前**に:

```cpp
  if (finfo.enabled) {
    // The virtual-bond parity ledger is mutable state (the simple update
    // rewrites it through svd_trunc), so it has to travel with the tensors:
    // reloading with a stale ledger changes the measured energy without any
    // error message.
    //
    // The CTM environment tensors saved above carry no information in fermion
    // mode: save_tensors() runs before measure() (main.cpp, run_groundstate)
    // and the fermionic CTM is rebuilt from scratch inside measure().
    save_fermion_parity(save_dir);
  }
```

  そして新しい関数を `save_tensors` の直後に置く:

```cpp
template <class ptensor>
void iTPS<ptensor>::save_fermion_parity(std::string const &save_dir) const {
  if (mpirank != 0) {
    return;
  }
  std::string filename = save_dir + "/fermion.dat";
  std::ofstream ofs(filename.c_str());
  constexpr int fermion_format_version = 1;
  ofs << fermion_format_version << " # Fermion_Format_Version\n";
  ofs << N_UNIT << " # N_UNIT\n";
  ofs << lattice.LX << " " << lattice.LY << " # L_sub\n";
  ofs << lattice.skew << " # skew\n";
  auto write_parity = [&ofs](tenes::fermion::parity_vector const &p) {
    for (std::size_t i = 0; i < p.size(); ++i) {
      ofs << (p[i] ? 1 : 0) << " ";
    }
  };
  for (int i = 0; i < N_UNIT; ++i) {
    write_parity(finfo.phys[i]);
    ofs << "# parity of the physical leg of Tn[" << i << "]\n";
    for (int leg = 0; leg < nleg; ++leg) {
      write_parity(finfo.virt[i][leg]);
      ofs << "# parity of the virtual leg " << leg << " of Tn[" << i << "]\n";
    }
  }
}
```

- [ ] **Step 4(Codex):** 同ファイル — `load_tensors()` を次のように書き換える
  (既存の `isdir` 検査とバージョン読み出しはそのまま):

```cpp
  bcast(tensor_format_version, 0, comm);

  load_fermion_ledger(load_dir);

  if (tensor_format_version == 0) {
    load_tensors_v0();
  } else if (tensor_format_version == 1) {
    load_tensors_v1();
  } else {
    // (既存のまま)
  }

  validate_loaded_fermion_tensors();
}
```

  そして `load_fermion_ledger` を実装する。**ファイルは rank 0 が生文字列として読み、
  そのまま bcast し、全ランクが同じ構文解析と検証を行う**(片方のランクだけが例外を
  投げると MPI がデッドロックする):

```cpp
template <class ptensor>
void iTPS<ptensor>::load_fermion_ledger(std::string const &load_dir) {
  const std::string filename = load_dir + "/fermion.dat";
  int exists = 0;
  if (mpirank == 0) {
    exists = util::path_exists(filename) ? 1 : 0;
  }
  bcast(exists, 0, comm);

  if (!finfo.enabled) {
    if (exists != 0) {
      throw tenes::load_error(
          "ERROR: " + filename +
          " exists, i.e. the saved tensors come from a fermionic run, but "
          "parameter.general.fermion is false.\n"
          "HINT: set fermion = true, or load tensors saved by a "
          "non-fermionic run.");
    }
    return;
  }
  if (exists == 0) {
    throw tenes::load_error(
        "ERROR: cannot find " + filename + ".\n"
        "The saved tensors do not carry the fermionic parity ledger of the "
        "virtual bonds, which fermion mode needs in order to interpret them.\n"
        "HINT: they were saved by a non-fermionic run, or by a version of "
        "TeNeS that could not save fermionic tensors.");
  }

  std::string content;
  if (mpirank == 0) {
    std::ifstream ifs(filename.c_str());
    std::stringstream ss;
    ss << ifs.rdbuf();
    content = ss.str();
  }
  bcast(content, 0, comm);

  std::istringstream iss(content);
  std::string line;
  auto next_line = [&iss, &line, &filename]() -> std::string {
    if (!std::getline(iss, line)) {
      throw tenes::load_error("ERROR: " + filename + " ends unexpectedly");
    }
    return util::drop_comment(line);
  };
  auto next_ints = [&next_line, &filename]() -> std::vector<int> {
    const auto words = util::split(next_line());
    std::vector<int> ret;
    ret.reserve(words.size());
    for (const auto &w : words) {
      ret.push_back(std::stoi(w));
    }
    return ret;
  };

  const int version = next_ints().at(0);
  if (version != 1) {
    std::stringstream ss;
    ss << "ERROR: " << filename << " has fermion format version " << version
       << " but this version of TeNeS supports only version 1";
    throw tenes::load_error(ss.str());
  }
  const int loaded_N_UNIT = next_ints().at(0);
  if (loaded_N_UNIT != N_UNIT) {
    std::stringstream ss;
    ss << "ERROR: N_UNIT is " << N_UNIT << " but " << filename << " has "
       << loaded_N_UNIT;
    throw tenes::load_error(ss.str());
  }
  const auto lsub = next_ints();
  const int loaded_skew = next_ints().at(0);
  if (lsub.size() != 2 || lsub[0] != lattice.LX || lsub[1] != lattice.LY ||
      loaded_skew != lattice.skew) {
    std::stringstream ss;
    ss << "ERROR: the unit cell of the saved tensors (L_sub = [";
    for (std::size_t i = 0; i < lsub.size(); ++i) {
      ss << (i == 0 ? "" : ", ") << lsub[i];
    }
    ss << "], skew = " << loaded_skew << ") differs from the input (L_sub = ["
       << lattice.LX << ", " << lattice.LY << "], skew = " << lattice.skew
       << ").\n"
       << "HINT: the parity ledger is indexed by (site, leg), so it only means "
          "the same thing on the same lattice.";
    throw tenes::load_error(ss.str());
  }

  auto to_parity = [](std::vector<int> const &v) {
    tenes::fermion::parity_vector p(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
      p[i] = (v[i] != 0);
    }
    return p;
  };

  std::vector<std::array<tenes::fermion::parity_vector, 4>> virt(N_UNIT);
  for (int i = 0; i < N_UNIT; ++i) {
    const auto phys = to_parity(next_ints());
    if (phys != finfo.phys[i]) {
      std::stringstream ss;
      ss << "ERROR: the physical parity of the tensor " << i << " in "
         << filename << " differs from tensor.unitcell.parity in the input";
      throw tenes::load_error(ss.str());
    }
    for (int leg = 0; leg < nleg; ++leg) {
      auto p = to_parity(next_ints());
      if (static_cast<int>(p.size()) != lattice.virtual_dims[i][leg]) {
        std::stringstream ss;
        ss << "ERROR: the virtual dimension of the leg " << leg
           << " of the tensor " << i << " is " << lattice.virtual_dims[i][leg]
           << " but the saved tensors have " << p.size() << ".\n"
           << "HINT: fermion mode cannot change virtual_dim on restart, "
              "because resizing a leg does not preserve its even/odd blocks. "
              "Keep virtual_dim as it was, or start without tensor_load.";
        throw tenes::load_error(ss.str());
      }
      virt[i][leg] = p;
    }
  }
  for (int i = 0; i < N_UNIT; ++i) {
    finfo.virt[i] = virt[i];
  }
  tenes::fermion::validate_neighbor_consistency(finfo, lattice);
}
```

- [ ] **Step 5(Codex):** 同ファイル — 後段の検証:

```cpp
template <class ptensor>
void iTPS<ptensor>::validate_loaded_fermion_tensors() const {
  if (!finfo.enabled) {
    return;
  }
  for (int i = 0; i < N_UNIT; ++i) {
    const auto ft = tenes::fermion::wrap_Tn(Tn[i], finfo, i);
    // parity_violation and max_abs scan the process-local slice only.
    std::vector<double> reduced{tenes::fermion::parity_violation(ft),
                                tenes::fermion::max_abs(ft)};
    tenes::allreduce_max(reduced, comm);
    const double violation = reduced[0];
    const double scale = std::max(reduced[1], 1.0e-300);
    if (violation > 1.0e-12 * scale) {
      std::stringstream ss;
      ss << "ERROR: the loaded tensor " << i
         << " breaks fermion parity under the loaded parity ledger "
         << "(max violating amplitude " << violation << ", max amplitude "
         << reduced[1] << ").\n"
         << "HINT: the tensors and " << peps_parameters.tensor_load_dir
         << "/fermion.dat do not belong together.";
      throw tenes::load_error(ss.str());
    }
  }
}
```

  必要な `#include`(`<sstream>`、`<array>`、`"../fermion/fops.hpp"`、
  `"../fermion/fermion_info.hpp"`)を確認して足す。

- [ ] **Step 5b(Codex):** 同ファイル — `save_tensors()` の λ 書き出しを厳密往復にする。
  現状は `std::ofstream` の既定精度(有効 6 桁)なので、保存・読み込みで λ が丸められる。
  `density.cpp:122` と同じ書式指定を足す:

```cpp
      std::ofstream ofs(save_dir + "/lambda_" + std::to_string(i) + ".dat");
      // max_digits10 round-trips a double exactly; the default 6 significant
      // digits silently truncated the Schmidt weights on every checkpoint.
      ofs << std::scientific
          << std::setprecision(std::numeric_limits<double>::max_digits10);
```

  `<iomanip>` と `<limits>` の include を確認して足す。これはボソン経路にも効くが、
  保存ファイルの**精度が上がるだけ**で読み出し側(`ifs >> temp`)は変更不要、
  古い保存ファイルもそのまま読める。`restart` ctest は 2 回の実行を rtol 1e-3 で
  比べているので、精度向上で悪化することはない。

- [ ] **Step 6(Codex):** `src/iTPS/load_toml.cpp` — 次のブロックを削除する:

```cpp
  if (!peps_parameters.tensor_load_dir.empty() ||
      !peps_parameters.tensor_save_dir.empty()) {
    throw_fermion_guard("tensor load/save directories");
  }
```

- [ ] **Step 7(Codex):** `cmake --build out-gcc/build -j 8` → `out-gcc/build/test/test_fermion_layer`
  と `out-gcc/build/test/test_input` が全件緑。報告を `work/sl-task1-report.md` に

- [ ] **Step 8(Claude):** sha256 でテスト不改変を確認。`ctest --test-dir out-gcc/build` **全件**
  (`restart` を含む)。`clang-format -i` を変更した src とテストに適用。
  変異テスト(scratchpad のコピーで): (a) `finfo.virt[i] = virt[i]` の代入を消す → 層1・層2 が赤、
  (b) V7 の閾値比較を `false` に固定 → 層3 の V7 が赤、(c) `load_fermion_ledger` の呼び出しを
  テンソル読み込みの**後**に移す → 層3 の V5(resize されていないこと)が赤。
  コミット `"Persist the fermionic parity ledger with the saved tensors"`

---

## Task 2: E2E テスト、ドキュメント、NEWS

**Files:**
- Create: `test/fermion/free_fermion_saveload.py.in`
- Modify: `test/CMakeLists.txt`(`FreeFermionMF` の登録の直後)
- Modify: `docs/sphinx/ja/file_specification/parameter_section.rst`
- Modify: `docs/sphinx/en/file_specification/parameter_section.rst`
- Modify: `NEWS.md`

### 契約書(E2E、テスト作成者に渡す散文)

`test/fermion/free_fermion_mf.py.in` を雛形に `test/fermion/free_fermion_saveload.py.in` を書く
(必要な関数はコピーしてよい。共有モジュール化はしない)。spinless free fermion、μ = 0、
D = 2、`meanfield_env = true`(CTM を作らないぶん速い)、`seed = 11`。

1. **run1**: `num_step = 200`、`tensor_save = "save1"` で実行。終了コード 0。
2. `save1/fermion.dat` が存在し、次を満たすこと:
   - 1 行目の形式バージョンが 1、2 行目が `N_UNIT`(= 4)
   - `L_sub` 行が `2 2`、`skew` 行が `0`
   - 続く 4×5 = 20 行のパリティ行があり、物理脚の行は長さ 2、仮想脚の行は長さ 2
   - 各行の要素が `0` か `1` のみ
   **`even_first_parity` と違う分割であることを要求してはいけない**(既定の初期化では
   一致するのが正常。設計書 §2.4)。
3. **run2**: 同じ入力から `num_step = 0`、`tensor_load = "save1"` で実行。終了コード 0。
   `output/density.dat` が run1 と `np.isclose(rtol=1e-10, atol=1e-12)` で一致すること。
4. **run3**(スキャン運用の煙テスト): run1 の保存を読み、μ を −1 に変えたハミルトニアンと
   発展演算子で `num_step = 200` 実行。終了コード 0、`Energy` が有限で負、
   `n` が 0 < n < 1 であること。
5. **run4**(D 変更の拒否): run1 の保存を `virtual_dim = [3,3,3,3]` の入力で読み、
   **終了コードが 0 でない**こと、標準出力または標準エラーに `virtual_dim` を含む
   メッセージが出ること。
6. golden ファイルは使わない(run1 と run2 の比較が自己完結しているため)。
   run1 の `Energy` / `n` を実行ログに印字する。

`test/CMakeLists.txt`: `configure_file` + `add_test(NAME FreeFermionSaveLoad ...)` +
`set_tests_properties(FreeFermionSaveLoad PROPERTIES TIMEOUT 1800)` を
`FreeFermionMF` の登録の直後に置く。

### ドキュメント・NEWS(Codex に渡す。`test/` 配下は変更禁止)

- [ ] **Step 1(Claude):** E2E をテスト作成者に書かせ、`ctest -R FreeFermionSaveLoad` が
  緑になることを確認(Task 1 が入っているので RED にはならない。ここは回帰網)
- [ ] **Step 2(Codex):** `docs/sphinx/ja/file_specification/parameter_section.rst` の
  フェルミオン非対応一覧から `tensor_save` / `tensor_load` を外し、直後に次の項を足す:
  「``tensor_save`` / ``tensor_load`` はフェルミオン模式でも使えます。仮想ボンドの
  偶奇台帳が保存先の ``fermion.dat`` に書き出され、読み込み時に物理脚のパリティ・
  ``virtual_dim``・``L_sub``・``skew``・テンソル自身のパリティが検証されます。
  ``virtual_dim`` を変えての読み込みは非対応です(偶奇ブロックの構造が保てないため)」。
  英語版も同内容
- [ ] **Step 3(Codex):** `NEWS.md` のフェルミオン非対応一覧から `tensor_save`/`tensor_load` を
  外し、項目を足す: "`tensor_save` / `tensor_load` now work in fermion mode, so a parameter
  scan can restart from a converged state. The parity ledger of the virtual bonds is saved
  alongside the tensors as `fermion.dat` and validated on load (physical parity, `virtual_dim`,
  `L_sub`, `skew`, and the parity of the tensors themselves); restarting with a different
  `virtual_dim` is not supported."。報告を `work/sl-task2-report.md` に
- [ ] **Step 4(Claude):** `ctest --test-dir out-gcc/build` 全件緑。docs の RST 構文確認。
  コミット `"Add the fermionic save/load end-to-end test and document it"`

---

## 最終レビュー

- [ ] レビューサブエージェントにブランチ全体をレビューさせる。観点: (1) ボソン経路の
  到達条件が変わっていないか、(2) MPI で片ランクだけが例外を投げる経路が無いか、
  (3) 検証の前後分割が実際に resize より前になっているか、(4) 変異テストの結果、
  (5) エラーメッセージが利用者に何をすればよいか伝えているか
