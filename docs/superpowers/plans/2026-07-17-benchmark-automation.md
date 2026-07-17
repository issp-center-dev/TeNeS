# ベンチマーク自動化 実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 変更前後の性能を A/B 比較できる汎用計測基盤(C++ タイマーレジストリ + Python ベンチマークハーネス)を作る。

**Architecture:** C++ 側に階層名(`contract/itps_ctm/2x2` 等)で引ける累積タイマーレジストリを追加し、実行終了時に `output/timers.json` を書く。Python 側は `benchmark/bench.py` が TOML スイート定義からケースを展開して `tenes_simple → tenes_std → tenes` を反復実行し、結果を名前マッチングで比較して Markdown レポートを生成する。

**Tech Stack:** C++17 / doctest / CMake / Python 3(numpy、tomllib または toml)/ pytest

**Spec:** `docs/superpowers/specs/2026-07-17-benchmark-automation-design.md`

## Global Constraints

- C++17(`CMAKE_CXX_STANDARD 17`)。C++20 機能は使わない。
- すべての C++ 実装は namespace `tenes` 内(ソルバー核は `tenes::itps`)。
- 新規ソースファイルの先頭には GPL ヘッダを付ける。C++ 用:

```cpp
/* TeNeS - Massively parallel tensor network solver /
/ Copyright (C) 2019- The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify /
/ it under the terms of the GNU General Public License as published by /
/ the Free Software Foundation, either version 3 of the License, or /
/ (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, /
/ but WITHOUT ANY WARRANTY; without even the implied warranty of /
/ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the /
/ GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License /
/ along with this program. If not, see http://www.gnu.org/licenses/. */
```

Python 用(既存 `test/python/test_tenes_std.py` と同一):

```python
# TeNeS - Massively parallel tensor network solver
# Copyright (C) 2019- The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses
```

以下のコードブロック中の `/* GPLヘッダ */` および `# GPLヘッダ` は、上記ヘッダをそのまま貼ることを意味する(この参照のみ例外的に許すプレースホルダ)。

- コミット前フォーマット: C++ は `clang-format -i <file>`、Python は `python3 -m black <files>`(line-length 88 = black デフォルト)。
- ビルドは既存の `build/` ディレクトリを使う: `cmake --build build -j`。存在しなければ `mkdir build && cd build && cmake ../ && cd ..`。
- テスト実行は `(cd build && ctest -R <name>)`。
- タイマー名の規約: スラッシュ区切りの階層名。`total`、`phase/<名前>`、`contract/<バックエンド>/<NxM>`。名前に使う文字は `[a-z0-9_/x]` のみ。
- ベンチ結果ディレクトリ `benchmark/results/` は git 管理外(Task 9 で .gitignore に追加)。

---

### Task 1: TimerRegistry と ScopedTimer

**Files:**
- Modify: `src/timer.hpp`(`Timer` クラスの下に追加)
- Create: `test/timer_registry.cpp`
- Modify: `test/CMakeLists.txt:9-10`(foreach リストに `timer_registry` を追加)

**Interfaces:**
- Consumes: 既存 `tenes::Timer<>`(`src/timer.hpp`、`elapsed()` が秒を返す)
- Produces: `tenes::TimerRegistry`(`void add(std::string const&, double)`、`std::map<std::string, Entry> const& entries() const`、`static TimerRegistry& instance()`、`Entry = {std::size_t count; double sum;}`)、`tenes::ScopedTimer(std::string name, TimerRegistry& = instance())`

- [ ] **Step 1: 失敗するテストを書く**

`test/timer_registry.cpp` を作成:

```cpp
/* GPLヘッダ */

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <string>

#include "../src/mpi.hpp"
#include "../src/timer.hpp"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  doctest::Context context(argc, argv);
  const int res = context.run();
  MPI_Finalize();
  return res;
}

TEST_CASE("TimerRegistry::add accumulates count and sum") {
  tenes::TimerRegistry registry;
  registry.add("contract/itps_ctm/2x2", 1.5);
  registry.add("contract/itps_ctm/2x2", 0.5);
  registry.add("phase/simple_update", 2.0);

  auto const &entries = registry.entries();
  REQUIRE(entries.size() == 2);
  CHECK(entries.at("contract/itps_ctm/2x2").count == 2);
  CHECK(entries.at("contract/itps_ctm/2x2").sum == doctest::Approx(2.0));
  CHECK(entries.at("phase/simple_update").count == 1);
  CHECK(entries.at("phase/simple_update").sum == doctest::Approx(2.0));
}

TEST_CASE("ScopedTimer records on destruction") {
  tenes::TimerRegistry registry;
  {
    tenes::ScopedTimer timer("scope/test", registry);
  }
  auto const &entries = registry.entries();
  REQUIRE(entries.count("scope/test") == 1);
  CHECK(entries.at("scope/test").count == 1);
  CHECK(entries.at("scope/test").sum >= 0.0);
}

TEST_CASE("TimerRegistry::instance returns a singleton") {
  auto &a = tenes::TimerRegistry::instance();
  auto &b = tenes::TimerRegistry::instance();
  CHECK(&a == &b);
}
```

注意: 既存 `test/util.cpp` の冒頭(doctest の include と `main`)と同じ流儀にする。`util.cpp` の `#include "doctest.h"` より前にある `#define`(`DOCTEST_CONFIG_IMPLEMENT` のはず)を確認し、異なればそちらに合わせる。

`test/CMakeLists.txt` の foreach リストを変更:

```cmake
foreach(basename input simple_update full_update tensor_util util
        correlation_length saveload arnoldi timer_registry)
```

- [ ] **Step 2: テストが失敗する(コンパイルエラーになる)ことを確認**

Run: `cmake --build build -j 2>&1 | tail -20`
Expected: `timer_registry.cpp` のコンパイルエラー(`TimerRegistry is not a member of 'tenes'` 相当)。
(初回は `cd build && cmake ../` で再構成が必要な場合がある)

- [ ] **Step 3: 実装を書く**

`src/timer.hpp` の include を拡張し、`Timer` クラス定義の直後(namespace `tenes` 内)に追加:

```cpp
// includes(ファイル先頭、<chrono> の行を置き換え):
#include <chrono>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
```

```cpp
/*! @brief Process-wide accumulator of named timers
 *
 * Names are hierarchical, slash-separated (e.g. "contract/itps_ctm/2x2").
 */
class TimerRegistry {
 public:
  struct Entry {
    std::size_t count = 0;
    double sum = 0.0;  //!< accumulated time in seconds
  };

  void add(std::string const &name, double seconds) {
    Entry &e = entries_[name];
    e.count += 1;
    e.sum += seconds;
  }

  std::map<std::string, Entry> const &entries() const { return entries_; }

  //! Registry used by instrumentation points across the process.
  static TimerRegistry &instance() {
    static TimerRegistry registry;
    return registry;
  }

 private:
  std::map<std::string, Entry> entries_;
};

//! RAII helper: measures from construction to destruction.
class ScopedTimer {
 public:
  explicit ScopedTimer(std::string name,
                       TimerRegistry &registry = TimerRegistry::instance())
      : name_(std::move(name)), registry_(registry) {}
  ~ScopedTimer() { registry_.add(name_, timer_.elapsed()); }
  ScopedTimer(ScopedTimer const &) = delete;
  ScopedTimer &operator=(ScopedTimer const &) = delete;

 private:
  std::string name_;
  TimerRegistry &registry_;
  Timer<> timer_;
};
```

- [ ] **Step 4: テストが通ることを確認**

Run: `cmake --build build -j && (cd build && ctest -R timer_registry --output-on-failure)`
Expected: `100% tests passed, 0 tests failed out of 1`

- [ ] **Step 5: フォーマットとコミット**

```bash
clang-format -i src/timer.hpp test/timer_registry.cpp
cmake --build build -j && (cd build && ctest -R timer_registry --output-on-failure)
git add src/timer.hpp test/timer_registry.cpp test/CMakeLists.txt
git commit -m "Add TimerRegistry and ScopedTimer for named cumulative timing"
```

---

### Task 2: MPI 集計と JSON 出力

**Files:**
- Modify: `src/mpi.hpp`(`allreduce_sum` 群の直後に `allreduce_max` / `allreduce_min` を追加)
- Modify: `src/timer.hpp`(宣言追加)
- Create: `src/timer.cpp`
- Modify: `src/CMakeLists.txt`(`add_library(tenes_mpi STATIC mpi.cpp)` → `add_library(tenes_mpi STATIC mpi.cpp timer.cpp)`)
- Modify: `test/timer_registry.cpp`(テスト追加)

**Interfaces:**
- Consumes: Task 1 の `TimerRegistry`、既存 `tenes::allreduce_sum`(`src/mpi.hpp:150` 付近のパターン)、`MPI_Comm`(`src/mpi.hpp` が `_NO_MPI` スタブを提供)
- Produces:
  - `tenes::TimerAggregate {std::size_t count; double sum; double max_rank; double min_rank;}`
  - `std::map<std::string, TimerAggregate> tenes::aggregate_timers(TimerRegistry const&, MPI_Comm)`(集団呼び出し)
  - `std::string tenes::timers_to_json(std::map<std::string, TimerAggregate> const&, std::string const& tenes_version, int mpi_size, int omp_threads)`
  - `template<class T> int tenes::allreduce_max(std::vector<T>&, MPI_Comm)` / `allreduce_min`(`_NO_MPI` では恒等変換)

- [ ] **Step 1: 失敗するテストを書く**

`test/timer_registry.cpp` の末尾に追加:

```cpp
TEST_CASE("aggregate_timers with a single process") {
  tenes::TimerRegistry registry;
  registry.add("a", 1.0);
  registry.add("a", 2.0);
  registry.add("b", 4.0);

  auto agg = tenes::aggregate_timers(registry, MPI_COMM_WORLD);
  REQUIRE(agg.size() == 2);
  CHECK(agg.at("a").count == 2);
  CHECK(agg.at("a").sum == doctest::Approx(3.0));
  CHECK(agg.at("a").max_rank == doctest::Approx(3.0));
  CHECK(agg.at("a").min_rank == doctest::Approx(3.0));
  CHECK(agg.at("b").sum == doctest::Approx(4.0));
}

TEST_CASE("timers_to_json renders the expected fields") {
  std::map<std::string, tenes::TimerAggregate> timers;
  timers["total"] = tenes::TimerAggregate{1, 1.5, 1.5, 1.5};
  timers["contract/itps_ctm/2x2"] = tenes::TimerAggregate{48, 4.5, 5.0, 4.0};
  const std::string json = tenes::timers_to_json(timers, "2.2-dev", 1, 8);
  CHECK(json.find("\"tenes_version\": \"2.2-dev\"") != std::string::npos);
  CHECK(json.find("\"mpi_size\": 1") != std::string::npos);
  CHECK(json.find("\"omp_threads\": 8") != std::string::npos);
  CHECK(json.find("\"total\": {\"count\": 1, \"sum\": 1.5") !=
        std::string::npos);
  CHECK(json.find("\"contract/itps_ctm/2x2\": {\"count\": 48") !=
        std::string::npos);
}
```

- [ ] **Step 2: コンパイルエラーで失敗することを確認**

Run: `cmake --build build -j 2>&1 | tail -10`
Expected: `aggregate_timers` / `TimerAggregate` 未定義のコンパイルエラー。

- [ ] **Step 3: 実装を書く**

`src/mpi.hpp` の `allreduce_sum(std::vector<T>&, MPI_Comm)` の直後に追加(同じパターン):

```cpp
template <class T>
int allreduce_max(std::vector<T> &val, MPI_Comm comm) {
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  int N = val.size();
  std::vector<T> recv(val);
  int ret = MPI_Allreduce(val.data(), recv.data(), N, datatype, MPI_MAX, comm);
  if (ret != 0) {
    return ret;
  }
  val.assign(recv.begin(), recv.end());
#endif
  return 0;
}

template <class T>
int allreduce_min(std::vector<T> &val, MPI_Comm comm) {
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  int N = val.size();
  std::vector<T> recv(val);
  int ret = MPI_Allreduce(val.data(), recv.data(), N, datatype, MPI_MIN, comm);
  if (ret != 0) {
    return ret;
  }
  val.assign(recv.begin(), recv.end());
#endif
  return 0;
}
```

`src/timer.hpp` の `ScopedTimer` の後(namespace `tenes` 内)に宣言を追加。`timer.hpp` の include に `#include "mpi.hpp"` を足す:

```cpp
//! Cross-rank view of one timer (times in seconds).
struct TimerAggregate {
  std::size_t count = 0;
  double sum = 0.0;       //!< value on the calling rank
  double max_rank = 0.0;  //!< maximum over MPI ranks
  double min_rank = 0.0;  //!< minimum over MPI ranks
};

/*! @brief Aggregate timer sums across MPI ranks (collective call)
 *
 * Assumes every rank recorded the same set of names
 * (tensor operations are collective).
 */
std::map<std::string, TimerAggregate> aggregate_timers(
    TimerRegistry const &registry, MPI_Comm comm);

//! Render aggregated timers and metadata as a JSON document.
std::string timers_to_json(std::map<std::string, TimerAggregate> const &timers,
                           std::string const &tenes_version, int mpi_size,
                           int omp_threads);
```

`src/timer.cpp` を作成:

```cpp
/* GPLヘッダ */

#include "timer.hpp"

#include <iomanip>
#include <sstream>
#include <vector>

#include "mpi.hpp"

namespace tenes {

std::map<std::string, TimerAggregate> aggregate_timers(
    TimerRegistry const &registry, MPI_Comm comm) {
  std::vector<double> sums;
  sums.reserve(registry.entries().size());
  for (auto const &kv : registry.entries()) {
    sums.push_back(kv.second.sum);
  }
  std::vector<double> maxs = sums;
  std::vector<double> mins = sums;
  allreduce_max(maxs, comm);
  allreduce_min(mins, comm);

  std::map<std::string, TimerAggregate> result;
  std::size_t i = 0;
  for (auto const &kv : registry.entries()) {
    result[kv.first] =
        TimerAggregate{kv.second.count, sums[i], maxs[i], mins[i]};
    ++i;
  }
  return result;
}

namespace {
std::string escape_json(std::string const &s) {
  std::string ret;
  for (char c : s) {
    if (c == '"' || c == '\\') {
      ret += '\\';
    }
    ret += c;
  }
  return ret;
}
}  // namespace

std::string timers_to_json(std::map<std::string, TimerAggregate> const &timers,
                           std::string const &tenes_version, int mpi_size,
                           int omp_threads) {
  std::ostringstream oss;
  oss << std::setprecision(17);
  oss << "{\n";
  oss << "  \"meta\": {\n";
  oss << "    \"tenes_version\": \"" << escape_json(tenes_version) << "\",\n";
  oss << "    \"mpi_size\": " << mpi_size << ",\n";
  oss << "    \"omp_threads\": " << omp_threads << "\n";
  oss << "  },\n";
  oss << "  \"timers\": {\n";
  bool first = true;
  for (auto const &kv : timers) {
    if (!first) {
      oss << ",\n";
    }
    first = false;
    oss << "    \"" << escape_json(kv.first) << "\": {"
        << "\"count\": " << kv.second.count << ", \"sum\": " << kv.second.sum
        << ", \"max_rank\": " << kv.second.max_rank
        << ", \"min_rank\": " << kv.second.min_rank << "}";
  }
  oss << "\n  }\n}\n";
  return oss.str();
}

}  // namespace tenes
```

`src/CMakeLists.txt` の変更:

```cmake
add_library(tenes_mpi STATIC mpi.cpp timer.cpp)
```

- [ ] **Step 4: テストが通ることを確認**

Run: `cmake --build build -j && (cd build && ctest -R timer_registry --output-on-failure)`
Expected: PASS

- [ ] **Step 5: フォーマットとコミット**

```bash
clang-format -i src/timer.hpp src/timer.cpp src/mpi.hpp test/timer_registry.cpp
cmake --build build -j && (cd build && ctest -R timer_registry --output-on-failure)
git add src/timer.hpp src/timer.cpp src/mpi.hpp src/CMakeLists.txt test/timer_registry.cpp
git commit -m "Add MPI aggregation and JSON serialization for TimerRegistry"
```

---

### Task 3: iTPS への統合 — timers.json の出力

**Files:**
- Modify: `src/iTPS/measure.cpp`(`iTPS<ptensor>::summary()`、80行目付近)

**Interfaces:**
- Consumes: Task 1-2 の `TimerRegistry::instance()` / `aggregate_timers` / `timers_to_json`、iTPS のメンバー `timer_all`, `time_simple_update`, `time_full_update`, `time_environment`, `time_observable`(`src/iTPS/iTPS.hpp:325-329`)、`mpirank`, `outdir`、`TENES_VERSION`(`src/version.hpp`)
- Produces: `tenes` 実行終了時に `<outdir>/timers.json`(`total`, `phase/*` エントリと meta を含む)

- [ ] **Step 1: iTPS のメンバー名を確認**

Run: `grep -n "MPI_Comm\|mpisize\|mpirank\|outdir" src/iTPS/iTPS.hpp | head -10`
Expected: `MPI_Comm comm`、`int mpisize`、`int mpirank`、`std::string outdir` 相当のメンバーが見つかる。名前が違う場合は以降のコードをその名前に合わせる。

- [ ] **Step 2: summary() を変更**

`src/iTPS/measure.cpp` の include 部に追加:

```cpp
#include "../timer.hpp"
#include "../version.hpp"
#ifndef _NO_OMP
#include <omp.h>
#endif
```

`summary()`(measure.cpp 80行目付近)を次のように変更。`if (mpirank == 0)` の**前**にレジストリ登録と集計(集団呼び出しのため全 rank で実行)を入れ、rank 0 ブロック内の time.dat 書き出しの後に timers.json 書き出しを追加する:

```cpp
template <class ptensor>
void iTPS<ptensor>::summary() const {
  auto &registry = TimerRegistry::instance();
  registry.add("total", timer_all.elapsed());
  registry.add("phase/simple_update", time_simple_update);
  registry.add("phase/full_update", time_full_update);
  registry.add("phase/environment", time_environment);
  registry.add("phase/observable", time_observable);
#ifndef _NO_OMP
  const int omp_threads = omp_get_max_threads();
#else
  const int omp_threads = 1;
#endif
  const auto aggregated = aggregate_timers(registry, comm);
  if (mpirank == 0) {
    const double time_all = timer_all.elapsed();
    {
      std::string filename = outdir + "/time.dat";
      std::ofstream ofs(filename.c_str());
      ofs << "time all           = " << time_all << std::endl;
      ofs << "time simple update = " << time_simple_update << std::endl;
      ofs << "time full update   = " << time_full_update << std::endl;
      ofs << "time environment  = " << time_environment << std::endl;
      ofs << "time observable    = " << time_observable << std::endl;
      if (peps_parameters.print_level >= PrintLevel::info) {
        std::cout << "    Save elapsed times to " << filename << std::endl;
      }
    }
    {
      std::string filename = outdir + "/timers.json";
      std::ofstream ofs(filename.c_str());
      ofs << timers_to_json(aggregated, TENES_VERSION, mpisize, omp_threads);
      if (peps_parameters.print_level >= PrintLevel::info) {
        std::cout << "    Save timers to " << filename << std::endl;
      }
    }
    if (peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "Wall times [sec.]:" << std::endl;
      std::cout << "  all           = " << time_all << std::endl;
      std::cout << "  simple update = " << time_simple_update << std::endl;
      std::cout << "  full update   = " << time_full_update << std::endl;
      std::cout << "  environment  = " << time_environment << std::endl;
      std::cout << "  observable    = " << time_observable << std::endl;
      std::cout << std::endl << "Done." << std::endl;
    }
  }
}
```

注: `phase/*` の `count` は 1 になる(summary で1回だけ add するため)。意味を持つのは `sum` のみ。この仕様は Task 10 の README に明記する。

- [ ] **Step 3: ビルドして E2E で確認**

Run:
```bash
cmake --build build -j
(cd build && ctest -R Honeycomb --output-on-failure)
find build/test -name timers.json | head -3
python3 -m json.tool "$(find build/test -name timers.json | head -1)" | head -20
```
Expected: ctest PASS。timers.json が存在し、JSON としてパースでき、`"total"` と `"phase/simple_update"` などのキーを含む。

- [ ] **Step 4: フォーマットとコミット**

```bash
clang-format -i src/iTPS/measure.cpp
cmake --build build -j && (cd build && ctest -R Honeycomb --output-on-failure)
git add src/iTPS/measure.cpp
git commit -m "Write output/timers.json with total and phase timings"
```

---

### Task 4: 縮約カーネルディスパッチャの計測

**Files:**
- Modify: `src/iTPS/core/contract_itps_ctm/ctm.cpp:44-85`(`Contract_iTPS_CTM`)
- Modify: `src/iTPS/core/contract_itps_mf/mf.cpp:44` 付近(`Contract_iTPS_MF`)
- Modify: `src/iTPS/core/contract_density_ctm/ctm.cpp:44` 付近(`Contract_density_CTM`)

**Interfaces:**
- Consumes: Task 1 の `tenes::ScopedTimer`(ディスパッチャは `tenes::itps::core` 内なので非修飾で参照可)
- Produces: timers.json に `contract/itps_ctm/<N>x<M>`、`contract/itps_mf/<N>x<M>`、`contract/density_ctm/<N>x<M>` エントリ

- [ ] **Step 1: 計測点が逐次コンテキストであることを確認**

Run: `grep -rn "pragma omp" src/iTPS/*.cpp src/iTPS/core/ | grep -v Binary`
Expected: `src/iTPS/transfer_matrix.cpp` の2箇所(要素代入ループ)のみ。縮約ディスパッチャの呼び出し元(`onesite_obs.cpp` / `twosite_obs.cpp` / `multisite_obs.cpp` 等)に OpenMP 並列ループはない。もし新たな並列呼び出しが見つかった場合は実装を止めてユーザーに報告する。

- [ ] **Step 2: 3つのディスパッチャに ScopedTimer を入れる**

`src/iTPS/core/contract_itps_ctm/ctm.cpp` — include 部に追加:

```cpp
#include <string>

#include "../../../timer.hpp"
```

`Contract_iTPS_CTM` の `nrow` / `ncol` 定義の直後(`#define CALL_CONTRACT` の前)に挿入:

```cpp
  const size_t nrow = Tn.size();
  const size_t ncol = Tn[0].size();
  ScopedTimer scoped_timer("contract/itps_ctm/" + std::to_string(nrow) + "x" +
                           std::to_string(ncol));
```

`src/iTPS/core/contract_itps_mf/mf.cpp` の `Contract_iTPS_MF` にも同様に(名前は `"contract/itps_mf/"`)、`src/iTPS/core/contract_density_ctm/ctm.cpp` の `Contract_density_CTM` にも同様に(名前は `"contract/density_ctm/"`)挿入する。3ファイルとも `nrow`/`ncol` を `Tn` から計算する同じ構造をしている(異なる場合はその関数の行数変数に合わせる)。

- [ ] **Step 3: ビルドして 3 バックエンドの計測を確認**

Run:
```bash
cmake --build build -j
(cd build && ctest -R "Honeycomb$|AntiferroHeisenberg_mf|FT_Kitaev" --output-on-failure)
for f in $(find build/test -name timers.json); do echo "== $f"; grep -o '"contract/[a-z_]*' "$f" | sort -u; done
```
Expected: ctest PASS(3件)。Honeycomb の timers.json に `"contract/itps_ctm`、AntiferroHeisenberg_mf に `"contract/itps_mf`、FT_Kitaev に `"contract/density_ctm` が現れる。

- [ ] **Step 4: 回帰がないことを全テストで確認**

Run: `(cd build && ctest --output-on-failure)`
Expected: 全テスト PASS(このブランチ作成時点で 23 テスト)。

- [ ] **Step 5: フォーマットとコミット**

```bash
clang-format -i src/iTPS/core/contract_itps_ctm/ctm.cpp src/iTPS/core/contract_itps_mf/mf.cpp src/iTPS/core/contract_density_ctm/ctm.cpp
cmake --build build -j
git add src/iTPS/core/contract_itps_ctm/ctm.cpp src/iTPS/core/contract_itps_mf/mf.cpp src/iTPS/core/contract_density_ctm/ctm.cpp
git commit -m "Instrument contraction dispatchers with per-cluster timers"
```

---

### Task 5: benchlib.suite — スイート定義の読み込みとスイープ展開

**Files:**
- Create: `benchmark/benchlib/__init__.py`
- Create: `benchmark/benchlib/suite.py`
- Create: `test/python/test_benchmark_suite.py`

**Interfaces:**
- Consumes: なし(純粋 Python)
- Produces(`benchlib.suite`):
  - `load_toml(path) -> dict`(tomllib、なければ toml にフォールバック)
  - `@dataclass Case(name: str, kind: str, source: Path, params: dict)`(kind は `"template"` か `"input"`)
  - `@dataclass Suite(name: str, repeat: int, cases: list)`
  - `expand_sweep(sweep: dict) -> list[dict]`(直積展開)
  - `derive_params(params: dict) -> dict`(`Lsub=[L,W]` → `L`,`W`、`chi_ratio=r` → `chi=r*D*D`)
  - `format_name(pattern: str, params: dict) -> str`(`string.Template` の `${var}` 置換)
  - `render_template(text: str, params: dict) -> str`(同上)
  - `load_suite(path) -> Suite`

- [ ] **Step 1: 失敗するテストを書く**

`test/python/test_benchmark_suite.py` を作成:

```python
# GPLヘッダ

import os
import sys

import pytest

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmark"),
)

from benchlib import suite


def test_expand_sweep_product():
    got = suite.expand_sweep({"D": [2, 4], "chi_ratio": [1, 2]})
    assert len(got) == 4
    assert {"D": 2, "chi_ratio": 1} in got
    assert {"D": 4, "chi_ratio": 2} in got


def test_expand_sweep_with_list_values():
    got = suite.expand_sweep({"Lsub": [[1, 1], [2, 2]]})
    assert got == [{"Lsub": [1, 1]}, {"Lsub": [2, 2]}]


def test_derive_params_lsub_and_chi():
    got = suite.derive_params({"Lsub": [2, 3], "D": 4, "chi_ratio": 2})
    assert got["L"] == 2
    assert got["W"] == 3
    assert got["chi"] == 32
    assert "Lsub" not in got
    assert "chi_ratio" not in got


def test_format_name():
    name = suite.format_name(
        "afh_${L}x${W}_D${D}_chi${chi}", {"L": 2, "W": 2, "D": 4, "chi": 32}
    )
    assert name == "afh_2x2_D4_chi32"


def test_load_suite_expands_template_cases(tmp_path):
    (tmp_path / "tpl.toml").write_text("D = ${D}\nsteps = ${steps}\n")
    suite_toml = tmp_path / "suite.toml"
    suite_toml.write_text(
        """
[suite]
name = "s"
repeat = 2

[[case]]
name = "c_D${D}"
template = "tpl.toml"
params = { steps = 10 }
sweep = { D = [2, 4] }
"""
    )
    s = suite.load_suite(suite_toml)
    assert s.name == "s"
    assert s.repeat == 2
    assert [c.name for c in s.cases] == ["c_D2", "c_D4"]
    assert s.cases[0].kind == "template"
    assert s.cases[0].params == {"steps": 10, "D": 2}
    assert s.cases[0].source == tmp_path / "tpl.toml"


def test_load_suite_input_case(tmp_path):
    (tmp_path / "in.toml").write_text("")
    suite_toml = tmp_path / "suite.toml"
    suite_toml.write_text(
        """
[suite]
name = "s"

[[case]]
name = "existing"
input = "in.toml"
"""
    )
    s = suite.load_suite(suite_toml)
    assert s.repeat == 1
    assert s.cases[0].kind == "input"
    assert s.cases[0].source == tmp_path / "in.toml"


def test_load_suite_rejects_duplicate_names(tmp_path):
    (tmp_path / "tpl.toml").write_text("")
    suite_toml = tmp_path / "suite.toml"
    suite_toml.write_text(
        """
[suite]
name = "s"

[[case]]
name = "same"
template = "tpl.toml"

[[case]]
name = "same"
template = "tpl.toml"
"""
    )
    with pytest.raises(ValueError):
        suite.load_suite(suite_toml)


def test_render_template():
    assert suite.render_template("L = ${L}\n", {"L": 2}) == "L = 2\n"
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `python3 -m pytest test/python/test_benchmark_suite.py -v 2>&1 | tail -5`
Expected: `ModuleNotFoundError: No module named 'benchlib'` で collection エラー。

- [ ] **Step 3: 実装を書く**

`benchmark/benchlib/__init__.py`:

```python
# GPLヘッダ
```

(ライセンスヘッダのみの空モジュール)

`benchmark/benchlib/suite.py`:

```python
# GPLヘッダ

import itertools
import string
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List


def load_toml(path):
    """Load a TOML file with tomllib (Python >= 3.11) or toml as fallback."""
    try:
        import tomllib
    except ImportError:
        import toml

        return toml.load(str(path))
    with open(path, "rb") as f:
        return tomllib.load(f)


@dataclass
class Case:
    name: str
    kind: str  # "template" or "input"
    source: Path
    params: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Suite:
    name: str
    repeat: int
    cases: List[Case]


def expand_sweep(sweep: Dict[str, List[Any]]) -> List[Dict[str, Any]]:
    """Expand {key: [values...]} into the direct product as a list of dicts."""
    keys = list(sweep.keys())
    values = [sweep[k] for k in keys]
    return [dict(zip(keys, comb)) for comb in itertools.product(*values)]


def derive_params(params: Dict[str, Any]) -> Dict[str, Any]:
    """Resolve derived parameters.

    - Lsub = [L, W]   -> L, W
    - chi_ratio = r   -> chi = r * D * D (requires D)
    """
    ret = dict(params)
    if "Lsub" in ret:
        ret["L"], ret["W"] = ret.pop("Lsub")
    if "chi_ratio" in ret:
        ret["chi"] = ret.pop("chi_ratio") * ret["D"] ** 2
    return ret


def _substitute(text: str, params: Dict[str, Any]) -> str:
    return string.Template(text).substitute({k: str(v) for k, v in params.items()})


def format_name(pattern: str, params: Dict[str, Any]) -> str:
    return _substitute(pattern, params)


def render_template(text: str, params: Dict[str, Any]) -> str:
    return _substitute(text, params)


def load_suite(path) -> Suite:
    path = Path(path)
    data = load_toml(path)
    suite_table = data["suite"]
    cases: List[Case] = []
    for case_table in data.get("case", []):
        if "template" in case_table:
            base = dict(case_table.get("params", {}))
            sweep = case_table.get("sweep", {})
            param_sets = expand_sweep(sweep) if sweep else [{}]
            for sweep_params in param_sets:
                params = derive_params({**base, **sweep_params})
                cases.append(
                    Case(
                        name=format_name(case_table["name"], params),
                        kind="template",
                        source=path.parent / case_table["template"],
                        params=params,
                    )
                )
        elif "input" in case_table:
            cases.append(
                Case(
                    name=case_table["name"],
                    kind="input",
                    source=path.parent / case_table["input"],
                )
            )
        else:
            raise ValueError("each [[case]] must have either 'template' or 'input'")
    names = [c.name for c in cases]
    if len(names) != len(set(names)):
        raise ValueError("case names must be unique after sweep expansion")
    return Suite(
        name=suite_table["name"],
        repeat=int(suite_table.get("repeat", 1)),
        cases=cases,
    )
```

- [ ] **Step 4: テストが通ることを確認**

Run: `python3 -m pytest test/python/test_benchmark_suite.py -v`
Expected: 8 passed

- [ ] **Step 5: フォーマットとコミット**

```bash
python3 -m black benchmark/benchlib test/python/test_benchmark_suite.py
python3 -m pytest test/python/test_benchmark_suite.py -q
git add benchmark/benchlib/__init__.py benchmark/benchlib/suite.py test/python/test_benchmark_suite.py
git commit -m "Add benchmark suite definition loader with sweep expansion"
```

---

### Task 6: benchlib.meta と benchlib.runner — メタデータ収集と実行

**Files:**
- Create: `benchmark/benchlib/meta.py`
- Create: `benchmark/benchlib/runner.py`
- Create: `test/python/test_benchmark_runner.py`

**Interfaces:**
- Consumes: Task 5 の `benchlib.suite`(`Case`, `Suite`, `render_template`)
- Produces:
  - `benchlib.meta.collect_meta(tenes_bin, repo_dir, launcher=None) -> dict`(キー: `date`, `hostname`, `omp_num_threads`, `launcher`, `git_commit`, `git_dirty`, `tenes_version`)
  - `benchlib.runner.RunContext(tenes_dir, tool_dir, results_dir, launcher=None, repeat=None)`
  - `benchlib.runner.prepare_input(case, workdir, ctx) -> Path`(workdir に input.toml を生成)
  - `benchlib.runner.run_case(case, repeat, ctx)`(`results_dir/<case>/run_<i>/` に成果物を保存)
  - `benchlib.runner.run_suite(suite, ctx)`(meta.json を書き、全ケースを実行)

- [ ] **Step 1: 失敗するテストを書く**

`test/python/test_benchmark_runner.py` を作成:

```python
# GPLヘッダ

import json
import os
import sys

import pytest

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmark"),
)

from benchlib import meta, runner, suite


def _make_stub(path, body):
    path.write_text("#!/bin/sh\n" + body + "\n")
    path.chmod(0o755)


@pytest.fixture
def stub_tools(tmp_path):
    """tenes_simple/tenes_std/tenes の代わりになるスタブ実行ファイル群"""
    tool_dir = tmp_path / "tool"
    tool_dir.mkdir()
    _make_stub(tool_dir / "tenes_simple", "touch std.toml")
    _make_stub(tool_dir / "tenes_std", "touch input.toml")
    tenes_dir = tmp_path / "src"
    tenes_dir.mkdir()
    _make_stub(
        tenes_dir / "tenes",
        'mkdir -p output && echo "{}" > output/timers.json'
        ' && echo "0.5" > output/density.dat',
    )
    return tenes_dir, tool_dir


def test_collect_meta_in_git_repo(tmp_path):
    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
    got = meta.collect_meta(tmp_path / "no_such_tenes", repo_dir=repo)
    assert len(got["git_commit"]) == 40
    assert isinstance(got["git_dirty"], bool)
    assert got["tenes_version"] is None  # 存在しないバイナリは None
    assert got["hostname"]


def test_prepare_input_renders_template(tmp_path, stub_tools):
    tenes_dir, tool_dir = stub_tools
    tpl = tmp_path / "tpl.toml"
    tpl.write_text("D = ${D}\n")
    case = suite.Case(name="c", kind="template", source=tpl, params={"D": 2})
    ctx = runner.RunContext(
        tenes_dir=tenes_dir, tool_dir=tool_dir, results_dir=tmp_path / "results"
    )
    workdir = tmp_path / "results" / "c" / "work"
    input_toml = runner.prepare_input(case, workdir, ctx)
    assert (workdir / "simple.toml").read_text() == "D = 2\n"
    assert input_toml.exists()


def test_run_case_stores_outputs_per_run(tmp_path, stub_tools):
    tenes_dir, tool_dir = stub_tools
    tpl = tmp_path / "tpl.toml"
    tpl.write_text("D = ${D}\n")
    case = suite.Case(name="c_D2", kind="template", source=tpl, params={"D": 2})
    ctx = runner.RunContext(
        tenes_dir=tenes_dir, tool_dir=tool_dir, results_dir=tmp_path / "results"
    )
    runner.run_case(case, repeat=2, ctx=ctx)
    for i in range(2):
        rundir = tmp_path / "results" / "c_D2" / "run_{}".format(i)
        assert (rundir / "timers.json").exists()
        assert (rundir / "density.dat").exists()


def test_run_suite_writes_meta(tmp_path, stub_tools):
    tenes_dir, tool_dir = stub_tools
    tpl = tmp_path / "tpl.toml"
    tpl.write_text("")
    s = suite.Suite(
        name="s",
        repeat=1,
        cases=[suite.Case(name="c", kind="template", source=tpl, params={})],
    )
    ctx = runner.RunContext(
        tenes_dir=tenes_dir, tool_dir=tool_dir, results_dir=tmp_path / "results"
    )
    runner.run_suite(s, ctx)
    got = json.loads((tmp_path / "results" / "meta.json").read_text())
    assert got["suite"] == "s"
    assert "git_commit" in got
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `python3 -m pytest test/python/test_benchmark_runner.py -v 2>&1 | tail -5`
Expected: `ImportError: cannot import name 'meta'`(または runner)で collection エラー。

- [ ] **Step 3: 実装を書く**

`benchmark/benchlib/meta.py`:

```python
# GPLヘッダ

import os
import socket
import subprocess
from datetime import datetime, timezone


def _run(cmd, cwd=None):
    return subprocess.run(
        cmd, cwd=cwd, check=True, capture_output=True, text=True
    ).stdout.strip()


def collect_meta(tenes_bin, repo_dir, launcher=None):
    """Collect provenance metadata for a benchmark run."""
    result = {
        "date": datetime.now(timezone.utc).astimezone().isoformat(),
        "hostname": socket.gethostname(),
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
        "launcher": launcher,
    }
    try:
        result["git_commit"] = _run(["git", "rev-parse", "HEAD"], cwd=repo_dir)
        result["git_dirty"] = bool(
            _run(["git", "status", "--porcelain"], cwd=repo_dir)
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        result["git_commit"] = None
        result["git_dirty"] = None
    try:
        result["tenes_version"] = _run([str(tenes_bin), "--version"])
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        result["tenes_version"] = None
    return result
```

`benchmark/benchlib/runner.py`:

```python
# GPLヘッダ

import json
import shutil
import subprocess
from pathlib import Path

from . import meta as meta_mod
from . import suite as suite_mod


class RunContext:
    def __init__(self, tenes_dir, tool_dir, results_dir, launcher=None, repeat=None):
        self.tenes_dir = Path(tenes_dir)
        self.tool_dir = Path(tool_dir)
        self.results_dir = Path(results_dir)
        self.launcher = launcher
        self.repeat = repeat  # None -> suite default


def _call(cmd, cwd):
    cmd = [str(c) for c in cmd]
    print("  $ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def prepare_input(case, workdir, ctx):
    """Generate input.toml in workdir via tenes_simple -> tenes_std."""
    workdir.mkdir(parents=True, exist_ok=True)
    simple_toml = workdir / "simple.toml"
    if case.kind == "template":
        simple_toml.write_text(
            suite_mod.render_template(case.source.read_text(), case.params)
        )
    else:
        shutil.copy(case.source, simple_toml)
    _call([ctx.tool_dir / "tenes_simple", "simple.toml"], cwd=workdir)
    _call([ctx.tool_dir / "tenes_std", "std.toml"], cwd=workdir)
    return workdir / "input.toml"


def run_case(case, repeat, ctx):
    """Run one case `repeat` times and store outputs under run_<i>/."""
    casedir = ctx.results_dir / case.name
    workdir = casedir / "work"
    prepare_input(case, workdir, ctx)
    tenes_cmd = []
    if ctx.launcher:
        tenes_cmd += ctx.launcher.split()
    tenes_cmd += [ctx.tenes_dir / "tenes", "input.toml"]
    for i in range(repeat):
        _call(tenes_cmd, cwd=workdir)
        rundir = casedir / "run_{}".format(i)
        rundir.mkdir(parents=True, exist_ok=True)
        outdir = workdir / "output"
        for f in sorted(outdir.glob("*.dat")) + sorted(outdir.glob("*.json")):
            shutil.copy(f, rundir / f.name)


def run_suite(s, ctx):
    ctx.results_dir.mkdir(parents=True, exist_ok=True)
    repo_dir = Path(__file__).resolve().parents[2]
    result = meta_mod.collect_meta(
        ctx.tenes_dir / "tenes", repo_dir=repo_dir, launcher=ctx.launcher
    )
    result["suite"] = s.name
    (ctx.results_dir / "meta.json").write_text(json.dumps(result, indent=2))
    repeat = ctx.repeat if ctx.repeat is not None else s.repeat
    for case in s.cases:
        print("== case: {} (repeat={}) ==".format(case.name, repeat), flush=True)
        run_case(case, repeat, ctx)
```

- [ ] **Step 4: テストが通ることを確認**

Run: `python3 -m pytest test/python/test_benchmark_runner.py -v`
Expected: 4 passed

- [ ] **Step 5: フォーマットとコミット**

```bash
python3 -m black benchmark/benchlib test/python/test_benchmark_runner.py
python3 -m pytest test/python/test_benchmark_runner.py test/python/test_benchmark_suite.py -q
git add benchmark/benchlib/meta.py benchmark/benchlib/runner.py test/python/test_benchmark_runner.py
git commit -m "Add benchmark runner: pipeline execution, repeats, metadata"
```

---

### Task 7: benchlib.stats と benchlib.obscheck — 統計集計と物理量チェック

**Files:**
- Create: `benchmark/benchlib/stats.py`
- Create: `benchmark/benchlib/obscheck.py`
- Create: `test/python/test_benchmark_stats.py`

**Interfaces:**
- Consumes: なし(純粋 Python + numpy)
- Produces:
  - `benchlib.stats.aggregate_runs(run_timers: list[dict]) -> dict`。入力は各反復の timers 辞書(`name -> {"count": int, "sum": float, ...}`)のリスト。出力は `name -> {"count": int, "count_varies": bool, "median": float, "min": float, "max": float}`
  - `benchlib.obscheck.extract_numbers(text: str) -> list[float]`(`#` 始まり行を無視して数値トークンを抽出)
  - `benchlib.obscheck.compare_dat_dirs(dir_a, dir_b, rtol=1e-3, atol=1e-4, exclude=("time.dat", "parameters.dat")) -> list[str]`(警告メッセージのリスト、一致なら空)

- [ ] **Step 1: 失敗するテストを書く**

`test/python/test_benchmark_stats.py` を作成:

```python
# GPLヘッダ

import os
import sys

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmark"),
)

from benchlib import obscheck, stats


def _t(count, s):
    return {"count": count, "sum": s}


def test_aggregate_runs_median_min_max():
    runs = [
        {"a": _t(10, 3.0)},
        {"a": _t(10, 1.0)},
        {"a": _t(10, 2.0)},
    ]
    got = stats.aggregate_runs(runs)
    assert got["a"]["median"] == 2.0
    assert got["a"]["min"] == 1.0
    assert got["a"]["max"] == 3.0
    assert got["a"]["count"] == 10
    assert got["a"]["count_varies"] is False


def test_aggregate_runs_missing_name_and_count_change():
    runs = [
        {"a": _t(10, 1.0), "b": _t(1, 5.0)},
        {"a": _t(20, 2.0)},
    ]
    got = stats.aggregate_runs(runs)
    assert got["a"]["count_varies"] is True
    assert got["b"]["median"] == 5.0


def test_extract_numbers_skips_comments():
    text = "# header 123\n1.0 2.5e-3\n-4 hello5\n"
    assert obscheck.extract_numbers(text) == [1.0, 2.5e-3, -4.0, 5.0]


def test_compare_dat_dirs(tmp_path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir()
    b.mkdir()
    (a / "density.dat").write_text("1.0 2.0\n")
    (b / "density.dat").write_text("1.0001 2.0\n")
    (a / "time.dat").write_text("999\n")
    (b / "time.dat").write_text("1\n")
    assert obscheck.compare_dat_dirs(a, b) == []

    (b / "density.dat").write_text("1.5 2.0\n")
    warnings = obscheck.compare_dat_dirs(a, b)
    assert len(warnings) == 1
    assert "density.dat" in warnings[0]

    (b / "extra.dat").write_text("0\n")
    warnings = obscheck.compare_dat_dirs(a, b)
    assert any("extra.dat" in w for w in warnings)
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `python3 -m pytest test/python/test_benchmark_stats.py -v 2>&1 | tail -5`
Expected: `ImportError` で collection エラー。

- [ ] **Step 3: 実装を書く**

`benchmark/benchlib/stats.py`:

```python
# GPLヘッダ

import statistics


def aggregate_runs(run_timers):
    """Aggregate repeated runs of one case.

    run_timers: list of timer dicts (name -> {"count": int, "sum": float, ...}).
    Returns name -> {"count", "count_varies", "median", "min", "max"}.
    Names missing from some runs are aggregated over the runs that have them.
    """
    names = []
    for timers in run_timers:
        for name in timers:
            if name not in names:
                names.append(name)
    result = {}
    for name in names:
        sums = [t[name]["sum"] for t in run_timers if name in t]
        counts = [t[name]["count"] for t in run_timers if name in t]
        result[name] = {
            "count": counts[0],
            "count_varies": len(set(counts)) > 1 or len(sums) != len(run_timers),
            "median": statistics.median(sums),
            "min": min(sums),
            "max": max(sums),
        }
    return result
```

`benchmark/benchlib/obscheck.py`:

```python
# GPLヘッダ

import re
from pathlib import Path

import numpy as np

_FLOAT_RE = re.compile(r"[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?")


def extract_numbers(text):
    """Extract numeric tokens, ignoring comment lines starting with '#'."""
    numbers = []
    for line in text.splitlines():
        if line.lstrip().startswith("#"):
            continue
        for m in _FLOAT_RE.finditer(line):
            numbers.append(float(m.group(0)))
    return numbers


def compare_dat_dirs(
    dir_a, dir_b, rtol=1e-3, atol=1e-4, exclude=("time.dat", "parameters.dat")
):
    """Numerically compare common .dat files. Returns a list of warnings."""
    dir_a = Path(dir_a)
    dir_b = Path(dir_b)
    files_a = {p.name for p in dir_a.glob("*.dat")} - set(exclude)
    files_b = {p.name for p in dir_b.glob("*.dat")} - set(exclude)
    warnings = []
    for name in sorted(files_a & files_b):
        num_a = extract_numbers((dir_a / name).read_text())
        num_b = extract_numbers((dir_b / name).read_text())
        if len(num_a) != len(num_b):
            warnings.append(
                "{}: number of values differs ({} vs {})".format(
                    name, len(num_a), len(num_b)
                )
            )
        elif not np.allclose(num_a, num_b, rtol=rtol, atol=atol):
            warnings.append("{}: values differ beyond tolerance".format(name))
    for name in sorted(files_a ^ files_b):
        warnings.append("{}: missing on one side".format(name))
    return warnings
```

- [ ] **Step 4: テストが通ることを確認**

Run: `python3 -m pytest test/python/test_benchmark_stats.py -v`
Expected: 4 passed

- [ ] **Step 5: フォーマットとコミット**

```bash
python3 -m black benchmark/benchlib test/python/test_benchmark_stats.py
python3 -m pytest test/python/test_benchmark_stats.py -q
git add benchmark/benchlib/stats.py benchmark/benchlib/obscheck.py test/python/test_benchmark_stats.py
git commit -m "Add benchmark statistics aggregation and observable check"
```

---

### Task 8: benchlib.compare — 比較レポート生成

**Files:**
- Create: `benchmark/benchlib/compare.py`
- Create: `test/python/test_benchmark_compare.py`

**Interfaces:**
- Consumes: Task 7 の `stats.aggregate_runs` / `obscheck.compare_dat_dirs`
- Produces:
  - `benchlib.compare.load_label_dir(label_dir) -> (meta: dict, cases: dict)`。`cases` は `case_name -> {"agg": <aggregate_runs結果>, "run0": Path}`
  - `benchlib.compare.compare_timers(agg_a, agg_b) -> list[dict]`(行: `{"name", "a", "b", "ratio", "overlap", "count_mismatch"}`。片側欠損時は `a` か `b` が None)
  - `benchlib.compare.compare_results(dir_a, dir_b, rtol=1e-3, atol=1e-4) -> str`(Markdown レポート)

- [ ] **Step 1: 失敗するテストを書く**

`test/python/test_benchmark_compare.py` を作成:

```python
# GPLヘッダ

import json
import os
import sys

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmark"),
)

from benchlib import compare


def _write_label(root, label, sums, count=10):
    """テスト用のラベルディレクトリ(1ケース・2反復)を作る"""
    labeldir = root / label
    for i, s in enumerate(sums):
        rundir = labeldir / "case1" / "run_{}".format(i)
        rundir.mkdir(parents=True)
        timers = {
            "total": {"count": 1, "sum": 10 * s, "max_rank": 10 * s, "min_rank": 10 * s},
            "phase/environment": {"count": 1, "sum": 5 * s, "max_rank": 5 * s, "min_rank": 5 * s},
            "contract/itps_ctm/2x2": {"count": count, "sum": s, "max_rank": s, "min_rank": s},
        }
        (rundir / "timers.json").write_text(
            json.dumps({"meta": {}, "timers": timers})
        )
        (rundir / "density.dat").write_text("1.0\n")
    (labeldir / "meta.json").write_text(
        json.dumps({"git_commit": label + "0" * 34, "hostname": "h"})
    )
    return labeldir


def test_compare_timers_ratio_and_overlap():
    agg_a = {"x": {"count": 1, "count_varies": False, "median": 2.0, "min": 1.9, "max": 2.1}}
    agg_b = {"x": {"count": 1, "count_varies": False, "median": 1.0, "min": 0.9, "max": 1.1}}
    rows = compare.compare_timers(agg_a, agg_b)
    assert rows[0]["ratio"] == 0.5
    assert rows[0]["overlap"] is False
    assert rows[0]["count_mismatch"] is False


def test_compare_timers_one_side_only():
    rows = compare.compare_timers({"only_a": {"count": 1, "count_varies": False, "median": 1.0, "min": 1.0, "max": 1.0}}, {})
    assert rows[0]["a"] is not None
    assert rows[0]["b"] is None


def test_compare_results_end_to_end(tmp_path):
    dir_a = _write_label(tmp_path, "base", [1.0, 1.2])
    dir_b = _write_label(tmp_path, "new", [0.5, 0.6])
    report = compare.compare_results(dir_a, dir_b)
    assert "# Benchmark comparison" in report
    assert "total" in report
    assert "## contract" in report
    assert "case1" in report
    assert "WARNING" not in report


def test_compare_results_detects_observable_mismatch(tmp_path):
    dir_a = _write_label(tmp_path, "base", [1.0])
    dir_b = _write_label(tmp_path, "new", [1.0])
    rundir = dir_b / "case1" / "run_0"
    (rundir / "density.dat").write_text("2.0\n")
    report = compare.compare_results(dir_a, dir_b)
    assert "WARNING" in report
    assert "density.dat" in report


def test_compare_results_flags_count_mismatch(tmp_path):
    dir_a = _write_label(tmp_path, "base", [1.0], count=10)
    dir_b = _write_label(tmp_path, "new", [1.0], count=20)
    report = compare.compare_results(dir_a, dir_b)
    assert "count differs" in report
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `python3 -m pytest test/python/test_benchmark_compare.py -v 2>&1 | tail -5`
Expected: `ImportError` で collection エラー。

- [ ] **Step 3: 実装を書く**

`benchmark/benchlib/compare.py`:

```python
# GPLヘッダ

import json
from pathlib import Path

from . import obscheck, stats

_META_KEYS = (
    "git_commit",
    "git_dirty",
    "tenes_version",
    "hostname",
    "date",
    "omp_num_threads",
    "launcher",
    "suite",
)


def load_label_dir(label_dir):
    """Load meta.json and per-case aggregated timers from a label directory."""
    label_dir = Path(label_dir)
    meta_path = label_dir / "meta.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    cases = {}
    for casedir in sorted(p for p in label_dir.iterdir() if p.is_dir()):
        runs = []
        for rundir in sorted(casedir.glob("run_*")):
            timers_path = rundir / "timers.json"
            if timers_path.exists():
                runs.append(json.loads(timers_path.read_text())["timers"])
        if runs:
            cases[casedir.name] = {
                "agg": stats.aggregate_runs(runs),
                "run0": casedir / "run_0",
            }
    return meta, cases


def compare_timers(agg_a, agg_b):
    """Compare two aggregate dicts name-by-name."""
    rows = []
    for name in sorted(set(agg_a) | set(agg_b)):
        a = agg_a.get(name)
        b = agg_b.get(name)
        row = {"name": name, "a": a, "b": b}
        if a is not None and b is not None:
            row["ratio"] = b["median"] / a["median"] if a["median"] > 0 else None
            row["overlap"] = a["min"] <= b["max"] and b["min"] <= a["max"]
            row["count_mismatch"] = a["count"] != b["count"]
        rows.append(row)
    return rows


def _group_key(name):
    """Group timers: total and phase/* into 'summary', others by prefix."""
    if name == "total" or name.startswith("phase/"):
        return "summary"
    return name.split("/", 1)[0] if "/" in name else "other"


def _stat_cell(s):
    return "{:.4g} [{:.4g}, {:.4g}]".format(s["median"], s["min"], s["max"])


def _format_row(case, row):
    a, b = row["a"], row["b"]
    if a is None:
        return "| {} | {} | (absent) | {} | n/a |  | B only |".format(
            case, row["name"], _stat_cell(b)
        )
    if b is None:
        return "| {} | {} | {} | (absent) | n/a |  | A only |".format(
            case, row["name"], _stat_cell(a)
        )
    ratio = "{:.3f}".format(row["ratio"]) if row["ratio"] is not None else "n/a"
    notes = []
    if row["overlap"]:
        notes.append("no sig. diff")
    if row["count_mismatch"]:
        notes.append("count differs")
    if a.get("count_varies") or b.get("count_varies"):
        notes.append("count varies across runs")
    return "| {} | {} | {} | {} | {} | {}/{} | {} |".format(
        case,
        row["name"],
        _stat_cell(a),
        _stat_cell(b),
        ratio,
        a["count"],
        b["count"],
        ", ".join(notes),
    )


_TABLE_HEADER = (
    "| case | timer | A: median [min, max] | B: median [min, max]"
    " | ratio B/A | count A/B | note |\n|---|---|---|---|---|---|---|"
)


def compare_results(dir_a, dir_b, rtol=1e-3, atol=1e-4):
    """Render a Markdown report comparing two label directories."""
    meta_a, cases_a = load_label_dir(dir_a)
    meta_b, cases_b = load_label_dir(dir_b)
    lines = []
    lines.append(
        "# Benchmark comparison: {} (A) vs {} (B)".format(
            Path(dir_a).name, Path(dir_b).name
        )
    )
    lines.append("")
    lines.append("| meta | A | B |")
    lines.append("|---|---|---|")
    for key in _META_KEYS:
        lines.append("| {} | {} | {} |".format(key, meta_a.get(key), meta_b.get(key)))
    lines.append("")

    common = sorted(set(cases_a) & set(cases_b))
    only = sorted(set(cases_a) ^ set(cases_b))
    if only:
        lines.append("**WARNING: cases on one side only:** " + ", ".join(only))
        lines.append("")

    obs_warnings = []
    for case in common:
        for w in obscheck.compare_dat_dirs(
            cases_a[case]["run0"], cases_b[case]["run0"], rtol=rtol, atol=atol
        ):
            obs_warnings.append("{}: {}".format(case, w))
    if obs_warnings:
        lines.append("## WARNING: observable mismatches")
        lines.append("")
        for w in obs_warnings:
            lines.append("- " + w)
        lines.append("")

    groups = {}
    for case in common:
        rows = compare_timers(cases_a[case]["agg"], cases_b[case]["agg"])
        for row in rows:
            groups.setdefault(_group_key(row["name"]), []).append((case, row))

    group_names = sorted(groups, key=lambda g: (g != "summary", g))
    for group in group_names:
        title = "Summary (total / phase)" if group == "summary" else group + "/"
        lines.append("## " + title)
        lines.append("")
        lines.append(_TABLE_HEADER)
        for case, row in groups[group]:
            lines.append(_format_row(case, row))
        lines.append("")
    return "\n".join(lines)
```

- [ ] **Step 4: テストが通ることを確認**

Run: `python3 -m pytest test/python/test_benchmark_compare.py -v`
Expected: 5 passed

- [ ] **Step 5: フォーマットとコミット**

```bash
python3 -m black benchmark/benchlib test/python/test_benchmark_compare.py
python3 -m pytest test/python/test_benchmark_compare.py -q
git add benchmark/benchlib/compare.py test/python/test_benchmark_compare.py
git commit -m "Add benchmark comparison report generator"
```

---

### Task 9: bench.py CLI とスイート定義ファイル

**Files:**
- Create: `benchmark/bench.py`
- Create: `benchmark/templates/afh_square.toml`
- Create: `benchmark/suites/contraction.toml`
- Create: `benchmark/suites/e2e.toml`
- Create: `benchmark/suites/smoke.toml`
- Create: `benchmark/suites/ci.toml`
- Modify: `.gitignore`(`benchmark/results/` を追加。ファイルが無ければ作成)

**Interfaces:**
- Consumes: Task 5-8 の `benchlib.{suite,runner,compare}`
- Produces: CLI `bench.py run <suite.toml> --label L [--tenes-dir D] [--tool-dir D] [--results-dir D] [--launcher CMD] [--repeat N]` / `bench.py compare <dirA> <dirB> [--output FILE] [--rtol R] [--atol A]`

- [ ] **Step 1: bench.py を書く**

`benchmark/bench.py`:

```python
#!/usr/bin/env python3
# GPLヘッダ

"""TeNeS benchmark harness.

Run a suite and store results:
    bench.py run suites/contraction.toml --label baseline

Compare two result sets:
    bench.py compare results/baseline results/cotengra
"""

import argparse
import sys
from pathlib import Path

BENCH_DIR = Path(__file__).resolve().parent
REPO_DIR = BENCH_DIR.parent
sys.path.insert(0, str(BENCH_DIR))

from benchlib import compare as compare_mod
from benchlib import runner, suite


def cmd_run(args):
    s = suite.load_suite(args.suite)
    results_dir = Path(args.results_dir) / args.label
    if results_dir.exists():
        sys.exit("error: {} already exists; use a new label".format(results_dir))
    ctx = runner.RunContext(
        tenes_dir=args.tenes_dir,
        tool_dir=args.tool_dir,
        results_dir=results_dir,
        launcher=args.launcher,
        repeat=args.repeat,
    )
    runner.run_suite(s, ctx)
    print("results saved to {}".format(results_dir))


def cmd_compare(args):
    report = compare_mod.compare_results(
        args.dir_a, args.dir_b, rtol=args.rtol, atol=args.atol
    )
    if args.output:
        Path(args.output).write_text(report)
        print("report saved to {}".format(args.output))
    else:
        print(report)


def main():
    parser = argparse.ArgumentParser(description="TeNeS benchmark harness")
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser("run", help="run a benchmark suite")
    p_run.add_argument("suite", help="path to a suite TOML file")
    p_run.add_argument("--label", required=True, help="identifier of this run")
    p_run.add_argument("--tenes-dir", default=str(REPO_DIR / "build" / "src"))
    p_run.add_argument("--tool-dir", default=str(REPO_DIR / "build" / "tool"))
    p_run.add_argument("--results-dir", default=str(BENCH_DIR / "results"))
    p_run.add_argument("--launcher", default=None, help='e.g. "mpirun -np 4"')
    p_run.add_argument("--repeat", type=int, default=None, help="override suite repeat")
    p_run.set_defaults(func=cmd_run)

    p_cmp = sub.add_parser("compare", help="compare two result directories")
    p_cmp.add_argument("dir_a")
    p_cmp.add_argument("dir_b")
    p_cmp.add_argument("--output", default=None, help="write report to a file")
    p_cmp.add_argument("--rtol", type=float, default=1e-3)
    p_cmp.add_argument("--atol", type=float, default=1e-4)
    p_cmp.set_defaults(func=cmd_compare)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
```

作成後: `chmod +x benchmark/bench.py`

- [ ] **Step 2: テンプレートとスイート定義を書く**

`benchmark/templates/afh_square.toml`(`sample/AFH_square/simple.toml` がベース。χ は `parameter.ctm.dimension`):

```toml
# Antiferromagnetic Heisenberg model on the square lattice.
# Placeholders: ${simple_steps} ${chi} ${L} ${W} ${D}
[parameter]
[parameter.general]
output = "output"

[parameter.simple_update]
num_step = ${simple_steps}
tau = 0.01

[parameter.full_update]
num_step = 0
tau = 0.01

[parameter.ctm]
dimension = ${chi}
iteration_max = 10
meanfield_env = false

[lattice]
type = "square lattice"
L = ${L}
W = ${W}
initial = "antiferro"
virtual_dim = ${D}

[model]
type = "spin"
J = 1.0
```

`benchmark/suites/contraction.toml`:

```toml
# Unit-cell size x D x chi sweep for contraction-kernel benchmarking.
[suite]
name = "contraction"
repeat = 3

[[case]]
name = "afh_square_${L}x${W}_D${D}_chi${chi}"
template = "../templates/afh_square.toml"
params = { simple_steps = 100 }
sweep = { Lsub = [[1, 1], [2, 2], [3, 3], [4, 4]], D = [2, 3], chi_ratio = [2] }
```

`benchmark/suites/e2e.toml`:

```toml
# Representative end-to-end workloads based on existing inputs.
[suite]
name = "e2e"
repeat = 3

[[case]]
name = "honeycomb"
input = "../../test/data/Honeycomb.toml"

[[case]]
name = "afh_square_sample"
input = "../../sample/AFH_square/simple.toml"
```

`benchmark/suites/smoke.toml`:

```toml
# Minimal suite for CI smoke testing (a few seconds).
[suite]
name = "smoke"
repeat = 1

[[case]]
name = "smoke_${L}x${W}_D${D}_chi${chi}"
template = "../templates/afh_square.toml"
params = { simple_steps = 10 }
sweep = { Lsub = [[1, 1]], D = [2], chi_ratio = [2] }
```

`benchmark/suites/ci.toml`:

```toml
# Reduced suite for GitHub Actions (reference values only).
[suite]
name = "ci"
repeat = 3

[[case]]
name = "afh_square_${L}x${W}_D${D}_chi${chi}"
template = "../templates/afh_square.toml"
params = { simple_steps = 50 }
sweep = { Lsub = [[1, 1], [2, 2]], D = [2], chi_ratio = [2] }
```

`.gitignore` に追記:

```
benchmark/results/
```

- [ ] **Step 3: 手動で動作確認(スモークスイート)**

Run:
```bash
python3 benchmark/bench.py run benchmark/suites/smoke.toml --label check_a
python3 benchmark/bench.py run benchmark/suites/smoke.toml --label check_b
python3 benchmark/bench.py compare benchmark/results/check_a benchmark/results/check_b | head -40
```
Expected: 両ラベルが正常終了し、レポートに meta 表・`Summary (total / phase)`・`contract/` セクションが出る。observable の WARNING は出ない(同一バイナリなので)。確認後 `rm -rf benchmark/results/check_a benchmark/results/check_b`。

- [ ] **Step 4: フォーマットとコミット**

```bash
python3 -m black benchmark/bench.py
python3 -m pytest test/python -q
git add benchmark/bench.py benchmark/templates benchmark/suites .gitignore
git commit -m "Add bench.py CLI, suite definitions, and templates"
```

---

### Task 10: ctest スモークテストと README

**Files:**
- Create: `test/benchmark_smoke.py.in`
- Modify: `test/CMakeLists.txt`(末尾に追加)
- Create: `benchmark/README.md`

**Interfaces:**
- Consumes: Task 9 の `bench.py` と `suites/smoke.toml`、CMake 変数 `TENES_PYTHON_EXECUTABLE` / `MPIEXEC` / `ENABLE_MPI`
- Produces: ctest エントリ `benchmark_smoke`

- [ ] **Step 1: スモークテストスクリプトを書く**

`test/benchmark_smoke.py.in`:

```python
# GPLヘッダ

import json
import os
import shutil
import subprocess
import sys
from os.path import join

SOURCE_DIR = "@CMAKE_SOURCE_DIR@"
BINARY_DIR = "@CMAKE_BINARY_DIR@"
LAUNCHER = "@BENCH_SMOKE_LAUNCHER@"

resdir = join("@CMAKE_CURRENT_BINARY_DIR@", "bench_smoke_results")
if os.path.exists(resdir):
    shutil.rmtree(resdir)

bench = join(SOURCE_DIR, "benchmark", "bench.py")
cmd = [
    sys.executable,
    bench,
    "run",
    join(SOURCE_DIR, "benchmark", "suites", "smoke.toml"),
    "--label",
    "smoke",
    "--tenes-dir",
    join(BINARY_DIR, "src"),
    "--tool-dir",
    join(BINARY_DIR, "tool"),
    "--results-dir",
    resdir,
]
if LAUNCHER:
    cmd += ["--launcher", LAUNCHER]
subprocess.run(cmd, check=True)

labeldir = join(resdir, "smoke")
casedirs = [
    d
    for d in os.listdir(labeldir)
    if os.path.isdir(join(labeldir, d))
]
assert casedirs, "no case directory produced"
timers_path = join(labeldir, casedirs[0], "run_0", "timers.json")
with open(timers_path) as f:
    data = json.load(f)
timers = data["timers"]
assert "total" in timers, "total timer missing"
assert any(
    name.startswith("contract/") for name in timers
), "no contraction timers recorded"

report = subprocess.run(
    [sys.executable, bench, "compare", labeldir, labeldir],
    check=True,
    capture_output=True,
    text=True,
).stdout
assert "total" in report
assert "WARNING" not in report, "self-comparison must not warn"

print("benchmark smoke test passed")
```

`test/CMakeLists.txt` の末尾に追加:

```cmake
if(ENABLE_MPI)
  set(BENCH_SMOKE_LAUNCHER
      "${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_POSTFLAGS}")
else()
  set(BENCH_SMOKE_LAUNCHER "")
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/benchmark_smoke.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/benchmark_smoke.py @ONLY)
add_test(NAME benchmark_smoke
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/benchmark_smoke.py)
```

- [ ] **Step 2: スモークテストが通ることを確認**

Run:
```bash
(cd build && cmake ../ > /dev/null && ctest -R benchmark_smoke --output-on-failure)
```
Expected: `1/1 Test #N: benchmark_smoke .... Passed`(数十秒以内)

- [ ] **Step 3: README を書く**

`benchmark/README.md`:

```markdown
# TeNeS benchmark harness

A/B performance comparison harness for TeNeS. It runs the
`tenes_simple -> tenes_std -> tenes` pipeline over parameterized cases,
collects per-timer results from `output/timers.json`, and renders a
Markdown comparison report.

## Quick start (A/B comparison)

```bash
# 1. Build variant A (e.g. develop) and run the suite
python3 benchmark/bench.py run benchmark/suites/contraction.toml \
    --label baseline --tenes-dir build/src --tool-dir build/tool

# 2. Switch to variant B (e.g. cotengra branch), rebuild, and run again
python3 benchmark/bench.py run benchmark/suites/contraction.toml \
    --label cotengra --tenes-dir build/src --tool-dir build/tool

# 3. Compare
python3 benchmark/bench.py compare \
    benchmark/results/baseline benchmark/results/cotengra \
    --output report.md
```

Alternatively keep two build directories and switch `--tenes-dir`.

## Timers

`tenes` always writes `output/timers.json`: cumulative wall times keyed by
hierarchical names.

- `total` — whole process (as `time.dat`'s "time all")
- `phase/*` — phase times (same values as `time.dat`); `count` is
  meaningless (always 1) for these entries
- `contract/<backend>/<N>x<M>` — contraction kernels per unit-cell shape;
  backends: `itps_ctm`, `itps_mf`, `density_ctm`
- `max_rank` / `min_rank` — per-rank extrema under MPI (load imbalance)

To add a new instrumentation point, place
`ScopedTimer scoped_timer("your/name");` (declared in `src/timer.hpp`)
inside the scope to measure. The harness and report pick it up
automatically; names are grouped by their first path component.

## Suites

Suite files (TOML) live in `benchmark/suites/`:

- `contraction.toml` — unit-cell x D x chi sweep (main suite)
- `e2e.toml` — representative end-to-end workloads
- `smoke.toml` — minimal case for CI smoke test
- `ci.toml` — reduced suite for GitHub Actions

A `[[case]]` either renders a `template` (a `simple.toml` with `${var}`
placeholders, parameters from `params` plus the direct product of `sweep`)
or copies an existing `input` file (a `simple.toml`). Derived parameters:
`Lsub = [L, W]` expands to `L`/`W`; `chi_ratio = r` yields
`chi = r * D * D`. Case names must be unique after expansion.

## Statistics and caveats

- Each case runs `repeat` times; the report shows median [min, max] and
  the ratio of medians (B/A). If the min-max intervals overlap, the row
  is marked "no sig. diff".
- If call counts differ between A and B ("count differs"), the algorithm
  changed and time ratios should be interpreted accordingly.
- Observables (`output/*.dat` except `time.dat` / `parameters.dat`) of
  `run_0` are compared with `np.isclose` (rtol 1e-3, atol 1e-4); any
  mismatch is reported as a WARNING at the top: faster-but-wrong must
  not go unnoticed.

## HPC (MPI) usage

Pass the launcher prefix and run both labels inside the *same* job to
avoid node-placement noise:

```bash
python3 benchmark/bench.py run benchmark/suites/contraction.toml \
    --label baseline --launcher "srun -n 64" ...
python3 benchmark/bench.py run benchmark/suites/contraction.toml \
    --label cotengra --launcher "srun -n 64" ...
```

`OMP_NUM_THREADS` and the launcher string are recorded in `meta.json`.

## GitHub Actions

The `Benchmark` workflow (manual trigger: workflow_dispatch) runs
`suites/ci.toml` on a shared runner and uploads the results directory as
an artifact. Shared-runner numbers are noisy; treat them as reference
values only and compare locally.
```

- [ ] **Step 4: コミット**

```bash
python3 -m pytest test/python -q
git add test/benchmark_smoke.py.in test/CMakeLists.txt benchmark/README.md
git commit -m "Add benchmark smoke test to ctest and harness README"
```

---

### Task 11: GitHub Actions ワークフローと総仕上げ

**Files:**
- Create: `.github/workflows/benchmark.yml`
- Modify: なし

**Interfaces:**
- Consumes: Task 9 の `bench.py` / `suites/ci.toml`。actions のバージョンとハッシュは `.github/workflows/linux.yml` の pin に合わせる(実装時に linux.yml から checkout / setup-python / upload-artifact の `uses:` 行をコピーする)。
- Produces: 手動トリガーの `Benchmark` ワークフロー(結果 JSON を artifact 化)

- [ ] **Step 1: ワークフローを書く**

`.github/workflows/benchmark.yml`(`uses:` のハッシュ pin は linux.yml の該当行をコピーして合わせること):

```yaml
name: Benchmark

on:
  workflow_dispatch:
    inputs:
      suite:
        description: "suite file under benchmark/suites/"
        default: "ci.toml"
      label:
        description: "label for this run"
        default: "ci"

jobs:
  benchmark:
    runs-on: ubuntu-24.04
    timeout-minutes: 60
    steps:
      - uses: actions/checkout@de0fac2e4500dabe0009e67214ff5f5447ce83dd # v6.0.2

      - name: Setup Python
        uses: actions/setup-python@a309ff8b426b58ec0e2a45f0f869d46889d02405 # v6.2.0
        with:
          python-version: "3.13"

      - name: apt
        run: |
          sudo apt update
          sudo apt install libblas-dev liblapack-dev

      - name: pip install
        run: |
          python -m pip install -U pip
          python -m pip install numpy scipy toml

      - name: cmake
        run: |
          cmake -E make_directory build
          cmake -B build -DENABLE_MPI=OFF -DTesting=OFF .

      - name: build
        run: cmake --build build -j4

      - name: run benchmark
        run: |
          python benchmark/bench.py run \
            "benchmark/suites/${{ github.event.inputs.suite }}" \
            --label "${{ github.event.inputs.label }}" \
            --tenes-dir build/src --tool-dir build/tool

      - name: upload results
        uses: actions/upload-artifact@v4
        with:
          name: benchmark-results
          path: benchmark/results/
```

注: `upload-artifact` は linux.yml に無い場合があるので、その場合は GitHub 公式の最新 v4 ハッシュを `gh api repos/actions/upload-artifact/git/refs/tags/v4` 等で確認して pin する(できなければ `@v4` のままにして PR レビューで指摘を仰ぐ)。

- [ ] **Step 2: YAML の構文チェック**

Run: `python3 -c "import yaml, sys; yaml.safe_load(open('.github/workflows/benchmark.yml'))" && echo OK`
(pyyaml が無ければ `python3 -m pip install --user pyyaml` するか、`gh workflow list` 用に push 後確認とする)
Expected: `OK`

- [ ] **Step 3: 全テストで回帰確認**

Run:
```bash
cmake --build build -j
(cd build && ctest --output-on-failure)
```
Expected: 全テスト PASS(既存 23 + timer_registry + benchmark_smoke = 25)。

- [ ] **Step 4: コミット**

```bash
git add .github/workflows/benchmark.yml
git commit -m "Add manually-triggered benchmark workflow for GitHub Actions"
```

---

## 完了条件

- `(cd build && ctest --output-on-failure)` が全件 PASS。
- `python3 benchmark/bench.py run benchmark/suites/smoke.toml --label x` → 同 `--label y` → `compare` が警告なしのレポートを出す。
- `output/timers.json` に `total` / `phase/*` / `contract/*` が記録される(CTM・MF・density の3バックエンドとも)。
