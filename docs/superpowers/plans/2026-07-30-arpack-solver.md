# ARPACK-NG 転送行列固有値ソルバ 実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 相関長計算の転送行列固有値計算を ARPACK-NG で実行できるようにする(オプショナル依存+自前 Arnoldi フォールバック、入力 toml で切替)。

**Architecture:** 全ランクが同一のシリアル ARPACK reverse communication ループを冗長実行し、matvec のみ既存の分散テンソル縮約を使う。新規 `src/arpack_solver.{hpp,cpp}` に ARPACK 呼び出しを隔離し、`TransferMatrix::eigenvalues`(`src/iTPS/transfer_matrix.cpp`)の Arnoldi 分岐に第 3 の選択肢として組み込む。

**Tech Stack:** C++17 / ARPACK-NG(Fortran シンボル直接呼び出し、ICB 不使用)/ mptensor / doctest / CMake

**Spec:** `docs/superpowers/specs/2026-07-30-arpack-ng-design.md`

## Global Constraints

- C++17。すべての C++ ソースは namespace `tenes`(iTPS 層は `tenes::itps`)。
- 新規ソースファイルには GPL v3 ヘッダ(既存ファイル冒頭 15 行と同一)を付ける。
- C++ 整形は `clang-format`(リポジトリの `.clang-format`)。コミット前に新規・変更ファイルへ適用する。
- ビルド: `cmake --preset gcc && cmake --build --preset gcc`、テスト: `ctest --preset gcc`(macOS ローカル。Debug, g++-16, MPI OFF)。ビルドディレクトリはプリセット定義により `out-gcc/build/`(`build/` ではない)。Python は `venv/bin/python3`(numpy/scipy/toml 導入済み)。
- スクラッチ作業は `/private/tmp/claude-501/-Users-yomichi-source-github-com-issp-center-dev-TeNeS/dc8e166b-4e2e-4acd-9258-19d3d1f4f54b/scratchpad`(以下 `$SCRATCH`)で行う。
- コンパイル定義 `TENES_USE_ARPACK` は ARPACK あり時のみ定義される。CMake 変数 `TENES_ARPACK_FOUND` が ON/OFF を保持する。
- コミットメッセージ末尾に以下を付ける:
  ```
  Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01XqkDpnELWYMqSqSJXkGyMu
  ```

---

### Task 1: CMake — ARPACK-NG の検出(ENABLE_ARPACK)

**Files:**
- Modify: `CMakeLists.txt`(トップレベル。`find_package(LAPACK REQUIRED)` の直後、および `TENES_COMPILE_OPTIONS` 組み立てブロック)

**Interfaces:**
- Produces: CMake 変数 `TENES_ARPACK_FOUND`(ON/OFF)、`TENES_ARPACK_LIBRARIES`(リンク対象)、コンパイル定義 `-DTENES_USE_ARPACK`(あり時のみ)。後続タスクの `src/CMakeLists.txt`・`test/CMakeLists.txt` がこの 2 変数を参照する。

- [ ] **Step 1: ローカルに ARPACK-NG があるか確認し、なければ入れる**

Run: `brew list arpack >/dev/null 2>&1 && echo present || brew install arpack`
Expected: `present` または brew によるインストール完了。
(brew の arpack は `$(brew --prefix arpack)` 配下に lib/libarpack.dylib を持つ)

- [ ] **Step 2: トップレベル CMakeLists.txt に検出ブロックを追加**

`find_package(LAPACK REQUIRED)` の直後(`if(ENABLE_MPI)` の前)に挿入:

```cmake
# ARPACK-NG (optional): eigensolver for the transfer matrix.
# AUTO: use if found / ON: required / OFF: do not search
set(ENABLE_ARPACK "AUTO" CACHE STRING
    "Use ARPACK-NG for the transfer-matrix eigensolver (AUTO/ON/OFF)")
set_property(CACHE ENABLE_ARPACK PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${ENABLE_ARPACK}" _enable_arpack)
set(TENES_ARPACK_FOUND OFF)
set(TENES_ARPACK_LIBRARIES "")
if(NOT _enable_arpack STREQUAL "OFF")
    find_package(arpackng QUIET CONFIG NAMES arpackng arpack-ng
                 HINTS ${ARPACK_ROOT})
    if(TARGET ARPACK::ARPACK)
        set(TENES_ARPACK_FOUND ON)
        set(TENES_ARPACK_LIBRARIES ARPACK::ARPACK)
    else()
        find_library(ARPACK_LIBRARY NAMES arpack
                     HINTS ${ARPACK_ROOT} PATH_SUFFIXES lib lib64)
        if(ARPACK_LIBRARY)
            set(TENES_ARPACK_FOUND ON)
            set(TENES_ARPACK_LIBRARIES ${ARPACK_LIBRARY})
        endif()
    endif()
    if(TENES_ARPACK_FOUND)
        message(STATUS "ARPACK-NG: ${TENES_ARPACK_LIBRARIES}")
    elseif(_enable_arpack STREQUAL "ON")
        message(FATAL_ERROR
            "ENABLE_ARPACK=ON but ARPACK-NG was not found "
            "(hint: -DARPACK_ROOT=<prefix>)")
    else()
        message(STATUS "ARPACK-NG: not found (builtin Arnoldi only)")
    endif()
endif()
```

- [ ] **Step 3: コンパイル定義を追加**

`TENES_COMPILE_OPTIONS` を組み立てているブロック(`if(NOT ENABLE_MPI)` で `-D_NO_MPI` を足している箇所)の直後に追加:

```cmake
if(TENES_ARPACK_FOUND)
    set(TENES_COMPILE_OPTIONS -DTENES_USE_ARPACK ${TENES_COMPILE_OPTIONS})
endif()
```

(`TENES_COMPILE_OPTIONS` は `tenes_mpi` に PUBLIC で付与され全ターゲットに伝播するので、
src/ と test/ の両方で `TENES_USE_ARPACK` が一貫して見える)

- [ ] **Step 4: 検出ありの configure を確認**

Run: `cmake --preset gcc 2>&1 | grep -i arpack`
Expected: `ARPACK-NG: <brewのlibarpackパス or ARPACK::ARPACK>` と `TENES_COMPILE_OPTIONS` に `-DTENES_USE_ARPACK`

- [ ] **Step 5: OFF と AUTO 未検出の configure を確認**

Run: `cmake --preset gcc -DENABLE_ARPACK=OFF 2>&1 | grep -i arpack; cmake --preset gcc -DENABLE_ARPACK=AUTO -DARPACK_ROOT=/nonexistent 2>&1 | grep -i arpack`
Expected: 1 回目は ARPACK の status 行なし(または not found 表示なしで正常終了)、2 回目は… AUTO でも brew の arpack がデフォルト検索パスで見つかる場合があるので、見つかっても失敗しないことだけ確認。最後に `cmake --preset gcc` で AUTO(検出あり)に戻す。
Expected: いずれも configure が成功する。

- [ ] **Step 6: ビルドが従来どおり通ることを確認してコミット**

Run: `cmake --build --preset gcc -j 2>&1 | tail -3`
Expected: エラーなし(この時点ではソース変更なしなので挙動不変)

```bash
git add CMakeLists.txt
git commit -m "Add optional ARPACK-NG detection (ENABLE_ARPACK=AUTO/ON/OFF)"
```

---

### Task 2: mpi.hpp — 複素ベクトルの allreduce_sum

**Files:**
- Modify: `src/mpi.hpp`(現在 `unimplemented_error` を投げている `allreduce_sum(std::vector<std::complex<T>>&, MPI_Comm)` を実装に置き換え)
- Test: `test/util.cpp`

**Interfaces:**
- Produces: `int tenes::allreduce_sum(std::vector<std::complex<T>> &val, MPI_Comm comm)` — 実部・虚部を分解して総和。Task 3 の gather ヘルパが使う。
- 注意: スカラー版 `allreduce_sum(std::complex<T>, MPI_Comm)` は使わないので throw のまま残す(YAGNI)。

- [ ] **Step 1: 失敗するテストを書く**

`test/util.cpp` の doctest main の有無を確認(`grep -n DOCTEST_CONFIG_IMPLEMENT test/util.cpp`)。
`#define DOCTEST_CONFIG_IMPLEMENT` + `main` があるファイル(なければ test/arnoldi.cpp:17-34 と同じ main を持つ別テストファイルに追記)の末尾に追加:

```cpp
TEST_CASE("allreduce_sum for complex vectors") {
  std::vector<std::complex<double>> v = {{1.0, 2.0}, {-3.0, 4.5}};
  int mpisize = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  tenes::allreduce_sum(v, MPI_COMM_WORLD);
  const double s = static_cast<double>(mpisize);
  CHECK(v[0].real() == doctest::Approx(1.0 * s));
  CHECK(v[0].imag() == doctest::Approx(2.0 * s));
  CHECK(v[1].real() == doctest::Approx(-3.0 * s));
  CHECK(v[1].imag() == doctest::Approx(4.5 * s));
}
```

`#include <complex>` と `#include "../src/mpi.hpp"` が無ければ足す。

- [ ] **Step 2: テストが失敗することを確認**

Run: `cmake --build --preset gcc -j --target test_util && ctest --preset gcc -R test_util -V | tail -5`
Expected: MPI ビルドでは `unimplemented_error` 送出で FAIL。
(注: `_NO_MPI` ビルド(gcc プリセット)では `#ifndef _NO_MPI` の外側が no-op なので throw に到達せず PASS してしまう。その場合はこのステップの「失敗確認」は形式的に PASS でよい — 実装の意味があるのは MPI ビルドで、コードパスは Step 3 の実装で MPI 有無両対応になる)

- [ ] **Step 3: 実装する**

`src/mpi.hpp` の throw している `allreduce_sum(std::vector<std::complex<T>>...)` を置き換え:

```cpp
template <class T>
int allreduce_sum(std::vector<std::complex<T>> &val, MPI_Comm comm) {
  int ret = 0;
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  const int sz = static_cast<int>(val.size());
  std::vector<T> reim(2 * sz);
  for (int i = 0; i < sz; ++i) {
    reim[2 * i] = val[i].real();
    reim[2 * i + 1] = val[i].imag();
  }
  std::vector<T> recv(reim);
  ret = MPI_Allreduce(reim.data(), recv.data(), 2 * sz, datatype, MPI_SUM,
                      comm);
  if (ret != 0) {
    return ret;
  }
  for (int i = 0; i < sz; ++i) {
    val[i] = std::complex<T>(recv[2 * i], recv[2 * i + 1]);
  }
#endif
  return ret;
}
```

- [ ] **Step 4: テストが通ることを確認**

Run: `cmake --build --preset gcc -j --target test_util && ctest --preset gcc -R test_util`
Expected: PASS

- [ ] **Step 5: コミット**

```bash
git add src/mpi.hpp test/util.cpp
git commit -m "Implement allreduce_sum for complex vectors"
```

---

### Task 3: arpack_solver — ARPACK 呼び出しの本体

**Files:**
- Create: `src/arpack_solver.hpp`
- Create: `src/arpack_solver.cpp`
- Create: `test/arpack.cpp`
- Modify: `src/CMakeLists.txt`(`tensor` ライブラリにソース追加とリンク)
- Modify: `test/CMakeLists.txt`(ARPACK あり時のみ test_arpack を登録)

**Interfaces:**
- Consumes: `tenes::allreduce_sum(std::vector<std::complex<T>>&, MPI_Comm)`(Task 2)
- Produces:
  - `constexpr bool tenes::arpack_available()`(ヘッダ、ビルド時定数)
  - `template <class ptensor> std::vector<std::complex<double>> tenes::arpack_eigenvalues(std::function<void(ptensor&, ptensor const&)> A, ptensor const& initial, std::size_t nev, int ncv, int maxiter, double tol)` — |λ| 降順の固有値を size nev で返す。`real_tensor`/`complex_tensor` で明示的実体化。ARPACK なしビルドでは呼ぶと `std::logic_error`。

- [ ] **Step 1: ヘッダを書く**

`src/arpack_solver.hpp`(GPL ヘッダ 15 行を先頭に付ける):

```cpp
#ifndef TENES_SRC_ARPACK_SOLVER_HPP_
#define TENES_SRC_ARPACK_SOLVER_HPP_

#include <complex>
#include <cstddef>
#include <functional>
#include <vector>      // IWYU pragma: export
#include "tensor.hpp"  // IWYU pragma: export

namespace tenes {

//! @return true if TeNeS was built with ARPACK-NG
constexpr bool arpack_available() {
#ifdef TENES_USE_ARPACK
  return true;
#else
  return false;
#endif
}

/*! @brief Largest-magnitude eigenvalues via ARPACK-NG (dnaupd/znaupd)
 *
 * Every rank redundantly runs the same serial ARPACK iteration; only the
 * matrix-vector product runs on the distributed tensor.
 *
 * @param[in] A "Matrix" as a function taking a "vector" `x` and returning
 *              another, `Ax`. The first [out] argument of `A` is for `Ax`
 *              and the second [in] argument is for `x` (same as Arnoldi).
 * @param[in] initial Initial (distributed) vector; also fixes the problem
 *                    size N via its shape.
 * @param[in] nev Number of eigenvalues to be calculated
 * @param[in] ncv Dimension of the Krylov subspace (clamped into
 *                [nev+2, N] internally)
 * @param[in] maxiter Maximum number of implicit restarts
 * @param[in] tol Relative tolerance on the Ritz values (<= 0 means
 *                machine epsilon)
 * @return Eigenvalues sorted by decreasing magnitude, size nev. If fewer
 *         than nev eigenvalues converge, the tail is NaN and a warning is
 *         printed on rank 0.
 *
 * Throws tenes::runtime_error on ARPACK failures and std::logic_error if
 * TeNeS was built without ARPACK-NG.
 */
template <class ptensor>
std::vector<std::complex<double>> arpack_eigenvalues(
    std::function<void(ptensor &, ptensor const &)> A, ptensor const &initial,
    std::size_t nev, int ncv, int maxiter, double tol);

}  // end of namespace tenes

#endif  // TENES_SRC_ARPACK_SOLVER_HPP_
```

- [ ] **Step 2: 失敗するテストを書く**

`test/arpack.cpp`(GPL ヘッダ付き。`test/arnoldi.cpp` と同じ構造):

```cpp
#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <cmath>
#include <complex>
#include <cstddef>

#include "../src/arpack_solver.hpp"
#include "../src/arnoldi.hpp"
#include "../src/tensor.hpp"
#include "../src/mpi.hpp"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  doctest::Context context(argc, argv);
  const int res = context.run();
  MPI_Finalize();
  return res;
}

namespace {

using rtensor = tenes::real_tensor;
using ctensor = tenes::complex_tensor;

// A x for the diagonal matrix A = diag(2^0, 2^-1, ..., 2^-(N-1)), scaled by 16
void matvec_real(rtensor &out, rtensor const &in, std::size_t N) {
  out = rtensor(mptensor::Shape(N));
  for (std::size_t i = 0; i < N; ++i) {
    double v;
    in.get_value({i}, v);
    out.set_value({i}, 16.0 * std::pow(2.0, -static_cast<double>(i)) * v);
  }
}

// the same spectrum rotated to the imaginary axis: A = diag(16i * 2^-i)
void matvec_complex(ctensor &out, ctensor const &in, std::size_t N) {
  out = ctensor(mptensor::Shape(N));
  const std::complex<double> I(0.0, 1.0);
  for (std::size_t i = 0; i < N; ++i) {
    std::complex<double> v;
    in.get_value({i}, v);
    out.set_value({i},
                  16.0 * I * std::pow(2.0, -static_cast<double>(i)) * v);
  }
}

rtensor ones_real(std::size_t N) {
  rtensor v{mptensor::Shape(N)};
  for (std::size_t i = 0; i < N; ++i) {
    v.set_value({i}, 1.0);
  }
  return v;
}

ctensor ones_complex(std::size_t N) {
  ctensor v{mptensor::Shape(N)};
  for (std::size_t i = 0; i < N; ++i) {
    v.set_value({i}, std::complex<double>(1.0, 0.0));
  }
  return v;
}

}  // namespace

TEST_CASE("arpack_eigenvalues, real") {
  const std::size_t N = 20;
  const std::size_t nev = 2;
  auto A = [](rtensor &out, rtensor const &in) { matvec_real(out, in, N); };

  auto ev = tenes::arpack_eigenvalues<rtensor>(A, ones_real(N), nev, 10, 10,
                                               1.0e-10);
  REQUIRE(ev.size() == nev);
  CHECK(std::abs(ev[0]) == doctest::Approx(16.0).epsilon(1.0e-8));
  CHECK(std::abs(ev[1]) == doctest::Approx(8.0).epsilon(1.0e-8));
  CHECK(ev[0].imag() == doctest::Approx(0.0).epsilon(1.0e-8));
}

TEST_CASE("arpack_eigenvalues, complex") {
  const std::size_t N = 20;
  const std::size_t nev = 2;
  auto A = [](ctensor &out, ctensor const &in) { matvec_complex(out, in, N); };

  auto ev = tenes::arpack_eigenvalues<ctensor>(A, ones_complex(N), nev, 10, 10,
                                               1.0e-10);
  REQUIRE(ev.size() == nev);
  // eigenvalues are purely imaginary: 16i, 8i, ...
  CHECK(ev[0].real() == doctest::Approx(0.0).epsilon(1.0e-8));
  CHECK(ev[0].imag() == doctest::Approx(16.0).epsilon(1.0e-8));
  CHECK(std::abs(ev[1]) == doctest::Approx(8.0).epsilon(1.0e-8));
}

TEST_CASE("arpack agrees with the builtin Arnoldi") {
  const std::size_t N = 20;
  const std::size_t nev = 2;
  auto A = [](rtensor &out, rtensor const &in) { matvec_real(out, in, N); };

  auto ev_arpack = tenes::arpack_eigenvalues<rtensor>(A, ones_real(N), nev, 10,
                                                      10, 1.0e-10);

  tenes::Arnoldi<rtensor> arnoldi(N, 10);
  arnoldi.initialize(ones_real(N));
  arnoldi.run(A, nev, 5, 10, 1.0e-8);
  auto ev_builtin = arnoldi.eigenvalues();

  for (std::size_t i = 0; i < nev; ++i) {
    CHECK(std::abs(ev_arpack[i]) ==
          doctest::Approx(std::abs(ev_builtin[i])).epsilon(1.0e-6));
  }
}
```

- [ ] **Step 3: CMake にソースとテストを登録**

`src/CMakeLists.txt`:

```cmake
add_library(tensor STATIC tensor.cpp arnoldi.cpp arpack_solver.cpp)
target_include_directories(tensor PUBLIC ${MPTENSOR_INCLUDE_DIR})
target_link_libraries(tensor INTERFACE ${SCALAPACK_LIBRARIES}
                                       ${LAPACK_LIBRARIES}
                                       ${TENES_ARPACK_LIBRARIES} mptensor)
```

`test/CMakeLists.txt` の foreach を変数化して条件追加:

```cmake
set(tenes_unit_tests input simple_update full_update tensor_util util
    correlation_length saveload arnoldi timer_registry)
if(TENES_ARPACK_FOUND)
  list(APPEND tenes_unit_tests arpack)
endif()
foreach(basename IN LISTS tenes_unit_tests)
```

(foreach の中身は変更しない)

- [ ] **Step 4: ビルドが失敗することを確認**

Run: `cmake --preset gcc && cmake --build --preset gcc -j 2>&1 | tail -5`
Expected: `arpack_solver.cpp` が無いので FAIL(または作成済みヘッダのみで cpp 欠如のエラー)

- [ ] **Step 5: 実装を書く**

`src/arpack_solver.cpp`(GPL ヘッダ付き)。全文:

```cpp
#include "arpack_solver.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "exception.hpp"
#include "mpi.hpp"
#include "util/abs.hpp"

namespace tenes {

#ifdef TENES_USE_ARPACK

// Direct declarations of the four ARPACK-NG routines (arpack++ style).
// The trailing std::size_t arguments are the hidden lengths of the
// CHARACTER dummy arguments appended by Fortran compilers (gfortran etc.).
extern "C" {
void dnaupd_(int *ido, const char *bmat, const int *n, const char *which,
             const int *nev, const double *tol, double *resid, const int *ncv,
             double *v, const int *ldv, int *iparam, int *ipntr, double *workd,
             double *workl, const int *lworkl, int *info, std::size_t bmat_len,
             std::size_t which_len);
void dneupd_(const int *rvec, const char *howmny, const int *select,
             double *dr, double *di, double *z, const int *ldz,
             const double *sigmar, const double *sigmai, double *workev,
             const char *bmat, const int *n, const char *which, const int *nev,
             const double *tol, double *resid, const int *ncv, double *v,
             const int *ldv, int *iparam, int *ipntr, double *workd,
             double *workl, const int *lworkl, int *info,
             std::size_t howmny_len, std::size_t bmat_len,
             std::size_t which_len);
void znaupd_(int *ido, const char *bmat, const int *n, const char *which,
             const int *nev, const double *tol, std::complex<double> *resid,
             const int *ncv, std::complex<double> *v, const int *ldv,
             int *iparam, int *ipntr, std::complex<double> *workd,
             std::complex<double> *workl, const int *lworkl, double *rwork,
             int *info, std::size_t bmat_len, std::size_t which_len);
void zneupd_(const int *rvec, const char *howmny, const int *select,
             std::complex<double> *d, std::complex<double> *z, const int *ldz,
             const std::complex<double> *sigma, std::complex<double> *workev,
             const char *bmat, const int *n, const char *which, const int *nev,
             const double *tol, std::complex<double> *resid, const int *ncv,
             std::complex<double> *v, const int *ldv, int *iparam, int *ipntr,
             std::complex<double> *workd, std::complex<double> *workl,
             const int *lworkl, double *rwork, int *info,
             std::size_t howmny_len, std::size_t bmat_len,
             std::size_t which_len);
}

namespace {

using dcomplex = std::complex<double>;

[[noreturn]] void throw_arpack_error(const char *routine, int info) {
  std::stringstream ss;
  ss << "ARPACK routine " << routine << " failed with info = " << info
     << " (see the ARPACK-NG documentation of " << routine << ")";
  if (info == 3 || info == -3) {
    ss << "; hint: increase correlation_length.arnoldi_maxdim";
  }
  throw tenes::runtime_error(ss.str());
}

// gather a distributed rank-1 tensor into a full local copy on every rank
template <class ptensor>
std::vector<typename ptensor::value_type> gather_vector(ptensor const &t) {
  using value_type = typename ptensor::value_type;
  const std::size_t N = t.shape()[0];
  std::vector<value_type> buf(N, 0.0);
  for (std::size_t i = 0; i < t.local_size(); ++i) {
    const auto index = t.global_index(i);
    value_type v;
    if (t.get_value(index, v)) {
      buf[index[0]] = v;
    }
  }
  // every element has exactly one owner, so the sum only adds zeros and
  // the result is bitwise identical on all ranks
  allreduce_sum(buf, t.get_comm());
  return buf;
}

// scatter a full local copy into the distributed rank-1 tensor
template <class ptensor>
void scatter_vector(std::vector<typename ptensor::value_type> const &buf,
                    ptensor &t) {
  for (std::size_t i = 0; i < t.local_size(); ++i) {
    const auto index = t.global_index(i);
    t.set_value(index, buf[index[0]]);
  }
}

using serial_matvec_real =
    std::function<void(std::vector<double> &, std::vector<double> const &)>;
using serial_matvec_complex =
    std::function<void(std::vector<dcomplex> &, std::vector<dcomplex> const &)>;

std::vector<dcomplex> run_arpack(serial_matvec_real const &av,
                                 std::vector<double> &resid, int nev, int ncv,
                                 int maxiter, double tol, bool print_warn) {
  const int n = static_cast<int>(resid.size());
  int ido = 0;
  int info = 1;  // resid contains the initial vector
  int iparam[11] = {};
  int ipntr[14] = {};
  iparam[0] = 1;        // exact shifts
  iparam[2] = maxiter;  // max implicit restarts
  iparam[6] = 1;        // mode 1: standard eigenvalue problem
  const int lworkl = 3 * ncv * ncv + 6 * ncv;
  std::vector<double> v(static_cast<std::size_t>(n) * ncv), workd(3 * n),
      workl(lworkl);
  std::vector<double> x(n), y(n);

  while (true) {
    dnaupd_(&ido, "I", &n, "LM", &nev, &tol, resid.data(), &ncv, v.data(), &n,
            iparam, ipntr, workd.data(), workl.data(), &lworkl, &info, 1, 2);
    if (ido != 1 && ido != -1) {
      break;
    }
    std::copy_n(workd.data() + ipntr[0] - 1, n, x.begin());
    av(y, x);
    std::copy_n(y.data(), n, workd.data() + ipntr[1] - 1);
  }
  if (info < 0 || info == 3) {
    throw_arpack_error("dnaupd", info);
  }
  if (info == 1 && print_warn) {
    std::cerr << "WARNING: ARPACK (dnaupd) reached the maximum number of "
                 "restarts ("
              << maxiter << ") before all eigenvalues converged" << std::endl;
  }

  const int rvec = 0;  // eigenvalues only
  std::vector<int> select(ncv);
  std::vector<double> dr(nev + 1), di(nev + 1), workev(3 * ncv);
  const double sigmar = 0.0, sigmai = 0.0;
  int ierr = 0;
  dneupd_(&rvec, "A", select.data(), dr.data(), di.data(), v.data(), &n,
          &sigmar, &sigmai, workev.data(), "I", &n, "LM", &nev, &tol,
          resid.data(), &ncv, v.data(), &n, iparam, ipntr, workd.data(),
          workl.data(), &lworkl, &ierr, 1, 1, 2);
  if (ierr != 0) {
    throw_arpack_error("dneupd", ierr);
  }

  const int nconv = std::min(iparam[4], nev);
  std::vector<dcomplex> ev;
  ev.reserve(nconv);
  for (int i = 0; i < nconv; ++i) {
    ev.emplace_back(dr[i], di[i]);
  }
  return ev;
}

std::vector<dcomplex> run_arpack(serial_matvec_complex const &av,
                                 std::vector<dcomplex> &resid, int nev,
                                 int ncv, int maxiter, double tol,
                                 bool print_warn) {
  const int n = static_cast<int>(resid.size());
  int ido = 0;
  int info = 1;  // resid contains the initial vector
  int iparam[11] = {};
  int ipntr[14] = {};
  iparam[0] = 1;        // exact shifts
  iparam[2] = maxiter;  // max implicit restarts
  iparam[6] = 1;        // mode 1: standard eigenvalue problem
  const int lworkl = 3 * ncv * ncv + 5 * ncv;
  std::vector<dcomplex> v(static_cast<std::size_t>(n) * ncv), workd(3 * n),
      workl(lworkl);
  std::vector<double> rwork(ncv);
  std::vector<dcomplex> x(n), y(n);

  while (true) {
    znaupd_(&ido, "I", &n, "LM", &nev, &tol, resid.data(), &ncv, v.data(), &n,
            iparam, ipntr, workd.data(), workl.data(), &lworkl, rwork.data(),
            &info, 1, 2);
    if (ido != 1 && ido != -1) {
      break;
    }
    std::copy_n(workd.data() + ipntr[0] - 1, n, x.begin());
    av(y, x);
    std::copy_n(y.data(), n, workd.data() + ipntr[1] - 1);
  }
  if (info < 0 || info == 3) {
    throw_arpack_error("znaupd", info);
  }
  if (info == 1 && print_warn) {
    std::cerr << "WARNING: ARPACK (znaupd) reached the maximum number of "
                 "restarts ("
              << maxiter << ") before all eigenvalues converged" << std::endl;
  }

  const int rvec = 0;  // eigenvalues only
  std::vector<int> select(ncv);
  std::vector<dcomplex> d(nev + 1), workev(2 * ncv);
  const dcomplex sigma(0.0, 0.0);
  int ierr = 0;
  zneupd_(&rvec, "A", select.data(), d.data(), v.data(), &n, &sigma,
          workev.data(), "I", &n, "LM", &nev, &tol, resid.data(), &ncv,
          v.data(), &n, iparam, ipntr, workd.data(), workl.data(), &lworkl,
          rwork.data(), &ierr, 1, 1, 2);
  if (ierr != 0) {
    throw_arpack_error("zneupd", ierr);
  }

  const int nconv = std::min(iparam[4], nev);
  return std::vector<dcomplex>(d.begin(), d.begin() + nconv);
}

}  // namespace

template <class ptensor>
std::vector<std::complex<double>> arpack_eigenvalues(
    std::function<void(ptensor &, ptensor const &)> A, ptensor const &initial,
    std::size_t nev, int ncv, int maxiter, double tol) {
  using value_type = typename ptensor::value_type;
  const int N = static_cast<int>(initial.shape()[0]);
  const int nev_ = static_cast<int>(nev);
  const bool print_warn = initial.get_comm_rank() == 0;

  int ncv_ = std::min(ncv, N);
  if (ncv_ < nev_ + 2) {
    ncv_ = std::min(nev_ + 2, N);
  }
  if (nev_ >= ncv_ - 1) {
    std::stringstream ss;
    ss << "ARPACK requires num_eigvals + 2 <= min(arnoldi_maxdim, N) but "
          "num_eigvals = "
       << nev_ << ", N = " << N
       << "; hint: increase correlation_length.maxdim_dense_eigensolver "
          "above N to use the dense eigensolver instead";
    throw tenes::input_error(ss.str());
  }

  const MPI_Comm comm = initial.get_comm();
  std::vector<value_type> resid = gather_vector(initial);
  auto serial_matvec = [&](std::vector<value_type> &out,
                           std::vector<value_type> const &in) {
    ptensor x(comm, mptensor::Shape(N));
    scatter_vector(in, x);
    ptensor y;
    A(y, x);
    out = gather_vector(y);
  };

  std::vector<dcomplex> ev =
      run_arpack(serial_matvec, resid, nev_, ncv_, maxiter, tol, print_warn);

  std::sort(ev.begin(), ev.end(), [](dcomplex const &a, dcomplex const &b) {
    return util::abs2(a) > util::abs2(b);
  });
  if (ev.size() < nev) {
    if (print_warn) {
      std::cerr << "WARNING: ARPACK converged only " << ev.size() << " of "
                << nev << " eigenvalues; the rest are reported as NaN"
                << std::endl;
    }
    ev.resize(nev, dcomplex(std::numeric_limits<double>::quiet_NaN(),
                            std::numeric_limits<double>::quiet_NaN()));
  }
  return ev;
}

#else  // TENES_USE_ARPACK

template <class ptensor>
std::vector<std::complex<double>> arpack_eigenvalues(
    std::function<void(ptensor &, ptensor const &)> /* A */,
    ptensor const & /* initial */, std::size_t /* nev */, int /* ncv */,
    int /* maxiter */, double /* tol */) {
  throw std::logic_error(
      "internal error: arpack_eigenvalues is called, but TeNeS was built "
      "without ARPACK-NG");
}

#endif  // TENES_USE_ARPACK

template std::vector<std::complex<double>> arpack_eigenvalues<real_tensor>(
    std::function<void(real_tensor &, real_tensor const &)>,
    real_tensor const &, std::size_t, int, int, double);
template std::vector<std::complex<double>> arpack_eigenvalues<complex_tensor>(
    std::function<void(complex_tensor &, complex_tensor const &)>,
    complex_tensor const &, std::size_t, int, int, double);

}  // end of namespace tenes
```

実装上の注意:
- `tenes::runtime_error` / `tenes::input_error` は `src/exception.hpp:38,44` に定義済み(確認済み)。
- `util::abs2` は `src/util/abs.hpp`(arnoldi.cpp が同じものを使用)。
- mptensor のローカル要素アクセスは `local_size()` / `global_index(i)` / `get_value` / `set_value`(`src/iTPS/tensors.cpp:114-118` と同じパターン)。
- `ptensor(comm, mptensor::Shape(N))` コンストラクタは `src/iTPS/transfer_matrix.cpp:153` で使用実績あり(`_NO_MPI` 両対応)。

- [ ] **Step 6: ビルドしてテストが通ることを確認**

Run: `cmake --build --preset gcc -j && ctest --preset gcc -R test_arpack -V | tail -10`
Expected: 3 つの TEST_CASE がすべて PASS

- [ ] **Step 7: ARPACK なしビルドでもコンパイルが通ることを確認**

Run:
```bash
SCRATCH=/private/tmp/claude-501/-Users-yomichi-source-github-com-issp-center-dev-TeNeS/dc8e166b-4e2e-4acd-9258-19d3d1f4f54b/scratchpad
cmake -S . -B $SCRATCH/tenes-noarpack -DENABLE_ARPACK=OFF -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_COMPILER=g++-16 -DTENES_PYTHON_EXECUTABLE=$PWD/venv/bin/python3 -DTesting=ON
SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk \
  cmake --build $SCRATCH/tenes-noarpack -j --target tensor 2>&1 | tail -3
```
Expected: `tensor` ターゲットのビルド成功(スタブ側がコンパイルされる)

- [ ] **Step 8: コミット**

```bash
git add src/arpack_solver.hpp src/arpack_solver.cpp test/arpack.cpp src/CMakeLists.txt test/CMakeLists.txt
git commit -m "Add ARPACK-NG eigensolver for distributed matvec"
```

---

### Task 4: 入力パラメータ correlation_length.eigensolver

**Files:**
- Modify: `src/iTPS/transfer_matrix.hpp`(enum とフィールド追加)
- Modify: `src/iTPS/transfer_matrix.cpp`(`Bcast` に追加)
- Modify: `src/iTPS/load_toml.cpp`(`gen_transfer_matrix_parameter`)
- Test: `test/input.cpp`

**Interfaces:**
- Consumes: `tenes::arpack_available()`(Task 3)
- Produces:
  - `enum class tenes::itps::TransferMatrixEigensolver : int { automatic = 0, arpack = 1, builtin = 2 };`(`transfer_matrix.hpp`)
  - `TransferMatrix_Parameters::eigensolver`(既定 `automatic`)。Task 5 が分岐に使う。

- [ ] **Step 1: 失敗するテストを書く**

`test/input.cpp` に SUBCASE を追加(`TEST_CASE("input")` 内、既存 SUBCASE の並び):

```cpp
SUBCASE("correlation_length eigensolver") {
  INFO("correlation_length eigensolver");
  auto toml_default = parse_str(R"([correlation_length])");
  auto p_default = gen_transfer_matrix_parameter(
      toml_default.at("correlation_length"), "correlation_length");
  CHECK(p_default.eigensolver == TransferMatrixEigensolver::automatic);

  auto toml_builtin = parse_str(R"(
[correlation_length]
eigensolver = "builtin"
)");
  auto p_builtin = gen_transfer_matrix_parameter(
      toml_builtin.at("correlation_length"), "correlation_length");
  CHECK(p_builtin.eigensolver == TransferMatrixEigensolver::builtin);

  auto toml_arpack = parse_str(R"(
[correlation_length]
eigensolver = "arpack"
)");
  if (tenes::arpack_available()) {
    auto p_arpack = gen_transfer_matrix_parameter(
        toml_arpack.at("correlation_length"), "correlation_length");
    CHECK(p_arpack.eigensolver == TransferMatrixEigensolver::arpack);
  } else {
    CHECK_THROWS_AS(
        gen_transfer_matrix_parameter(toml_arpack.at("correlation_length"),
                                      "correlation_length"),
        tenes::input_error);
  }

  auto toml_bad = parse_str(R"(
[correlation_length]
eigensolver = "lapack"
)");
  CHECK_THROWS_AS(
      gen_transfer_matrix_parameter(toml_bad.at("correlation_length"),
                                    "correlation_length"),
      tenes::input_error);
}
```

`test/input.cpp` の冒頭に `#include "../src/arpack_solver.hpp"` と
`#include "../src/iTPS/transfer_matrix.hpp"` が無ければ足す。
`gen_transfer_matrix_parameter` が `load_toml.hpp` から公開されていることを確認
(`src/iTPS/load_toml.hpp:46` に宣言あり)。

- [ ] **Step 2: テストが失敗することを確認**

Run: `cmake --build --preset gcc -j --target test_input 2>&1 | tail -5`
Expected: `TransferMatrixEigensolver` が未定義でコンパイルエラー

- [ ] **Step 3: enum・フィールド・Bcast・パースを実装**

`src/iTPS/transfer_matrix.hpp` の `TransferMatrix_Parameters` の直前に:

```cpp
//! How to solve the transfer-matrix eigenproblem (when larger than
//! maxdim_dense_eigensolver)
enum class TransferMatrixEigensolver : int {
  automatic = 0,  //!< ARPACK-NG if built in, otherwise builtin
  arpack = 1,     //!< ARPACK-NG (input error if not built in)
  builtin = 2,    //!< builtin implicit-restart Arnoldi
};
```

`TransferMatrix_Parameters` にフィールドとコンストラクタ初期化を追加:

```cpp
  TransferMatrixEigensolver eigensolver;
```

```cpp
        arnoldi_rtol(1.0e-10),
        eigensolver(TransferMatrixEigensolver::automatic) {}
```

`src/iTPS/transfer_matrix.cpp` の `Bcast` の `PARAMS_INT_INDEX` に
`I_eigensolver` を追加し、SAVE/LOAD 両ブロックに
`SAVE_PARAM(eigensolver, int);` / `LOAD_PARAM(eigensolver, int);` を追加
(`static_cast` は enum class ↔ int を扱えるのでマクロはそのまま使える)。

`src/iTPS/load_toml.cpp` の `gen_transfer_matrix_parameter` の
`load_if(clength.arnoldi_rtol, ...)` の後に:

```cpp
  std::string eigensolver_name = "auto";
  load_if(eigensolver_name, toml, "eigensolver");
  if (eigensolver_name == "auto") {
    clength.eigensolver = TransferMatrixEigensolver::automatic;
  } else if (eigensolver_name == "arpack") {
    if (!arpack_available()) {
      throw tenes::input_error(
          "correlation_length.eigensolver = \"arpack\", but this TeNeS "
          "binary is built without ARPACK-NG "
          "(configure with -DENABLE_ARPACK=ON)");
    }
    clength.eigensolver = TransferMatrixEigensolver::arpack;
  } else if (eigensolver_name == "builtin") {
    clength.eigensolver = TransferMatrixEigensolver::builtin;
  } else {
    throw tenes::input_error(
        "correlation_length.eigensolver must be \"auto\", \"arpack\", or "
        "\"builtin\", but is \"" +
        eigensolver_name + "\"");
  }
```

`load_toml.cpp` の冒頭に `#include "../arpack_solver.hpp"` を追加。
`load_if` が `std::string` を扱えるか確認(`grep -n "load_if" src/iTPS/load_toml.cpp | head`
で既存の string 用例を探す。無ければ toml11 の `toml::find_or` 相当の既存流儀に合わせる)。

- [ ] **Step 4: テストが通ることを確認**

Run: `cmake --build --preset gcc -j --target test_input && ctest --preset gcc -R test_input`
Expected: PASS

- [ ] **Step 5: コミット**

```bash
git add src/iTPS/transfer_matrix.hpp src/iTPS/transfer_matrix.cpp src/iTPS/load_toml.cpp test/input.cpp
git commit -m "Add correlation_length.eigensolver input parameter"
```

---

### Task 5: TransferMatrix::eigenvalues への統合とタイマー

**Files:**
- Modify: `src/iTPS/transfer_matrix.cpp:114-138`(Arnoldi 分岐)
- Modify: `src/iTPS/correlation_length.cpp`(タイマー追加)

**Interfaces:**
- Consumes: `arpack_eigenvalues<ptensor>`(Task 3)、`TransferMatrixEigensolver`(Task 4)、`ScopedTimer`(`src/timer.hpp`)
- Produces: 実行時のソルバ選択。タイマーキー `"measure/correlation_length"`(Task 7 のベンチが参照)。

- [ ] **Step 1: eigenvalues の分岐を書き換える**

`src/iTPS/transfer_matrix.cpp` の `} else {  // use Arnoldi` ブロックを:

```cpp
  } else {  // use an iterative eigensolver (ARPACK-NG or builtin Arnoldi)
    auto maxvec = params.arnoldi_maxdim;
    auto maxiter = params.arnoldi_maxiter;
    if (N < static_cast<size_t>(maxvec)) {
      maxvec = N;
      maxiter = 1;
    }

    ptensor initial_vec = initial_vector(dir, fixed_coord, rng);
    std::function<void(ptensor &, ptensor const &)> matvec;
    if (dir == 0) {
      matvec = [&](ptensor &out, ptensor const &in) {
        matvec_horizontal(out, in, fixed_coord);
      };
    } else {
      matvec = [&](ptensor &out, ptensor const &in) {
        matvec_vertical(out, in, fixed_coord);
      };
    }

    const bool use_arpack =
        params.eigensolver == TransferMatrixEigensolver::arpack ||
        (params.eigensolver == TransferMatrixEigensolver::automatic &&
         arpack_available());
    if (use_arpack) {
      eigvals = arpack_eigenvalues<ptensor>(matvec, initial_vec, nev, maxvec,
                                            maxiter, params.arnoldi_rtol);
    } else {
      Arnoldi<ptensor> arnoldi(N, maxvec);
      arnoldi.initialize(initial_vec);
      arnoldi.run(matvec, nev, params.arnoldi_restartdim, maxiter,
                  params.arnoldi_rtol);
      eigvals = arnoldi.eigenvalues();
    }
  }
```

`transfer_matrix.cpp` の冒頭に `#include "../arpack_solver.hpp"` を追加。

- [ ] **Step 2: タイマーを追加する**

`src/iTPS/correlation_length.cpp` の
`iTPS<ptensor>::measure_transfer_matrix_eigenvalues()` の先頭
(`std::vector<transfer_matrix_eigenvalues_type> res;` の前)に:

```cpp
  ScopedTimer scoped_timer("measure/correlation_length");
```

冒頭に `#include "../timer.hpp"` を追加。

- [ ] **Step 3: ビルドして単体テストと既存ゴールデンテストを回す**

Run: `cmake --build --preset gcc -j && ctest --preset gcc -R 'test_|AntiferroHeisenberg|Kitaev' --output-on-failure | tail -15`
Expected: すべて PASS(AFH 系は χ が小さく密ソルバ経路なので数値は不変。単体テストで ARPACK 経路を担保)

- [ ] **Step 4: 反復経路の実機 A/B 検証(手動)**

ゴールデンテストは χ が小さく密ソルバ経路しか通らないため、
`maxdim_dense_eigensolver = 1` で反復経路を強制して builtin と ARPACK の一致を確認する。
`[correlation_length]` は simple.toml では**トップレベル**セクション
(`test/data/AntiferroHeisenberg_real.toml` には存在しないので追記でよい。
`[correlation]` とは別セクションなことに注意):

```bash
SCRATCH=/private/tmp/claude-501/-Users-yomichi-source-github-com-issp-center-dev-TeNeS/dc8e166b-4e2e-4acd-9258-19d3d1f4f54b/scratchpad
REPO=/Users/yomichi/source/github.com/issp-center-dev/TeNeS
PY=$REPO/venv/bin/python3
TENES=$REPO/out-gcc/build/src/tenes
for solver in builtin arpack; do
  mkdir -p $SCRATCH/ab/$solver && cd $SCRATCH/ab/$solver
  { cat $REPO/test/data/AntiferroHeisenberg_real.toml
    printf '\n[correlation_length]\nmeasure = true\nmaxdim_dense_eigensolver = 1\neigensolver = "%s"\n' $solver
  } > simple.toml
  $PY $REPO/tool/tenes_simple.py simple.toml -o std.toml
  $PY $REPO/tool/tenes_std.py std.toml -o input.toml
  $TENES input.toml
done
paste $SCRATCH/ab/builtin/output/correlation_length.dat \
      $SCRATCH/ab/arpack/output/correlation_length.dat | grep -v '^#' | head
```

Expected: 両者の `correlation_length.dat` の相関長・固有値列が相対誤差 1e-6 程度で一致。

- [ ] **Step 5: timers.json にキーが出ることを確認**

Run: `grep -o 'measure/correlation_length' $SCRATCH/ab/arpack/output/timers.json`
Expected: `measure/correlation_length` が出力される

- [ ] **Step 6: コミット**

```bash
cd /Users/yomichi/source/github.com/issp-center-dev/TeNeS
git add src/iTPS/transfer_matrix.cpp src/iTPS/correlation_length.cpp
git commit -m "Use ARPACK-NG for transfer-matrix eigenvalues when available"
```

---

### Task 6: ドキュメント更新

**Files:**
- Modify: `docs/sphinx/ja/file_specification/correlation_length_section.rst`
- Modify: `docs/sphinx/en/file_specification/correlation_length_section.rst`

**Interfaces:**
- Consumes: Task 4 の入力仕様(`eigensolver` の値と意味)

- [ ] **Step 1: 日本語版を更新**

csv-table に行を追加(`arnoldi_rtol` の行の後):

```
   ``eigensolver``,              "転送行列の固有値ソルバ(auto / arpack / builtin)",   文字列, \"auto\"
```

末尾の説明文の後に追記:

```rst
行列サイズが ``maxdim_dense_eigensolver`` より大きい場合に使う反復固有値ソルバは
``eigensolver`` で選べます。

- ``"auto"`` (デフォルト): ARPACK-NG 付きでビルドされていれば ARPACK-NG を、
  そうでなければ組み込みの IRA 法を使います。
- ``"arpack"``: ARPACK-NG を使います。ARPACK-NG なしでビルドされたバイナリでは
  入力エラーになります(CMake オプション ``-DENABLE_ARPACK=ON`` でビルドしてください)。
- ``"builtin"``: 組み込みの IRA 法を使います。

ARPACK-NG を使う場合、 ``arnoldi_maxdim`` は Krylov 部分空間の次元 (ncv)、
``arnoldi_maxiterations`` は implicit restart の最大回数、 ``arnoldi_rtol`` は
Ritz 値の相対残差の許容値として引き継がれます。
``arnoldi_restartdim`` は組み込みソルバでのみ意味を持ちます。
```

- [ ] **Step 2: 英語版を更新**

英語版 csv-table に行を追加:

```
   ``eigensolver``,              "Eigensolver for the transfer matrix (auto / arpack / builtin)", String, \"auto\"
```

説明文(日本語版と同内容):

```rst
The iterative eigensolver used when the matrix size exceeds
``maxdim_dense_eigensolver`` can be chosen by ``eigensolver``.

- ``"auto"`` (default): use ARPACK-NG if TeNeS is built with it,
  otherwise the builtin IRA solver.
- ``"arpack"``: use ARPACK-NG. This is an input error for binaries built
  without ARPACK-NG (configure with ``-DENABLE_ARPACK=ON``).
- ``"builtin"``: use the builtin IRA solver.

With ARPACK-NG, ``arnoldi_maxdim`` is used as the dimension of the Krylov
subspace (ncv), ``arnoldi_maxiterations`` as the maximum number of implicit
restarts, and ``arnoldi_rtol`` as the relative tolerance on the Ritz values.
``arnoldi_restartdim`` only affects the builtin solver.
```

(英語版の csv-table の実際の列構成・既存文面を確認し、表現を合わせること)

- [ ] **Step 3: コミット**

```bash
git add docs/sphinx/ja/file_specification/correlation_length_section.rst docs/sphinx/en/file_specification/correlation_length_section.rst
git commit -m "Document correlation_length.eigensolver"
```

---

### Task 7: ベンチマークスイート correlation_length

**Files:**
- Create: `benchmark/templates/afh_square_clength.toml`
- Create: `benchmark/suites/correlation_length.toml`

**Interfaces:**
- Consumes: タイマーキー `"measure/correlation_length"`(Task 5)、`${eigensolver}` プレースホルダ(benchlib の string.Template 置換)

- [ ] **Step 1: テンプレートを作る**

`benchmark/templates/afh_square_clength.toml`
(`afh_square.toml` のコピーに correlation_length セクションを追加):

```toml
# Antiferromagnetic Heisenberg model on the square lattice,
# with the correlation-length measurement enabled.
# Placeholders: ${simple_steps} ${chi} ${L} ${W} ${D} ${eigensolver}
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

[correlation_length]
measure = true
# force the iterative (Arnoldi/ARPACK) path even for moderate chi
maxdim_dense_eigensolver = 50
eigensolver = "${eigensolver}"

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

- [ ] **Step 2: スイートを作る**

`benchmark/suites/correlation_length.toml`:

```toml
# A/B comparison of the transfer-matrix eigensolver
# (builtin Arnoldi vs ARPACK-NG).
# The "arpack" cases require a binary built with ARPACK-NG (ENABLE_ARPACK);
# without it they fail at input loading.
# Compare the measure/correlation_length timer between the two cases of the
# SAME run; `bench.py compare` is for commit-to-commit comparison.
[suite]
name = "correlation_length"
repeat = 3

[[case]]
name = "clength_${eigensolver}_D${D}_chi${chi}"
template = "../templates/afh_square_clength.toml"
params = { simple_steps = 100 }
sweep = { Lsub = [[2, 1]], D = [3], chi_ratio = [3], eigensolver = ["builtin", "arpack"] }
```

(chi = chi_ratio × D² = 27 → 転送行列の次元 N = 729 > 50 で反復経路に入る)

- [ ] **Step 3: スイートを実行して確認**

Run:
```bash
./venv/bin/python3 benchmark/bench.py run benchmark/suites/correlation_length.toml \
  --label arpack-dev --tenes-dir out-gcc/build/src --tool-dir out-gcc/build/tool --force
grep -rl 'measure/correlation_length' benchmark/results/arpack-dev/
grep -rh 'measure/correlation_length' benchmark/results/arpack-dev/ | head -4
```
(results ディレクトリの詳細レイアウトは `benchmark/README.md` を参照)
Expected: builtin / arpack 両ケースの timers.json がヒットし、どちらにも
`measure/correlation_length` がある。arpack ケースの時間が builtin と同程度以下。
`benchmark/results/` は git 管理外であることを確認してからコミットする
(`git status --short benchmark/` に results が出ないこと)。

- [ ] **Step 4: コミット**

```bash
git add benchmark/templates/afh_square_clength.toml benchmark/suites/correlation_length.toml
git commit -m "Add correlation_length benchmark suite (builtin vs ARPACK A/B)"
```

---

### Task 8: 仕上げ — NEWS・整形・両ビルド検証

**Files:**
- Modify: `NEWS.md`
- Modify: (clang-format の差分があれば該当 C++ ファイル)

- [ ] **Step 1: NEWS.md にエントリを追加**

`### Changes` の `tenes` 項に追加(PR 番号は PR 作成時に付けるので今は付けない):

```markdown
  - The transfer-matrix eigensolver for the correlation length can now use ARPACK-NG; TeNeS links it automatically when found (CMake option `ENABLE_ARPACK=AUTO/ON/OFF`, hint `ARPACK_ROOT`), and the solver can be chosen per run by `correlation_length.eigensolver = "auto" / "arpack" / "builtin"`
```

`### Development` に追加:

```markdown
- Added a `correlation_length` benchmark suite comparing the builtin Arnoldi and ARPACK-NG eigensolvers within one run
```

- [ ] **Step 2: clang-format を適用**

Run:
```bash
clang-format -i src/arpack_solver.hpp src/arpack_solver.cpp src/mpi.hpp src/iTPS/transfer_matrix.hpp src/iTPS/transfer_matrix.cpp src/iTPS/load_toml.cpp src/iTPS/correlation_length.cpp test/arpack.cpp test/util.cpp test/input.cpp
git diff --stat
```
Expected: 差分ゼロ、または整形差分のみ(あればそのままコミットに含める)

- [ ] **Step 3: ARPACK ありビルドで全 ctest**

Run: `cmake --preset gcc && cmake --build --preset gcc -j && ctest --preset gcc --output-on-failure | tail -5`
Expected: 100% tests passed

- [ ] **Step 4: ARPACK なしビルドで全 ctest**

Run:
```bash
SCRATCH=/private/tmp/claude-501/-Users-yomichi-source-github-com-issp-center-dev-TeNeS/dc8e166b-4e2e-4acd-9258-19d3d1f4f54b/scratchpad
cmake -S . -B $SCRATCH/tenes-noarpack -DENABLE_ARPACK=OFF -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_COMPILER=g++-16 -DTENES_PYTHON_EXECUTABLE=$PWD/venv/bin/python3 -DTesting=ON
SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk \
  cmake --build $SCRATCH/tenes-noarpack -j
OMP_NUM_THREADS=1 ctest --test-dir $SCRATCH/tenes-noarpack --output-on-failure | tail -5
```
Expected: 100% tests passed(test_arpack は登録されない。フォールバック経路の挙動不変を確認)

- [ ] **Step 5: コミット**

```bash
git add NEWS.md
git commit -m "Add NEWS entries for the ARPACK-NG eigensolver"
```

---

## 実行後(タスク外)

- PR 作成は `superpowers:finishing-a-development-branch` に従う(ブランチ `arpack_solver` → base `develop`)。PR 番号確定後に NEWS.md へ参照を追記する(リポジトリの慣行)。
- ベンチ結果の要約(builtin vs arpack の `measure/correlation_length` 比較)を PR 本文に載せる。
