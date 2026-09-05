# TeNeS Release Notes

## Changes between develop branch (2026-07-11) and v2.1.3

### New features

- `tenes`
  - Experimental support for fermionic models (`[parameter.general] fermion = true` with per-site `[[tensor.unitcell]] parity`): site tensors become Z2-graded and all bond updates and measurements generate the fermionic exchange signs automatically; verified against exact free-fermion results (energy within 0.7% at D=4, chi=32 for half-filled spinless fermions on the square lattice)
  - The full update is available in fermion mode with the CTM environment (`fastfullupdate = true`, the default, falls back to the plain full update with a warning; `meanfield_env = true` with a full update is rejected). The two-site environment is built by an open-channel fold of the graded QR factors and the bosonic ALS core is reused on masked variables; the parity violation of that environment is checked against `max(1e-8, convergence_epsilon)` of the CTM. The fermionic CTM environment carries a window-dependent overall phase, so the RDM convergence check is made phase-invariant in fermion mode and the full-update metric, the one-site RDM output and the norm diagnostics fix their own phase. Note that at small bond dimensions (e.g. `D = 2` for free fermions) the full update started from a converged simple-update state can raise the energy; `D >= 3` is recommended
  - Limitations of the current fermion mode: ground-state mode only (simple update with the CTM or mean-field environment, full update with the CTM environment); `Use_RSVD`, `Simple_Gauge_Fix`, finite temperature, real-time evolution, multi-site observables, `ops`-form two-site observables, two-site observables beyond nearest neighbors, skewed unit cells, one-wide cells (`LX < 2` or `LY < 2`), correlation functions, and correlation length are rejected at input time (the correlation length is forcibly disabled); away from half filling the simple update needs considerably more imaginary time to converge
  - `tensor_save` / `tensor_load` now work in fermion mode, so a parameter scan can restart from a converged state. The parity ledger of the virtual bonds is saved alongside the tensors as `fermion.dat` and validated on load (physical parity, `virtual_dim`, `L_sub`, `skew`, and the parity of the tensors themselves); restarting with a different `virtual_dim` is not supported. The lambda files are now written with full double precision instead of the default six significant digits.
  - `meanfield_env = true` is now accepted in fermion mode; two-site observables are then evaluated by a single-layer graded contraction (cost D^6 d^2 instead of the D^12 d^4 reduced-pair blob of the CTM path), which makes D > 2 and the Hubbard model (d = 4) measurable at simple-update accuracy
  - `tenes` adds an identity simple-update gate on every virtual bond with dimension > 1 that no evolution operator acts on (reported as an INFO line). Such a bond was never touched by the simple update, kept its random initial tensors, and the CTM then often failed to converge on it (in fermion mode this showed up as negative norms and energies below the variational bound); the identity gate is physically a no-op but canonicalises the bond every step
- `tenes_simple` / `tenes_std`
  - `tenes_simple` now supports two fermionic models on the square lattice with nearest-neighbor bonds: `type = "spinless fermion"` (spinless fermions; parameters `t`, `v`, `mu`) and `type = "hubbard"` (the fermionic Hubbard model; parameters `t`, `u`, `v`, `mu`, `h`; the Bose-Hubbard model remains `type = "boson"`), with initial states `vacuum` (both), `full` and `cdw` (Hubbard); `tenes_std` carries the `parity` metadata through to `input.toml` and rejects inputs the fermion mode cannot handle (long-range bonds, `ops`-form two-site observables, skewed unit cells)

### Changes

- Building TeNeS now requires a C++17 compiler and mptensor v0.5.0 or later ([#104][], [#109][])
- The minimum required CMake version is raised from 3.6 to 3.8, so that the language standard is honored in configure-time checks (CMP0067); note that the bundled toml11 requires CMake 3.16 or later anyway ([#112][])
- On macOS, CMake now locates a Homebrew libomp instead of assuming `/usr/local`: it queries `brew --prefix libomp` and falls back to Homebrew's usual prefixes, then lets `find_package(OpenMP)` drive Apple clang. Building with Apple clang -- including `-DCMAKE_CXX_COMPILER=mpicxx`, since Homebrew's `mpicxx` wraps it -- therefore requires CMake 3.12 or later, which is when `OpenMP_ROOT` started to steer the library search ([#113][])
- `tenes`
  - Replaced the TOML parser cpptoml (archived, TOML v0.5.0) with toml11 v4.4.0 (TOML v1.0.0); input errors now report the file name, the line number, and the offending value ([#108][])
  - Writing `[observable.onesite]` and similar sections as a single table instead of an array of tables (`[[...]]`) is now an input error instead of being silently ignored ([#108][])
  - CMake option `CPPTOML_ROOT` is renamed to `TOML11_ROOT` ([#108][])
  - Now always writes `output/timers.json`: cumulative wall times per hierarchical timer name (total / phase/* / contract/*), with per-MPI-rank max/min ([#110][])
  - The transfer-matrix eigensolver for the correlation length can now use ARPACK-NG; TeNeS links it automatically when found (CMake option `ENABLE_ARPACK=AUTO/ON/OFF`, hint `ARPACK_ROOT`), and the solver can be chosen per run by `correlation_length.eigensolver = "auto" / "arpack" / "builtin"` ([#111][])
  - `correlation_length.arnoldi_maxdim`, `arnoldi_restartdim`, and `arnoldi_maxiterations` now default to 0, meaning a solver-dependent automatic value scaled by `num_eigvals` (ARPACK-NG: small subspace with restarts; builtin: a large single sweep reproducing the previous defaults); explicit values are used as-is ([#111][])
  - When ARPACK-NG fails to converge with an automatic `arnoldi_maxdim`, the Krylov subspace is doubled (up to the matrix size) and the computation retried; with an explicit `arnoldi_maxdim`, unconverged eigenvalues and the resulting correlation length are reported as NaN instead of a misleading finite value ([#111][])

### Bug fixes

- `tenes`
  - Fixed a link error (undefined reference to `Make_single_tensor_density`) with some compilers, e.g. GCC 13 on Linux ([#107][])
  - Fixed Inf/NaN handling with the Intel icpx compiler: the default `-fp-model=fast` broke the detection of divergent correlation lengths and could leak `nan` rows of unmeasured sites into `onesite_obs.dat` / `density.dat`; `-fp-model=precise` is now enforced for icpx builds ([#107][])
  - Fixed an undefined behavior (null-pointer dereference) when the input file has no `[[evolution.simple]]` / `[[evolution.full]]` sections; they are now treated as an empty list of evolution operators ([#108][])
  - Fixed an out-of-bounds access in mptensor (`make_l2g_map`) that aborted finite-temperature calculations in Debug builds, by updating mptensor to v0.5.0 ([#104][])
  - Fixed an out-of-bounds access in mptensor (`&v[0]` on an empty `std::vector`) that aborted every MPI run with two or more processes in Debug builds. Empty local blocks are routine with the ScaLAPACK backend, and GCC 15 and later enable `_GLIBCXX_ASSERTIONS` by default when compiling without optimization, turning the previously harmless undefined behavior into an abort ([#113][])
  - `tenes --version` / `--help` no longer initialize MPI, fixing a hang of MPI-linked binaries where the launcher infrastructure is unavailable, e.g. on cluster login nodes ([#111][])

### Development

- Modernized the C++ core to C++17: removed the vendored boost headers (`boost::optional` → `std::optional`), reimplemented `util/file` with `std::filesystem`, and replaced SFINAE-based type traits with `if constexpr` ([#109][])
- CMake now probes at configure time whether `std::filesystem` needs a separate library (`stdc++fs` / `c++fs`) and links it automatically; this fixes a link error (undefined reference to `std::filesystem::status` etc.) with icpx or clang picking up GCC 8's libstdc++, e.g. on RHEL 8 without a gcc-toolset ([#112][])
- Updated the bundled mptensor to v0.5.0 ([#104][])
- Repaired the macOS CI: install `libomp` for mptensor's OpenMP requirement ([#104][])
- Added a `benchmark/` harness for A/B performance comparison (`bench.py run` / `bench.py compare`); see `benchmark/README.md` ([#110][])
- Added a `correlation_length` benchmark suite comparing the builtin Arnoldi and ARPACK-NG eigensolvers within one run ([#111][])
- The benchmark harness now applies a 10-second timeout when collecting `tenes --version` for provenance, so a non-returning binary cannot wedge the whole run ([#111][])
- Added `bench.py show <label-dir>` rendering the results of a single run, pairing within-run A/B cases (builtin vs arpack) side by side with a ratio column ([#111][])

### Documentation and samples

- Updated the requirements in README and the install guides: C++17, mptensor >= v0.5.0, dependency list, and notes for Intel compilers ([#104][], [#107][], [#108][])
- Install guides: how to set up libomp for Apple clang on macOS, a note that Homebrew's `open-mpi` + `scalapack` works with `-DENABLE_MPI=ON` (that ScaLAPACK links OpenBLAS, not Accelerate), and a recommendation to run with `OMP_NUM_THREADS=1` on macOS, where every OpenMP barrier costs a system call ([#113][])

## Changes between v2.1.3 and v2.1.2

### Changes

- `tenes`
  - `group` keyword in `evolution.simple`/`evolution.full` sections is now mandatory ([#101][])

### Bug fixes

- `tenes`
  - Fixed an out-of-bounds access in the convergence check of CTMRG ([#103][])
  - Fixed `tau` of the full update being read from the `parameter.simple_update` section instead of `parameter.full_update` ([#106][])
  - Fixed a broken identity matrix generation that disabled the shift in the Arnoldi restart ([#106][])
  - Fixed a hang in the Arnoldi solver when `Lanczos_restartdim` exceeds the transfer-matrix dimension ([#106][])
  - Fixed crashes and silent data corruption when loading saved tensors from missing or truncated files ([#106][])
  - Fixed a crash when no evolution operators are defined ([#106][])
  - Fixed correlation length output for degenerate transfer-matrix eigenvalues ([#106][])
- `tenes_simple`
  - Fixed the onsite repulsion `U` of the Bose-Hubbard model being silently dropped when `--use-site-hamiltonian` is specified ([#105][])
  - Fixed the names of the two-site (hopping) operators of the Bose-Hubbard model ([#106][])
  - Fixed an `AttributeError` when a unit-cell site is left undefined ([#106][])
- `tenes_std`
  - Optional `parameter` subsections are no longer required in the input file ([#106][])
  - Fixed merging of `parameter` subsections when multiple input files are given ([#106][])

### Development

- Updated and pinned versions of GitHub Actions ([#102][])
- Added unit tests for the fixes in [#106][] (doctest suites and a pytest suite, registered in CTest)

### Documentation and samples

- Added new paper information

## Changes between v2.1.2 and v2.1.1

### Bug fixes

- Fixed missing group parameter in path-decomposed evolution operators ([#100][])

### Development

- Updated CI environments ([#99][])

### Documentation and samples

- Updated README
- Updated samples of QMC in finite temperature calculations
- Fixed plot scripts in finite temperature
- Updated plot scripts for finite temperature calculation
- Updated sample for finite temperature calculation

## Changes between v2.1.1 and v2.1.0

### Bug fixes

- `tenes_std`
  - Fixed a bug in the calculation of the evolutionary tensor for complex Hamiltonians ([#97][])

## Changes between v2.1.0 and v2.0.0

### New features

- `tenes`
  - Enabled to save tensors in the real-time evolution and the finite-temperature calculation ([#88][])
  - Enabled to specify coefficient of observables ([#91][])
- `tenes_std`
  - Enabled to read multiple input files ([#92][])

### Bug fixes

- `tenes`
  - Fixed a bug of combination of RSVD and CTMRG method with shrinkage of chi ([#86][], [#87][])

[#86]: https://github.com/issp-center-dev/TeNeS/pull/86
[#87]: https://github.com/issp-center-dev/TeNeS/pull/87
[#88]: https://github.com/issp-center-dev/TeNeS/pull/88
[#91]: https://github.com/issp-center-dev/TeNeS/pull/91
[#92]: https://github.com/issp-center-dev/TeNeS/pull/92
[#97]: https://github.com/issp-center-dev/TeNeS/pull/97
[#99]: https://github.com/issp-center-dev/TeNeS/pull/99
[#100]: https://github.com/issp-center-dev/TeNeS/pull/100
[#101]: https://github.com/issp-center-dev/TeNeS/pull/101
[#102]: https://github.com/issp-center-dev/TeNeS/pull/102
[#103]: https://github.com/issp-center-dev/TeNeS/pull/103
[#104]: https://github.com/issp-center-dev/TeNeS/pull/104
[#105]: https://github.com/issp-center-dev/TeNeS/pull/105
[#106]: https://github.com/issp-center-dev/TeNeS/pull/106
[#107]: https://github.com/issp-center-dev/TeNeS/pull/107
[#108]: https://github.com/issp-center-dev/TeNeS/pull/108
[#109]: https://github.com/issp-center-dev/TeNeS/pull/109
[#110]: https://github.com/issp-center-dev/TeNeS/pull/110
[#111]: https://github.com/issp-center-dev/TeNeS/pull/111
[#112]: https://github.com/issp-center-dev/TeNeS/pull/112
[#113]: https://github.com/issp-center-dev/TeNeS/pull/113
