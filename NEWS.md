# TeNeS v2.1.x Release Notes

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
[#105]: https://github.com/issp-center-dev/TeNeS/pull/105
[#106]: https://github.com/issp-center-dev/TeNeS/pull/106
