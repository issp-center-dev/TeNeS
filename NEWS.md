# TeNeS v2.1.x Release Notes

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
