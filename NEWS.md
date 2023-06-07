# TeNeS v1.3.x Release Notes

## Changes between v1.3.2 and v1.3.1

### Bugfixes

- `tenes_simple`
    - Antiferromagnetic initial state for S>1/2 is wrong (#61 #62)
    - Bilinear biquadratic term $B S \cdot S$ is wrong (#63)
    - `tenes_simple` fails to treat longer interaction and field simultaneously unless using `--use-site-hamiltonian` option (#66)

## Changes between v1.3.1 and v1.3.0

- Update bundled doctest.h
    - New version supports Intel compiler 19

## Changes between v1.3.0 and v1.2

### New Features

- `tenes_simple`
    - Enable to generate site hamiltonian (#56)
- `tenes_std`
    - Introduce site hamiltonian (#56)
    - Copy `hamiltonian` as `observable` with `group=0` when `observable` with `group=0` is not defined (#53)
- `tenes`
    - Introduce site imaginary time evolutionary operator (#56)

### Fixes

- `tenes_simple`
    - Fix a bug that an initial state always becomes "antiferro"  in the case of honeycomb lattice (#59)
- `tenes`
    - Undefined observables (skipped groups) will no longer be written in `density.dat` (#60)

### Documents and tutorials

- Update (#52, #55, #57, #58)

### For developers

- Fix bugs that compilation fails in some environments (#47, #51)
- CMake can pass `MPIEXEC_PREFLAGS` and `MPIEXEC_POSTFLAGS` options to `mpiexec` commands in test (#48)
