![TeNeS logo](docs/sphinx/img/TeNeS_logo_banner.png)

| Branch | Build status | Documentation |
| :-: | :-: | :-: |
| master (latest stable) | [![master][ci/master/badge]][ci/master/uri] | [![doc_en][doc/en/badge]][doc/master/en/uri] [![doc_ja][doc/ja/badge]][doc/master/ja/uri] |
| develop (latest) | [![develop][ci/develop/badge]][ci/develop/uri] | [![doc_en][doc/en/badge]][doc/develop/en/uri] [![doc_ja][doc/ja/badge]][doc/develop/ja/uri] |

# TeNeS

TeNeS (**Te**nsor **Ne**twork **S**olver) is a solver for 2D quantum lattice system based on a PEPS wave function and the CTM method.
TeNeS can make use of many CPU/nodes through an OpenMP/MPI hybirid parallel tensor operation library, [mptensor](https://github.com/smorita/mptensor).

## Online manual

- develop (Latest, UNSTABLE)
  - [English][doc/develop/en/uri]
  - [日本語][doc/develop/ja/uri]
- master (Latest, stable)
  - [English][doc/master/en/uri]
  - [日本語][doc/master/ja/uri]

## Getting started

- [Prerequisites and dependencies](#prerequisites-and-dependencies)
- [Install](#install)
  - [Simplest way to build](#simplest-way-to-build)
  - [Install binaries and samples](#install-binaries-and-samples)
  - [Specify compiler](#specify-compiler)
  - [Disable MPI/ScaLAPACK parallelization](#disable-mpiscalapack-parallelization)
  - [Specify ScaLAPACK](#specify-scalapack)
  - [Use the pre-built mptensor](#use-the-pre-built-mptensor)
  - [Specify Python interpreter](#specify-python-interpreter)
- [Usage](#usage)
  - [Use pre-defined model and lattice](#use-pre-defined-model-and-lattice)
  - [Calculate imaginary time evolution operators](#calculate-imaginary-time-evolution-operators)
  - [Perform](#perform)
- [Question or comment](#question-or-comment)
- [Contibution](#contibution)
- [License](#license)
- [Acknowledgement](#acknowledgement)

## Prerequisites and dependencies

The following tools are required for building TeNeS.

- C++11 compiler
- CMake (>=3.6.0)
- BLAS/LAPACK

TeNeS depends on the following libraries, but these are downloaded automatically through the build process.

- [mptensor](https://github.com/smorita/mptensor)
- [cpptoml](https://github.com/skystrife/cpptoml)

TeNeS can be parallerized by using MPI and ScaLAPACK.

TeNeS tools (`tenes_simple`, `tenes_std`) are written in Python3.
The following external packages are required:

- numpy
- scipy
- toml

## Install

### Simplest way to build

``` bash
mkdir build
cd build
cmake ../
make
```

(NOTE: Some system (e.g. CentOS) provides CMake 3 as `cmake3`)

The above commands makes an exectutable file `tenes` in the `build/src` directory.

### Install binaries and samples

``` bash
cmake -DCMAKE_INSTALL_PREFIX=<path to install to> ../
make
make install
```

Noted that the parallel building `make -j <num_parallel>` can reduce the time to build.

The `make install` command installs `tenes`, `tenes_std`, and `tenes_simple` into the `<path to install to>/bin` .
Samples will be also installed into the `<path to install to>/share/tenes/<VERSION>/sample` .
The default value of the `<path to install to>` is `/usr/local` .

### Specify compiler

CMake detects your compiler automatically but sometimes this is not what you want.
In this case, you can specify the compiler by the following way,

``` bash
cmake -DCMAKE_CXX_COMPILER=<path to your compiler> ../
```

### Disable MPI/ScaLAPACK parallelization

To disable parallelization, pass the `-DENABLE_MPI=OFF` option to `cmake` commands.

If you use macos, MPI/ScaLAPACK parallelization is disabled by default because the combination of Apple Accelerate BLAS/LAPACK library with ScaLAPACK seems to have some troubles.

### Specify ScaLAPACK

TeNeS finds ScaLAPACK automatically, but may fail.
In such a case, `-DSCALAPACK_ROOT=<path>` option specifies the path to the ScaLAPACK library file, `<path>/lib/libscalapack.so`.

### Use the pre-built mptensor

TeNeS is based on the parallerized tensor library, [mptensor](https://github.com/smorita/mptensor) (>= v0.3).
The build system of TeNeS installs this automatically, but you can use the extra pre-built mptensor by the following way.

``` bash
cmake -DMPTENSOR_ROOT=<path to mptensor> ../
```

### Specify Python interpreter

TeNeS tools `tenes_simple` and `tenes_std` use `python3` which can be found in `PATH` via `/usr/bin/env python3`.
Please make sure that `python3` command invokes Python3 interpreter, for example, by using `type python3` .

If you want to fix the interpreter to be used (or `/usr/bin/env` does not exist), you can specify it by the following way,

``` bash
cmake -DTENES_PYTHON_EXECUTABLE=<path to your interpreter> ../
```

## Usage

### Use pre-defined model and lattice

For example, the following file `simple.toml` represents the transverse field Ising model on the square lattice ([sample/01_transverse_field_ising/simple.toml](sample/01_transverse_field_ising/simple.toml)).

``` toml
[parameter]
[parameter.general]
is_real = true  # Limit tensor elements in real (not complex)

[parameter.simple_update]
num_step = 1000 # Number of steps in simple update
tau = 0.01      # Imaginary time slice

[parameter.full_update]
num_step = 0    # Number of steps in full update
tau = 0.01      # Imaginary time slice

[parameter.ctm]
meanfield_env = false # Use meanfield environment to contract iTNS
iteration_max = 10    # Maximum number of iterations in CTMRG
dimension = 10        # Bond dimension of corner transfer matrix

[lattice]
type = "square lattice" # Type of lattice
L = 2                   # X length of unit cell
W = 2                   # Y length of unit cell
virtual_dim = 2         # Bond dimension of bulk tensors
initial = "ferro"       # Initial condition

[model]
type = "spin" # Type of model
Jz = -1.0     # Jz SzSz
Jx = 0.0      # Jx SxSx
Jy = 0.0      # Jy SySy
hx = 0.0      # hx Sx

```

`tenes_simple` is a utility tool for converting this file to another file, `std.toml`, denoting the operator tensors including bond hamiltonian.

``` bash
tenes_simple simple.toml
```

### Calculate imaginary time evolution operators

`tenes_std` is another utility tool for calculating imaginary time evolution operators and converting `std.toml` to the input file of `tenes`, `input.toml`.

``` bash
tenes_std std.toml
```

By editing `std.toml`, users can perform other models and lattices as ones like.

### Perform

To perform simulation, pass `input.toml` to `tenes` as the following

``` bash
tenes input.toml
```

Results can be found in `output` directory.
For example, expectation values of operators per site are stored in `output/densities.dat` as the following,

``` text
Energy           = -5.00000000000000000e-01  0.00000000000000000e+00
Sz               =  5.00000000000000000e-01  0.00000000000000000e+00
Sx               = -1.28526262481784176e-13  0.00000000000000000e+00
bond_hamiltonian = -5.00000000000000000e-01  0.00000000000000000e+00
SzSz             =  5.00000000000000000e-01  0.00000000000000000e+00
SxSx             = -1.73749199820346390e-18  0.00000000000000000e+00
SySy             =  1.73749202732501919e-18  0.00000000000000000e+00

```

The file format of input/output files is described in the manual page.

## Question or comment

Feel free to ask any question through an issue (public) or an e-mail (private) (`tenes-dev__at__issp.u-tokyo.ac.jp`, `__at__ -> @`).

## Contibution

Pull request is welcome (even for a small typo, of course!).
Before send a PR, please make sure the following:

- Rebase (or merge) `develop` branch into your feature branch
- Check `make` and `ctest` processes pass
- Format Codes by using `clang-format` (C++) and `black` (Python)

(Incomplete) developer's document written in doxygen is available.

1. Move to `docs/doxygen`
2. Invoke `doxygen`
3. Open `doxygen_out/html/index.html` in your browser

## License

TeNeS is available under the GNU GPL v3.

## Paper

When you publish the results by using TeNeS, we would appreciate if you cite the following paper:

- [Y. Motoyama, Tsuyoshi Okubo, Kazuyoshi Yoshimi, Satoshi Morita, Takeo Kato, and Naoki Kawashima, "TeNeS: Tensor Network Solver for Quantum Lattice Systems", Comput. Phys. Commun. **279**, 108437 (2022)][paper/tenes_1.2]
- [Y. Motoyama, Tsuyoshi Okubo, Kazuyoshi Yoshimi, Satoshi Morita, Tatsumi Aoyama, Takeo Kato, and Naoki Kawashima, "TeNeS-v2: Enhancement for real-time and finite temperature simulations of quantum many-body systems", Comput. Phys. Commun. **315**, 109692 (2025)][paper/tenes_2.1]

## Acknowledgement

TeNeS was supported by MEXT as "Exploratory Challenge on Post-K computer" (Frontiers of Basic Science: Challenging the Limits) and "Priority Issue on Post-K computer" (Creation of New Functional Devices and High-Performance Materials to Support Next-Generation Industries).
We also would also like to express our thanks for the support of the "Project for advancement of software usability in materials science" of The Institute for Solid State Physics, The University of Tokyo, for the development of TeNeS.

[ci/master/badge]: https://github.com/issp-center-dev/TeNeS/actions/workflows/linux.yml/badge.svg?branch=master
[ci/master/uri]: https://github.com/issp-center-dev/TeNeS/actions?query=workflow%3ATest+branch%3Amaster
[ci/develop/badge]: https://github.com/issp-center-dev/TeNeS/actions/workflows/linux.yml/badge.svg?branch=develop
[ci/develop/uri]: https://github.com/issp-center-dev/TeNeS/actions?query=workflow%3ATest+branch%3Adevelop

[doc/en/badge]: https://img.shields.io/badge/doc-English-blue.svg
[doc/ja/badge]: https://img.shields.io/badge/doc-Japanese-blue.svg
[doc/master/en/uri]: https://issp-center-dev.github.io/TeNeS/manual/master/en/html/index.html
[doc/master/ja/uri]: https://issp-center-dev.github.io/TeNeS/manual/master/ja/html/index.html
[doc/develop/en/uri]: https://issp-center-dev.github.io/TeNeS/manual/develop/en/html/index.html
[doc/develop/ja/uri]: https://issp-center-dev.github.io/TeNeS/manual/develop/ja/html/index.html

[paper/tenes_1.2]: https://www.sciencedirect.com/science/article/pii/S0010465522001564
[paper/tenes_2.1]: https://www.sciencedirect.com/science/article/pii/S0010465525001948
