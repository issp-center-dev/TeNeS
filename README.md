![TeNeS logo](docs/sphinx/img/TeNeS_logo_banner.png)

| Branch | Build status | Documentation |
| :-: | :-: | :-: |
| master (latest stable) | [![master](https://travis-ci.org/issp-center-dev/TeNeS.svg?branch=master)](https://travis-ci.org/issp-center-dev/TeNeS) | [![doc_en](https://img.shields.io/badge/doc-English-blue.svg)](https://issp-center-dev.github.io/TeNeS/manual/master/en/html/index.html) [![doc_ja](https://img.shields.io/badge/doc-Japanese-blue.svg)](https://issp-center-dev.github.io/TeNeS/manual/master/ja/html/index.html) |
| develop (latest) | [![develop/Test](https://github.com/issp-center-dev/TeNeS/workflows/Test/badge.svg?branch=develop)](https://github.com/issp-center-dev/TeNeS/actions?query=workflow%3ATest+branch%3Adevelop) | [![doc_en](https://img.shields.io/badge/doc-English-blue.svg)](https://issp-center-dev.github.io/TeNeS/manual/develop/en/html/index.html) [![doc_ja](https://img.shields.io/badge/doc-Japanese-blue.svg)](https://issp-center-dev.github.io/TeNeS/manual/develop/ja/html/index.html) |

# TeNeS

TeNeS (**Te**nsor **Ne**twork **S**olver) is a solver for 2D quantum lattice system based on a PEPS wave function and the CTM method.
TeNeS can make use of many CPU/nodes through an OpenMP/MPI hybirid parallel tensor operation library, [mptensor](https://github.com/smorita/mptensor).

## Online manual

- develop (Latest, UNSTABLE)
    - [English](https://issp-center-dev.github.io/TeNeS/manual/develop/en/html/index.html)
    - [日本語](https://issp-center-dev.github.io/TeNeS/manual/develop/ja/html/index.html)
- master (Latest, stable)
    - [English](https://issp-center-dev.github.io/TeNeS/manual/master/en/html/index.html)
    - [日本語](https://issp-center-dev.github.io/TeNeS/manual/master/ja/html/index.html)

## Getting started

- [Prerequisites and dependencies](#prerequisites-and-dependencies)
- [Install](#install)
    - [Simplest way to build](#simplest-way-to-build)
    - [Install binaries and samples](#install-binaries-and-samples)
    - [Specify compiler](#specify-compiler)
    - [Disable MPI/ScaLAPACK parallelization](#disable-mpi/scalapack-parallelization)
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

TeNeS depends on the following libraries, but these are downloaded automatically through the build process.

- [mptensor](https://github.com/smorita/mptensor)
- [cpptoml](https://github.com/skystrife/cpptoml)

TeNeS can be parallerized by using MPI and ScaLAPACK.

TeNeS tools (`tenes_simple`, `tenes_std`) are written in Python3.
The following external packages are required:

- numpy
- scipy
- toml
- typing (mandatory for python < 3.5)

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

For example, the following file `simple.toml` represents the transverse field Ising model on the square lattice.

```
[parameter]
[parameter.general]
is_real = true

[parameter.simple_update]
num_step = 1000
tau = 0.01

[parameter.full_update]
num_step = 0
tau = 0.01

[parameter.ctm]
iteration_max = 10
dimension = 10

[lattice]
type = "square lattice"
L = 2
W = 2
virtual_dim = 2
initial = "ferro"

[model]
type = "spin"
Jz = -1.0 # negative for FM interaction
Jx = 0.0
Jy = 0.0
hx = 1.0   # transverse field
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

```
Sz          =  2.97866964051826333e-01  0.00000000000000000e+00
Sx          =  3.86024172907023511e-01  0.00000000000000000e+00
hamiltonian = -7.57303058659582140e-01  0.00000000000000000e+00
SzSz        =  2.16869216589772901e-01  0.00000000000000000e+00
SxSx        =  3.19350111777505108e-01  0.00000000000000000e+00
SySy        = -4.77650003168152704e-02  0.00000000000000000e+00
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

## Acknowledgement

TeNeS was supported by MEXT as "Exploratory Challenge on Post-K computer" (Frontiers of Basic Science: Challenging the Limits) and "Priority Issue on Post-K computer" (Creation of New Functional Devices and High-Performance Materials to Support Next-Generation Industries).
We also would also like to express our thanks for the support of the "Project for advancement of software usability in materials science" of The Institute for Solid State Physics, The University of Tokyo, for the development of TeNeS.
