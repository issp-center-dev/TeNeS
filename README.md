| Build status | Documentation (latest) | Documentation (latest stable) |
| :-: | :-: | :-: |
| [![Build status (master)](https://travis-ci.org/issp-center-dev/TeNeS.svg?branch=master)](https://travis-ci.org/issp-center-dev/TeNeS) | [![doc latest en](https://img.shields.io/badge/doc--en-master-blue.svg)](https://issp-center-dev.github.io/TeNeS/manual/master/en/html/index.html) [![doc latest ja](https://img.shields.io/badge/doc--ja-master-blue.svg)](https://issp-center-dev.github.io/TeNeS/manual/master/ja/html/index.html) | [![doc latest_stable en](https://img.shields.io/badge/doc--en-v0.1.0-blue.svg)](https://issp-center-dev.github.io/TeNeS/manual/v0.1.0/en/html/index.html) [![doc latest_stable ja](https://img.shields.io/badge/doc--ja-v0.1.0-blue.svg)](https://issp-center-dev.github.io/TeNeS/manual/v0.1.0/ja/html/index.html) |


# TeNeS

TeNeS (**Te**nsor **Ne**twork **S**olver) is a solver for 2D quantum lattice system based on a PEPS wave function and the CTM method.

## Online manual
- master (Latest, UNSTABLE)
    - [English](https://issp-center-dev.github.io/TeNeS/manual/master/en/html/index.html)
    - [日本語](https://issp-center-dev.github.io/TeNeS/manual/master/ja/html/index.html)
- v0.1 (Latest stable)
    - [English](https://issp-center-dev.github.io/TeNeS/manual/v0.1.0/en/html/index.html)
    - [日本語](https://issp-center-dev.github.io/TeNeS/manual/v0.1.0/ja/html/index.html)


## Getting started

- [Prerequisites](#prerequisites)
- [Install](#install)
    - [Simplest way to build](#simplest-way-to-build)
    - [Install](#install)
    - [Specify compiler](#specify-compiler)
    - [Use the pre-built mptensor](#use-the-pre-built-mptensor)
- [Usage](#usage)
- [Question or comment](#question-or-comment)
- [License](#license)
- [Acknowledgement](#acknowledgement)


## Prerequisites
The following tools are required for building TeNeS.

- C++11 compiler
- CMake (>=2.8.14)
- Python (>=2.7)
    - numpy
    - toml

TeNeS depends on the following libraries, but these are downloaded automatically through the build process.

- [mptensor](https://github.com/smorita/mptensor)
- [cpptoml](https://github.com/skystrife/cpptoml)
- [sanitizers-cmake](https://github.com/arsenm/sanitizers-cmake)

TeNeS can be parallerized by using MPI and ScaLAPACK.

## Install

### Simplest way to build

``` bash
$ mkdir build
$ cd build
$ cmake ../
$ make
```

The above commands makes an exectutable file `teses` in the `build/src` directory.

### Install

``` bash
$ cmake -DCMAKE_INSTALL_PREFIX=<path to install to> ../
$ make
$ make install
```

The above installs `tenes` into the `<path to install to>/bin` .
The default value of the `<path to install to>` is `/usr/local` .

### Specify compiler

CMake detects your compiler automatically but sometimes this is not what you want.
In this case, you can specify the compiler by the following way,

``` bash
$ cmake -DCMAKE_CXX_COMPILER=<path to your compiler> ../
```

### Specify Python interpreter

CMake detects also python interpreter automatically but sometimes this is not what you want.
In this case, you can specify the python interpreter by the following way,

``` bash
$ cmake -DPYTHON_EXECUTABLE=<path to your interpreter> ../
```

### Disable MPI/ScaLAPACK parallelization

To disable parallelization, pass the `-DENABLE_MPI=OFF` option to `cmake` commands.

### Use the pre-built mptensor

TeNeS is based on the parallerized tensor library, [mptensor](https://github.com/smorita/mptensor) .
The build system of TeNeS installs this automatically, but you can use the extra pre-built mptensor by the following way.

``` bash
$ cmake -DMPTENSOR_ROOT=<path to mptensor> ../
```

## Usage

``` bash
$ tenes input.toml
```

The file format of an input file is described in the manual page.

## Question or comment

Feel free to ask any question through an issue (public) or an e-mail (private) (`tenes-dev__at__issp.u-tokyo.ac.jp`, `__at__ -> @`).

Pull request is welcome (even for a small typo, of course!).

## License
TeNeS is available under the GNU GPL v3.

## Acknowledgement
TeNeS was supported by MEXT as “Exploratory Challenge on Post-K computer” (Frontiers of Basic Science: Challenging the Limits) and “Priority Issue on Post-K computer” (Creation of New Functional Devices and High-Performance Materials to Support Next-Generation Industries).
We also would also like to express our thanks for the support of the “Project for advancement of software usability in materials science” of The Institute for Solid State Physics, The University of Tokyo, for the development of TeNeS.
