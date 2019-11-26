[![Build status (master)](https://travis-ci.org/issp-center-dev/TeNeS.svg?branch=master)](https://travis-ci.org/issp-center-dev/TeNeS)

# TeNeS

TeNeS (**Te**nsor **Ne**twork **S**olver) is a solver for 2D quantum lattice system based on a PEPS wave function and the CTM method.

## Table of Contents

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
- MPI and ScaLAPACK

TeNeS depends on the following libraries, but these are downloaded automatically through the build process.

- [mptensor](https://github.com/smorita/mptensor)
- [cpptoml](https://github.com/skystrife/cpptoml)
- [sanitizers-cmake](https://github.com/arsenm/sanitizers-cmake)

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

The file format of an input file is described in the manual page (under the construction).

## Question or comment

Feel free to ask any question through an issue (public) or an e-mail (private) (`tenes-dev__at__issp.u-tokyo.ac.jp`, `__at__ -> @`).

Pull request is welcome (even for a small typo, of course!).

## License
TeNeS is available under the GNU GPL v3.

## Acknowledgement
TeNeS is developed under the support of "Project for advancement of software usability in materials science" in fiscal year 2019 by The Institute for Solid State Physics, The University of Tokyo.

