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

#include <iostream>
#include <string>

#include "version.hpp"
#include "mpi.hpp"
#include "exception.hpp"
#include "printlevel.hpp"

#include "tensor.hpp"
#include "iTPS/main.hpp"

int main2(int argc, char **argv) {
  int mpirank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int status = 0;

  std::string usage =
      R"(TeNeS: TEnsor NEtwork Solver for 2D quantum lattice system
  
  Usage:
    tenes [--quiet] <input_toml>
    tenes --help
    tenes --version

  Options:
    -h --help       Show this help message.
    -v --version    Show the version.
    -q --quiet      Do not print any messages.
  )";

  if (argc == 1) {
    if (mpirank == 0) std::cout << usage << std::endl;
    return 0;
  }

  for (int i = 1; i < argc; ++i) {
    std::string opt = argv[i];
    if (opt == "-h" || opt == "--help") {
      if (mpirank == 0) std::cout << usage << std::endl;
      return 0;
    }
  }

  for (int i = 1; i < argc; ++i) {
    std::string opt = argv[i];
    if (opt == "-v" || opt == "--version") {
      if (mpirank == 0) std::cout << "TeNeS v" << TENES_VERSION << std::endl;
      return 0;
    }
  }

  using PrintLevel = tenes::PrintLevel;

  PrintLevel print_level = PrintLevel::info;
  std::string input_filename;
  for (int i = 1; i < argc; ++i) {
    std::string opt = argv[i];
    if (opt == "-q" || opt == "--quiet") {
      print_level = PrintLevel::none;
    } else {
      input_filename = opt;
    }
  }

  try {
    status =
        tenes::itps::itps_main(input_filename, MPI_COMM_WORLD, print_level);
  } catch (const tenes::input_error e) {
    if (mpirank == 0) {
      std::cerr << "[INPUT ERROR]" << std::endl;
      std::cerr << e.what() << std::endl;
    }
    status = 1;
  } catch (const tenes::load_error e) {
    if (mpirank == 0) {
      std::cerr << "[TENSOR LOAD ERROR]" << std::endl;
      std::cerr << e.what() << std::endl;
    }
    status = 1;
  } catch (const tenes::runtime_error e) {
    if (mpirank == 0) {
      std::cerr << "[ERROR]" << std::endl;
      std::cerr << e.what() << std::endl;
    }
    status = 1;
  }

  return status;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  tenes::initialize_mptensor();
  int status = main2(argc, argv);
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Comm parent;
  MPI_Comm_get_parent(&parent);
  if (parent != MPI_COMM_NULL) {
    int mpirank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    if (mpirank == 0) {
      MPI_Send(&status, 1, MPI_INT, 0, 0, parent);
    }
    MPI_Comm_disconnect(&parent);
  }
  MPI_Finalize();
  return status;
}

/*! \mainpage notitle
 *
 * ## About TeNeS
 *
 * [TeNeS (TEnsor NEtwork Solver)](https://github.com/issp-center-dev/TeNeS) is
 * a solver for 2D quantum lattice system based on a PEPS wave function and the
 * CTM method. TeNeS can make use of many CPU/nodes through an OpenMP/MPI
 * hybirid parallel tensor operation library,
 * [mptensor](https://github.com/smorita/mptensor).
 *
 * ## About this document
 *
 * This document, for developers, describes classes and functions in TeNeS.
 * For users, see the [user's
 * manual](https://issp-center-dev.github.io/TeNeS/manual/develop/en/html/index.html).
 *
 * ## Architecture
 *
 * \ref tenes is the main program written in C++.
 *
 * \ref tenes_simple and \ref tenes_std are the utility tools for generating
 * input files written in Python.
 *
 */
