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

#ifndef TENES_SRC_ITPS_MAIN_HPP_
#define TENES_SRC_ITPS_MAIN_HPP_

#include <string>

#include "../printlevel.hpp"
#include "../mpi.hpp"

extern "C" {
int tenes_itps_main(const char* input_filename, MPI_Comm comm, int print_level);
}

namespace tenes {
namespace itps {
int itps_main(const char* input_filename, MPI_Comm comm,
              PrintLevel print_level = PrintLevel::info);
int itps_main(std::string input_filename, MPI_Comm comm,
              PrintLevel print_level = PrintLevel::info);
}  // namespace itps
}  // namespace tenes

#endif  // TENES_SRC_ITPS_MAIN_HPP_
