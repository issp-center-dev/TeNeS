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

#include "../../../tensor.hpp"

#ifndef INSTANTIATE_CONTRACT
#define INSTANTIATE_CONTRACT(tensor_type, NROW, NCOL)                \
  template tensor_type::value_type Contract_##NROW##x##NCOL##_iTPS_MF( \
      const std::vector<std::vector<const tensor_type *>> &Tn,       \
      const std::vector<std::vector<const tensor_type *>> &op)
#endif  // INSTANTIATE_CONTRACT
