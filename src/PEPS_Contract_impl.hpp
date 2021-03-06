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

#ifndef SRC_PEPS_CONTRACT_IMPL_HPP_
#define SRC_PEPS_CONTRACT_IMPL_HPP_

#include <vector>
#include <mptensor/tensor.hpp>
#include <mptensor/complex.hpp>

namespace tenes {
using namespace mptensor;

/*
Declare template functions like

template <class tensor>
typename tensor::value_type
Contract_1x1 (
  const std::vector<const tensor*> &C,
  const std::vector<const tensor*> &eTt,
  const std::vector<const tensor*> &eTr,
  const std::vector<const tensor*> &eTb,
  const std::vector<const tensor*> &eTl,
  const std::vector<std::vector<const tensor*>> &Tn,
  const std::vector<std::vector<const tensor*>> &op
)
*/

#define DECLARE_CONTRACT(NROW, NCOL)                      \
  template <class tensor>                                 \
  typename tensor::value_type Contract_##NROW##x##NCOL(   \
      const std::vector<const tensor *> &C,               \
      const std::vector<const tensor *> &eTt,             \
      const std::vector<const tensor *> &eTr,             \
      const std::vector<const tensor *> &eTb,             \
      const std::vector<const tensor *> &eTl,             \
      const std::vector<std::vector<const tensor *>> &Tn, \
      const std::vector<std::vector<const tensor *>> &op)

DECLARE_CONTRACT(1, 1);
DECLARE_CONTRACT(2, 1);
DECLARE_CONTRACT(3, 1);
DECLARE_CONTRACT(4, 1);
DECLARE_CONTRACT(1, 2);
DECLARE_CONTRACT(2, 2);
DECLARE_CONTRACT(3, 2);
DECLARE_CONTRACT(4, 2);
DECLARE_CONTRACT(1, 3);
DECLARE_CONTRACT(2, 3);
DECLARE_CONTRACT(3, 3);
DECLARE_CONTRACT(4, 3);
DECLARE_CONTRACT(1, 4);
DECLARE_CONTRACT(2, 4);
DECLARE_CONTRACT(3, 4);
DECLARE_CONTRACT(4, 4);

#undef DECLARE_CONTRACT

}  // end of namespace tenes

#endif  // SRC_PEPS_CONTRACT_IMPL_HPP_
