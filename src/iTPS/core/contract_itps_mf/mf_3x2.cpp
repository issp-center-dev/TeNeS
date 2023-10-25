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

#include "../contract_itps_mf.hpp"
#include "instantiate.hpp"

namespace tenes {
namespace itps {
namespace core {
using mptensor::Axes;

template <class tensor>
typename tensor::value_type Contract_3x2_iTPS_MF(
    const std::vector<std::vector<const tensor *>> &Tn,
    const std::vector<std::vector<const tensor *>> &op) {
#ifndef NDEBUG
  const size_t nrow = Tn.size();
  const size_t ncol = Tn[0].size();
  for (const auto &r : Tn) assert(r.size() == ncol);
  assert(op.size() == nrow);
  for (const auto &r : op) assert(r.size() == ncol);
#endif
  ////////////////////////////////////////////////////////////
  // 3_2_mf.dat
  ////////////////////////////////////////////////////////////
  // (*(op[0][0])*(*(Tn[0][0])*(conj(*(Tn[0][0]))*((*(Tn[0][1])*(conj(*(Tn[0][1]))**(op[0][1])))*((*(Tn[1][0])*(conj(*(Tn[1][0]))**(op[1][0])))*((*(Tn[1][1])*(conj(*(Tn[1][1]))**(op[1][1])))*((*(Tn[2][0])*(conj(*(Tn[2][0]))**(op[2][0])))*(*(Tn[2][1])*(conj(*(Tn[2][1]))**(op[2][1]))))))))))
  // cpu_cost= 7.75947e+07  memory= 921728
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      *(op[0][0]),
      tensordot(
          *(Tn[0][0]),
          tensordot(
              conj(*(Tn[0][0])),
              tensordot(
                  tensordot(*(Tn[0][1]),
                            tensordot(conj(*(Tn[0][1])), *(op[0][1]), Axes(4),
                                      Axes(1)),
                            Axes(1, 2, 4), Axes(1, 2, 4)),
                  tensordot(
                      tensordot(*(Tn[1][0]),
                                tensordot(conj(*(Tn[1][0])), *(op[1][0]),
                                          Axes(4), Axes(1)),
                                Axes(0, 4), Axes(0, 4)),
                      tensordot(
                          tensordot(*(Tn[1][1]),
                                    tensordot(conj(*(Tn[1][1])), *(op[1][1]),
                                              Axes(4), Axes(1)),
                                    Axes(2, 4), Axes(2, 4)),
                          tensordot(tensordot(*(Tn[2][0]),
                                              tensordot(conj(*(Tn[2][0])),
                                                        *(op[2][0]), Axes(4),
                                                        Axes(1)),
                                              Axes(0, 3, 4), Axes(0, 3, 4)),
                                    tensordot(*(Tn[2][1]),
                                              tensordot(conj(*(Tn[2][1])),
                                                        *(op[2][1]), Axes(4),
                                                        Axes(1)),
                                              Axes(2, 3, 4), Axes(2, 3, 4)),
                                    Axes(1, 3), Axes(0, 2)),
                          Axes(2, 5), Axes(2, 3)),
                      Axes(1, 2, 4, 5), Axes(0, 4, 2, 5)),
                  Axes(1, 3), Axes(2, 3)),
              Axes(2, 3), Axes(1, 3)),
          Axes(0, 1, 2, 3), Axes(0, 1, 3, 4)),
      Axes(0, 1), Axes(0, 1));
}

//! @cond
INSTANTIATE_CONTRACT(real_tensor, 3, 2);
INSTANTIATE_CONTRACT(complex_tensor, 3, 2);
//! @endcond

}  // namespace core
}  // namespace itps
}  // namespace tenes
