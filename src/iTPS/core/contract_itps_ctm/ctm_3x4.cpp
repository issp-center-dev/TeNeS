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

#include "../contract_itps_ctm.hpp"
#include "instantiate.hpp"

namespace tenes {
namespace itps {
namespace core {
using mptensor::Axes;

template <class tensor>
typename tensor::value_type Contract_3x4_iTPS_CTM(
    const std::vector<const tensor *> &C,
    const std::vector<const tensor *> &eTt,
    const std::vector<const tensor *> &eTr,
    const std::vector<const tensor *> &eTb,
    const std::vector<const tensor *> &eTl,
    const std::vector<std::vector<const tensor *>> &Tn,
    const std::vector<std::vector<const tensor *>> &op) {
#ifndef NDEBUG
  const size_t nrow = Tn.size();
  const size_t ncol = Tn[0].size();
  for (const auto &r : Tn) assert(r.size() == ncol);
  assert(C.size() == 4);
  assert(eTt.size() == ncol);
  assert(eTr.size() == nrow);
  assert(eTb.size() == ncol);
  assert(eTl.size() == nrow);
  assert(op.size() == nrow);
  for (const auto &r : op) assert(r.size() == ncol);
#endif
  ////////////////////////////////////////////////////////////
  // 3_4.dat
  ////////////////////////////////////////////////////////////
  // (*(op[0][0])*(*(Tn[0][0])*((*(eTt[0])*(*(C[0])**(eTl[0])))*(conj(*(Tn[0][0]))*(*(eTl[1])*(*(Tn[1][0])*((conj(*(Tn[1][0]))**(op[1][0]))*((*(Tn[2][0])*((conj(*(Tn[2][0]))**(op[2][0]))*(*(eTb[0])*(*(C[3])**(eTl[2])))))*(*(eTt[1])*(*(Tn[0][1])*((conj(*(Tn[0][1]))**(op[0][1]))*(*(Tn[1][1])*((conj(*(Tn[1][1]))**(op[1][1]))*(*(Tn[2][1])*((conj(*(Tn[2][1]))**(op[2][1]))*(*(eTb[1])*(*(eTt[2])*(*(Tn[0][2])*((conj(*(Tn[0][2]))**(op[0][2]))*(*(Tn[1][2])*((conj(*(Tn[1][2]))**(op[1][2]))*(*(Tn[2][2])*((conj(*(Tn[2][2]))**(op[2][2]))*(*(eTb[2])*((*(Tn[0][3])*((conj(*(Tn[0][3]))**(op[0][3]))*(*(eTt[3])*(*(C[1])**(eTr[0])))))*(*(Tn[1][3])*((conj(*(Tn[1][3]))**(op[1][3]))*(*(eTr[1])*(*(Tn[2][3])*((conj(*(Tn[2][3]))**(op[2][3]))*(*(eTr[2])*(*(C[2])**(eTb[3])))))))))))))))))))))))))))))))))
  // cpu_cost= 1.48356e+10  memory= 8.405e+07
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(*(op[0][0]),
               tensordot(
                   *(Tn[0][0]),
                   tensordot(
                       tensordot(
                           *(eTt[0]),
                           tensordot(*(C[0]), *(eTl[0]), Axes(0), Axes(1)),
                           Axes(0), Axes(0)),
                       tensordot(
                           conj(*(Tn[0][0])),
                           tensordot(
                               *(eTl[1]),
                               tensordot(
                                   *(Tn[1][0]),
                                   tensordot(
                                       tensordot(conj(*(Tn[1][0])), *(op[1][0]),
                                                 Axes(4), Axes(1)),
                                       tensordot(
                                           tensordot(
                                               *(Tn[2][0]),
                                               tensordot(
                                                   tensordot(conj(*(Tn[2][0])),
                                                             *(op[2][0]),
                                                             Axes(4), Axes(1)),
                                                   tensordot(
                                                       *(eTb[0]),
                                                       tensordot(*(C[3]), *(eTl[2]),
                                                                 Axes(1), Axes(0)),
                                                       Axes(1), Axes(0)),
                                                   Axes(0, 3), Axes(5, 2)),
                                               Axes(0, 3, 4), Axes(6, 4, 2)),
                                           tensordot(
                                               *(eTt[1]),
                                               tensordot(
                                                   *(Tn[0][1]),
                                                   tensordot(
                                                       tensordot(conj(
                                                                     *(Tn[0][1])),
                                                                 *(op[0][1]),
                                                                 Axes(4), Axes(1)),
                                                       tensordot(
                                                           *(Tn[1][1]),
                                                           tensordot(
                                                               tensordot(
                                                                   conj(
                                                                       *(Tn[1][1])),
                                                                   *(op[1][1]),
                                                                   Axes(4), Axes(1)),
                                                               tensordot(
                                                                   *(Tn[2][1]),
                                                                   tensordot(
                                                                       tensordot(
                                                                           conj(
                                                                               *(Tn[2][1])),
                                                                           *(op[2]
                                                                               [1]),
                                                                           Axes(4), Axes(1)),
                                                                       tensordot(*(eTb[1]),
                                                                                 tensordot(*(eTt[2]),
                                                                                           tensordot(*(Tn[0]
                                                                                                         [2]),
                                                                                                     tensordot(tensordot(conj(
                                                                                                                             *(Tn[0][2])),
                                                                                                                         *(op[0]
                                                                                                                             [2]),
                                                                                                                         Axes(4), Axes(1)),
                                                                                                               tensordot(*(Tn[1]
                                                                                                                             [2]),
                                                                                                                         tensordot(tensordot(conj(
                                                                                                                                                 *(Tn
                                                                                                                                                       [1][2])),
                                                                                                                                             *(op[1]
                                                                                                                                                 [2]),
                                                                                                                                             Axes(
                                                                                                                                                 4),
                                                                                                                                             Axes(
                                                                                                                                                 1)),
                                                                                                                                   tensordot(
                                                                                                                                       *(Tn[2]
                                                                                                                                           [2]),
                                                                                                                                       tensordot(
                                                                                                                                           tensordot(
                                                                                                                                               conj(*(Tn[2][2])), *(op[2][2]),
                                                                                                                                               Axes(
                                                                                                                                                   4),
                                                                                                                                               Axes(
                                                                                                                                                   1)),
                                                                                                                                           tensordot(*(eTb[2]),
                                                                                                                                                     tensordot(tensordot(*(Tn[0]
                                                                                                                                                                             [3]),
                                                                                                                                                                         tensordot(tensordot(conj(*(Tn[0][3])), *(op[0][3]), Axes(4), Axes(1)), tensordot(*(eTt[3]), tensordot(*(C[1]), *(eTr[0]), Axes(1), Axes(0)), Axes(1), Axes(0)), Axes(1, 2), Axes(2, 5)),
                                                                                                                                                                         Axes(1,
                                                                                                                                                                              2,
                                                                                                                                                                              4),
                                                                                                                                                                         Axes(4,
                                                                                                                                                                              6,
                                                                                                                                                                              2)),
                                                                                                                                                               tensordot(*(Tn[1][3]), tensordot(tensordot(conj(*(Tn[1][3])), *(op[1][3]), Axes(4), Axes(1)), tensordot(*(eTr[1]), tensordot(*(Tn[2][3]), tensordot(tensordot(conj(*(Tn[2][3])), *(op[2][3]), Axes(4), Axes(1)), tensordot(*(eTr[2]), tensordot(*(C[2]), *(eTb[3]), Axes(1), Axes(0)), Axes(1), Axes(0)), Axes(2, 3), Axes(2, 5)), Axes(2, 3, 4), Axes(4, 6, 2)), Axes(1), Axes(4)), Axes(2, 3), Axes(2, 6)), Axes(2, 3, 4), Axes(4, 6, 2)),
                                                                                                                                                               Axes(1, 3, 5),
                                                                                                                                                               Axes(1, 3, 4)),
                                                                                                                                                     Axes(
                                                                                                                                                         0),
                                                                                                                                                     Axes(
                                                                                                                                                         7)),
                                                                                                                                           Axes(
                                                                                                                                               2,
                                                                                                                                               3),
                                                                                                                                           Axes(
                                                                                                                                               9,
                                                                                                                                               2)),
                                                                                                                                       Axes(2,
                                                                                                                                            3,
                                                                                                                                            4),
                                                                                                                                       Axes(
                                                                                                                                           10,
                                                                                                                                           4,
                                                                                                                                           2)),
                                                                                                                                   Axes(
                                                                                                                                       2,
                                                                                                                                       3),
                                                                                                                                   Axes(
                                                                                                                                       9,
                                                                                                                                       3)),
                                                                                                                         Axes(
                                                                                                                             2,
                                                                                                                             3,
                                                                                                                             4),
                                                                                                                         Axes(
                                                                                                                             10,
                                                                                                                             4,
                                                                                                                             2)),
                                                                                                               Axes(
                                                                                                                   2,
                                                                                                                   3),
                                                                                                               Axes(
                                                                                                                   8,
                                                                                                                   3)),
                                                                                                     Axes(
                                                                                                         2,
                                                                                                         3,
                                                                                                         4),
                                                                                                     Axes(
                                                                                                         9,
                                                                                                         4,
                                                                                                         2)),
                                                                                           Axes(
                                                                                               1,
                                                                                               2,
                                                                                               3),
                                                                                           Axes(
                                                                                               9,
                                                                                               1,
                                                                                               3)),
                                                                                 Axes(
                                                                                     0),
                                                                                 Axes(
                                                                                     7)),
                                                                       Axes(2,
                                                                            3),
                                                                       Axes(9,
                                                                            2)),
                                                                   Axes(2, 3,
                                                                        4),
                                                                   Axes(10, 4,
                                                                        2)),
                                                               Axes(2, 3),
                                                               Axes(9, 3)),
                                                           Axes(2, 3, 4),
                                                           Axes(10, 4, 2)),
                                                       Axes(2, 3), Axes(9, 3)),
                                                   Axes(2, 3, 4),
                                                   Axes(10, 4, 2)),
                                               Axes(1, 2, 3), Axes(9, 1, 3)),
                                           Axes(1, 3, 4), Axes(5, 6, 7)),
                                       Axes(2, 3), Axes(7, 1)),
                                   Axes(2, 3, 4), Axes(8, 3, 2)),
                               Axes(0, 2, 3), Axes(4, 0, 2)),
                           Axes(2, 3), Axes(5, 2)),
                       Axes(0, 2, 3, 5), Axes(5, 1, 3, 0)),
                   Axes(0, 1, 2, 3), Axes(1, 0, 4, 3)),
               Axes(0, 1), Axes(0, 1));
}

//! @cond
INSTANTIATE_CONTRACT_ITPS_CTM(real_tensor, 3, 4);
INSTANTIATE_CONTRACT_ITPS_CTM(complex_tensor, 3, 4);
//! @endcond

}  // namespace core
}  // namespace itps
}  // namespace tenes
