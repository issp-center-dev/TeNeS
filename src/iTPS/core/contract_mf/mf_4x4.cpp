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

#include "../contract_mf.hpp"
#include "instantiate.hpp"

namespace tenes {
using mptensor::Axes;

template <class tensor>
typename tensor::value_type Contract_MF_4x4(
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
  // 4_4_mf.dat
  ////////////////////////////////////////////////////////////
  // (*(op[0][0])*(*(Tn[0][0])*(conj(*(Tn[0][0]))*((*(Tn[0][1])*(conj(*(Tn[0][1]))**(op[0][1])))*((*(Tn[1][0])*(conj(*(Tn[1][0]))**(op[1][0])))*(*(Tn[1][1])*((conj(*(Tn[1][1]))**(op[1][1]))*((*(Tn[1][2])*((conj(*(Tn[1][2]))**(op[1][2]))*((*(Tn[0][2])*(conj(*(Tn[0][2]))**(op[0][2])))*((*(Tn[0][3])*(conj(*(Tn[0][3]))**(op[0][3])))*(*(Tn[1][3])*(conj(*(Tn[1][3]))**(op[1][3])))))))*((*(Tn[2][1])*((conj(*(Tn[2][1]))**(op[2][1]))*((*(Tn[2][0])*(conj(*(Tn[2][0]))**(op[2][0])))*((*(Tn[3][0])*(conj(*(Tn[3][0]))**(op[3][0])))*(*(Tn[3][1])*(conj(*(Tn[3][1]))**(op[3][1])))))))*(*(Tn[2][2])*((conj(*(Tn[2][2]))**(op[2][2]))*((*(Tn[2][3])*(conj(*(Tn[2][3]))**(op[2][3])))*((*(Tn[3][3])*(conj(*(Tn[3][3]))**(op[3][3])))*(*(Tn[3][2])*(conj(*(Tn[3][2]))**(op[3][2]))))))))))))))))
  // cpu_cost= 2.10667e+11  memory= 1.8527e+08
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
                            Axes(1, 4), Axes(1, 4)),
                  tensordot(
                      tensordot(*(Tn[1][0]),
                                tensordot(conj(*(Tn[1][0])), *(op[1][0]),
                                          Axes(4), Axes(1)),
                                Axes(0, 4), Axes(0, 4)),
                      tensordot(
                          *(Tn[1][1]),
                          tensordot(
                              tensordot(conj(*(Tn[1][1])), *(op[1][1]), Axes(4),
                                        Axes(1)),
                              tensordot(
                                  tensordot(
                                      *(Tn[1][2]),
                                      tensordot(
                                          tensordot(conj(*(Tn[1][2])),
                                                    *(op[1][2]), Axes(4),
                                                    Axes(1)),
                                          tensordot(
                                              tensordot(
                                                  *(Tn[0][2]),
                                                  tensordot(conj(*(Tn[0][2])),
                                                            *(op[0][2]),
                                                            Axes(4), Axes(1)),
                                                  Axes(1, 4), Axes(1, 4)),
                                              tensordot(
                                                  tensordot(
                                                      *(Tn[0][3]),
                                                      tensordot(
                                                          conj(*(Tn[0][3])),
                                                          *(op[0][3]), Axes(4),
                                                          Axes(1)),
                                                      Axes(1, 2, 4),
                                                      Axes(1, 2, 4)),
                                                  tensordot(
                                                      *(Tn[1][3]),
                                                      tensordot(
                                                          conj(*(Tn[1][3])),
                                                          *(op[1][3]), Axes(4),
                                                          Axes(1)),
                                                      Axes(2, 4), Axes(2, 4)),
                                                  Axes(1, 3), Axes(1, 4)),
                                              Axes(1, 4), Axes(0, 1)),
                                          Axes(1, 2), Axes(3, 6)),
                                      Axes(1, 2, 4), Axes(4, 6, 2)),
                                  tensordot(
                                      tensordot(
                                          *(Tn[2][1]),
                                          tensordot(
                                              tensordot(conj(*(Tn[2][1])),
                                                        *(op[2][1]), Axes(4),
                                                        Axes(1)),
                                              tensordot(
                                                  tensordot(
                                                      *(Tn[2][0]),
                                                      tensordot(
                                                          conj(*(Tn[2][0])),
                                                          *(op[2][0]), Axes(4),
                                                          Axes(1)),
                                                      Axes(0, 4), Axes(0, 4)),
                                                  tensordot(
                                                      tensordot(
                                                          *(Tn[3][0]),
                                                          tensordot(
                                                              conj(*(Tn[3][0])),
                                                              *(op[3][0]),
                                                              Axes(4), Axes(1)),
                                                          Axes(0, 3, 4),
                                                          Axes(0, 3, 4)),
                                                      tensordot(
                                                          *(Tn[3][1]),
                                                          tensordot(
                                                              conj(*(Tn[3][1])),
                                                              *(op[3][1]),
                                                              Axes(4), Axes(1)),
                                                          Axes(3, 4),
                                                          Axes(3, 4)),
                                                      Axes(1, 3), Axes(0, 3)),
                                                  Axes(2, 5), Axes(0, 1)),
                                              Axes(0, 3), Axes(3, 6)),
                                          Axes(0, 3, 4), Axes(4, 6, 2)),
                                      tensordot(
                                          *(Tn[2][2]),
                                          tensordot(
                                              tensordot(conj(*(Tn[2][2])),
                                                        *(op[2][2]), Axes(4),
                                                        Axes(1)),
                                              tensordot(
                                                  tensordot(
                                                      *(Tn[2][3]),
                                                      tensordot(
                                                          conj(*(Tn[2][3])),
                                                          *(op[2][3]), Axes(4),
                                                          Axes(1)),
                                                      Axes(2, 4), Axes(2, 4)),
                                                  tensordot(
                                                      tensordot(
                                                          *(Tn[3][3]),
                                                          tensordot(
                                                              conj(*(Tn[3][3])),
                                                              *(op[3][3]),
                                                              Axes(4), Axes(1)),
                                                          Axes(2, 3, 4),
                                                          Axes(2, 3, 4)),
                                                      tensordot(
                                                          *(Tn[3][2]),
                                                          tensordot(
                                                              conj(*(Tn[3][2])),
                                                              *(op[3][2]),
                                                              Axes(4), Axes(1)),
                                                          Axes(3, 4),
                                                          Axes(3, 4)),
                                                      Axes(0, 2), Axes(2, 5)),
                                                  Axes(2, 5), Axes(0, 1)),
                                              Axes(2, 3), Axes(2, 7)),
                                          Axes(2, 3, 4), Axes(3, 7, 2)),
                                      Axes(1, 3, 6, 7), Axes(0, 2, 6, 7)),
                                  Axes(1, 3, 6, 7), Axes(4, 5, 6, 7)),
                              Axes(2, 3), Axes(1, 5)),
                          Axes(2, 3, 4), Axes(3, 6, 2)),
                      Axes(1, 2, 4, 5), Axes(0, 6, 2, 7)),
                  Axes(1, 2, 4, 5), Axes(4, 2, 5, 3)),
              Axes(2, 3), Axes(1, 3)),
          Axes(0, 1, 2, 3), Axes(0, 1, 3, 4)),
      Axes(0, 1), Axes(0, 1));
}

INSTANTIATE_CONTRACT(real_tensor, 4, 4);
INSTANTIATE_CONTRACT(complex_tensor, 4, 4);

}  // end of namespace tenes
