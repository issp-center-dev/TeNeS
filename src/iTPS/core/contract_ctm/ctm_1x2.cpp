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

#include "../contract_ctm.hpp"
#include "instantiate.hpp"

namespace tenes {
namespace itps {
namespace core {
using mptensor::Axes;

template <class tensor>
typename tensor::value_type Contract_1x2(
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
  // 1_2.dat
  ////////////////////////////////////////////////////////////
  // (*(op[0][0])*(*(Tn[0][0])*((*(eTt[0])*(*(C[0])**(eTl[0])))*(conj(*(Tn[0][0]))*((*(C[3])**(eTb[0]))*(*(eTt[1])*(*(Tn[0][1])*((conj(*(Tn[0][1]))**(op[0][1]))*(*(eTb[1])*(*(C[1])*(*(C[2])**(eTr[0]))))))))))))
  // cpu_cost= 2.12419e+07  memory= 405776
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      *(op[0][0]),
      tensordot(
          *(Tn[0][0]),
          tensordot(
              tensordot(*(eTt[0]),
                        tensordot(*(C[0]), *(eTl[0]), Axes(0), Axes(1)),
                        Axes(0), Axes(0)),
              tensordot(
                  conj(*(Tn[0][0])),
                  tensordot(
                      tensordot(*(C[3]), *(eTb[0]), Axes(0), Axes(1)),
                      tensordot(
                          *(eTt[1]),
                          tensordot(
                              *(Tn[0][1]),
                              tensordot(
                                  tensordot(conj(*(Tn[0][1])), *(op[0][1]),
                                            Axes(4), Axes(1)),
                                  tensordot(
                                      *(eTb[1]),
                                      tensordot(*(C[1]),
                                                tensordot(*(C[2]), *(eTr[0]),
                                                          Axes(0), Axes(1)),
                                                Axes(1), Axes(1)),
                                      Axes(0), Axes(1)),
                                  Axes(2, 3), Axes(5, 2)),
                              Axes(2, 3, 4), Axes(6, 4, 2)),
                          Axes(1, 2, 3), Axes(5, 1, 3)),
                      Axes(1), Axes(3)),
                  Axes(2, 3), Axes(5, 2)),
              Axes(0, 2, 3, 5), Axes(5, 1, 3, 0)),
          Axes(0, 1, 2, 3), Axes(1, 0, 4, 3)),
      Axes(0, 1), Axes(0, 1));
}

//! @cond
INSTANTIATE_CONTRACT(real_tensor, 1, 2);
INSTANTIATE_CONTRACT(complex_tensor, 1, 2);
//! @endcond

}  // namespace core
}  // namespace itps
}  // namespace tenes
