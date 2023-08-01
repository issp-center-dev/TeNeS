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

/*
 *
 Basic routines independent on unit cell structures.
 Using mptensor libraries
 (Test version)
 2015 Dec.  Tsuyoshi Okubo
*/

#include "../contract_density_ctm.hpp"
#include "../../../tensor.hpp"

namespace tenes {
namespace itps {
namespace core {

using mptensor::Axes;

template <class tensor>
void StartCorrelation_density_CTM(tensor &A, const tensor &C1, const tensor &C4,
                                  const tensor &eT1, const tensor &eT3,
                                  const tensor &eT4, const tensor &Tn1,
                                  const tensor &op) {
  ////////////////////////////////////////////////////////////
  // ((C1*eT1)*((Tn1*op)*(eT3*(C4*eT4))))
  // cpu_cost= 144384  memory= 6560
  // final_bond_order (eT1_r, eT3_r, Tn1_r)
  ////////////////////////////////////////////////////////////
  A = transpose(
      tensordot(tensordot(C1, eT1, Axes(1), Axes(0)),
                tensordot(tensordot(Tn1, op, Axes(4, 5), Axes(0, 1)),
                          tensordot(eT3, tensordot(C4, eT4, Axes(1), Axes(0)),
                                    Axes(1), Axes(0)),
                          Axes(0, 3), Axes(3, 1)),
                Axes(0, 2), Axes(3, 0)),
      Axes(0, 2, 1));
}

template <class tensor>
void Transfer_density_CTM(tensor &A, const tensor &eT1, const tensor &eT3,
                          const tensor &Tn1) {
  ////////////////////////////////////////////////////////////
  // ./transfer_den.dat
  ////////////////////////////////////////////////////////////
  // (eT1*(Tn1*(A*eT3)))
  // cpu_cost= 196608  memory= 9472
  // final_bond_order (eT1_r, eT3_r, Tn1_r)
  ////////////////////////////////////////////////////////////
  A = transpose(tensordot(eT1,
                          tensordot(contract(Tn1, Axes(4), Axes(5)),
                                    tensordot(A, eT3, Axes(1), Axes(1)),
                                    Axes(0, 3), Axes(1, 3)),
                          Axes(0, 2), Axes(2, 0)),
                Axes(0, 2, 1));
}

template <class tensor>
typename tensor::value_type FinishCorrelation_density_CTM(
    const tensor &A, const tensor &C2, const tensor &C3, const tensor &eT1,
    const tensor &eT2, const tensor &eT3, const tensor &Tn1, const tensor &op) {
  ////////////////////////////////////////////////////////////
  // ./finishcorrelation_den.dat
  ////////////////////////////////////////////////////////////
  // (op*(Tn1*((A*eT1)*(eT3*(C2*(C3*eT2))))))
  // cpu_cost= 230404  memory= 11268
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op,
      tensordot(
          Tn1,
          tensordot(
              tensordot(A, eT1, Axes(0), Axes(0)),
              tensordot(eT3,
                        tensordot(C2, tensordot(C3, eT2, Axes(0), Axes(1)),
                                  Axes(1), Axes(1)),
                        Axes(0), Axes(1)),
              Axes(0, 2), Axes(0, 2)),
          Axes(0, 1, 2, 3), Axes(0, 1, 3, 2)),
      Axes(0, 1), Axes(0, 1));
}

// start of explicit instantiation

template void StartCorrelation_density_CTM(
    real_tensor &A, const real_tensor &C1, const real_tensor &C4,
    const real_tensor &eT1, const real_tensor &eT3, const real_tensor &eT4,
    const real_tensor &Tn1, const real_tensor &op);
template void StartCorrelation_density_CTM(
    complex_tensor &A, const complex_tensor &C1, const complex_tensor &C4,
    const complex_tensor &eT1, const complex_tensor &eT3,
    const complex_tensor &eT4, const complex_tensor &Tn1,
    const complex_tensor &op);

template void Transfer_density_CTM(real_tensor &A, const real_tensor &eT1,
                                   const real_tensor &eT3,
                                   const real_tensor &Tn1);
template void Transfer_density_CTM(complex_tensor &A, const complex_tensor &eT1,
                                   const complex_tensor &eT3,
                                   const complex_tensor &Tn1);

template typename real_tensor::value_type FinishCorrelation_density_CTM(
    const real_tensor &A, const real_tensor &C2, const real_tensor &C3,
    const real_tensor &eT1, const real_tensor &eT2, const real_tensor &eT3,
    const real_tensor &Tn1, const real_tensor &op);
template typename complex_tensor::value_type FinishCorrelation_density_CTM(
    const complex_tensor &A, const complex_tensor &C2, const complex_tensor &C3,
    const complex_tensor &eT1, const complex_tensor &eT2,
    const complex_tensor &eT3, const complex_tensor &Tn1,
    const complex_tensor &op);

// end of explicit instantiation

}  // end of namespace core
}  // namespace itps
}  // namespace tenes
