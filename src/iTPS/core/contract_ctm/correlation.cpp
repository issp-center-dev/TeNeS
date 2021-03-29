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

#include "../contract_ctm.hpp"
#include "../../../tensor.hpp"

namespace tenes {

using mptensor::Axes;

template <class tensor>
void StartCorrelation(tensor &A, const tensor &C1, const tensor &C4,
                      const tensor &eT1, const tensor &eT3, const tensor &eT4,
                      const tensor &Tn1, const tensor &op) {
  ////////////////////////////////////////////////////////////
  // (eT1*(Tn1*((Tn2*op)*(eT3*(C1*(C4*eT4))))))
  // cpu_cost= 7.5525e+06  memory= 192500
  // final_bond_order (e1r, e3r, n1r, n2r)
  ////////////////////////////////////////////////////////////
  A = transpose(
      tensordot(
          eT1,
          tensordot(
              Tn1,
              tensordot(
                  tensordot(conj(Tn1), op, Axes(4), Axes(1)),
                  tensordot(eT3,
                            tensordot(C1, tensordot(C4, eT4, Axes(1), Axes(0)),
                                      Axes(0), Axes(1)),
                            Axes(1), Axes(1)),
                  Axes(0, 3), Axes(5, 2)),
              Axes(0, 3, 4), Axes(6, 4, 2)),
          Axes(0, 2, 3), Axes(5, 0, 2)),
      Axes(0, 3, 1, 2));
}

template <class tensor>
void Transfer(tensor &A, const tensor &eT1, const tensor &eT3,
              const tensor &Tn1) {
  ////////////////////////////////////////////////////////////
  // (eT1*(Tn1*(Tn2*(A*eT3))))
  // cpu_cost= 7.5e+06  memory= 192500
  // final_bond_order (e1r, e3r, n1r, n2r)
  ////////////////////////////////////////////////////////////
  A = transpose(
      tensordot(
          eT1,
          tensordot(Tn1,
                    tensordot(conj(Tn1), tensordot(A, eT3, Axes(1), Axes(1)),
                              Axes(0, 3), Axes(2, 5)),
                    Axes(0, 3, 4), Axes(4, 6, 2)),
          Axes(0, 2, 3), Axes(4, 0, 2)),
      Axes(0, 3, 1, 2));
}

template <class tensor>
typename tensor::value_type FinishCorrelation(
    const tensor &A, const tensor &C2, const tensor &C3, const tensor &eT1,
    const tensor &eT2, const tensor &eT3, const tensor &Tn1, const tensor &op) {
  ////////////////////////////////////////////////////////////
  // (op*(Tn1*((A*eT1)*(Tn2*(eT3*(C2*(C3*eT2)))))))
  // cpu_cost= 7.5525e+06  memory= 252504
  // final_bond_order ()
  ////////////////////////////////////////////////////////////
  return trace(
      op,
      tensordot(
          Tn1,
          tensordot(
              tensordot(A, eT1, Axes(0), Axes(0)),
              tensordot(
                  conj(Tn1),
                  tensordot(eT3,
                            tensordot(C2, tensordot(C3, eT2, Axes(0), Axes(1)),
                                      Axes(1), Axes(1)),
                            Axes(0), Axes(1)),
                  Axes(2, 3), Axes(5, 2)),
              Axes(0, 2, 3, 5), Axes(3, 0, 5, 1)),
          Axes(0, 1, 2, 3), Axes(0, 1, 4, 3)),
      Axes(0, 1), Axes(0, 1));
}

template void StartCorrelation(real_tensor &A, const real_tensor &C1,
                               const real_tensor &C4, const real_tensor &eT1,
                               const real_tensor &eT3, const real_tensor &eT4,
                               const real_tensor &Tn1, const real_tensor &op);
template void StartCorrelation(complex_tensor &A, const complex_tensor &C1,
                               const complex_tensor &C4,
                               const complex_tensor &eT1,
                               const complex_tensor &eT3,
                               const complex_tensor &eT4,
                               const complex_tensor &Tn1,
                               const complex_tensor &op);

template void Transfer(real_tensor &A, const real_tensor &eT1,
                       const real_tensor &eT3, const real_tensor &Tn1);
template void Transfer(complex_tensor &A, const complex_tensor &eT1,
                       const complex_tensor &eT3, const complex_tensor &Tn1);

template typename real_tensor::value_type FinishCorrelation(
    const real_tensor &A, const real_tensor &C2, const real_tensor &C3,
    const real_tensor &eT1, const real_tensor &eT2, const real_tensor &eT3,
    const real_tensor &Tn1, const real_tensor &op);
template typename complex_tensor::value_type FinishCorrelation(
    const complex_tensor &A, const complex_tensor &C2, const complex_tensor &C3,
    const complex_tensor &eT1, const complex_tensor &eT2,
    const complex_tensor &eT3, const complex_tensor &Tn1,
    const complex_tensor &op);

}  // end of namespace tenes
