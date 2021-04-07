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

#include "../contract_mf.hpp"
#include "../../../tensor.hpp"

namespace tenes {
namespace itps {
namespace core {

using mptensor::Axes;

/* Lambda tensors should be absorbed before entering */
template <class tensor>
void StartCorrelation_MF(tensor &A, const tensor &Tn, const tensor &op,
                         size_t direction) {
  Axes axes;
  for (size_t dir = 0; dir < direction; ++dir) {
    axes.push(dir);
  }
  for (size_t dir = direction + 1; dir < 4; ++dir) {
    axes.push(dir);
  }
  axes.push(4);

  A = tensordot(Tn, tensordot(conj(Tn), op, Axes(4), Axes(1)), axes, axes);
}

/* Lambda tensors should be absorbed before entering */
template <class tensor>
void Transfer_MF(tensor &A, const tensor &Tn, size_t direction) {
  size_t direction2 = (direction + 2) % 4;

  Axes axes;
  for (size_t d = 0; d < 4; ++d) {
    if (d == direction || d == direction2) continue;
    axes.push(d);
  }
  axes.push(4);

  tensor transfer = tensordot(Tn, conj(Tn), axes, axes);
  size_t leg = direction < direction2 ? 1 : 0;
  A = tensordot(A, transfer, Axes(0, 1), Axes(leg, leg + 2));
}

/* Lambda tensors should be absorbed before entering */
template <class tensor>
typename tensor::value_type FinishCorrelation_MF(const tensor &A,
                                                 const tensor &Tn,
                                                 const tensor &op,
                                                 size_t direction) {
  direction += 2;
  direction %= 4;

  Axes axes;
  for (size_t dir = 0; dir < direction; ++dir) {
    axes.push(dir);
  }
  for (size_t dir = direction + 1; dir < 4; ++dir) {
    axes.push(dir);
  }
  axes.push(4);

  tensor B =
      tensordot(Tn, tensordot(conj(Tn), op, Axes(4), Axes(1)), axes, axes);
  return trace(tensordot(A, B, Axes(0, 1), Axes(0, 1)));
}

template void StartCorrelation_MF(real_tensor &A, const real_tensor &Tn,
                                  const real_tensor &op, size_t direction);
template void StartCorrelation_MF(complex_tensor &A, const complex_tensor &Tn,
                                  const complex_tensor &op, size_t direction);

template void Transfer_MF(real_tensor &A, const real_tensor &Tn,
                          size_t direction);
template void Transfer_MF(complex_tensor &A, const complex_tensor &Tn,
                          size_t direction);

template typename real_tensor::value_type FinishCorrelation_MF(
    const real_tensor &A, const real_tensor &Tn, const real_tensor &op,
    size_t direction);
template typename complex_tensor::value_type FinishCorrelation_MF(
    const complex_tensor &A, const complex_tensor &Tn, const complex_tensor &op,
    size_t direction);

}  // namespace core
}  // namespace itps
}  // namespace tenes
