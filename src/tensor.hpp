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

#ifndef TENES_SRC_TENSOR_HPP_
#define TENES_SRC_TENSOR_HPP_

#include <cstddef>
#include <vector>
#include <complex>

#include "mptensor/tensor.hpp"  // IWYU pragma: export

namespace tenes {

#ifdef _NO_MPI
template <class T>
using mptensor_matrix_type = mptensor::lapack::Matrix<T>;
#else
template <class T>
using mptensor_matrix_type = mptensor::scalapack::Matrix<T>;
#endif  // USE_MPI

template <class T>
using mptensor_tensor_type = mptensor::Tensor<mptensor_matrix_type, T>;

template <class T>
using small_tensor = mptensor::Tensor<mptensor::lapack::Matrix, T>;

using real_tensor = mptensor_tensor_type<double>;
using complex_tensor = mptensor_tensor_type<std::complex<double>>;

template <class tensor>
tensor resize_tensor(tensor const& src, mptensor::Shape target_shape);

template <class T>
void eigen(small_tensor<T> const& A, std::vector<std::complex<double>>& eigvals,
           std::vector<std::complex<double>>& eigvecs_last, int nev);

template <class T>
small_tensor<T> identity(std::size_t k, T v) {
  small_tensor<T> ret{mptensor::Shape(k, k)};
  for (std::size_t i = 0; i < k; ++i) {
    ret.set_value({k, k}, v);
  }
  return ret;
}

}  // end of namespace tenes

#endif  // TENES_SRC_TENSOR_HPP_
