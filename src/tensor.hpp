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

#ifndef TENSOR_HPP
#define TENSOR_HPP

#include <sstream>

#include <mptensor/tensor.hpp>

#include "exception.hpp"

namespace tenes {

#ifdef _NO_MPI
template <class T> using mptensor_matrix_type = mptensor::lapack::Matrix<T>;
#else
template <class T> using mptensor_matrix_type = mptensor::scalapack::Matrix<T>;
#endif // USE_MPI

template <class T> using mptensor_tensor_type = mptensor::Tensor<mptensor_matrix_type, T>;
using real_tensor = mptensor_tensor_type<double>;
using complex_tensor = mptensor_tensor_type<std::complex<double>>;

template <class T>
mptensor_tensor_type<T> resize_tensor(mptensor_tensor_type<T> const& src, mptensor::Shape target_shape){
  mptensor::Shape shape = src.shape();
  const size_t ndim = shape.size();
  if(target_shape.size() != ndim){
    std::stringstream ss;
    ss << "dimension mismatch in resize_tensor: source = " << ndim << ", target = " << target_shape.size();
    tenes::logic_error(ss.str());
  }
  mptensor::Shape zero = shape;

  bool to_extend = false;
  bool to_shrink = false;

  for(size_t i=0; i<ndim; ++i){
    if(shape[i] < target_shape[i]){
      to_extend = true;
      shape[i] = target_shape[i];
    }else if(shape[i] > target_shape[i]){
      to_shrink = true;
    }
    zero[i] = 0;
  }
  mptensor_tensor_type<T> A;
  if(to_extend){
    A = mptensor::extend(src, shape);
  }else{
    A = src;
  }
  if(to_shrink){
    return mptensor::slice(A, zero, target_shape);
  }else{
    return A;
  }
}

} // end of namespace tenes

#endif // TENSOR_HPP
