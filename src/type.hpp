#ifndef TYPES_HPP
#define TYPES_HPP

#include <mptensor/complex.hpp>
#include <mptensor/tensor.hpp>

#ifdef _NO_MPI
template <class T>
using mptensor_matrix_type = mptensor::lapack::Matrix<T>;
#else
template <class T>
using mptensor_matrix_type = mptensor::scalapack::Matrix<T>;
#endif // USE_MPI

#endif // TYPES_HPP
