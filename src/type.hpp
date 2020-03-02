#ifndef TYPES_HPP
#define TYPES_HPP

#include <mptensor/tensor.hpp>

namespace tenes {

#ifdef _NO_MPI
template <class T> using mptensor_matrix_type = mptensor::lapack::Matrix<T>;
#else
template <class T> using mptensor_matrix_type = mptensor::scalapack::Matrix<T>;
#endif // USE_MPI

} // end of namespace tenes

#endif // TYPES_HPP
