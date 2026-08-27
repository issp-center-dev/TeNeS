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

/*! @file
 *  @brief MPI wrapper: typed bcast / allreduce helpers, and stub
 *         definitions of the used MPI symbols so _NO_MPI builds compile
 *         without an MPI library.
 */

#ifndef TENES_SRC_MPI_HPP_
#define TENES_SRC_MPI_HPP_

#include <type_traits>
#include <cstddef>        // for size_t
#include <array>          // for array
#include <complex>        // for complex
#include <string>         // for string
#include <vector>         // for vector
#include "exception.hpp"  // for unimplemented_error

#ifdef _NO_MPI

// Mocks

using MPI_Comm = int;
using MPI_Datatype = int;

constexpr MPI_Comm MPI_COMM_WORLD = 0;
constexpr MPI_Comm MPI_COMM_NULL = 0;
constexpr MPI_Datatype MPI_BYTE = 0;
constexpr MPI_Datatype MPI_INT = 0;
constexpr MPI_Datatype MPI_DOUBLE = 0;

// Do nothing
int MPI_Init(int *, char ***);
int MPI_Barrier(MPI_Comm);
int MPI_Finalize();
int MPI_Bcast(void *, int, MPI_Datatype, int, MPI_Comm);
int MPI_Send(const void *, int, MPI_Datatype, int, int, MPI_Comm);
int MPI_Comm_disconnect(MPI_Comm *);

int MPI_Comm_size(MPI_Comm, int *);   // return 1 as size
int MPI_Comm_rank(MPI_Comm, int *);   // return 0 as rank
int MPI_Comm_get_parent(MPI_Comm *);  // return MPI_COMM_NULL

#else

#include <mpi.h>  // IWYU pragma: export

#endif  // _NO_MPI

namespace tenes {

//! The MPI datatype tag of T (bool travels as MPI_INT); throws
//! tenes::unimplemented_error for unsupported types.
template <class T>
MPI_Datatype get_MPI_Datatype() {
  if constexpr (std::is_same_v<T, int>) {
    return MPI_INT;
  } else if constexpr (std::is_same_v<T, bool>) {
    return MPI_INT;
  } else if constexpr (std::is_same_v<T, double>) {
    return MPI_DOUBLE;
  } else {
    throw tenes::unimplemented_error("");
  }
}

/*! @brief Broadcast a value from root to every rank (no-op without MPI).
 *
 *  The overloads cover scalars, strings, complex values, and vectors of
 *  those; vectors are resized on the receiving ranks. All return the MPI
 *  error code (0 on success).
 */
template <class T>
int bcast(T &val, int root, MPI_Comm comm) {
  int ret = 0;
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  ret = MPI_Bcast(&val, 1, datatype, root, comm);
#endif
  return ret;
}

//! bcast() overload for bool.
int bcast(bool &val, int root, MPI_Comm comm);
//! bcast() overload for strings.
int bcast(std::string &val, int root, MPI_Comm comm);
//! bcast() overload for vectors of strings.
int bcast(std::vector<std::string> &val, int root, MPI_Comm comm);

//! bcast() overload for complex scalars.
template <class T>
int bcast(std::complex<T> &val, int root, MPI_Comm comm) {
  int ret = 0;
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  std::array<T, 2> reim = {val.real(), val.imag()};
  ret = MPI_Bcast(&(reim[0]), 2, datatype, root, comm);
  val = std::complex<T>(reim[0], reim[1]);
#endif
  return ret;
}

//! bcast() overload for vectors; receiving ranks are resized.
template <class T>
int bcast(std::vector<T> &val, int root, MPI_Comm comm) {
  int ret = 0;
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  int rank;
  MPI_Comm_rank(comm, &rank);
  int sz = val.size();
  ret = MPI_Bcast(&sz, 1, MPI_INT, root, comm);
  if (ret != 0) {
    return ret;
  }
  if (rank != root) {
    val.resize(sz);
  }
  ret = MPI_Bcast(val.data(), sz, datatype, root, comm);
#endif
  return ret;
}

//! bcast() overload for vectors of complex values.
template <class T>
int bcast(std::vector<std::complex<T>> &val, int root, MPI_Comm comm) {
  int ret = 0;
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  int rank;
  MPI_Comm_rank(comm, &rank);
  int sz = val.size();
  ret = MPI_Bcast(&sz, 1, MPI_INT, root, comm);
  if (ret != 0) {
    return ret;
  }
  if (rank != root) {
    val.resize(sz);
  }
  std::vector<T> reim(2 * sz);
  for (int i = 0; i < sz; ++i) {
    reim[2 * i] = val[i].real();
    reim[2 * i + 1] = val[i].imag();
  }
  ret = MPI_Bcast(reim.data(), 2 * sz, datatype, root, comm);
  for (int i = 0; i < sz; ++i) {
    val[i] = std::complex<T>(reim[2 * i], reim[2 * i + 1]);
  }
#endif
  return ret;
}

/*! @brief In-place sum of val over all ranks (no-op without MPI).
 *
 *  The overloads cover scalars and vectors, real and complex (the plain
 *  complex scalar overload is unimplemented and throws). All return the
 *  MPI error code (0 on success).
 */
template <class T>
int allreduce_sum(T &val, MPI_Comm comm) {
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  T recv;
  int ret = MPI_Allreduce(&val, &recv, 1, datatype, MPI_SUM, comm);
  if (ret != 0) {
    return ret;
  }
  val = recv;
#endif
  return 0;
}

//! allreduce_sum() overload for vectors (elementwise sum).
template <class T>
int allreduce_sum(std::vector<T> &val, MPI_Comm comm) {
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  int N = val.size();
  std::vector<T> recv(val);
  int ret = MPI_Allreduce(val.data(), recv.data(), N, datatype, MPI_SUM, comm);
  if (ret != 0) {
    return ret;
  }
  val.assign(recv.begin(), recv.end());
#endif
  return 0;
}

//! Unimplemented complex-scalar overload: always throws.
template <class T>
int allreduce_sum(std::complex<T> /* &val */, MPI_Comm /* comm */) {
  throw tenes::unimplemented_error("allreduce for complex is not implemented");
}

//! allreduce_sum() overload for vectors of complex values.
template <class T>
int allreduce_sum(std::vector<std::complex<T>> &val, MPI_Comm comm) {
  int ret = 0;
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  const int sz = static_cast<int>(val.size());
  std::vector<T> reim(2 * sz);
  for (int i = 0; i < sz; ++i) {
    reim[2 * i] = val[i].real();
    reim[2 * i + 1] = val[i].imag();
  }
  std::vector<T> recv(reim);
  ret =
      MPI_Allreduce(reim.data(), recv.data(), 2 * sz, datatype, MPI_SUM, comm);
  if (ret != 0) {
    return ret;
  }
  for (int i = 0; i < sz; ++i) {
    val[i] = std::complex<T>(recv[2 * i], recv[2 * i + 1]);
  }
#endif
  return ret;
}

//! In-place elementwise maximum of val over all ranks (no-op without
//! MPI); returns the MPI error code.
template <class T>
int allreduce_max(std::vector<T> &val, MPI_Comm comm) {
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  int N = val.size();
  std::vector<T> recv(val);
  int ret = MPI_Allreduce(val.data(), recv.data(), N, datatype, MPI_MAX, comm);
  if (ret != 0) {
    return ret;
  }
  val.assign(recv.begin(), recv.end());
#endif
  return 0;
}

//! In-place elementwise minimum of val over all ranks (no-op without
//! MPI); returns the MPI error code.
template <class T>
int allreduce_min(std::vector<T> &val, MPI_Comm comm) {
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  int N = val.size();
  std::vector<T> recv(val);
  int ret = MPI_Allreduce(val.data(), recv.data(), N, datatype, MPI_MIN, comm);
  if (ret != 0) {
    return ret;
  }
  val.assign(recv.begin(), recv.end());
#endif
  return 0;
}

}  // namespace tenes

#endif  // TENES_SRC_MPI_HPP_
