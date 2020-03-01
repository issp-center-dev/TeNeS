#ifndef TENES_MPI_H
#define TENES_MPI_H

#include <array>
#include <complex>
#include <vector>

#include "exception.hpp"

#ifdef _NO_MPI

using MPI_Comm = int;
using MPI_Datatype = int;

constexpr MPI_Comm MPI_COMM_WORLD = 0;
constexpr MPI_Datatype MPI_BYTE = 0;
constexpr MPI_Datatype MPI_INT = 0;
constexpr MPI_Datatype MPI_DOUBLE = 0;

int MPI_Init(int*, char***);
int MPI_Barrier(MPI_Comm);
int MPI_Finalize();

int MPI_Comm_size(MPI_Comm, int*);
int MPI_Comm_rank(MPI_Comm, int*);

int MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm);

#else

#include <mpi.h>

#endif  //_NO_MPI

namespace tenes{

template <class T>
MPI_Datatype get_MPI_Datatype();
template <> MPI_Datatype get_MPI_Datatype<int>();
template <> MPI_Datatype get_MPI_Datatype<double>();

template <class T>
int bcast(T &val, int root, MPI_Comm comm){
  int ret=0;
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  ret=MPI_Bcast(&val, 1, datatype, root, comm);
#endif
  return ret;
}

template <class T>
int bcast(std::complex<T> &val, int root, MPI_Comm comm){
  int ret=0;
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  std::array<T,2> reim = {val.real(), val.imag()};
  ret=MPI_Bcast(&(reim[0]), 2, datatype, root, comm);
  val = std::complex<T>(reim[0], reim[1]);
#endif
  return ret;
}

template <class T>
int bcast(std::vector<T> &val, int root, MPI_Comm comm){
  int ret=0;
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  int rank;
  MPI_Comm_rank(com, &rank);
  int sz = val.size();
  int ret=MPI_Bcast(&sz, 1, MPI_INT, root, comm);
  if(ret!=0){
    return ret;
  }
  if(rank != root){
    val.resize(sz);
  }
  int ret=MPI_Bcast(&(val[0]), sz, datatype, root, comm);
#endif
  return ret;
}

template <class T>
int bcast(std::vector<std::complex<T>> &val, int root, MPI_Comm comm){
  int ret=0;
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  int rank;
  MPI_Comm_rank(com, &rank);
  int sz = val.size();
  ret=MPI_Bcast(&sz, 1, MPI_INT, root, comm);
  if(ret!=0){
    return ret;
  }
  if(rank != root){
    val.resize(sz);
  }
  std::vector<T> reim(2*sz);
  for(size_t i=0; i<sz; ++i){
    reim[2*i] = val[i].real();
    reim[2*i+1] = val[i].imag();
  }
  ret=MPI_Bcast(&(reim[0]), 2*sz, datatype, root, comm);
  for(size_t i=0; i<sz; ++i){
    val[i] = std::complex<T>(reim[2*i], reim[2*i+1]);
  }
#endif
  return ret;
}

template <class T>
int allreduce_sum(T &val, MPI_Comm comm){
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  T recv;
  int ret = MPI_Allreduce(&val, &recv, 1, datatype, MPI_SUM, comm);
  if(ret!=0){
    return ret;
  }
  val = recv;
#endif
  return 0;
}

template <class T>
int allreduce_sum(std::vector<T> &val, MPI_Comm comm){
#ifndef _NO_MPI
  const MPI_Datatype datatype = get_MPI_Datatype<T>();
  int N = val.size();
  std::vector<T> recv(val);
  int ret = MPI_Allreduce(&(val[0]), &(recv[0]), N, datatype, MPI_SUM, comm);
  if(ret!=0){
    return ret;
  }
  val.assign(recv.begin(), recv.end());
#endif
  return 0;
}

template <class T>
int allreduce_sum(std::complex<T> &val, MPI_Comm comm){
  throw tenes::unimplemented_error("allreduce for complex is not implemented");
}

template <class T>
int allreduce_sum(std::vector<std::complex<T>> &val, MPI_Comm comm){
  throw tenes::unimplemented_error("allreduce for complex is not implemented");
}

}


#endif  // TENES_MPI_H
