#include "mpi.hpp"

#ifdef _NO_MPI
int MPI_Init(int*, char***) { return 0; }
int MPI_Barrier(MPI_Comm) { return 0; }
int MPI_Finalize() { return 0; }

int MPI_Comm_size(MPI_Comm, int* size) {
  *size = 1;
  return 0;
}
int MPI_Comm_rank(MPI_Comm, int* rank) {
  *rank = 0;
  return 0;
}

int MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm) { return 0; }
#endif

namespace tenes{

template <>
MPI_Datatype get_MPI_Datatype<int>() { return MPI_INT; }
template <>
MPI_Datatype get_MPI_Datatype<double>() { return MPI_DOUBLE; }
template <>
MPI_Datatype get_MPI_Datatype<bool>() { return MPI_INT; }

int bcast(bool &val, int root, MPI_Comm comm){
  int ret=0;
#ifndef _NO_MPI
  int v = val?0:1;
  ret=MPI_Bcast(&v, 1, MPI_INT, root, comm);
  val = v==0;
#endif
  return ret;
}

}
