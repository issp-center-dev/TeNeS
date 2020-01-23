#ifdef _NO_MPI

#include "mpi.hpp"

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
