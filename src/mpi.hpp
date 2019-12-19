#ifndef TENES_MPI_H
#define TENES_MPI_H

#ifdef _NO_MPI

using MPI_Comm = int;
using MPI_Datatype = int;

constexpr MPI_Comm MPI_COMM_WORLD=0;
constexpr MPI_Datatype MPI_INT=0;
constexpr MPI_Datatype MPI_DOUBLE=0;

int MPI_Init(int*, char***);
int MPI_Barrier(MPI_Comm);
int MPI_Finalize();

int MPI_Comm_size(MPI_Comm, int*);
int MPI_Comm_rank(MPI_Comm, int*);

int MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm);

#else

#include <mpi.h>

#endif  //_NO_MPI

#endif  // TENES_MPI_H
