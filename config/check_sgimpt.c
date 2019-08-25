/* This file is originally included in mptensor */
/* https://github.com/smorita/mptensor */

#include <mpi.h>

#ifdef OMPI_MPI_H // OpenMPI
#error
#endif

#ifdef MPICH2 // MPICH2 or Intel MPI
#error
#endif

#ifdef MVAPICH2_VERSION //MVAPICH2
#error
#endif

int main(int argc, char** argv) {
  return 0;
}

