/* This file is originally included in mptensor */
/* https://github.com/smorita/mptensor */

#include <mpi.h>

#ifdef OMPI_MPI_H // OpenMPI
int main(int argc, char** argv) {
  return 0;
}
#else
#error
#endif
