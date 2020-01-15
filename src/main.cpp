#include "mpi.hpp"

namespace tenes{
int main_impl(int argc, char **argv);
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int status = tenes::main_impl(argc, argv);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return status;
}
