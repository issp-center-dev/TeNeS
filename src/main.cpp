#include <iostream>

#include "mpi.hpp"

#include "exception.hpp"

namespace tenes {
int main_impl(int argc, char **argv);
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int status = 0;

  try{
    status = tenes::main_impl(argc, argv);
  }catch(const tenes::input_error e){
    std::cerr << "[INPUT ERROR]" << std::endl;
    std::cerr << e.what() << std::endl;;
    status = 1;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return status;
}
