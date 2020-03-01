#include <iostream>

#include "version.hpp"

#include "mpi.hpp"

#include "exception.hpp"
#include "printlevel.hpp"

namespace tenes {
int main_impl(std::string input_filename, MPI_Comm com, PrintLevel print_level);
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int mpirank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int status = 0;
  try{
    std::string usage =
        R"(TeNeS: TEnsor NEtwork Solver for 2D quantum lattice system
    
    Usage:
      tenes [--quiet] <input_toml>
      tenes --help
      tenes --version

    Options:
      -h --help       Show this help message.
      -v --version    Show the version.
      -q --quiet      Do not print any messages.
    )";

    if (argc == 1) {
      if (mpirank == 0)
        std::cout << usage << std::endl;
      return 0;
    }

    for (int i = 1; i < argc; ++i) {
      std::string opt = argv[i];
      if (opt == "-h" || opt == "--help") {
        if (mpirank == 0)
          std::cout << usage << std::endl;
        return 0;
      }
    }

    for (int i = 1; i < argc; ++i) {
      std::string opt = argv[i];
      if (opt == "-v" || opt == "--version") {
        if (mpirank == 0)
          std::cout << "TeNeS v" << TENES_VERSION << std::endl;
        return 0;
      }
    }

    using PrintLevel = tenes::PrintLevel;

    PrintLevel print_level = PrintLevel::info;
    std::string input_filename;
    for (int i = 1; i < argc; ++i) {
      std::string opt = argv[i];
      if (opt == "-q" || opt == "--quiet") {
        print_level = PrintLevel::none;
      } else {
        input_filename = opt;
      }
    }

    status = tenes::main_impl(input_filename, MPI_COMM_WORLD, print_level);
  }catch(const tenes::input_error e){
    if(mpirank==0){
      std::cerr << "[INPUT ERROR]" << std::endl;
      std::cerr << e.what() << std::endl;;
    }
    status = 1;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return status;
}
