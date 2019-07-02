#include <vector>
#include <fstream>

#include <mpi.h>
#include <toml11/toml.hpp>

#include "util/find_or.hpp"
#include "Parameters.hpp"

Parameters::Parameters(){
  tau_simple = 0.0/0.0;
  num_simple_step = -1;

  tau_full = 0.0/0.0;
  num_full_step = 0;
}
Parameters::Parameters(const char *filename): Parameters(toml::parse(filename)){}
Parameters::Parameters(toml::Table data){
  set(data);
}

void Parameters::set(const char *filename){ set(toml::parse(filename)); }
void Parameters::set(toml::Table data){
  toml::Table param = toml::find<toml::Table>(data, "parameter");

  tau_simple = toml::find<double>(param, "tau_simple");
  num_simple_step = toml::find<int>(param, "num_simple_step");

  tau_full = toml::find<double>(param, "tau_full");
  num_full_step = util::find_or(param, "num_full_step", 0);
}

void Parameters::Bcast(MPI_Comm comm, int root) {
  int irank;
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);

  std::vector<double> params_double(4);
  std::vector<int> params_int(4);

  if (irank == root) {
    params_int[0] = num_simple_step;
    params_int[1] = num_full_step;

    params_double[0] = tau_simple;
    params_double[1] = tau_full;

    MPI_Bcast(&params_int.front(), 2, MPI_INT, 0, comm);
    MPI_Bcast(&params_double.front(), 2, MPI_DOUBLE, 0, comm);
  } else {
    MPI_Bcast(&params_int.front(), 2, MPI_INT, 0, comm);
    MPI_Bcast(&params_double.front(), 2, MPI_DOUBLE, 0, comm);

    num_simple_step = params_int[0];
    num_full_step = params_int[1];

    tau_simple = params_double[0];
    tau_full = params_double[1];
  }
}

void Parameters::save(const char *filename, const bool append) {
  std::ofstream ofs;
  if (append) {
    ofs.open(filename, std::ios::out | std::ios::app);
  } else {
    ofs.open(filename, std::ios::out);
  }

  // Tensor
  ofs << "tau_simple " << tau_simple << std::endl;
  ofs << "num_simple_step " << num_simple_step << std::endl;

  ofs << "tau_full " << tau_full << std::endl;
  ofs << "num_full_step " << num_full_step << std::endl;
  ofs.close();
}
