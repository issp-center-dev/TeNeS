#ifndef PARAMETERS_HPP

#include <toml11/toml.hpp>
#include "util/find_or.hpp"

struct Parameters{

  double tau_simple;
  int tau_simple_step;

  double tau_full;
  int tau_full_step;

  Parameters(){
    tau_simple = 0.0/0.0;
    tau_simple_step = -1;

    tau_full = 0.0/0.0;
    tau_full_step = 0;
  }

  explicit Parameters(const char *filename): Parameters(toml::parse(filename)){}
  Parameters(toml::table data){
    set(data);
  }

  void set(const char *filename){ set(toml::parse(filename)); }
  void set(toml::table data){
    toml::table param = toml::find<toml::table>(data, "parameter");

    tau_simple = toml::find<double>(param, "tau_simple");
    tau_simple_step = toml::find<int>(param, "tau_simple_step");

    tau_full = toml::find<double>(param, "tau_full");
    tau_full_step = util::find_or(param, "tau_full_step", 0);
  }

  void output_parameters(const char *filename, const bool append) {
    std::ofstream ofs;
    if (append) {
      ofs.open(filename, std::ios::out | std::ios::app);
    } else {
      ofs.open(filename, std::ios::out);
    }

    // Tensor
    ofs << "tau_simple " << tau_simple << std::endl;
    ofs << "tau_simple_step " << tau_simple_step << std::endl;

    ofs << "tau_full " << tau_full << std::endl;
    ofs << "tau_full_step " << tau_full_step << std::endl;
    ofs.close();
  }

  void output_parameters(const char *filename) {
    output_parameters(filename, false);
  }

  void output_parameters_append(const char *filename) {
    output_parameters(filename, true);
  }
  void Bcast_parameters(MPI_Comm comm) {
    int irank;
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);

    std::vector<double> params_double(4);
    std::vector<int> params_int(4);

    if (irank == 0) {
      params_int[0] = tau_simple_step;
      params_int[1] = tau_full_step;

      params_double[0] = tau_simple;
      params_double[1] = tau_full;

      MPI_Bcast(&params_int.front(), 2, MPI_INT, 0, comm);
      MPI_Bcast(&params_double.front(), 2, MPI_DOUBLE, 0, comm);
    } else {
      MPI_Bcast(&params_int.front(), 2, MPI_INT, 0, comm);
      MPI_Bcast(&params_double.front(), 2, MPI_DOUBLE, 0, comm);

      tau_simple_step = params_int[0];
      tau_full_step = params_int[1];

      tau_simple = params_double[0];
      tau_full = params_double[1];
    }
  }

};

#define PARAMETERS_HPP
#endif // PARAMETERS_HPP
