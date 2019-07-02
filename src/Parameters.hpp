#ifndef PARAMETERS_HPP

#include <mpi.h>

#include <toml11/toml/types.hpp>

struct Parameters{

  double tau_simple;
  int num_simple_step;

  double tau_full;
  int num_full_step;

  Parameters();
  explicit Parameters(const char *filename);
  explicit Parameters(toml::Table data);

  void set(const char *filename);
  void set(toml::Table data);

  void save(const char *filename, const bool append=false);

  void save_append(const char *filename) {
    save(filename, true);
  }

  void Bcast(MPI_Comm comm, int root=0);
};

#define PARAMETERS_HPP
#endif // PARAMETERS_HPP
