#ifndef _PEPS_PARAMETERS_HPP_
#define _PEPS_PARAMETERS_HPP_

#include <mpi.h>

//#include <toml11/toml/value.hpp>
#include <cpptoml.h>


class PEPS_Parameters {
public:
  // Tensor
  int D;    // Bond dimension for central tensor
  int CHI;  // Bond dimension for environment tensor

  // Debug
  bool Debug_flag;
  bool Warning_flag;

  // Simple update
  int num_simple_step;
  double Inverse_lambda_cut;

  // Environment
  double Inverse_projector_cut;
  double CTM_Convergence_Epsilon;
  int Max_CTM_Iteration;
  bool CTM_Projector_corner;
  bool Use_RSVD;
  int RSVD_Oversampling_factor;

  // Full update
  int num_full_step;
  double Inverse_Env_cut;
  double Full_Inverse_precision;
  double Full_Convergence_Epsilon;
  int Full_max_iteration;
  bool Full_Gauge_Fix;
  bool Full_Use_FastFullUpdate; // Fast Full Update

  // observable
  int Lcor;

  PEPS_Parameters();

  explicit PEPS_Parameters(const char *filename): PEPS_Parameters(){ set(filename); }
  PEPS_Parameters(decltype(cpptoml::parse_file("")) toml);
  void set(const char *filename);
  void set(decltype(cpptoml::parse_file("")) toml);

  void save(const char *filename, bool append=false);
  void save_append(const char *filename){save(filename, true);}

  void Bcast(MPI_Comm comm, int root=0);
};
#endif // _PEPS_PARAMETERS_HPP_
