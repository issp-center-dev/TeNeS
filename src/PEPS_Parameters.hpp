#ifndef _PEPS_PARAMETERS_HPP_
#define _PEPS_PARAMETERS_HPP_

#include <string>

#include "mpi.hpp"
#include "printlevel.hpp"

namespace tenes {

class PEPS_Parameters {
 public:
  // Tensor
  int D;    // Bond dimension for central tensor
  int CHI;  // Bond dimension for environment tensor

  PrintLevel print_level;

  // Simple update
  int num_simple_step;
  double Inverse_lambda_cut;

  // Environment
  double Inverse_projector_cut;
  double CTM_Convergence_Epsilon;
  int Max_CTM_Iteration;
  bool CTM_Projector_corner;
  bool Use_RSVD;
  double RSVD_Oversampling_factor;

  // Full update
  int num_full_step;
  double Inverse_Env_cut;
  double Full_Inverse_precision;
  double Full_Convergence_Epsilon;
  int Full_max_iteration;
  bool Full_Gauge_Fix;
  bool Full_Use_FastFullUpdate;  // Fast Full Update

  // observable
  int Lcor;

  // random
  int seed;

  // general
  bool is_real;
  double iszero_tol;
  std::string tensor_load_dir;
  std::string tensor_save_dir;
  std::string outdir;

  PEPS_Parameters();

  void save(const char *filename, bool append = false);
  void save_append(const char *filename) { save(filename, true); }

  void Bcast(MPI_Comm comm, int root = 0);
};

}  // end of namespace tenes
#endif  // _PEPS_PARAMETERS_HPP_
