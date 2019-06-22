#include <fstream>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

#include <toml11/toml.hpp>

#include "util/find_or.hpp"

#include "PEPS_Parameters.hpp"

PEPS_Parameters::PEPS_Parameters(toml::Table data): PEPS_Parameters() { set(data); }
void PEPS_Parameters::set(const char *filename){ set(toml::parse(filename)); }
void PEPS_Parameters::set(toml::Table data){
  toml::Table param = toml::find<toml::Table>(data, "parameter");

  // Tensor
  D = util::find_or(param, "D", 2);
  CHI = util::find_or(param, "CHI", 4);


  // Debug
  Debug_flag = util::find_or(param, "Debug", false);
  Warning_flag = util::find_or(param, "Warning", true);

  // Simple update
  Inverse_lambda_cut = util::find_or(param, "inverse_lambda_cutoff", 1e-12);

  // Environment
  Inverse_projector_cut = util::find_or(param, "inverse_projector_cutoff", 1e-12);
  CTM_Convergence_Epsilon = util::find_or(param, "ctm_convergence_epsilon", 1e-10);
  Max_CTM_Iteration = util::find_or(param, "ctm_iteration_max", 100);
  CTM_Projector_corner = util::find_or(param, "ctm_projector_corner", false);
  Use_RSVD = util::find_or(param, "use_rsvd", false);
  RSVD_Oversampling_factor = util::find_or(param, "rsvd_oversampling_factor", 2);

  // Full update
  Inverse_Env_cut = util::find_or(param, "inverse_projector_cutoff", 1e-12);
  Full_Inverse_precision = util::find_or(param, "full_inverse_precision", 1e-12);
  Full_Convergence_Epsilon = util::find_or(param, "full_convergence_epsilon", 1e-12);
  Full_max_iteration = util::find_or(param, "full_max_iteration", 1000);
  Full_Gauge_Fix = util::find_or(param, "full_gauge_fix", true);
  Full_Use_FFU = util::find_or(param, "full_use_ffu", true);
}
