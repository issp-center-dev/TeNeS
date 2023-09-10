/* TeNeS - Massively parallel tensor network solver /
/ Copyright (C) 2019- The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify /
/ it under the terms of the GNU General Public License as published by /
/ the Free Software Foundation, either version 3 of the License, or /
/ (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, /
/ but WITHOUT ANY WARRANTY; without even the implied warranty of /
/ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the /
/ GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License /
/ along with this program. If not, see http://www.gnu.org/licenses/. */

#include "PEPS_Parameters.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../exception.hpp"

namespace tenes {
namespace itps {

PEPS_Parameters::PEPS_Parameters() {
  // Tensor
  CHI = 2;

  // Debug
  print_level = PrintLevel::info;

  // Simple update
  num_simple_step = std::vector<int>{0};
  tau_simple_step = std::vector<double>{0.0};
  measure_interval = std::vector<int>{10};
  Inverse_lambda_cut = 1e-12;
  Simple_Gauge_Fix = false;
  Simple_Gauge_maxiter = 100;
  Simple_Gauge_Convergence_Epsilon = 1.0e-2;

  // Environment
  Inverse_projector_cut = 1e-12;
  CTM_Convergence_Epsilon = 1e-6;
  Max_CTM_Iteration = 100;
  CTM_Projector_corner = true;
  Use_RSVD = false;
  RSVD_Oversampling_factor = 2.0;
  MeanField_Env = false;

  // Full update
  num_full_step = std::vector<int>{0};
  tau_full_step = std::vector<double>{0.0};
  Inverse_Env_cut = 1e-12;
  Full_Inverse_precision = 1e-12;
  Full_Convergence_Epsilon = 1e-6;
  Full_max_iteration = 100;
  Full_Gauge_Fix = true;
  Full_Use_FastFullUpdate = true;

  // observable
  corder_opt = ContractionOrderOptimization::automatic;
  contraction_order_file = "";

  // random
  seed = 11;

  // general
  is_real = false;
  iszero_tol = 0.0;
  to_measure = true;
  tensor_load_dir = "";
  tensor_save_dir = "";
  outdir = "output";
  calcmode = CalculationMode::ground_state;
}

#define SAVE_PARAM(name, type) params_##type[I_##name] = static_cast<type>(name)
#define LOAD_PARAM(name, type) \
  name = static_cast<decltype(name)>(params_##type[I_##name])

void PEPS_Parameters::Bcast(MPI_Comm comm, int root) {
  using std::string;

  enum PARAMS_INT_INDEX {
    I_CHI,
    I_print_level,
    I_num_simple_step,
    I_Simple_Gauge_Fix,
    I_Simple_Gauge_maxiter,
    I_Max_CTM_Iteration,
    I_CTM_Projector_corner,
    I_Use_RSVD,
    I_MeanField_Env,
    I_num_full_step,
    I_Full_max_iteration,
    I_Full_Gauge_Fix,
    I_Full_Use_FastFullUpdate,
    I_seed,
    I_is_real,
    I_to_measure,

    N_PARAMS_INT_INDEX,
  };
  enum PARAMS_DOUBLE_INDEX {
    I_Inverse_lambda_cut,
    I_Simple_Gauge_Convergence_Epsilon,
    I_Inverse_projector_cut,
    I_CTM_Convergence_Epsilon,
    I_Inverse_Env_cut,
    I_Full_Inverse_precision,
    I_Full_Convergence_Epsilon,
    I_RSVD_Oversampling_factor,
    I_iszero_tol,

    N_PARAMS_DOUBLE_INDEX,
  };

  enum PARAMS_STRING_INDEX {
    I_tensor_load_dir,
    I_tensor_save_dir,
    I_outdir,

    N_PARAMS_STRING_INDEX,
  };

  int irank;
  MPI_Comm_rank(comm, &irank);

  std::vector<int> params_int(N_PARAMS_INT_INDEX);
  std::vector<double> params_double(N_PARAMS_DOUBLE_INDEX);
  std::vector<std::string> params_string(N_PARAMS_STRING_INDEX);

  if (irank == root) {
    SAVE_PARAM(CHI, int);
    SAVE_PARAM(print_level, int);
    // SAVE_PARAM(num_simple_step, int);
    SAVE_PARAM(Simple_Gauge_Fix, int);
    SAVE_PARAM(Simple_Gauge_maxiter, int);
    SAVE_PARAM(Max_CTM_Iteration, int);
    SAVE_PARAM(CTM_Projector_corner, int);
    SAVE_PARAM(Use_RSVD, int);
    SAVE_PARAM(MeanField_Env, int);
    // SAVE_PARAM(num_full_step, int);
    SAVE_PARAM(Full_max_iteration, int);
    SAVE_PARAM(Full_Gauge_Fix, int);
    SAVE_PARAM(Full_Use_FastFullUpdate, int);
    SAVE_PARAM(seed, int);

    SAVE_PARAM(Inverse_lambda_cut, double);
    SAVE_PARAM(Simple_Gauge_Convergence_Epsilon, double);
    SAVE_PARAM(Inverse_projector_cut, double);
    SAVE_PARAM(CTM_Convergence_Epsilon, double);
    SAVE_PARAM(Inverse_Env_cut, double);
    SAVE_PARAM(Full_Inverse_precision, double);
    SAVE_PARAM(Full_Convergence_Epsilon, double);
    SAVE_PARAM(RSVD_Oversampling_factor, double);

    SAVE_PARAM(is_real, int);
    SAVE_PARAM(iszero_tol, double);
    SAVE_PARAM(to_measure, int);
    SAVE_PARAM(tensor_load_dir, string);
    SAVE_PARAM(tensor_save_dir, string);
    SAVE_PARAM(outdir, string);

    bcast(params_int, 0, comm);
    bcast(params_double, 0, comm);
    bcast(params_string, 0, comm);
    // MPI_Bcast(&params_int.front(), N_PARAMS_INT_INDEX, MPI_INT, 0, comm);
    // MPI_Bcast(&params_double.front(), N_PARAMS_DOUBLE_INDEX, MPI_DOUBLE, 0,
    //           comm);
  } else {
    bcast(params_int, 0, comm);
    bcast(params_double, 0, comm);
    bcast(params_string, 0, comm);
    // MPI_Bcast(&params_int.front(), N_PARAMS_INT_INDEX, MPI_INT, 0, comm);
    // MPI_Bcast(&params_double.front(), N_PARAMS_DOUBLE_INDEX, MPI_DOUBLE, 0,
    //           comm);

    LOAD_PARAM(CHI, int);
    LOAD_PARAM(print_level, int);
    // LOAD_PARAM(num_simple_step, int);
    LOAD_PARAM(Simple_Gauge_Fix, int);
    LOAD_PARAM(Simple_Gauge_maxiter, int);
    LOAD_PARAM(Max_CTM_Iteration, int);
    LOAD_PARAM(CTM_Projector_corner, int);
    LOAD_PARAM(Use_RSVD, int);
    LOAD_PARAM(MeanField_Env, int);
    // LOAD_PARAM(num_full_step, int);
    LOAD_PARAM(Full_max_iteration, int);
    LOAD_PARAM(Full_Gauge_Fix, int);
    LOAD_PARAM(Full_Use_FastFullUpdate, int);
    LOAD_PARAM(seed, int);

    LOAD_PARAM(Inverse_lambda_cut, double);
    LOAD_PARAM(Simple_Gauge_Convergence_Epsilon, double);
    LOAD_PARAM(Inverse_projector_cut, double);
    LOAD_PARAM(CTM_Convergence_Epsilon, double);
    LOAD_PARAM(Inverse_Env_cut, double);
    LOAD_PARAM(Full_Inverse_precision, double);
    LOAD_PARAM(Full_Convergence_Epsilon, double);
    LOAD_PARAM(RSVD_Oversampling_factor, double);

    LOAD_PARAM(is_real, int);
    LOAD_PARAM(iszero_tol, double);
    LOAD_PARAM(to_measure, int);
    LOAD_PARAM(tensor_load_dir, string);
    LOAD_PARAM(tensor_save_dir, string);
    LOAD_PARAM(outdir, string);
  }
}

#undef SAVE_PARAM
#undef LOAD_PARAM

void PEPS_Parameters::save(const char *filename, bool append) {
  std::ofstream ofs;
  if (append) {
    ofs.open(filename, std::ios::out | std::ios::app);
  } else {
    ofs.open(filename, std::ios::out);
  }

  // Simple update
  ofs << "simple_num_step = [";
  for (int i = 0; i < num_simple_step.size(); ++i) {
    if (i != 0) ofs << ", ";
    ofs << num_simple_step[i];
  }
  ofs << "]" << std::endl;
  ofs << "simple_tau = [";
  for (int i = 0; i < tau_simple_step.size(); ++i) {
    if (i != 0) ofs << ", ";
    ofs << tau_simple_step[i];
  }
  ofs << "]" << std::endl;
  ofs << "simple_inverse_lambda_cutoff = " << Inverse_lambda_cut << std::endl;
  ofs << "simple_gauge_fix = " << Simple_Gauge_Fix << std::endl;
  ofs << "simple_gauge_maxiter = " << Simple_Gauge_maxiter << std::endl;
  ofs << "simple_gauge_convergence_epsilon = "
      << Simple_Gauge_Convergence_Epsilon << std::endl;

  ofs << std::endl;

  // Full update
  ofs << "full_num_step = [";
  for (int i = 0; i < num_full_step.size(); ++i) {
    if (i != 0) ofs << ", ";
    ofs << num_full_step[i];
  }
  ofs << "]" << std::endl;
  ofs << "full_inverse_projector_cutoff = " << Inverse_projector_cut
      << std::endl;
  ofs << "full_inverse_precision = " << Full_Inverse_precision << std::endl;
  ofs << "full_convergence_epsilon = " << Full_Convergence_Epsilon << std::endl;
  ofs << "full_iteration_max = " << Full_max_iteration << std::endl;
  ofs << "full_gauge_fix = " << (Full_Gauge_Fix ? "true" : "false")
      << std::endl;
  ofs << "full_fastfullupdate = "
      << (Full_Use_FastFullUpdate ? "true" : "false") << std::endl;

  ofs << std::endl;

  // Environment
  ofs << "ctm_dimension = " << CHI << std::endl;
  ofs << "ctm_inverse_projector_cutoff = " << Inverse_Env_cut << std::endl;
  ofs << "ctm_convergence_epsilon = " << CTM_Convergence_Epsilon << std::endl;
  ofs << "ctm_iteration_max = " << Max_CTM_Iteration << std::endl;
  ofs << "ctm_projector_corner = " << (CTM_Projector_corner ? "true" : "false")
      << std::endl;
  ofs << "use_rsvd = " << (Use_RSVD ? "true" : "false") << std::endl;
  ofs << "rsvd_oversampling_factor = " << RSVD_Oversampling_factor << std::endl;
  ofs << "meanfield_env = " << (MeanField_Env ? "true" : "false") << std::endl;

  ofs << std::endl;

  ofs << "mode = ";
  switch(calcmode){
    case PEPS_Parameters::CalculationMode::ground_state:
      ofs << "ground state" << std::endl;
      break;
    case PEPS_Parameters::CalculationMode::time_evolution:
      ofs << "time evolution" << std::endl;
      break;
    case PEPS_Parameters::CalculationMode::finite_temperature:
      ofs << "finite temperature" << std::endl;
      break;
    default:
      break;
  }
  ofs << "simple" << std::endl;
  switch(corder_opt){
    case PEPS_Parameters::ContractionOrderOptimization::automatic:
      ofs << "contraction_optimize = automatic" << std::endl;
      break;
    case PEPS_Parameters::ContractionOrderOptimization::never:
      ofs << "contraction_optimize = never" << std::endl;
      break;
    case PEPS_Parameters::ContractionOrderOptimization::always:
      ofs << "contraction_optimize = always" << std::endl;
      break;
    case PEPS_Parameters::ContractionOrderOptimization::old:
      ofs << "contraction_mode = old" << std::endl;
      break;
    default:
      break;
  }
  ofs << "contraction_order_file = " << contraction_order_file << std::endl;
  ofs << "seed = " << seed << std::endl;
  ofs << "is_real = " << is_real << std::endl;
  ofs << "iszero_tol = " << iszero_tol << std::endl;
  ofs << "measure = " << to_measure << std::endl;
  ofs << "tensor_load_dir = " << tensor_load_dir << std::endl;
  ofs << "tensor_save_dir = " << tensor_save_dir << std::endl;
  ofs << "outdir = " << outdir << std::endl;

  ofs.close();
}

void PEPS_Parameters::check() const {
  if (MeanField_Env && num_full_step[0] > 0) {
    std::stringstream ss;
    ss << "ERROR: Cannot enable full update and mean field environment "
          "simultaneously";
    throw tenes::input_error(ss.str());
  }
}

}  // namespace itps
}  // namespace tenes
