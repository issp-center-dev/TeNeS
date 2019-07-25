#include <fstream>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

#include <toml11/toml.hpp>

#include "util/find_or.hpp"

#include "PEPS_Parameters.hpp"

PEPS_Parameters::PEPS_Parameters() {
  S = 0.5;

  // Tensor
  D = 2;
  CHI = D * D;

  // Debug
  Debug_flag = false;
  Warning_flag = true;

  // Simple update
  tau_simple = 0.0;
  num_simple_step = 0;
  Inverse_lambda_cut = 1e-12;

  // Environment
  Inverse_projector_cut = 1e-12;
  CTM_Convergence_Epsilon = 1e-10;
  Max_CTM_Iteration = 100;
  CTM_Projector_corner = false;
  Use_RSVD = false;
  RSVD_Oversampling_factor = 2;

  // Full update
  tau_full = 0.0;
  num_full_step = 0;
  Inverse_Env_cut = 1e-12;
  Full_Inverse_precision = 1e-12;
  Full_Convergence_Epsilon = 1e-12;
  Full_max_iteration = 1000;
  Full_Gauge_Fix = true;
  Full_Use_FFU = true;
}

PEPS_Parameters::PEPS_Parameters(toml::Table data): PEPS_Parameters() { set(data); }
void PEPS_Parameters::set(const char *filename){ set(toml::parse(filename)); }
void PEPS_Parameters::set(toml::Table data){
  toml::Table param = toml::find<toml::Table>(data, "parameter");

  S = util::find_or(param, "S", 0.5);
  if(2*S!=static_cast<int>(2*S)){
    std::cerr << "S should be an integer or a half-integer" << std::endl;
    std::exit(1);
  }

  // Tensor
  D = util::find_or(param, "D", 2);
  CHI = util::find_or(param, "CHI", 4);


  // Debug
  Debug_flag = util::find_or(param, "Debug", false);
  Warning_flag = util::find_or(param, "Warning", true);

  // Simple update
  tau_simple = util::find_or(param, "tau_simple", 0.0);
  num_simple_step = util::find_or(param, "num_simple_step", 0);
  Inverse_lambda_cut = util::find_or(param, "inverse_lambda_cutoff", 1e-12);

  // Environment
  Inverse_projector_cut = util::find_or(param, "inverse_projector_cutoff", 1e-12);
  CTM_Convergence_Epsilon = util::find_or(param, "ctm_convergence_epsilon", 1e-10);
  Max_CTM_Iteration = util::find_or(param, "ctm_iteration_max", 100);
  CTM_Projector_corner = util::find_or(param, "ctm_projector_corner", false);
  Use_RSVD = util::find_or(param, "use_rsvd", false);
  RSVD_Oversampling_factor = util::find_or(param, "rsvd_oversampling_factor", 2);

  // Full update
  tau_full = util::find_or(param, "tau_full", 0.0);
  num_full_step = util::find_or(param, "num_full_step", 0);
  Inverse_Env_cut = util::find_or(param, "inverse_projector_cutoff", 1e-12);
  Full_Inverse_precision = util::find_or(param, "full_inverse_precision", 1e-12);
  Full_Convergence_Epsilon = util::find_or(param, "full_convergence_epsilon", 1e-12);
  Full_max_iteration = util::find_or(param, "full_max_iteration", 1000);
  Full_Gauge_Fix = util::find_or(param, "full_gauge_fix", true);
  Full_Use_FFU = util::find_or(param, "full_use_ffu", true);
}

#define SAVE_PARAM(name, type) params_##type [I_##name] = name
#define LOAD_PARAM(name, type) name = params_##type [I_##name]

void PEPS_Parameters::Bcast(MPI_Comm comm, int root) {
  enum PARAMS_INT_INDEX{ 
    I_D,
    I_CHI,
    I_Debug_flag,
    I_Warning_flag,
    I_num_simple_step,
    I_Max_CTM_Iteration,
    I_CTM_Projector_corner,
    I_Use_RSVD,
    I_RSVD_Oversampling_factor,
    I_Full_max_iteration,
    I_Full_Gauge_Fix,
    I_Full_Use_FFU,

    N_PARAMS_INT_INDEX,
  };
  enum PARAMS_DOUBLE_INDEX{
    I_S,
    I_tau_simple,
    I_Inverse_lambda_cut,
    I_Inverse_projector_cut,
    I_CTM_Convergence_Epsilon,
    I_Inverse_Env_cut,
    I_tau_full,
    I_Full_Inverse_precision,
    I_Full_Convergence_Epsilon,

    N_PARAMS_DOUBLE_INDEX,
  };

  int irank;
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);

  std::vector<int> params_int(N_PARAMS_INT_INDEX);
  std::vector<double> params_double(N_PARAMS_DOUBLE_INDEX);

  if (irank == root) {
    SAVE_PARAM(D, int);
    SAVE_PARAM(CHI, int);
    SAVE_PARAM(Debug_flag, int);
    SAVE_PARAM(Warning_flag, int);
    SAVE_PARAM(num_simple_step, int);
    SAVE_PARAM(Max_CTM_Iteration, int);
    SAVE_PARAM(CTM_Projector_corner, int);
    SAVE_PARAM(Use_RSVD, int);
    SAVE_PARAM(RSVD_Oversampling_factor, int);
    SAVE_PARAM(Full_max_iteration, int);
    SAVE_PARAM(Full_Gauge_Fix, int);
    SAVE_PARAM(Full_Use_FFU, int);

    SAVE_PARAM(S, double);
    SAVE_PARAM(tau_simple, double);
    SAVE_PARAM(Inverse_lambda_cut, double);
    SAVE_PARAM(Inverse_projector_cut, double);
    SAVE_PARAM(CTM_Convergence_Epsilon, double);
    SAVE_PARAM(Inverse_Env_cut, double);
    SAVE_PARAM(tau_full, double);
    SAVE_PARAM(Full_Inverse_precision, double);
    SAVE_PARAM(Full_Convergence_Epsilon, double);


    MPI_Bcast(&params_int.front(), N_PARAMS_INT_INDEX, MPI_INT, 0, comm);
    MPI_Bcast(&params_double.front(), N_PARAMS_DOUBLE_INDEX, MPI_DOUBLE, 0, comm);
  } else {
    MPI_Bcast(&params_int.front(), N_PARAMS_INT_INDEX, MPI_INT, 0, comm);
    MPI_Bcast(&params_double.front(), N_PARAMS_DOUBLE_INDEX, MPI_DOUBLE, 0, comm);

    LOAD_PARAM(D, int);
    LOAD_PARAM(CHI, int);
    LOAD_PARAM(Debug_flag, int);
    LOAD_PARAM(Warning_flag, int);
    LOAD_PARAM(num_simple_step, int);
    LOAD_PARAM(Max_CTM_Iteration, int);
    LOAD_PARAM(CTM_Projector_corner, int);
    LOAD_PARAM(Use_RSVD, int);
    LOAD_PARAM(RSVD_Oversampling_factor, int);
    LOAD_PARAM(Full_max_iteration, int);
    LOAD_PARAM(Full_Gauge_Fix, int);
    LOAD_PARAM(Full_Use_FFU, int);

    LOAD_PARAM(S, double);
    LOAD_PARAM(tau_simple, double);
    LOAD_PARAM(Inverse_lambda_cut, double);
    LOAD_PARAM(Inverse_projector_cut, double);
    LOAD_PARAM(CTM_Convergence_Epsilon, double);
    LOAD_PARAM(Inverse_Env_cut, double);
    LOAD_PARAM(tau_full, double);
    LOAD_PARAM(Full_Inverse_precision, double);
    LOAD_PARAM(Full_Convergence_Epsilon, double);

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
  ofs << "S " << S << std::endl;
  // Tensor
  ofs << "D " << D << std::endl;
  ofs << "CHI " << CHI << std::endl;

  // Debug
  ofs << "Debug_flag " << Debug_flag << std::endl;
  ofs << "Warning_flag " << Warning_flag << std::endl;

  // Simple update
  ofs << "Inverse_lambda_cut " << Inverse_lambda_cut << std::endl;

  // Environment
  ofs << "Inverse_projector_cut " << Inverse_projector_cut << std::endl;
  ofs << "CTM_Convergence_Epsilon " << CTM_Convergence_Epsilon << std::endl;
  ofs << "Max_CTM_Iteration " << Max_CTM_Iteration << std::endl;
  ofs << "CTM_Projector_corner " << CTM_Projector_corner << std::endl;
  ofs << "Use_RSVD " << Use_RSVD << std::endl;
  ofs << "RSVD_Oversampling_factor " << RSVD_Oversampling_factor << std::endl;

  // Full update
  ofs << "Inverse_Env_cut " << Inverse_Env_cut << std::endl;
  ofs << "Full_Inverse_precision " << Full_Inverse_precision << std::endl;
  ofs << "Full_Convergence_Epsilon " << Full_Convergence_Epsilon << std::endl;
  ofs << "Full_max_iteration " << Full_max_iteration << std::endl;
  ofs << "Full_Gauge_Fix " << Full_Gauge_Fix << std::endl;
  ofs << "Full_Use_FFU " << Full_Use_FFU << std::endl;

  ofs.close();
}
