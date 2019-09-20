#include <fstream>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

#include "PEPS_Parameters.hpp"

PEPS_Parameters::PEPS_Parameters() {
  // Tensor
  D = 2;
  CHI = D * D;

  // Debug
  print_level = PrintLevel::info;

  // Simple update
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
  num_full_step = 0;
  Inverse_Env_cut = 1e-12;
  Full_Inverse_precision = 1e-12;
  Full_Convergence_Epsilon = 1e-12;
  Full_max_iteration = 1000;
  Full_Gauge_Fix = true;
  Full_Use_FastFullUpdate = true;

  Lcor = 0;
}

#define SAVE_PARAM(name, type) params_##type[I_##name] = static_cast<type>(name)
#define LOAD_PARAM(name, type) name = static_cast<decltype(name)>(params_##type[I_##name])

void PEPS_Parameters::Bcast(MPI_Comm comm, int root) {
  enum PARAMS_INT_INDEX {
    I_D,
    I_CHI,
    I_print_level,
    I_num_simple_step,
    I_Max_CTM_Iteration,
    I_CTM_Projector_corner,
    I_Use_RSVD,
    I_RSVD_Oversampling_factor,
    I_Full_max_iteration,
    I_Full_Gauge_Fix,
    I_Full_Use_FastFullUpdate,
    I_Lcor,

    N_PARAMS_INT_INDEX,
  };
  enum PARAMS_DOUBLE_INDEX {
    I_Inverse_lambda_cut,
    I_Inverse_projector_cut,
    I_CTM_Convergence_Epsilon,
    I_Inverse_Env_cut,
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
    SAVE_PARAM(print_level, int);
    SAVE_PARAM(num_simple_step, int);
    SAVE_PARAM(Max_CTM_Iteration, int);
    SAVE_PARAM(CTM_Projector_corner, int);
    SAVE_PARAM(Use_RSVD, int);
    SAVE_PARAM(RSVD_Oversampling_factor, int);
    SAVE_PARAM(Full_max_iteration, int);
    SAVE_PARAM(Full_Gauge_Fix, int);
    SAVE_PARAM(Full_Use_FastFullUpdate, int);
    SAVE_PARAM(Lcor, int);

    SAVE_PARAM(Inverse_lambda_cut, double);
    SAVE_PARAM(Inverse_projector_cut, double);
    SAVE_PARAM(CTM_Convergence_Epsilon, double);
    SAVE_PARAM(Inverse_Env_cut, double);
    SAVE_PARAM(Full_Inverse_precision, double);
    SAVE_PARAM(Full_Convergence_Epsilon, double);

    MPI_Bcast(&params_int.front(), N_PARAMS_INT_INDEX, MPI_INT, 0, comm);
    MPI_Bcast(&params_double.front(), N_PARAMS_DOUBLE_INDEX, MPI_DOUBLE, 0,
              comm);
  } else {
    MPI_Bcast(&params_int.front(), N_PARAMS_INT_INDEX, MPI_INT, 0, comm);
    MPI_Bcast(&params_double.front(), N_PARAMS_DOUBLE_INDEX, MPI_DOUBLE, 0,
              comm);

    LOAD_PARAM(D, int);
    LOAD_PARAM(CHI, int);
    LOAD_PARAM(print_level, int);
    LOAD_PARAM(num_simple_step, int);
    LOAD_PARAM(Max_CTM_Iteration, int);
    LOAD_PARAM(CTM_Projector_corner, int);
    LOAD_PARAM(Use_RSVD, int);
    LOAD_PARAM(RSVD_Oversampling_factor, int);
    LOAD_PARAM(Full_max_iteration, int);
    LOAD_PARAM(Full_Gauge_Fix, int);
    LOAD_PARAM(Full_Use_FastFullUpdate, int);
    LOAD_PARAM(Lcor, int);

    LOAD_PARAM(Inverse_lambda_cut, double);
    LOAD_PARAM(Inverse_projector_cut, double);
    LOAD_PARAM(CTM_Convergence_Epsilon, double);
    LOAD_PARAM(Inverse_Env_cut, double);
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
  // Tensor
  ofs << "D = " << D << std::endl;
  ofs << "CHI = " << CHI << std::endl;

  // Simple update
  ofs << "simple_num_step = " << num_simple_step << std::endl;
  ofs << "simple_inverse_lambda_cutoff = " << Inverse_lambda_cut << std::endl;

  // Full update
  ofs << "full_num_step = " << num_full_step << std::endl;
  ofs << "full_inverse_projector_cutoff = " << Inverse_projector_cut
      << std::endl;
  ofs << "full_inverse_precision = " << Full_Inverse_precision << std::endl;
  ofs << "full_convergence_epsilon = " << Full_Convergence_Epsilon << std::endl;
  ofs << "full_iteration_max = " << Full_max_iteration << std::endl;
  ofs << "full_gauge_fix = " << (Full_Gauge_Fix ? "true" : "false")
      << std::endl;
  ofs << "full_fastfullupdate = "
      << (Full_Use_FastFullUpdate ? "true" : "false") << std::endl;

  // Environment
  ofs << "ctm_inverse_projector_cutoff = " << Inverse_Env_cut << std::endl;
  ofs << "ctm_convergence_epsilon = " << CTM_Convergence_Epsilon << std::endl;
  ofs << "ctm_iteration_max = " << Max_CTM_Iteration << std::endl;
  ofs << "ctm_projector_corner = " << (CTM_Projector_corner ? "true" : "false")
      << std::endl;
  ofs << "use_rsvd = " << (Use_RSVD ? "true" : "false") << std::endl;
  ofs << "rsvd_oversampling_factor = " << RSVD_Oversampling_factor << std::endl;

  ofs.close();
}
