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
  bool to_measure;
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
