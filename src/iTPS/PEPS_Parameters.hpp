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

/*! @file
 *  @brief Runtime parameters of the solver
 *         (::tenes::itps::PEPS_Parameters).
 */

#ifndef TENES_SRC_PEPS_PARAMETERS_HPP_
#define TENES_SRC_PEPS_PARAMETERS_HPP_

#include <string>
#include <vector>

#include "../mpi.hpp"
#include "../printlevel.hpp"

namespace tenes::itps {

/*! @brief Runtime parameters of the solver, read from the [parameter]
 *         table of input.toml.
 *
 *  Each member names the input.toml key it comes from (see the file
 *  specification in the manual); the constructor sets the documented
 *  defaults.
 */
struct PEPS_Parameters {
  // Tensor
  int D;    //!< Virtual bond dimension of the center tensors
  int CHI;  //!< Bond dimension of the CTM environment (ctm.dimension)

  PrintLevel print_level;  //!< Verbosity of the standard output

  /*! @name Simple update (parameter.simple_update) */
  //!@{
  //! Number of simple-update steps, per evolution-operator group
  //! (num_step).
  std::vector<int> num_simple_step;
  //! (Imaginary) time step of each evolution-operator group (tau).
  std::vector<double> tau_simple_step;
  //! Steps between measurements in time-evolution / finite-temperature
  //! runs (measure_interval).
  std::vector<int> measure_interval;
  //! Mean-field weights below this are treated as zero when inverted
  //! (lambda_cutoff).
  double Inverse_lambda_cut;
  //! Fix the tensor gauge after each step (gauge_fix).
  bool Simple_Gauge_Fix;
  //! Maximum number of gauge-fixing iterations (gauge_maxiter).
  int Simple_Gauge_maxiter;
  //! Convergence criterion of the gauge-fixing iteration
  //! (gauge_converge_epsilon).
  double Simple_Gauge_Convergence_Epsilon;
  //!@}

  /*! @name CTM environment (parameter.ctm) */
  //!@{
  //! Singular values below this are treated as zero when computing CTM
  //! projectors (projector_cutoff).
  double Inverse_projector_cut;
  //! CTM convergence criterion (convergence_epsilon).
  double CTM_Convergence_Epsilon;
  //! Maximum number of CTM iterations (iteration_max).
  int Max_CTM_Iteration;
  //! Use only the 1/4 corner tensor in the projector calculation
  //! (projector_corner).
  bool CTM_Projector_corner;
  //! Replace SVD with randomized SVD (use_rsvd).
  bool Use_RSVD;
  //! Include one-site RDM distance in the CTM convergence check
  //! (use_onesite_rdm_convergence).
  bool CTM_Convergence_Onesite_RDM;
  //! Oversampling ratio of the randomized SVD
  //! (rsvd_oversampling_factor).
  double RSVD_Oversampling_factor;
  //! Measure with the mean-field environment from the simple update
  //! instead of the CTM (meanfield_env).
  bool MeanField_Env;
  //!@}

  /*! @name Full update (parameter.full_update) */
  //!@{
  //! Number of full-update steps, per evolution-operator group
  //! (num_step).
  std::vector<int> num_full_step;
  //! (Imaginary) time step of each evolution-operator group (tau).
  std::vector<double> tau_full_step;
  //! Singular values below this are treated as zero when computing the
  //! environment (env_cutoff).
  double Inverse_Env_cut;
  //! Cutoff for the pseudoinverse in the truncation optimization
  //! (inverse_precision).
  double Full_Inverse_precision;
  //! Convergence criterion of the truncation optimization
  //! (convergence_epsilon).
  double Full_Convergence_Epsilon;
  //! Maximum number of truncation-optimization iterations
  //! (iteration_max).
  int Full_max_iteration;
  //! Fix the tensor gauge before truncating (gauge_fix).
  bool Full_Gauge_Fix;
  //! Use the fast full update (fastfullupdate).
  bool Full_Use_FastFullUpdate;
  //!@}

  // observable
  //! Legacy field kept in the params.dat save format; the correlation
  //! distance actually used comes from CorrelationParameter::r_max.
  int Lcor;

  // random
  //! Seed of the pseudo-random number generator used to initialize the
  //! tensors (random.seed); every rank seeds identically and picks its
  //! locally stored elements from the same sequence.
  int seed;

  /*! @name General (parameter.general) */
  //!@{
  bool is_real;       //!< Restrict all tensors to real values (is_real)
  double iszero_tol;  //!< Operator elements below this read as zero
                      //!< (iszero_tol)
  bool to_measure;    //!< Calculate and save observables (measure)
  std::string tensor_load_dir;  //!< Directory to load initial tensors
                                //!< from (tensor_load); empty: none
  std::string tensor_save_dir;  //!< Directory to save optimized tensors
                                //!< into (tensor_save); empty: none
  std::string outdir;  //!< Directory the results are written to (output)
  //!@}

  //! What the run computes (parameter.general.mode).
  enum CalculationMode {
    ground_state,        //!< imaginary-time evolution to the ground state
    time_evolution,      //!< real-time evolution of an initial state
    finite_temperature,  //!< finite-temperature (density matrix) run
  };
  CalculationMode calcmode;  //!< see CalculationMode

  //! Set the documented default of every parameter.
  PEPS_Parameters();

  //! Append the parameter values to a file (rank 0 only).
  void save(const char *filename, bool append = false);
  //! save() with append = true.
  void save_append(const char *filename) { save(filename, true); }

  //! Broadcast all parameters from root to every rank.
  void Bcast(MPI_Comm comm, int root = 0);

  //! Validate cross-parameter consistency; may throw tenes::input_error.
  void check() const;
};

}  // namespace tenes::itps
#endif  // TENES_SRC_PEPS_PARAMETERS_HPP_
