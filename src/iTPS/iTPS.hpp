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
 *  @brief The solver class ::tenes::itps::iTPS.
 */

#ifndef TENES_SRC_ITPS_ITPS_HPP_
#define TENES_SRC_ITPS_ITPS_HPP_

#include <cstdlib>
#include <complex>
#include <map>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>
#include <set>

#include <optional>

// IWYU pragma begin_exports
#include "../mpi.hpp"
#include "../tensor.hpp"

#include "../fermion/fermion_info.hpp"
#include "../util/type_traits.hpp"

#include "../timer.hpp"
#include "../SquareLattice.hpp"
#include "../operator.hpp"

#include "transfer_matrix.hpp"
#include "correlation_function.hpp"

#include "PEPS_Parameters.hpp"
// IWYU pragma end_exports

namespace tenes::itps {

struct iTPSTestAccessor;

//! A pair of sites identified by the source site and the displacement to
//! the other one; the key type of the two-site measurement results.
struct Bond {
  int source_site;  //!< site index within the unit cell
  int dx;           //!< x displacement to the other site
  int dy;           //!< y displacement to the other site
};
//! Lexicographic order so Bond can key a std::map.
inline bool operator<(const Bond &a, const Bond &b) {
  return std::tie(a.source_site, a.dx, a.dy) <
         std::tie(b.source_site, b.dx, b.dy);
}

//! A cluster of sites identified by the source site and one displacement
//! per further site; the key type of the multi-site measurement results.
struct Multisites {
  int source_site;      //!< site index within the unit cell
  std::vector<int> dx;  //!< x displacements to the other sites
  std::vector<int> dy;  //!< y displacements to the other sites
};
//! Lexicographic order so Multisites can key a std::map.
inline bool operator<(const Multisites &a, const Multisites &b) {
  return std::tie(a.source_site, a.dx, a.dy) <
         std::tie(b.source_site, b.dx, b.dy);
}

/*! @brief The solver: an iTPS wave function (or density matrix), its
 *         environment, and the algorithms acting on them.
 *
 *  One instance holds the whole state of a run: the center tensors, the
 *  CTM and mean-field environments, the evolution operators, and the
 *  observables. main.cpp builds it from input.toml via load_toml.cpp and
 *  then calls optimize() / time_evolution() / finite_temperature()
 *  according to the calculation mode.
 *
 *  The *_density members are the finite-temperature (density matrix /
 *  purification) variants: there the center tensors carry two physical
 *  legs (ket and bra) and the CTM is built from single-layer tensors.
 *
 *  @tparam tensor ::tenes::real_tensor or ::tenes::complex_tensor.
 */
template <class tensor>
class iTPS {
  friend struct iTPSTestAccessor;

 public:
  //! Scalar type of the tensor elements (double or complex).
  using tensor_type = typename tensor::value_type;
  //! True for the real_tensor instantiation.
  static constexpr bool is_tensor_real = std::is_floating_point_v<tensor_type>;

  //! One transfer-matrix result: direction, fixed coordinate, and the
  //! leading eigenvalues.
  using transfer_matrix_eigenvalues_type =
      std::tuple<int, int, std::vector<std::complex<double>>>;

  /*! @brief constructor
   *
   *  @param[in] comm_ communicator all tensors live on
   *  @param[in] peps_parameters_ runtime parameters
   *  @param[in] lattice_ unit-cell geometry
   *  @param[in] simple_updates_ time evolution operators for simple updates
   *  @param[in] full_updates_ time evolution operators for full updates
   *  @param[in] onesite_operators_ onesite operators to be measured
   *  @param[in] twosite_operators_ twosite operators to be measured
   *  @param[in] multisite_operators_ multisite operators to be measured
   *  @param[in] corparam_  parameters for measuring correlation functions
   *  @param[in] tmatrix_param_  parameters for measuring eigenvalues of
   * transfer matrix
   */
  iTPS(MPI_Comm comm_, PEPS_Parameters peps_parameters_, SquareLattice lattice_,
       EvolutionOperators<tensor> simple_updates_,
       EvolutionOperators<tensor> full_updates_,
       Operators<tensor> onesite_operators_,
       Operators<tensor> twosite_operators_,
       Operators<tensor> multisite_operators_, CorrelationParameter corparam_,
       TransferMatrix_Parameters tmatrix_param_);

  //! Allocate and initialize all tensors (random or loaded initial state).
  void initialize_tensors();
  //! Finite-temperature variant of initialize_tensors(): the initial
  //! density matrix is the (infinite-temperature) identity.
  void initialize_tensors_density();
  //! Trace out the bra-side physical legs, yielding one single-layer
  //! tensor per site for the finite-temperature CTM.
  std::vector<tensor> make_single_tensor_density();

  //! Converge the corner transfer matrices for the current state.
  void update_CTM();
  //! Finite-temperature variant of update_CTM().
  void update_CTM_density();

  //! Run all simple-update steps of every operator group.
  void simple_update();
  //! Apply one simple-update gate (one bond or one site).
  void simple_update(EvolutionOperator<tensor> const &up);
  //! Iteratively fix the local gauge of the bonds (simple_update.gauge_fix).
  void fix_local_gauge();
  //! Finite-temperature variant of simple_update().
  void simple_update_density();
  //! Finite-temperature variant of simple_update(up): the gate acts on the
  //! ket side and its conjugate on the bra side.
  void simple_update_density(EvolutionOperator<tensor> const &up);
  //! Purification variant of simple_update_density().
  void simple_update_density_purification();
  //! Purification variant of simple_update_density(up): the state is a
  //! purified wave function, so the gate acts on one layer of the doubled
  //! physical legs only and the other (ancilla) layer is left untouched.
  void simple_update_density_purification(EvolutionOperator<tensor> const &up);
  //! Finite-temperature variant of fix_local_gauge().
  void fix_local_gauge_density();

  //! Run all full-update steps of every operator group (converges the CTM
  //! between gates).
  void full_update();
  //! Apply one full-update gate.
  void full_update(EvolutionOperator<tensor> const &up);

  //! Optimize the state in ground-state mode: simple updates, then full
  //! updates.
  void optimize();
  // void optimize_density();

  /*! @brief Measure and save all observables for the current state.
   *
   *  Converges the environment (CTM or mean field), evaluates the one-,
   *  two-, and multi-site observables, the short-range correlation
   *  functions, and the transfer-matrix correlation length, and writes
   *  the *.dat files into the output directory.
   *
   *  @param[in] time evolution time (or inverse temperature) stamped onto
   *             each output row; no stamp when std::nullopt
   *  @param[in] filename_prefix prefix of the output file names
   */
  void measure(std::optional<double> time = std::nullopt,
               std::string filename_prefix = "");
  //! Finite-temperature variant of measure() at inverse temperature beta.
  void measure_density(double beta, std::string filename_prefix = "FT_");

  //! Print the elapsed-time summary and write time.dat.
  void summary() const;

  //! Run time-evolution mode: apply the real-time gates, measuring every
  //! measure_interval steps.
  void time_evolution();

  //! Run finite-temperature mode: evolve the density matrix in inverse
  //! temperature, measuring every measure_interval steps.
  void finite_temperature();

  //! Measure the one-site observables; result indexed by [operator][site].
  std::vector<std::vector<tensor_type>> measure_onesite();
  //! Finite-temperature variant of measure_onesite().
  std::vector<std::vector<tensor_type>> measure_onesite_density();

  //! Measure the two-site observables; result indexed by [operator][bond].
  std::vector<std::map<Bond, tensor_type>> measure_twosite();
  //! Finite-temperature variant of measure_twosite().
  std::vector<std::map<Bond, tensor_type>> measure_twosite_density();

  //! Measure the multi-site observables; result indexed by
  //! [operator][cluster].
  std::vector<std::map<Multisites, tensor_type>> measure_multisite();
  //! Finite-temperature variant of measure_multisite().
  std::vector<std::map<Multisites, tensor_type>> measure_multisite_density();

  //! Measure the correlation functions up to distance
  //! correlation.r_max along the x and y axes.
  std::vector<Correlation> measure_correlation();

  //! Compute the leading transfer-matrix eigenvalues for every row and
  //! column (the correlation lengths follow from their ratios).
  std::vector<transfer_matrix_eigenvalues_type>
  measure_transfer_matrix_eigenvalues();

  /*! @brief write measured onesite observables
   *
   *  @param[in] onesite_obs
   *  @param[in] time
   *  @param[in] filename_prefix
   */
  void save_onesite(std::vector<std::vector<tensor_type>> const &onesite_obs,
                    std::optional<double> time = std::nullopt,
                    std::string filename_prefix = "");

  /*! @brief write measured twosite observables
   *
   *  @param[in] twosite_obs
   *  @param[in] time
   *  @param[in] filename_prefix
   */
  void save_twosite(std::vector<std::map<Bond, tensor_type>> const &twosite_obs,
                    std::optional<double> time = std::nullopt,
                    std::string filename_prefix = "");

  /*! @brief write measured multisite observables
   *
   *  @param[in] multisite_obs
   *  @param[in] time
   *  @param[in] filename_prefix
   */
  void save_multisite(
      std::vector<std::map<Multisites, tensor_type>> const &multisite_obs,
      std::optional<double> time = std::nullopt,
      std::string filename_prefix = "");

  /*! @brief write measured correlation functions
   *
   *  @param[in] correlations
   *  @param[in] time
   *  @param[in] filename_prefix
   */
  void save_correlation(std::vector<Correlation> const &correlations,
                        std::optional<double> time = std::nullopt,
                        std::string filename_prefix = "");

  /*! @brief calculate and write correlation length
   *
   *  @param[in] eigvals
   *  @param[in] time
   *  @param[in] filename_prefix
   */
  void save_correlation_length(
      std::vector<transfer_matrix_eigenvalues_type> const &eigvals,
      std::optional<double> time = std::nullopt,
      std::string filename_prefix = "");

  /*! @brief write the measured energy and other densities (per-site
   *         averages of the observables)
   *
   *  @param[in] onesite_obs measured one-site observables
   *  @param[in] twosite_obs measured two-site observables
   *  @param[in] multisite_obs measured multi-site observables
   *  @param[in] time time (or inverse temperature) stamped onto each row
   *  @param[in] filename_prefix prefix of the output file name
   */
  void save_density(
      std::vector<std::vector<tensor_type>> const &onesite_obs,
      std::vector<std::map<Bond, tensor_type>> const &twosite_obs,
      std::vector<std::map<Multisites, tensor_type>> const &multisite_obs,
      std::optional<double> time = std::nullopt,
      std::string filename_prefix = "");

  //! save optimized tensors into files
  void save_tensors() const;

  //! load tensors from files
  void load_tensors();

 private:
  //! Reject measurement requests unsupported by the fermionic CTM path.
  void validate_fermion_ctm_measurement() const;

  //! Index of the group-th one-site operator acting on site.
  int siteoperator_index(int site, int group) const {
    return site_ops_indices[site][group];
  }

  //! Convert a scalar to tensor_type, dropping the imaginary part for
  //! real runs.
  template <class T>
  tensor_type to_tensor_type(T const &v) const {
    return convert_complex<tensor_type>(v);
  }

  //! Read tensors in the current save format, whose params.dat records
  //! the format version, unit cell, and per-site shapes.
  void load_tensors_v1();
  //! Read tensors in the legacy save format (bare tensor files, no
  //! params.dat metadata).
  void load_tensors_v0();
  /*!
   * @brief Save the fermionic parity ledger to fermion.dat (rank 0 only).
   *
   * The virtual-bond ledger is mutable state (the simple update rewrites
   * it through the graded svd_trunc), so it must travel with the saved
   * tensors: reloading them under a stale ledger silently changes the
   * measured energy. The file records a format version, the unit-cell
   * geometry, and the physical and virtual parity vectors of every site.
   */
  void save_fermion_parity(std::string const &save_dir) const;
  /*!
   * @brief Load and validate fermion.dat before reading the tensors.
   *
   * Restores finfo from the saved ledger. Throws tenes::load_error when
   * the saved run and the current one disagree about fermion mode
   * (fermion.dat present but fermion = false, or vice versa), or when the
   * file does not match the current geometry.
   */
  void load_fermion_ledger(std::string const &load_dir);
  //! check the loaded tensors against the restored ledger (after reading them)
  void validate_loaded_fermion_tensors() const;

  //! measure_correlation() with the CTM environment.
  std::vector<Correlation> measure_correlation_ctm();
  //! measure_correlation() with the mean-field environment.
  std::vector<Correlation> measure_correlation_mf();

  static constexpr int nleg = 4;  //!< virtual legs per center tensor

  MPI_Comm comm;  //!< communicator all tensors live on
  int mpisize;    //!< size of comm
  int mpirank;    //!< rank within comm

  PEPS_Parameters peps_parameters;  //!< runtime parameters
  SquareLattice lattice;            //!< unit-cell geometry
  //! Physical-leg parities copied from PEPS_Parameters::phys_parity at
  //! initialize_tensors(); the immutable half of what seeds finfo.
  std::vector<std::vector<bool>> phys_parity;

  //! Simple-update gates, all groups interleaved.
  EvolutionOperators<tensor> simple_updates;
  //! Number of simple-update operator groups.
  int num_simple_update_groups;
  //! Full-update gates, all groups interleaved.
  EvolutionOperators<tensor> full_updates;
  //! Number of full-update operator groups.
  int num_full_update_groups;
  Operators<tensor> onesite_operators;    //!< observables on one site
  Operators<tensor> twosite_operators;    //!< observables on two sites
  Operators<tensor> multisite_operators;  //!< observables on 3+ sites
  //! site_ops_indices[site][group]: index into onesite_operators.
  std::vector<std::vector<int>> site_ops_indices;
  int num_onesite_operators;    //!< number of one-site operator groups
  int num_twosite_operators;    //!< number of two-site operator groups
  int num_multisite_operators;  //!< number of multi-site operator groups
  //! Name of each one-site operator group (for the output files).
  std::vector<std::string> onesite_operator_names;
  //! Name of each two-site operator group.
  std::vector<std::string> twosite_operator_names;
  //! Name of each multi-site operator group.
  std::vector<std::string> multisite_operator_names;
  //! Number of terms in each one-site operator group.
  std::vector<int> onesite_operator_counts;
  //! Number of terms in each two-site operator group.
  std::vector<int> twosite_operator_counts;
  //! Number of terms in each multi-site operator group.
  std::vector<int> multisite_operator_counts;
  //! Number of sites each multi-site operator group acts on.
  std::vector<int> multisite_operator_nsites;
  //! Distinct multi-site cluster sizes present.
  std::set<int> multisite_operator_nsites_set;

  //! Identity operator on the physical leg, per site.
  std::vector<tensor> op_identity;

  CorrelationParameter corparam;            //!< correlation-function settings
  TransferMatrix_Parameters tmatrix_param;  //!< correlation-length settings

  /*! @name Tensors
   *  @brief Tensors of an iTPS
   *
   *  An index of vector is corresponding with the index of the site in the
   * unitcell.
   *
   *  The ordering of tensors is as following:
   *
   *  C1 eTt C2
   *
   *  eTl Tn eTr
   *
   *  C4 eTb C3
   */
  //!@{

  /*! @brief Center tensors
   *
   *  Tn[i] has 4 virtual legs and 1 physical leg;
   *  Tn[i][left, up, right, bottom, physical]
   *
   */
  std::vector<tensor> Tn;
  std::vector<tensor> eTt;  //!< Top edge tensor for each center
  std::vector<tensor> eTr;  //!< Right edge tensor for each center
  std::vector<tensor> eTb;  //!< Bottom edge tensor for each center
  std::vector<tensor> eTl;  //!< Left edge tensor for each center
  std::vector<tensor> C1;   //!< Left-top CTM for each center
  std::vector<tensor> C2;   //!< Right-top CTM for each center
  std::vector<tensor> C3;   //!< Right-bottom CTM for each center
  std::vector<tensor> C4;   //!< Left-bottom CTM for each center
  //!@}
  //! Parity ledgers of every site tensor (fermion mode). Mutable state:
  //! the simple update rewrites the virtual-bond ledgers through the
  //! graded svd_trunc, so it is saved to and loaded from fermion.dat
  //! alongside the tensors. Disabled (enabled = false) in bosonic runs.
  tenes::fermion::FermionInfo finfo;
  std::vector<std::vector<std::vector<double>>>
      lambda_tensor;  //!< Meanfield environments

  std::size_t CHI;  //!< Bond dimension of corner transfer matrices
  int LX;           //!< Length of a unitcell along with X axes
  int LY;           //!< Length of a unitcell along with Y axes
  int N_UNIT;       //!< The number of sites in a unitcell

  std::string outdir;  //!< path to the directory where results will be written

  Timer<> timer_all;          //!< wall clock of the whole run
  double time_simple_update;  //!< seconds spent in simple updates
  double time_full_update;    //!< seconds spent in full updates
  double time_environment;    //!< seconds spent converging environments
  double time_observable;     //!< seconds spent measuring observables
};

}  // namespace tenes::itps

#endif  //  TENES_SRC_ITPS_ITPS_HPP_
