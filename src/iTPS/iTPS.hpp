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

#include <boost/optional.hpp>

// IWYU pragma begin_exports
#include "../mpi.hpp"
#include "../tensor.hpp"

#include "../util/type_traits.hpp"

#include "../timer.hpp"
#include "../SquareLattice.hpp"
#include "../operator.hpp"

#include "transfer_matrix.hpp"
#include "correlation_function.hpp"

#include "PEPS_Parameters.hpp"
// IWYU pragma end_exports

namespace tenes {
namespace itps {

struct Bond {
  int source_site;
  int dx;
  int dy;
};
// std::map<Bond, T> requires operator<
inline bool operator<(const Bond &a, const Bond &b) {
  return std::tie(a.source_site, a.dx, a.dy) <
         std::tie(b.source_site, b.dx, b.dy);
}

struct Multisites {
  int source_site;
  std::vector<int> dx;
  std::vector<int> dy;
};
// std::map<Multisites, T> requires operator<
inline bool operator<(const Multisites &a, const Multisites &b) {
  return std::tie(a.source_site, a.dx, a.dy) <
         std::tie(b.source_site, b.dx, b.dy);
}

//! Solver main class
template <class tensor>
class iTPS {
 public:
  using tensor_type = typename tensor::value_type;
  static constexpr bool is_tensor_real =
      std::is_floating_point<tensor_type>::value;

  using transfer_matrix_eigenvalues_type =
      std::tuple<int, int, std::vector<std::complex<double>>>;

  /*! @brief constructor
   *
   *  @param[in] comm_
   *  @param[in] peps_parameters
   *  @param[in] lattice_
   *  @param[in] simple_updates Time Evolution operators for simple updates
   *  @param[in] full_updates Time Evolution operators for full updates
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

  //! initialize tensors
  void initialize_tensors();
  void initialize_tensors_density();
  std::vector<tensor> make_single_tensor_density();
  
  //! update corner transfer matrices
  void update_CTM();
  void update_CTM_density();

  //! perform simple update
  void simple_update();
  void simple_update(EvolutionOperator<tensor> const &up);
  void fix_local_gauge();
  void simple_update_density();
  void simple_update_density(EvolutionOperator<tensor> const &up);
  void simple_update_density_purification();
  void simple_update_density_purification(EvolutionOperator<tensor> const &up);
  void fix_local_gauge_density();

  //! perform full update
  void full_update();
  void full_update(EvolutionOperator<tensor> const &up);

  //! optimize tensors
  void optimize();
  //void optimize_density();

  //! measure expectation value of observables
  void measure(boost::optional<double> time = boost::optional<double>(),
               std::string filename_prefix = "");
  void measure_density(double beta, std::string filename_prefix = "FT_");

  //! print elapsed time
  void summary() const;

  void time_evolution();

  void finite_temperature();

  //! measure expectation value of onesite observables
  std::vector<std::vector<tensor_type>> measure_onesite();
  std::vector<std::vector<tensor_type>> measure_onesite_density();

  //! measure expectation value of twosite observables
  std::vector<std::map<Bond, tensor_type>> measure_twosite();
  std::vector<std::map<Bond, tensor_type>> measure_twosite_density();

  //! measure expectation value of multisite observables
  std::vector<std::map<Multisites, tensor_type>> measure_multisite();
  std::vector<std::map<Multisites, tensor_type>> measure_multisite_density();

  //! measure correlation functions
  std::vector<Correlation> measure_correlation();

  //! measure eigenvalues of transfer matrix
  std::vector<transfer_matrix_eigenvalues_type>
  measure_transfer_matrix_eigenvalues();

  /*! @brief write measured onesite observables
   *
   *  @param[in] onesite_obs
   *  @param[in] time
   *  @param[in] filename_prefix
   */
  void save_onesite(std::vector<std::vector<tensor_type>> const &onesite_obs,
                    boost::optional<double> time = boost::optional<double>(),
                    std::string filename_prefix = "");

  /*! @brief write measured twosite observables
   *
   *  @param[in] twosite_obs
   *  @param[in] time
   *  @param[in] filename_prefix
   */
  void save_twosite(std::vector<std::map<Bond, tensor_type>> const &twosite_obs,
                    boost::optional<double> time = boost::optional<double>(),
                    std::string filename_prefix = "");

  /*! @brief write measured multisite observables
   *
   *  @param[in] multisite_obs
   *  @param[in] time
   *  @param[in] filename_prefix
   */
  void save_multisite(
      std::vector<std::map<Multisites, tensor_type>> const &multisite_obs,
      boost::optional<double> time = boost::optional<double>(),
      std::string filename_prefix = "");

  /*! @brief write measured correlation functions
   *
   *  @param[in] correlations
   *  @param[in] time
   *  @param[in] filename_prefix
   */
  void save_correlation(
      std::vector<Correlation> const &correlations,
      boost::optional<double> time = boost::optional<double>(),
      std::string filename_prefix = "");

  /*! @brief calculate and write correlation length
   *
   *  @param[in] eigvals
   *  @param[in] time
   *  @param[in] filename_prefix
   */
  void save_correlation_length(
      std::vector<transfer_matrix_eigenvalues_type> const &eigvals,
      boost::optional<double> time = boost::optional<double>(),
      std::string filename_prefix = "");

  /*! @brief calculate and write correlation length
   *
   *  @param[in] onesite_obs
   *  @param[in] twosite_obs,
   *  @param[in] multisite_obs,
   *  @param[in] time
   *  @param[in] filename_prefix
   */
  void save_density(std::vector<std::vector<tensor_type>> const &onesite_obs,
                    std::vector<std::map<Bond, tensor_type>> const &twosite_obs,
                    std::vector<std::map<Multisites, tensor_type>> const &multisite_obs,
                    boost::optional<double> time = boost::optional<double>(),
                    std::string filename_prefix = "");

  //! save optimized tensors into files
  void save_tensors() const;

  //! load tensors from files
  void load_tensors();

 private:
  int siteoperator_index(int site, int group) const {
    return site_ops_indices[site][group];
  }

  template <class T>
  tensor_type to_tensor_type(T const &v) const {
    return convert_complex<tensor_type>(v);
  }

  void load_tensors_v1();
  void load_tensors_v0();

  std::vector<Correlation> measure_correlation_ctm();
  std::vector<Correlation> measure_correlation_mf();

  static constexpr int nleg = 4;

  MPI_Comm comm;
  int mpisize, mpirank;

  PEPS_Parameters peps_parameters;
  SquareLattice lattice;

  EvolutionOperators<tensor> simple_updates;
  int num_simple_update_groups;
  EvolutionOperators<tensor> full_updates;
  int num_full_update_groups;
  Operators<tensor> onesite_operators;
  Operators<tensor> twosite_operators;
  Operators<tensor> multisite_operators;
  std::vector<std::vector<int>> site_ops_indices;
  int num_onesite_operators;
  int num_twosite_operators;
  int num_multisite_operators;
  std::vector<std::string> onesite_operator_names;
  std::vector<std::string> twosite_operator_names;
  std::vector<std::string> multisite_operator_names;
  std::vector<int> onesite_operator_counts;
  std::vector<int> twosite_operator_counts;
  std::vector<int> multisite_operator_counts;
  std::vector<int> multisite_operator_nsites;
  std::set<int> multisite_operator_nsites_set;

  std::vector<tensor> op_identity;

  CorrelationParameter corparam;
  TransferMatrix_Parameters tmatrix_param;

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
  std::vector<std::vector<std::vector<double>>>
      lambda_tensor;  //!< Meanfield environments

  std::size_t CHI;  //!< Bond dimension of corner transfer matrices
  int LX;           //!< Length of a unitcell along with X axes
  int LY;           //!< Length of a unitcell along with Y axes
  int N_UNIT;       //!< The number of sites in a unitcell

  std::string outdir;  //!< path to the directory where results will be written

  Timer<> timer_all;
  double time_simple_update;
  double time_full_update;
  double time_environment;
  double time_observable;
};

}  // namespace itps
}  // namespace tenes

#endif  //  TENES_SRC_ITPS_ITPS_HPP_
