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

#include <complex>
#include <ctime>
#include <map>
#include <string>
#include <tuple>
#include <type_traits>
#include <array>
#include <functional>
#include <vector>

#include "../tensor.hpp"

#include "../Lattice.hpp"
#include "../operator.hpp"
#include "../PEPS_Parameters.hpp"
#include "../correlation.hpp"
#include "../printlevel.hpp"
#include "../timer.hpp"
#include "../util/file.hpp"
#include "../util/string.hpp"
#include "../util/type_traits.hpp"
#include "../correlation_length.hpp"

namespace tenes {

struct Bond {
  int source_site;
  int dx;
  int dy;
};

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
   *  @param[in] simple_updates ITE operators for simple updates
   *  @param[in] full_updates ITE operators for full updates
   *  @param[in] onesite_operators
   *  @param[in] twosite_operators
   *  @param[in] corparam_
   */
  iTPS(MPI_Comm comm_, PEPS_Parameters peps_parameters_, Lattice lattice_,
       NNOperators<tensor> simple_updates_, NNOperators<tensor> full_updates_,
       Operators<tensor> onesite_operators,
       Operators<tensor> twosite_operators, CorrelationParameter corparam_,
       CorrelationLengthCalculator_Parameters clength_param_);

  void initialize_tensors();
  void update_CTM();
  void simple_update();
  void full_update();

  //! optimize tensors
  void optimize();

  //! measure expectation value of observables
  void measure();

  //! print elapsed time
  void summary() const;

  //! measure expectation value of onesite observables
  std::vector<std::vector<tensor_type>> measure_onesite();

  //! measure expectation value of twosite observables
  std::vector<std::map<Bond, tensor_type>> measure_twosite();

  //! measure correlation functions
  std::vector<Correlation> measure_correlation();

  std::vector<transfer_matrix_eigenvalues_type>
  measure_transfer_matrix_eigenvalues();

  /*! @brief write measured onesite observables
   *
   *  @params[in] onesite_obs
   */
  void save_onesite(std::vector<std::vector<tensor_type>> const &onesite_obs);

  /*! @brief write measured twosite observables
   *
   *  @params[in] twosite_obs
   */
  void save_twosite(
      std::vector<std::map<Bond, tensor_type>> const &twosite_obs);

  /*! @brief write measured correlation functions
   *
   *  @params[in] correlations
   */
  void save_correlation(std::vector<Correlation> const &correlations);

  /*! @brief write measured correlation length
   *
   *  @params[in] xi
   */
  void save_correlation_length(
      std::vector<transfer_matrix_eigenvalues_type> const &xi);

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
  Lattice lattice;

  NNOperators<tensor> simple_updates;
  NNOperators<tensor> full_updates;
  Operators<tensor> onesite_operators;
  Operators<tensor> twosite_operators;
  std::vector<std::vector<int>> site_ops_indices;
  int num_onesite_operators;
  int num_twosite_operators;
  std::vector<std::string> onesite_operator_names;
  std::vector<std::string> twosite_operator_names;

  std::vector<tensor> op_identity;

  CorrelationParameter corparam;
  CorrelationLengthCalculator_Parameters clength_param;

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

  size_t CHI;  //!< Bond dimension of corner transfer matrices
  int LX;      //!< Length of a unitcell along with X axes
  int LY;      //!< Length of a unitcell along with Y axes
  int N_UNIT;  //!< The number of sites in a unitcell

  std::string outdir;  //!< path to the directory where results will be written

  Timer<> timer_all;
  double time_simple_update;
  double time_full_update;
  double time_environment;
  double time_observable;
};

}  // end of namespace tenes

#endif  //  TENES_SRC_ITPS_ITPS_HPP_
