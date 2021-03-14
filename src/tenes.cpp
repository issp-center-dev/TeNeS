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

#define _USE_MATH_DEFINES
#include <sys/stat.h>
#include <algorithm>
#include <complex>
#include <ctime>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <type_traits>
#include <array>
#include <functional>

#ifdef _NO_OMP
int omp_get_max_threads() { return 1; }
#else
#include <omp.h>
#endif

#include <mptensor/complex.hpp>
#include <mptensor/rsvd.hpp>
#include <mptensor/tensor.hpp>

#include "tensor.hpp"

#include "Lattice.hpp"
#include "PEPS_Basics.hpp"
#include "PEPS_Parameters.hpp"
#include "Square_lattice_CTM.hpp"
#include "correlation.hpp"
#include "printlevel.hpp"
#include "timer.hpp"
#include "util/file.hpp"
#include "util/string.hpp"
#include "util/type_traits.hpp"

#include "tenes.hpp"

namespace tenes {
// using namespace mptensor;

struct Bond {
  int source_site;
  int dx;
  int dy;
};

bool operator<(const Bond &a, const Bond &b) {
  return std::tie(a.source_site, a.dx, a.dy) <
         std::tie(b.source_site, b.dx, b.dy);
}

std::string datetime(std::time_t const &t) {
  std::tm *lt = std::localtime(&t);
  char dt[64];
  std::strftime(dt, 64, "%FT%T", lt);
  char tz[7];
  std::strftime(tz, 7, "%z", lt);
  if (tz[0] != 'Z' && tz[0] != 'z') {
    for (size_t i = 6; i > 3; --i) {
      tz[i] = tz[i - 1];
    }
    tz[3] = ':';
  }
  std::string ret = dt;
  return ret + tz;
}
std::string datetime() {
  std::time_t t = std::time(nullptr);
  return datetime(t);
}

//! Solver main class
template <class ptensor>
class TeNeS {
 public:
  using tensor_type = typename ptensor::value_type;
  static constexpr bool is_tensor_real =
      std::is_floating_point<tensor_type>::value;

  using transfer_matrix_eigenvalues_type =
      std::tuple<int, int, double, double, double>;

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
  TeNeS(MPI_Comm comm_, PEPS_Parameters peps_parameters_, Lattice lattice_,
        NNOperators<ptensor> simple_updates_,
        NNOperators<ptensor> full_updates_,
        Operators<ptensor> onesite_operators,
        Operators<ptensor> twosite_operators, CorrelationParameter corparam_);

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
  void transfer_matvec_horizontal(ptensor &vec, int y);
  void transfer_matvec_vertical(ptensor &vec, int x);
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

  NNOperators<ptensor> simple_updates;
  NNOperators<ptensor> full_updates;
  Operators<ptensor> onesite_operators;
  Operators<ptensor> twosite_operators;
  std::vector<std::vector<int>> site_ops_indices;
  int num_onesite_operators;
  int num_twosite_operators;
  std::vector<std::string> onesite_operator_names;
  std::vector<std::string> twosite_operator_names;

  std::vector<ptensor> op_identity;

  CorrelationParameter corparam;

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
  std::vector<ptensor> Tn;
  std::vector<ptensor> eTt;  //!< Top edge tensor for each center
  std::vector<ptensor> eTr;  //!< Right edge tensor for each center
  std::vector<ptensor> eTb;  //!< Bottom edge tensor for each center
  std::vector<ptensor> eTl;  //!< Left edge tensor for each center
  std::vector<ptensor> C1;   //!< Left-top CTM for each center
  std::vector<ptensor> C2;   //!< Right-top CTM for each center
  std::vector<ptensor> C3;   //!< Right-bottom CTM for each center
  std::vector<ptensor> C4;   //!< Left-bottom CTM for each center
  //!@}
  std::vector<std::vector<std::vector<double>>>
      lambda_tensor;  //!< Meanfield environments

  int CHI;     //!< Bond dimension of corner transfer matrices
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

template <class ptensor>
TeNeS<ptensor>::TeNeS(MPI_Comm comm_, PEPS_Parameters peps_parameters_,
                      Lattice lattice_, NNOperators<ptensor> simple_updates_,
                      NNOperators<ptensor> full_updates_,
                      Operators<ptensor> onesite_operators_,
                      Operators<ptensor> twosite_operators_,
                      CorrelationParameter corparam_)
    : comm(comm_),
      peps_parameters(peps_parameters_),
      lattice(lattice_),
      simple_updates(simple_updates_),
      full_updates(full_updates_),
      onesite_operators(onesite_operators_),
      twosite_operators(twosite_operators_),
      corparam(corparam_),
      outdir("output"),
      timer_all(),
      time_simple_update(),
      time_full_update(),
      time_environment(),
      time_observable() {
  MPI_Comm_size(comm, &mpisize);
  MPI_Comm_rank(comm, &mpirank);

  peps_parameters.check();

  peps_parameters.Bcast(comm);

  // output debug or warning info only from process 0
  if (mpirank != 0) {
    peps_parameters.print_level = PrintLevel::none;
  }

  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "Number of Processes: " << mpisize << std::endl;
    std::cout << "Number of Threads / Process: " << omp_get_max_threads()
              << std::endl;

    if (peps_parameters.is_real) {
      std::cout << "Tensor type: real" << std::endl;
    } else {
      std::cout << "Tensor type: complex" << std::endl;
    }
  }

  CHI = peps_parameters.CHI;

  lattice.Bcast(comm);

  LX = lattice.LX;
  LY = lattice.LY;
  N_UNIT = lattice.N_UNIT;

  // set seed for randomized svd
  int seed = peps_parameters.seed;
  random_tensor::set_seed(seed + mpirank);

  outdir = peps_parameters.outdir;

  bool is_ok = true;

  if (mpirank == 0) {
    if (!util::isdir(outdir)) {
      is_ok = is_ok && util::mkdir(outdir);
    }

    std::string param_file = outdir + "/parameters.dat";

    peps_parameters.save(param_file.c_str());
    lattice.save_append(param_file.c_str());

    std::fstream ofs(param_file, std::ios::out | std::ios::app);
    ofs << std::endl;
    ofs << "start_datetime =  " << datetime() << std::endl;
  }
  bcast(is_ok, 0, comm);
  if (!is_ok) {
    std::stringstream ss;
    ss << "Cannot mkdir " << outdir;
    throw tenes::runtime_error(ss.str());
  }

  std::string savedir = peps_parameters_.tensor_save_dir;
  if (mpirank == 0) {
    if (!savedir.empty()) {
      if (!util::isdir(savedir)) {
        is_ok = is_ok && util::mkdir(savedir);
      }
    }
  }
  bcast(is_ok, 0, comm);
  if (!is_ok) {
    std::stringstream ss;
    ss << "Cannot mkdir " << savedir;
    throw tenes::runtime_error(ss.str());
  }

  initialize_tensors();

  int maxops = -1;
  size_t maxlength = 0;
  for (auto const &op : onesite_operators) {
    maxops = std::max(op.group, maxops);
  }
  num_onesite_operators = maxops + 1;
  onesite_operator_names.resize(num_onesite_operators);
  for (auto const &op : onesite_operators) {
    if (op.name.empty()) {
      std::stringstream ss;
      ss << "onesite[" << op.group << "]";
      onesite_operator_names[op.group] = ss.str();
    } else {
      onesite_operator_names[op.group] = op.name;
    }
  }
  for (auto const &s : onesite_operator_names) {
    maxlength = std::max(s.size(), maxlength);
  }

  maxops = -1;
  for (auto const &op : twosite_operators) {
    maxops = std::max(op.group, maxops);
  }
  num_twosite_operators = maxops + 1;
  twosite_operator_names.resize(num_twosite_operators);
  for (auto const &op : twosite_operators) {
    if (op.name.empty()) {
      std::stringstream ss;
      ss << "twosite[" << op.group << "]";
      twosite_operator_names[op.group] = ss.str();
    } else {
      twosite_operator_names[op.group] = op.name;
    }
  }
  for (auto const &s : twosite_operator_names) {
    maxlength = std::max(s.size(), maxlength);
  }

  for (auto &s : onesite_operator_names) {
    const auto l = maxlength - s.size();
    for (size_t i = 0; i < l; ++i) {
      s += " ";
    }
  }
  for (auto &s : twosite_operator_names) {
    const auto l = maxlength - s.size();
    for (size_t i = 0; i < l; ++i) {
      s += " ";
    }
  }

  site_ops_indices.resize(N_UNIT, std::vector<int>(num_onesite_operators, -1));
  for (int i = 0; i < onesite_operators.size(); ++i) {
    auto const &op = onesite_operators[i];
    site_ops_indices[op.source_site][op.group] = i;
  }
}

template <class ptensor>
void TeNeS<ptensor>::initialize_tensors() {
  Tn.clear();
  eTt.clear();
  eTr.clear();
  eTb.clear();
  eTl.clear();
  C1.clear();
  C2.clear();
  C3.clear();
  C4.clear();
  lambda_tensor.clear();

  for (int i = 0; i < N_UNIT; ++i) {
    const auto pdim = lattice.physical_dims[i];
    const auto vdim = lattice.virtual_dims[i];

    Tn.push_back(ptensor(Shape(vdim[0], vdim[1], vdim[2], vdim[3], pdim)));
    eTt.push_back(ptensor(Shape(CHI, CHI, vdim[1], vdim[1])));
    eTr.push_back(ptensor(Shape(CHI, CHI, vdim[2], vdim[2])));
    eTb.push_back(ptensor(Shape(CHI, CHI, vdim[3], vdim[3])));
    eTl.push_back(ptensor(Shape(CHI, CHI, vdim[0], vdim[0])));
    C1.push_back(ptensor(Shape(CHI, CHI)));
    C2.push_back(ptensor(Shape(CHI, CHI)));
    C3.push_back(ptensor(Shape(CHI, CHI)));
    C4.push_back(ptensor(Shape(CHI, CHI)));

    std::vector<std::vector<double>> lambda(nleg);
    for (int j = 0; j < nleg; ++j) {
      lambda[j] = std::vector<double>(vdim[j], 1.0);
    }
    lambda_tensor.push_back(lambda);

    ptensor id(mptensor::Shape(pdim, pdim));
    for (int j = 0; j < pdim; ++j) {
      for (int k = 0; k < pdim; ++k) {
        id.set_value(mptensor::Index(j, k), (j == k ? 1.0 : 0.0));
      }
    }
    op_identity.push_back(id);
  }

  std::mt19937 gen(peps_parameters.seed);
  // use another rng for backward compatibility
  std::mt19937 gen_im(peps_parameters.seed * 11 + 137);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  int nr;

  if (peps_parameters.tensor_load_dir.empty()) {
    mptensor::Index index;
    for (int i = 0; i < lattice.N_UNIT; ++i) {
      const auto pdim = lattice.physical_dims[i];
      const auto vdim = lattice.virtual_dims[i];

      const size_t ndim = vdim[0] * vdim[1] * vdim[2] * vdim[3] * pdim;
      std::vector<double> ran_re(ndim);
      std::vector<double> ran_im(ndim);

      for (int j = 0; j < ndim; j++) {
        ran_re[j] = dist(gen);
        ran_im[j] = dist(gen_im);
      }
      auto &dir = lattice.initial_dirs[i];
      std::vector<double> dir_im(pdim);
      if (std::all_of(dir.begin(), dir.end(),
                      [=](double x) { return x == 0.0; })) {
        // random
        dir.resize(pdim);
        for (int j = 0; j < pdim; ++j) {
          dir[j] = dist(gen);
          dir_im[j] = dist(gen_im);
        }
      }

      for (int n = 0; n < Tn[i].local_size(); ++n) {
        index = Tn[i].global_index(n);
        if (index[0] == 0 && index[1] == 0 && index[2] == 0 && index[3] == 0) {
          auto v = std::complex<double>(dir[index[4]], dir_im[index[4]]);
          Tn[i].set_value(index, to_tensor_type(v));
        } else {
          nr = index[0] + index[1] * vdim[0] + index[2] * vdim[0] * vdim[1] +
               index[3] * vdim[0] * vdim[1] * vdim[2] +
               index[4] * vdim[0] * vdim[1] * vdim[2] * vdim[3];
          auto v =
              lattice.noises[i] * std::complex<double>(ran_re[nr], ran_im[nr]);
          Tn[i].set_value(index, to_tensor_type(v));
        }
      }
    }
  } else {
    load_tensors();
    if (peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "Tensors loaded from " << peps_parameters.tensor_load_dir
                << std::endl;
    }
  }  // end of else part of if(load_dir.empty())
}

template <class ptensor>
inline void TeNeS<ptensor>::update_CTM() {
  Timer<> timer;
  Calc_CTM_Environment(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, peps_parameters,
                       lattice);
  time_environment += timer.elapsed();
}

template <class ptensor>
void TeNeS<ptensor>::simple_update() {
  Timer<> timer;
  ptensor Tn1_new;
  ptensor Tn2_new;
  std::vector<double> lambda_c;
  const int nsteps = peps_parameters.num_simple_step;
  double next_report = 10.0;

  for (int int_tau = 0; int_tau < nsteps; ++int_tau) {
    for (auto up : simple_updates) {
      const int source = up.source_site;
      const int source_leg = up.source_leg;
      const int target = lattice.neighbor(source, source_leg);
      const int target_leg = (source_leg + 2) % 4;
      Simple_update_bond(Tn[source], Tn[target], lambda_tensor[source],
                         lambda_tensor[target], up.op, source_leg,
                         peps_parameters, Tn1_new, Tn2_new, lambda_c);
      lambda_tensor[source][source_leg] = lambda_c;
      lambda_tensor[target][target_leg] = lambda_c;
      Tn[source] = Tn1_new;
      Tn[target] = Tn2_new;
    }

    if (peps_parameters.print_level >= PrintLevel::info) {
      double r_tau = 100.0 * (int_tau + 1) / nsteps;
      if (r_tau >= next_report) {
        while (r_tau >= next_report) {
          next_report += 10.0;
        }
        std::cout << "  " << next_report - 10.0
                  << "% "
                     "["
                  << int_tau + 1 << "/" << nsteps << "] done" << std::endl;
      }
    }
  }
  time_simple_update += timer.elapsed();
}

template <class ptensor>
void TeNeS<ptensor>::full_update() {
  if (peps_parameters.num_full_step > 0) {
    update_CTM();
  }

  Timer<> timer;
  ptensor Tn1_new, Tn2_new;
  const int nsteps = peps_parameters.num_full_step;
  double next_report = 10.0;

  for (int int_tau = 0; int_tau < nsteps; ++int_tau) {
    for (auto up : full_updates) {
      const int source = up.source_site;
      const int source_leg = up.source_leg;
      const int target = lattice.neighbor(source, source_leg);
      // const int target_leg = (source_leg + 2) % 4;

      if (source_leg == 0) {
        /*
         *  C1' t' t C3
         *  l'  T' T r
         *  C2' b' b C4
         *
         *   |
         *   | rotate
         *   V
         *
         *  C4 b b' C2'
         *  r  T T' l'
         *  C3 t t' C1'
         */
        Full_update_bond(C4[source], C2[target], C1[target], C3[source],
                         eTb[source], eTb[target], eTl[target], eTt[target],
                         eTt[source], eTr[source], Tn[source], Tn[target],
                         up.op, source_leg, peps_parameters, Tn1_new, Tn2_new);
      } else if (source_leg == 1) {
        /*
         * C1' t' C2'
         *  l' T'  r'
         *  l  T   r
         * C4  b  C3
         *
         *   |
         *   | rotate
         *   V
         *
         *  C4 l l' C1'
         *  b  T T' t'
         *  C3 r r' C2'
         */
        Full_update_bond(C4[source], C1[target], C2[target], C3[source],
                         eTl[source], eTl[target], eTt[target], eTr[target],
                         eTr[source], eTb[source], Tn[source], Tn[target],
                         up.op, source_leg, peps_parameters, Tn1_new, Tn2_new);
      } else if (source_leg == 2) {
        /*
         *  C1 t t' C2'
         *  l  T T' r'
         *  C4 b b' C3'
         */
        Full_update_bond(C1[source], C2[target], C3[target], C4[source],
                         eTt[source], eTt[target], eTr[target],  // t  t' r'
                         eTb[target], eTb[source], eTl[source],  // b' b  l
                         Tn[source], Tn[target], up.op, source_leg,
                         peps_parameters, Tn1_new, Tn2_new);
      } else {
        /*
         * C1  t C2
         *  l  T  r
         *  l' T' r'
         * C4' b C3'
         *
         *   |
         *   | rotate
         *   V
         *
         *  C2 r r' C3'
         *  t  T T' b'
         *  C1 l l' C4'
         */
        Full_update_bond(C2[source], C3[target], C4[target], C1[source],
                         eTr[source], eTr[target], eTb[target], eTl[target],
                         eTl[source], eTt[source], Tn[source], Tn[target],
                         up.op, source_leg, peps_parameters, Tn1_new, Tn2_new);
      }
      Tn[source] = Tn1_new;
      Tn[target] = Tn2_new;

      if (peps_parameters.Full_Use_FastFullUpdate) {
        if (source_leg == 0) {
          const int source_x = source % LX;
          const int target_x = target % LX;
          Right_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, source_x,
                     peps_parameters, lattice);
          Left_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, target_x,
                    peps_parameters, lattice);
        } else if (source_leg == 1) {
          const int source_y = source / LX;
          const int target_y = target / LX;
          Bottom_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, source_y,
                      peps_parameters, lattice);
          Top_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, target_y,
                   peps_parameters, lattice);
        } else if (source_leg == 2) {
          const int source_x = source % LX;
          const int target_x = target % LX;
          Left_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, source_x,
                    peps_parameters, lattice);
          Right_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, target_x,
                     peps_parameters, lattice);
        } else {
          const int source_y = source / LX;
          const int target_y = target / LX;
          Top_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, source_y,
                   peps_parameters, lattice);
          Bottom_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, target_y,
                      peps_parameters, lattice);
        }
      } else {
        update_CTM();
      }
    }

    if (peps_parameters.print_level >= PrintLevel::info) {
      double r_tau = 100.0 * (int_tau + 1) / nsteps;
      if (r_tau >= next_report) {
        while (r_tau >= next_report) {
          next_report += 10.0;
        }
        std::cout << "  " << next_report - 10.0
                  << "% "
                     "["
                  << int_tau + 1 << "/" << nsteps << "] done" << std::endl;
      }
    }
  }
  time_full_update += timer.elapsed();
}

template <class ptensor>
void TeNeS<ptensor>::optimize() {
  // for measure time
  std::cout << std::setprecision(12);

  ptensor Tn1_new, Tn2_new;
  std::vector<double> lambda_c;

  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "Start simple update" << std::endl;
  }
  simple_update();

  if (peps_parameters.num_full_step > 0) {
    if (peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "Start full update" << std::endl;
    }
    full_update();
  }
}

template <class ptensor>
auto TeNeS<ptensor>::measure_onesite()
    -> std::vector<std::vector<typename TeNeS<ptensor>::tensor_type>> {
  Timer<> timer;
  const int nlops = num_onesite_operators;
  std::vector<std::vector<tensor_type>> local_obs(
      nlops, std::vector<tensor_type>(
                 N_UNIT, std::numeric_limits<double>::quiet_NaN()));

  if (peps_parameters.MeanField_Env) {
    std::vector<ptensor> Tn_(Tn);
    for (int i = 0; i < N_UNIT; ++i) {
      for (int leg = 0; leg < nleg; ++leg) {
        const std::vector<double> mf = lambda_tensor[i][leg];
        Tn_[i].multiply_vector(mf, leg);
      }
    }

    std::vector<double> norm(N_UNIT);
    for (int i = 0; i < N_UNIT; ++i) {
      const auto n = Contract_one_site_MF(Tn_[i], op_identity[i]);
      norm[i] = std::real(n);
    }

    for (auto const &op : onesite_operators) {
      const int i = op.source_site;
      const auto val = Contract_one_site_MF(Tn_[i], op.op);
      local_obs[op.group][i] = val / norm[i];
    }
  } else {
    std::vector<double> norm(N_UNIT);
    for (int i = 0; i < N_UNIT; ++i) {
      const auto n =
          Contract_one_site(C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i], eTb[i],
                            eTl[i], Tn[i], op_identity[i]);
      norm[i] = std::real(n);
    }
    for (auto const &op : onesite_operators) {
      const int i = op.source_site;
      const auto val = Contract_one_site(C1[i], C2[i], C3[i], C4[i], eTt[i],
                                         eTr[i], eTb[i], eTl[i], Tn[i], op.op);
      local_obs[op.group][i] = val / norm[i];
    }
  }

  time_observable += timer.elapsed();
  return local_obs;
}

template <class ptensor>
void TeNeS<ptensor>::save_onesite(
    std::vector<std::vector<typename TeNeS<ptensor>::tensor_type>> const
        &onesite_obs) {
  if (mpirank != 0) {
    return;
  }

  const int nlops = num_onesite_operators;
  std::string filename = outdir + "/onesite_obs.dat";
  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "    Save onesite observables to " << filename << std::endl;
  }
  std::ofstream ofs(filename.c_str());
  ofs << std::scientific
      << std::setprecision(std::numeric_limits<double>::max_digits10);
  ofs << "# $1: op_group\n";
  ofs << "# $2: site_index\n";
  ofs << "# $3: real\n";
  ofs << "# $4: imag\n";
  ofs << std::endl;

  for (int ilops = 0; ilops < nlops; ++ilops) {
    int num = 0;
    tensor_type sum = 0.0;
    for (int i = 0; i < N_UNIT; ++i) {
      if (std::isnan(std::real(onesite_obs[ilops][i]))) {
        continue;
      }
      num += 1;
      const auto v = onesite_obs[ilops][i];
      sum += v;
      ofs << ilops << " " << i << " " << std::real(v) << " " << std::imag(v)
          << std::endl;
    }
  }
}

template <class ptensor>
auto TeNeS<ptensor>::measure_twosite()
    -> std::vector<std::map<Bond, typename TeNeS<ptensor>::tensor_type>> {
  Timer<> timer;

  const int nlops = num_twosite_operators;
  std::vector<std::map<Bond, tensor_type>> ret(nlops);

  constexpr int nmax = 4;

  std::map<std::tuple<int, int, int>, double> norms;

  for (const auto &op : twosite_operators) {
    const int source = op.source_site;
    const int dx = op.dx[0];
    const int dy = op.dy[0];

    const int ncol = std::abs(dx) + 1;
    const int nrow = std::abs(dy) + 1;
    if (ncol > nmax || nrow > nmax) {
      std::cerr << "Warning: now version of TeNeS does not support too long "
                   "operator"
                << std::endl;
      std::cerr << "group = " << op.group << " (dx = " << dx << ", dy = " << dy
                << ")" << std::endl;
      continue;
    }

    std::vector<const ptensor *> C_(4, nullptr);
    std::vector<const ptensor *> eTt_(ncol, nullptr);
    std::vector<const ptensor *> eTr_(nrow, nullptr);
    std::vector<const ptensor *> eTb_(ncol, nullptr);
    std::vector<const ptensor *> eTl_(nrow, nullptr);

    /*
     * Caution: orders of tensors in unitcell and Contract_* function are
     * different
     *
     * Lattice:
     *
     *    y
     *    ^
     *    |
     *    0--> x
     *
     * Contract_*:
     *    0-->col
     *    |
     *    v
     *    row
     *
     */
    std::vector<std::vector<const ptensor *>> Tn_(
        nrow, std::vector<const ptensor *>(ncol, nullptr));
    std::vector<std::vector<const ptensor *>> op_(
        nrow, std::vector<const ptensor *>(ncol, nullptr));

    std::vector<std::vector<int>> indices(nrow, std::vector<int>(ncol));
    std::vector<ptensor> boundaries;

    int source_col, source_row, target_col, target_row;

    if (dx >= 0) {
      source_col = 0;
      target_col = ncol - 1;
    } else {
      source_col = ncol - 1;
      target_col = 0;
    }
    if (dy >= 0) {
      source_row = nrow - 1;
      target_row = 0;
    } else {
      source_row = 0;
      target_row = nrow - 1;
    }

    if (peps_parameters.MeanField_Env) {
      int iboundary = 0;
      const int nboundary = 2 * (ncol + nrow - 2);
      boundaries.reserve(nboundary);

      for (int row = 0; row < nrow; ++row) {
        for (int col = 0; col < ncol; ++col) {
          const int index =
              lattice.other(source, col - source_col, source_row - row);
          indices[row][col] = index;
          op_[row][col] = &(op_identity[index]);
          if ((0 < row && row < nrow - 1) && (0 < col && col < ncol - 1)) {
            Tn_[row][col] = &(Tn[index]);
          } else {
            boundaries.push_back(Tn[index]);
            Tn_[row][col] = &(boundaries[iboundary++]);
          }
        }
      }
      assert(boundaries.size() == nboundary);

      // absorb MF ENV into center tensors on boundary
      for (int row = 0; row < nrow; ++row) {
        const_cast<ptensor *>(Tn_[row][0])
            ->multiply_vector(lambda_tensor[indices[row][0]][0], 0);
        const_cast<ptensor *>(Tn_[row][ncol - 1])
            ->multiply_vector(lambda_tensor[indices[row][ncol - 1]][2], 2);
      }
      for (int col = 0; col < ncol; ++col) {
        const_cast<ptensor *>(Tn_[0][col])
            ->multiply_vector(lambda_tensor[indices[0][col]][1], 1);
        const_cast<ptensor *>(Tn_[nrow - 1][col])
            ->multiply_vector(lambda_tensor[indices[nrow - 1][col]][3], 3);
      }
    } else {  // Use CTM
      for (int row = 0; row < nrow; ++row) {
        for (int col = 0; col < ncol; ++col) {
          const int index =
              lattice.other(source, col - source_col, source_row - row);
          indices[row][col] = index;
          op_[row][col] = &(op_identity[index]);
          Tn_[row][col] = &(Tn[index]);
        }
        eTl_[row] = &(eTl[indices[row][0]]);
        eTr_[row] = &(eTr[indices[row][ncol - 1]]);
      }
      for (int col = 0; col < ncol; ++col) {
        eTt_[col] = &(eTt[indices[0][col]]);
        eTb_[col] = &(eTb[indices[nrow - 1][col]]);
      }
      C_[0] = &(C1[indices[0][0]]);
      C_[1] = &(C2[indices[0][ncol - 1]]);
      C_[2] = &(C3[indices[nrow - 1][ncol - 1]]);
      C_[3] = &(C4[indices[nrow - 1][0]]);
    }

    const auto norm_key = std::make_tuple(indices[0][0], nrow, ncol);
    auto norm =
        (norms.count(norm_key) ? norms[norm_key]
                               : std::numeric_limits<double>::quiet_NaN());
    if (std::isnan(norm)) {
      if (peps_parameters.MeanField_Env) {
        norm = std::real(Contract_MF(Tn_, op_));
      } else {
        norm = std::real(Contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_));
      }
      norms[norm_key] = norm;
    }

    tensor_type value = 0.0;
    if (op.ops_indices.empty()) {
      if (nrow * ncol == 2) {
        if (nrow == 2) {
          const int top = indices[0][0];
          const int bottom = indices[1][0];
          ptensor o =
              (top == source ? op.op
                             : mptensor::transpose(op.op, {1, 0, 3, 2}));
          value = peps_parameters.MeanField_Env
                      ? Contract_two_sites_vertical_op12_MF(*(Tn_[0][0]),
                                                            *(Tn_[1][0]), o)
                      : Contract_two_sites_vertical_op12(
                            C1[top], C2[top], C3[bottom], C4[bottom], eTt[top],
                            eTr[top], eTr[bottom], eTb[bottom], eTl[bottom],
                            eTl[top], Tn[top], Tn[bottom], o);
        } else {  // ncol == 2
          const int left = indices[0][0];
          const int right = indices[0][1];
          ptensor o =
              (left == source ? op.op
                              : mptensor::transpose(op.op, {1, 0, 3, 2}));
          value = peps_parameters.MeanField_Env
                      ? Contract_two_sites_horizontal_op12_MF(*(Tn_[0][0]),
                                                              *(Tn_[0][1]), o)
                      : Contract_two_sites_horizontal_op12(
                            C1[left], C2[right], C3[right], C4[left], eTt[left],
                            eTt[right], eTr[right], eTb[right], eTb[left],
                            eTl[left], Tn[left], Tn[right], o);
        }
      } else {
        ptensor U, VT;
        std::vector<double> s;
        mptensor::svd(op.op, {0, 2}, {1, 3}, U, s, VT);
        const int ns = s.size();
        for (int is = 0; is < ns; ++is) {
          ptensor source_op =
              reshape(slice(U, 2, is, is + 1), {U.shape()[0], U.shape()[0]});
          op_[source_row][source_col] = &source_op;
          ptensor target_op =
              reshape(slice(VT, 0, is, is + 1), {VT.shape()[1], VT.shape()[1]});
          op_[target_row][target_col] = &target_op;
          auto localvalue =
              peps_parameters.MeanField_Env
                  ? Contract_MF(Tn_, op_)
                  : Contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_);
          value += localvalue * s[is];
        }
      }
    } else {
      op_[source_row][source_col] =
          &(onesite_operators[siteoperator_index(op.source_site,
                                                 op.ops_indices[0])]
                .op);
      const int target_site = lattice.other(op.source_site, dx, dy);
      op_[target_row][target_col] = &(
          onesite_operators[siteoperator_index(target_site, op.ops_indices[1])]
              .op);
      auto localvalue = peps_parameters.MeanField_Env
                            ? Contract_MF(Tn_, op_)
                            : Contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_);
      value += localvalue;
    }
    ret[op.group][{op.source_site, op.dx[0], op.dy[0]}] = value / norm;
  }

  time_observable += timer.elapsed();
  return ret;
}

template <class ptensor>
void TeNeS<ptensor>::save_twosite(
    std::vector<std::map<Bond, typename TeNeS<ptensor>::tensor_type>> const
        &twosite_obs) {
  if (mpirank != 0) {
    return;
  }

  const int nlops = num_twosite_operators;
  std::string filename = outdir + "/twosite_obs.dat";
  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "    Save twosite observables to " << filename << std::endl;
  }
  std::ofstream ofs(filename.c_str());
  ofs << std::scientific
      << std::setprecision(std::numeric_limits<double>::max_digits10);
  ofs << "# $1: op_group\n";
  ofs << "# $2: source_site\n";
  ofs << "# $3: dx\n";
  ofs << "# $4: dy\n";
  ofs << "# $5: real\n";
  ofs << "# $6: imag\n";
  ofs << std::endl;
  for (int ilops = 0; ilops < nlops; ++ilops) {
    tensor_type sum = 0.0;
    int num = 0;
    for (const auto &r : twosite_obs[ilops]) {
      auto bond = r.first;
      auto value = r.second;
      sum += value;
      num += 1;
      ofs << ilops << " " << bond.source_site << " " << bond.dx << " "
          << bond.dy << " " << std::real(value) << " " << std::imag(value)
          << std::endl;
    }
  }
}

template <class ptensor>
std::vector<Correlation> TeNeS<ptensor>::measure_correlation() {
  if (peps_parameters.MeanField_Env) {
    return measure_correlation_mf();
  } else {
    return measure_correlation_ctm();
  }
}

template <class ptensor>
std::vector<Correlation> TeNeS<ptensor>::measure_correlation_ctm() {
  Timer<> timer;

  const int nlops = num_onesite_operators;
  const int r_max = corparam.r_max;
  std::vector<std::vector<int>> r_ops(nlops);
  for (auto ops : corparam.operators) {
    r_ops[std::get<0>(ops)].push_back(std::get<1>(ops));
  }

  std::vector<Correlation> correlations;
  for (int left_index = 0; left_index < N_UNIT; ++left_index) {
    const auto vdim = lattice.virtual_dims[left_index];
    ptensor correlation_T(Shape(CHI, CHI, vdim[0], vdim[0]));
    ptensor correlation_norm(Shape(CHI, CHI, vdim[0], vdim[0]));
    for (int left_ilop = 0; left_ilop < nlops; ++left_ilop) {
      if (r_ops[left_ilop].empty()) {
        continue;
      }

      {  // horizontal
        int left_op_index = siteoperator_index(left_index, left_ilop);
        if (left_op_index < 0) {
          continue;
        }
        const auto left_op = onesite_operators[left_op_index].op;
        StartCorrelation(correlation_T, C1[left_index], C4[left_index],
                         eTt[left_index], eTb[left_index], eTl[left_index],
                         Tn[left_index], left_op);
        StartCorrelation(correlation_norm, C1[left_index], C4[left_index],
                         eTt[left_index], eTb[left_index], eTl[left_index],
                         Tn[left_index], op_identity[left_index]);

        int right_index = left_index;
        for (int r = 0; r < r_max; ++r) {
          right_index = lattice.right(right_index);
          double norm = std::real(FinishCorrelation(
              correlation_norm, C2[right_index], C3[right_index],
              eTt[right_index], eTr[right_index], eTb[right_index],
              Tn[right_index], op_identity[right_index]));
          for (auto right_ilop : r_ops[left_ilop]) {
            int right_op_index = siteoperator_index(right_index, right_ilop);
            if (right_op_index < 0) {
              continue;
            }
            const auto right_op = onesite_operators[right_op_index].op;
            auto val = FinishCorrelation(correlation_T, C2[right_index],
                                         C3[right_index], eTt[right_index],
                                         eTr[right_index], eTb[right_index],
                                         Tn[right_index], right_op) /
                       norm;
            correlations.push_back(Correlation{left_index, r + 1, 0, left_ilop,
                                               right_ilop, std::real(val),
                                               std::imag(val)});
          }

          Transfer(correlation_T, eTt[right_index], eTb[right_index],
                   Tn[right_index]);
          Transfer(correlation_norm, eTt[right_index], eTb[right_index],
                   Tn[right_index]);
        }
      }
      {  // vertical
        int left_op_index = siteoperator_index(left_index, left_ilop);
        if (left_op_index < 0) {
          continue;
        }
        const auto left_op = onesite_operators[left_op_index].op;
        ptensor tn = transpose(Tn[left_index], mptensor::Axes(3, 0, 1, 2, 4));
        StartCorrelation(correlation_T, C4[left_index], C3[left_index],
                         eTl[left_index], eTr[left_index], eTb[left_index], tn,
                         left_op);
        StartCorrelation(correlation_norm, C4[left_index], C3[left_index],
                         eTl[left_index], eTr[left_index], eTb[left_index], tn,
                         op_identity[left_index]);

        int right_index = left_index;
        for (int r = 0; r < r_max; ++r) {
          right_index = lattice.top(right_index);
          tn = transpose(Tn[right_index], mptensor::Axes(3, 0, 1, 2, 4));
          double norm = std::real(FinishCorrelation(
              correlation_norm, C1[right_index], C2[right_index],
              eTl[right_index], eTt[right_index], eTr[right_index], tn,
              op_identity[right_index]));
          for (auto right_ilop : r_ops[left_ilop]) {
            int right_op_index = siteoperator_index(right_index, right_ilop);
            if (right_op_index < 0) {
              continue;
            }
            const auto right_op = onesite_operators[right_op_index].op;
            auto val = FinishCorrelation(correlation_T, C1[right_index],
                                         C2[right_index], eTl[right_index],
                                         eTt[right_index], eTr[right_index], tn,
                                         right_op) /
                       norm;
            correlations.push_back(Correlation{left_index, 0, r + 1, left_ilop,
                                               right_ilop, std::real(val),
                                               std::imag(val)});
          }

          Transfer(correlation_T, eTl[right_index], eTr[right_index], tn);
          Transfer(correlation_norm, eTl[right_index], eTr[right_index], tn);
        }
      }
    }
  }

  time_observable += timer.elapsed();
  return correlations;
}

template <class ptensor>
void TeNeS<ptensor>::transfer_matvec_horizontal(ptensor &vec, int y) {
  for (int x = 0; lattice.LX; ++x) {
    const int site = x + lattice.LX * y;
    ////////////////////////////////////////////////////////////
    // ./transfer_matvec.dat
    ////////////////////////////////////////////////////////////
    // (eTt[site]*(vec*eTb[site]))
    // cpu_cost= 131072  memory= 12544
    // final_bond_order (et_r, eb_r)
    ////////////////////////////////////////////////////////////
    vec = tensordot(eTt[site], tensordot(vec, eTb[site], Axes(1), Axes(1)),
                    Axes(0, 2, 3), Axes(0, 2, 3));
  }
}

template <class ptensor>
void TeNeS<ptensor>::transfer_matvec_vertical(ptensor &vec, int x) {
  for (int y = 0; lattice.LY; ++y) {
    const int site = x + lattice.LX * y;
    ////////////////////////////////////////////////////////////
    // ./transfer_matvec.dat
    ////////////////////////////////////////////////////////////
    // (eTt[site]*(vec*eTb[site]))
    // cpu_cost= 131072  memory= 12544
    // final_bond_order (et_r, eb_r)
    ////////////////////////////////////////////////////////////
    vec = tensordot(eTl[site], tensordot(vec, eTr[site], Axes(1), Axes(1)),
                    Axes(0, 2, 3), Axes(0, 2, 3));
  }
}

template <class ptensor>
std::vector<typename TeNeS<ptensor>::transfer_matrix_eigenvalues_type>
TeNeS<ptensor>::measure_transfer_matrix_eigenvalues() {
  // res[id][0]: direction
  // res[id][1]: coord
  // res[id][2]: value
  // res[id][3]: eigval 0
  // res[id][4]: eigval 1
  // res[id][5]: eigval 2
  std::vector<transfer_matrix_eigenvalues_type> res;

  std::array<std::function<const ptensor(int, int)>, 2> topedge{
      [&](int run, int fixed) { return eTt[lattice.index(run, fixed)]; },
      [&](int run, int fixed) { return eTl[lattice.index(fixed, run)]; },
  };

  std::array<std::function<const ptensor(int, int)>, 2> bottomedge{
      [&](int run, int fixed) { return eTb[lattice.index(run, fixed - 1)]; },
      [&](int run, int fixed) { return eTr[lattice.index(fixed + 1, run)]; },
  };

  ptensor T;
  std::vector<double> lambda;
  for (int dir = 0; dir < 2; ++dir) {
    int L = dir == 0 ? LX : LY;
    int W = dir == 0 ? LY : LX;
    for (int fixed = 0; fixed < W; ++fixed) {
      int run = 0;
      T = tensordot(topedge[dir](run, fixed), bottomedge[dir](run, fixed),
                    Axes(2, 3), Axes(2, 3));
      for (run = 1; run < L; ++run) {
        T = tensordot(
            T,
            tensordot(topedge[dir](run, fixed), bottomedge[dir](run, fixed),
                      Axes(2, 3), Axes(2, 3)),
            Axes(1, 2), Axes(0, 3));
      }
      psvd(T, Axes(0, 3), Axes(1, 2), lambda, 3);
      for (int i = 0; i < 3; ++i) {
        lambda[i] = std::sqrt(lambda[i]);
      }
      res.push_back(std::make_tuple(0, fixed, lambda[0], lambda[1], lambda[2]));
    }
  }

  return res;
}

template <class ptensor>
void TeNeS<ptensor>::save_correlation_length(
    std::vector<typename TeNeS<ptensor>::transfer_matrix_eigenvalues_type> const
        &lambdas) {
  if (mpirank != 0) {
    return;
  }
  std::string filename = outdir + "/correlation_length.dat";
  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "    Save correlation length to " << filename << std::endl;
  }
  std::ofstream ofs(filename.c_str());
  ofs << std::scientific
      << std::setprecision(std::numeric_limits<double>::max_digits10);
  ofs << "# $1: direction\n";
  ofs << "# $2: col or row index\n";
  ofs << "# $3: value\n";
  ofs << "# $4: correction x (log(l1/l2))\n";
  ofs << std::endl;

  for (const auto &lambda : lambdas) {
    int dir, x;
    double l0, l1, l2;
    std::tie(dir, x, l0, l1, l2) = lambda;

    int L = dir == 0 ? LX : LY;
    double correlation_length = L / std::log(l0 / l1);
    double correction_x = std::log(l1 / l2);

    ofs << dir << " ";
    ofs << x << " ";
    ofs << correlation_length << " ";
    ofs << correction_x << " ";
    ofs << std::endl;
  }
}

template <class ptensor>
std::vector<Correlation> TeNeS<ptensor>::measure_correlation_mf() {
  Timer<> timer;

  const int nlops = num_onesite_operators;
  const int r_max = corparam.r_max;
  std::vector<std::vector<int>> r_ops(nlops);
  for (auto ops : corparam.operators) {
    r_ops[std::get<0>(ops)].push_back(std::get<1>(ops));
  }

  std::vector<ptensor> Tn_horizontal(Tn.begin(), Tn.end());
  std::vector<ptensor> Tn_vertical(Tn.begin(), Tn.end());
  for (int index = 0; index < N_UNIT; ++index) {
    std::vector<std::vector<double>> const &lambda = lambda_tensor[index];
    Tn_horizontal[index].multiply_vector(lambda[1], 1, lambda[3], 3);
    Tn_vertical[index].multiply_vector(lambda[0], 0, lambda[2], 2);
  }

  std::vector<Correlation> correlations;
  for (int left_index = 0; left_index < N_UNIT; ++left_index) {
    const auto vdim = lattice.virtual_dims[left_index];
    ptensor correlation_T(Shape(vdim[0], vdim[0]));
    ptensor correlation_norm(Shape(vdim[0], vdim[0]));
    for (int left_ilop = 0; left_ilop < nlops; ++left_ilop) {
      if (r_ops[left_ilop].empty()) {
        continue;
      }

      {  // horizontal
        int left_op_index = siteoperator_index(left_index, left_ilop);
        if (left_op_index < 0) {
          continue;
        }
        ptensor T = Tn_horizontal[left_index];
        T.multiply_vector(lambda_tensor[left_index][0], 0);
        const auto left_op = onesite_operators[left_op_index].op;
        StartCorrelation_MF(correlation_T, T, left_op, 2);
        StartCorrelation_MF(correlation_norm, T, op_identity[left_index], 2);

        int right_index = left_index;
        for (int r = 0; r < r_max; ++r) {
          right_index = lattice.right(right_index);
          T = Tn_horizontal[right_index];
          T.multiply_vector(lambda_tensor[right_index][2], 2);
          double norm = std::real(FinishCorrelation_MF(
              correlation_norm, T, op_identity[right_index], 2));
          for (auto right_ilop : r_ops[left_ilop]) {
            int right_op_index = siteoperator_index(right_index, right_ilop);
            if (right_op_index < 0) {
              continue;
            }
            const auto right_op = onesite_operators[right_op_index].op;
            auto val =
                FinishCorrelation_MF(correlation_T, T, right_op, 2) / norm;
            correlations.push_back(Correlation{left_index, r + 1, 0, left_ilop,
                                               right_ilop, std::real(val),
                                               std::imag(val)});
          }

          Transfer_MF(correlation_T, Tn_horizontal[right_index], 2);
          Transfer_MF(correlation_norm, Tn_horizontal[right_index], 2);
        }
      }
      {  // vertical
        int left_op_index = siteoperator_index(left_index, left_ilop);
        if (left_op_index < 0) {
          continue;
        }
        ptensor T = Tn_vertical[left_index];
        T.multiply_vector(lambda_tensor[left_index][3], 3);
        const auto left_op = onesite_operators[left_op_index].op;
        StartCorrelation_MF(correlation_T, T, left_op, 1);
        StartCorrelation_MF(correlation_norm, T, op_identity[left_index], 1);

        int right_index = left_index;
        for (int r = 0; r < r_max; ++r) {
          right_index = lattice.top(right_index);
          T = Tn_vertical[right_index];
          T.multiply_vector(lambda_tensor[right_index][1], 1);
          double norm = std::real(FinishCorrelation_MF(
              correlation_norm, T, op_identity[right_index], 1));
          for (auto right_ilop : r_ops[left_ilop]) {
            int right_op_index = siteoperator_index(right_index, right_ilop);
            if (right_op_index < 0) {
              continue;
            }
            const auto right_op = onesite_operators[right_op_index].op;
            auto val =
                FinishCorrelation_MF(correlation_T, T, right_op, 1) / norm;
            correlations.push_back(Correlation{left_index, 0, r + 1, left_ilop,
                                               right_ilop, std::real(val),
                                               std::imag(val)});
          }

          Transfer_MF(correlation_T, Tn_vertical[right_index], 1);
          Transfer_MF(correlation_norm, Tn_vertical[right_index], 1);
        }
      }
    }
  }

  time_observable += timer.elapsed();
  return correlations;
}

template <class ptensor>
void TeNeS<ptensor>::save_correlation(
    std::vector<Correlation> const &correlations) {
  if (mpirank != 0) {
    return;
  }
  std::string filename = outdir + "/correlation.dat";
  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "    Save long-range correlations to " << filename
              << std::endl;
  }
  std::ofstream ofs(filename.c_str());
  ofs << std::scientific
      << std::setprecision(std::numeric_limits<double>::max_digits10);
  ofs << "# $1: left_op\n";
  ofs << "# $2: left_site\n";
  ofs << "# $3: right_op\n";
  ofs << "# $4: right_dx\n";
  ofs << "# $5: right_dy\n";
  ofs << "# $6: real\n";
  ofs << "# $7: imag\n";
  ofs << std::endl;
  for (auto const &cor : correlations) {
    ofs << cor.left_op << " " << cor.left_index << " " << cor.right_op << " "
        << cor.right_dx << " " << cor.right_dy << " " << cor.real << " "
        << cor.imag << " " << std::endl;
  }
}

template <class ptensor>
void TeNeS<ptensor>::measure() {
  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "Start calculating observables" << std::endl;
    std::cout << "  Start updating environment" << std::endl;
  }
  if (!peps_parameters.MeanField_Env) {
    update_CTM();
  }

  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "  Start calculating onesite operators" << std::endl;
  }
  auto onesite_obs = measure_onesite();
  save_onesite(onesite_obs);

  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "  Start calculating twosite operators" << std::endl;
  }
  auto twosite_obs = measure_twosite();
  save_twosite(twosite_obs);

  if (corparam.r_max > 0) {
    if (peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "  Start calculating long range correlation" << std::endl;
    }
    auto correlations = measure_correlation();
    save_correlation(correlations);
  }

  auto correlation_length = measure_transfer_matrix_eigenvalues();
  save_correlation_length(correlation_length);

  if (mpirank == 0) {
    std::vector<tensor_type> loc_obs(num_onesite_operators);
    int numsites = 0;
    for (int i = 0; i < N_UNIT; ++i) {
      if (lattice.physical_dims[i] > 1) {
        ++numsites;
        for (int ilops = 0; ilops < num_onesite_operators; ++ilops) {
          loc_obs[ilops] += onesite_obs[ilops][i];
        }
      }
    }
    std::vector<tensor_type> two_obs(num_twosite_operators);
    for (int iops = 0; iops < num_twosite_operators; ++iops) {
      for (const auto &obs : twosite_obs[iops]) {
        two_obs[iops] += obs.second;
      }
    }

    auto energy = 0.0;
    for (const auto &obs : twosite_obs[0]) {
      energy += std::real(obs.second);
    }

    {
      const double invV = 1.0 / numsites;
      std::string filename = outdir + "/density.dat";
      std::ofstream ofs(filename.c_str());
      ofs << std::scientific
          << std::setprecision(std::numeric_limits<double>::max_digits10);

      if (mpirank == 0) {
        for (int ilops = 0; ilops < num_onesite_operators; ++ilops) {
          const auto v = loc_obs[ilops] * invV;
          ofs << onesite_operator_names[ilops] << " = ";
          if (std::real(v) >= 0.0) {
            ofs << " ";
          }
          ofs << std::real(v) << " ";
          if (std::imag(v) >= 0.0) {
            ofs << " ";
          }
          ofs << std::imag(v) << std::endl;
        }

        for (int ilops = 0; ilops < num_twosite_operators; ++ilops) {
          const auto v = two_obs[ilops] * invV;
          ofs << twosite_operator_names[ilops] << " = ";
          if (std::real(v) >= 0.0) {
            ofs << " ";
          }
          ofs << std::real(v) << " ";
          if (std::imag(v) >= 0.0) {
            ofs << " ";
          }
          ofs << std::imag(v) << std::endl;
        }
        std::cout << "    Save observable densities to " << filename
                  << std::endl;
      }
    }
    {
      std::string filename = outdir + "/parameters.dat";
      std::ofstream ofs(filename.c_str(), std::ios::out | std::ios::app);
      ofs << "finish_datetime = " << datetime() << std::endl;
    }

    if (peps_parameters.print_level >= PrintLevel::info) {
      const double invV = 1.0 / numsites;
      std::cout << std::endl;

      std::cout << "Onesite observables per site:" << std::endl;
      for (int ilops = 0; ilops < num_onesite_operators; ++ilops) {
        const auto v = loc_obs[ilops] * invV;
        std::cout << "  " << onesite_operator_names[ilops] << " = "
                  << std::real(v) << " " << std::imag(v) << std::endl;
      }

      std::cout << "Twosite observables per site:" << std::endl;
      for (int ilops = 0; ilops < num_twosite_operators; ++ilops) {
        const auto v = two_obs[ilops] * invV;
        std::cout << "  " << twosite_operator_names[ilops] << " = "
                  << std::real(v) << " " << std::imag(v) << std::endl;
      }
    }
  }  // end of if(mpirank == 0)
}

template <class ptensor>
void TeNeS<ptensor>::summary() const {
  if (mpirank == 0) {
    const double time_all = timer_all.elapsed();
    {
      std::string filename = outdir + "/time.dat";
      std::ofstream ofs(filename.c_str());
      ofs << "time all           = " << time_all << std::endl;
      ofs << "time simple update = " << time_simple_update << std::endl;
      ofs << "time full update   = " << time_full_update << std::endl;
      ofs << "time environmnent  = " << time_environment << std::endl;
      ofs << "time observable    = " << time_observable << std::endl;
      if (peps_parameters.print_level >= PrintLevel::info) {
        std::cout << "    Save elapsed times to " << filename << std::endl;
      }
    }
    if (peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "Wall times [sec.]:" << std::endl;
      std::cout << "  all           = " << time_all << std::endl;
      std::cout << "  simple update = " << time_simple_update << std::endl;
      std::cout << "  full update   = " << time_full_update << std::endl;
      std::cout << "  environmnent  = " << time_environment << std::endl;
      std::cout << "  observable    = " << time_observable << std::endl;
      std::cout << std::endl << "Done." << std::endl;
    }
  }
}

template <class ptensor>
void TeNeS<ptensor>::save_tensors() const {
  std::string const &save_dir = peps_parameters.tensor_save_dir;
  if (save_dir.empty()) {
    return;
  }
  {
    // metadata
    std::string filename = save_dir + "/params.dat";
    std::ofstream ofs(filename.c_str());

    constexpr int tensor_format_version = 1;
    ofs << tensor_format_version << " # Format_Version\n";
    ofs << N_UNIT << " # N_UNIT\n";
    ofs << CHI << " # CHI\n";
    for (int i = 0; i < N_UNIT; ++i) {
      for (int j = 0; j < nleg; ++j) {
        ofs << lattice.virtual_dims[i][j] << " ";
      }
      ofs << lattice.physical_dims[i] << " # Shape of Tn[" << i << "]\n";
    }
  }
  for (int i = 0; i < N_UNIT; ++i) {
    std::string filename = save_dir + "/";
    std::string suffix = "_" + std::to_string(i) + ".dat";
    Tn[i].save((filename + "T" + suffix).c_str());
    eTt[i].save((filename + "Et" + suffix).c_str());
    eTr[i].save((filename + "Er" + suffix).c_str());
    eTb[i].save((filename + "Eb" + suffix).c_str());
    eTl[i].save((filename + "El" + suffix).c_str());
    C1[i].save((filename + "C1" + suffix).c_str());
    C2[i].save((filename + "C2" + suffix).c_str());
    C3[i].save((filename + "C3" + suffix).c_str());
    C4[i].save((filename + "C4" + suffix).c_str());
  }
  if (mpirank == 0) {
    for (int i = 0; i < N_UNIT; ++i) {
      std::ofstream ofs(save_dir + "/lambda_" + std::to_string(i) + ".dat");
      for (int j = 0; j < nleg; ++j) {
        for (int k = 0; k < lattice.virtual_dims[i][j]; ++k) {
          ofs << lambda_tensor[i][j][k] << std::endl;
        }
      }
    }
  }
  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "Tensors saved in " << save_dir << std::endl;
  }
}

template <class ptensor>
void TeNeS<ptensor>::load_tensors() {
  std::string const &load_dir = peps_parameters.tensor_load_dir;

  if (!util::isdir(load_dir)) {
    std::string msg = load_dir + " does not exists.";
    throw tenes::load_error(msg);
  }

  int tensor_format_version = 0;
  if (mpirank == 0) {
    std::string filename = load_dir + "/params.dat";
    std::string line;
    if (util::path_exists(filename)) {
      std::ifstream ifs(filename.c_str());
      std::getline(ifs, line);
      tensor_format_version = std::stoi(util::drop_comment(line));
    }
  }
  bcast(tensor_format_version, 0, comm);
  if (tensor_format_version == 0) {
    load_tensors_v0();
  } else if (tensor_format_version == 1) {
    load_tensors_v1();
  } else {
    std::stringstream ss;
    ss << "ERROR: Unknown checkpoint format version: " << tensor_format_version;
    throw tenes::load_error(ss.str());
  }
}

template <class ptensor>
void TeNeS<ptensor>::load_tensors_v1() {
  std::string const &load_dir = peps_parameters.tensor_load_dir;

  int loaded_CHI = 1;
  std::vector<std::vector<int>> loaded_shape(N_UNIT,
                                             std::vector<int>(nleg + 1));
  if (mpirank == 0) {
    std::string filename = load_dir + "/params.dat";
    std::string line;
    std::ifstream ifs(filename.c_str());
    std::getline(ifs, line);

    std::getline(ifs, line);
    const int loaded_N_UNIT = std::stoi(util::drop_comment(line));
    if (N_UNIT != loaded_N_UNIT) {
      std::stringstream ss;
      ss << "ERROR: N_UNIT is " << N_UNIT << " but loaded N_UNIT has "
         << loaded_N_UNIT << std::endl;
      throw tenes::load_error(ss.str());
    }

    std::getline(ifs, line);
    loaded_CHI = std::stoi(util::drop_comment(line));
    if (CHI != loaded_CHI) {
      if (peps_parameters.print_level >= PrintLevel::info) {
        std::cout << "WARNING: parameters.ctm.dimension is " << CHI
                  << " but loaded tensors have CHI = " << loaded_CHI
                  << std::endl;
      }
    }

    for (int i = 0; i < N_UNIT; ++i) {
      std::getline(ifs, line);
      const auto shape = util::split(util::drop_comment(line));
      for (int j = 0; j < nleg; ++j) {
        loaded_shape[i][j] = std::stoi(shape[j]);
        const int vd_param = lattice.virtual_dims[i][j];
        if (vd_param != loaded_shape[i][j]) {
          if (peps_parameters.print_level >= PrintLevel::info) {
            std::cout << "WARNING: virtual dimension of the leg " << j
                      << " of the tensor " << i << " is " << vd_param
                      << " but loaded tensor has " << loaded_shape[i][j]
                      << std::endl;
          }
        }
      }
      loaded_shape[i][nleg] = std::stoi(shape[nleg]);
      const int pdim = lattice.physical_dims[i];
      if (pdim != loaded_shape[i][nleg]) {
        std::stringstream ss;
        ss << "ERROR: dimension of the physical bond of the tensor " << i
           << " is " << pdim << " but loaded tensor has "
           << loaded_shape[i][nleg] << std::endl;
        throw tenes::load_error(ss.str());
      }
    }
  }
  for (int i = 0; i < N_UNIT; ++i) {
    bcast(loaded_shape[i], 0, comm);
  }

#define LOAD_TENSOR_(A, name)                      \
  do {                                             \
    ptensor temp;                                  \
    temp.load((filename + name + suffix).c_str()); \
    A = resize_tensor(temp, A.shape());            \
  } while (false)

  for (int i = 0; i < N_UNIT; ++i) {
    std::string filename = load_dir + "/";
    std::string suffix = "_" + std::to_string(i) + ".dat";

    LOAD_TENSOR_(Tn[i], "T");
    LOAD_TENSOR_(eTl[i], "El");
    LOAD_TENSOR_(eTt[i], "Et");
    LOAD_TENSOR_(eTr[i], "Er");
    LOAD_TENSOR_(eTb[i], "Eb");
    LOAD_TENSOR_(C1[i], "C1");
    LOAD_TENSOR_(C2[i], "C2");
    LOAD_TENSOR_(C3[i], "C3");
    LOAD_TENSOR_(C4[i], "C4");
  }
#undef LOAD_TENSOR_

  std::vector<double> ls;
  if (mpirank == 0) {
    for (int i = 0; i < N_UNIT; ++i) {
      std::ifstream ifs(load_dir + "/lambda_" + std::to_string(i) + ".dat");
      for (int j = 0; j < nleg; ++j) {
        for (int k = 0; k < loaded_shape[i][j]; ++k) {
          double temp = 0.0;
          ifs >> temp;
          ls.push_back(temp);
        }
      }
    }
  }
  bcast(ls, 0, comm);
  int index = 0;
  for (int i = 0; i < N_UNIT; ++i) {
    const auto vdim = lattice.virtual_dims[i];
    for (int j = 0; j < nleg; ++j) {
      lambda_tensor[i][j].clear();
      for (int k = 0; k < loaded_shape[i][j]; ++k) {
        lambda_tensor[i][j].push_back(ls[index]);
        ++index;
      }
      lambda_tensor[i][j].resize(vdim[j]);
    }
  }
}

template <class ptensor>
void TeNeS<ptensor>::load_tensors_v0() {
  std::string const &load_dir = peps_parameters.tensor_load_dir;

  // load from the checkpoint
  if (!util::isdir(load_dir)) {
    std::string msg = load_dir + " does not exists.";
    throw tenes::load_error(msg);
  }
  for (int i = 0; i < N_UNIT; ++i) {
    std::string filename = load_dir + "/";
    std::string suffix = "_" + std::to_string(i) + ".dat";
    Tn[i].load((filename + "T" + suffix).c_str());
    eTt[i].load((filename + "Et" + suffix).c_str());
    eTr[i].load((filename + "Er" + suffix).c_str());
    eTb[i].load((filename + "Eb" + suffix).c_str());
    eTl[i].load((filename + "El" + suffix).c_str());
    C1[i].load((filename + "C1" + suffix).c_str());
    C2[i].load((filename + "C2" + suffix).c_str());
    C3[i].load((filename + "C3" + suffix).c_str());
    C4[i].load((filename + "C4" + suffix).c_str());
  }
  std::vector<double> ls;
  if (mpirank == 0) {
    for (int i = 0; i < N_UNIT; ++i) {
      const auto vdim = lattice.virtual_dims[i];
      std::ifstream ifs(load_dir + "/lambda_" + std::to_string(i) + ".dat");
      for (int j = 0; j < nleg; ++j) {
        for (int k = 0; k < vdim[j]; ++k) {
          double temp = 0.0;
          ifs >> temp;
          ls.push_back(temp);
        }
      }
    }
  }
  bcast(ls, 0, comm);
  int index = 0;
  for (int i = 0; i < N_UNIT; ++i) {
    const auto vdim = lattice.virtual_dims[i];
    for (int j = 0; j < nleg; ++j) {
      for (int k = 0; k < vdim[j]; ++k) {
        lambda_tensor[i][j][k] = ls[index];
        ++index;
      }
    }
  }

  // overwrite dimensions
  const Shape Cshape = C1[0].shape();
  if (CHI != Cshape[0]) {
    if (peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "WARNING: parameters.ctm.dimension is " << CHI
                << " but loaded tensors have CHI = " << Cshape[0] << std::endl;
    }
  }
  for (int i = 0; i < N_UNIT; ++i) {
    const Shape Tshape = Tn[i].shape();
    const int pdim = lattice.physical_dims[i];
    if (pdim != Tshape[4]) {
      std::stringstream ss;
      ss << "ERROR: dimension of the physical bond of the tensor " << i
         << " is " << pdim << " but loaded tensor has " << Tshape[4]
         << std::endl;
      throw tenes::input_error(ss.str());
    }

    for (int l = 0; l < nleg; ++l) {
      const int vd_param = lattice.virtual_dims[i][l];
      const int vd_loaded = Tshape[l];
      if (vd_param != vd_loaded) {
        if (peps_parameters.print_level >= PrintLevel::info) {
          std::cout << "WARNING: virtual dimension of the leg " << l
                    << " of the tensor " << i << " is " << vd_param
                    << " but loaded tensor has " << vd_loaded << std::endl;
        }
      }
    }
  }
}

template <class tensor>
int tenes(MPI_Comm comm, PEPS_Parameters peps_parameters, Lattice lattice,
          NNOperators<tensor> simple_updates, NNOperators<tensor> full_updates,
          Operators<tensor> onesite_operators,
          Operators<tensor> twosite_operators, CorrelationParameter corparam) {
  TeNeS<tensor> tns(comm, peps_parameters, lattice, simple_updates,
                    full_updates, onesite_operators, twosite_operators,
                    corparam);
  tns.optimize();
  tns.save_tensors();
  if (peps_parameters.to_measure) {
    tns.measure();
  }
  tns.summary();
  return 0;
}

// template specialization
template int tenes<real_tensor>(MPI_Comm comm, PEPS_Parameters peps_parameters,
                                Lattice lattice,
                                NNOperators<real_tensor> simple_updates,
                                NNOperators<real_tensor> full_updates,
                                Operators<real_tensor> onesite_operators,
                                Operators<real_tensor> twosite_operators,
                                CorrelationParameter corparam);

template int tenes<complex_tensor>(MPI_Comm comm,
                                   PEPS_Parameters peps_parameters,
                                   Lattice lattice,
                                   NNOperators<complex_tensor> simple_updates,
                                   NNOperators<complex_tensor> full_updates,
                                   Operators<complex_tensor> onesite_operators,
                                   Operators<complex_tensor> twosite_operators,
                                   CorrelationParameter corparam);

}  // end of namespace tenes
