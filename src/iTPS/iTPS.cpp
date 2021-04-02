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

#include "iTPS.hpp"
#include "transfer_matrix.hpp"

#define _USE_MATH_DEFINES
#include <sys/stat.h>
#include <algorithm>
#include <complex>
#include <ctime>
#include <limits>
#include <map>
#include <memory>
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

#include <mptensor/rsvd.hpp>

#include "../tensor.hpp"

#include "../operator.hpp"
#include "../printlevel.hpp"
#include "../timer.hpp"
#include "../util/file.hpp"
#include "../util/string.hpp"
#include "../util/type_traits.hpp"
#include "../util/datetime.hpp"

#include "PEPS_Parameters.hpp"
#include "core/ctm.hpp"

namespace tenes {

template <class tensor>
iTPS<tensor>::iTPS(MPI_Comm comm_, PEPS_Parameters peps_parameters_,
                      SquareLattice lattice_, NNOperators<tensor> simple_updates_,
                      NNOperators<tensor> full_updates_,
                      Operators<tensor> onesite_operators_,
                      Operators<tensor> twosite_operators_,
                      CorrelationParameter corparam_,
                      TransferMatrix_Parameters clength_param_)
    : comm(comm_),
      peps_parameters(peps_parameters_),
      lattice(lattice_),
      simple_updates(simple_updates_),
      full_updates(full_updates_),
      onesite_operators(onesite_operators_),
      twosite_operators(twosite_operators_),
      corparam(corparam_),
      tmatrix_param(clength_param_),
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
  mptensor::random_tensor::set_seed(seed + mpirank);

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
    ofs << "start_datetime =  " << util::datetime() << std::endl;
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
void iTPS<ptensor>::update_CTM() {
  Timer<> timer;
  Calc_CTM_Environment(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, peps_parameters,
                       lattice);
  time_environment += timer.elapsed();
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // end of namespace tenes
