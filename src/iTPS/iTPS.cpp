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

#ifdef _NO_OMP
int omp_get_max_threads() { return 1; }
#else
#include <omp.h>
#endif

#define _USE_MATH_DEFINES
#include <algorithm>
#include <ctime>
#include <iostream>
#include <string>
#include <set>

#include <mptensor/rsvd.hpp>

#include "../tensor.hpp"

#include "../operator.hpp"
#include "../printlevel.hpp"
#include "../timer.hpp"
#include "../util/file.hpp"
#include "../util/datetime.hpp"

#include "PEPS_Parameters.hpp"
#include "transfer_matrix.hpp"
#include "core/ctm.hpp"

namespace tenes {
namespace itps {

template <class tensor>
iTPS<tensor>::iTPS(MPI_Comm comm_, PEPS_Parameters peps_parameters_,
                   SquareLattice lattice_,
                   EvolutionOperators<tensor> simple_updates_,
                   EvolutionOperators<tensor> full_updates_,
                   Operators<tensor> onesite_operators_,
                   Operators<tensor> twosite_operators_,
                   Operators<tensor> multisite_operators_,
                   CorrelationParameter corparam_,
                   TransferMatrix_Parameters clength_param_)
    : comm(comm_),
      peps_parameters(peps_parameters_),
      lattice(lattice_),
      simple_updates(simple_updates_),
      full_updates(full_updates_),
      onesite_operators(onesite_operators_),
      twosite_operators(twosite_operators_),
      multisite_operators(multisite_operators_),
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
    int Dmax = 0;
    for (int site = 0; site < N_UNIT; ++site) {
      for (int leg = 0; leg < 4; ++leg) {
        Dmax = std::max(lattice.virtual_dims[site][leg], Dmax);
      }
    }
    if (peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "Bond dimensions:\n";
      std::cout << "   D  (Bulk): " << Dmax << "\n";
      std::cout << "   chi (CTM): " << CHI << std::endl;
      if (peps_parameters.calcmode == PEPS_Parameters::CalculationMode::finite_temperature){
	if (CHI < Dmax) {
	  std::cerr << "WARNING: CTM may be too small (chi < D) for finite_temperature mode" << std::endl;
	}
	else {
	  if (CHI < Dmax * Dmax) {
	    std::cerr << "WARNING: CTM may be too small (chi < D*D)" << std::endl;
	  }
	}
      }
    }

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

  if (peps_parameters.calcmode == PEPS_Parameters::CalculationMode::finite_temperature){
    initialize_tensors_density();
  } else {
    initialize_tensors();
  }

  std::set<int> simple_update_groups;
  bool notwarned = true;
  for (auto const &op : simple_updates) {
    auto g = op.group;
    if (g < 0) {
      std::stringstream ss;
      ss << "ERROR: a simple update has negative group number " << g;
      throw std::runtime_error(ss.str());
    }
    if (notwarned && g > 0) {
      std::cerr << "WARNING: a simple update has nonzero group number " << g
                << std::endl;
      std::cerr << "         This feature is reserved for future use."
                << std::endl;
      std::cerr << "         Currently, all simple updates with nonzero group "
                   "number are ignored."
                << std::endl;
      notwarned = false;
    }
    simple_update_groups.insert(g);
  }
  num_simple_update_groups = *simple_update_groups.rbegin() + 1;
  for (int g = 0; g < num_simple_update_groups; ++g) {
    auto it = simple_update_groups.find(g);
    if (it == simple_update_groups.end()) {
      std::stringstream ss;
      ss << "ERROR: no simple update with group number " << g;
      throw std::runtime_error(ss.str());
    }
  }

  notwarned = true;

  std::set<int> full_update_groups;
  for (auto const &op : full_updates) {
    auto g = op.group;
    if (g < 0) {
      std::stringstream ss;
      ss << "ERROR: a full update has negative group number " << g;
      throw std::runtime_error(ss.str());
    }
    if (notwarned && g > 0) {
      std::cerr << "WARNING: a full update has nonzero group number " << g
                << std::endl;
      std::cerr << "         This feature is reserved for future use."
                << std::endl;
      std::cerr << "         Currently, all full updates with nonzero group "
                   "number are ignored."
                << std::endl;
      notwarned = false;
    }
    full_update_groups.insert(g);
  }
  num_full_update_groups = *full_update_groups.rbegin() + 1;
  for (int g = 0; g < num_full_update_groups; ++g) {
    auto it = full_update_groups.find(g);
    if (it == full_update_groups.end()) {
      std::stringstream ss;
      ss << "ERROR: no full update with group number " << g;
      throw std::runtime_error(ss.str());
    }
  }

  if (peps_parameters.calcmode !=
      PEPS_Parameters::CalculationMode::ground_state) {
    if (peps_parameters.num_simple_step[0] > 0 &&
        peps_parameters.num_full_step[0]) {
      std::string msg =
          "ERROR: While mode is not ground state calculation,\n"
          "       simple update and full update are both required.";
      throw std::runtime_error(msg);
    }
    if (peps_parameters.num_simple_step[0] > 0) {
      int nss_size = peps_parameters.num_simple_step.size();
      if (nss_size != num_simple_update_groups) {
        if (nss_size == 1) {
          peps_parameters.num_simple_step.resize(
              num_simple_update_groups, peps_parameters.num_simple_step[0]);
        } else {
          std::stringstream msg;
          msg << "ERROR: size of num_step " << nss_size
              << " is not equal to the number of simple update groups.";
          throw std::runtime_error(msg.str());
        }
      }

      int tss_size = peps_parameters.tau_simple_step.size();
      if (tss_size != num_simple_update_groups) {
        if (tss_size == 1) {
          peps_parameters.tau_simple_step.resize(
              num_simple_update_groups, peps_parameters.tau_simple_step[0]);
        } else {
          std::stringstream msg;
          msg << "ERROR: size of tau " << nss_size
              << " is not equal to the number of simple update groups.";
          throw std::runtime_error(msg.str());
        }
      }
    } else {
      int nss_size = peps_parameters.num_full_step.size();
      if (nss_size != num_full_update_groups) {
        if (nss_size == 1) {
          peps_parameters.num_full_step.resize(
              num_full_update_groups, peps_parameters.num_full_step[0]);
        } else {
          std::stringstream msg;
          msg << "ERROR: size of num_step " << nss_size
              << " is not equal to the number of full update groups.";
          throw std::runtime_error(msg.str());
        }
      }
      int tss_size = peps_parameters.tau_full_step.size();
      if (tss_size != num_full_update_groups) {
        if (tss_size == 1) {
          peps_parameters.tau_full_step.resize(
              num_full_update_groups, peps_parameters.tau_full_step[0]);
        } else {
          std::stringstream msg;
          msg << "ERROR: size of tau " << tss_size
              << " is not equal to the number of full update groups.";
          throw std::runtime_error(msg.str());
        }
      }
    }
  } else {
  }

  size_t maxlength = std::string("Energy").size();
  int maxops = -1;
  for (auto const &op : onesite_operators) {
    maxops = std::max(op.group, maxops);
  }
  num_onesite_operators = maxops + 1;
  onesite_operator_names.resize(num_onesite_operators);
  onesite_operator_counts.resize(num_onesite_operators);
  for (auto const &op : onesite_operators) {
    ++onesite_operator_counts[op.group];
    std::string name = op.name;
    if (name.empty()) {
      std::stringstream ss;
      ss << "onesite[" << op.group << "]";
      name = ss.str();
    }
    if (onesite_operator_names[op.group].empty()) {
      onesite_operator_names[op.group] = op.name;
    } else if (onesite_operator_names[op.group] != op.name) {
      std::stringstream ss;
      ss << "ERROR: onesite operator group " << op.group
         << " has two different names: " << name << " and "
         << onesite_operator_names[op.group] << "\n";
      ss << "       onesite operators with the same group must "
             "have the same name.\n";
      throw std::runtime_error(ss.str());
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
  twosite_operator_counts.resize(num_twosite_operators);
  for (auto const &op : twosite_operators) {
    ++twosite_operator_counts[op.group];
    std::string name = op.name;
    if (name.empty()) {
      std::stringstream ss;
      ss << "twosite[" << op.group << "]";
      name = ss.str();
    }
    if (twosite_operator_names[op.group].empty()) {
      twosite_operator_names[op.group] = op.name;
    } else if (twosite_operator_names[op.group] != name) {
      std::stringstream ss;
      ss << "ERROR: twosite operator group " << op.group
         << " has two names: " << twosite_operator_names[op.group] << " and "
         << name << "\n";
      ss << "       twosite operators with the same group must "
             "have the same name.\n";
      throw std::runtime_error(ss.str());
    }
  }
  for (auto const &s : twosite_operator_names) {
    maxlength = std::max(s.size(), maxlength);
  }

  maxops = -1;
  for (auto const &op : multisite_operators) {
    maxops = std::max(op.group, maxops);
  }
  num_multisite_operators = maxops + 1;
  multisite_operator_names.resize(num_multisite_operators);
  multisite_operator_counts.resize(num_multisite_operators);
  multisite_operator_nsites.resize(num_multisite_operators);
  for (auto const &op : multisite_operators) {
    ++multisite_operator_counts[op.group];
    std::string name = op.name;
    if (name.empty()) {
      std::stringstream ss;
      ss << "multisite[" << op.group << "]";
      name = ss.str();
    }
    if (multisite_operator_names[op.group].empty()) {
      multisite_operator_names[op.group] = name;
    } else if (multisite_operator_names[op.group] != name) {
      std::stringstream msg;
      msg << "ERROR: multisite operator " << name << " has group " << op.group
          << " but was previously set to " << multisite_operator_names[op.group]
          << "\n";
      msg << "       multisite operators with the same group must "
             "have the same name.\n";
      throw std::runtime_error(msg.str());
    }

    if (multisite_operator_nsites[op.group] == 0) {
      multisite_operator_nsites[op.group] = op.nsites();
    } else if (multisite_operator_nsites[op.group] != op.nsites()) {
      std::stringstream msg;
      msg << "ERROR: multisite operator " << op.name << " has nsites "
          << op.nsites() << " but was previously set to "
          << multisite_operator_nsites[op.group] << "\n";
      msg << "       Multiple multisite operators with the same group must "
             "have the same number of sites.\n";
      throw std::runtime_error(msg.str());
    }
    multisite_operator_nsites_set.insert(op.nsites());
  }
  for (auto const &s : multisite_operator_names) {
    maxlength = std::max(s.size(), maxlength);
  }

  // add right margin to operator names
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
  for (auto &s : multisite_operator_names) {
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
  core::Calc_CTM_Environment(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
                             peps_parameters, lattice);
  time_environment += timer.elapsed();
}

template <class ptensor>
void iTPS<ptensor>::update_CTM_density() {
  Timer<> timer;
  core::Calc_CTM_Environment_density(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
                             peps_parameters, lattice);
  time_environment += timer.elapsed();
}

template <class ptensor>
std::vector<ptensor> iTPS<ptensor>::make_single_tensor_density(){
  return core::Make_single_tensor_density(Tn);
}

  
// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
