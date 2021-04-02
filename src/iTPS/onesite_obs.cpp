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

#include "iTPS.hpp"

#include "../tensor.hpp"

#include "core/contract_ctm.hpp"
#include "core/contract_mf.hpp"

namespace tenes {

template <class tensor>
auto iTPS<tensor>::measure_onesite()
    -> std::vector<std::vector<typename iTPS<tensor>::tensor_type>> {
  Timer<> timer;
  const int nlops = num_onesite_operators;
  std::vector<std::vector<tensor_type>> local_obs(
      nlops, std::vector<tensor_type>(
                 N_UNIT, std::numeric_limits<double>::quiet_NaN()));

  if (peps_parameters.MeanField_Env) {
    std::vector<tensor> Tn_(Tn);
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
void iTPS<ptensor>::save_onesite(
    std::vector<std::vector<typename iTPS<ptensor>::tensor_type>> const
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

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // end of namespace tenes
