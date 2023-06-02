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
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "iTPS.hpp"

#include "../tensor.hpp"

#include "core/contract_ctm.hpp"
#include "core/contract_mf.hpp"

namespace tenes {
namespace itps {

template <class tensor>
auto iTPS<tensor>::measure_onesite()
    -> std::vector<std::vector<typename iTPS<tensor>::tensor_type>> {
  Timer<> timer;
  const int nlops = num_onesite_operators;
  std::vector<std::vector<tensor_type>> local_obs(
      nlops, std::vector<tensor_type>(
                 N_UNIT, std::numeric_limits<double>::quiet_NaN()));
  std::vector<tensor_type> norm(N_UNIT);

  if (peps_parameters.MeanField_Env) {
    std::vector<tensor> Tn_(Tn);
    for (int i = 0; i < N_UNIT; ++i) {
      for (int leg = 0; leg < nleg; ++leg) {
        const std::vector<double> mf = lambda_tensor[i][leg];
        Tn_[i].multiply_vector(mf, leg);
      }
    }

    for (int i = 0; i < N_UNIT; ++i) {
      norm[i] = core::Contract_one_site_MF(Tn_[i], op_identity[i]);
    }

    for (auto const &op : onesite_operators) {
      const int i = op.source_site;
      const auto val = core::Contract_one_site_MF(Tn_[i], op.op);
      local_obs[op.group][i] = val / norm[i];
    }
  } else {
    for (int i = 0; i < N_UNIT; ++i) {
      norm[i] =
          core::Contract_one_site(C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i],
                                  eTb[i], eTl[i], Tn[i], op_identity[i]);
    }
    for (auto const &op : onesite_operators) {
      const int i = op.source_site;
      const auto val =
          core::Contract_one_site(C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i],
                                  eTb[i], eTl[i], Tn[i], op.op);
      local_obs[op.group][i] = val / norm[i];
    }
  }
  double norm_real_min = 1e100;
  double norm_imag_abs_max = 0.0;
  for (int i = 0; i < N_UNIT; ++i) {
    norm_real_min = std::min(std::real(norm[i]), norm_real_min);
    norm_imag_abs_max =
        std::max(std::abs(std::imag(norm[i])), norm_imag_abs_max);
  }
  if (mpirank == 0) {
    if (norm_real_min <= 0.0) {
      std::cerr << "WARNING: Norm is negative [min(real(NORM)) = "
                << norm_real_min << "].\n";
      std::cerr << "HINT: Increase the bond dimension of CTM." << std::endl;
    }
    if (norm_imag_abs_max > 1e-6) {
      std::cerr << "WARNING: Norm is not real [max(abs(imag(NORM))) = "
                << norm_imag_abs_max << " > 1e-6].\n";
      std::cerr << "HINT: Increase the bond dimension of CTM." << std::endl;
    }
  }
  local_obs.push_back(norm);

  time_observable += timer.elapsed();
  return local_obs;
}

template <class ptensor>
void iTPS<ptensor>::save_onesite(
    std::vector<std::vector<typename iTPS<ptensor>::tensor_type>> const
        &onesite_obs,
    boost::optional<double> time, std::string filename_prefix) {
  if (mpirank != 0) {
    return;
  }

  const int nlops = num_onesite_operators;
  std::string filepath = outdir + "/" + filename_prefix + "onesite_obs.dat";
  if (!time && peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "    Save onesite observables to " << filepath << std::endl;
  }

  static bool first_time = true;
  if(first_time){
    std::ofstream ofs(filepath.c_str());
    int index = 1;
    if (time) {
      ofs << "# $" << index++ << ": (imaginary) time\n";
    }
    ofs << "# $" << index++ << ": op_group\n";
    ofs << "# $" << index++ << ": site_index\n";
    ofs << "# $" << index++ << ": real\n";
    ofs << "# $" << index++ << ": imag\n";
    ofs << std::endl;
    first_time = false;
  }

  std::ofstream ofs(filepath.c_str(), std::ios::out | std::ios::app);
  ofs << std::scientific
      << std::setprecision(std::numeric_limits<double>::max_digits10);

  for (int ilops = 0; ilops < nlops; ++ilops) {
    for (int i = 0; i < N_UNIT; ++i) {
      const auto v = onesite_obs[ilops][i];
      if (std::isnan(std::real(v))) {
        continue;
      }
      if (time) {
        ofs << time.get() << " ";
      }
      ofs << ilops << " " << i << " " << std::real(v) << " " << std::imag(v)
          << std::endl;
    }
  }
  if (onesite_obs.size() == nlops + 1) {
    // includes norm
    for (int i = 0; i < N_UNIT; ++i) {
      const auto v = onesite_obs[nlops][i];
      if (std::isnan(std::real(v))) {
        continue;
      }
      if (time) {
        ofs << time.get() << " ";
      }
      ofs << "-1 " << i << " " << std::real(v) << " " << std::imag(v)
          << std::endl;
    }
  }
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
