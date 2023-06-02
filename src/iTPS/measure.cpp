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

#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "iTPS.hpp"
#include "../util/datetime.hpp"

namespace tenes {
namespace itps {

template <class ptensor>
void iTPS<ptensor>::measure() {
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

  if (tmatrix_param.to_calculate) {
    if (peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "  Start calculating correlation length" << std::endl;
    }
    auto correlation_length = measure_transfer_matrix_eigenvalues();
    save_correlation_length(correlation_length);
  }

  save_density(onesite_obs, twosite_obs);
}

template <class ptensor>
void iTPS<ptensor>::measure(double time, std::string filename_prefix) {
  if (!peps_parameters.MeanField_Env) {
    update_CTM();
  }

  auto onesite_obs = measure_onesite();
  save_onesite(onesite_obs, time, filename_prefix);

  auto twosite_obs = measure_twosite();
  save_twosite(twosite_obs, time, filename_prefix);

  if (corparam.r_max > 0) {
    auto correlations = measure_correlation();
    save_correlation(correlations, time, filename_prefix);
  }

  if (tmatrix_param.to_calculate) {
    auto correlation_length = measure_transfer_matrix_eigenvalues();
    save_correlation_length(correlation_length, time, filename_prefix);
  }

  save_density(onesite_obs, twosite_obs, time, filename_prefix);
}

template <class ptensor>
void iTPS<ptensor>::summary(std::string filename_prefix) const {
  if (mpirank == 0) {
    const double time_all = timer_all.elapsed();
    {
      std::string filename = outdir + "/" + filename_prefix + "time.dat";
      std::ofstream ofs(filename.c_str());
      ofs << "time all           = " << time_all << std::endl;
      ofs << "time simple update = " << time_simple_update << std::endl;
      ofs << "time full update   = " << time_full_update << std::endl;
      ofs << "time environment  = " << time_environment << std::endl;
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
      std::cout << "  environment  = " << time_environment << std::endl;
      std::cout << "  observable    = " << time_observable << std::endl;
      std::cout << std::endl << "Done." << std::endl;
    }
  }
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
