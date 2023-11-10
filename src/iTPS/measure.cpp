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
void iTPS<ptensor>::measure(boost::optional<double> time,
                            std::string filename_prefix) {
  if (!time && peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "Start calculating observables" << std::endl;
    std::cout << "  Start updating environment" << std::endl;
  }

  if (!peps_parameters.MeanField_Env) {
    update_CTM();
  }

  if (!time && peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "  Start calculating onesite operators" << std::endl;
  }
  auto onesite_obs = measure_onesite();
  save_onesite(onesite_obs, time, filename_prefix);

  if (!time && peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "  Start calculating twosite operators" << std::endl;
  }
  auto twosite_obs = measure_twosite();
  save_twosite(twosite_obs, time, filename_prefix);

  if (!time && peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "  Start calculating multisite operators" << std::endl;
  }
  auto multisite_obs = measure_multisite();
  if (multisite_operators.size() > 0) {
    save_multisite(multisite_obs, time, filename_prefix);
  }

  if (corparam.r_max > 0) {
    if (!time && peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "  Start calculating long range correlation" << std::endl;
    }
    auto correlations = measure_correlation();
    save_correlation(correlations, time, filename_prefix);
  }

  if (tmatrix_param.to_calculate) {
    if (!time && peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "  Start calculating correlation length" << std::endl;
    }
    auto correlation_length = measure_transfer_matrix_eigenvalues();
    save_correlation_length(correlation_length, time, filename_prefix);
  }

  save_density(onesite_obs, twosite_obs, multisite_obs, time, filename_prefix);
}

template <class ptensor>
void iTPS<ptensor>::summary() const {
  if (mpirank == 0) {
    const double time_all = timer_all.elapsed();
    {
      std::string filename = outdir + "/time.dat";
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


template <class ptensor>
void iTPS<ptensor>::measure_density(double beta, std::string filename_prefix) {
  //if (!peps_parameters.MeanField_Env) {
  update_CTM_density();
  //  }

  // DEBUG
#define DUMP_TENSOR(NAME) \
  { \
    int site = 0; \
    for (int n=0; n < NAME[site].local_size(); ++n){ \
      typename ptensor::value_type v; \
      mptensor::Index idx = NAME[site].global_index(n); \
      NAME[site].get_value(idx, v); \
      std::cerr << "DEBUG " #NAME " " << beta << " "; \
      for(int d=0; d<idx.size(); ++d){ \
        std::cerr << idx[d] << " "; \
      } \
      std::cerr << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << std::real(v) << " " << std::imag(v) << std::endl; \
    } \
  }

  DUMP_TENSOR(C1);
  DUMP_TENSOR(C2);
  DUMP_TENSOR(C3);
  DUMP_TENSOR(C4);
  DUMP_TENSOR(eTl);
  DUMP_TENSOR(eTt);
  DUMP_TENSOR(eTl);
  DUMP_TENSOR(eTb);

#undef DUMP_TENSOR

  // auto onesite_obs = measure_onesite_density();
  auto onesite_obs = measure_onesite();
  save_onesite(onesite_obs, beta, filename_prefix);

  // auto twosite_obs = measure_twosite_density();
  auto twosite_obs = measure_twosite();
  save_twosite(twosite_obs, beta, filename_prefix);

  // In finite temperature simplation, multisite operators are not supported so far.
  // auto multisite_obs = measure_multisite_density();
  auto multisite_obs = measure_multisite();
  if (multisite_operators.size() > 0) {
    save_multisite(multisite_obs, beta, filename_prefix);
  }

  if (corparam.r_max > 0) {
    auto correlations = measure_correlation();
    save_correlation(correlations, beta, filename_prefix);
  }

  if (tmatrix_param.to_calculate) {
    if(beta > 0.0){
      // the method is unstable at beta = 0.0, so we skip it.
      auto correlation_length = measure_transfer_matrix_eigenvalues();
      save_correlation_length(correlation_length, beta, filename_prefix);
    }
  }
  save_density(onesite_obs, twosite_obs, multisite_obs, beta, filename_prefix);
}
  
// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
