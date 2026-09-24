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
#include "../exception.hpp"
#include "../fermion/fops.hpp"
#include "../fermion/reduced_measure.hpp"
#include "../util/datetime.hpp"
#include "../timer.hpp"
#include "../version.hpp"
#include "core/ctm.hpp"
#ifndef _NO_OMP
#include <omp.h>
#endif

namespace tenes::itps {

template <class ptensor>
void iTPS<ptensor>::validate_fermion_ctm_measurement() const {
  if (!finfo.enabled || peps_parameters.MeanField_Env) {
    return;
  }

  bool has_non_nearest_twosite = false;
  for (const auto &op : twosite_operators) {
    const int abs_dx = std::abs(op.dx[0]);
    const int abs_dy = std::abs(op.dy[0]);
    const bool is_nearest_neighbor =
        (abs_dx == 1 && abs_dy == 0) || (abs_dx == 0 && abs_dy == 1);
    // A same-site pair falls through to the raw-Tn path in measure_twosite(),
    // which would mix it with the reduced CTM environment and silently return
    // an incorrect result.
    if (!is_nearest_neighbor) {
      has_non_nearest_twosite = true;
      break;
    }
  }
  if (has_non_nearest_twosite || !multisite_operators.empty() ||
      corparam.r_max > 0) {
    throw tenes::input_error(
        "fermion CTM measurement supports nearest-neighbor two-site "
        "observables only");
  }
}

template <class ptensor>
void iTPS<ptensor>::measure(std::optional<double> time,
                            std::string filename_prefix) {
  validate_fermion_ctm_measurement();

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
  auto onesite_rdm = measure_onesite_rdm();
  save_onesite_rdm(onesite_rdm, time, filename_prefix);

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

  if (finfo.enabled && tmatrix_param.to_calculate) {
    if (mpirank == 0) {
      std::cerr << "WARNING: fermion mode disables correlation_length.measure "
                   "because the transfer-matrix correlation length is not "
                   "fermion-aware in this version"
                << std::endl;
    }
    tmatrix_param.to_calculate = false;
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
  auto &registry = TimerRegistry::instance();
  registry.add("total", timer_all.elapsed());
  registry.add("phase/simple_update", time_simple_update);
  registry.add("phase/full_update", time_full_update);
  registry.add("phase/environment", time_environment);
  registry.add("phase/observable", time_observable);
#ifndef _NO_OMP
  const int omp_threads = omp_get_max_threads();
#else
  const int omp_threads = 1;
#endif
  const auto aggregated = aggregate_timers(registry, comm);
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
    {
      std::string filename = outdir + "/timers.json";
      std::ofstream ofs(filename.c_str());
      ofs << timers_to_json(aggregated, TENES_VERSION, mpisize, omp_threads);
      if (peps_parameters.print_level >= PrintLevel::info) {
        std::cout << "    Save timers to " << filename << std::endl;
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
  // if (!peps_parameters.MeanField_Env) {
  update_CTM_density();
  //  }

  // auto onesite_obs = measure_onesite_density();
  auto onesite_obs = measure_onesite();
  save_onesite(onesite_obs, beta, filename_prefix);
  auto onesite_rdm = measure_onesite_rdm();
  save_onesite_rdm(onesite_rdm, beta, filename_prefix);

  // auto twosite_obs = measure_twosite_density();
  auto twosite_obs = measure_twosite();
  save_twosite(twosite_obs, beta, filename_prefix);

  // In finite temperature simplation, multisite operators are not supported so
  // far. auto multisite_obs = measure_multisite_density();
  auto multisite_obs = measure_multisite();
  if (multisite_operators.size() > 0) {
    save_multisite(multisite_obs, beta, filename_prefix);
  }

  if (corparam.r_max > 0) {
    auto correlations = measure_correlation();
    save_correlation(correlations, beta, filename_prefix);
  }

  if (tmatrix_param.to_calculate) {
    if (beta > 0.0) {
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

}  // namespace tenes::itps
