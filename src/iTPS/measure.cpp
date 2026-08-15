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
#include "../fermion/fops.hpp"
#include "../util/datetime.hpp"
#include "../timer.hpp"
#include "../version.hpp"
#include "core/ctm.hpp"
#ifndef _NO_OMP
#include <omp.h>
#endif

namespace tenes::itps {

template <class ptensor>
void iTPS<ptensor>::measure(std::optional<double> time,
                            std::string filename_prefix) {
  if (!time && peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "Start calculating observables" << std::endl;
    std::cout << "  Start updating environment" << std::endl;
  }

  if (!peps_parameters.MeanField_Env && finfo.enabled) {
    Timer<> timer;
    using ftensor = tenes::fermion::ftensor<ptensor>;
    std::vector<ftensor> fC1, fC2, fC3, fC4, feTt, feTr, feTb, feTl, fTn;
    fC1.reserve(N_UNIT);
    fC2.reserve(N_UNIT);
    fC3.reserve(N_UNIT);
    fC4.reserve(N_UNIT);
    feTt.reserve(N_UNIT);
    feTr.reserve(N_UNIT);
    feTb.reserve(N_UNIT);
    feTl.reserve(N_UNIT);
    fTn.reserve(N_UNIT);
    for (int site = 0; site < N_UNIT; ++site) {
      fC1.push_back(tenes::fermion::wrap_C(C1[site], finfo, 0, site));
      fC2.push_back(tenes::fermion::wrap_C(C2[site], finfo, 1, site));
      fC3.push_back(tenes::fermion::wrap_C(C3[site], finfo, 2, site));
      fC4.push_back(tenes::fermion::wrap_C(C4[site], finfo, 3, site));
      feTt.push_back(tenes::fermion::wrap_eT(eTt[site], finfo, 0, site));
      feTr.push_back(tenes::fermion::wrap_eT(eTr[site], finfo, 1, site));
      feTb.push_back(tenes::fermion::wrap_eT(eTb[site], finfo, 2, site));
      feTl.push_back(tenes::fermion::wrap_eT(eTl[site], finfo, 3, site));
      fTn.push_back(tenes::fermion::wrap_Tn(Tn[site], finfo, site));
    }
    core::Calc_CTM_Environment(fC1, fC2, fC3, fC4, feTt, feTr, feTb, feTl, fTn,
                               peps_parameters, lattice);
    for (int site = 0; site < N_UNIT; ++site) {
      tenes::fermion::unwrap_C(fC1[site], C1[site], finfo, 0, site);
      tenes::fermion::unwrap_C(fC2[site], C2[site], finfo, 1, site);
      tenes::fermion::unwrap_C(fC3[site], C3[site], finfo, 2, site);
      tenes::fermion::unwrap_C(fC4[site], C4[site], finfo, 3, site);
      tenes::fermion::unwrap_eT(feTt[site], eTt[site], finfo, 0, site);
      tenes::fermion::unwrap_eT(feTr[site], eTr[site], finfo, 1, site);
      tenes::fermion::unwrap_eT(feTb[site], eTb[site], finfo, 2, site);
      tenes::fermion::unwrap_eT(feTl[site], eTl[site], finfo, 3, site);
    }
    time_environment += timer.elapsed();
  } else if (!peps_parameters.MeanField_Env) {
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

  if (finfo.enabled && tmatrix_param.to_calculate) {
    if (mpirank == 0) {
      std::cerr << "WARNING: fermion mode disables correlation_length.measure "
                   "because the transfer-matrix correlation length is not "
                   "fermion-aware in M1"
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
