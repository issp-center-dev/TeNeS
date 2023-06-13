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
          if (onesite_operator_counts[ilops] == 0){
            continue;
          }
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
          if (twosite_operator_counts[ilops] == 0){
            continue;
          }
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
      ofs << "finish_datetime = " << util::datetime() << std::endl;
    }

    if (peps_parameters.print_level >= PrintLevel::info) {
      const double invV = 1.0 / numsites;
      std::cout << std::endl;

      std::cout << "Onesite observables per site:" << std::endl;
      for (int ilops = 0; ilops < num_onesite_operators; ++ilops) {
        if ( onesite_operator_counts[ilops] == 0) {
          continue;
        }
        const auto v = loc_obs[ilops] * invV;
        std::cout << "  " << onesite_operator_names[ilops] << " = "
                  << std::real(v) << " " << std::imag(v) << std::endl;
      }

      std::cout << "Twosite observables per site:" << std::endl;
      for (int ilops = 0; ilops < num_twosite_operators; ++ilops) {
        if ( twosite_operator_counts[ilops] == 0) {
          continue;
        }
        const auto v = two_obs[ilops] * invV;
        std::cout << "  " << twosite_operator_names[ilops] << " = "
                  << std::real(v) << " " << std::imag(v) << std::endl;
      }
    }
  }  // end of if(mpirank == 0)
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

/*
template <class ptensor>
void iTPS<ptensor>::measure_density() {
  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "Start calculating observables" << std::endl;
    std::cout << "  Start updating environment" << std::endl;
  }
  if (!peps_parameters.MeanField_Env) {
    update_CTM_density();
  }

  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "  Start calculating onesite operators" << std::endl;
  }
  auto onesite_obs = measure_onesite_density();
  save_onesite(onesite_obs);

  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "  Start calculating twosite operators" << std::endl;
  }
  auto twosite_obs = measure_twosite_density();
  save_twosite(twosite_obs);

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
          if (onesite_operator_counts[ilops] == 0){
            continue;
          }
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
          if (twosite_operator_counts[ilops] == 0){
            continue;
          }
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
      ofs << "finish_datetime = " << util::datetime() << std::endl;
    }

    if (peps_parameters.print_level >= PrintLevel::info) {
      const double invV = 1.0 / numsites;
      std::cout << std::endl;

      std::cout << "Onesite observables per site:" << std::endl;
      for (int ilops = 0; ilops < num_onesite_operators; ++ilops) {
        if ( onesite_operator_counts[ilops] == 0) {
          continue;
        }
        const auto v = loc_obs[ilops] * invV;
        std::cout << "  " << onesite_operator_names[ilops] << " = "
                  << std::real(v) << " " << std::imag(v) << std::endl;
      }

      std::cout << "Twosite observables per site:" << std::endl;
      for (int ilops = 0; ilops < num_twosite_operators; ++ilops) {
        if ( twosite_operator_counts[ilops] == 0) {
          continue;
        }
        const auto v = two_obs[ilops] * invV;
        std::cout << "  " << twosite_operator_names[ilops] << " = "
                  << std::real(v) << " " << std::imag(v) << std::endl;
      }
    }
  }  // end of if(mpirank == 0)
}
 */
  
// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
