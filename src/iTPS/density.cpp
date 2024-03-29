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
void iTPS<ptensor>::save_density(
    std::vector<std::vector<tensor_type>> const &onesite_obs,
    std::vector<std::map<Bond, tensor_type>> const &twosite_obs,
    std::vector<std::map<Multisites, tensor_type>> const &multisite_obs,
    boost::optional<double> time, std::string filename_prefix) {
  if (mpirank != 0) {
    return;
  }
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
  std::vector<tensor_type> multi_obs(num_multisite_operators);
  for (int iops = 0; iops < num_multisite_operators; ++iops) {
    for (const auto &obs : multisite_obs[iops]) {
      multi_obs[iops] += obs.second;
    }
  }

  auto energy = 0.0;
  for (int iops=0; iops<num_onesite_operators; ++iops) {
    if (onesite_operator_names[iops] == "site_hamiltonian"){
      energy += std::real(loc_obs[iops]);
    }
  }
  for (int iops=0; iops<num_twosite_operators; ++iops) {
    if (twosite_operator_names[iops] == "bond_hamiltonian"){
      energy += std::real(two_obs[iops]);
    }
  }

  const double invV = 1.0 / numsites;
  std::string filename = outdir + "/" + filename_prefix + "density.dat";
  energy *= invV;

  std::ios::openmode mode = std::ios::out;
  if (time) {
    static bool first_time = true;
    if (first_time) {
      first_time = false;
      std::ofstream ofs(filename.c_str());
      ofs << "# The meaning of each column is the following: \n";
      if (peps_parameters.calcmode ==
          PEPS_Parameters::CalculationMode::time_evolution) {
        ofs << "# $1: time\n";
      } else if (peps_parameters.calcmode ==
                 PEPS_Parameters::CalculationMode::finite_temperature) {
        ofs << "# $1: inverse temperature\n";
      }
      ofs << "# $2: observable ID\n";
      ofs << "# $3: real\n";
      ofs << "# $4: imag\n";

      int index = 0;
      ofs << "# The meaning of observable IDs are the following: \n";
      ofs << "# " << index++ << ": Energy\n";
      for (int ilops = 0; ilops < num_onesite_operators; ++ilops) {
        if (onesite_operator_counts[ilops] == 0) {
          continue;
        }
        ofs << "# " << index++ << ": " << onesite_operator_names[ilops] << "\n";
      }
      for (int ilops = 0; ilops < num_twosite_operators; ++ilops) {
        if (twosite_operator_counts[ilops] == 0) {
          continue;
        }
        ofs << "# " << index++ << ": " << twosite_operator_names[ilops] << "\n";
      }
      for (int ilops = 0; ilops < num_multisite_operators; ++ilops) {
        if (multisite_operator_counts[ilops] == 0) {
          continue;
        }
        ofs << "# " << index++ << ": " << multisite_operator_names[ilops] << "\n";
      }
      ofs << std::endl;
    }
    mode |= std::ios::app;
  }
  std::ofstream ofs(filename.c_str(), mode);
  ofs << std::scientific
      << std::setprecision(std::numeric_limits<double>::max_digits10);

  if (time) {
    ofs << time.get() << " 0 "; // 0 means the index of "Energy"
  } else {
    size_t namesize = onesite_operator_names[0].size();
    std::string s = "Energy";
    const auto l = namesize - s.size();
    for (size_t i = 0; i < l; ++i) {
      s += " ";
    }
    ofs << s << " = ";
  }
  if (energy >= 0.0) {
    ofs << " ";
  }
  ofs << energy << "  " << 0.0 << std::endl;

  int index = 1;
  for (int ilops = 0; ilops < num_onesite_operators; ++ilops) {
    if (onesite_operator_counts[ilops] == 0) {
      continue;
    }
    const auto v = loc_obs[ilops] * invV;
    if (time) {
      ofs << time.get() << " " << index++ << " ";
    } else {
      ofs << onesite_operator_names[ilops] << " = ";
    }
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
    if (twosite_operator_counts[ilops] == 0) {
      continue;
    }
    if (time) {
      ofs << time.get() << " " << index++ << " ";
    } else {
      ofs << twosite_operator_names[ilops] << " = ";
    }
    const auto v = two_obs[ilops] * invV;
    if (std::real(v) >= 0.0) {
      ofs << " ";
    }
    ofs << std::real(v) << " ";
    if (std::imag(v) >= 0.0) {
      ofs << " ";
    }
    ofs << std::imag(v) << std::endl;
  }

  for (int ilops = 0; ilops < num_multisite_operators; ++ilops) {
    if (multisite_operator_counts[ilops] == 0) {
      continue;
    }
    if (time) {
      ofs << time.get() << " " << index++ << " ";
    } else {
      ofs << multisite_operator_names[ilops] << " = ";
    }
    const auto v = multi_obs[ilops] * invV;
    if (std::real(v) >= 0.0) {
      ofs << " ";
    }
    ofs << std::real(v) << " ";
    if (std::imag(v) >= 0.0) {
      ofs << " ";
    }
    ofs << std::imag(v) << std::endl;
  }

  if (!time) {
    std::cout << "    Save observable densities to " << filename << std::endl;
    if (peps_parameters.print_level >= PrintLevel::info) {
      std::cout << std::endl;

      std::cout << "Onesite observables per site:" << std::endl;
      for (int ilops = 0; ilops < num_onesite_operators; ++ilops) {
        if (onesite_operator_counts[ilops] == 0) {
          continue;
        }
        const auto v = loc_obs[ilops] * invV;
        std::cout << "  " << onesite_operator_names[ilops] << " = "
                  << std::real(v) << " " << std::imag(v) << std::endl;
      }

      std::cout << "Twosite observables per site:" << std::endl;
      for (int ilops = 0; ilops < num_twosite_operators; ++ilops) {
        if (twosite_operator_counts[ilops] == 0) {
          continue;
        }
        const auto v = two_obs[ilops] * invV;
        std::cout << "  " << twosite_operator_names[ilops] << " = "
                  << std::real(v) << " " << std::imag(v) << std::endl;
      }
    }
  }
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
