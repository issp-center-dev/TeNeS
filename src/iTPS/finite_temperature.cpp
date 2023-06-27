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

namespace tenes {
namespace itps {

template <class tensor>
void iTPS<tensor>::finite_temperature() {
  Timer<> timer;
  double time_measure = 0.0;
  double next_report = 10.0;
  double beta = 0.0;
  measure_density(beta, "FT_");

  const int num_groups = num_simple_update_groups;

  for (int group = 0; group < num_groups; ++group) {
    const int nsteps = peps_parameters.num_simple_step[group];
    const double dt = peps_parameters.tau_simple_step[group];
    int measure_interval = peps_parameters.measure_interval[group];

    for (int int_tau = 0; int_tau < nsteps; ++int_tau) {
      auto const& evols = simple_updates;
      for (auto up : evols) {
        if (up.group != group) {
          continue;
        }
	simple_update_density(up);
      }
      beta += dt * 2.0;
      // local gauge fixing
      if (peps_parameters.Simple_Gauge_Fix) {
        fix_local_gauge_density();
      }
      if ((int_tau + 1) % measure_interval == 0) {
        Timer<> timer_m;
        measure_density(beta, "FT_");
        time_measure += timer_m.elapsed();
      }
      if (peps_parameters.print_level >= PrintLevel::info) {
        double r_tau = 100.0 * (int_tau + 1) / nsteps;
        if (r_tau >= next_report) {
          while (r_tau >= next_report) {
            next_report += 10.0;
          }
          std::cout << "  " << next_report - 10.0 << "% [" << int_tau + 1 << "/"
                    << nsteps << "] done" << std::endl;
        }
      }
    }  // end of for (int_tau)
    if (nsteps % measure_interval != 0) {
      Timer<> timer_m;
      measure_density(beta, "TE_");
      time_measure += timer_m.elapsed();
    }
    time_simple_update += timer.elapsed() - time_measure;
  }
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
