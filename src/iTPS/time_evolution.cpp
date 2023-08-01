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
void iTPS<tensor>::time_evolution() {
  Timer<> timer;
  double t = 0.0;
  measure(t, "TE_");

  const bool su = peps_parameters.num_simple_step[0] > 0;
  const int num_groups = su ? num_simple_update_groups : num_full_update_groups;

  if (peps_parameters.measure_interval.size() < num_groups){
    while(peps_parameters.measure_interval.size() < num_groups){
      peps_parameters.measure_interval.push_back(*peps_parameters.measure_interval.rbegin());
    }
  }

  for (int group = 0; group < num_groups; ++group) {
    const int nsteps = su ? peps_parameters.num_simple_step[group]
                          : peps_parameters.num_full_step[group];
    const double dt = su ? peps_parameters.tau_simple_step[group]
                         : peps_parameters.tau_full_step[group];
    int measure_interval = peps_parameters.measure_interval[group];
    double time_measure = 0.0;
    double next_report = 10.0;

    for (int int_tau = 0; int_tau < nsteps; ++int_tau) {
      auto const& evols = su ? simple_updates : full_updates;
      for (auto up : evols) {
        if (up.group != group) {
          continue;
        }
        if(su){
          simple_update(up);
        }else{
          full_update(up);
        }
      }
      t += dt;
      // local gauge fixing
      if (su && peps_parameters.Simple_Gauge_Fix) {
        fix_local_gauge();
      }
      if ((int_tau + 1) % measure_interval == 0) {
        Timer<> timer_m;
        measure(t, "TE_");
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
      measure(t, "TE_");
      time_measure += timer_m.elapsed();
    }
    if(su){
      time_simple_update += timer.elapsed() - time_measure;
    }else{
      time_full_update += timer.elapsed() - time_measure;
    }
  }
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
