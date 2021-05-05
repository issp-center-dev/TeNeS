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

#include "core/simple_update.hpp"
#include "core/local_gauge.hpp"

namespace tenes {
namespace itps {

template <class tensor>
void iTPS<tensor>::simple_update() {
  Timer<> timer;
  tensor Tn1_new(comm);
  tensor Tn2_new(comm);
  std::vector<double> lambda_c;
  const int nsteps = peps_parameters.num_simple_step;
  double next_report = 10.0;

  for (int int_tau = 0; int_tau < nsteps; ++int_tau) {
    for (auto up : simple_updates) {
      const int source = up.source_site;
      const int source_leg = up.source_leg;
      const int target = lattice.neighbor(source, source_leg);
      const int target_leg = (source_leg + 2) % 4;
      core::Simple_update_bond(Tn[source], Tn[target], lambda_tensor[source],
                               lambda_tensor[target], up.op, source_leg,
                               peps_parameters, Tn1_new, Tn2_new, lambda_c);
      lambda_tensor[source][source_leg] = lambda_c;
      lambda_tensor[target][target_leg] = lambda_c;
      Tn[source] = Tn1_new;
      Tn[target] = Tn2_new;
    }

    // local gauge fixing
    const int maxiter_gauge = peps_parameters.Simple_Gauge_maxiter;
    const double conv_tol_gauge =
        peps_parameters.Simple_Gauge_Convergence_Epsilon;
    if (peps_parameters.Simple_Gauge_Fix) {
      int iter_gauge = 0;
      for (iter_gauge = 0; iter_gauge < maxiter_gauge; ++iter_gauge) {
        for (int parity : {1, -1}) {
          for (int source_leg : {2, 1}) {
            int target_leg = (source_leg + 2) % 4;
            for (int source = 0; source < N_UNIT; ++source) {
              if (lattice.parity(source) != parity) {
                continue;
              }
              if (lattice.virtual_dims[source][source_leg] <= 1) {
                continue;
              }
              int target = lattice.neighbor(source, source_leg);
              core::fix_local_gauge(
                  Tn[source], Tn[target], lambda_tensor[source],
                  lambda_tensor[target], source_leg, peps_parameters, Tn1_new,
                  Tn2_new, lambda_c);
              lambda_tensor[source][source_leg] = lambda_c;
              lambda_tensor[target][target_leg] = lambda_c;
              Tn[source] = Tn1_new;
              Tn[target] = Tn2_new;
            }
          }
        }

        // convergence check
        double score = 0.0;
        for (int site = 0; site < N_UNIT; ++site) {
          for (int leg = 0; leg < nleg; ++leg) {
            if (lattice.virtual_dims[site][leg] <= 1) {
              continue;
            }
            auto M = core::boundary_tensor(Tn[site], lambda_tensor[site], leg,
                                           peps_parameters);
            tensor U;
            std::vector<double> D;
            eigh(M, mptensor::Axes(0), mptensor::Axes(1), D, U);
            for (auto d : D) {
              score = std::max(score, std::abs(d - 1.0));
            }
          }
        }  // end of for (source)
        // std::cout << int_tau << " " << iter_gauge << " " << score
        //           << std::endl;
        if (score < conv_tol_gauge) {
          // std::cout << int_tau << " " << iter_gauge << " " << score
          //           << std::endl;
          break;
        }
      }  // end of for (iter_gauge)
    }    // end of if (Simple_Gauge_Fix)

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
  }
  time_simple_update += timer.elapsed();
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
