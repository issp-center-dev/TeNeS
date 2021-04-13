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

#include "../printlevel.hpp"
#include "../timer.hpp"

#include "core/full_update.hpp"
#include "core/ctm.hpp"
#include "../tensor.hpp"

namespace tenes {
namespace itps {

template <class ptensor>
void iTPS<ptensor>::full_update() {
  if (peps_parameters.num_full_step > 0) {
    update_CTM();
  }

  Timer<> timer;
  ptensor Tn1_new(comm), Tn2_new(comm);
  const int nsteps = peps_parameters.num_full_step;
  double next_report = 10.0;

  for (int int_tau = 0; int_tau < nsteps; ++int_tau) {
    for (auto up : full_updates) {
      const int source = up.source_site;
      const int source_leg = up.source_leg;
      const int target = lattice.neighbor(source, source_leg);
      // const int target_leg = (source_leg + 2) % 4;

      if (source_leg == 0) {
        /*
         *  C1' t' t C3
         *  l'  T' T r
         *  C2' b' b C4
         *
         *   |
         *   | rotate
         *   V
         *
         *  C4 b b' C2'
         *  r  T T' l'
         *  C3 t t' C1'
         */
        core::Full_update_bond(C4[source], C2[target], C1[target], C3[source],
                               eTb[source], eTb[target], eTl[target],
                               eTt[target], eTt[source], eTr[source],
                               Tn[source], Tn[target], up.op, source_leg,
                               peps_parameters, Tn1_new, Tn2_new);
      } else if (source_leg == 1) {
        /*
         * C1' t' C2'
         *  l' T'  r'
         *  l  T   r
         * C4  b  C3
         *
         *   |
         *   | rotate
         *   V
         *
         *  C4 l l' C1'
         *  b  T T' t'
         *  C3 r r' C2'
         */
        core::Full_update_bond(C4[source], C1[target], C2[target], C3[source],
                               eTl[source], eTl[target], eTt[target],
                               eTr[target], eTr[source], eTb[source],
                               Tn[source], Tn[target], up.op, source_leg,
                               peps_parameters, Tn1_new, Tn2_new);
      } else if (source_leg == 2) {
        /*
         *  C1 t t' C2'
         *  l  T T' r'
         *  C4 b b' C3'
         */
        core::Full_update_bond(C1[source], C2[target], C3[target], C4[source],
                               eTt[source], eTt[target],
                               eTr[target],  // t  t' r'
                               eTb[target], eTb[source], eTl[source],  // b' b l
                               Tn[source], Tn[target], up.op, source_leg,
                               peps_parameters, Tn1_new, Tn2_new);
      } else {
        /*
         * C1  t C2
         *  l  T  r
         *  l' T' r'
         * C4' b C3'
         *
         *   |
         *   | rotate
         *   V
         *
         *  C2 r r' C3'
         *  t  T T' b'
         *  C1 l l' C4'
         */
        core::Full_update_bond(C2[source], C3[target], C4[target], C1[source],
                               eTr[source], eTr[target], eTb[target],
                               eTl[target], eTl[source], eTt[source],
                               Tn[source], Tn[target], up.op, source_leg,
                               peps_parameters, Tn1_new, Tn2_new);
      }
      Tn[source] = Tn1_new;
      Tn[target] = Tn2_new;

      if (peps_parameters.Full_Use_FastFullUpdate) {
        if (source_leg == 0) {
          const int source_x = source % LX;
          const int target_x = target % LX;
          core::Right_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, source_x,
                           peps_parameters, lattice);
          core::Left_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, target_x,
                          peps_parameters, lattice);
        } else if (source_leg == 1) {
          const int source_y = source / LX;
          const int target_y = target / LX;
          core::Bottom_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, source_y,
                            peps_parameters, lattice);
          core::Top_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, target_y,
                         peps_parameters, lattice);
        } else if (source_leg == 2) {
          const int source_x = source % LX;
          const int target_x = target % LX;
          core::Left_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, source_x,
                          peps_parameters, lattice);
          core::Right_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, target_x,
                           peps_parameters, lattice);
        } else {
          const int source_y = source / LX;
          const int target_y = target / LX;
          core::Top_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, source_y,
                         peps_parameters, lattice);
          core::Bottom_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, target_y,
                            peps_parameters, lattice);
        }
      } else {
        update_CTM();
      }
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
  }
  time_full_update += timer.elapsed();
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
