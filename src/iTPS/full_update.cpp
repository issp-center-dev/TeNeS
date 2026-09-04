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
#include "core/full_update_fermion.hpp"
#include "core/ctm.hpp"
#include "../fermion/fops.hpp"
#include "../tensor.hpp"

namespace tenes::itps {

template <class tensor>
void iTPS<tensor>::full_update(EvolutionOperator<tensor> const &up) {
  if (up.is_onesite()) {
    const int source = up.source_site;
    if (finfo.enabled) {
      apply_onesite_gate_fermion(up);
      return;
    }
    Tn[source] =
        tensordot(Tn[source], up.op, mptensor::Axes(4), mptensor::Axes(0));
  } else {
    tensor Tn1_work(comm), Tn2_work(comm);
    const int source = up.source_site;
    const int source_leg = up.source_leg;
    const int target = lattice.neighbor(source, source_leg);
    const int target_leg = (source_leg + 2) % 4;

    if (finfo.enabled) {
      int s1 = source;
      int s2 = target;
      int s1_leg = source_leg;
      int s2_leg = target_leg;
      auto fop = tenes::fermion::wrap_twosite_gate(up.op, finfo.phys[source],
                                                   finfo.phys[target]);
      if (source_leg == 0 || source_leg == 1) {
        std::swap(s1, s2);
        std::swap(s1_leg, s2_leg);
        fop = tenes::fermion::transpose(fop, mptensor::Axes(1, 0, 3, 2));
      }
      auto fTn1 = tenes::fermion::wrap_Tn(Tn[s1], finfo, s1);
      auto fTn2 = tenes::fermion::wrap_Tn(Tn[s2], finfo, s2);
      tenes::fermion::ftensor<tensor> fTn1_work, fTn2_work;
      if (s1_leg == 2) {
        core::Full_update_bond_fermion(
            C1[s1], C2[s2], C3[s2], C4[s1], eTt[s1], eTt[s2], eTr[s2],
            eTb[s2], eTb[s1], eTl[s1], fTn1, fTn2, fop,
            tenes::fermion::reduced_pair_direction::horizontal,
            peps_parameters, fTn1_work, fTn2_work);
      } else {
        core::Full_update_bond_fermion(
            C1[s1], C2[s1], C3[s2], C4[s2], eTt[s1], eTr[s1], eTr[s2],
            eTb[s2], eTl[s2], eTl[s1], fTn1, fTn2, fop,
            tenes::fermion::reduced_pair_direction::vertical,
            peps_parameters, fTn1_work, fTn2_work);
      }
      finfo.virt[s1][s1_leg] = fTn1_work.parity[s1_leg];
      finfo.virt[s2][s2_leg] = fTn2_work.parity[s2_leg];
      tenes::fermion::unwrap_Tn(fTn1_work, Tn[s1], finfo, s1);
      tenes::fermion::unwrap_Tn(fTn2_work, Tn[s2], finfo, s2);
      tenes::fermion::validate_neighbor_consistency(finfo, lattice);
      update_CTM();
      return;
    }

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
                             eTb[source], eTb[target], eTl[target], eTt[target],
                             eTt[source], eTr[source], Tn[source], Tn[target],
                             up.op, source_leg, peps_parameters, Tn1_work,
                             Tn2_work);
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
                             eTl[source], eTl[target], eTt[target], eTr[target],
                             eTr[source], eTb[source], Tn[source], Tn[target],
                             up.op, source_leg, peps_parameters, Tn1_work,
                             Tn2_work);
    } else if (source_leg == 2) {
      /*
       *  C1 t t' C2'
       *  l  T T' r'
       *  C4 b b' C3'
       */
      core::Full_update_bond(C1[source], C2[target], C3[target], C4[source],
                             eTt[source], eTt[target],
                             eTr[target],  // t  t' r'
                             eTb[target], eTb[source],
                             eTl[source],  // b' b l
                             Tn[source], Tn[target], up.op, source_leg,
                             peps_parameters, Tn1_work, Tn2_work);
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
                             eTr[source], eTr[target], eTb[target], eTl[target],
                             eTl[source], eTt[source], Tn[source], Tn[target],
                             up.op, source_leg, peps_parameters, Tn1_work,
                             Tn2_work);
    }
    Tn[source] = Tn1_work;
    Tn[target] = Tn2_work;

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
}

template <class ptensor>
void iTPS<ptensor>::full_update() {
  const int group = 0;
  if (peps_parameters.num_full_step[group] > 0) {
    update_CTM();
  }

  Timer<> timer;
  ptensor Tn1_work(comm), Tn2_work(comm);
  const int nsteps = peps_parameters.num_full_step[group];
  double next_report = 10.0;

  for (int int_tau = 0; int_tau < nsteps; ++int_tau) {
    for (auto up : full_updates) {
      if (up.group != group) {
        continue;
      }
      full_update(up);
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

}  // namespace tenes::itps
