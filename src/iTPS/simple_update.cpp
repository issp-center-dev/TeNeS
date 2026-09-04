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

#include <cstdlib>
#include <iostream>

#include "../fermion/fops.hpp"
#include "core/simple_update.hpp"
#include "core/local_gauge.hpp"

namespace tenes::itps {

namespace {
// Diagnostic knob: setting TENES_FERMION_SECTOR_LOG (any value) dumps the
// even/odd dimension split of every virtual-bond ledger to stderr after
// each simple-update step, to watch how the truncation distributes the
// bond dimension over the two parity sectors.
bool fermion_sector_log_enabled() {
  return std::getenv("TENES_FERMION_SECTOR_LOG") != nullptr;
}

void log_fermion_sector_dimensions(const tenes::fermion::FermionInfo& finfo,
                                   const SquareLattice& lattice, int n_unit,
                                   int mpirank, int step) {
  if (!fermion_sector_log_enabled() || !finfo.enabled || mpirank != 0) {
    return;
  }
  for (int source = 0; source < n_unit; ++source) {
    for (int leg : {2, 1}) {
      const int target = lattice.neighbor(source, leg);
      const auto& parity = finfo.virt[source][leg];
      int even = 0;
      int odd = 0;
      for (bool p : parity) {
        if (p) {
          ++odd;
        } else {
          ++even;
        }
      }
      std::cerr << "TENES_FERMION_SECTOR step=" << step << " source=" << source
                << " leg=" << leg << " target=" << target << " Deven=" << even
                << " Dodd=" << odd << "\n";
    }
  }
}
}  // namespace

template <class tensor>
void iTPS<tensor>::apply_onesite_gate_fermion(
    EvolutionOperator<tensor> const& up) {
  const int source = up.source_site;
  auto fTn = tenes::fermion::wrap_Tn(Tn[source], finfo, source);
  tenes::fermion::ftensor<tensor> fop{
      up.op, {finfo.phys[source], finfo.phys[source]}};
  auto updated = tenes::fermion::tensordot(fTn, fop, mptensor::Axes(4),
                                           mptensor::Axes(0));
  tenes::fermion::unwrap_Tn(updated, Tn[source], finfo, source);
}

template <class tensor>
void iTPS<tensor>::simple_update(EvolutionOperator<tensor> const& up) {
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
    std::vector<double> lambda_work;
    const int source = up.source_site;
    const int source_leg = up.source_leg;
    const int target = lattice.neighbor(source, source_leg);
    const int target_leg = (source_leg + 2) % 4;
    if (finfo.enabled) {
      // Canonical fermionic bond orientation is raster order (left-to-right,
      // top-to-bottom): the source of the graded gate must be the JW-earlier
      // site, i.e. the left site of a horizontal bond (source_leg == 2) or
      // the upper site of a vertical bond (source_leg == 3). Updates given
      // in the opposite orientation are normalized by swapping the site
      // roles; feeding the gate from the JW-later side applies wrong Koszul
      // masks (verified by the vertical-chain lambda trajectory diagnostic).
      int s1 = source;
      int s2 = target;
      int s1_leg = source_leg;
      int s2_leg = target_leg;
      auto fop = tenes::fermion::wrap_twosite_gate(up.op, finfo.phys[source],
                                                   finfo.phys[target]);
      if (source_leg == 0 || source_leg == 1) {
        std::swap(s1, s2);
        std::swap(s1_leg, s2_leg);
        // Graded transpose: exchanging the two sites of a two-site operator
        // conjugates it by S(|a> x |b>) = (-1)^{|a||b|} |b> x |a>, so it
        // carries the Fock reordering sign on both leg pairs, exactly as the
        // measurement path does (twosite_obs.cpp). A plain transpose drops
        // (-1)^{p(j1)p(j2) + p(k1)p(k2)}, which is -1 on the
        // (odd,odd) <-> (even,even) channels: empty for a
        // particle-number-conserving spinless model, but the doublon-holon
        // channel of an electron model at first order in the hopping.
        fop = tenes::fermion::transpose(fop, mptensor::Axes(1, 0, 3, 2));
      }
      auto fTn1 = tenes::fermion::wrap_Tn(Tn[s1], finfo, s1);
      auto fTn2 = tenes::fermion::wrap_Tn(Tn[s2], finfo, s2);
      tenes::fermion::ftensor<tensor> fTn1_work, fTn2_work;
      core::Simple_update_bond(fTn1, fTn2, lambda_tensor[s1], lambda_tensor[s2],
                               fop, s1_leg, peps_parameters, fTn1_work,
                               fTn2_work, lambda_work);
      lambda_tensor[s1][s1_leg] = lambda_work;
      lambda_tensor[s2][s2_leg] = lambda_work;
      finfo.virt[s1][s1_leg] = fTn1_work.parity[s1_leg];
      finfo.virt[s2][s2_leg] = fTn2_work.parity[s2_leg];
      tenes::fermion::unwrap_Tn(fTn1_work, Tn[s1], finfo, s1);
      tenes::fermion::unwrap_Tn(fTn2_work, Tn[s2], finfo, s2);
      return;
    }
    core::Simple_update_bond(Tn[source], Tn[target], lambda_tensor[source],
                             lambda_tensor[target], up.op, source_leg,
                             peps_parameters, Tn1_work, Tn2_work, lambda_work);
    lambda_tensor[source][source_leg] = lambda_work;
    lambda_tensor[target][target_leg] = lambda_work;
    Tn[source] = Tn1_work;
    Tn[target] = Tn2_work;
  }
}

template <class tensor>
void iTPS<tensor>::fix_local_gauge() {
  tensor Tn1_work(comm), Tn2_work(comm);
  std::vector<double> lambda_work;
  const int maxiter_gauge = peps_parameters.Simple_Gauge_maxiter;
  const double conv_tol_gauge =
      peps_parameters.Simple_Gauge_Convergence_Epsilon;
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
          core::fix_local_gauge(Tn[source], Tn[target], lambda_tensor[source],
                                lambda_tensor[target], source_leg,
                                peps_parameters, Tn1_work, Tn2_work,
                                lambda_work);
          lambda_tensor[source][source_leg] = lambda_work;
          lambda_tensor[target][target_leg] = lambda_work;
          Tn[source] = Tn1_work;
          Tn[target] = Tn2_work;
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
    if (score < conv_tol_gauge) {
      break;
    }
  }  // end of for (iter_gauge)
}

template <class tensor>
void iTPS<tensor>::simple_update() {
  const int group = 0;
  Timer<> timer;
  tensor Tn1_new(comm);
  tensor Tn2_new(comm);
  std::vector<double> lambda_c;
  const int nsteps = peps_parameters.num_simple_step[group];
  double next_report = 10.0;

  log_fermion_sector_dimensions(finfo, lattice, N_UNIT, mpirank, 0);
  for (int int_tau = 0; int_tau < nsteps; ++int_tau) {
    for (auto up : simple_updates) {
      if (up.group != group) {
        continue;
      }
      simple_update(up);
    }

    // local gauge fixing
    if (peps_parameters.Simple_Gauge_Fix) {
      fix_local_gauge();
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
    const int step = int_tau + 1;
    if (step == 10 || step == 50 || step == 100 || step == 300) {
      log_fermion_sector_dimensions(finfo, lattice, N_UNIT, mpirank, step);
    }
  }  // end of for (int_tau)
  time_simple_update += timer.elapsed();
}

template <class tensor>
void iTPS<tensor>::simple_update_density(EvolutionOperator<tensor> const& up) {
  if (up.is_onesite()) {
    const int source = up.source_site;
    Tn[source] =
        tensordot(Tn[source], up.op, mptensor::Axes(4), mptensor::Axes(0));
    Tn[source] = tensordot(Tn[source], conj(up.op), mptensor::Axes(4),
                           mptensor::Axes(1));
  } else {
    tensor Tn1_work(comm), Tn2_work(comm);
    std::vector<double> lambda_work;
    const int source = up.source_site;
    const int source_leg = up.source_leg;
    const int target = lattice.neighbor(source, source_leg);
    const int target_leg = (source_leg + 2) % 4;
    core::Simple_update_bond(
        reshape(Tn[source],
                mptensor::Shape(Tn[source].shape()[0], Tn[source].shape()[1],
                                Tn[source].shape()[2], Tn[source].shape()[3],
                                Tn[source].shape()[4] * Tn[source].shape()[5])),
        reshape(Tn[target],
                mptensor::Shape(Tn[target].shape()[0], Tn[target].shape()[1],
                                Tn[target].shape()[2], Tn[target].shape()[3],
                                Tn[target].shape()[4] * Tn[target].shape()[5])),
        lambda_tensor[source], lambda_tensor[target],
        mptensor::kron(up.op, conj(up.op)), source_leg, peps_parameters,
        Tn1_work, Tn2_work, lambda_work);
    lambda_tensor[source][source_leg] = lambda_work;
    lambda_tensor[target][target_leg] = lambda_work;
    Tn[source] = reshape(
        Tn1_work, mptensor::Shape(Tn[source].shape()[0], Tn[source].shape()[1],
                                  Tn[source].shape()[2], Tn[source].shape()[3],
                                  up.op.shape()[2], up.op.shape()[2]));
    Tn[target] = reshape(
        Tn2_work, mptensor::Shape(Tn[target].shape()[0], Tn[target].shape()[1],
                                  Tn[target].shape()[2], Tn[target].shape()[3],
                                  up.op.shape()[3], up.op.shape()[3]));
  }
}

template <class tensor>
void iTPS<tensor>::fix_local_gauge_density() {
  tensor Tn1_work(comm), Tn2_work(comm);
  std::vector<double> lambda_work;
  const int maxiter_gauge = peps_parameters.Simple_Gauge_maxiter;
  const double conv_tol_gauge =
      peps_parameters.Simple_Gauge_Convergence_Epsilon;
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
              reshape(Tn[source],
                      mptensor::Shape(
                          Tn[source].shape()[0], Tn[source].shape()[1],
                          Tn[source].shape()[2], Tn[source].shape()[3],
                          Tn[source].shape()[4] * Tn[source].shape()[5])),
              reshape(Tn[target],
                      mptensor::Shape(
                          Tn[target].shape()[0], Tn[target].shape()[1],
                          Tn[target].shape()[2], Tn[target].shape()[3],
                          Tn[target].shape()[4] * Tn[target].shape()[5])),
              lambda_tensor[source], lambda_tensor[target], source_leg,
              peps_parameters, Tn1_work, Tn2_work, lambda_work);
          lambda_tensor[source][source_leg] = lambda_work;
          lambda_tensor[target][target_leg] = lambda_work;
          Tn[source] = reshape(
              Tn1_work,
              mptensor::Shape(Tn[source].shape()[0], Tn[source].shape()[1],
                              Tn[source].shape()[2], Tn[source].shape()[3],
                              Tn[source].shape()[4], Tn[source].shape()[5]));
          Tn[target] = reshape(
              Tn2_work,
              mptensor::Shape(Tn[target].shape()[0], Tn[target].shape()[1],
                              Tn[target].shape()[2], Tn[target].shape()[3],
                              Tn[target].shape()[4], Tn[target].shape()[5]));
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
        auto M = core::boundary_tensor(
            reshape(Tn[site],
                    mptensor::Shape(Tn[site].shape()[0], Tn[site].shape()[1],
                                    Tn[site].shape()[2], Tn[site].shape()[3],
                                    Tn[site].shape()[4] * Tn[site].shape()[5])),
            lambda_tensor[site], leg, peps_parameters);

        tensor U;
        std::vector<double> D;
        eigh(M, mptensor::Axes(0), mptensor::Axes(1), D, U);
        for (auto d : D) {
          score = std::max(score, std::abs(d - 1.0));
        }
      }
    }  // end of for (source)
    if (score < conv_tol_gauge) {
      break;
    }
  }  // end of for (iter_gauge)
}

template <class tensor>
void iTPS<tensor>::simple_update_density() {
  const int group = 0;
  Timer<> timer;
  const int nsteps = peps_parameters.num_simple_step[group];
  double next_report = 10.0;

  for (int int_tau = 0; int_tau < nsteps; ++int_tau) {
    for (auto up : simple_updates) {
      if (up.group != group) {
        continue;
      }
      simple_update_density(up);
    }

    // local gauge fixing
    if (peps_parameters.Simple_Gauge_Fix) {
      fix_local_gauge_density();
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
  time_simple_update += timer.elapsed();
}

template <class tensor>
void iTPS<tensor>::simple_update_density_purification(
    EvolutionOperator<tensor> const& up) {
  if (up.is_onesite()) {
    const int source = up.source_site;
    Tn[source] =
        tensordot(Tn[source], up.op, mptensor::Axes(5), mptensor::Axes(0));
  } else {
    tensor Tn1_work(comm), Tn2_work(comm);
    std::vector<double> lambda_work;
    const int source = up.source_site;
    const int source_leg = up.source_leg;
    const int target = lattice.neighbor(source, source_leg);
    const int target_leg = (source_leg + 2) % 4;

    tensor identity(
        comm, mptensor::Shape(Tn[source].shape()[5], Tn[target].shape()[5],
                              Tn[source].shape()[5], Tn[target].shape()[5]));
    for (size_t j1 = 0; j1 < Tn[source].shape()[5]; ++j1) {
      for (size_t k1 = 0; k1 < Tn[target].shape()[5]; ++k1) {
        for (size_t j2 = 0; j2 < Tn[source].shape()[5]; ++j2) {
          for (size_t k2 = 0; k2 < Tn[target].shape()[5]; ++k2) {
            identity.set_value(mptensor::Index(j1, k1, j2, k2),
                               (j1 == j2 && k1 == k2 ? 1.0 : 0.0));
          }
        }
      }
    }

    core::Simple_update_bond(
        reshape(Tn[source],
                mptensor::Shape(Tn[source].shape()[0], Tn[source].shape()[1],
                                Tn[source].shape()[2], Tn[source].shape()[3],
                                Tn[source].shape()[4] * Tn[source].shape()[5])),
        reshape(Tn[target],
                mptensor::Shape(Tn[target].shape()[0], Tn[target].shape()[1],
                                Tn[target].shape()[2], Tn[target].shape()[3],
                                Tn[target].shape()[4] * Tn[target].shape()[5])),
        lambda_tensor[source], lambda_tensor[target],
        mptensor::kron(up.op, identity), source_leg, peps_parameters, Tn1_work,
        Tn2_work, lambda_work);
    lambda_tensor[source][source_leg] = lambda_work;
    lambda_tensor[target][target_leg] = lambda_work;
    Tn[source] = reshape(
        Tn1_work, mptensor::Shape(Tn[source].shape()[0], Tn[source].shape()[1],
                                  Tn[source].shape()[2], Tn[source].shape()[3],
                                  up.op.shape()[2], Tn[source].shape()[5]));
    Tn[target] = reshape(
        Tn2_work, mptensor::Shape(Tn[target].shape()[0], Tn[target].shape()[1],
                                  Tn[target].shape()[2], Tn[target].shape()[3],
                                  up.op.shape()[3], Tn[target].shape()[5]));
  }
}

template <class tensor>
void iTPS<tensor>::simple_update_density_purification() {
  const int group = 0;
  Timer<> timer;
  const int nsteps = peps_parameters.num_simple_step[group];
  double next_report = 10.0;

  for (int int_tau = 0; int_tau < nsteps; ++int_tau) {
    for (auto up : simple_updates) {
      if (up.group != group) {
        continue;
      }
      simple_update_density_purification(up);
    }

    // local gauge fixing
    if (peps_parameters.Simple_Gauge_Fix) {
      fix_local_gauge_density();
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
  time_simple_update += timer.elapsed();
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace tenes::itps
