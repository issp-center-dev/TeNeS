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

#include "correlation_function.hpp"

#include <algorithm>
#include <iomanip>
#include <string>
#include <type_traits>

#include "../fermion/relay.hpp"
#include "iTPS.hpp"

#include "core/contract.hpp"

namespace tenes::itps {

using mptensor::Shape;

namespace {

template <class ptensor>
ptensor fermion_correlation_site(
    const tenes::fermion::ftensor<ptensor> &Tn, tenes::fermion::relay_role role,
    int entry, int exit, const tenes::fermion::relay_channel<ptensor> &channel,
    bool vertical) {
  ptensor site =
      tenes::fermion::build_relay_site(Tn, role, entry, exit, channel);
  if (vertical) {
    site = transpose(site, mptensor::Axes(3, 0, 1, 2, 4, 5));
  }
  return site;
}

}  // namespace

template <class ptensor>
std::vector<Correlation> iTPS<ptensor>::measure_correlation() {
  validate_fermion_ctm_measurement();

  if (peps_parameters.MeanField_Env) {
    return measure_correlation_mf();
  } else {
    return measure_correlation_ctm();
  }
}

template <class ptensor>
std::vector<Correlation> iTPS<ptensor>::measure_correlation_fermion(
    bool meanfield_env) {
  Timer<> timer;
  ScopedTimer scoped_timer("measure/correlation");

  const int nlops = num_onesite_operators;
  const int r_max = corparam.r_max;
  std::vector<std::vector<int>> r_ops(nlops);
  for (auto [left_op, right_op] : corparam.operators) {
    r_ops[left_op].push_back(right_op);
  }

  std::vector<tenes::fermion::ftensor<ptensor>> fTn;
  std::vector<ptensor> reduced;
  std::vector<ptensor> reduced_vertical;
  std::vector<ptensor> Tn_horizontal;
  std::vector<ptensor> Tn_vertical;
  if (meanfield_env) {
    Tn_horizontal.assign(Tn.begin(), Tn.end());
    Tn_vertical.assign(Tn.begin(), Tn.end());
    for (int index = 0; index < N_UNIT; ++index) {
      std::vector<std::vector<double>> const &lambda = lambda_tensor[index];
      Tn_horizontal[index].multiply_vector(lambda[1], 1, lambda[3], 3);
      Tn_vertical[index].multiply_vector(lambda[0], 0, lambda[2], 2);
    }
  } else {
    fTn.reserve(N_UNIT);
    reduced.reserve(N_UNIT);
    reduced_vertical.reserve(N_UNIT);
    for (int index = 0; index < N_UNIT; ++index) {
      fTn.push_back(tenes::fermion::wrap_Tn(Tn[index], finfo, index));
      reduced.push_back(tenes::fermion::build_reduced_op(fTn.back()));
      reduced_vertical.push_back(
          transpose(reduced.back(), mptensor::Axes(3, 0, 1, 2, 4, 5)));
    }
  }

  const auto delta_corner = tenes::fermion::make_delta_corner<ptensor>(comm);
  const auto delta_edge = [&](int D) {
    return tenes::fermion::make_delta_edge<ptensor>(D, comm);
  };

  std::vector<std::vector<ptensor>> relay_middle_horizontal(
      2, std::vector<ptensor>(N_UNIT));
  std::vector<std::vector<ptensor>> relay_middle_vertical(
      2, std::vector<ptensor>(N_UNIT));
  std::vector<std::vector<bool>> relay_middle_horizontal_ready(
      2, std::vector<bool>(N_UNIT, false));
  std::vector<std::vector<bool>> relay_middle_vertical_ready(
      2, std::vector<bool>(N_UNIT, false));
  std::vector<ptensor> reduced_horizontal(N_UNIT);
  std::vector<ptensor> reduced_vertical_mf(N_UNIT);
  std::vector<bool> reduced_horizontal_ready(N_UNIT, false);
  std::vector<bool> reduced_vertical_mf_ready(N_UNIT, false);

  const auto horizontal_reduced = [&](int index) -> const ptensor & {
    if (meanfield_env) {
      if (!reduced_horizontal_ready[index]) {
        const auto fT =
            tenes::fermion::wrap_Tn(Tn_horizontal[index], finfo, index);
        reduced_horizontal[index] = tenes::fermion::build_reduced_op(fT);
        reduced_horizontal_ready[index] = true;
      }
      return reduced_horizontal[index];
    }
    return reduced[index];
  };
  const auto vertical_reduced = [&](int index) -> const ptensor & {
    if (meanfield_env) {
      if (!reduced_vertical_mf_ready[index]) {
        const auto fT =
            tenes::fermion::wrap_Tn(Tn_vertical[index], finfo, index);
        reduced_vertical_mf[index] =
            transpose(tenes::fermion::build_reduced_op(fT),
                      mptensor::Axes(3, 0, 1, 2, 4, 5));
        reduced_vertical_mf_ready[index] = true;
      }
      return reduced_vertical_mf[index];
    }
    return reduced_vertical[index];
  };
  const auto horizontal_middle =
      [&](int index, const tenes::fermion::relay_channel<ptensor> &channel)
      -> const ptensor & {
    const int parity = channel.u.parity[2][0] ? 1 : 0;
    if (!relay_middle_horizontal_ready[parity][index]) {
      const auto fT = meanfield_env ? tenes::fermion::wrap_Tn(
                                          Tn_horizontal[index], finfo, index)
                                    : fTn[index];
      relay_middle_horizontal[parity][index] = tenes::fermion::build_relay_site(
          fT, tenes::fermion::relay_role::middle, 0, 2, channel);
      relay_middle_horizontal_ready[parity][index] = true;
    }
    return relay_middle_horizontal[parity][index];
  };
  const auto vertical_middle =
      [&](int index, const tenes::fermion::relay_channel<ptensor> &channel)
      -> const ptensor & {
    const int parity = channel.u.parity[2][0] ? 1 : 0;
    if (!relay_middle_vertical_ready[parity][index]) {
      const auto fT = meanfield_env ? tenes::fermion::wrap_Tn(
                                          Tn_vertical[index], finfo, index)
                                    : fTn[index];
      relay_middle_vertical[parity][index] = fermion_correlation_site(
          fT, tenes::fermion::relay_role::middle, 3, 1, channel, true);
      relay_middle_vertical_ready[parity][index] = true;
    }
    return relay_middle_vertical[parity][index];
  };

  std::vector<Correlation> correlations;
  for (int left_index = 0; left_index < N_UNIT; ++left_index) {
    const auto vdim = lattice.virtual_dims[left_index];
    ptensor correlation_T =
        meanfield_env ? ptensor{} : ptensor(comm, Shape(CHI, CHI, vdim[0]));
    ptensor correlation_norm =
        meanfield_env ? ptensor{} : ptensor(comm, Shape(CHI, CHI, vdim[0]));
    for (int left_ilop = 0; left_ilop < nlops; ++left_ilop) {
      if (r_ops[left_ilop].empty()) {
        continue;
      }

      {  // horizontal
        int left_op_index = siteoperator_index(left_index, left_ilop);
        if (left_op_index < 0) {
          continue;
        }
        tenes::fermion::ftensor<ptensor> fSource;
        if (meanfield_env) {
          ptensor T = Tn_horizontal[left_index];
          T.multiply_vector(lambda_tensor[left_index][0], 0);
          fSource = tenes::fermion::wrap_Tn(T, finfo, left_index);
        } else {
          fSource = fTn[left_index];
        }
        const auto left_op = onesite_operators[left_op_index].op;
        tenes::fermion::relay_channel<ptensor> channel;
        const bool left_odd =
            onesite_parity[left_op_index] == tenes::fermion::op_parity::odd;
        channel.u = tenes::fermion::relay_product_source(
            left_op, finfo.phys[left_index], left_odd);
        channel.vt = tenes::fermion::relay_product_target(
            left_op, finfo.phys[left_index], left_odd);
        auto source = tenes::fermion::build_relay_site(
            fSource, tenes::fermion::relay_role::source, -1, 2, channel);
        if (meanfield_env) {
          auto source_reduced = tenes::fermion::build_reduced_op(fSource);
          auto eTt_source = delta_edge(static_cast<int>(fSource.shape()[1]));
          auto eTb_source = delta_edge(static_cast<int>(fSource.shape()[3]));
          auto eTl_source = delta_edge(static_cast<int>(fSource.shape()[0]));
          core::StartCorrelation_density_CTM(
              correlation_T, delta_corner, delta_corner, eTt_source, eTb_source,
              eTl_source, source, op_identity[left_index]);
          core::StartCorrelation_density_CTM(
              correlation_norm, delta_corner, delta_corner, eTt_source,
              eTb_source, eTl_source, source_reduced, op_identity[left_index]);
        } else {
          core::StartCorrelation_density_CTM(correlation_T, C1[left_index],
                                             C4[left_index], eTt[left_index],
                                             eTb[left_index], eTl[left_index],
                                             source, op_identity[left_index]);
          core::StartCorrelation_density_CTM(
              correlation_norm, C1[left_index], C4[left_index], eTt[left_index],
              eTb[left_index], eTl[left_index], reduced[left_index],
              op_identity[left_index]);
        }

        int right_index = left_index;
        for (int r = 0; r < r_max; ++r) {
          right_index = lattice.right(right_index);
          tenes::fermion::ftensor<ptensor> fTarget;
          tensor_type norm = 0.0;
          ptensor eTt_target, eTr_target, eTb_target;
          if (meanfield_env) {
            ptensor Ttarget = Tn_horizontal[right_index];
            Ttarget.multiply_vector(lambda_tensor[right_index][2], 2);
            fTarget = tenes::fermion::wrap_Tn(Ttarget, finfo, right_index);
            auto target_reduced = tenes::fermion::build_reduced_op(fTarget);
            eTt_target = delta_edge(static_cast<int>(fTarget.shape()[1]));
            eTr_target = delta_edge(static_cast<int>(fTarget.shape()[2]));
            eTb_target = delta_edge(static_cast<int>(fTarget.shape()[3]));
            norm = core::FinishCorrelation_density_CTM(
                correlation_norm, delta_corner, delta_corner, eTt_target,
                eTr_target, eTb_target, target_reduced,
                op_identity[right_index]);
          } else {
            fTarget = fTn[right_index];
            norm = core::FinishCorrelation_density_CTM(
                correlation_norm, C2[right_index], C3[right_index],
                eTt[right_index], eTr[right_index], eTb[right_index],
                reduced[right_index], op_identity[right_index]);
          }
          for (auto right_ilop : r_ops[left_ilop]) {
            int right_op_index = siteoperator_index(right_index, right_ilop);
            if (right_op_index < 0) {
              continue;
            }
            const auto right_op = onesite_operators[right_op_index].op;
            tensor_type val = 0.0;
            const auto pA = onesite_parity[left_op_index];
            const auto pB = onesite_parity[right_op_index];
            if (pA == pB) {
              channel.vt = tenes::fermion::relay_product_target(
                  right_op, finfo.phys[right_index], left_odd);
              auto target = tenes::fermion::build_relay_site(
                  fTarget, tenes::fermion::relay_role::target, 0, -1, channel);
              if (meanfield_env) {
                val = core::FinishCorrelation_density_CTM(
                          correlation_T, delta_corner, delta_corner, eTt_target,
                          eTr_target, eTb_target, target,
                          op_identity[right_index]) /
                      norm;
              } else {
                val = core::FinishCorrelation_density_CTM(
                          correlation_T, C2[right_index], C3[right_index],
                          eTt[right_index], eTr[right_index], eTb[right_index],
                          target, op_identity[right_index]) /
                      norm;
              }
            }
            correlations.push_back(Correlation{left_index, r + 1, 0, left_ilop,
                                               right_ilop, std::real(val),
                                               std::imag(val)});
          }

          if (!meanfield_env) {
            core::Transfer_density_CTM(correlation_T, eTt[right_index],
                                       eTb[right_index],
                                       horizontal_middle(right_index, channel));
            core::Transfer_density_CTM(correlation_norm, eTt[right_index],
                                       eTb[right_index],
                                       horizontal_reduced(right_index));
          } else {
            auto eTt_middle = delta_edge(
                static_cast<int>(Tn_horizontal[right_index].shape()[1]));
            auto eTb_middle = delta_edge(
                static_cast<int>(Tn_horizontal[right_index].shape()[3]));
            core::Transfer_density_CTM(correlation_T, eTt_middle, eTb_middle,
                                       horizontal_middle(right_index, channel));
            core::Transfer_density_CTM(correlation_norm, eTt_middle, eTb_middle,
                                       horizontal_reduced(right_index));
          }
        }
      }
      {  // vertical
        int left_op_index = siteoperator_index(left_index, left_ilop);
        if (left_op_index < 0) {
          continue;
        }
        const auto left_op = onesite_operators[left_op_index].op;
        tenes::fermion::ftensor<ptensor> fSource;
        if (meanfield_env) {
          ptensor T = Tn_vertical[left_index];
          T.multiply_vector(lambda_tensor[left_index][3], 3);
          fSource = tenes::fermion::wrap_Tn(T, finfo, left_index);
        } else {
          fSource = fTn[left_index];
        }
        tenes::fermion::relay_channel<ptensor> channel;
        const bool left_odd =
            onesite_parity[left_op_index] == tenes::fermion::op_parity::odd;
        channel.u = tenes::fermion::relay_product_source(
            left_op, finfo.phys[left_index], left_odd);
        channel.vt = tenes::fermion::relay_product_target(
            left_op, finfo.phys[left_index], left_odd);
        auto source = fermion_correlation_site(
            fSource, tenes::fermion::relay_role::source, -1, 1, channel, true);
        if (meanfield_env) {
          auto source_reduced =
              transpose(tenes::fermion::build_reduced_op(fSource),
                        mptensor::Axes(3, 0, 1, 2, 4, 5));
          auto eTl_source = delta_edge(static_cast<int>(fSource.shape()[0]));
          auto eTr_source = delta_edge(static_cast<int>(fSource.shape()[2]));
          auto eTb_source = delta_edge(static_cast<int>(fSource.shape()[3]));
          core::StartCorrelation_density_CTM(
              correlation_T, delta_corner, delta_corner, eTl_source, eTr_source,
              eTb_source, source, op_identity[left_index]);
          core::StartCorrelation_density_CTM(
              correlation_norm, delta_corner, delta_corner, eTl_source,
              eTr_source, eTb_source, source_reduced, op_identity[left_index]);
        } else {
          core::StartCorrelation_density_CTM(correlation_T, C4[left_index],
                                             C3[left_index], eTl[left_index],
                                             eTr[left_index], eTb[left_index],
                                             source, op_identity[left_index]);
          core::StartCorrelation_density_CTM(
              correlation_norm, C4[left_index], C3[left_index], eTl[left_index],
              eTr[left_index], eTb[left_index], reduced_vertical[left_index],
              op_identity[left_index]);
        }

        int right_index = left_index;
        for (int r = 0; r < r_max; ++r) {
          right_index = lattice.top(right_index);
          tenes::fermion::ftensor<ptensor> fTarget;
          tensor_type norm = 0.0;
          ptensor eTl_target, eTt_target, eTr_target;
          if (meanfield_env) {
            ptensor Ttarget = Tn_vertical[right_index];
            Ttarget.multiply_vector(lambda_tensor[right_index][1], 1);
            fTarget = tenes::fermion::wrap_Tn(Ttarget, finfo, right_index);
            auto target_reduced =
                transpose(tenes::fermion::build_reduced_op(fTarget),
                          mptensor::Axes(3, 0, 1, 2, 4, 5));
            eTl_target = delta_edge(static_cast<int>(fTarget.shape()[0]));
            eTt_target = delta_edge(static_cast<int>(fTarget.shape()[1]));
            eTr_target = delta_edge(static_cast<int>(fTarget.shape()[2]));
            norm = core::FinishCorrelation_density_CTM(
                correlation_norm, delta_corner, delta_corner, eTl_target,
                eTt_target, eTr_target, target_reduced,
                op_identity[right_index]);
          } else {
            fTarget = fTn[right_index];
            norm = core::FinishCorrelation_density_CTM(
                correlation_norm, C1[right_index], C2[right_index],
                eTl[right_index], eTt[right_index], eTr[right_index],
                reduced_vertical[right_index], op_identity[right_index]);
          }
          for (auto right_ilop : r_ops[left_ilop]) {
            int right_op_index = siteoperator_index(right_index, right_ilop);
            if (right_op_index < 0) {
              continue;
            }
            const auto right_op = onesite_operators[right_op_index].op;
            tensor_type val = 0.0;
            const auto pA = onesite_parity[left_op_index];
            const auto pB = onesite_parity[right_op_index];
            if (pA == pB) {
              channel.vt = tenes::fermion::relay_product_target(
                  right_op, finfo.phys[right_index], left_odd);
              auto target = fermion_correlation_site(
                  fTarget, tenes::fermion::relay_role::target, 3, -1, channel,
                  true);
              if (meanfield_env) {
                val = core::FinishCorrelation_density_CTM(
                          correlation_T, delta_corner, delta_corner, eTl_target,
                          eTt_target, eTr_target, target,
                          op_identity[right_index]) /
                      norm;
              } else {
                val = core::FinishCorrelation_density_CTM(
                          correlation_T, C1[right_index], C2[right_index],
                          eTl[right_index], eTt[right_index], eTr[right_index],
                          target, op_identity[right_index]) /
                      norm;
              }
            }
            correlations.push_back(Correlation{left_index, 0, r + 1, left_ilop,
                                               right_ilop, std::real(val),
                                               std::imag(val)});
          }

          if (!meanfield_env) {
            core::Transfer_density_CTM(correlation_T, eTl[right_index],
                                       eTr[right_index],
                                       vertical_middle(right_index, channel));
            core::Transfer_density_CTM(correlation_norm, eTl[right_index],
                                       eTr[right_index],
                                       reduced_vertical[right_index]);
          } else {
            auto eTl_middle = delta_edge(
                static_cast<int>(Tn_vertical[right_index].shape()[0]));
            auto eTr_middle = delta_edge(
                static_cast<int>(Tn_vertical[right_index].shape()[2]));
            core::Transfer_density_CTM(correlation_T, eTl_middle, eTr_middle,
                                       vertical_middle(right_index, channel));
            core::Transfer_density_CTM(correlation_norm, eTl_middle, eTr_middle,
                                       vertical_reduced(right_index));
          }
        }
      }
    }
  }

  time_observable += timer.elapsed();
  return correlations;
}

template <class ptensor>
std::vector<Correlation> iTPS<ptensor>::measure_correlation_ctm() {
  const bool is_tpo = peps_parameters.calcmode ==
                      PEPS_Parameters::CalculationMode::finite_temperature;
  if (finfo.enabled && !is_tpo) {
    return measure_correlation_fermion(false);
  }

  Timer<> timer;
  ScopedTimer scoped_timer("measure/correlation");

  const int nlops = num_onesite_operators;
  const int r_max = corparam.r_max;
  std::vector<std::vector<int>> r_ops(nlops);
  for (auto [left_op, right_op] : corparam.operators) {
    r_ops[left_op].push_back(right_op);
  }

  std::vector<Correlation> correlations;
  for (int left_index = 0; left_index < N_UNIT; ++left_index) {
    const auto vdim = lattice.virtual_dims[left_index];
    ptensor correlation_T, correlation_norm;
    if (is_tpo) {
      correlation_T = ptensor(comm, Shape(CHI, CHI, vdim[0]));
      correlation_norm = ptensor(comm, Shape(CHI, CHI, vdim[0]));
    } else {
      correlation_T = ptensor(comm, Shape(CHI, CHI, vdim[0], vdim[0]));
      correlation_norm = ptensor(comm, Shape(CHI, CHI, vdim[0], vdim[0]));
    }
    for (int left_ilop = 0; left_ilop < nlops; ++left_ilop) {
      if (r_ops[left_ilop].empty()) {
        continue;
      }

      {  // horizontal
        int left_op_index = siteoperator_index(left_index, left_ilop);
        if (left_op_index < 0) {
          continue;
        }
        const auto left_op = onesite_operators[left_op_index].op;
        if (is_tpo) {
          core::StartCorrelation_density_CTM(
              correlation_T, C1[left_index], C4[left_index], eTt[left_index],
              eTb[left_index], eTl[left_index], Tn[left_index], left_op);
          core::StartCorrelation_density_CTM(
              correlation_norm, C1[left_index], C4[left_index], eTt[left_index],
              eTb[left_index], eTl[left_index], Tn[left_index],
              op_identity[left_index]);
        } else {
          core::StartCorrelation_iTPS_CTM(
              correlation_T, C1[left_index], C4[left_index], eTt[left_index],
              eTb[left_index], eTl[left_index], Tn[left_index], left_op);
          core::StartCorrelation_iTPS_CTM(
              correlation_norm, C1[left_index], C4[left_index], eTt[left_index],
              eTb[left_index], eTl[left_index], Tn[left_index],
              op_identity[left_index]);
        }

        int right_index = left_index;
        for (int r = 0; r < r_max; ++r) {
          right_index = lattice.right(right_index);
          tensor_type norm =
              is_tpo ? core::FinishCorrelation_density_CTM(
                           correlation_norm, C2[right_index], C3[right_index],
                           eTt[right_index], eTr[right_index], eTb[right_index],
                           Tn[right_index], op_identity[right_index])
                     : core::FinishCorrelation_iTPS_CTM(
                           correlation_norm, C2[right_index], C3[right_index],
                           eTt[right_index], eTr[right_index], eTb[right_index],
                           Tn[right_index], op_identity[right_index]);
          for (auto right_ilop : r_ops[left_ilop]) {
            int right_op_index = siteoperator_index(right_index, right_ilop);
            if (right_op_index < 0) {
              continue;
            }
            const auto right_op = onesite_operators[right_op_index].op;
            auto val =
                is_tpo ? core::FinishCorrelation_density_CTM(
                             correlation_T, C2[right_index], C3[right_index],
                             eTt[right_index], eTr[right_index],
                             eTb[right_index], Tn[right_index], right_op) /
                             norm
                       : core::FinishCorrelation_iTPS_CTM(
                             correlation_T, C2[right_index], C3[right_index],
                             eTt[right_index], eTr[right_index],
                             eTb[right_index], Tn[right_index], right_op) /
                             norm;
            correlations.push_back(Correlation{left_index, r + 1, 0, left_ilop,
                                               right_ilop, std::real(val),
                                               std::imag(val)});
          }

          if (is_tpo) {
            core::Transfer_density_CTM(correlation_T, eTt[right_index],
                                       eTb[right_index], Tn[right_index]);
            core::Transfer_density_CTM(correlation_norm, eTt[right_index],
                                       eTb[right_index], Tn[right_index]);
          } else {
            core::Transfer_iTPS_CTM(correlation_T, eTt[right_index],
                                    eTb[right_index], Tn[right_index]);
            core::Transfer_iTPS_CTM(correlation_norm, eTt[right_index],
                                    eTb[right_index], Tn[right_index]);
          }
        }
      }
      {  // vertical
        int left_op_index = siteoperator_index(left_index, left_ilop);
        if (left_op_index < 0) {
          continue;
        }
        const auto left_op = onesite_operators[left_op_index].op;
        ptensor tn =
            is_tpo ? transpose(Tn[left_index], mptensor::Axes(3, 0, 1, 2, 4, 5))
                   : transpose(Tn[left_index], mptensor::Axes(3, 0, 1, 2, 4));
        if (is_tpo) {
          core::StartCorrelation_density_CTM(
              correlation_T, C4[left_index], C3[left_index], eTl[left_index],
              eTr[left_index], eTb[left_index], tn, left_op);
          core::StartCorrelation_density_CTM(
              correlation_norm, C4[left_index], C3[left_index], eTl[left_index],
              eTr[left_index], eTb[left_index], tn, op_identity[left_index]);
        } else {
          core::StartCorrelation_iTPS_CTM(
              correlation_T, C4[left_index], C3[left_index], eTl[left_index],
              eTr[left_index], eTb[left_index], tn, left_op);
          core::StartCorrelation_iTPS_CTM(
              correlation_norm, C4[left_index], C3[left_index], eTl[left_index],
              eTr[left_index], eTb[left_index], tn, op_identity[left_index]);
        }

        int right_index = left_index;
        for (int r = 0; r < r_max; ++r) {
          right_index = lattice.top(right_index);
          tn =
              is_tpo
                  ? transpose(Tn[right_index], mptensor::Axes(3, 0, 1, 2, 4, 5))
                  : transpose(Tn[right_index], mptensor::Axes(3, 0, 1, 2, 4));
          tensor_type norm =
              is_tpo ? core::FinishCorrelation_density_CTM(
                           correlation_norm, C1[right_index], C2[right_index],
                           eTl[right_index], eTt[right_index], eTr[right_index],
                           tn, op_identity[right_index])
                     : core::FinishCorrelation_iTPS_CTM(
                           correlation_norm, C1[right_index], C2[right_index],
                           eTl[right_index], eTt[right_index], eTr[right_index],
                           tn, op_identity[right_index]);
          for (auto right_ilop : r_ops[left_ilop]) {
            int right_op_index = siteoperator_index(right_index, right_ilop);
            if (right_op_index < 0) {
              continue;
            }
            const auto right_op = onesite_operators[right_op_index].op;
            auto val =
                is_tpo ? core::FinishCorrelation_density_CTM(
                             correlation_T, C1[right_index], C2[right_index],
                             eTl[right_index], eTt[right_index],
                             eTr[right_index], tn, right_op) /
                             norm
                       : core::FinishCorrelation_iTPS_CTM(
                             correlation_T, C1[right_index], C2[right_index],
                             eTl[right_index], eTt[right_index],
                             eTr[right_index], tn, right_op) /
                             norm;
            correlations.push_back(Correlation{left_index, 0, r + 1, left_ilop,
                                               right_ilop, std::real(val),
                                               std::imag(val)});
          }

          if (is_tpo) {
            core::Transfer_density_CTM(correlation_T, eTl[right_index],
                                       eTr[right_index], tn);
            core::Transfer_density_CTM(correlation_norm, eTl[right_index],
                                       eTr[right_index], tn);
          } else {
            core::Transfer_iTPS_CTM(correlation_T, eTl[right_index],
                                    eTr[right_index], tn);
            core::Transfer_iTPS_CTM(correlation_norm, eTl[right_index],
                                    eTr[right_index], tn);
          }
        }
      }
    }
  }

  time_observable += timer.elapsed();
  return correlations;
}

template <class ptensor>
std::vector<Correlation> iTPS<ptensor>::measure_correlation_mf() {
  const bool is_tpo = peps_parameters.calcmode ==
                      PEPS_Parameters::CalculationMode::finite_temperature;
  if (is_tpo) {
    throw std::runtime_error(
        "iTPS::measure_correlation_mf() is not implemented for finite "
        "temperature");
  }
  if (finfo.enabled) {
    return measure_correlation_fermion(true);
  }

  ScopedTimer scoped_timer("measure/correlation");
  Timer<> timer;

  const int nlops = num_onesite_operators;
  const int r_max = corparam.r_max;
  std::vector<std::vector<int>> r_ops(nlops);
  for (auto [left_op, right_op] : corparam.operators) {
    r_ops[left_op].push_back(right_op);
  }

  std::vector<ptensor> Tn_horizontal(Tn.begin(), Tn.end());
  std::vector<ptensor> Tn_vertical(Tn.begin(), Tn.end());
  for (int index = 0; index < N_UNIT; ++index) {
    std::vector<std::vector<double>> const &lambda = lambda_tensor[index];
    Tn_horizontal[index].multiply_vector(lambda[1], 1, lambda[3], 3);
    Tn_vertical[index].multiply_vector(lambda[0], 0, lambda[2], 2);
  }

  std::vector<Correlation> correlations;
  for (int left_index = 0; left_index < N_UNIT; ++left_index) {
    const auto vdim = lattice.virtual_dims[left_index];
    ptensor correlation_T(comm, Shape(vdim[0], vdim[0]));
    ptensor correlation_norm(comm, Shape(vdim[0], vdim[0]));
    for (int left_ilop = 0; left_ilop < nlops; ++left_ilop) {
      if (r_ops[left_ilop].empty()) {
        continue;
      }

      {
        const int direction = 2;  // right
        int left_op_index = siteoperator_index(left_index, left_ilop);
        if (left_op_index < 0) {
          continue;
        }
        ptensor T = Tn_horizontal[left_index];
        T.multiply_vector(lambda_tensor[left_index][0], 0);  // 0 means left
        const auto left_op = onesite_operators[left_op_index].op;
        core::StartCorrelation_iTPS_MF(correlation_T, T, left_op, direction);
        core::StartCorrelation_iTPS_MF(correlation_norm, T,
                                       op_identity[left_index], direction);

        int right_index = left_index;
        for (int r = 0; r < r_max; ++r) {
          right_index = lattice.right(right_index);
          T = Tn_horizontal[right_index];
          T.multiply_vector(lambda_tensor[right_index][2], direction);
          tensor_type norm = core::FinishCorrelation_iTPS_MF(
              correlation_norm, T, op_identity[right_index], direction);
          for (auto right_ilop : r_ops[left_ilop]) {
            int right_op_index = siteoperator_index(right_index, right_ilop);
            if (right_op_index < 0) {
              continue;
            }
            const auto right_op = onesite_operators[right_op_index].op;
            auto val = core::FinishCorrelation_iTPS_MF(correlation_T, T,
                                                       right_op, direction) /
                       norm;
            correlations.push_back(Correlation{left_index, r + 1, 0, left_ilop,
                                               right_ilop, std::real(val),
                                               std::imag(val)});
          }

          core::Transfer_iTPS_MF(correlation_T, Tn_horizontal[right_index],
                                 direction);
          core::Transfer_iTPS_MF(correlation_norm, Tn_horizontal[right_index],
                                 direction);
        }
      }
      {                     // vertical
        int direction = 1;  // top
        int left_op_index = siteoperator_index(left_index, left_ilop);
        if (left_op_index < 0) {
          continue;
        }
        ptensor T = Tn_vertical[left_index];
        T.multiply_vector(lambda_tensor[left_index][3], 3);  // 3 means bottom
        const auto left_op = onesite_operators[left_op_index].op;
        core::StartCorrelation_iTPS_MF(correlation_T, T, left_op, direction);
        core::StartCorrelation_iTPS_MF(correlation_norm, T,
                                       op_identity[left_index], direction);

        int right_index = left_index;
        for (int r = 0; r < r_max; ++r) {
          right_index = lattice.top(right_index);
          T = Tn_vertical[right_index];
          T.multiply_vector(lambda_tensor[right_index][1], direction);
          tensor_type norm = core::FinishCorrelation_iTPS_MF(
              correlation_norm, T, op_identity[right_index], direction);
          for (auto right_ilop : r_ops[left_ilop]) {
            int right_op_index = siteoperator_index(right_index, right_ilop);
            if (right_op_index < 0) {
              continue;
            }
            const auto right_op = onesite_operators[right_op_index].op;
            auto val = core::FinishCorrelation_iTPS_MF(correlation_T, T,
                                                       right_op, direction) /
                       norm;
            correlations.push_back(Correlation{left_index, 0, r + 1, left_ilop,
                                               right_ilop, std::real(val),
                                               std::imag(val)});
          }

          core::Transfer_iTPS_MF(correlation_T, Tn_vertical[right_index],
                                 direction);
          core::Transfer_iTPS_MF(correlation_norm, Tn_vertical[right_index],
                                 direction);
        }
      }
    }
  }

  time_observable += timer.elapsed();
  return correlations;
}

template <class ptensor>
void iTPS<ptensor>::save_correlation(
    std::vector<Correlation> const &correlations, std::optional<double> time,
    std::string filename_prefix) {
  if (mpirank != 0) {
    return;
  }
  std::string filepath = outdir + "/" + filename_prefix + "correlation.dat";
  if (!time && peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "    Save long-range correlations to " << filepath
              << std::endl;
  }

  static bool first_time = true;
  if (first_time) {
    first_time = false;
    std::ofstream ofs(filepath.c_str());
    ofs << "# The meaning of each column is the following: \n";
    int index = 1;
    if (time) {
      if (peps_parameters.calcmode ==
          PEPS_Parameters::CalculationMode::time_evolution) {
        ofs << "# $" << index++ << ": time\n";
      } else if (peps_parameters.calcmode ==
                 PEPS_Parameters::CalculationMode::finite_temperature) {
        ofs << "# $" << index++ << ": inverse temperature\n";
      }
    }
    ofs << "# $" << index++ << ": left_op\n";
    ofs << "# $" << index++ << ": left_site\n";
    ofs << "# $" << index++ << ": right_op\n";
    ofs << "# $" << index++ << ": right_dx\n";
    ofs << "# $" << index++ << ": right_dy\n";
    ofs << "# $" << index++ << ": real\n";
    ofs << "# $" << index++ << ": imag\n";

    ofs << "# The names of operators are the following: \n";
    for (int ilops = 0; ilops < num_onesite_operators; ++ilops) {
      ofs << "# " << ilops << ": " << onesite_operator_names[ilops] << "\n";
    }
    ofs << std::endl;
  }
  std::ofstream ofs(filepath.c_str(), std::ios::out | std::ios::app);
  ofs << std::scientific
      << std::setprecision(std::numeric_limits<double>::max_digits10);
  for (auto const &cor : correlations) {
    if (time) {
      ofs << (*time) << " ";
    }
    ofs << cor.left_op << " " << cor.left_index << " " << cor.right_op << " "
        << cor.right_dx << " " << cor.right_dy << " " << cor.real << " "
        << cor.imag << " " << std::endl;
  }
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace tenes::itps
