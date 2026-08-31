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

#define _USE_MATH_DEFINES
#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "iTPS.hpp"

#include "../tensor.hpp"

#include "core/contract.hpp"

namespace tenes::itps {

namespace {

template <class tensor>
std::vector<typename tensor::value_type> gather_rank2_tensor(tensor const &t) {
  using value_type = typename tensor::value_type;
  const size_t d0 = t.shape()[0];
  const size_t d1 = t.shape()[1];
  std::vector<value_type> buf(d0 * d1, 0.0);
  for (size_t local = 0; local < t.local_size(); ++local) {
    const auto index = t.global_index(local);
    value_type v;
    if (t.get_value(index, v)) {
      buf[index[0] * d1 + index[1]] = v;
    }
  }
  allreduce_sum(buf, t.get_comm());
  return buf;
}

}  // namespace

template <class tensor>
auto iTPS<tensor>::measure_onesite()
    -> std::vector<std::vector<typename iTPS<tensor>::tensor_type>> {
  Timer<> timer;
  const bool is_meanfield = peps_parameters.MeanField_Env;
  const bool is_density = peps_parameters.calcmode ==
                          PEPS_Parameters::CalculationMode::finite_temperature;

  const int nlops = num_onesite_operators;
  std::vector<std::vector<tensor_type>> local_obs(
      nlops, std::vector<tensor_type>(
                 N_UNIT, std::numeric_limits<double>::quiet_NaN()));
  std::vector<tensor_type> norm(N_UNIT);

  if (is_meanfield) {
    if (is_density) {
      throw std::runtime_error(
          "Mean field calculation is not implemented for finite temperature.");
    }
    std::vector<tensor> Tn_(Tn);
    for (int i = 0; i < N_UNIT; ++i) {
      for (int leg = 0; leg < nleg; ++leg) {
        const std::vector<double> mf = lambda_tensor[i][leg];
        Tn_[i].multiply_vector(mf, leg);
      }
    }

    for (int i = 0; i < N_UNIT; ++i) {
      norm[i] = core::Contract_one_site_iTPS_MF(Tn_[i], op_identity[i]);
    }

    for (auto const &op : onesite_operators) {
      const int i = op.source_site;
      const auto val = core::Contract_one_site_iTPS_MF(Tn_[i], op.op);
      local_obs[op.group][i] = op.coeff * val / norm[i];
    }
  } else {  // CTM
    if (is_density) {
      for (int i = 0; i < N_UNIT; ++i) {
        norm[i] = core::Contract_one_site_density_CTM(
            C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i], eTb[i], eTl[i], Tn[i],
            op_identity[i]);
      }
      for (auto const &op : onesite_operators) {
        const int i = op.source_site;
        const auto val = core::Contract_one_site_density_CTM(
            C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i], eTb[i], eTl[i], Tn[i],
            op.op);
        local_obs[op.group][i] = op.coeff * val / norm[i];
      }
    } else {
      for (int i = 0; i < N_UNIT; ++i) {
        norm[i] = core::Contract_one_site_iTPS_CTM(
            C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i], eTb[i], eTl[i], Tn[i],
            op_identity[i]);
      }
      for (auto const &op : onesite_operators) {
        const int i = op.source_site;
        const auto val = core::Contract_one_site_iTPS_CTM(
            C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i], eTb[i], eTl[i], Tn[i],
            op.op);
        local_obs[op.group][i] = op.coeff * val / norm[i];
      }
    }
  }
  double norm_real_min = 1e100;
  double norm_imag_abs_max = 0.0;
  for (int i = 0; i < N_UNIT; ++i) {
    norm_real_min = std::min(std::real(norm[i]), norm_real_min);
    norm_imag_abs_max =
        std::max(std::abs(std::imag(norm[i])), norm_imag_abs_max);
  }
  if (mpirank == 0) {
    if (norm_real_min <= 0.0) {
      std::cerr << "WARNING: Norm is negative [min(real(NORM)) = "
                << norm_real_min << "].\n";
      std::cerr << "HINT: Increase the bond dimension of CTM." << std::endl;
    }
    if (norm_imag_abs_max > 1e-6) {
      std::cerr << "WARNING: Norm is not real [max(abs(imag(NORM))) = "
                << norm_imag_abs_max << " > 1e-6].\n";
      std::cerr << "HINT: Increase the bond dimension of CTM." << std::endl;
    }
  }
  local_obs.push_back(norm);

  time_observable += timer.elapsed();
  return local_obs;
}

template <class tensor>
auto iTPS<tensor>::measure_onesite_rdm()
    -> std::vector<small_tensor<typename iTPS<tensor>::tensor_type>> {
  const bool is_meanfield = peps_parameters.MeanField_Env;
  const bool is_density = peps_parameters.calcmode ==
                          PEPS_Parameters::CalculationMode::finite_temperature;

  std::vector<small_tensor<tensor_type>> rdm_all;
  rdm_all.reserve(N_UNIT);

  std::vector<tensor> Tn_mf;
  if (is_meanfield) {
    if (is_density) {
      throw std::runtime_error(
          "Mean field calculation is not implemented for finite temperature.");
    }
    Tn_mf = Tn;
    for (int i = 0; i < N_UNIT; ++i) {
      for (int leg = 0; leg < nleg; ++leg) {
        const std::vector<double> mf = lambda_tensor[i][leg];
        Tn_mf[i].multiply_vector(mf, leg);
      }
    }
  }

  for (int i = 0; i < N_UNIT; ++i) {
    tensor rdm = is_meanfield ? core::Contract_one_site_RDM_iTPS_MF(Tn_mf[i])
                 : is_density ? core::Contract_one_site_RDM_density_CTM(
                                    C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i],
                                    eTb[i], eTl[i], Tn[i])
                              : core::Contract_one_site_RDM_iTPS_CTM(
                                    C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i],
                                    eTb[i], eTl[i], Tn[i]);
    const size_t d0 = rdm.shape()[0];
    const size_t d1 = rdm.shape()[1];
    assert(d0 == d1);
    const auto buf = gather_rank2_tensor(rdm);
    small_tensor<tensor_type> rdm_local{mptensor::Shape(d0, d1)};
    tensor_type tr = 0.0;

    for (size_t a = 0; a < d0; ++a) {
      tr += buf[a * d1 + a];
    }

    for (size_t row = 0; row < d0; ++row) {
      for (size_t col = 0; col < d1; ++col) {
        const tensor_type v = (is_meanfield || is_density)
                                  ? buf[col * d1 + row]
                                  : buf[row * d1 + col];
        rdm_local.set_value({row, col}, v / tr);
      }
    }
    rdm_all.push_back(rdm_local);
  }

  return rdm_all;
}

template <class ptensor>
void iTPS<ptensor>::save_onesite(
    std::vector<std::vector<typename iTPS<ptensor>::tensor_type>> const
        &onesite_obs,
    std::optional<double> time, std::string filename_prefix) {
  if (mpirank != 0) {
    return;
  }

  const int nlops = num_onesite_operators;
  std::string filepath = outdir + "/" + filename_prefix + "onesite_obs.dat";
  if (!time && peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "    Save onesite observables to " << filepath << std::endl;
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
    ofs << "# $" << index++ << ": op_group\n";
    ofs << "# $" << index++ << ": site_index\n";
    ofs << "# $" << index++ << ": real\n";
    ofs << "# $" << index++ << ": imag\n";

    ofs << "# The names of op_group are the following: \n";
    for (int ilops = 0; ilops < num_onesite_operators; ++ilops) {
      ofs << "# " << ilops << ": " << onesite_operator_names[ilops] << "\n";
    }
    if (onesite_obs.size() == static_cast<std::size_t>(nlops) + 1) {
      ofs << "# -1: norm\n";
    }
    ofs << std::endl;
  }

  std::ofstream ofs(filepath.c_str(), std::ios::out | std::ios::app);
  ofs << std::scientific
      << std::setprecision(std::numeric_limits<double>::max_digits10);

  for (int ilops = 0; ilops < nlops; ++ilops) {
    for (int i = 0; i < N_UNIT; ++i) {
      const auto v = onesite_obs[ilops][i];
      if (std::isnan(std::real(v))) {
        continue;
      }
      if (time) {
        ofs << (*time) << " ";
      }
      ofs << ilops << " " << i << " " << std::real(v) << " " << std::imag(v)
          << std::endl;
    }
  }
  if (onesite_obs.size() == static_cast<std::size_t>(nlops) + 1) {
    // includes norm
    for (int i = 0; i < N_UNIT; ++i) {
      const auto v = onesite_obs[nlops][i];
      if (std::isnan(std::real(v))) {
        continue;
      }
      if (time) {
        ofs << (*time) << " ";
      }
      ofs << "-1 " << i << " " << std::real(v) << " " << std::imag(v)
          << std::endl;
    }
  }
}

template <class ptensor>
void iTPS<ptensor>::save_onesite_rdm(
    std::vector<small_tensor<typename iTPS<ptensor>::tensor_type>> const
        &onesite_rdm,
    std::optional<double> time, std::string filename_prefix) {
  if (mpirank != 0) {
    return;
  }

  std::string filepath =
      outdir + "/" + filename_prefix + "onesite_density_matrix.dat";
  if (!time && peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "    Save one-site RDMs to " << filepath << std::endl;
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
    ofs << "# $" << index++ << ": site\n";
    ofs << "# $" << index++ << ": row index (bra)\n";
    ofs << "# $" << index++ << ": col index (ket)\n";
    ofs << "# $" << index++ << ": real part\n";
    ofs << "# $" << index++ << ": imag part\n";
    ofs << std::endl;
  }

  std::ofstream ofs(filepath.c_str(), std::ios::out | std::ios::app);
  ofs << std::scientific
      << std::setprecision(std::numeric_limits<double>::max_digits10);

  for (int i = 0; i < N_UNIT; ++i) {
    const size_t d0 = onesite_rdm[i].shape()[0];
    const size_t d1 = onesite_rdm[i].shape()[1];
    for (size_t row = 0; row < d0; ++row) {
      for (size_t col = 0; col < d1; ++col) {
        tensor_type v;
        onesite_rdm[i].get_value({row, col}, v);
        if (time) {
          ofs << (*time) << " ";
        }
        ofs << i << " " << row << " " << col << " " << std::real(v) << " "
            << std::imag(v) << std::endl;
      }
    }
  }
}

template <class tensor>
auto iTPS<tensor>::measure_onesite_density()
    -> std::vector<std::vector<typename iTPS<tensor>::tensor_type>> {
  Timer<> timer;
  const int nlops = num_onesite_operators;
  std::vector<std::vector<tensor_type>> local_obs(
      nlops, std::vector<tensor_type>(
                 N_UNIT, std::numeric_limits<double>::quiet_NaN()));
  std::vector<tensor_type> norm(N_UNIT);

  /*
  if (peps_parameters.MeanField_Env) {
    std::vector<tensor> Tn_(Tn);
    for (int i = 0; i < N_UNIT; ++i) {
      for (int leg = 0; leg < nleg; ++leg) {
        const std::vector<double> mf = lambda_tensor[i][leg];
        Tn_[i].multiply_vector(mf, leg);
      }
    }

    for (int i = 0; i < N_UNIT; ++i) {
      norm[i] = core::Contract_one_site_MF_density(Tn_[i], op_identity[i]);
    }

    for (auto const &op : onesite_operators) {
      const int i = op.source_site;
      const auto val = core::Contract_one_site_MF_density(Tn_[i], op.op);
      local_obs[op.group][i] = val / norm[i];
    }
    } else {*/
  for (int i = 0; i < N_UNIT; ++i) {
    norm[i] = core::Contract_one_site_density_CTM(
        C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i], eTb[i], eTl[i], Tn[i],
        op_identity[i]);
  }
  for (auto const &op : onesite_operators) {
    const int i = op.source_site;
    const auto val = core::Contract_one_site_density_CTM(
        C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i], eTb[i], eTl[i], Tn[i],
        op.op);
    local_obs[op.group][i] = op.coeff * val / norm[i];
  }
  //  }
  double norm_real_min = 1e100;
  double norm_imag_abs_max = 0.0;
  for (int i = 0; i < N_UNIT; ++i) {
    norm_real_min = std::min(std::real(norm[i]), norm_real_min);
    norm_imag_abs_max =
        std::max(std::abs(std::imag(norm[i])), norm_imag_abs_max);
  }
  if (mpirank == 0) {
    if (norm_real_min <= 0.0) {
      std::cerr << "WARNING: Norm is negative [min(real(NORM)) = "
                << norm_real_min << "].\n";
      std::cerr << "HINT: Increase the bond dimension of CTM." << std::endl;
    }
    if (norm_imag_abs_max > 1e-6) {
      std::cerr << "WARNING: Norm is not real [max(abs(imag(NORM))) = "
                << norm_imag_abs_max << " > 1e-6].\n";
      std::cerr << "HINT: Increase the bond dimension of CTM." << std::endl;
    }
  }
  local_obs.push_back(norm);

  time_observable += timer.elapsed();
  return local_obs;
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace tenes::itps
