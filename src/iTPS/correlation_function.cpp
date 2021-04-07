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
#include <string>
#include <type_traits>

#include "iTPS.hpp"

#include "core/contract_ctm.hpp"
#include "core/contract_mf.hpp"

namespace tenes {
namespace itps {

using mptensor::Shape;

template <class ptensor>
std::vector<Correlation> iTPS<ptensor>::measure_correlation() {
  if (peps_parameters.MeanField_Env) {
    return measure_correlation_mf();
  } else {
    return measure_correlation_ctm();
  }
}

template <class ptensor>
std::vector<Correlation> iTPS<ptensor>::measure_correlation_ctm() {
  Timer<> timer;

  const int nlops = num_onesite_operators;
  const int r_max = corparam.r_max;
  std::vector<std::vector<int>> r_ops(nlops);
  for (auto ops : corparam.operators) {
    r_ops[std::get<0>(ops)].push_back(std::get<1>(ops));
  }

  std::vector<Correlation> correlations;
  for (int left_index = 0; left_index < N_UNIT; ++left_index) {
    const auto vdim = lattice.virtual_dims[left_index];
    ptensor correlation_T(Shape(CHI, CHI, vdim[0], vdim[0]));
    ptensor correlation_norm(Shape(CHI, CHI, vdim[0], vdim[0]));
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
        core::StartCorrelation(correlation_T, C1[left_index], C4[left_index],
                               eTt[left_index], eTb[left_index],
                               eTl[left_index], Tn[left_index], left_op);
        core::StartCorrelation(correlation_norm, C1[left_index], C4[left_index],
                               eTt[left_index], eTb[left_index],
                               eTl[left_index], Tn[left_index],
                               op_identity[left_index]);

        int right_index = left_index;
        for (int r = 0; r < r_max; ++r) {
          right_index = lattice.right(right_index);
          double norm = std::real(core::FinishCorrelation(
              correlation_norm, C2[right_index], C3[right_index],
              eTt[right_index], eTr[right_index], eTb[right_index],
              Tn[right_index], op_identity[right_index]));
          for (auto right_ilop : r_ops[left_ilop]) {
            int right_op_index = siteoperator_index(right_index, right_ilop);
            if (right_op_index < 0) {
              continue;
            }
            const auto right_op = onesite_operators[right_op_index].op;
            auto val = core::FinishCorrelation(
                           correlation_T, C2[right_index], C3[right_index],
                           eTt[right_index], eTr[right_index], eTb[right_index],
                           Tn[right_index], right_op) /
                       norm;
            correlations.push_back(Correlation{left_index, r + 1, 0, left_ilop,
                                               right_ilop, std::real(val),
                                               std::imag(val)});
          }

          core::Transfer(correlation_T, eTt[right_index], eTb[right_index],
                         Tn[right_index]);
          core::Transfer(correlation_norm, eTt[right_index], eTb[right_index],
                         Tn[right_index]);
        }
      }
      {  // vertical
        int left_op_index = siteoperator_index(left_index, left_ilop);
        if (left_op_index < 0) {
          continue;
        }
        const auto left_op = onesite_operators[left_op_index].op;
        ptensor tn = transpose(Tn[left_index], mptensor::Axes(3, 0, 1, 2, 4));
        core::StartCorrelation(correlation_T, C4[left_index], C3[left_index],
                               eTl[left_index], eTr[left_index],
                               eTb[left_index], tn, left_op);
        core::StartCorrelation(correlation_norm, C4[left_index], C3[left_index],
                               eTl[left_index], eTr[left_index],
                               eTb[left_index], tn, op_identity[left_index]);

        int right_index = left_index;
        for (int r = 0; r < r_max; ++r) {
          right_index = lattice.top(right_index);
          tn = transpose(Tn[right_index], mptensor::Axes(3, 0, 1, 2, 4));
          double norm = std::real(core::FinishCorrelation(
              correlation_norm, C1[right_index], C2[right_index],
              eTl[right_index], eTt[right_index], eTr[right_index], tn,
              op_identity[right_index]));
          for (auto right_ilop : r_ops[left_ilop]) {
            int right_op_index = siteoperator_index(right_index, right_ilop);
            if (right_op_index < 0) {
              continue;
            }
            const auto right_op = onesite_operators[right_op_index].op;
            auto val = core::FinishCorrelation(
                           correlation_T, C1[right_index], C2[right_index],
                           eTl[right_index], eTt[right_index], eTr[right_index],
                           tn, right_op) /
                       norm;
            correlations.push_back(Correlation{left_index, 0, r + 1, left_ilop,
                                               right_ilop, std::real(val),
                                               std::imag(val)});
          }

          core::Transfer(correlation_T, eTl[right_index], eTr[right_index], tn);
          core::Transfer(correlation_norm, eTl[right_index], eTr[right_index],
                         tn);
        }
      }
    }
  }

  time_observable += timer.elapsed();
  return correlations;
}

template <class ptensor>
std::vector<Correlation> iTPS<ptensor>::measure_correlation_mf() {
  Timer<> timer;

  const int nlops = num_onesite_operators;
  const int r_max = corparam.r_max;
  std::vector<std::vector<int>> r_ops(nlops);
  for (auto ops : corparam.operators) {
    r_ops[std::get<0>(ops)].push_back(std::get<1>(ops));
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
    ptensor correlation_T(Shape(vdim[0], vdim[0]));
    ptensor correlation_norm(Shape(vdim[0], vdim[0]));
    for (int left_ilop = 0; left_ilop < nlops; ++left_ilop) {
      if (r_ops[left_ilop].empty()) {
        continue;
      }

      {  // horizontal
        int left_op_index = siteoperator_index(left_index, left_ilop);
        if (left_op_index < 0) {
          continue;
        }
        ptensor T = Tn_horizontal[left_index];
        T.multiply_vector(lambda_tensor[left_index][0], 0);
        const auto left_op = onesite_operators[left_op_index].op;
        core::StartCorrelation_MF(correlation_T, T, left_op, 2);
        core::StartCorrelation_MF(correlation_norm, T, op_identity[left_index],
                                  2);

        int right_index = left_index;
        for (int r = 0; r < r_max; ++r) {
          right_index = lattice.right(right_index);
          T = Tn_horizontal[right_index];
          T.multiply_vector(lambda_tensor[right_index][2], 2);
          double norm = std::real(core::FinishCorrelation_MF(
              correlation_norm, T, op_identity[right_index], 2));
          for (auto right_ilop : r_ops[left_ilop]) {
            int right_op_index = siteoperator_index(right_index, right_ilop);
            if (right_op_index < 0) {
              continue;
            }
            const auto right_op = onesite_operators[right_op_index].op;
            auto val =
                core::FinishCorrelation_MF(correlation_T, T, right_op, 2) /
                norm;
            correlations.push_back(Correlation{left_index, r + 1, 0, left_ilop,
                                               right_ilop, std::real(val),
                                               std::imag(val)});
          }

          core::Transfer_MF(correlation_T, Tn_horizontal[right_index], 2);
          core::Transfer_MF(correlation_norm, Tn_horizontal[right_index], 2);
        }
      }
      {  // vertical
        int left_op_index = siteoperator_index(left_index, left_ilop);
        if (left_op_index < 0) {
          continue;
        }
        ptensor T = Tn_vertical[left_index];
        T.multiply_vector(lambda_tensor[left_index][3], 3);
        const auto left_op = onesite_operators[left_op_index].op;
        core::StartCorrelation_MF(correlation_T, T, left_op, 1);
        core::StartCorrelation_MF(correlation_norm, T, op_identity[left_index],
                                  1);

        int right_index = left_index;
        for (int r = 0; r < r_max; ++r) {
          right_index = lattice.top(right_index);
          T = Tn_vertical[right_index];
          T.multiply_vector(lambda_tensor[right_index][1], 1);
          double norm = std::real(core::FinishCorrelation_MF(
              correlation_norm, T, op_identity[right_index], 1));
          for (auto right_ilop : r_ops[left_ilop]) {
            int right_op_index = siteoperator_index(right_index, right_ilop);
            if (right_op_index < 0) {
              continue;
            }
            const auto right_op = onesite_operators[right_op_index].op;
            auto val =
                core::FinishCorrelation_MF(correlation_T, T, right_op, 1) /
                norm;
            correlations.push_back(Correlation{left_index, 0, r + 1, left_ilop,
                                               right_ilop, std::real(val),
                                               std::imag(val)});
          }

          core::Transfer_MF(correlation_T, Tn_vertical[right_index], 1);
          core::Transfer_MF(correlation_norm, Tn_vertical[right_index], 1);
        }
      }
    }
  }

  time_observable += timer.elapsed();
  return correlations;
}

template <class ptensor>
void iTPS<ptensor>::save_correlation(
    std::vector<Correlation> const &correlations) {
  if (mpirank != 0) {
    return;
  }
  std::string filename = outdir + "/correlation.dat";
  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "    Save long-range correlations to " << filename
              << std::endl;
  }
  std::ofstream ofs(filename.c_str());
  ofs << std::scientific
      << std::setprecision(std::numeric_limits<double>::max_digits10);
  ofs << "# $1: left_op\n";
  ofs << "# $2: left_site\n";
  ofs << "# $3: right_op\n";
  ofs << "# $4: right_dx\n";
  ofs << "# $5: right_dy\n";
  ofs << "# $6: real\n";
  ofs << "# $7: imag\n";
  ofs << std::endl;
  for (auto const &cor : correlations) {
    ofs << cor.left_op << " " << cor.left_index << " " << cor.right_op << " "
        << cor.right_dx << " " << cor.right_dy << " " << cor.real << " "
        << cor.imag << " " << std::endl;
  }
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
