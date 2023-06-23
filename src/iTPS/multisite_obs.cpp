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

#include <cassert>

#include "iTPS.hpp"

#include "../tensor.hpp"

#include "../printlevel.hpp"
#include "../timer.hpp"

#include "core/contract_ctm.hpp"
#include "core/contract_mf.hpp"

namespace tenes {
namespace itps {

template <class ptensor>
auto iTPS<ptensor>::measure_multisite()
    -> std::vector<std::map<Multisites, typename iTPS<ptensor>::tensor_type>> {
  Timer<> timer;

  const int nlops = num_multisite_operators;
  std::vector<std::map<Multisites, tensor_type>> ret(nlops);
  if (nlops == 0) {
    return ret;
  }

  constexpr int nmax = 4;

  std::map<Bond, tensor_type> norms;

  for (const auto &op : multisite_operators) {
    const int nothers = op.dx.size();
    const int nsites = nothers + 1;
    const int source = op.source_site;
    int mindx = 0, maxdx = 0, mindy = 0, maxdy = 0;
    for (auto dx : op.dx) {
      mindx = std::min(mindx, dx);
      maxdx = std::max(maxdx, dx);
    }
    for (auto dy : op.dy) {
      mindy = std::min(mindy, dy);
      maxdy = std::max(maxdy, dy);
    }
    const int ncol = maxdx - mindx + 1;
    const int nrow = maxdy - mindy + 1;

    if (ncol > nmax || nrow > nmax) {
      std::cerr
          << "Warning: now version of TeNeS does not support too long-ranged "
             "operator"
          << std::endl;
      continue;
    }

    std::vector<const ptensor *> C_(4, nullptr);
    std::vector<const ptensor *> eTt_(ncol, nullptr);
    std::vector<const ptensor *> eTr_(nrow, nullptr);
    std::vector<const ptensor *> eTb_(ncol, nullptr);
    std::vector<const ptensor *> eTl_(nrow, nullptr);

    /*
     * Caution: orders of tensors in unitcell and Contract_* function are
     * different
     *
     * Lattice:
     *
     *    y
     *    ^
     *    |
     *    0--> x
     *
     * Contract_*:
     *    0-->col
     *    |
     *    v
     *    row
     *
     */
    std::vector<std::vector<const ptensor *>> Tn_(
        nrow, std::vector<const ptensor *>(ncol, nullptr));
    std::vector<std::vector<const ptensor *>> op_(
        nrow, std::vector<const ptensor *>(ncol, nullptr));

    std::vector<std::vector<int>> indices(nrow, std::vector<int>(ncol));
    std::vector<ptensor> boundaries;

    const int source_col = -mindx;
    const int source_row = maxdy;

    if (peps_parameters.MeanField_Env) {
      int iboundary = 0;
      const int nboundary = 2 * (ncol + nrow - 2);
      boundaries.reserve(nboundary);

      for (int row = 0; row < nrow; ++row) {
        for (int col = 0; col < ncol; ++col) {
          const int index =
              lattice.other(source, col - source_col, source_row - row);
          indices[row][col] = index;
          op_[row][col] = &(op_identity[index]);
          if ((0 < row && row < nrow - 1) && (0 < col && col < ncol - 1)) {
            Tn_[row][col] = &(Tn[index]);
          } else {
            boundaries.push_back(Tn[index]);
            Tn_[row][col] = &(boundaries[iboundary++]);
          }
        }
      }
      assert(boundaries.size() == nboundary);

      // absorb MF ENV into center tensors on boundary
      for (int row = 0; row < nrow; ++row) {
        const_cast<ptensor *>(Tn_[row][0])
            ->multiply_vector(lambda_tensor[indices[row][0]][0], 0);
        const_cast<ptensor *>(Tn_[row][ncol - 1])
            ->multiply_vector(lambda_tensor[indices[row][ncol - 1]][2], 2);
      }
      for (int col = 0; col < ncol; ++col) {
        const_cast<ptensor *>(Tn_[0][col])
            ->multiply_vector(lambda_tensor[indices[0][col]][1], 1);
        const_cast<ptensor *>(Tn_[nrow - 1][col])
            ->multiply_vector(lambda_tensor[indices[nrow - 1][col]][3], 3);
      }
    } else {  // Use CTM
      for (int row = 0; row < nrow; ++row) {
        for (int col = 0; col < ncol; ++col) {
          const int index =
              lattice.other(source, col - source_col, source_row - row);
          indices[row][col] = index;
          op_[row][col] = &(op_identity[index]);
          Tn_[row][col] = &(Tn[index]);
        }
        eTl_[row] = &(eTl[indices[row][0]]);
        eTr_[row] = &(eTr[indices[row][ncol - 1]]);
      }
      for (int col = 0; col < ncol; ++col) {
        eTt_[col] = &(eTt[indices[0][col]]);
        eTb_[col] = &(eTb[indices[nrow - 1][col]]);
      }
      C_[0] = &(C1[indices[0][0]]);
      C_[1] = &(C2[indices[0][ncol - 1]]);
      C_[2] = &(C3[indices[nrow - 1][ncol - 1]]);
      C_[3] = &(C4[indices[nrow - 1][0]]);
    }

    const auto norm_key = Bond{indices[nrow - 1][0], nrow - 1, ncol - 1};
    if (norms.count(norm_key) == 0) {
      if (peps_parameters.MeanField_Env) {
        norms[norm_key] = core::Contract_MF(Tn_, op_);
      } else {
        norms[norm_key] = core::Contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_);
      }
    }
    auto norm = norms[norm_key];

    tensor_type value = 0.0;
    if (op.ops_indices.empty()) {
      throw std::runtime_error(
          "Empty op.ops_indices is not supported for multisites observable");
    } else {
      op_[source_row][source_col] =
          &(onesite_operators[siteoperator_index(op.source_site,
                                                 op.ops_indices[0])]
                .op);
      for (int i = 0; i < nothers; ++i) {
        const int dx = op.dx[i];
        const int dy = op.dy[i];
        const int target_site = lattice.other(op.source_site, dx, dy);
        op_[source_row - dy][source_col + dx] =
            &(onesite_operators[siteoperator_index(target_site,
                                                   op.ops_indices[i + 1])]
                  .op);
      }
      auto localvalue =
          peps_parameters.MeanField_Env
              ? core::Contract_MF(Tn_, op_)
              : core::Contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_);
      value += localvalue;
    }
    ret[op.group][{op.source_site, op.dx, op.dy}] = value / norm;
  }
  // ret.push_back(norms);

  double norm_real_min = 1e100;
  double norm_imag_abs_max = 0.0;
  for (auto &r : norms) {
    double norm_re = std::real(r.second);
    double norm_im = std::imag(r.second);
    norm_real_min = std::min(norm_re, norm_real_min);
    norm_imag_abs_max = std::max(std::abs(norm_im), norm_imag_abs_max);
  }
  if (mpirank == 0) {
    if (norm_real_min < 0.0) {
      std::cerr << "WARNING: Norm is negative [min(real(NORM)) = "
                << norm_real_min << "].\n";
      std::cerr << "HINT: Increase the bond dimension of CTM." << std::endl;
    }
    if (norm_imag_abs_max > 1.0e-6) {
      std::cerr << "WARNING: Norm is not real [max(abs(imag(NORM))) = "
                << norm_imag_abs_max << " > 1e-6].\n";
      std::cerr << "HINT: Increase the bond dimension of CTM." << std::endl;
    }
  }

  time_observable += timer.elapsed();
  return ret;
}

template <class ptensor>
void iTPS<ptensor>::save_multisite(
    std::vector<std::map<Multisites, typename iTPS<ptensor>::tensor_type>> const
        &multisite_obs,
    boost::optional<double> time, std::string filename_prefix) {
  if (mpirank != 0) {
    return;
  }

  const int nlops = num_multisite_operators;
  std::map<int, std::string> filepath;
  for (auto nsites : multisite_operator_nsites_set) {
    std::stringstream ss;
    ss << outdir << "/" << filename_prefix << "multisite_obs_" << nsites
       << ".dat";
    filepath[nsites] = ss.str();
    if (!time && peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "    Save " << nsites << "-site observables to " << ss.str() <<
      std::endl;
    }
  }

  static bool first_time = true;
  if (first_time) {
    first_time = false;
    for (auto nsites : multisite_operator_nsites_set) {
      std::ofstream ofs(filepath[nsites].c_str());
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
      ofs << "# $" << index++ << ": source_site\n";
      for (int other = 1; other < nsites; ++other) {
        ofs << "# $" << index++ << ": dx[" << other << "]\n";
        ofs << "# $" << index++ << ": dy[" << other << "]\n";
      }
      ofs << "# $" << index++ << ": real\n";
      ofs << "# $" << index++ << ": imag\n";

      ofs << "# The names of op_group are the following: \n";
      for (int ilops = 0; ilops < num_multisite_operators; ++ilops) {
        if (multisite_operator_nsites[ilops] == nsites) {
          ofs << "# " << ilops << ": " << multisite_operator_names[ilops]
              << "\n";
        }
      }
      ofs << std::endl;
    }
  }

  std::map<int, std::ofstream> ofs;
  for (auto nsites : multisite_operator_nsites_set) {
    ofs[nsites] =
        std::ofstream(filepath[nsites].c_str(), std::ios::out | std::ios::app);
    ofs[nsites] << std::scientific
                << std::setprecision(std::numeric_limits<double>::max_digits10);
  }

  for (int ilops = 0; ilops < nlops; ++ilops) {
    const int nsites = multisite_operator_nsites[ilops];
    for (const auto &r : multisite_obs[ilops]) {
      auto multi = r.first;
      auto value = r.second;
      if (time) {
        ofs[nsites] << time.get() << " ";
      }
      ofs[nsites] << ilops << " " << multi.source_site;
      for (int other = 0; other < nsites-1; ++other) {
        ofs[nsites] << " " << multi.dx[other] << " " << multi.dy[other];
      }
      ofs[nsites] << " " << std::real(value) << " " << std::imag(value)
                  << std::endl;
    }
  }
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
