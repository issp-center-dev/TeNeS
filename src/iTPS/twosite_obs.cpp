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

#include "core/contract.hpp"

namespace tenes {
namespace itps {

template <class ptensor>
auto iTPS<ptensor>::measure_twosite()
    -> std::vector<std::map<Bond, typename iTPS<ptensor>::tensor_type>> {
  Timer<> timer;

  const bool is_TPO = peps_parameters.calcmode ==
                      PEPS_Parameters::CalculationMode::finite_temperature;
  const bool is_mf = peps_parameters.MeanField_Env;
  const bool use_old =
      peps_parameters.contraction_mode == PEPS_Parameters::ContractionMode::old;
  const int nlops = num_twosite_operators;
  std::vector<std::map<Bond, tensor_type>> ret(nlops);

  const int nmax = is_mf ? 5 : 4;

  std::map<Bond, tensor_type> norms;

  for (const auto &op : twosite_operators) {
    const int source = op.source_site;
    const int dx = op.dx[0];
    const int dy = op.dy[0];

    const int ncol = std::abs(dx) + 1;
    const int nrow = std::abs(dy) + 1;
    if (ncol > nmax || nrow > nmax) {
      if (peps_parameters.contraction_mode ==
              PEPS_Parameters::ContractionMode::force_static ||
          use_old) {
        std::cerr << "Warning: When contraction_mode = force_static or old, "
                     "too long-ranged operator does not supported."
                  << std::endl;
        std::cerr << "Skipped: group = " << op.group << " (dx = " << dx
                  << ", dy = " << dy << ")" << std::endl;
        continue;
      }
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
    std::vector<int> shape_types(nrow * ncol);
    std::vector<ptensor> boundaries;

    int source_col, source_row, target_col, target_row;

    if (dx >= 0) {
      source_col = 0;
      target_col = ncol - 1;
    } else {
      source_col = ncol - 1;
      target_col = 0;
    }
    if (dy >= 0) {
      source_row = nrow - 1;
      target_row = 0;
    } else {
      source_row = 0;
      target_row = nrow - 1;
    }

    if (is_mf) {
      int iboundary = 0;
      const int nboundary = 2 * (ncol + nrow - 2);
      boundaries.reserve(nboundary);

      for (int row = 0; row < nrow; ++row) {
        for (int col = 0; col < ncol; ++col) {
          const int index =
              lattice.other(source, col - source_col, source_row - row);
          indices[row][col] = index;
          shape_types[row * ncol + col] = tensor_shape_types[index];
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

    TensorNetworkContractor<ptensor> tnc;
    if (!use_old) {
      auto tnc_it =
          contraction_paths.find(std::forward_as_tuple(nrow, ncol, shape_types));
      if (tnc_it == contraction_paths.end()) {
        auto res = contraction_paths.emplace(
            std::piecewise_construct, std::forward_as_tuple(nrow, ncol, shape_types),
            std::forward_as_tuple(nrow, ncol, is_TPO, is_mf));
        tnc_it = res.first;
        if (peps_parameters.contraction_mode ==
            PEPS_Parameters::ContractionMode::force_dynamic) {
          tnc_it->second.optimize(shape_types, tensor_shape_dims, peps_parameters.CHI);
        } else if (peps_parameters.contraction_mode ==
                   PEPS_Parameters::ContractionMode::force_static) {
          std::vector<int> path = default_path(nrow, ncol, is_TPO, is_mf);
          assert(!path.empty());
          tnc_it->second.load_path(path);
        } else {
          std::vector<int> path = default_path(nrow, ncol, is_TPO, is_mf);
          if (path.empty()) {
            tnc_it->second.optimize(shape_types, tensor_shape_dims, peps_parameters.CHI);
          } else {
            bool use_default_path = true;
            const int first_type = shape_types[0];
            for (size_t d = 1; d < 4; ++d) {
              if (tensor_shape_dims[first_type][d] !=
                  tensor_shape_dims[first_type][0]) {
                use_default_path = false;
                break;
              }
            }
            if (use_default_path) {
              for (const auto &stype : shape_types) {
                if (stype != first_type) {
                  use_default_path = false;
                  break;
                }
              }
            }
            if (use_default_path) {
              tnc_it->second.load_path(path);
            } else {
              tnc_it->second.optimize(shape_types, tensor_shape_dims, peps_parameters.CHI);
            }
          }  // end use_default_path
        }    // end contraction_mode
      }      // end tnc_it == contraction_paths.end()
      tnc = tnc_it->second;
    }  // end !use_old

    const auto norm_key = Bond{indices[nrow - 1][0], nrow - 1, ncol - 1};
    if (norms.count(norm_key) == 0) {
      if (use_old) {
        norms[norm_key] =
            core::Contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_, is_TPO, is_mf);
      } else {
        norms[norm_key] = tnc.contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_);
      }
    }
    auto norm = norms[norm_key];

    tensor_type value = 0.0;
    if (op.ops_indices.empty()) {
      if (nrow * ncol == 2) {
        if (nrow == 2) {
          const int top = indices[0][0];
          const int bottom = indices[1][0];
          ptensor o =
              (top == source ? op.op
                             : mptensor::transpose(op.op, {1, 0, 3, 2}));
          value = core::Contract_two_sites_vertical_op12(
              C1[top], C2[top], C3[bottom], C4[bottom], eTt[top], eTr[top],
              eTr[bottom], eTb[bottom], eTl[bottom], eTl[top], *(Tn_[0][0]),
              *(Tn_[1][0]), o, is_TPO, is_mf);
        } else {  // ncol == 2
          const int left = indices[0][0];
          const int right = indices[0][1];
          ptensor o =
              (left == source ? op.op
                              : mptensor::transpose(op.op, {1, 0, 3, 2}));

          value = core::Contract_two_sites_horizontal_op12(
              C1[left], C2[right], C3[right], C4[left], eTt[left], eTt[right],
              eTr[right], eTb[right], eTb[left], eTl[left], *(Tn_[0][0]),
              *(Tn_[0][1]), o, is_TPO, is_mf);
        }
      } else {
        ptensor U, VT;
        std::vector<double> s;
        mptensor::svd(op.op, {0, 2}, {1, 3}, U, s, VT);
        const int ns = s.size();
        for (int is = 0; is < ns; ++is) {
          ptensor source_op =
              reshape(slice(U, 2, is, is + 1), {U.shape()[0], U.shape()[0]});
          ptensor target_op =
              reshape(slice(VT, 0, is, is + 1), {VT.shape()[1], VT.shape()[1]});
          if (!use_old) {
            op_[source_row][source_col] = &source_op;
            op_[target_row][target_col] = &target_op;
            auto localvalue = core::Contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_,
                                             op_, is_TPO, is_mf);
            value += localvalue * s[is];
          } else {
            std::map<std::tuple<int, int>, ptensor const *> ops;
            ops[std::forward_as_tuple(source_row, source_col)] = &source_op;
            ops[std::forward_as_tuple(target_row, target_col)] = &target_op;
            auto localvalue =
                tnc.contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_, ops);
            value += localvalue * s[is];
          }
        }
      }
    } else {  // op.ops_indices.size() == 2
      ptensor const *op_source =
          &(onesite_operators[siteoperator_index(op.source_site,
                                                 op.ops_indices[0])]
                .op);
      const int target_site = lattice.other(op.source_site, dx, dy);
      ptensor const *op_target = &(
          onesite_operators[siteoperator_index(target_site, op.ops_indices[1])]
              .op);
      if (use_old) {
        op_[source_row][source_col] = op_source;
        op_[target_row][target_col] = op_target;
        auto localvalue =
            core::Contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_, is_TPO, is_mf);
        value += localvalue;
      } else {
        std::map<std::tuple<int, int>, ptensor const *> ops;
        ops[std::forward_as_tuple(source_row, source_col)] = op_source;
        ops[std::forward_as_tuple(target_row, target_col)] = op_target;
        auto localvalue = tnc.contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_, ops);
        value += localvalue;
      }
    }
    ret[op.group][{op.source_site, op.dx[0], op.dy[0]}] = value / norm;
  }  // end of loop over twosite_operators
  ret.push_back(norms);

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
void iTPS<ptensor>::save_twosite(
    std::vector<std::map<Bond, typename iTPS<ptensor>::tensor_type>> const
        &twosite_obs,
    boost::optional<double> time, std::string filename_prefix) {
  if (mpirank != 0) {
    return;
  }

  const int nlops = num_twosite_operators;
  std::string filepath = outdir + "/" + filename_prefix + "twosite_obs.dat";
  if (!time && peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "    Save twosite observables to " << filepath << std::endl;
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
    ofs << "# $" << index++ << ": source_site\n";
    ofs << "# $" << index++ << ": dx\n";
    ofs << "# $" << index++ << ": dy\n";
    ofs << "# $" << index++ << ": real\n";
    ofs << "# $" << index++ << ": imag\n";

    ofs << "# The names of op_group are the following: \n";
    for (int ilops = 0; ilops < num_twosite_operators; ++ilops) {
      ofs << "# " << ilops << ": " << twosite_operator_names[ilops] << "\n";
    }
    if (twosite_obs.size() == nlops + 1) {
      ofs << "# -1: norm\n";
    }
    ofs << std::endl;
  }
  std::ofstream ofs(filepath.c_str(), std::ios::out | std::ios::app);
  ofs << std::scientific
      << std::setprecision(std::numeric_limits<double>::max_digits10);

  for (int ilops = 0; ilops < nlops; ++ilops) {
    for (const auto &r : twosite_obs[ilops]) {
      auto bond = r.first;
      auto value = r.second;
      if (time) {
        ofs << time.get() << " ";
      }
      ofs << ilops << " " << bond.source_site << " " << bond.dx << " "
          << bond.dy << " " << std::real(value) << " " << std::imag(value)
          << std::endl;
    }
  }

  if (twosite_obs.size() == nlops + 1) {
    // includes norm
    for (const auto &r : twosite_obs[nlops]) {
      auto bond = r.first;
      auto value = r.second;
      if (time) {
        ofs << time.get() << " ";
      }
      ofs << "-1 " << bond.source_site << " " << bond.dx << " " << bond.dy
          << " " << std::real(value) << " " << std::imag(value) << std::endl;
    }
  }
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
