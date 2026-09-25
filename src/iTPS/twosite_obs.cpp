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

#include <algorithm>
#include <cassert>
#include <iomanip>

#include "iTPS.hpp"

#include "../fermion/fops.hpp"
#include "../fermion/relay.hpp"
#include "../fermion/reduced_measure.hpp"
#include "../tensor.hpp"

#include "../printlevel.hpp"
#include "../timer.hpp"

#include "core/contract.hpp"

namespace tenes::itps {

template <class ptensor>
auto iTPS<ptensor>::measure_twosite()
    -> std::vector<std::map<Bond, typename iTPS<ptensor>::tensor_type>> {
  validate_fermion_ctm_measurement();

  Timer<> timer;
  ScopedTimer scoped_timer("measure/twosite");

  const bool is_TPO = peps_parameters.calcmode ==
                      PEPS_Parameters::CalculationMode::finite_temperature;
  const bool is_mf = peps_parameters.MeanField_Env;
  const int nlops = num_twosite_operators;
  std::vector<std::map<Bond, tensor_type>> ret(nlops);

  constexpr int nmax = 4;

  std::map<Bond, tensor_type> norms;

  for (const auto &op : twosite_operators) {
    const int source = op.source_site;
    const int dx = op.dx[0];
    const int dy = op.dy[0];

    const int ncol = std::abs(dx) + 1;
    const int nrow = std::abs(dy) + 1;
    if (ncol > nmax || nrow > nmax) {
      std::cerr
          << "Warning: now version of TeNeS does not support too long-ranged "
             "operator"
          << std::endl;
      std::cerr << "group = " << op.group << " (dx = " << dx << ", dy = " << dy
                << ")" << std::endl;
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

    if (peps_parameters.MeanField_Env) {
      int iboundary = 0;
      const int nboundary =
          nrow * ncol - std::max(nrow - 2, 0) * std::max(ncol - 2, 0);
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

    const bool is_fermion_longrange_density =
        finfo.enabled && !is_TPO && nrow * ncol != 2;
    std::vector<std::vector<tenes::fermion::ftensor<ptensor>>> fTn;
    std::vector<std::vector<ptensor>> reduced;
    std::vector<std::vector<const ptensor *>> reduced_ptr;
    std::vector<ptensor> delta_C;
    std::vector<ptensor> delta_eTt;
    std::vector<ptensor> delta_eTr;
    std::vector<ptensor> delta_eTb;
    std::vector<ptensor> delta_eTl;
    if (is_fermion_longrange_density) {
      fTn.resize(nrow);
      reduced.resize(nrow);
      reduced_ptr.assign(nrow, std::vector<const ptensor *>(ncol, nullptr));
      for (int row = 0; row < nrow; ++row) {
        fTn[row].reserve(ncol);
        reduced[row].reserve(ncol);
        for (int col = 0; col < ncol; ++col) {
          fTn[row].push_back(tenes::fermion::wrap_Tn(*(Tn_[row][col]), finfo,
                                                     indices[row][col]));
          reduced[row].push_back(
              tenes::fermion::build_reduced_op(fTn[row].back()));
          reduced_ptr[row][col] = &reduced[row].back();
        }
      }
      if (is_mf) {
        delta_C.assign(4, tenes::fermion::make_delta_corner<ptensor>(comm));
        delta_eTt.reserve(ncol);
        delta_eTb.reserve(ncol);
        for (int col = 0; col < ncol; ++col) {
          delta_eTt.push_back(tenes::fermion::make_delta_edge<ptensor>(
              static_cast<int>(fTn[0][col].shape()[1]), comm));
          delta_eTb.push_back(tenes::fermion::make_delta_edge<ptensor>(
              static_cast<int>(fTn[nrow - 1][col].shape()[3]), comm));
          eTt_[col] = &delta_eTt.back();
          eTb_[col] = &delta_eTb.back();
        }
        delta_eTl.reserve(nrow);
        delta_eTr.reserve(nrow);
        for (int row = 0; row < nrow; ++row) {
          delta_eTl.push_back(tenes::fermion::make_delta_edge<ptensor>(
              static_cast<int>(fTn[row][0].shape()[0]), comm));
          delta_eTr.push_back(tenes::fermion::make_delta_edge<ptensor>(
              static_cast<int>(fTn[row][ncol - 1].shape()[2]), comm));
          eTl_[row] = &delta_eTl.back();
          eTr_[row] = &delta_eTr.back();
        }
        for (int i = 0; i < 4; ++i) {
          C_[i] = &delta_C[i];
        }
      }
    }

    const auto norm_key = Bond{indices[nrow - 1][0], nrow - 1, ncol - 1};
    if (norms.count(norm_key) == 0) {
      if (finfo.enabled && !is_TPO && nrow * ncol == 2) {
        if (nrow == 2) {
          const int top = indices[0][0];
          const int bottom = indices[1][0];
          const auto fTop = tenes::fermion::wrap_Tn(*(Tn_[0][0]), finfo, top);
          const auto fBottom =
              tenes::fermion::wrap_Tn(*(Tn_[1][0]), finfo, bottom);
          if (is_mf) {
            // Mean field: Tn_ are the lambda-dressed boundary copies, so the
            // single-layer graded contraction already closes the window.
            norms[norm_key] = tenes::fermion::contract_pair_MF(
                tenes::fermion::build_pair_state(
                    fTop, fBottom,
                    tenes::fermion::reduced_pair_direction::vertical));
          } else {
            const auto halves = tenes::fermion::build_reduced_identity_halves(
                fTop, fBottom,
                tenes::fermion::reduced_pair_direction::vertical);
            norms[norm_key] =
                tenes::fermion::contract_reduced_pair_halves_density_CTM(
                    C1[top], C2[top], C3[bottom], C4[bottom], eTt[top],
                    eTr[top], eTr[bottom], eTb[bottom], eTl[bottom], eTl[top],
                    halves);
          }
        } else {
          const int left = indices[0][0];
          const int right = indices[0][1];
          const auto fLeft = tenes::fermion::wrap_Tn(*(Tn_[0][0]), finfo, left);
          const auto fRight =
              tenes::fermion::wrap_Tn(*(Tn_[0][1]), finfo, right);
          if (is_mf) {
            norms[norm_key] = tenes::fermion::contract_pair_MF(
                tenes::fermion::build_pair_state(
                    fLeft, fRight,
                    tenes::fermion::reduced_pair_direction::horizontal));
          } else {
            const auto halves = tenes::fermion::build_reduced_identity_halves(
                fLeft, fRight,
                tenes::fermion::reduced_pair_direction::horizontal);
            norms[norm_key] =
                tenes::fermion::contract_reduced_pair_halves_density_CTM(
                    C1[left], C2[right], C3[right], C4[left], eTt[left],
                    eTt[right], eTr[right], eTb[right], eTb[left], eTl[left],
                    halves);
          }
        }
      } else if (is_fermion_longrange_density) {
        norms[norm_key] = core::Contract_density_CTM(C_, eTt_, eTr_, eTb_, eTl_,
                                                     reduced_ptr, op_);
      } else if (is_mf) {
        norms[norm_key] = core::Contract_iTPS_MF(Tn_, op_);
      } else if (is_TPO) {
        norms[norm_key] =
            core::Contract_density_CTM(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_);
      } else {
        norms[norm_key] =
            core::Contract_iTPS_CTM(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_);
      }
    }
    auto norm = norms[norm_key];

    const int target_site = lattice.other(op.source_site, dx, dy);
    ptensor fermion_product_op;
    const ptensor *op12 = &op.op;
    bool use_fermion_product_op = false;
    if (finfo.enabled && !is_TPO && !op.ops_indices.empty()) {
      const int opA_index =
          siteoperator_index(op.source_site, op.ops_indices[0]);
      const int opB_index = siteoperator_index(target_site, op.ops_indices[1]);
      const auto pA = onesite_parity[opA_index];
      const auto pB = onesite_parity[opB_index];
      if (pA == tenes::fermion::op_parity::mixed ||
          pB == tenes::fermion::op_parity::mixed) {
        throw std::runtime_error("fermion ops form contains mixed parity");
      }
      if (pA != pB) {
        throw std::runtime_error(
            "fermion ops form combines one-site operators of different parity");
      }
      fermion_product_op = tenes::fermion::product_twosite_op(
          onesite_operators[opA_index].op, onesite_operators[opB_index].op,
          finfo.phys[op.source_site], finfo.phys[target_site],
          pB == tenes::fermion::op_parity::odd);
      op12 = &fermion_product_op;
      use_fermion_product_op = true;
    }

    tensor_type value = 0.0;
    if (op.ops_indices.empty() || use_fermion_product_op) {
      if (nrow * ncol == 2) {
        if (nrow == 2) {
          const int top = indices[0][0];
          const int bottom = indices[1][0];
          if (finfo.enabled && !is_TPO) {
            const int target = top == source ? bottom : top;
            auto o = tenes::fermion::wrap_twosite_gate(
                *op12, finfo.phys[source], finfo.phys[target]);
            if (top != source) {
              // Graded transpose: carries the Fock reordering sign
              // |n_B n_A> = (-1)^{n_A n_B} |n_A n_B> on both leg pairs.
              o = tenes::fermion::transpose(o, mptensor::Axes(1, 0, 3, 2));
            }
            const auto fTop = tenes::fermion::wrap_Tn(*(Tn_[0][0]), finfo, top);
            const auto fBottom =
                tenes::fermion::wrap_Tn(*(Tn_[1][0]), finfo, bottom);
            if (is_mf) {
              value = tenes::fermion::contract_pair_MF(
                  tenes::fermion::build_pair_state(
                      fTop, fBottom,
                      tenes::fermion::reduced_pair_direction::vertical),
                  o);
            } else {
              // The window rows already carry geometric roles (row 0 =
              // upper), and the infinite lattice has no boundary, so bonds
              // whose target wraps around the unit cell need no special
              // ordering.
              const auto halves = tenes::fermion::build_reduced_pair_halves(
                  fTop, fBottom, o,
                  tenes::fermion::reduced_pair_direction::vertical);
              value = tenes::fermion::contract_reduced_pair_halves_density_CTM(
                  C1[top], C2[top], C3[bottom], C4[bottom], eTt[top], eTr[top],
                  eTr[bottom], eTb[bottom], eTl[bottom], eTl[top], halves);
            }
          } else {
            ptensor o =
                (top == source ? *op12
                               : mptensor::transpose(*op12, {1, 0, 3, 2}));
            value = core::Contract_two_sites_vertical_op12(
                C1[top], C2[top], C3[bottom], C4[bottom], eTt[top], eTr[top],
                eTr[bottom], eTb[bottom], eTl[bottom], eTl[top], *(Tn_[0][0]),
                *(Tn_[1][0]), o, is_TPO, is_mf);
          }
          // value = peps_parameters.MeanField_Env
          //             ? core::Contract_two_sites_vertical_op12_MF(
          //                   *(Tn_[0][0]), *(Tn_[1][0]), o)
          //             : core::Contract_two_sites_vertical_op12(
          //                   C1[top], C2[top], C3[bottom], C4[bottom],
          //                   eTt[top], eTr[top], eTr[bottom], eTb[bottom],
          //                   eTl[bottom], eTl[top], Tn[top], Tn[bottom], o);
        } else {  // ncol == 2
          const int left = indices[0][0];
          const int right = indices[0][1];
          if (finfo.enabled && !is_TPO) {
            const int target = left == source ? right : left;
            auto o = tenes::fermion::wrap_twosite_gate(
                *op12, finfo.phys[source], finfo.phys[target]);
            if (left != source) {
              // Graded transpose; see the vertical branch.
              o = tenes::fermion::transpose(o, mptensor::Axes(1, 0, 3, 2));
            }
            const auto fLeft =
                tenes::fermion::wrap_Tn(*(Tn_[0][0]), finfo, left);
            const auto fRight =
                tenes::fermion::wrap_Tn(*(Tn_[0][1]), finfo, right);
            if (is_mf) {
              value = tenes::fermion::contract_pair_MF(
                  tenes::fermion::build_pair_state(
                      fLeft, fRight,
                      tenes::fermion::reduced_pair_direction::horizontal),
                  o);
            } else {
              // Window columns already carry geometric roles (col 0 = left);
              // see the vertical branch.
              const auto halves = tenes::fermion::build_reduced_pair_halves(
                  fLeft, fRight, o,
                  tenes::fermion::reduced_pair_direction::horizontal);
              value = tenes::fermion::contract_reduced_pair_halves_density_CTM(
                  C1[left], C2[right], C3[right], C4[left], eTt[left],
                  eTt[right], eTr[right], eTb[right], eTb[left], eTl[left],
                  halves);
            }
          } else {
            ptensor o =
                (left == source ? *op12
                                : mptensor::transpose(*op12, {1, 0, 3, 2}));

            value = core::Contract_two_sites_horizontal_op12(
                C1[left], C2[right], C3[right], C4[left], eTt[left], eTt[right],
                eTr[right], eTb[right], eTb[left], eTl[left], *(Tn_[0][0]),
                *(Tn_[0][1]), o, is_TPO, is_mf);
          }
          // value = peps_parameters.MeanField_Env
          //             ? core::Contract_two_sites_horizontal_op12_MF(
          //                   *(Tn_[0][0]), *(Tn_[0][1]), o)
          //             : core::Contract_two_sites_horizontal_op12(
          //                   C1[left], C2[right], C3[right], C4[left],
          //                   eTt[left], eTt[right], eTr[right], eTb[right],
          //                   eTb[left], eTl[left], Tn[left], Tn[right], o);
        }
      } else {
        if (is_fermion_longrange_density) {
          const auto wrapped_op = tenes::fermion::wrap_twosite_gate(
              *op12, finfo.phys[source], finfo.phys[target_site]);
          const auto channels = tenes::fermion::relay_channels(wrapped_op);
          const auto path = tenes::fermion::relay_path(
              tenes::fermion::window_cell{source_row, source_col},
              tenes::fermion::window_cell{target_row, target_col},
              tenes::fermion::relay_order::x_first);
          for (const auto &channel : channels) {
            auto relay_Tn = reduced_ptr;
            std::vector<ptensor> path_tensors;
            path_tensors.reserve(path.size());
            for (std::size_t i = 0; i < path.size(); ++i) {
              const auto cell = path[i];
              tenes::fermion::relay_role role =
                  tenes::fermion::relay_role::middle;
              int entry = -1;
              int exit = -1;
              if (i == 0) {
                role = tenes::fermion::relay_role::source;
                exit = tenes::fermion::relay_leg(path[i], path[i + 1]);
              } else if (i + 1 == path.size()) {
                role = tenes::fermion::relay_role::target;
                entry = tenes::fermion::relay_leg(path[i], path[i - 1]);
              } else {
                entry = tenes::fermion::relay_leg(path[i], path[i - 1]);
                exit = tenes::fermion::relay_leg(path[i], path[i + 1]);
              }
              path_tensors.push_back(tenes::fermion::build_relay_site(
                  fTn[cell.row][cell.col], role, entry, exit, channel));
              relay_Tn[cell.row][cell.col] = &path_tensors.back();
            }
            value += core::Contract_density_CTM(C_, eTt_, eTr_, eTb_, eTl_,
                                                relay_Tn, op_);
          }
        } else {
          ptensor U, VT;
          std::vector<double> s;
          mptensor::svd(*op12, {0, 2}, {1, 3}, U, s, VT);
          const int ns = s.size();
          for (int is = 0; is < ns; ++is) {
            ptensor source_op =
                reshape(slice(U, 2, is, is + 1), {U.shape()[0], U.shape()[0]});
            op_[source_row][source_col] = &source_op;
            ptensor target_op = reshape(slice(VT, 0, is, is + 1),
                                        {VT.shape()[1], VT.shape()[1]});
            op_[target_row][target_col] = &target_op;
            auto localvalue = core::Contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_,
                                             op_, is_TPO, is_mf);
            // auto localvalue =
            //     peps_parameters.MeanField_Env
            //         ? core::Contract_MF(Tn_, op_)
            //         : core::Contract_CTM(C_, eTt_, eTr_, eTb_, eTl_, Tn_,
            //         op_);
            value += localvalue * s[is];
          }
        }
      }
    } else {
      op_[source_row][source_col] =
          &(onesite_operators[siteoperator_index(op.source_site,
                                                 op.ops_indices[0])]
                .op);
      op_[target_row][target_col] = &(
          onesite_operators[siteoperator_index(target_site, op.ops_indices[1])]
              .op);
      auto localvalue =
          core::Contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_, is_TPO, is_mf);
      // auto localvalue =
      //     peps_parameters.MeanField_Env
      //         ? core::Contract_MF(Tn_, op_)
      //         : core::Contract_CTM(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_);
      value += localvalue;
    }
    ret[op.group][{op.source_site, op.dx[0], op.dy[0]}] =
        op.coeff * value / norm;
  }
  ret.push_back(norms);

  double norm_real_min = 1e100;
  double norm_imag_abs_max = 0.0;
  for (const auto &[bond, norm] : norms) {
    tensor_type diagnostic_norm = norm;
    if (finfo.enabled && !is_mf) {
      if (std::isfinite(std::real(norm)) && std::isfinite(std::imag(norm)) &&
          std::abs(norm) > 0.0) {
        diagnostic_norm = tensor_type(std::abs(norm));
      } else {
        norm_real_min = -std::numeric_limits<double>::infinity();
        norm_imag_abs_max = std::numeric_limits<double>::infinity();
        continue;
      }
    }
    double norm_re = std::real(diagnostic_norm);
    double norm_im = std::imag(diagnostic_norm);
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
    std::optional<double> time, std::string filename_prefix) {
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
    if (twosite_obs.size() == static_cast<std::size_t>(nlops) + 1) {
      ofs << "# -1: norm\n";
    }
    ofs << std::endl;
  }
  std::ofstream ofs(filepath.c_str(), std::ios::out | std::ios::app);
  ofs << std::scientific
      << std::setprecision(std::numeric_limits<double>::max_digits10);

  for (int ilops = 0; ilops < nlops; ++ilops) {
    for (const auto &[bond, value] : twosite_obs[ilops]) {
      if (time) {
        ofs << (*time) << " ";
      }
      ofs << ilops << " " << bond.source_site << " " << bond.dx << " "
          << bond.dy << " " << std::real(value) << " " << std::imag(value)
          << std::endl;
    }
  }

  if (twosite_obs.size() == static_cast<std::size_t>(nlops) + 1) {
    // includes norm
    for (const auto &[bond, value] : twosite_obs[nlops]) {
      if (time) {
        ofs << (*time) << " ";
      }
      ofs << "-1 " << bond.source_site << " " << bond.dx << " " << bond.dy
          << " " << std::real(value) << " " << std::imag(value) << std::endl;
    }
  }
}

template <class ptensor>
auto iTPS<ptensor>::measure_twosite_density()
    -> std::vector<std::map<Bond, typename iTPS<ptensor>::tensor_type>> {
  Timer<> timer;
  ScopedTimer scoped_timer("measure/twosite");

  const int nlops = num_twosite_operators;
  std::vector<std::map<Bond, tensor_type>> ret(nlops);

  // constexpr int nmax = 4;

  std::map<Bond, tensor_type> norms;

  for (const auto &op : twosite_operators) {
    const int dx = op.dx[0];
    const int dy = op.dy[0];

    const int ncol = std::abs(dx) + 1;
    const int nrow = std::abs(dy) + 1;
  }

  for (const auto &op : twosite_operators) {
    const int source = op.source_site;
    const int dx = op.dx[0];
    const int dy = op.dy[0];

    const int ncol = std::abs(dx) + 1;
    const int nrow = std::abs(dy) + 1;
    if (ncol * nrow != 2) {
      std::cerr
          << "Warning: now version of TeNeS does not support too long-ranged "
             "operator"
          << std::endl;
      std::cerr << "group = " << op.group << " (dx = " << dx << ", dy = " << dy
                << ")" << std::endl;
      continue;
    }

    std::vector<const ptensor *> C_(4, nullptr);
    std::vector<const ptensor *> eTt_(ncol, nullptr);
    std::vector<const ptensor *> eTr_(nrow, nullptr);
    std::vector<const ptensor *> eTb_(ncol, nullptr);
    std::vector<const ptensor *> eTl_(nrow, nullptr);

    std::vector<std::vector<const ptensor *>> Tn_(
        nrow, std::vector<const ptensor *>(ncol, nullptr));
    std::vector<std::vector<const ptensor *>> op_(
        nrow, std::vector<const ptensor *>(ncol, nullptr));

    std::vector<std::vector<int>> indices(nrow, std::vector<int>(ncol));
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

    const auto norm_key = Bond{indices[nrow - 1][0], nrow - 1, ncol - 1};
    /*
    if (norms.count(norm_key) == 0) {
      if (peps_parameters.MeanField_Env) {
        norms[norm_key] = core::Contract_MF_density(Tn_, op_);
      } else {
        norms[norm_key] = core::Contract_density(C_, eTt_, eTr_, eTb_,
        eTl_, Tn_, op_);
      }
    }
    */

    if (norms.count(norm_key) == 0) {
      if (nrow == 2) {
        const int top = indices[0][0];
        const int bottom = indices[1][0];
        norms[norm_key] = core::Contract_two_sites_vertical_density_CTM(
            C1[top], C2[top], C3[bottom], C4[bottom], eTt[top], eTr[top],
            eTr[bottom], eTb[bottom], eTl[bottom], eTl[top], Tn[top],
            Tn[bottom], op_identity[top], op_identity[bottom]);
      } else {
        const int left = indices[0][0];
        const int right = indices[0][1];
        norms[norm_key] = core::Contract_two_sites_horizontal_density_CTM(
            C1[left], C2[right], C3[right], C4[left], eTt[left], eTt[right],
            eTr[right], eTb[right], eTb[left], eTl[left], Tn[left], Tn[right],
            op_identity[left], op_identity[right]);
      }
    }
    auto norm = norms[norm_key];

    tensor_type value = 0.0;
    if (nrow * ncol == 2) {
      if (op.ops_indices.empty()) {
        if (nrow == 2) {
          const int top = indices[0][0];
          const int bottom = indices[1][0];
          ptensor o =
              (top == source ? op.op
                             : mptensor::transpose(op.op, {1, 0, 3, 2}));
          /*
          value = peps_parameters.MeanField_Env
                      ? core::Contract_two_sites_vertical_op12_MF_density(
                            *(Tn_[0][0]), *(Tn_[1][0]), o)
                      : core::Contract_two_sites_vertical_op12_density(
                            C1[top], C2[top], C3[bottom], C4[bottom], eTt[top],
                            eTr[top], eTr[bottom], eTb[bottom], eTl[bottom],
                            eTl[top], Tn[top], Tn[bottom], o);
          */
          value = core::Contract_two_sites_vertical_op12_density_CTM(
              C1[top], C2[top], C3[bottom], C4[bottom], eTt[top], eTr[top],
              eTr[bottom], eTb[bottom], eTl[bottom], eTl[top], Tn[top],
              Tn[bottom], o);

        } else {  // ncol == 2
          const int left = indices[0][0];
          const int right = indices[0][1];
          ptensor o =
              (left == source ? op.op
                              : mptensor::transpose(op.op, {1, 0, 3, 2}));
          /*
            value = peps_parameters.MeanField_Env
            ? core::Contract_two_sites_horizontal_op12_MF_density(
            *(Tn_[0][0]), *(Tn_[0][1]), o)
            : core::Contract_two_sites_horizontal_op12_density(
            C1[left], C2[right], C3[right], C4[left], eTt[left],
            eTt[right], eTr[right], eTb[right], eTb[left],
            eTl[left], Tn[left], Tn[right], o);
          */
          value = core::Contract_two_sites_horizontal_op12_density_CTM(
              C1[left], C2[right], C3[right], C4[left], eTt[left], eTt[right],
              eTr[right], eTb[right], eTb[left], eTl[left], Tn[left], Tn[right],
              o);
        }
      } else {
        if (nrow == 2) {
          const int top = indices[0][0];
          const int bottom = indices[1][0];
          const int target_site = lattice.other(op.source_site, dx, dy);

          ptensor op_t, op_b;

          if (top == source) {
            op_t = onesite_operators[siteoperator_index(op.source_site,
                                                        op.ops_indices[0])]
                       .op;
            op_b = onesite_operators[siteoperator_index(target_site,
                                                        op.ops_indices[1])]
                       .op;
          } else {
            op_t = onesite_operators[siteoperator_index(target_site,
                                                        op.ops_indices[1])]
                       .op;
            op_b = onesite_operators[siteoperator_index(op.source_site,
                                                        op.ops_indices[0])]
                       .op;
          }

          value = core::Contract_two_sites_vertical_density_CTM(
              C1[top], C2[top], C3[bottom], C4[bottom], eTt[top], eTr[top],
              eTr[bottom], eTb[bottom], eTl[bottom], eTl[top], Tn[top],
              Tn[bottom], op_t, op_b);

        } else {  // ncol == 2
          const int left = indices[0][0];
          const int right = indices[0][1];
          const int target_site = lattice.other(op.source_site, dx, dy);

          ptensor op_l, op_r;

          if (left == source) {
            op_l = onesite_operators[siteoperator_index(op.source_site,
                                                        op.ops_indices[0])]
                       .op;
            op_r = onesite_operators[siteoperator_index(target_site,
                                                        op.ops_indices[1])]
                       .op;
          } else {
            op_l = onesite_operators[siteoperator_index(target_site,
                                                        op.ops_indices[1])]
                       .op;
            op_r = onesite_operators[siteoperator_index(op.source_site,
                                                        op.ops_indices[0])]
                       .op;
          }
          value = core::Contract_two_sites_horizontal_density_CTM(
              C1[left], C2[right], C3[right], C4[left], eTt[left], eTt[right],
              eTr[right], eTb[right], eTb[left], eTl[left], Tn[left], Tn[right],
              op_l, op_r);
        }
      }
    }
    ret[op.group][{op.source_site, op.dx[0], op.dy[0]}] =
        op.coeff * value / norm;
  }
  ret.push_back(norms);

  double norm_real_min = 1e100;
  double norm_imag_abs_max = 0.0;
  for (const auto &[bond, norm] : norms) {
    double norm_re = std::real(norm);
    double norm_im = std::imag(norm);
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

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace tenes::itps
