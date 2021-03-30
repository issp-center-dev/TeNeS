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

#ifndef TENES_SRC_ITPS_CORE_CTM_HPP_
#define TENES_SRC_ITPS_CORE_CTM_HPP_

#include <vector>

namespace tenes {

class PEPS_Parameters;
class SquareLattice;

template <class tensor>
void Calc_projector_left_block(
    const tensor &C1, const tensor &C4,
    const tensor &eT1, const tensor &eT6,
    const tensor &eT7, const tensor &eT8,
    const tensor &Tn1, const tensor &Tn4,
    const PEPS_Parameters peps_parameters, tensor &PU,
    tensor &PL);

template <class tensor>
void Calc_projector_updown_blocks(
    const tensor &C1, const tensor &C2,
    const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2,
    const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6,
    const tensor &eT7, const tensor &eT8,
    const tensor &Tn1, const tensor &Tn2,
    const tensor &Tn3, const tensor &Tn4,
    const PEPS_Parameters peps_parameters, tensor &PU,
    tensor &PL);

template <class tensor>
void Calc_Next_CTM(const tensor &C1, const tensor &C4,
                   const tensor &eT1, const tensor &eT6,
                   const tensor &PU, const tensor &PL,
                   tensor &C1_out, tensor &C4_out);

template <class tensor>
void Calc_Next_eT(const tensor &eT8, const tensor &Tn1,
                  const tensor &PU, const tensor &PL,
                  tensor &eT_out);

template <class tensor>
void Left_move(std::vector<tensor> &C1,
               const std::vector<tensor> &C2,
               const std::vector<tensor> &C3,
               std::vector<tensor> &C4,
               const std::vector<tensor> &eTt,
               const std::vector<tensor> &eTr,
               const std::vector<tensor> &eTb,
               std::vector<tensor> &eTl,
               const std::vector<tensor> &Tn, const int ix,
               const PEPS_Parameters peps_parameters, const SquareLattice lattice);

template <class tensor>
void Right_move(const std::vector<tensor> &C1,
                std::vector<tensor> &C2,
                std::vector<tensor> &C3,
                const std::vector<tensor> &C4,
                const std::vector<tensor> &eTt,
                std::vector<tensor> &eTr,
                const std::vector<tensor> &eTb,
                const std::vector<tensor> &eTl,
                const std::vector<tensor> &Tn, const int ix,
                const PEPS_Parameters peps_parameters, const SquareLattice lattice);

template <class tensor>
void Top_move(std::vector<tensor> &C1,
              std::vector<tensor> &C2,
              const std::vector<tensor> &C3,
              const std::vector<tensor> &C4,
              std::vector<tensor> &eTt,
              const std::vector<tensor> &eTr,
              const std::vector<tensor> &eTb,
              const std::vector<tensor> &eTl,
              const std::vector<tensor> &Tn, const int iy,
              const PEPS_Parameters peps_parameters, const SquareLattice lattice);

template <class tensor>
void Bottom_move(const std::vector<tensor> &C1,
                 const std::vector<tensor> &C2,
                 std::vector<tensor> &C3,
                 std::vector<tensor> &C4,
                 const std::vector<tensor> &eTt,
                 const std::vector<tensor> &eTr,
                 std::vector<tensor> &eTb,
                 const std::vector<tensor> &eTl,
                 const std::vector<tensor> &Tn, const int iy,
                 const PEPS_Parameters peps_parameters, const SquareLattice lattice);

template <class tensor>
bool Check_Convergence_CTM(const std::vector<tensor> &C1,
                           const std::vector<tensor> &C2,
                           const std::vector<tensor> &C3,
                           const std::vector<tensor> &C4,
                           const std::vector<tensor> &C1_old,
                           const std::vector<tensor> &C2_old,
                           const std::vector<tensor> &C3_old,
                           const std::vector<tensor> &C4_old,
                           const PEPS_Parameters peps_parameters,
                           const SquareLattice lattice, double &sig_max);

template <class tensor>
int Calc_CTM_Environment(
    std::vector<tensor> &C1, std::vector<tensor> &C2,
    std::vector<tensor> &C3, std::vector<tensor> &C4,
    std::vector<tensor> &eTt, std::vector<tensor> &eTr,
    std::vector<tensor> &eTb, std::vector<tensor> &eTl,
    const std::vector<tensor> &Tn,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice,
    bool initialize = true);

}  // end of namespace tenes

#endif  // TENES_SRC_ITPS_CORE_CTM_HPP_
