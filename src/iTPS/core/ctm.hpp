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

#include <complex>
#include <vector>

#include "../../tensor.hpp"

namespace tenes {
class SquareLattice;

namespace itps {
class PEPS_Parameters;

namespace core {

template <class tensor>
void Calc_projector_left_block(const tensor &C1, const tensor &C4,
                               const tensor &eT1, const tensor &eT6,
                               const tensor &eT7, const tensor &eT8,
                               const tensor &Tn1, const tensor &Tn4,
                               const PEPS_Parameters peps_parameters,
                               tensor &PU, tensor &PL);
template <class tensor>
void Calc_projector_left_block_single(const tensor &C1, const tensor &C4,
                               const tensor &eT1, const tensor &eT6,
                               const tensor &eT7, const tensor &eT8,
                               const tensor &Tn1, const tensor &Tn4,
                               const PEPS_Parameters peps_parameters,
                               tensor &PU, tensor &PL);

  
template <class tensor>
void Calc_projector_updown_blocks(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &eT7, const tensor &eT8,
    const tensor &Tn1, const tensor &Tn2, const tensor &Tn3, const tensor &Tn4,
    const PEPS_Parameters peps_parameters, tensor &PU, tensor &PL);

template <class tensor>
void Calc_projector_updown_blocks_single(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &eT7, const tensor &eT8,
    const tensor &Tn1, const tensor &Tn2, const tensor &Tn3, const tensor &Tn4,
    const PEPS_Parameters peps_parameters, tensor &PU, tensor &PL);

  
template <class tensor>
void Calc_Next_CTM(const tensor &C1, const tensor &C4, const tensor &eT1,
                   const tensor &eT6, const tensor &PU, const tensor &PL,
                   tensor &C1_out, tensor &C4_out);

template <class tensor>
void Calc_Next_eT(const tensor &eT8, const tensor &Tn1, const tensor &PU,
                  const tensor &PL, tensor &eT_out);

template <class tensor>
void Calc_Next_CTM_single(const tensor &C1, const tensor &C4, const tensor &eT1,
                   const tensor &eT6, const tensor &PU, const tensor &PL,
                   tensor &C1_out, tensor &C4_out);

template <class tensor>
void Calc_Next_eT_single(const tensor &eT8, const tensor &Tn1, const tensor &PU,
                  const tensor &PL, tensor &eT_out);

  
template <class tensor>
void Left_move(std::vector<tensor> &C1, const std::vector<tensor> &C2,
               const std::vector<tensor> &C3, std::vector<tensor> &C4,
               const std::vector<tensor> &eTt, const std::vector<tensor> &eTr,
               const std::vector<tensor> &eTb, std::vector<tensor> &eTl,
               const std::vector<tensor> &Tn, const int ix,
               const PEPS_Parameters peps_parameters,
               const SquareLattice lattice);

template <class tensor>
void Right_move(const std::vector<tensor> &C1, std::vector<tensor> &C2,
                std::vector<tensor> &C3, const std::vector<tensor> &C4,
                const std::vector<tensor> &eTt, std::vector<tensor> &eTr,
                const std::vector<tensor> &eTb, const std::vector<tensor> &eTl,
                const std::vector<tensor> &Tn, const int ix,
                const PEPS_Parameters peps_parameters,
                const SquareLattice lattice);

template <class tensor>
void Top_move(std::vector<tensor> &C1, std::vector<tensor> &C2,
              const std::vector<tensor> &C3, const std::vector<tensor> &C4,
              std::vector<tensor> &eTt, const std::vector<tensor> &eTr,
              const std::vector<tensor> &eTb, const std::vector<tensor> &eTl,
              const std::vector<tensor> &Tn, const int iy,
              const PEPS_Parameters peps_parameters,
              const SquareLattice lattice);

template <class tensor>
void Bottom_move(const std::vector<tensor> &C1, const std::vector<tensor> &C2,
                 std::vector<tensor> &C3, std::vector<tensor> &C4,
                 const std::vector<tensor> &eTt, const std::vector<tensor> &eTr,
                 std::vector<tensor> &eTb, const std::vector<tensor> &eTl,
                 const std::vector<tensor> &Tn, const int iy,
                 const PEPS_Parameters peps_parameters,
                 const SquareLattice lattice);

template <class tensor>
bool Check_Convergence_CTM(
    const std::vector<tensor> &C1, const std::vector<tensor> &C2,
    const std::vector<tensor> &C3, const std::vector<tensor> &C4,
    const std::vector<tensor> &C1_old, const std::vector<tensor> &C2_old,
    const std::vector<tensor> &C3_old, const std::vector<tensor> &C4_old,
    const PEPS_Parameters peps_parameters, const SquareLattice lattice,
    double &sig_max);

/**
 * @brief Check CTM convergence from one-site RDMs.
 *
 * For every site, this routine contracts the one-site reduced density matrix
 * (RDM), gathers its distributed elements over MPI, and computes the maximum
 * elementwise distance from the previous iteration divided by the current
 * trace magnitude.  The function returns true only when that scaled distance
 * is smaller than @p epsilon.  The first iteration always returns false because
 * no previous RDM exists.
 *
 * Every path rejects non-finite elements and a non-finite or zero trace as
 * unhealthy.  With @p phase_invariant false, a non-positive real trace or a
 * sizable imaginary trace also makes the iteration unconverged.  With it
 * true, each RDM is divided by its trace phase before comparison and no
 * real-positive trace condition is imposed.  When the distance is not
 * evaluated, the output is left as NaN.  The latest RDMs are still stored in
 * @p rdm_old for the next iteration.  All MPI ranks follow the same site, row,
 * and column order and use an allreduce sum, so @p rdm_dist and the return
 * value are identical on every rank.
 *
 * @param epsilon Threshold for the trace-scaled one-site RDM distance.
 * @param rdm_dist Output distance; NaN when the distance is not evaluated.
 * @param phase_invariant Whether to remove the trace phase before comparison.
 */
template <class tensor>
bool Check_Convergence_CTM_RDM(
    const std::vector<tensor> &C1, const std::vector<tensor> &C2,
    const std::vector<tensor> &C3, const std::vector<tensor> &C4,
    const std::vector<tensor> &eTt, const std::vector<tensor> &eTr,
    const std::vector<tensor> &eTb, const std::vector<tensor> &eTl,
    const std::vector<tensor> &Tn, const SquareLattice lattice,
    std::vector<small_tensor<typename tensor::value_type>> &rdm_old,
    bool &has_rdm_old, const double epsilon, const bool is_density,
    double &rdm_dist, const bool phase_invariant = false);

/*! @brief Divide an RDM by the phase of its trace.
 *  @param[in,out] rdm RDM to phase-normalize; it must have finite elements and
 *                 a finite, nonzero trace.
 *  @return Unit complex trace phase removed from @p rdm.
 *  @throw std::runtime_error If @p rdm is not square, if an element or the
 *         trace is non-finite, or if the trace is zero.
 */
template <class value_type>
std::complex<double> normalize_rdm_phase(small_tensor<value_type> &rdm);

template <class tensor>
int Calc_CTM_Environment(std::vector<tensor> &C1, std::vector<tensor> &C2,
                         std::vector<tensor> &C3, std::vector<tensor> &C4,
                         std::vector<tensor> &eTt, std::vector<tensor> &eTr,
                         std::vector<tensor> &eTb, std::vector<tensor> &eTl,
                         const std::vector<tensor> &Tn,
                         const PEPS_Parameters peps_parameters,
                         const SquareLattice lattice, bool initialize = true);



template <class tensor>
void Left_move_single(std::vector<tensor> &C1, const std::vector<tensor> &C2,
               const std::vector<tensor> &C3, std::vector<tensor> &C4,
               const std::vector<tensor> &eTt, const std::vector<tensor> &eTr,
               const std::vector<tensor> &eTb, std::vector<tensor> &eTl,
               const std::vector<tensor> &Tn, const int ix,
               const PEPS_Parameters peps_parameters,
               const SquareLattice lattice);

template <class tensor>
void Right_move_single(const std::vector<tensor> &C1, std::vector<tensor> &C2,
                std::vector<tensor> &C3, const std::vector<tensor> &C4,
                const std::vector<tensor> &eTt, std::vector<tensor> &eTr,
                const std::vector<tensor> &eTb, const std::vector<tensor> &eTl,
                const std::vector<tensor> &Tn, const int ix,
                const PEPS_Parameters peps_parameters,
                const SquareLattice lattice);

template <class tensor>
void Top_move_single(std::vector<tensor> &C1, std::vector<tensor> &C2,
              const std::vector<tensor> &C3, const std::vector<tensor> &C4,
              std::vector<tensor> &eTt, const std::vector<tensor> &eTr,
              const std::vector<tensor> &eTb, const std::vector<tensor> &eTl,
              const std::vector<tensor> &Tn, const int iy,
              const PEPS_Parameters peps_parameters,
              const SquareLattice lattice);

template <class tensor>
void Bottom_move_single(const std::vector<tensor> &C1, const std::vector<tensor> &C2,
                 std::vector<tensor> &C3, std::vector<tensor> &C4,
                 const std::vector<tensor> &eTt, const std::vector<tensor> &eTr,
                 std::vector<tensor> &eTb, const std::vector<tensor> &eTl,
                 const std::vector<tensor> &Tn, const int iy,
                 const PEPS_Parameters peps_parameters,
                 const SquareLattice lattice);

/*! @brief Calculate a CTM environment for density tensors.
 *
 *  Every convergence path checks for finite RDM elements and a finite,
 *  nonzero trace.  When @p phase_invariant is true, convergence comparison
 *  divides each RDM by its trace phase and does not require the trace to be
 *  real and positive.
 *
 *  @param[in,out] C1,C2,C3,C4 Corner transfer matrices.
 *  @param[in,out] eTt,eTr,eTb,eTl Edge tensors.
 *  @param[in] Tn Density tensors for the unit cell.
 *  @param[in] peps_parameters CTM hyperparameters.
 *  @param[in] lattice Unit-cell geometry.
 *  @param[in] initialize Whether to initialize rather than reuse the input CTM.
 *  @param[in] phase_invariant Whether to compare RDMs modulo trace phase.
 *  @return Number of CTM iterations performed.
 */
template <class tensor>
int Calc_CTM_Environment_density(
    std::vector<tensor> &C1, std::vector<tensor> &C2, std::vector<tensor> &C3,
    std::vector<tensor> &C4, std::vector<tensor> &eTt, std::vector<tensor> &eTr,
    std::vector<tensor> &eTb, std::vector<tensor> &eTl,
    const std::vector<tensor> &Tn, const PEPS_Parameters peps_parameters,
    const SquareLattice lattice, bool initialize = true,
    bool phase_invariant = false);

template <class tensor>
std::vector<tensor> Make_single_tensor_density(const std::vector<tensor> &Tn);
}  // end of namespace core
}  // namespace itps
}  // namespace tenes

#endif  // TENES_SRC_ITPS_CORE_CTM_HPP_
