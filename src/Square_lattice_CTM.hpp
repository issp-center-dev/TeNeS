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

#ifndef TENES_SRC_SQUARE_LATTICE_HPP_
#define TENES_SRC_SQUARE_LATTICE_HPP_

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mptensor/complex.hpp>
#include <mptensor/tensor.hpp>
#include <numeric>
#include <vector>

#include "Lattice.hpp"
#include "PEPS_Basics.hpp"
#include "PEPS_Parameters.hpp"
#include "mpi.hpp"
#include "printlevel.hpp"

namespace tenes {

/*
 * corner and edge tensor
 *
 * C1 t C2
 * l  .  r
 * C4 b C3
 */

using namespace mptensor;
template <template <typename> class Matrix, typename C>
void Left_move(std::vector<Tensor<Matrix, C>> &C1,
               const std::vector<Tensor<Matrix, C>> &C2,
               const std::vector<Tensor<Matrix, C>> &C3,
               std::vector<Tensor<Matrix, C>> &C4,
               const std::vector<Tensor<Matrix, C>> &eTt,
               const std::vector<Tensor<Matrix, C>> &eTr,
               const std::vector<Tensor<Matrix, C>> &eTb,
               std::vector<Tensor<Matrix, C>> &eTl,
               const std::vector<Tensor<Matrix, C>> &Tn, const int ix,
               const PEPS_Parameters peps_parameters, const Lattice lattice) {
  /* Do one step left move absoving X=ix column
     part of C1, C4, eTl will be modified */

  std::vector<Tensor<Matrix, C>> PUs, PLs;
  PUs.resize(lattice.LY);
  PLs.resize(lattice.LY);
  int i, j, k, l;
  for (int iy = 0; iy < lattice.LY; ++iy) {
    i = lattice.index(ix, iy);
    j = lattice.right(i);
    k = lattice.bottom(j);
    l = lattice.left(k);

    if (peps_parameters.CTM_Projector_corner) {
      Calc_projector_left_block(C1[i], C4[l], eTt[i], eTb[l], eTl[l], eTl[i],
                                Tn[i], Tn[l], peps_parameters, PUs[iy],
                                PLs[iy]);
    } else {
      Calc_projector_updown_blocks(C1[i], C2[j], C3[k], C4[l], eTt[i], eTt[j],
                                   eTr[j], eTr[k], eTb[k], eTb[l], eTl[l],
                                   eTl[i], Tn[i], Tn[j], Tn[k], Tn[l],
                                   peps_parameters, PUs[iy], PLs[iy]);
    }
  }
  // update
  std::vector<Tensor<Matrix, C>> C1_bak(lattice.N_UNIT), C4_bak(lattice.N_UNIT),
      eTl_bak(lattice.N_UNIT);
  for (int num = 0; num < lattice.N_UNIT; num++) {
    C1_bak[num] = C1[num];
    C4_bak[num] = C4[num];
    eTl_bak[num] = eTl[num];
  }
  int iy_up, iy_down;
  for (int iy = 0; iy < lattice.LY; ++iy) {
    i = lattice.index(ix, iy);
    j = lattice.right(i);
    k = lattice.bottom(j);
    l = lattice.left(k);
    iy_up = (iy + 1) % lattice.LY;
    iy_down = (iy - 1 + lattice.LY) % lattice.LY;

    Calc_Next_CTM(C1_bak[i], C4_bak[l], eTt[i], eTb[l], PUs[iy_up],
                  PLs[iy_down], C1[j], C4[k]);
    Calc_Next_eT(eTl_bak[i], Tn[i], PUs[iy], PLs[iy_up], eTl[j]);
    Calc_Next_eT(eTl_bak[l], Tn[l], PUs[iy_down], PLs[iy], eTl[k]);
  }
}

template <template <typename> class Matrix, typename C>
void Right_move(const std::vector<Tensor<Matrix, C>> &C1,
                std::vector<Tensor<Matrix, C>> &C2,
                std::vector<Tensor<Matrix, C>> &C3,
                const std::vector<Tensor<Matrix, C>> &C4,
                const std::vector<Tensor<Matrix, C>> &eTt,
                std::vector<Tensor<Matrix, C>> &eTr,
                const std::vector<Tensor<Matrix, C>> &eTb,
                const std::vector<Tensor<Matrix, C>> &eTl,
                const std::vector<Tensor<Matrix, C>> &Tn, const int ix,
                const PEPS_Parameters peps_parameters, const Lattice lattice) {
  /*
    Do one step right move absorbing X=ix column
    part of C2, C3, eTr will be modified
  */
  std::vector<Tensor<Matrix, C>> PUs, PLs;
  PUs.resize(lattice.LY);
  PLs.resize(lattice.LY);
  int i, j, k, l;
  for (int iy = 0; iy < lattice.LY; ++iy) {
    k = lattice.index(ix, iy);
    l = lattice.left(k);
    i = lattice.top(l);
    j = lattice.right(i);

    if (peps_parameters.CTM_Projector_corner) {
      Calc_projector_left_block(C3[k], C2[j], eTb[k], eTt[j], eTr[j], eTr[k],
                                transpose(Tn[k], Axes(2, 3, 0, 1, 4)),
                                transpose(Tn[j], Axes(2, 3, 0, 1, 4)),
                                peps_parameters, PUs[iy], PLs[iy]);
    } else {
      Calc_projector_updown_blocks(
          C3[k], C4[l], C1[i], C2[j], eTb[k], eTb[l], eTl[l], eTl[i], eTt[i],
          eTt[j], eTr[j], eTr[k], transpose(Tn[k], Axes(2, 3, 0, 1, 4)),
          transpose(Tn[l], Axes(2, 3, 0, 1, 4)),
          transpose(Tn[i], Axes(2, 3, 0, 1, 4)),
          transpose(Tn[j], Axes(2, 3, 0, 1, 4)), peps_parameters, PUs[iy],
          PLs[iy]);
    }
  }
  // update
  std::vector<Tensor<Matrix, C>> C2_bak(lattice.N_UNIT), C3_bak(lattice.N_UNIT),
      eTr_bak(lattice.N_UNIT);
  for (int num = 0; num < lattice.N_UNIT; num++) {
    C2_bak[num] = C2[num];
    C3_bak[num] = C3[num];
    eTr_bak[num] = eTr[num];
  }
  int iy_up, iy_down;
  for (int iy = 0; iy < lattice.LY; ++iy) {
    k = lattice.index(ix, iy);
    l = lattice.left(k);
    i = lattice.top(l);
    j = lattice.right(i);

    iy_up = (iy + 1) % lattice.LY;
    iy_down = (iy - 1 + lattice.LY) % lattice.LY;

    Calc_Next_CTM(C3_bak[k], C2_bak[j], eTb[k], eTt[j], PUs[iy_down],
                  PLs[iy_up], C3[l], C2[i]);

    Calc_Next_eT(eTr_bak[k], transpose(Tn[k], Axes(2, 3, 0, 1, 4)), PUs[iy],
                 PLs[iy_down], eTr[l]);
    Calc_Next_eT(eTr_bak[j], transpose(Tn[j], Axes(2, 3, 0, 1, 4)), PUs[iy_up],
                 PLs[iy], eTr[i]);
  }
}

template <template <typename> class Matrix, typename C>
void Top_move(std::vector<Tensor<Matrix, C>> &C1,
              std::vector<Tensor<Matrix, C>> &C2,
              const std::vector<Tensor<Matrix, C>> &C3,
              const std::vector<Tensor<Matrix, C>> &C4,
              std::vector<Tensor<Matrix, C>> &eTt,
              const std::vector<Tensor<Matrix, C>> &eTr,
              const std::vector<Tensor<Matrix, C>> &eTb,
              const std::vector<Tensor<Matrix, C>> &eTl,
              const std::vector<Tensor<Matrix, C>> &Tn, const int iy,
              const PEPS_Parameters peps_parameters, const Lattice lattice) {
  /*
    ## Do one step top move absorbing Y=iy row
    ## part of C1, C2, eTt will be modified
  */
  std::vector<Tensor<Matrix, C>> PUs, PLs;
  PUs.resize(lattice.LX);
  PLs.resize(lattice.LX);
  int i, j, k, l;
  for (int ix = 0; ix < lattice.LX; ++ix) {
    j = lattice.index(ix, iy);
    k = lattice.bottom(j);
    l = lattice.left(k);
    i = lattice.top(l);

    if (peps_parameters.CTM_Projector_corner) {
      Calc_projector_left_block(C2[j], C1[i], eTr[j], eTl[i], eTt[i], eTt[j],
                                transpose(Tn[j], Axes(1, 2, 3, 0, 4)),
                                transpose(Tn[i], Axes(1, 2, 3, 0, 4)),
                                peps_parameters, PUs[ix], PLs[ix]);
    } else {
      Calc_projector_updown_blocks(
          C2[j], C3[k], C4[l], C1[i], eTr[j], eTr[k], eTb[k], eTb[l], eTl[l],
          eTl[i], eTt[i], eTt[j], transpose(Tn[j], Axes(1, 2, 3, 0, 4)),
          transpose(Tn[k], Axes(1, 2, 3, 0, 4)),
          transpose(Tn[l], Axes(1, 2, 3, 0, 4)),
          transpose(Tn[i], Axes(1, 2, 3, 0, 4)), peps_parameters, PUs[ix],
          PLs[ix]);
    }
  }
  // update
  std::vector<Tensor<Matrix, C>> C1_bak(lattice.N_UNIT), C2_bak(lattice.N_UNIT),
      eTt_bak(lattice.N_UNIT);
  for (int num = 0; num < lattice.N_UNIT; num++) {
    C1_bak[num] = C1[num];
    C2_bak[num] = C2[num];
    eTt_bak[num] = eTt[num];
  }
  int ix_right, ix_left;
  for (int ix = 0; ix < lattice.LX; ++ix) {
    j = lattice.index(ix, iy);
    k = lattice.bottom(j);
    l = lattice.left(k);
    i = lattice.top(l);

    ix_right = (ix + 1) % lattice.LX;
    ix_left = (ix - 1 + lattice.LX) % lattice.LX;

    Calc_Next_CTM(C2_bak[j], C1_bak[i], eTr[j], eTl[i], PUs[ix_right],
                  PLs[ix_left], C2[k], C1[l]);

    Calc_Next_eT(eTt_bak[j], transpose(Tn[j], Axes(1, 2, 3, 0, 4)), PUs[ix],
                 PLs[ix_right], eTt[k]);
    Calc_Next_eT(eTt_bak[i], transpose(Tn[i], Axes(1, 2, 3, 0, 4)),
                 PUs[ix_left], PLs[ix], eTt[l]);
  }
}
template <template <typename> class Matrix, typename C>
void Bottom_move(const std::vector<Tensor<Matrix, C>> &C1,
                 const std::vector<Tensor<Matrix, C>> &C2,
                 std::vector<Tensor<Matrix, C>> &C3,
                 std::vector<Tensor<Matrix, C>> &C4,
                 const std::vector<Tensor<Matrix, C>> &eTt,
                 const std::vector<Tensor<Matrix, C>> &eTr,
                 std::vector<Tensor<Matrix, C>> &eTb,
                 const std::vector<Tensor<Matrix, C>> &eTl,
                 const std::vector<Tensor<Matrix, C>> &Tn, const int iy,
                 const PEPS_Parameters peps_parameters, const Lattice lattice) {
  /*
    ## Do one step bottom move absorbing Y=iy row
    ## part of C3, C4, eTb will be modified
  */

  std::vector<Tensor<Matrix, C>> PUs, PLs;
  PUs.resize(lattice.LX);
  PLs.resize(lattice.LX);
  int i, j, k, l;
  for (int ix = 0; ix < lattice.LX; ++ix) {
    l = lattice.index(ix, iy);
    i = lattice.top(l);
    j = lattice.right(i);
    k = lattice.bottom(j);

    if (peps_parameters.CTM_Projector_corner) {
      Calc_projector_left_block(C4[l], C3[k], eTl[l], eTr[k], eTb[k], eTb[l],
                                transpose(Tn[l], Axes(3, 0, 1, 2, 4)),
                                transpose(Tn[k], Axes(3, 0, 1, 2, 4)),
                                peps_parameters, PUs[ix], PLs[ix]);
    } else {
      Calc_projector_updown_blocks(
          C4[l], C1[i], C2[j], C3[k], eTl[l], eTl[i], eTt[i], eTt[j], eTr[j],
          eTr[k], eTb[k], eTb[l], transpose(Tn[l], Axes(3, 0, 1, 2, 4)),
          transpose(Tn[i], Axes(3, 0, 1, 2, 4)),
          transpose(Tn[j], Axes(3, 0, 1, 2, 4)),
          transpose(Tn[k], Axes(3, 0, 1, 2, 4)), peps_parameters, PUs[ix],
          PLs[ix]);
    }
  }

  // update
  std::vector<Tensor<Matrix, C>> C3_bak(lattice.N_UNIT), C4_bak(lattice.N_UNIT),
      eTb_bak(lattice.N_UNIT);
  for (int num = 0; num < lattice.N_UNIT; num++) {
    C3_bak[num] = C3[num];
    C4_bak[num] = C4[num];
    eTb_bak[num] = eTb[num];
  }
  int ix_left, ix_right;
  for (int ix = 0; ix < lattice.LX; ++ix) {
    l = lattice.index(ix, iy);
    i = lattice.top(l);
    j = lattice.right(i);
    k = lattice.bottom(j);

    ix_right = (ix + 1) % lattice.LX;
    ix_left = (ix - 1 + lattice.LX) % lattice.LX;

    Calc_Next_CTM(C4_bak[l], C3_bak[k], eTl[l], eTr[k], PUs[ix_left],
                  PLs[ix_right], C4[i], C3[j]);

    Calc_Next_eT(eTb_bak[l], transpose(Tn[l], Axes(3, 0, 1, 2, 4)), PUs[ix],
                 PLs[ix_left], eTb[i]);
    Calc_Next_eT(eTb_bak[k], transpose(Tn[k], Axes(3, 0, 1, 2, 4)),
                 PUs[ix_right], PLs[ix], eTb[j]);
  }
}

template <template <typename> class Matrix, typename C>
bool Check_Convergence_CTM(const std::vector<Tensor<Matrix, C>> &C1,
                           const std::vector<Tensor<Matrix, C>> &C2,
                           const std::vector<Tensor<Matrix, C>> &C3,
                           const std::vector<Tensor<Matrix, C>> &C4,
                           const std::vector<Tensor<Matrix, C>> &C1_old,
                           const std::vector<Tensor<Matrix, C>> &C2_old,
                           const std::vector<Tensor<Matrix, C>> &C3_old,
                           const std::vector<Tensor<Matrix, C>> &C4_old,
                           const PEPS_Parameters peps_parameters,
                           const Lattice lattice, double &sig_max) {
  sig_max = 0.0;
  bool convergence = true;
  double sig, norm;
  std::vector<double> lam_new, lam_old;

  for (int i = 0; i < lattice.N_UNIT; ++i) {
    // C1
    svd(C1[i], lam_new);
    norm = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      norm += lam_new[j] * lam_new[j];
    }
    norm = sqrt(norm);
    for (int k = 0; k < lam_new.size(); ++k) {
      lam_new[k] /= norm;
    }
    svd(C1_old[i], lam_old);
    norm = 0.0;
    for (int j = 0; j < lam_old.size(); ++j) {
      norm += lam_old[j] * lam_old[j];
    }
    norm = sqrt(norm);
    for (int k = 0; k < lam_old.size(); ++k) {
      lam_old[k] /= norm;
    }

    sig = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      sig += (lam_new[j] - lam_old[j]) * (lam_new[j] - lam_old[j]);
    }
    sig = sqrt(sig);

    if (sig > peps_parameters.CTM_Convergence_Epsilon) {
      sig_max = sig;
      convergence = false;
      break;
    } else if (sig > sig_max) {
      sig_max = sig;
    }

    // C2
    svd(C2[i], lam_new);
    norm = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      norm += lam_new[j] * lam_new[j];
    }
    norm = sqrt(norm);
    for (int k = 0; k < lam_new.size(); ++k) {
      lam_new[k] /= norm;
    }

    svd(C2_old[i], lam_old);
    norm = 0.0;
    for (int j = 0; j < lam_old.size(); ++j) {
      norm += lam_old[j] * lam_old[j];
    }

    norm = sqrt(norm);
    for (int k = 0; k < lam_old.size(); ++k) {
      lam_old[k] /= norm;
    }

    sig = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      sig += (lam_new[j] - lam_old[j]) * (lam_new[j] - lam_old[j]);
    }
    sig = sqrt(sig);

    if (sig > peps_parameters.CTM_Convergence_Epsilon) {
      sig_max = sig;
      convergence = false;
      break;
    } else if (sig > sig_max) {
      sig_max = sig;
    }

    // C3
    svd(C3[i], lam_new);
    norm = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      norm += lam_new[j] * lam_new[j];
    }
    norm = sqrt(norm);
    for (int k = 0; k < lam_new.size(); ++k) {
      lam_new[k] /= norm;
    }
    svd(C3_old[i], lam_old);
    norm = 0.0;
    for (int j = 0; j < lam_old.size(); ++j) {
      norm += lam_old[j] * lam_old[j];
    }

    norm = sqrt(norm);
    for (int k = 0; k < lam_old.size(); ++k) {
      lam_old[k] /= norm;
    }

    sig = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      sig += (lam_new[j] - lam_old[j]) * (lam_new[j] - lam_old[j]);
    }
    sig = sqrt(sig);

    if (sig > peps_parameters.CTM_Convergence_Epsilon) {
      sig_max = sig;
      convergence = false;
      break;
    } else if (sig > sig_max) {
      sig_max = sig;
    }

    // C4
    svd(C4[i], lam_new);
    norm = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      norm += lam_new[j] * lam_new[j];
    }
    norm = sqrt(norm);
    for (int k = 0; k < lam_new.size(); ++k) {
      lam_new[k] /= norm;
    }
    svd(C4_old[i], lam_old);
    norm = 0.0;
    for (int j = 0; j < lam_old.size(); ++j) {
      norm += lam_old[j] * lam_old[j];
    }
    norm = sqrt(norm);
    for (int k = 0; k < lam_old.size(); ++k) {
      lam_old[k] /= norm;
    }

    sig = 0.0;
    for (int j = 0; j < lam_new.size(); ++j) {
      sig += (lam_new[j] - lam_old[j]) * (lam_new[j] - lam_old[j]);
    }
    sig = sqrt(sig);

    if (sig > peps_parameters.CTM_Convergence_Epsilon) {
      sig_max = sig;
      convergence = false;
      break;
    } else if (sig > sig_max) {
      sig_max = sig;
    }
  }
  return convergence;
}

template <template <typename> class Matrix, typename C>
int Calc_CTM_Environment(
    std::vector<Tensor<Matrix, C>> &C1, std::vector<Tensor<Matrix, C>> &C2,
    std::vector<Tensor<Matrix, C>> &C3, std::vector<Tensor<Matrix, C>> &C4,
    std::vector<Tensor<Matrix, C>> &eTt, std::vector<Tensor<Matrix, C>> &eTr,
    std::vector<Tensor<Matrix, C>> &eTb, std::vector<Tensor<Matrix, C>> &eTl,
    const std::vector<Tensor<Matrix, C>> &Tn,
    const PEPS_Parameters peps_parameters, const Lattice lattice,
    bool initialize = true) {
  /*
    ## Calc environment tensors
    ## C1,C2,C3,C4 and eTt,eTl,eTr,eTb will be modified
  */
  // Initialize
  if (initialize) {
    int num, d1, d2, d34;
    Index index;
    Tensor<Matrix, C> Projector;

    for (int i = 0; i < lattice.N_UNIT; ++i) {
      num = lattice.top(lattice.left(i));
      d1 = Tn[num].shape()[3] * Tn[num].shape()[3];
      d2 = Tn[num].shape()[2] * Tn[num].shape()[2];
      C1[i] = reshape(
          tensordot(Tn[num], conj(Tn[num]), Axes(0, 1, 4), Axes(0, 1, 4))
              .transpose(Axes(1, 3, 0, 2)),
          (Shape(d1, d2)));

      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        C1[i] = extend(C1[i], Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        C1[i] = slice(slice(C1[i], 0, 0, peps_parameters.CHI), 1, 0,
                      peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        C1[i] = extend(slice(C1[i], 1, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else {
        // d1 >= CHI
        // d2 < CHI
        C1[i] = extend(slice(C1[i], 0, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      }

      num = lattice.top(lattice.right(i));
      d1 = Tn[num].shape()[0] * Tn[num].shape()[0];
      d2 = Tn[num].shape()[3] * Tn[num].shape()[3];
      C2[i] = reshape(
          tensordot(Tn[num], conj(Tn[num]), Axes(1, 2, 4), Axes(1, 2, 4))
              .transpose(Axes(0, 2, 1, 3)),
          (Shape(d1, d2)));

      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        C2[i] = extend(C2[i], Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        C2[i] = slice(slice(C2[i], 0, 0, peps_parameters.CHI), 1, 0,
                      peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        C2[i] = extend(slice(C2[i], 1, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else {
        // d1 >= CHI
        // d2 < CHI
        C2[i] = extend(slice(C2[i], 0, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      }

      num = lattice.bottom(lattice.right(i));
      d1 = Tn[num].shape()[1] * Tn[num].shape()[1];
      d2 = Tn[num].shape()[0] * Tn[num].shape()[0];
      C3[i] = reshape(
          tensordot(Tn[num], conj(Tn[num]), Axes(2, 3, 4), Axes(2, 3, 4))
              .transpose(Axes(1, 3, 0, 2)),
          (Shape(d1, d2)));

      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        C3[i] = extend(C3[i], Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        C3[i] = slice(slice(C3[i], 0, 0, peps_parameters.CHI), 1, 0,
                      peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        C3[i] = extend(slice(C3[i], 1, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else {
        // d1 >= CHI
        // d2 < CHI
        C3[i] = extend(slice(C3[i], 0, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      }

      num = lattice.bottom(lattice.left(i));
      d1 = Tn[num].shape()[2] * Tn[num].shape()[2];
      d2 = Tn[num].shape()[1] * Tn[num].shape()[1];
      C4[i] = reshape(
          tensordot(Tn[num], conj(Tn[num]), Axes(0, 3, 4), Axes(0, 3, 4))
              .transpose(Axes(1, 3, 0, 2)),
          (Shape(d1, d2)));

      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        C4[i] = extend(C4[i], Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        C4[i] = slice(slice(C4[i], 0, 0, peps_parameters.CHI), 1, 0,
                      peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        C4[i] = extend(slice(C4[i], 1, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      } else {
        // d1 >= CHI
        // d2 < CHI
        C4[i] = extend(slice(C4[i], 0, 0, peps_parameters.CHI),
                       Shape(peps_parameters.CHI, peps_parameters.CHI));
      }

      num = lattice.top(i);
      d1 = Tn[num].shape()[0] * Tn[num].shape()[0];
      d2 = Tn[num].shape()[2] * Tn[num].shape()[2];
      d34 = Tn[num].shape()[3];
      eTt[i] = reshape(tensordot(Tn[num], conj(Tn[num]), Axes(1, 4), Axes(1, 4))
                           .transpose(Axes(0, 3, 1, 4, 2, 5)),
                       (Shape(d1, d2, d34, d34)));

      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        eTt[i] = extend(
            eTt[i], Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        eTt[i] = slice(slice(eTt[i], 0, 0, peps_parameters.CHI), 1, 0,
                       peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        eTt[i] =
            extend(slice(eTt[i], 1, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else {
        // d1 >= CHI
        // d2 < CHI
        eTt[i] =
            extend(slice(eTt[i], 0, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      }

      num = lattice.right(i);
      d1 = Tn[num].shape()[1] * Tn[num].shape()[1];
      d2 = Tn[num].shape()[3] * Tn[num].shape()[3];
      d34 = Tn[num].shape()[0];
      eTr[i] = reshape(tensordot(Tn[num], conj(Tn[num]), Axes(2, 4), Axes(2, 4))
                           .transpose(Axes(1, 4, 2, 5, 0, 3)),
                       (Shape(d1, d2, d34, d34)));

      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        eTr[i] = extend(
            eTr[i], Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        eTr[i] = slice(slice(eTr[i], 0, 0, peps_parameters.CHI), 1, 0,
                       peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        eTr[i] =
            extend(slice(eTr[i], 1, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else {
        // d1 >= CHI
        // d2 < CHI
        eTr[i] =
            extend(slice(eTr[i], 0, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      }

      num = lattice.bottom(i);
      d1 = Tn[num].shape()[2] * Tn[num].shape()[2];
      d2 = Tn[num].shape()[0] * Tn[num].shape()[0];
      d34 = Tn[num].shape()[1];
      eTb[i] = reshape(tensordot(Tn[num], conj(Tn[num]), Axes(3, 4), Axes(3, 4))
                           .transpose(Axes(2, 5, 0, 3, 1, 4)),
                       (Shape(d1, d2, d34, d34)));

      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        eTb[i] = extend(
            eTb[i], Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        eTb[i] = slice(slice(eTb[i], 0, 0, peps_parameters.CHI), 1, 0,
                       peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        eTb[i] =
            extend(slice(eTb[i], 1, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else {
        // d1 >= CHI
        // d2 < CHI
        eTb[i] =
            extend(slice(eTb[i], 0, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      }

      num = lattice.left(i);
      d1 = Tn[num].shape()[3] * Tn[num].shape()[3];
      d2 = Tn[num].shape()[1] * Tn[num].shape()[1];
      d34 = Tn[num].shape()[2];
      eTl[i] = reshape(tensordot(Tn[num], conj(Tn[num]), Axes(0, 4), Axes(0, 4))
                           .transpose(Axes(2, 5, 0, 3, 1, 4)),
                       (Shape(d1, d2, d34, d34)));
      if (d1 < peps_parameters.CHI && d2 < peps_parameters.CHI) {
        eTl[i] = extend(
            eTl[i], Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else if (d1 >= peps_parameters.CHI && d2 >= peps_parameters.CHI) {
        eTl[i] = slice(slice(eTl[i], 0, 0, peps_parameters.CHI), 1, 0,
                       peps_parameters.CHI);
      } else if (d1 < peps_parameters.CHI) {
        // d2 >= CHI
        eTl[i] =
            extend(slice(eTl[i], 1, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      } else {
        // d1 >= CHI
        // d2 < CHI
        eTl[i] =
            extend(slice(eTl[i], 0, 0, peps_parameters.CHI),
                   Shape(peps_parameters.CHI, peps_parameters.CHI, d34, d34));
      }
    }
  }
  // Initialize done

  bool convergence = false;
  int count = 0;
  std::vector<Tensor<Matrix, C>> C1_old = C1;
  std::vector<Tensor<Matrix, C>> C2_old = C2;
  std::vector<Tensor<Matrix, C>> C3_old = C3;
  std::vector<Tensor<Matrix, C>> C4_old = C4;

  double sig_max = 0.0;
  while ((!convergence) && (count < peps_parameters.Max_CTM_Iteration)) {
    // left move
    for (int ix = 0; ix < lattice.LX; ++ix) {
      Left_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, ix, peps_parameters,
                lattice);
    }

    // right move
    for (int ix = 0; ix > -lattice.LX; --ix) {
      Right_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
                 (ix + 1 + lattice.LX) % lattice.LX, peps_parameters, lattice);
    }

    // top move
    for (int iy = 0; iy > -lattice.LY; --iy) {
      Top_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
               (iy + 1 + lattice.LY) % lattice.LY, peps_parameters, lattice);
    }

    // bottom move

    for (int iy = 0; iy < lattice.LY; ++iy) {
      Bottom_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, iy, peps_parameters,
                  lattice);
    }

    convergence =
        Check_Convergence_CTM(C1, C2, C3, C4, C1_old, C2_old, C3_old, C4_old,
                              peps_parameters, lattice, sig_max);
    count += 1;

    C1_old = C1;
    C2_old = C2;
    C3_old = C3;
    C4_old = C4;
    if (peps_parameters.print_level >= PrintLevel::debug) {
      std::cout << "CTM: count, sig_max " << count << " " << sig_max
                << std::endl;
    }
  }

  if (!convergence && peps_parameters.print_level >= PrintLevel::warn) {
    std::cout << "Warning: CTM did not converge! count, sig_max = " << count
              << " " << sig_max << std::endl;
  }
  if (peps_parameters.print_level >= PrintLevel::debug) {
    std::cout << "CTM: count to convergence= " << count << std::endl;
  }
  return count;
}

}  // end of namespace tenes

#endif  // TENES_SRC_SQUARE_LATTICE_HPP_
