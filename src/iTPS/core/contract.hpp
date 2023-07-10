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

/*
 *
 Basic routines independent on unit cell structures.
 Using mptensor libraries
 (Test version)
 2015 Dec.  Tsuyoshi Okubo
*/

#ifndef TENES_SRC_ITPS_CORE_CONTRACT_HPP_
#define TENES_SRC_ITPS_CORE_CONTRACT_HPP_

#include <vector>
#include <cstddef>

#include "contract_itps_ctm.hpp"
#include "contract_itps_mf.hpp"
#include "contract_density_ctm.hpp"

namespace tenes {
namespace itps {
namespace core {

/*! @brief contract tensors with CTM
 *
 *  @param[in] C corner transfer matrix
 *  @param[in] eTt top edge tensors
 *  @param[in] eTr right edge tensors
 *  @param[in] eTb bottom edge tensors
 *  @param[in] eTl left edge tensors
 *  @param[in] Tn center tensors
 *  @param[in] op onesite operators
 */
template <class tensor>
typename tensor::value_type Contract(
    const std::vector<const tensor *> &C,
    const std::vector<const tensor *> &eTt,
    const std::vector<const tensor *> &eTr,
    const std::vector<const tensor *> &eTb,
    const std::vector<const tensor *> &eTl,
    const std::vector<std::vector<const tensor *>> &Tn,
    const std::vector<std::vector<const tensor *>> &op, bool is_density,
    bool is_meanfield) {
  if (is_density) {
    if (is_meanfield) {
      throw std::runtime_error("Not implemented");
      // return Contract_density_MF(Tn, op);
    } else {
      return Contract_density_CTM(C, eTt, eTr, eTb, eTl, Tn, op);
    }
  } else {
    if (is_meanfield) {
      return Contract_iTPS_MF(Tn, op);
    } else {
      return Contract_iTPS_CTM(C, eTt, eTr, eTb, eTl, Tn, op);
    }
  }
}

template <class tensor>
typename tensor::value_type Contract_one_site(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &Tn1, const tensor &op1, bool is_density, bool is_meanfield) {
  if (is_density) {
    if (is_meanfield) {
      throw std::runtime_error("Not implemented");
      // return Contract_one_site_density_MF(Tn1, op1);
    } else {
      return Contract_one_site_density_CTM(C1, C2, C3, C4, eT1, eT2, eT3, eT4,
                                           Tn1, op1);
    }
  } else {
    if (is_meanfield) {
      return Contract_one_site_iTPS_MF(Tn1, op1);
    } else {
      return Contract_one_site_iTPS_CTM(C1, C2, C3, C4, eT1, eT2, eT3, eT4, Tn1,
                                   op1);
    }
  }
}

template <class tensor>
typename tensor::value_type Contract_two_sites_horizontal(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op1, const tensor &op2, bool is_density, bool is_meanfield) {
  if (is_density) {
    if (is_meanfield) {
      throw std::runtime_error("Not implemented");
    } else {
      return Contract_two_site_horizontal_density_CTM(
          C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6, Tn1, Tn2, op1, op2);
    }
  } else {
    if (is_meanfield) {
      return Contract_two_sites_horizontal_iTPS_MF(Tn1, op1);
    } else {
      return Contract_two_sites_horizontal_iTPS_CTM(
          C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6, Tn1, Tn2, op1, op2);
    }
  }
}

template <class tensor>
typename tensor::value_type Contract_two_sites_vertical(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op1, const tensor &op2, bool is_density, bool is_meanfield) {
  if (is_density) {
    if (is_meanfield) {
      throw std::runtime_error("Not implemented");
    } else {
      return Contract_two_site_vertical_density_CTM(
          C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6, Tn1, Tn2, op1, op2);
    }
  } else {
    if (is_meanfield) {
      return Contract_two_sites_vertical_iTPS_MF(Tn1, op1);
    } else {
      return Contract_two_sites_vertical_iTPS_CTM(C1, C2, C3, C4, eT1, eT2, eT3, eT4,
                                             eT5, eT6, Tn1, Tn2, op1, op2);
    }
  }
}

template <class tensor>
typename tensor::value_type Contract_two_sites_horizontal_op12(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op12, bool is_density, bool is_meanfield) {
  if (is_density) {
    if (is_meanfield) {
      throw std::runtime_error("Not implemented");
    } else {
      return Contract_two_sites_horizontal_op12_density_CTM(
          C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6, Tn1, Tn2, op12);
    }
  } else {
    if (is_meanfield) {
      return Contract_two_sites_horizontal_op12_iTPS_MF(Tn1, Tn2, op12);
    } else {
      return Contract_two_sites_horizontal_op12_iTPS_CTM(
          C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6, Tn1, Tn2, op12);
    }
  }
}

template <class tensor>
typename tensor::value_type Contract_two_sites_vertical_op12(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &Tn1, const tensor &Tn2,
    const tensor &op12, bool is_density, bool is_meanfield) {
  if (is_density) {
    if (is_meanfield) {
      throw std::runtime_error("Not implemented");
    } else {
      return Contract_two_sites_vertical_op12_density_CTM(
          C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6, Tn1, Tn2, op12);
    }
  } else {
    if (is_meanfield) {
      return Contract_two_sites_vertical_op12_iTPS_MF(Tn1, Tn2, op12);
    } else {
      return Contract_two_sites_vertical_op12_iTPS_CTM(
          C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5, eT6, Tn1, Tn2, op12);
    }
  }
}

template <class tensor>
typename tensor::value_type Contract_four_sites(
    const tensor &C1, const tensor &C2, const tensor &C3, const tensor &C4,
    const tensor &eT1, const tensor &eT2, const tensor &eT3, const tensor &eT4,
    const tensor &eT5, const tensor &eT6, const tensor &eT7, const tensor &eT8,
    const tensor &Tn1, const tensor &Tn2, const tensor &Tn3, const tensor &Tn4,
    const tensor &op1, const tensor &op2, const tensor &op3, const tensor &op4,
    bool is_density, bool is_meanfield) {
  if (is_density) {
    if (is_meanfield) {
      throw std::runtime_error("Not implemented");
    } else {
      return Contract_four_sites_density_CTM(C1, C2, C3, C4, eT1, eT2, eT3, eT4,
                                             eT5, eT6, eT7, eT8, Tn1, Tn2, Tn3,
                                             Tn4, op1, op2, op3, op4);
    }
  } else {
    if (is_meanfield) {
      return Contract_four_sites_iTPS_MF(Tn1, Tn2, Tn3, Tn4, op1, op2, op3, op4);
    } else {
      return Contract_four_sites_iTPS_CTM(C1, C2, C3, C4, eT1, eT2, eT3, eT4, eT5,
                                     eT6, eT7, eT8, Tn1, Tn2, Tn3, Tn4, op1,
                                     op2, op3, op4);
    }
  }
}

template <class tensor>
void StartCorrelation(tensor &A, const tensor &C1, const tensor &C4,
                      const tensor &eT1, const tensor &eT3, const tensor &eT4,
                      const tensor &Tn1, const tensor &op, bool is_density,
                      bool is_meanfield) {
  if (is_density) {
    throw std::runtime_error("Not implemented");
    // if (is_meanfield) {
    //   throw std::runtime_error("Not implemented");
    // } else {
    //   StartCorrelation_density_CTM(A, C1, C4, eT1, eT3, eT4, Tn1, op);
    // }
  } else {
    if (is_meanfield) {
      StartCorrelation_iTPS_MF(A, Tn1, op);
    } else {
      StartCorrelation_iTPS_CTM(A, C1, C4, eT1, eT3, eT4, Tn1, op);
    }
  }
}

template <class tensor>
void Transfer(tensor &A, const tensor &eT1, const tensor &eT3,
              const tensor &Tn1, bool is_density, bool is_meanfield) {
  if (is_density) {
    throw std::runtime_error("Not implemented");
    // if (is_meanfield) {
    //   throw std::runtime_error("Not implemented");
    // } else {
    //   Transfer_density_CTM(A, eT1, eT3, Tn1);
    // }
  } else {
    if (is_meanfield) {
      Transfer_iTPS_MF(A, Tn1);
    } else {
      Transfer_iTPS_CTM(A, eT1, eT3, Tn1);
    }
  }
};

template <class tensor>
typename tensor::value_type FinishCorrelation(
    const tensor &A, const tensor &C2, const tensor &C3, const tensor &eT1,
    const tensor &eT2, const tensor &eT3, const tensor &Tn1, const tensor &op,
    bool is_density, bool is_meanfield) {
  if (is_density) {
    throw std::runtime_error("Not implemented");
    // if (is_meanfield) {
    //   throw std::runtime_error("Not implemented");
    // } else {
    //   return FinishCorrelation_density_CTM(A, C2, C3, eT1, eT2, eT3, Tn1,
    //   op);
    // }
  } else {
    if (is_meanfield) {
      return FinishCorrelation_iTPS_MF(A, Tn1, op);
    } else {
      return FinishCorrelation_iTPS_CTM(A, C2, C3, eT1, eT2, eT3, Tn1, op);
    }
  }
};

template <class tensor>
void TransferMatrix_MatVec(tensor &inoutvec, const tensor &eT1, bool is_density,
                           bool is_meanfield) {
  if (is_density) {
    throw std::runtime_error("Not implemented");
    // if (is_meanfield) {
    //   throw std::runtime_error("Not implemented");
    // } else {
    //   TransferMatrix_MatVec_density_CTM(inoutvec, eT1);
    // }
  } else {
    if (is_meanfield) {
      TransferMatrix_MatVec_iTPS_MF(inoutvec, eT1);
    } else {
      TransferMatrix_MatVec_iTPS_CTM(inoutvec, eT1);
    }
  }
};

}  // end of namespace core
}  // namespace itps
}  // namespace tenes

#endif  // TENES_SRC_ITPS_CORE_CONTRACT_CTM_HPP_
