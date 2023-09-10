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

#include "../optimize_contraction_order/optimize.hpp"

#include <stack>

#include "iTPS.hpp"
#include "contraction_order.hpp"

namespace {
inline int coord2index(int row, int col, int nrows, int ncols) {
  return row * ncols + col;
}
inline std::pair<int, int> index2coord(int index, int nrows, int ncols) {
  return std::make_pair(index / ncols, index % ncols);
}

template <class tensor>
struct Tensor_w_bondname {
  tensor t;
  tensor const* t_ptr;
  // std::vector<std::string> bonds;
  std::vector<int> bonds;
  Tensor_w_bondname(tensor const* t_ptr_, std::vector<int> const& b_)
      : t_ptr(t_ptr_), bonds(b_) {}
  Tensor_w_bondname(tensor const& t_, std::vector<int> const& b_)
      : t(t_), t_ptr(nullptr), bonds(b_) {}
};

template <class tensor>
Tensor_w_bondname<tensor> contract(Tensor_w_bondname<tensor> const& lhs,
                                   Tensor_w_bondname<tensor> const& rhs) {
  mptensor::Axes axes_lhs, axes_rhs;
  const size_t nb_lhs = lhs.bonds.size();
  const size_t nb_rhs = rhs.bonds.size();

  std::vector<int> bonds_lhs = lhs.bonds;
  std::vector<int> bonds_rhs = rhs.bonds;

  size_t nb_common = 0;
  for (size_t i = 0; i < nb_lhs; ++i) {
    const int bond = lhs.bonds[i];
    for (size_t j = 0; j < nb_rhs; ++j) {
      if (bond == rhs.bonds[j]) {
        axes_lhs.push(i);
        axes_rhs.push(j);
        bonds_lhs[i] = -1;
        bonds_rhs[j] = -1;
        ++nb_common;
        break;
      }
    }
  }
  std::vector<int> bonds;
  bonds.reserve(nb_lhs + nb_rhs - 2 * nb_common);
  for (size_t i = 0; i < nb_lhs; ++i) {
    if (bonds_lhs[i] != -1) {
      bonds.push_back(bonds_lhs[i]);
    }
  }
  for (size_t i = 0; i < nb_rhs; ++i) {
    if (bonds_rhs[i] != -1) {
      bonds.push_back(bonds_rhs[i]);
    }
  }

  tensor const* t_ptr_lhs = lhs.t_ptr;
  tensor const* t_ptr_rhs = rhs.t_ptr;
  if (t_ptr_lhs == nullptr) {
    t_ptr_lhs = &lhs.t;
  }
  if (t_ptr_rhs == nullptr) {
    t_ptr_rhs = &rhs.t;
  }

  return Tensor_w_bondname<tensor>(
      mptensor::tensordot(*t_ptr_lhs, *t_ptr_rhs, axes_lhs, axes_rhs), bonds);
}

}  // namespace

namespace tenes {
namespace itps {

template <class tensor>
TensorNetworkContractor<tensor>::TensorNetworkContractor(int nrows, int ncols,
                                                         bool is_tpo,
                                                         bool is_mf) {
  initialize(nrows, ncols, is_tpo, is_mf);
}

template <class tensor>
void TensorNetworkContractor<tensor>::initialize(int nrows, int ncols,
                                                 bool is_tpo, bool is_mf) {
  this->nrows = nrows;
  this->ncols = ncols;
  this->is_tpo = is_tpo;
  this->is_mf = is_mf;
  order.clear();

  if (is_mf) {
    initialize_mf();
  } else {
    initialize_ctm();
  }
}

template <class tensor>
std::pair<std::vector<int>, double> TensorNetworkContractor<tensor>::optimize(
    std::vector<int> const& shape_types,
    std::vector<std::vector<int>> const& shape_dims, int chi) {
  const int tn_num_bonds = (is_tpo ? 4 : 5);

  optimize_contraction_order::TensorNetwork tnw;
  tnw.set_default_bond_dim(chi);

  // add tensor
  for (const auto& t : C_name) {
    tnw.add_tensor(t.name, t.bonds);
  }
  for (const auto& t : eTt_name) {
    tnw.add_tensor(t.name, t.bonds);
  }
  for (const auto& t : eTr_name) {
    tnw.add_tensor(t.name, t.bonds);
  }
  for (const auto& t : eTb_name) {
    tnw.add_tensor(t.name, t.bonds);
  }
  for (const auto& t : eTl_name) {
    tnw.add_tensor(t.name, t.bonds);
  }
  for (const auto& t : Tn_name) {
    tnw.add_tensor(t.name, t.bonds);
  }

  const int N_UNIT = nrows * ncols;
  for (int i = 0; i < N_UNIT; ++i) {
    int type = shape_types[i];
    std::vector<int> const& dims = shape_dims[type];
    if (is_tpo) {
      for (int j = 0; j < tn_num_bonds; ++j) {
        tnw.set_bond_dim(Tn_name[i].bonds[j], dims[j]);
      }
    } else {
      for (int j = 0; j < tn_num_bonds; ++j) {
        tnw.set_bond_dim(Tn_name[i].bonds[j], dims[j]);
        tnw.set_bond_dim(Tn_name[i + N_UNIT].bonds[j], dims[j]);
      }
    }
  }

  auto res = tnw.optimize();
  order = std::get<0>(res);
  assert(check_order(order));
  return res;
}

template <class tensor>
bool TensorNetworkContractor<tensor>::check_order(
    std::vector<int> const& order) const {
  size_t num_tensors = 0;
  num_tensors += C_name.size();
  num_tensors += eTt_name.size();
  num_tensors += eTr_name.size();
  num_tensors += eTb_name.size();
  num_tensors += eTl_name.size();
  num_tensors += Tn_name.size();

  std::vector<int> t(num_tensors, 0);

  size_t num_contract = 0;
  for (auto p : order) {
    if (p >= 0) {
      ++t[p];
    } else {
      ++num_contract;
    }
  }
  bool res = true;
  if (num_contract != num_tensors - 1) {
    std::cerr << "ERROR: number of -1 is " << num_contract
              << " (expected: " << num_tensors - 1 << ")" << std::endl;
    res = false;
  }
  for (size_t i = 0; i < num_tensors; ++i) {
    int n = t[i];
    if (n != 1) {
      std::cerr << "ERROR: number of " << i << " is " << n << " (expected: 1)"
                << std::endl;
      res = false;
    }
  }
  std::flush(std::cerr);
  return res;
}

template <class tensor>
void TensorNetworkContractor<tensor>::set_order(std::vector<int> const& order) {
  if (!check_order(order)) {
    throw std::runtime_error("invalid contraction order");
  }
  this->order = order;
}

template <class tensor>
void TensorNetworkContractor<tensor>::initialize_mf() {
  if (is_tpo) {
    throw std::runtime_error("not implemented");
  }

  const int tn_num_bonds = 5;
  const int N_UNIT = ncols * nrows;
  const int num_tn = 2 * N_UNIT;
  Tn_name.resize(num_tn, TensorName(tn_num_bonds));

  int bond_index = 0;

  for (int i = 0; i < N_UNIT; ++i) {
    Tn_name[i].name = "Tn" + std::to_string(i);
    Tn_name[i + N_UNIT].name = "cTn" + std::to_string(i);
    for (int j = 0; j < tn_num_bonds; ++j) {
      Tn_name[i].bonds[j] = bond_index++;
      Tn_name[i + N_UNIT].bonds[j] = bond_index++;
      // Tn_name[i].bonds[j] = Tn_name[i].name + "_" + std::to_string(j);
      // Tn_name[i + N_UNIT].bonds[j] =
      //     Tn_name[i + N_UNIT].name + "_" + std::to_string(j);
    }
    Tn_name[i + N_UNIT].bonds[4] = Tn_name[i].bonds[4];
  }

  // connect bonds of Tn
  for (int r = 0; r < nrows; ++r) {
    for (int c = 0; c < ncols; ++c) {
      int index = coord2index(r, c, nrows, ncols);

      if (r == 0) {
        Tn_name[index + N_UNIT].bonds[1] = Tn_name[index].bonds[1];
      } else {
        int index_t = coord2index(r - 1, c, nrows, ncols);
        Tn_name[index].bonds[1] = Tn_name[index_t].bonds[3];
        Tn_name[index + N_UNIT].bonds[1] = Tn_name[index_t + N_UNIT].bonds[3];
      }
      if (r == nrows - 1) {
        Tn_name[index + N_UNIT].bonds[3] = Tn_name[index].bonds[3];
      }

      if (c == 0) {
        Tn_name[index + N_UNIT].bonds[0] = Tn_name[index].bonds[0];
      } else {
        int index_l = coord2index(r, c - 1, nrows, ncols);
        Tn_name[index].bonds[0] = Tn_name[index_l].bonds[2];
        Tn_name[index + N_UNIT].bonds[0] = Tn_name[index_l + N_UNIT].bonds[2];
      }
      if (c == ncols - 1) {
        Tn_name[index + N_UNIT].bonds[2] = Tn_name[index].bonds[2];
      }
    }
  }
}

template <class tensor>
void TensorNetworkContractor<tensor>::initialize_ctm() {
  const int edge_size = (is_tpo ? 3 : 4);
  const int tn_num_bonds = (is_tpo ? 4 : 5);
  const int N_UNIT = ncols * nrows;
  const int num_tn = (is_tpo ? nrows * ncols : nrows * ncols * 2);

  int bond_index = 0;

  C_name.resize(4);
  eTt_name.resize(ncols, TensorName(edge_size));
  eTr_name.resize(nrows, TensorName(edge_size));
  eTb_name.resize(ncols, TensorName(edge_size));
  eTl_name.resize(nrows, TensorName(edge_size));
  Tn_name.resize(num_tn, TensorName(tn_num_bonds));

  for (int i = 0; i < 4; ++i) {
    C_name[i].name = "C" + std::to_string(i);
    C_name[i].bonds.resize(2);
    C_name[i].bonds[0] = bond_index++;
    C_name[i].bonds[1] = bond_index++;
    // C_name[i].bonds[0] = C_name[i].name + "_0";
    // C_name[i].bonds[1] = C_name[i].name + "_1";
  }

  for (int i = 0; i < ncols; ++i) {
    eTt_name[i].name = "eTt" + std::to_string(i);
    eTb_name[i].name = "eTb" + std::to_string(i);
    for (int j = 0; j < edge_size; ++j) {
      eTt_name[i].bonds[j] = bond_index++;
      eTb_name[i].bonds[j] = bond_index++;
      // eTt_name[i].bonds[j] = eTt_name[i].name + "_" + std::to_string(j);
      // eTb_name[i].bonds[j] = eTb_name[i].name + "_" + std::to_string(j);
    }
  }

  for (int i = 0; i < nrows; ++i) {
    eTr_name[i].name = "eTr" + std::to_string(i);
    eTl_name[i].name = "eTl" + std::to_string(i);
    for (int j = 0; j < edge_size; ++j) {
      eTr_name[i].bonds[j] = bond_index++;
      eTl_name[i].bonds[j] = bond_index++;
      // eTr_name[i].bonds[j] = eTr_name[i].name + "_" + std::to_string(j);
      // eTl_name[i].bonds[j] = eTl_name[i].name + "_" + std::to_string(j);
    }
  }

  if (is_tpo) {
    for (int i = 0; i < N_UNIT; ++i) {
      Tn_name[i].name = "Tn" + std::to_string(i);
      for (int j = 0; j < tn_num_bonds; ++j) {
        Tn_name[i].bonds[j] = bond_index++;
        // Tn_name[i].bonds[j] = Tn_name[i].name + "_" + std::to_string(j);
      }
    }
  } else {
    for (int i = 0; i < N_UNIT; ++i) {
      Tn_name[i].name = "Tn" + std::to_string(i);
      Tn_name[i + N_UNIT].name = "cTn" + std::to_string(i);
      for (int j = 0; j < tn_num_bonds; ++j) {
        Tn_name[i].bonds[j] = bond_index++;
        Tn_name[i + N_UNIT].bonds[j] = bond_index++;
        // Tn_name[i].bonds[j] = Tn_name[i].name + "_" + std::to_string(j);
        // Tn_name[i + N_UNIT].bonds[j] =
        //     Tn_name[i + N_UNIT].name + "_" + std::to_string(j);
      }
      Tn_name[i + N_UNIT].bonds[4] = Tn_name[i].bonds[4];
    }
  }

  // connect bonds of C
  C_name[0].bonds[0] = eTl_name[0].bonds[1];
  C_name[0].bonds[1] = eTt_name[0].bonds[0];
  C_name[1].bonds[0] = eTt_name[ncols - 1].bonds[1];
  C_name[1].bonds[1] = eTr_name[0].bonds[0];
  C_name[2].bonds[0] = eTr_name[nrows - 1].bonds[1];
  C_name[2].bonds[1] = eTb_name[ncols - 1].bonds[0];
  C_name[3].bonds[0] = eTb_name[0].bonds[1];
  C_name[3].bonds[1] = eTl_name[nrows - 1].bonds[0];

  // connect bonds of edge
  for (int i = 0; i < ncols; ++i) {
    eTt_name[i].bonds[2] = Tn_name[coord2index(0, i, nrows, ncols)].bonds[1];
    eTb_name[i].bonds[2] =
        Tn_name[coord2index(nrows - 1, i, nrows, ncols)].bonds[3];
    if (!is_tpo) {
      eTt_name[i].bonds[3] =
          Tn_name[coord2index(0, i, nrows, ncols) + N_UNIT].bonds[1];
      eTb_name[i].bonds[3] =
          Tn_name[coord2index(nrows - 1, i, nrows, ncols) + N_UNIT].bonds[3];
    }
    if (i != 0) {
      eTt_name[i].bonds[0] = eTt_name[i - 1].bonds[1];
      eTb_name[i].bonds[1] = eTb_name[i - 1].bonds[0];
    }
  }
  for (int i = 0; i < nrows; ++i) {
    eTr_name[i].bonds[2] =
        Tn_name[coord2index(i, ncols - 1, nrows, ncols)].bonds[2];
    eTl_name[i].bonds[2] = Tn_name[coord2index(i, 0, nrows, ncols)].bonds[0];
    if (!is_tpo) {
      eTr_name[i].bonds[3] =
          Tn_name[coord2index(i, ncols - 1, nrows, ncols) + N_UNIT].bonds[2];
      eTl_name[i].bonds[3] =
          Tn_name[coord2index(i, 0, nrows, ncols) + N_UNIT].bonds[0];
    }
    if (i != 0) {
      eTr_name[i].bonds[0] = eTr_name[i - 1].bonds[1];
      eTl_name[i].bonds[1] = eTl_name[i - 1].bonds[0];
    }
  }

  // connect bonds of Tn
  for (int r = 0; r < nrows; ++r) {
    for (int c = 0; c < ncols; ++c) {
      int index = coord2index(r, c, nrows, ncols);
      if (r > 0) {
        int index_t = coord2index(r - 1, c, nrows, ncols);
        Tn_name[index].bonds[1] = Tn_name[index_t].bonds[3];
        if (!is_tpo) {
          Tn_name[index + N_UNIT].bonds[1] = Tn_name[index_t + N_UNIT].bonds[3];
        }
      }
      if (c > 0) {
        int index_l = coord2index(r, c - 1, nrows, ncols);
        Tn_name[index].bonds[0] = Tn_name[index_l].bonds[2];
        if (!is_tpo) {
          Tn_name[index + N_UNIT].bonds[0] = Tn_name[index_l + N_UNIT].bonds[2];
        }
      }
    }
  }
}  // end of initialize_ctm

template <class tensor>
typename tensor::value_type TensorNetworkContractor<tensor>::contract(
    std::vector<tensor const*> const& C, std::vector<tensor const*> const& eTt,
    std::vector<tensor const*> const& eTr,
    std::vector<tensor const*> const& eTb,
    std::vector<tensor const*> const& eTl,
    std::vector<std::vector<tensor const*>> const& Tn) const {
  std::map<std::tuple<int, int>, tensor const*> ops;
  return contract(C, eTt, eTr, eTb, eTl, Tn, ops);
}

template <class tensor>
typename tensor::value_type TensorNetworkContractor<tensor>::contract(
    std::vector<tensor const*> const& C, std::vector<tensor const*> const& eTt,
    std::vector<tensor const*> const& eTr,
    std::vector<tensor const*> const& eTb,
    std::vector<tensor const*> const& eTl,
    std::vector<std::vector<tensor const*>> const& Tn,
    std::map<std::tuple<int, int>, tensor const*> const& ops) const {
  using mptensor::Axes;
  std::stack<Tensor_w_bondname<tensor>> st;
  for (auto action : order) {
    if (action < 0) {
      // contract and push
      Tensor_w_bondname<tensor> rhs = std::move(st.top());
      st.pop();
      Tensor_w_bondname<tensor> lhs = std::move(st.top());
      st.pop();
      st.push(::contract(lhs, rhs));
    } else {
      // push original tensor
      if (!is_mf) {
        if (action < 4) {
          st.emplace(C[action], C_name[action].bonds);
          continue;
        }
        action -= 4;
        if (action < ncols) {
          st.emplace(eTt[action], eTt_name[action].bonds);
          continue;
        }
        action -= ncols;
        if (action < nrows) {
          st.emplace(eTr[action], eTr_name[action].bonds);
          continue;
        }
        action -= nrows;
        if (action < ncols) {
          st.emplace(eTb[action], eTb_name[action].bonds);
          continue;
        }
        action -= ncols;
        if (action < nrows) {
          st.emplace(eTl[action], eTl_name[action].bonds);
          continue;
        }
        action -= nrows;
        const int N_UNIT = nrows * ncols;
        if (action < N_UNIT) {
          int x, y;
          std::tie(x, y) = index2coord(action, nrows, ncols);
          if (!is_tpo) {
            st.emplace(Tn[x][y], Tn_name[action].bonds);
          } else {
            auto it = ops.find(std::forward_as_tuple(x, y));
            if (it == ops.end()) {
              st.emplace(mptensor::contract(*Tn[x][y], Axes(4), Axes(5)),
                         Tn_name[action].bonds);
            } else {
              st.emplace(mptensor::tensordot(*Tn[x][y], *(it->second),
                                             Axes(4, 5), Axes(0, 1)),
                         Tn_name[action].bonds);
            }
          }
          continue;
        }
        action -= N_UNIT;
        if (!is_tpo && action < N_UNIT) {
          int x, y;
          std::tie(x, y) = index2coord(action, nrows, ncols);
          auto it = ops.find(std::forward_as_tuple(x, y));
          if (it == ops.end()) {
            st.emplace(conj(*Tn[x][y]), Tn_name[action + N_UNIT].bonds);
          } else {
            st.emplace(mptensor::tensordot(conj(*Tn[x][y]), *(it->second),
                                           Axes(4), Axes(1)),
                       Tn_name[action + N_UNIT].bonds);
          }
          continue;
        }
        const int action_orig =
            action + 4 + ncols + nrows + ncols + nrows + N_UNIT;
        throw std::runtime_error("invalid action " +
                                 std::to_string(action_orig) +
                                 " in contraction order");
      } else {  // mean field
        if (is_tpo) {
          throw std::runtime_error("not implemented");
        }
        const int N_UNIT = nrows * ncols;
        if (action < N_UNIT) {
          // ket
          int x, y;
          std::tie(x, y) = index2coord(action, nrows, ncols);
          st.emplace(Tn[x][y], Tn_name[action].bonds);
          continue;
        }
        action -= N_UNIT;
        if (action < N_UNIT) {
          // bra
          int x, y;
          std::tie(x, y) = index2coord(action, nrows, ncols);
          auto it = ops.find(std::forward_as_tuple(x, y));
          if (it == ops.end()) {
            st.emplace(conj(*Tn[x][y]), Tn_name[action + N_UNIT].bonds);
          } else {
            st.emplace(mptensor::tensordot(conj(*Tn[x][y]), *(it->second),
                                           Axes(4), Axes(1)),
                       Tn_name[action + N_UNIT].bonds);
          }
          continue;
        }
        const int action_orig = action + N_UNIT;
        throw std::runtime_error("invalid action " +
                                 std::to_string(action_orig) +
                                 " in contraction order");
      }
    }
  }
  if (st.size() != 1) {
    throw std::runtime_error("invalid contraction order");
  }
  if (st.top().t.rank() != 0) {
    throw std::runtime_error("invalid contraction order");
  }

  return mptensor::trace(st.top().t, Axes(), Axes());
}

template <class tensor>
std::vector<int> TensorNetworkContractor<tensor>::get_order() const {
  return order;
}

template class TensorNetworkContractor<real_tensor>;
template class TensorNetworkContractor<complex_tensor>;


}  // namespace itps
}  // namespace tenes
