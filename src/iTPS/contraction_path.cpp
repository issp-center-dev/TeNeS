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

#include "../find_contraction_path/finder.hpp"

#include <stack>

#include "iTPS.hpp"
#include "contraction_path.hpp"

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
  path.clear();

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

  find_contraction_path::TensorNetwork tnw;
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
  path = std::get<0>(res);
  assert(check_path(path));
  return res;
}

template <class tensor>
bool TensorNetworkContractor<tensor>::check_path(
    std::vector<int> const& path) const {
  size_t num_tensors = 0;
  num_tensors += C_name.size();
  num_tensors += eTt_name.size();
  num_tensors += eTr_name.size();
  num_tensors += eTb_name.size();
  num_tensors += eTl_name.size();
  num_tensors += Tn_name.size();

  std::vector<int> t(num_tensors, 0);

  size_t num_contract = 0;
  for (auto p : path) {
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
void TensorNetworkContractor<tensor>::load_path(std::vector<int> const& path) {
  if (!check_path(path)) {
    throw std::runtime_error("invalid contraction path");
  }
  this->path = path;
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
  for (auto action : path) {
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
                                 " in contraction path");
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
                                 " in contraction path");
      }
    }
  }
  if (st.size() != 1) {
    throw std::runtime_error("invalid contraction path");
  }
  if (st.top().t.rank() != 0) {
    throw std::runtime_error("invalid contraction path");
  }

  return mptensor::trace(st.top().t, Axes(), Axes());
}

template <class tensor>
std::vector<int> TensorNetworkContractor<tensor>::get_path() const {
  return path;
}

template class TensorNetworkContractor<real_tensor>;
template class TensorNetworkContractor<complex_tensor>;

std::vector<std::vector<std::vector<int>>> make_default_path_itpo_ctm() {
  std::vector<std::vector<std::vector<int>>> path_table(
      4, std::vector<std::vector<int>>(4));

  // nrow = 1
  path_table[0][0] = {};
  path_table[0][1] = {
      10, 4, 0, 9, -1, -1, 3,  7,  -1, 5,  8,  11,
      -1, 1, 2, 6, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[0][2] = {
      12, 4, 0, 11, -1, -1, 3,  8,  -1, 5,  13, 9,  6,  10, 14,
      -1, 1, 2, 7,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[0][3] = {
      14, 4, 0,  13, -1, -1, 3,  9,  -1, 5,  15, 10, 6,  16, 11, 7,  12, 1,
      2,  8, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };

  // nrow = 2
  path_table[1][0] = {
      10, 4, 1, 5, -1, -1, 0,  8,  -1, 6,  9,  11,
      -1, 2, 3, 7, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[1][1] = {
      12, 4,  0,  10, -1, -1, 1,  5, -1, 6,  13, -1, -1, 3,  8,  -1,
      11, 14, -1, -1, 2,  7,  -1, 9, 15, -1, -1, -1, -1, -1, -1,
  };
  path_table[1][2] = {
      14, 4,  0,  12, -1, -1, 3,  9,  -1, 13, 17, -1, -1,
      5,  15, 18, 10, 1,  6,  -1, 7,  16, -1, -1, 2,  8,
      -1, 11, 19, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[1][3] = {
      16, 4,  0,  14, -1, -1, 3,  10, -1, 15, 20, -1, -1, 5,  17, 21,
      11, 6,  18, 22, 12, 1,  7,  -1, 8,  19, -1, -1, 2,  9,  -1, 13,
      23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };

  // nrow = 3
  path_table[2][0] = {
      12, 4, 1, 5, -1, -1, 0,  9,  -1, 6,  13, 10, 7,  11, 14,
      -1, 2, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[2][1] = {
      14, 4,  0,  11, -1, -1, 1,  5,  -1, 6,  15, -1, -1,
      7,  17, 16, 12, 3,  9,  -1, 13, 18, -1, -1, 2,  8,
      -1, 10, 19, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[2][2] = {
      16, 4,  0,  13, -1, -1, 5,  17, 6,  18, 1,  7,  -1, 8,  21, 20, 19,
      14, 2,  9,  -1, 3,  10, 15, 22, -1, -1, -1, 11, 23, -1, 12, 24, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[2][3] = {
      18, 4,  0,  15, -1, -1, 16, 22, 17, 26, 3,  11, -1, 5,  19,
      23, 27, 12, 6,  20, 24, 28, 13, 1,  7,  -1, 8,  21, -1, 9,
      25, -1, -1, 2,  10, 14, 29, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };

  // nrow = 4
  path_table[3][0] = {
      14, 4, 1,  5,  -1, -1, 0,  10, -1, 6,  15, 11, 7,  16, 12, 8,  13, 2,
      3,  9, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[3][1] = {
      16, 4,  0,  12, -1, -1, 1,  5,  -1, 6,  17, -1, -1, 7,  19, 18,
      13, 8,  21, 20, 14, 3,  10, -1, 15, 22, -1, -1, 2,  9,  -1, 11,
      23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[3][2] = {
      18, 4,  0,  14, -1, -1, 5,  19, 6,  20, 1,  7,  -1, 8,  23,
      22, 21, 15, 9,  26, 25, 24, 16, 2,  10, -1, 3,  11, 17, 27,
      -1, -1, -1, 12, 28, -1, 13, 29, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[3][3] = {
      20, 4,  0,  16, -1, -1, 5,  21, 6,  22, 1,  7,  8,  23, -1, -1, -1, 9,
      27, 26, 25, 24, 17, 10, 31, 30, 29, 28, 18, 11, 35, 32, 33, 34, -1, -1,
      14, 2,  15, -1, -1, 13, 12, 3,  19, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };

  return path_table;
}

std::vector<std::vector<std::vector<int>>> make_default_path_itps_ctm() {
  std::vector<std::vector<std::vector<int>>> path_table(
      4, std::vector<std::vector<int>>(4));

  // nrow = 1
  path_table[0][0] = {};
  path_table[0][1] = {
      0, 4,  9,  10, 12, 3,  7,  -1, 5,  11, 1,  6,  -1, 2,
      8, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[0][2] = {
      0,  4,  11, 12, 15, 3,  8,  -1, 5,  13, 16, 9,  14, 6,  1,  7,  -1, 2,
      10, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[0][3] = {
      0,  4,  13, 14, 18, 3,  9,  -1, 5,  15, 19, 10, 6,  16, 20,
      11, 17, 7,  1,  8,  -1, 2,  12, 21, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };

  // nrow = 2
  path_table[1][0] = {
      0,  4,  1,  5,  10, 12, 8,  6,  11, 2,  7,  -1, 3,  9,
      13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[1][1] = {
      0,  4,  10, 12, 16, 5,  13, 17, 1,  6,  -1, 3,  8,
      -1, 11, 14, 18, -1, -1, -1, 2,  7,  -1, 9,  15, 19,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[1][2] = {
      0,  4,  12, 14, 20, 13, 17, 23, 3,  9,  -1, 5,  15, 21, 18, 24, 10,
      1,  6,  -1, 7,  16, 22, -1, -1, -1, 2,  8,  -1, 11, 19, 25, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[1][3] = {
      5,  17, 25, 21, 29, 11, 0,  4,  -1, 14, 16, 24, -1, -1, -1, 3,
      10, -1, 15, 20, 28, -1, -1, -1, -1, -1, -1, -1, -1, -1, 6,  18,
      26, 22, 30, 12, 1,  7,  -1, 8,  19, 27, -1, -1, -1, 2,  9,  -1,
      13, 23, 31, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };

  // nrow = 3
  path_table[2][0] = {
      0,  4,  1,  5,  12, 15, 9,  6,  13, 16, 10, 14, 7,  2,  8,  -1, 3,  11,
      17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[2][1] = {
      0,  4,  11, 14, 20, 5,  15, 21, 1,  6,  -1, 7,  17, 23, 16, 22, 12,
      3,  9,  -1, 13, 18, 24, -1, -1, -1, 2,  8,  -1, 10, 19, 25, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[2][2] = {
      0,  4,  13, 16, 25, 5,  17, 26, 1,  6,  18, -1, 7,  27, -1, -1, -1,
      8,  21, 30, 20, 29, 19, 28, 14, 3,  10, 22, -1, 15, 31, -1, -1, -1,
      23, 32, 2,  11, 9,  24, -1, 12, 33, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[2][3] = {
      0,  4,  15, 18, 30, 16, 22, 34, 3,  11, 26, -1, 17, 38, -1, -1, -1,
      5,  19, 31, 23, 35, 27, 39, 12, 6,  20, 32, 24, 36, 28, 40, 13, 1,
      7,  21, -1, 8,  33, -1, -1, -1, 25, 37, 29, 2,  9,  10, 14, 41, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };

  // nrow = 4
  path_table[3][0] = {
      0,  4,  1,  5,  14, 18, 10, 6,  15, 19, 11, 7,  16, 20, 12,
      17, 8,  2,  9,  -1, 3,  13, 21, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[3][1] = {
      7,  19, 27, 18, 26, 13, 0,  4,  -1, 12, 16, 24, -1, -1, -1, 1,
      5,  -1, 6,  17, 25, -1, -1, -1, -1, -1, -1, -1, -1, -1, 8,  21,
      29, 20, 28, 14, 3,  10, -1, 15, 22, 30, -1, -1, -1, 2,  9,  -1,
      11, 23, 31, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[3][2] = {
      0,  4,  14, 18, 30, 5,  19, 31, 1,  6,  20, -1, 7,  32, -1, -1, -1,
      8,  23, 35, 22, 34, 21, 33, 15, 9,  26, 38, 25, 37, 24, 36, 16, 3,
      11, 27, -1, 17, 39, -1, -1, -1, 28, 40, 29, 2,  12, 10, 13, 41, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[3][3] = {
      0,  4,  16, 20, 36, 5,  21, 37, 6,  22, 38, 1,  7,  23, -1, 8,  39, -1,
      -1, -1, 9,  27, 43, 26, 42, 25, 41, 24, 40, 17, 10, 31, 47, 30, 46, 29,
      45, 28, 44, 18, 3,  12, 32, -1, 19, 48, -1, -1, -1, 33, 49, 35, 51, -1,
      11, 2,  15, -1, -1, 13, 14, -1, 34, 50, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };

  return path_table;
}

std::vector<std::vector<std::vector<int>>> make_default_path_itps_mf() {
  std::vector<std::vector<std::vector<int>>> path_table(
      5, std::vector<std::vector<int>>(5));

  // nrow = 1
  path_table[0][0] = {};
  path_table[0][1] = {
      0, 2, -1, 1, 3, -1, -1,
  };
  path_table[0][2] = {
      0, 3, -1, 1, 4, 2, 5, -1, -1, -1, -1,
  };
  path_table[0][3] = {
      0, 4, -1, 1, 5, -1, 3, 2, 6, 7, -1, -1, -1, -1, -1,
  };
  path_table[0][4] = {
      0, 5, -1, 4, 9, -1, 1, 6, -1, 2, 3, 7, 8, -1, -1, -1, -1, -1, -1,
  };

  // nrow = 2
  path_table[1][0] = {
      0, 2, -1, 1, 3, -1, -1,
  };
  path_table[1][1] = {
      0, 4, -1, 1, 5, -1, 2, 6, -1, 3, 7, -1, -1, -1, -1,
  };
  path_table[1][2] = {
      0, 6,  -1, 3,  9,  -1, 1,  7,  -1, 2,  8,  -1,
      4, 10, 5,  11, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[1][3] = {
      0,  8,  -1, 4, 12, -1, 1,  9,  -1, 5,  13, -1, 2,  10, -1, 6,
      14, -1, 3,  7, 11, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[1][4] = {
      0,  10, -1, 5,  15, -1, 1,  11, -1, 6,  16, -1, 2,
      12, -1, 7,  17, -1, 3,  8,  13, 18, -1, 4,  9,  14,
      19, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };

  // nrow = 3
  path_table[2][0] = {
      0, 3, -1, 1, 4, 2, 5, -1, -1, -1, -1,
  };
  path_table[2][1] = {
      0, 6, -1, 1,  7,  -1, 2,  8,  -1, 4,  10, -1,
      3, 9, 5,  11, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[2][2] = {
      0,  9,  -1, 1, 10, -1, 2,  11, -1, 3,  12, -1, 6,  15, -1, 4,  5,  13,
      14, -1, -1, 8, 7,  16, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[2][3] = {
      0,  12, -1, 4,  16, -1, 8,  20, -1, 1,  13, -1, 5,  17, 9,  21,
      -1, 2,  14, -1, 6,  10, 18, 22, -1, 3,  7,  -1, 11, 15, 19, 23,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[2][4] = {
      0,  15, -1, 5,  20, -1, 10, 25, -1, 1,  16, -1, 6,  21, 11,
      26, -1, 2,  17, -1, 7,  22, 12, 27, -1, 3,  18, -1, 8,  13,
      28, 23, 24, -1, 19, 29, 4,  9,  14, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };

  // nrow = 4
  path_table[3][0] = {
      0, 4, -1, 1, 5, -1, 3, 2, 6, 7, -1, -1, -1, -1, -1,
  };
  path_table[3][1] = {
      0,  8,  -1, 1, 9,  -1, 2,  10, -1, 3,  11, -1, 4,  12, -1, 5,
      13, -1, 6,  7, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[3][2] = {
      0,  12, -1, 1,  13, -1, 2,  14, -1, 3,  15, -1, 4,  16, 5,  17,
      -1, 6,  18, -1, 7,  8,  19, 20, -1, 9,  10, -1, 11, 21, 22, 23,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[3][3] = {
      0,  16, -1, 1,  17, -1, 4,  20, -1, 5,  21, 2,  18, -1, 6,  22,
      3,  7,  19, 23, -1, -1, -1, 8,  24, -1, 9,  25, 12, 13, 28, 29,
      -1, -1, -1, 10, 11, 27, 31, -1, -1, 14, 15, -1, 26, 30, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[3][4] = {
      0,  20, -1, 1,  21, -1, 5,  25, -1, 6,  26, 10, 30, -1, 11, 31,
      15, 16, 35, 36, -1, -1, -1, 2,  22, -1, 7,  27, 12, 32, 17, 37,
      -1, 3,  23, -1, 8,  28, 13, 18, 24, 29, -1, 4,  9,  14, -1, -1,
      -1, 19, 33, 34, 38, 39, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };

  // nrow = 5
  path_table[4][0] = {
      0, 5, -1, 4, 9, -1, 1, 6, -1, 2, 3, 7, 8, -1, -1, -1, -1, -1, -1,
  };
  path_table[4][1] = {
      0,  10, -1, 1,  11, -1, 2,  12, -1, 3,  13, -1, 4,
      14, -1, 5,  15, -1, 6,  7,  16, 17, -1, 8,  9,  18,
      19, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[4][2] = {
      0,  15, -1, 1,  16, -1, 2,  17, -1, 3,  18, -1, 4,  19, 5,
      20, -1, 6,  21, -1, 7,  22, 8,  23, -1, 9,  24, -1, 10, 11,
      26, 25, 28, -1, 27, 29, 12, 13, 14, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[4][3] = {
      0,  20, -1, 1,  21, -1, 4,  24, -1, 5,  25, 2,  22, -1, 6,  26,
      3,  7,  23, 27, -1, -1, -1, 8,  28, -1, 9,  29, 10, 30, 11, 31,
      -1, 12, 32, -1, 13, 33, 14, 15, 36, 37, -1, 16, 17, 18, -1, -1,
      -1, 19, 34, 35, 38, 39, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  path_table[4][4] = {
      0,  25, -1, 1,  26, -1, 5,  30, -1, 6,  31, 2,  27, -1, 7,  32, 3,
      28, -1, 8,  33, 4,  9,  29, 34, -1, -1, -1, 10, 35, -1, 11, 36, 12,
      37, 13, 38, 14, 39, -1, 15, 40, -1, 16, 41, 20, 21, 45, 46, -1, -1,
      -1, 17, 22, 44, 18, 19, -1, -1, 24, 43, 42, 47, -1, -1, 23, 48, 49,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };

  return path_table;
}

std::vector<int> default_path(int nrows, int ncols, bool is_tpo, bool is_mf) {
  static std::vector<std::vector<std::vector<int>>> path_table_itpo_ctm =
      make_default_path_itpo_ctm();
  static std::vector<std::vector<std::vector<int>>> path_table_itps_ctm =
      make_default_path_itps_ctm();
  static std::vector<std::vector<std::vector<int>>> path_table_itps_mf =
      make_default_path_itps_mf();
  if (is_tpo) {
    if (is_mf) {
      throw std::runtime_error("not implemented");
    } else {
      if (nrows > 4 || ncols > 4) {
        return {};
      }
      return path_table_itpo_ctm[nrows - 1][ncols - 1];
    }
  } else {
    if (is_mf) {
      if (nrows > 5 || ncols > 5) {
        return {};
      }
      return path_table_itps_mf[nrows - 1][ncols - 1];
    } else {
      if (nrows > 4 || ncols > 4) {
        return {};
      }
      return path_table_itps_ctm[nrows - 1][ncols - 1];
    }
  }
}

}  // namespace itps
}  // namespace tenes
