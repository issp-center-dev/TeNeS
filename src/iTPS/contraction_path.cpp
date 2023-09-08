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
  std::vector<std::string> bonds;
  Tensor_w_bondname(tensor const* t_ptr_, std::vector<std::string> const& b_)
      : t_ptr(t_ptr_), bonds(b_) {}
  Tensor_w_bondname(tensor const& t_, std::vector<std::string> const& b_)
      : t(t_), t_ptr(nullptr), bonds(b_) {}
};

template <class tensor>
Tensor_w_bondname<tensor> contract(Tensor_w_bondname<tensor> const& lhs,
                                   Tensor_w_bondname<tensor> const& rhs) {
  mptensor::Axes axes_lhs, axes_rhs;
  const size_t nb_lhs = lhs.bonds.size();
  const size_t nb_rhs = rhs.bonds.size();

  std::vector<std::string> bonds_lhs = lhs.bonds;
  std::vector<std::string> bonds_rhs = rhs.bonds;

  size_t nb_common = 0;
  for (size_t i = 0; i < nb_lhs; ++i) {
    const std::string& bond = lhs.bonds[i];
    for (size_t j = 0; j < nb_rhs; ++j) {
      if (bond == rhs.bonds[j]) {
        axes_lhs.push(i);
        axes_rhs.push(j);
        bonds_lhs[i] = "";
        bonds_rhs[j] = "";
        ++nb_common;
        break;
      }
    }
  }
  std::vector<std::string> bonds;
  bonds.reserve(nb_lhs + nb_rhs - 2 * nb_common);
  for (size_t i = 0; i < nb_lhs; ++i) {
    if (bonds_lhs[i] != "") {
      bonds.push_back(bonds_lhs[i]);
    }
  }
  for (size_t i = 0; i < nb_rhs; ++i) {
    if (bonds_rhs[i] != "") {
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
TensorNetworkContractor<tensor>::TensorNetworkContractor(
    int nrows, int ncols, std::vector<int> const& shape_types,
    std::vector<std::vector<int>> const& shape_dims, int chi, bool is_tpo,
    bool is_mf) {
  initialize(nrows, ncols, shape_types, shape_dims, chi, is_tpo, is_mf);
}

template <class tensor>
void TensorNetworkContractor<tensor>::initialize(
    int nrows, int ncols, std::vector<int> const& shape_types,
    std::vector<std::vector<int>> const& shape_dims, int chi, bool is_tpo,
    bool is_mf) {
  this->nrows = nrows;
  this->ncols = ncols;
  this->is_tpo = is_tpo;
  this->is_mf = is_mf;
  path.clear();
  tnw = find_contraction_path::TensorNetwork();

  if (is_mf) {
    initialize_mf(shape_types, shape_dims);
  } else {
    initialize_ctm(shape_types, shape_dims, chi);
  }

  const int tn_size = (is_tpo ? 4 : 5);

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
      for (int j = 0; j < tn_size; ++j) {
        tnw.set_bond_dim(Tn_name[i].bonds[j], dims[j]);
      }
    } else {
      for (int j = 0; j < tn_size; ++j) {
        tnw.set_bond_dim(Tn_name[i].bonds[j], dims[j]);
        tnw.set_bond_dim(Tn_name[i + N_UNIT].bonds[j], dims[j]);
      }
    }
  }
}

template <class tensor>
void TensorNetworkContractor<tensor>::optimize() {
  auto res = tnw.optimize();
  path = std::move(std::get<0>(res));
  assert(check_path(path));
  // double cost = std::get<1>(res);
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
    std::cerr << "ERROR: number of -1 is " << num_contract << " (expected: " << num_tensors-1 <<  ")" << std::endl;
    res = false;
  }
  for (size_t i=0; i < num_tensors; ++i) {
    int n = t[i];
    if (n != 1) {
      std::cerr << "ERROR: number of " << i << " is " << n << " (expected: 1)" << std::endl;
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
void TensorNetworkContractor<tensor>::initialize_mf(
    std::vector<int> const& shape_types,
    std::vector<std::vector<int>> const& shape_dims) {
  if (is_tpo) {
    throw std::runtime_error("not implemented");
  }

  const int tn_size = 5;
  const int N_UNIT = ncols * nrows;
  const int tn_num = 2 * N_UNIT;
  Tn_name.resize(tn_num, TensorName(tn_size));

  for (int i = 0; i < N_UNIT; ++i) {
    Tn_name[i].name = "Tn" + std::to_string(i);
    Tn_name[i + N_UNIT].name = "cTn" + std::to_string(i);
    for (int j = 0; j < tn_size; ++j) {
      Tn_name[i].bonds[j] = "Tn" + std::to_string(i) + "_" + std::to_string(j);
      Tn_name[i + N_UNIT].bonds[j] =
          "cTn" + std::to_string(i) + "_" + std::to_string(j);
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
        if (r == nrows - 1) {
          Tn_name[index + N_UNIT].bonds[3] = Tn_name[index].bonds[3];
        }
      }

      if (c == 0) {
        Tn_name[index + N_UNIT].bonds[0] = Tn_name[index].bonds[0];
      } else {
        int index_l = coord2index(r, c - 1, nrows, ncols);
        Tn_name[index].bonds[0] = Tn_name[index_l].bonds[2];
        Tn_name[index + N_UNIT].bonds[0] = Tn_name[index_l + N_UNIT].bonds[2];
        if (c == ncols - 1) {
          Tn_name[index + N_UNIT].bonds[2] = Tn_name[index].bonds[2];
        }
      }
    }
  }
}

template <class tensor>
void TensorNetworkContractor<tensor>::initialize_ctm(
    std::vector<int> const& shape_types,
    std::vector<std::vector<int>> const& shape_dims, int chi) {
  const int edge_size = (is_tpo ? 3 : 4);
  const int tn_size = (is_tpo ? 4 : 5);
  const int N_UNIT = ncols * nrows;
  const int tn_num = (is_tpo ? nrows * ncols : nrows * ncols * 2);

  C_name.resize(4);
  eTt_name.resize(ncols, TensorName(edge_size));
  eTr_name.resize(nrows, TensorName(edge_size));
  eTb_name.resize(ncols, TensorName(edge_size));
  eTl_name.resize(nrows, TensorName(edge_size));
  Tn_name.resize(tn_num, TensorName(tn_size));

  for (int i = 0; i < 4; ++i) {
    C_name[i].name = "C" + std::to_string(i);
    C_name[i].bonds.resize(2);
    C_name[i].bonds[0] = "C" + std::to_string(i) + "_0";
    C_name[i].bonds[1] = "C" + std::to_string(i) + "_1";
  }

  for (int i = 0; i < ncols; ++i) {
    eTt_name[i].name = "eTt" + std::to_string(i);
    eTb_name[i].name = "eTb" + std::to_string(i);
    for (int j = 0; j < edge_size; ++j) {
      eTt_name[i].bonds[j] =
          "eTt" + std::to_string(i) + "_" + std::to_string(j);
      eTb_name[i].bonds[j] =
          "eTb" + std::to_string(i) + "_" + std::to_string(j);
    }
  }

  for (int i = 0; i < nrows; ++i) {
    eTr_name[i].name = "eTr" + std::to_string(i);
    eTl_name[i].name = "eTl" + std::to_string(i);
    for (int j = 0; j < edge_size; ++j) {
      eTr_name[i].bonds[j] =
          "eTr" + std::to_string(i) + "_" + std::to_string(j);
      eTl_name[i].bonds[j] =
          "eTl" + std::to_string(i) + "_" + std::to_string(j);
    }
  }

  if (is_tpo) {
    for (int i = 0; i < N_UNIT; ++i) {
      Tn_name[i].name = "Tn" + std::to_string(i);
      for (int j = 0; j < tn_size; ++j) {
        Tn_name[i].bonds[j] =
            "Tn" + std::to_string(i) + "_" + std::to_string(j);
      }
    }
  } else {
    for (int i = 0; i < N_UNIT; ++i) {
      Tn_name[i].name = "Tn" + std::to_string(i);
      Tn_name[i + N_UNIT].name = "cTn" + std::to_string(i);
      for (int j = 0; j < tn_size; ++j) {
        Tn_name[i].bonds[j] =
            "Tn" + std::to_string(i) + "_" + std::to_string(j);
        Tn_name[i + N_UNIT].bonds[j] =
            "cTn" + std::to_string(i) + "_" + std::to_string(j);
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
              tensor T = mptensor::tensordot(*Tn[x][y], *(it->second), Axes(5),
                                             Axes(1));
              st.emplace(mptensor::contract(T, Axes(4), Axes(5)),
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
          int x, y;
          std::tie(x, y) = index2coord(action, nrows, ncols);
          st.emplace(Tn[x][y], Tn_name[action].bonds);
          continue;
        }
        action -= N_UNIT;
        if (action < N_UNIT) {
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
  return mptensor::trace(st.top().t);
}

template class TensorNetworkContractor<real_tensor>;
template class TensorNetworkContractor<complex_tensor>;

}  // namespace itps
}  // namespace tenes
