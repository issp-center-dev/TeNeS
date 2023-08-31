/*
 * finder.cpp and finder.hpp are translated from the original Python code,
 * Tensordot (https://github.com/smorita/Tensordot) written by Satoshi Morita.
 * The original license is the MIT License.
 */

#include "../icecream.cpp"

#include "finder.hpp"

#include <iostream>
#include <cassert>

#include <vector>
#include <set>
#include <map>
#include <string>
#include <algorithm>

#include <fstream>
#include <sstream>
#include <iterator>
#include <limits>

namespace tenes {
namespace find_contraction_path {

inline bool are_direct_product(const TensorFrame& t1, const TensorFrame& t2) {
  for (const auto& b : t1.bonds) {
    if (t2.bonds.find(b) != t2.bonds.end()) {
      return false;
    }
  }
  return true;
}

inline bool are_overlap(const TensorFrame& t1, const TensorFrame& t2) {
  return (t1.bits & t2.bits) > 0;
}

double get_contracting_cost(const TensorFrame& t1, const TensorFrame& t2,
                            const std::vector<int>& bond_dims) {
  double cost = 1.0;
  std::set<int> union_bonds(t1.bonds);
  union_bonds.insert(t2.bonds.begin(), t2.bonds.end());
  for (const auto& b : union_bonds) {
    cost *= bond_dims[b];
  }
  cost += t1.cost + t2.cost;
  return cost;
}

TensorFrame contract_tf(const TensorFrame& t1, const TensorFrame& t2,
                        const std::vector<int>& bond_dims) {
  assert(!are_direct_product(t1, t2));

  std::vector<int> rpn = t1.rpn;
  rpn.insert(rpn.end(), t2.rpn.begin(), t2.rpn.end());
  rpn.push_back(-1);

  // XOR bits and bonds
  boost::multiprecision::cpp_int bits = t1.bits ^ t2.bits;
  std::set<int> bonds = t1.bonds;
  for (const auto& b : t2.bonds) {
    if (bonds.find(b) == bonds.end()) {
      bonds.insert(b);
    } else {
      bonds.erase(b);
    }
  }

  double cost = get_contracting_cost(t1, t2, bond_dims);
  return TensorFrame(rpn, bits, bonds, cost);
}

void TensorNetwork::add_tensor(const std::string& t_name,
                               const std::vector<std::string>& b_names) {
  int t_index = tensors.size();
  std::vector<int> b_indices;
  for (const auto& b : b_names) {
    auto it = std::find(bond_names.begin(), bond_names.end(), b);
    int i = std::distance(bond_names.begin(), it);
    if (it == bond_names.end()) {
      bonds.push_back(Bond());
      bond_names.push_back(b);
      bond_dims.push_back(default_bond_dim);
    }
    bonds[i].connect(t_index);
    b_indices.push_back(i);
  }
  tensor_names.push_back(t_name);
  tensors.push_back(Tensor(t_index, b_indices));
}

TensorNetwork TensorNetwork::clone() const {
  TensorNetwork tn;
  tn.total_memory = total_memory;
  tn.max_memory = max_memory;
  tn.cpu_cost = cpu_cost;
  for (const auto& b : bonds) {
    tn.bonds.push_back(Bond(b.t0, b.t1));
  }
  for (const auto& t : tensors) {
    tn.tensors.push_back(Tensor(t.name, t.bonds));
  }
  return tn;
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>
TensorNetwork::find_bonds(int tensor_a, int tensor_b) const {
  const auto& bonds_a = tensors[tensor_a].bonds;
  const auto& bonds_b = tensors[tensor_b].bonds;

  std::vector<int> contracts;
  std::vector<int> replaced_a;
  std::vector<int> replaced_b;

  for (const auto& b : bonds_a) {
    if (std::find(bonds_b.begin(), bonds_b.end(), b) != bonds_b.end()) {
      contracts.push_back(b);
    } else {
      replaced_a.push_back(b);
    }
  }

  for (const auto& b : bonds_b) {
    if (std::find(bonds_a.begin(), bonds_a.end(), b) == bonds_a.end()) {
      replaced_b.push_back(b);
    }
  }

  return std::make_tuple(contracts, replaced_a, replaced_b);
}

TensorNetwork TensorNetwork::contract(int t0, int t1,
                                      const std::vector<int>& bc,
                                      const std::vector<int>& br0,
                                      const std::vector<int>& br1) {
  TensorNetwork tn = clone();

  // Create the contracted tensor
  Tensor& t_new = tn.tensors[t0];
  // Change names of tensors using Reverse Polish Notation
  t_new.name.insert(t_new.name.end(), tensors[t0].name.begin(),
                    tensors[t0].name.end());
  t_new.name.insert(t_new.name.end(), tensors[t1].name.begin(),
                    tensors[t1].name.end());
  t_new.name.push_back(-1);
  // Remove contracted bonds
  for (const auto& b : bc) {
    t_new.bonds.erase(std::remove(t_new.bonds.begin(), t_new.bonds.end(), b),
                      t_new.bonds.end());
  }
  // Add bonds from deleted tensor
  for (const auto& b : br1) {
    t_new.bonds.push_back(b);
  }

  // Clear the removed tensor
  tn.tensors[t1] = Tensor();

  // Update bonds
  std::vector<Bond>& bonds = tn.bonds;
  // Remove contracted bonds from the bond list
  for (const auto& b : bc) {
    bonds[b].t0 = bonds[b].t1 = -1;
  }
  // Change bond connections
  int old_idx = t1;
  int new_idx = t0;
  for (const auto& b : br1) {
    if (bonds[b].t0 == old_idx) {
      bonds[b].t0 = new_idx;
    } else if (bonds[b].t1 == old_idx) {
      bonds[b].t1 = new_idx;
    }
  }

  return tn;
}

double TensorNetwork::calc_memory(const std::vector<int>& rpn) const {
  TensorNetwork tn = clone();
  std::vector<std::tuple<int, double, double, double>> cost;

  for (const auto& item : rpn) {
    if (item == -1) {
      auto c1 = cost.back();
      cost.pop_back();
      auto c0 = cost.back();
      cost.pop_back();

      int index1 = std::get<0>(c1);
      int index0 = std::get<0>(c0);

      auto bc_br0_br1 = tn.find_bonds(index0, index1);
      const auto& bc = std::get<0>(bc_br0_br1);
      const auto& br0 = std::get<1>(bc_br0_br1);
      const auto& br1 = std::get<2>(bc_br0_br1);

      double mem_start = std::get<2>(c0) + std::get<2>(c1);
      double mem_end = 1.0;
      for (const auto& b : br0) {
        mem_end *= bond_dims[b];
      }
      for (const auto& b : br1) {
        mem_end *= bond_dims[b];
      }
      double mem_req = std::max(
          {std::get<1>(c0) + std::get<2>(c1), std::get<1>(c0) + std::get<3>(c1),
           std::get<2>(c0) + std::get<1>(c1), std::get<3>(c0) + std::get<1>(c1),
           mem_end + std::get<3>(c0) + std::get<3>(c1)});

      tn = tn.contract(index0, index1, bc, br0, br1);

      cost.push_back(std::make_tuple(index0, mem_req, mem_start, mem_end));
    } else {
      const auto& t = tn.tensors[item];
      double val = 1.0;
      for (const auto& b : t.bonds) {
        val *= bond_dims[b];
      }
      cost.push_back(std::make_tuple(item, val, val, val));
    }
  }

  return std::get<1>(cost[0]);
}

void TensorNetwork::set_bond_dim(const std::string& bond_name, int dim) {
  auto it = std::find(bond_names.begin(), bond_names.end(), bond_name);
  if (it != bond_names.end()) {
    int index = std::distance(bond_names.begin(), it);
    bond_dims[index] = dim;
  }
}

TensorNetwork TensorNetwork::from_file(const std::string& filename) {
  TensorNetwork tn;
  std::ifstream infile(filename);
  std::string line;

  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    std::vector<std::string> data((std::istream_iterator<std::string>(iss)),
                                  std::istream_iterator<std::string>());

    if (data.empty()) {
      continue;
    }

    std::string command = data[0];
    std::transform(command.begin(), command.end(), command.begin(), ::tolower);

    if (command == "default_dimension") {
      tn.default_bond_dim = std::stoi(data[1]);
    } else if (command == "tensor") {
      tn.add_tensor(data[1],
                    std::vector<std::string>(data.begin() + 2, data.end()));
    } else if (command == "bond") {
      for (size_t i = 1; i < data.size() - 1; ++i) {
        tn.set_bond_dim(data[i], std::stoi(data.back()));
      }
    } else if (command == "bond_dim") {
      for (size_t i = 2; i < data.size(); ++i) {
        tn.set_bond_dim(data[i], std::stoi(data[1]));
      }
    }
  }

  infile.close();
  return tn;
}

std::vector<std::map<boost::multiprecision::cpp_int, TensorFrame>>
TensorNetwork::init_tensordict_of_size() {
  namespace mp = boost::multiprecision;
  std::vector<std::map<mp::cpp_int, TensorFrame>> tensordict_of_size(
      tensors.size() + 1);

  for (const auto& t : tensors) {
    const auto& rpn = t.name;
    mp::cpp_int bits = mp::cpp_int(1) << rpn[0];
    std::set<int> bonds(t.bonds.begin(), t.bonds.end());
    double cost = 0.0;
    tensordict_of_size[1][bits] = TensorFrame(rpn, bits, bonds, cost);
  }

  return tensordict_of_size;
}

std::pair<std::vector<int>, double> TensorNetwork::optimize() {
  namespace mp = boost::multiprecision;
  auto tensordict_of_size = init_tensordict_of_size();

  int n = tensors.size();
  double xi_min = *std::min_element(bond_dims.begin(), bond_dims.end());
  double mu_cap = 1.0;
  double prev_mu_cap = 0.0;

  while (tensordict_of_size.back().size() < 1) {
    double next_mu_cap = std::numeric_limits<double>::max();

    for (int c = 2; c <= n; ++c) {
      for (int d1 = 1; d1 <= c / 2; ++d1) {
        int d2 = c - d1;

        for (const auto& t1 : tensordict_of_size[d1]) {
          for (const auto& t2 : tensordict_of_size[d2]) {
            if (are_overlap(t1.second, t2.second)) continue;
            if (are_direct_product(t1.second, t2.second)) continue;

            double cost = get_contracting_cost(t1.second, t2.second, bond_dims);
            mp::cpp_int bits = t1.first ^ t2.first;

            if (next_mu_cap <= cost){
              continue;
            } else if (mu_cap < cost){
              next_mu_cap = cost;
            } else if (t1.second.is_new || t2.second.is_new ||
                     prev_mu_cap < cost) {
              auto t_old_it = tensordict_of_size[c].find(bits);
              if (t_old_it == tensordict_of_size[c].end() ||
                  cost < t_old_it->second.cost) {
                tensordict_of_size[c][bits] =
                    contract_tf(t1.second, t2.second, bond_dims);
              }
            }
          }
        }
      }
    }

    prev_mu_cap = mu_cap;
    mu_cap = std::max(next_mu_cap, mu_cap * xi_min);
    for (auto& s : tensordict_of_size) {
      for (auto& t : s) {
        t.second.is_new = false;
      }
    }
  } // while

  auto t_final = tensordict_of_size.back().find((mp::cpp_int(1) << n) - 1);
  return {t_final->second.rpn, t_final->second.cost};
}

int test(std::string filename) {
  TensorNetwork tn = TensorNetwork::from_file(filename);

  auto result = tn.optimize();
  std::vector<int> rpn = result.first;
  double cpu = result.second;

  // calc_memory(tn, rpn);

  for (const auto& name : tn.get_tensor_names()) {
    std::cout << name << " ";
  }
  std::cout << std::endl;

  for (const auto& val : rpn) {
    std::cout << val << " ";
  }
  std::cout << std::endl;
  std::cout << cpu << std::endl;

  return 0;
}

}  // namespace find_contraction_path
}  // namespace tenes
