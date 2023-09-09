/*
 * finder.cpp and finder.hpp are translated from the original Python code,
 * Tensordot (https://github.com/smorita/Tensordot) written by Satoshi Morita.
 * The original license is the MIT License.
 */

#include "finder.hpp"

#include <iostream>
#include <cassert>

#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <algorithm>

#include <fstream>
#include <sstream>
#include <iterator>
#include <limits>

// #include <omp.h>

namespace tenes {
namespace find_contraction_path {

bool are_overlap(const TensorFrame& t1, const TensorFrame& t2) {
  return (t1.bits & t2.bits) > 0;
}

// returns the cost of contracting two tensors
// if the tensors are not connected, returns 0.0
double get_contracting_cost(const TensorFrame& t1, const TensorFrame& t2,
                            const std::vector<int>& bond_dims) {
  bool connected = false;
  double cost = 1.0;
  auto it1 = t1.bonds.begin();
  auto it2 = t2.bonds.begin();
  auto end1 = t1.bonds.end();
  auto end2 = t2.bonds.end();
  while (it1 != end1 || it2 != end2) {
    if (it2 == end2) {
      cost *= bond_dims[*it1];
      ++it1;
    } else if (it1 == end1) {
      cost *= bond_dims[*it2];
      ++it2;
    } else {
      int b1 = *it1;
      int b2 = *it2;
      if (b1 < b2) {
        cost *= bond_dims[b1];
        ++it1;
      } else if (b1 > b2) {
        cost *= bond_dims[b2];
        ++it2;
      } else {
        cost *= bond_dims[b1];
        ++it1;
        ++it2;
        connected = true;
      }
    }
  }
  if (!connected) {
    return 0.0;
  }
  cost += t1.cost + t2.cost;
  return cost;
}

TensorFrame contract_tf(const TensorFrame& t1, const TensorFrame& t2,
                        const std::vector<int>& bond_dims) {
  std::vector<int> rpn = t1.rpn;
  rpn.insert(rpn.end(), t2.rpn.begin(), t2.rpn.end());
  rpn.push_back(-1);

  // XOR bits and bonds
  boost::multiprecision::cpp_int bits = t1.bits ^ t2.bits;

  double cost = 1.0;
  std::vector<int> bonds;
  bonds.reserve(t1.bonds.size() + t2.bonds.size());
  auto it1 = t1.bonds.begin();
  auto it2 = t2.bonds.begin();
  auto end1 = t1.bonds.end();
  auto end2 = t2.bonds.end();
  while (it1 != end1 || it2 != end2) {
    if (it2 == end2) {
      bonds.push_back(*it1);
      cost *= bond_dims[*it1];
      ++it1;
    } else if (it1 == end1) {
      bonds.push_back(*it2);
      cost *= bond_dims[*it2];
      ++it2;
    } else {
      int b1 = *it1;
      int b2 = *it2;
      if (b1 < b2) {
        bonds.push_back(b1);
        cost *= bond_dims[b1];
        ++it1;
      } else if (b1 > b2) {
        bonds.push_back(b2);
        cost *= bond_dims[b2];
        ++it2;
      } else {
        cost *= bond_dims[b1];
        ++it1;
        ++it2;
      }
    }
  }
  cost += t1.cost + t2.cost;

  return TensorFrame(rpn, bits, bonds, cost);
}

void TensorNetwork::add_tensor(const std::string& t_name,
                               const std::vector<int>& b_names) {
  std::vector<std::string> b_names_str;
  b_names_str.reserve(b_names.size());
  for (const auto& b : b_names) {
    b_names_str.push_back(std::to_string(b));
  }
  add_tensor(t_name, b_names_str);
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

void TensorNetwork::set_bond_dim(int bond_name, int dim) {
  std::string name = std::to_string(bond_name);
  set_bond_dim(name, dim);
}

void TensorNetwork::set_bond_dim(const std::string& bond_name, int dim) {
  const size_t n = bond_names.size();
  for(size_t i=0; i<n; ++i){
    if(bond_names[i] == bond_name){
      bond_dims[i] = dim;
      break;
    }
  }
}

void TensorNetwork::set_default_bond_dim(int dim) {
  default_bond_dim = dim;
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
    std::vector<int> bonds = t.bonds;
    std::sort(bonds.begin(), bonds.end());
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
      std::vector<size_t> ns(c + 1);
      std::vector<std::vector<std::pair<mp::cpp_int, TensorFrame>>> t_list(c +
                                                                           1);
      for (int d = 1; d <= c; ++d) {
        t_list[d].assign(tensordict_of_size[d].begin(),
                         tensordict_of_size[d].end());
        ns[d] = t_list[d].size();
      }
      // bits -> d1, d2, i1, i2, cost
      std::map<mp::cpp_int, std::tuple<int, int, int, int, double>> winners;
      for (int d1 = 1; d1 <= c / 2; ++d1) {
        int d2 = c - d1;

        // const auto& t1_list = t_list[d1];
        // const auto& t2_list = t_list[d2];
        // std::vector<std::pair<mp::cpp_int, TensorFrame>> t1_list(
        //     tensordict_of_size[d1].begin(), tensordict_of_size[d1].end());
        // std::vector<std::pair<mp::cpp_int, TensorFrame>> t2_list(
        //     tensordict_of_size[d2].begin(), tensordict_of_size[d2].end());
        // const size_t n1 = t_list[d1].size();
        // const size_t n2 = t_list[d2].size();
        const size_t n1 = ns[d1];
        const size_t n2 = ns[d2];
        // const size_t n1 = t1_list.size();
        // const size_t n2 = t2_list.size();
        // const size_t n1 = tensordict_of_size[d1].size();
        // const size_t n2 = tensordict_of_size[d2].size();
        if (n1 * n2 == 0) continue;

        // for (const auto& t1 : tensordict_of_size[d1]) {
        // #pragma omp parallel for
        for (size_t i1 = 0; i1 < n1; ++i1) {
          const auto& t1 = t_list[d1][i1];
          // for (const auto& t2 : tensordict_of_size[d2]) {
          for (size_t i2 = 0; i2 < n2; ++i2) {
            const auto& t2 = t_list[d2][i2];
            if (are_overlap(t1.second, t2.second)) {
              continue;
            }

            double cost = get_contracting_cost(t1.second, t2.second, bond_dims);
            if (cost == 0.0) continue;  // disconnected
            mp::cpp_int bits = t1.first ^ t2.first;

            if (next_mu_cap <= cost) {
              continue;
            } else if (mu_cap < cost) {
              next_mu_cap = cost;
            } else if (t1.second.is_new || t2.second.is_new ||
                       prev_mu_cap < cost) {
              auto winner = winners.find(bits);
              if (winner == winners.end() ||
                  cost < std::get<4>(winner->second)) {
                winners[bits] = std::make_tuple(d1, d2, i1, i2, cost);
              }
              // auto t_old_it = tensordict_of_size[c].find(bits);
              // if (t_old_it == tensordict_of_size[c].end() ||
              //     cost < t_old_it->second.cost) {
              //   tensordict_of_size[c][bits] =
              //   contract_tf(t1.second, t2.second, bond_dims);
              // }
            }
          }  // end for i2
        }    // end for i1
      }      // end for d1
      for (auto const& winner : winners) {
        const int d1 = std::get<0>(winner.second);
        const int d2 = std::get<1>(winner.second);
        const int i1 = std::get<2>(winner.second);
        const int i2 = std::get<3>(winner.second);
        tensordict_of_size[c][winner.first] = contract_tf(
            t_list[d1][i1].second, t_list[d2][i2].second, bond_dims);
      }
    }  // end for c
    prev_mu_cap = mu_cap;
    mu_cap = std::max(next_mu_cap, mu_cap * xi_min);
    for (auto& s : tensordict_of_size) {
      for (auto& t : s) {
        t.second.is_new = false;
      }
    }
  }  // end while

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
