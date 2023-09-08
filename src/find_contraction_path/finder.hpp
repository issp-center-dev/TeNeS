#ifndef TENES_SRC_FINDCONTRACTIONPATH_FINDER_HPP_
#define TENES_SRC_FINDCONTRACTIONPATH_FINDER_HPP_

#include <cassert>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>

#include <boost/multiprecision/cpp_int.hpp>

namespace tenes {
namespace find_contraction_path {

struct Tensor {
  std::vector<int> name;
  std::vector<int> bonds;

  Tensor() : name(), bonds() {}

  Tensor(int name_value) : name(1, name_value), bonds() {}

  Tensor(const std::vector<int>& name_values) : name(name_values), bonds() {}

  Tensor(const int name_value, const std::vector<int>& bond_values)
      : name(1, name_value), bonds(bond_values) {}

  Tensor(const std::vector<int>& name_values,
         const std::vector<int>& bond_values)
      : name(name_values), bonds(bond_values) {}
};

struct Bond {
  int t0;
  int t1;

  Bond() : t0(-1), t1(-1) {}

  Bond(int t0_value, int t1_value) : t0(t0_value), t1(t1_value) {}

  bool isFree() const { return t0 < 0 || t1 < 0; }

  void connect(int tensor_index) {
    assert(isFree() && "edge already connected to two tensors");
    if (t0 < 0) {
      t0 = tensor_index;
    } else {
      assert(t0 != tensor_index && "edge connects to the same tensor");
      t1 = tensor_index;
    }
  }
};

struct TensorFrame {
  std::vector<int> rpn;
  boost::multiprecision::cpp_int bits;
  std::vector<int> bonds;
  double cost;
  bool is_new;

  TensorFrame(const std::vector<int>& rpn_ = {}, boost::multiprecision::cpp_int bits_ = 0,
              const std::vector<int>& bonds_ = {}, double cost_ = 0.0,
              bool is_new_ = true)
      : rpn(rpn_), bits(bits_), bonds(bonds_), cost(cost_), is_new(is_new_) {}
};

double get_contracting_cost(const TensorFrame& t1, const TensorFrame& t2,
                            const std::vector<int>& bond_dims);
TensorFrame contract(const TensorFrame& t1, const TensorFrame& t2);

class TensorNetwork {
 private:
  std::vector<Tensor> tensors;
  std::vector<Bond> bonds;
  std::vector<std::string> tensor_names;
  std::vector<std::string> bond_names;
  std::vector<int> bond_dims;
  double total_memory;
  double max_memory;
  double cpu_cost;
  int default_bond_dim;
  std::vector<int> result;

 public:
  TensorNetwork()
      : tensors(),
        bonds(),
        tensor_names(),
        bond_names(),
        bond_dims(),
        total_memory(0.0),
        max_memory(0.0),
        cpu_cost(0.0),
        default_bond_dim(10) {}

  std::vector<std::string> const& get_tensor_names() const {
    return tensor_names;
  }

  void add_tensor(const std::string& t_name,
                  const std::vector<std::string>& b_names);

  static TensorNetwork from_file(const std::string& filename);

  TensorNetwork clone() const;

  std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> find_bonds(
      int tensor_a, int tensor_b) const;

  TensorNetwork contract(int t0, int t1, const std::vector<int>& bc,
                         const std::vector<int>& br0,
                         const std::vector<int>& br1);
  double calc_memory(const std::vector<int>& rpn) const;
  void set_bond_dim(const std::string& bond_name, int dim);
  void set_default_bond_dim(int dim);

  std::vector<std::map<boost::multiprecision::cpp_int, TensorFrame>>
  init_tensordict_of_size();
  std::pair<std::vector<int>, double> optimize();
  std::vector<int> path() const { return result; }
};

}  // namespace find_contraction_path
}  // namespace tenes

#endif  // TENES_SRC_FINDCONTRACTIONPATH_FINDER_HPP_
