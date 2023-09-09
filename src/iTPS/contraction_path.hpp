#ifndef TENES_SRC_ITPS_CONTRACTION_PATH_HPP_
#define TENES_SRC_ITPS_CONTRACTION_PATH_HPP_

#include <mptensor/tensor.hpp>
#include <map>

namespace tenes {
namespace itps {

template <class tensor>
class TensorNetworkContractor {
 private:
  int nrows;
  int ncols;
  bool is_tpo;
  bool is_mf;
  std::vector<int> path;

  struct TensorName {
    std::string name;
    std::vector<int> bonds;
    TensorName() : name(), bonds() {}
    TensorName(int bond_size) : name(), bonds(bond_size) {}
  };

  std::vector<TensorName> C_name;
  std::vector<TensorName> eTt_name;
  std::vector<TensorName> eTr_name;
  std::vector<TensorName> eTb_name;
  std::vector<TensorName> eTl_name;
  std::vector<TensorName> Tn_name;

  void initialize_ctm();
  void initialize_mf();

  bool check_path(std::vector<int> const& path) const;

 public:
  TensorNetworkContractor(){};
  TensorNetworkContractor(int nrows, int ncols, bool is_tpo, bool is_mf);

  void initialize(int nrows, int ncols, bool is_tpo, bool is_mf);

  std::pair<std::vector<int>, double> optimize(
      std::vector<int> const& shape_types,
      std::vector<std::vector<int>> const& shape_dims, int chi);
  void load_path(std::vector<int> const& path);
  std::vector<int> get_path() const;

  typename tensor::value_type contract(
      std::vector<tensor const*> const& C,
      std::vector<tensor const*> const& eTt,
      std::vector<tensor const*> const& eTr,
      std::vector<tensor const*> const& eTb,
      std::vector<tensor const*> const& eTl,
      std::vector<std::vector<tensor const*>> const& Tn) const;

  typename tensor::value_type contract(
      std::vector<tensor const*> const& C,
      std::vector<tensor const*> const& eTt,
      std::vector<tensor const*> const& eTr,
      std::vector<tensor const*> const& eTb,
      std::vector<tensor const*> const& eTl,
      std::vector<std::vector<tensor const*>> const& Tn,
      std::map<std::tuple<int, int>, tensor const*> const& ops) const;
};

std::vector<int> default_path(int nrows, int ncols, bool is_tpo, bool is_mf);

}  // namespace itps
}  // namespace tenes
#endif  // TENES_SRC_ITPS_CONTRACTION_PATH_HPP_
