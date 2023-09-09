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

#include <complex>  // for complex
#include <cstdlib>  // for abs, size_t
#include <set>

#include "../version.hpp"
#include "../printlevel.hpp"
#include "../util/string.hpp"
#include "../util/file.hpp"
#include "../exception.hpp"
#include "../mpi.hpp"

#include "load_toml.hpp"
#include "contraction_path.hpp"

using std::size_t;

namespace tenes {
namespace itps {

struct entry {
  int nrow;
  int ncol;
  std::vector<int> shape_types;
  double cost;
  std::vector<int> path;

  void write_toml(std::ofstream& ofs) const {
    ofs << "[[contraction]]" << std::endl;
    ofs << "nrow = " << nrow << std::endl;
    ofs << "ncol = " << ncol << std::endl;
    ofs << "shape_types = [";
    for (const auto& type : shape_types) {
      ofs << type << ", ";
    }
    ofs << "]" << std::endl;
    ofs << "cost = " << cost << std::endl;
    ofs << "contraction_path = [" << std::endl;
    for (const auto& p : path) {
      ofs << p << ", ";
    }
    ofs << std::endl;
    ofs << "]" << std::endl;
    ofs << std::endl;
  }

  static entry recv(int src, MPI_Comm comm) {
    entry e;
    MPI_Status status;
    double buffer_double[2];
    MPI_Recv(buffer_double, 2, MPI_DOUBLE, src, src, comm, &status);
    e.cost = buffer_double[0];

    int size = buffer_double[1];
    std::vector<int> buffer(size);
    MPI_Recv(buffer.data(), size, MPI_INT, src, src, comm, &status);
    e.nrow = buffer[0];
    e.ncol = buffer[1];
    e.shape_types.assign(buffer.begin() + 2,
                         buffer.begin() + 2 + e.nrow * e.ncol);
    e.path.assign(buffer.begin() + 2 + e.nrow * e.ncol, buffer.end());

    return e;
  }

  void send(int tgt, MPI_Comm comm) {
    int mpirank = 0;
    MPI_Comm_rank(comm, &mpirank);
    double buffer_double[2];
    buffer_double[0] = cost;
    buffer_double[1] = shape_types.size() + path.size() + 2;
    MPI_Send(buffer_double, 2, MPI_DOUBLE, tgt, mpirank, comm);

    std::vector<int> buffer;
    buffer.push_back(nrow);
    buffer.push_back(ncol);
    buffer.insert(buffer.end(), shape_types.begin(), shape_types.end());
    buffer.insert(buffer.end(), path.begin(), path.end());
    MPI_Send(buffer.data(), buffer.size(), MPI_INT, tgt, mpirank, comm);
  }
};

template <class tensor>
void make_key(Operator<tensor> const op,
              std::set<std::tuple<int, int, std::vector<int>>>& contractions,
              std::vector<int>& tensor_shape_types,
              const SquareLattice& lattice) {
  const int source = op.source_site;
  int mindx = 0, maxdx = 0, mindy = 0, maxdy = 0;
  for (auto dx : op.dx) {
    mindx = std::min(mindx, dx);
    maxdx = std::max(maxdx, dx);
  }
  for (auto dy : op.dy) {
    mindy = std::min(mindy, dy);
    maxdy = std::max(maxdy, dy);
  }
  const int ncol = maxdx - mindx + 1;
  const int nrow = maxdy - mindy + 1;
  const int source_col = -mindx;
  const int source_row = maxdy;

  std::vector<int> shape_types(nrow * ncol);
  for (int row = 0; row < nrow; ++row) {
    for (int col = 0; col < ncol; ++col) {
      const int index =
          lattice.other(source, col - source_col, source_row - row);
      shape_types[row * ncol + col] = tensor_shape_types[index];
    }
  }
  contractions.insert(std::forward_as_tuple(nrow, ncol, shape_types));
}

int main_make_contraction(std::string input_filename, MPI_Comm comm,
                          PrintLevel print_level) {
  using tensor_complex =
      mptensor::Tensor<mptensor_matrix_type, std::complex<double>>;

  int mpisize = 0, mpirank = 0;
  MPI_Comm_rank(comm, &mpirank);
  MPI_Comm_size(comm, &mpisize);

  if (!util::path_exists(input_filename)) {
    std::stringstream ss;
    ss << "ERROR: cannot find the input file: " << input_filename << std::endl;
    throw tenes::input_error(ss.str());
  }

  auto input_toml = cpptoml::parse_file(input_filename);

  // Parameters
  auto toml_param = input_toml->get_table("parameter");
  PEPS_Parameters peps_parameters =
      (toml_param != nullptr ? gen_param(toml_param) : PEPS_Parameters());
  peps_parameters.print_level = print_level;
  peps_parameters.Bcast(comm);

  std::string output_filename = peps_parameters.contraction_path_file;
  if (output_filename.empty()) {
    throw tenes::input_error(
        "parameter.contraction.pathfile is not specified (this is used for "
        "output filename)");
  }

  const int CHI = peps_parameters.CHI;
  const bool is_TPO = peps_parameters.calcmode ==
                      PEPS_Parameters::CalculationMode::finite_temperature;
  const bool is_mf = peps_parameters.MeanField_Env;

  auto toml_lattice = input_toml->get_table("tensor");
  if (toml_lattice == nullptr) {
    throw tenes::input_error("[tensor] not found");
  }
  SquareLattice lattice = gen_lattice(toml_lattice);
  lattice.Bcast(comm);

  std::map<std::vector<int>, int> shape_type_map;
  std::vector<int> tensor_shape_types(lattice.N_UNIT);
  std::vector<std::vector<int>> tensor_shape_dims;
  for (int i = 0; i < lattice.N_UNIT; ++i) {
    std::vector<int> shape(lattice.virtual_dims[i].begin(),
                           lattice.virtual_dims[i].end());
    shape.push_back(lattice.physical_dims[i]);
    if (is_TPO) {
      shape.push_back(lattice.physical_dims[i]);
    }
    auto it = shape_type_map.find(shape);
    if (it == shape_type_map.end()) {
      shape_type_map[shape] = tensor_shape_types[i] = shape_type_map.size();
      tensor_shape_dims.push_back(shape);
    } else {
      tensor_shape_types[i] = it->second;
    }
  }

  // observable
  auto toml_observable = input_toml->get_table("observable");
  if (toml_observable == nullptr) {
    throw tenes::input_error("[observable] not found");
  }

  const double tol = 0.0;
  // const auto onesite_obs = load_operators<tensor_complex>(
  //     input_toml, comm, lattice.N_UNIT, 1, tol, "observable.onesite");
  const auto twosite_obs = load_operators<tensor_complex>(
      input_toml, comm, lattice.N_UNIT, 2, tol, "observable.twosite");
  const auto multisite_obs = load_operators<tensor_complex>(
      input_toml, comm, lattice.N_UNIT, 3, tol, "observable.multisite");

  using keytype = std::tuple<int, int, std::vector<int>>;

  std::set<keytype> keys_set;
  for (const auto& obs : twosite_obs) {
    make_key(obs, keys_set, tensor_shape_types, lattice);
  }
  for (const auto& obs : multisite_obs) {
    make_key(obs, keys_set, tensor_shape_types, lattice);
  }

  std::vector<keytype> keys(keys_set.begin(), keys_set.end());
  const size_t nkeys = keys.size();
  std::vector<int> indices(nkeys);
  for (int i = 0; i < nkeys; ++i) {
    indices[i] = i;
  }
  std::sort(indices.begin(), indices.end(), [&](int i_a, int i_b) {
    const auto& a = keys[i_a];
    const auto& b = keys[i_b];
    const int ns_a = (std::get<0>(a) + 2) * (std::get<1>(a) + 2);
    const int ns_b = (std::get<0>(b) + 2) * (std::get<1>(b) + 2);
    return ns_a >= ns_b;
  });
  std::vector<int> where(nkeys);
  for (int i = 0; i < nkeys; i++) {
    where[indices[i]] = i;
  }

  std::vector<int> who(nkeys, 0);
  if (mpisize > 1) {
    int r = 0;
    int dr = 1;
    for (int i = 0; i < nkeys; i++) {
      who[i] = r;
      if (r + dr < 0) {
        dr = 1;
      } else if (r + dr >= mpisize) {
        dr = -1;
      }
      r += dr;
    }
  }

  std::vector<entry> tasks(nkeys);
  for (int i = 0; i < nkeys; ++i) {
    if (who[i] == mpirank) {
      const auto& key = keys[indices[i]];
      const int nrow = std::get<0>(key);
      const int ncol = std::get<1>(key);
      const auto& shape_types = std::get<2>(key);
      TensorNetworkContractor<tensor_complex> tnc(nrow, ncol, is_TPO, is_mf);
      auto res = tnc.optimize(shape_types, tensor_shape_dims, CHI);
      tasks[i] = {nrow, ncol, shape_types, res.second, res.first};
      std::cout << "task " << i << " done on rank " << mpirank << std::endl;
    }
  }

  std::vector<entry> results(nkeys);
  if (mpisize == 1) {
    for (int i = 0; i < nkeys; ++i) {
      // results[i] = tasks[where[i]];
      results[indices[i]] = tasks[i];
    }
  } else {
    for (int i = 0; i < nkeys; ++i) {
      int r = who[i];
      if (mpirank == 0) {
        if (r == 0) {
          results[indices[i]] = tasks[i];
        } else {
          results[indices[i]] = entry::recv(r, comm);
        }
      } else {
        if (r == mpirank) {
          tasks[i].send(0, comm);
        }
      }
    }
  }

  if (mpirank == 0) {
    std::ofstream ofs(output_filename.c_str());
    ofs << "[params]" << std::endl;
    ofs << "chi = " << CHI << std::endl;
    ofs << "is_TPO = " << (is_TPO ? "true" : "false") << std::endl;
    ofs << "is_mf = " << (is_mf ? "true" : "false") << std::endl;
    ofs << std::endl;
    ofs << "[tensor]" << std::endl;
    ofs << "L_sub = [" << lattice.LX << ", " << lattice.LY << "]" << std::endl;
    ofs << "skew = " << lattice.skew << std::endl;
    ofs << "shape_types = [";
    for (auto type : tensor_shape_types) {
      ofs << type << ", ";
    }
    ofs << "]" << std::endl;
    ofs << "shape_dims = [" << std::endl;
    for (const auto& dims : tensor_shape_dims) {
      ofs << "[";
      for (const auto& dim : dims) {
        ofs << dim << ", ";
      }
      ofs << "]," << std::endl;
    }
    ofs << "]" << std::endl;
    ofs << std::endl;

    for (const auto& e : results) {
      e.write_toml(ofs);
    }
  }
  MPI_Barrier(comm);
  return 0;
}

}  // namespace itps
}  // namespace tenes

int main_impl(int argc, char** argv) {
  int mpisize = 1;
  int mpirank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int status = 0;

  std::string usage =
      R"(TeNeS: TEnsor NEtwork Solver for 2D quantum lattice system
  
  Usage:
    tenes [--quiet] <input_toml>
    tenes --help
    tenes --version

  Options:
    -h --help       Show this help message.
    -v --version    Show the version.
    -q --quiet      Do not print any messages.
  )";

  using PrintLevel = tenes::PrintLevel;

  PrintLevel print_level = PrintLevel::info;
  std::string input_filename;
  for (int i = 1; i < argc; ++i) {
    std::string opt = argv[i];
    if (opt == "-h" || opt == "--help") {
      if (mpirank == 0) std::cout << usage << std::endl;
      return 0;
    } else if (opt == "-v" || opt == "--version") {
      if (mpirank == 0) std::cout << "TeNeS v" << TENES_VERSION << std::endl;
      return 0;
    } else if (opt == "-q" || opt == "--quiet") {
      print_level = PrintLevel::none;
    } else {
      input_filename = opt;
    }
  }

  if (input_filename.empty()) {
    if (mpirank == 0) std::cout << usage << std::endl;
    return 0;
  }

  try {
    status = tenes::itps::main_make_contraction(input_filename, MPI_COMM_WORLD,
                                                print_level);
  } catch (const tenes::input_error e) {
    if (mpirank == 0) {
      std::cerr << "[INPUT ERROR]" << std::endl;
      std::cerr << e.what() << std::endl;
    }
    status = 1;
  } catch (const tenes::load_error e) {
    if (mpirank == 0) {
      std::cerr << "[TENSOR LOAD ERROR]" << std::endl;
      std::cerr << e.what() << std::endl;
    }
    status = 1;
  } catch (const tenes::runtime_error e) {
    if (mpirank == 0) {
      std::cerr << "[ERROR]" << std::endl;
      std::cerr << e.what() << std::endl;
    }
    status = 1;
  }

  return status;
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int status = main_impl(argc, argv);
  MPI_Finalize();
  return status;
}
