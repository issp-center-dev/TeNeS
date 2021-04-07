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

#include "../printlevel.hpp"
#include "../util/file.hpp"
#include "../exception.hpp"

#include "load_toml.hpp"
#include "iTPS.hpp"

using std::size_t;

namespace {
tenes::Operators<mptensor::Tensor<tenes::mptensor_matrix_type, double>> to_real(
    tenes::Operators<mptensor::Tensor<tenes::mptensor_matrix_type,
                                      std::complex<double>>> const& ops) {
  tenes::Operators<mptensor::Tensor<tenes::mptensor_matrix_type, double>> ret;
  for (auto const& op : ops) {
    if (op.is_onesite()) {
      mptensor::Tensor<tenes::mptensor_matrix_type, double> A(op.op.shape());
      for (size_t lindex = 0; lindex < op.op.local_size(); ++lindex) {
        A[lindex] = op.op[lindex].real();
      }
      ret.emplace_back(op.name, op.group, op.source_site, A);
    } else if (op.ops_indices.empty()) {
      mptensor::Tensor<tenes::mptensor_matrix_type, double> A(op.op.shape());
      for (size_t lindex = 0; lindex < op.op.local_size(); ++lindex) {
        A[lindex] = op.op[lindex].real();
      }
      ret.emplace_back(op.name, op.group, op.source_site, op.dx, op.dy, A);
    } else {
      ret.emplace_back(op.name, op.group, op.source_site, op.dx, op.dy,
                       op.ops_indices);
    }
  }
  return ret;
}

bool is_real(tenes::Operators<mptensor::Tensor<
                 tenes::mptensor_matrix_type, std::complex<double>>> const& ops,
             double tol) {
  for (auto const& op : ops) {
    if (!op.ops_indices.empty()) {
      continue;
    }
    int res = 0;
    for (size_t lindex = 0; lindex < op.op.local_size(); ++lindex) {
      if (std::abs(op.op[lindex].imag()) > tol) {
        res = 1;
        break;
      }
    }
    tenes::allreduce_sum(res, op.op.get_comm());
    if (res > 0) {
      return false;
    }
  }
  return true;
}

tenes::NNOperators<mptensor::Tensor<tenes::mptensor_matrix_type, double>>
to_real(tenes::NNOperators<mptensor::Tensor<tenes::mptensor_matrix_type,
                                            std::complex<double>>> const& ops) {
  tenes::NNOperators<mptensor::Tensor<tenes::mptensor_matrix_type, double>> ret;
  for (auto const& op : ops) {
    mptensor::Tensor<tenes::mptensor_matrix_type, double> A(op.op.shape());
    for (size_t lindex = 0; lindex < op.op.local_size(); ++lindex) {
      A[lindex] = op.op[lindex].real();
    }
    ret.emplace_back(op.source_site, op.source_leg, A);
  }
  return ret;
}

bool is_real(tenes::NNOperators<mptensor::Tensor<
                 tenes::mptensor_matrix_type, std::complex<double>>> const& ops,
             double tol) {
  for (auto const& op : ops) {
    int res = 0;
    for (size_t lindex = 0; lindex < op.op.local_size(); ++lindex) {
      if (std::abs(op.op[lindex].imag()) > tol) {
        res = 1;
        break;
      }
    }
    tenes::allreduce_sum(res, op.op.get_comm());
    if (res > 0) {
      return false;
    }
  }
  return true;
}

}  // end of unnamed namespace

namespace tenes {

template <class tensor>
int tenes(MPI_Comm comm, PEPS_Parameters peps_parameters, SquareLattice lattice,
          NNOperators<tensor> simple_updates, NNOperators<tensor> full_updates,
          Operators<tensor> onesite_operators,
          Operators<tensor> twosite_operators, CorrelationParameter corparam,
          TransferMatrix_Parameters clength_param) {
  iTPS<tensor> tns(comm, peps_parameters, lattice, simple_updates, full_updates,
                   onesite_operators, twosite_operators, corparam,
                   clength_param);
  tns.optimize();
  tns.save_tensors();
  if (peps_parameters.to_measure) {
    tns.measure();
  }
  tns.summary();
  return 0;
}

int main_impl(std::string input_filename, MPI_Comm com,
              PrintLevel print_level = PrintLevel::info) {
  using tensor_complex =
      mptensor::Tensor<mptensor_matrix_type, std::complex<double>>;

  int mpisize = 0, mpirank = 0;
  MPI_Comm_rank(com, &mpirank);
  MPI_Comm_size(com, &mpisize);

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
  peps_parameters.Bcast(MPI_COMM_WORLD);

  auto toml_lattice = input_toml->get_table("tensor");
  if (toml_lattice == nullptr) {
    throw tenes::input_error("[tensor] not found");
  }
  SquareLattice lattice = gen_lattice(toml_lattice);
  lattice.Bcast(MPI_COMM_WORLD);

  // time evolution
  auto toml_evolution = input_toml->get_table("evolution");
  if (toml_evolution == nullptr) {
    throw tenes::input_error("[evolution] not found");
  }

  const double tol = peps_parameters.iszero_tol;

  const auto simple_updates = load_simple_updates<tensor_complex>(input_toml);
  const auto full_updates = load_full_updates<tensor_complex>(input_toml);

  // observable
  auto toml_observable = input_toml->get_table("observable");
  if (toml_observable == nullptr) {
    throw tenes::input_error("[observable] not found");
  }

  // onesite observable
  const auto onesite_obs = load_operators<tensor_complex>(
      input_toml, lattice.N_UNIT, 1, tol, "observable.onesite");
  const auto twosite_obs = load_operators<tensor_complex>(
      input_toml, lattice.N_UNIT, 2, tol, "observable.twosite");

  // correlation
  auto toml_correlation = input_toml->get_table("correlation");
  const auto corparam = (toml_correlation != nullptr
                             ? gen_corparam(toml_correlation, "correlation")
                             : CorrelationParameter());

  // correlation length
  auto toml_clength = input_toml->get_table("correlation_length");
  const auto clength_param =
      (toml_clength != nullptr
           ? gen_transfer_matrix_parameter(toml_clength, "correlation_length")
           : TransferMatrix_Parameters());

  bool is_real = peps_parameters.is_real;
  is_real = is_real && ::is_real(simple_updates, tol);
  is_real = is_real && ::is_real(full_updates, tol);
  is_real = is_real && ::is_real(onesite_obs, tol);
  is_real = is_real && ::is_real(twosite_obs, tol);

  if (peps_parameters.is_real && !is_real) {
    std::stringstream ss;
    ss << "TeNeS invoked in real tensor mode (parameter.general.is_real = "
          "true) but some operators are complex.\n";
    ss << "Consider using larger parameter.general.iszero_tol (present: " << tol
       << ")";
    throw tenes::input_error(ss.str());
  }

  std::string outdir = peps_parameters.outdir;
  bool is_ok = true;
  if (mpirank == 0) {
    if (!util::isdir(outdir)) {
      is_ok = util::mkdir(outdir);
    }
  }
  bcast(is_ok, 0, com);
  if (!is_ok) {
    std::stringstream ss;
    ss << "Cannot mkdir " << outdir;
    throw tenes::runtime_error(ss.str());
  }

  if (mpirank == 0) {
    std::string basename = util::basename(input_filename);
    std::ifstream ifs(input_filename.c_str());
    std::string dst_filename = outdir + "/" + basename;
    std::ofstream ofs(dst_filename.c_str());
    std::string line;
    while (std::getline(ifs, line)) {
      ofs << line << std::endl;
    }
  }

  if (is_real) {
    return tenes(MPI_COMM_WORLD, peps_parameters, lattice,
                 to_real(simple_updates), to_real(full_updates),
                 to_real(onesite_obs), to_real(twosite_obs), corparam,
                 clength_param);
  } else {
    return tenes(MPI_COMM_WORLD, peps_parameters, lattice, simple_updates,
                 full_updates, onesite_obs, twosite_obs, corparam,
                 clength_param);
  }
}

}  // end of namespace tenes
