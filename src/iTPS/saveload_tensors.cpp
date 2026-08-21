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

#define _USE_MATH_DEFINES
#include <algorithm>
#include <array>
#include <complex>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <cstdlib>

#include "iTPS.hpp"

#include <mptensor/file_io/file_io.hpp>

#include "../tensor.hpp"

#include "../fermion/fermion_info.hpp"
#include "../fermion/fops.hpp"
#include "../printlevel.hpp"
#include "../util/file.hpp"
#include "../util/string.hpp"

using std::size_t;

namespace tenes::itps {

template <class ptensor>
void iTPS<ptensor>::save_tensors() const {
  std::string const &save_dir = peps_parameters.tensor_save_dir;
  if (save_dir.empty()) {
    return;
  }
  {
    // metadata
    std::string filename = save_dir + "/params.dat";
    std::ofstream ofs(filename.c_str());

    constexpr int tensor_format_version = 1;
    ofs << tensor_format_version << " # Format_Version\n";
    ofs << N_UNIT << " # N_UNIT\n";
    ofs << CHI << " # CHI\n";
    for (int i = 0; i < N_UNIT; ++i) {
      for (int j = 0; j < nleg; ++j) {
        ofs << lattice.virtual_dims[i][j] << " ";
      }
      ofs << lattice.physical_dims[i] << " # Shape of Tn[" << i << "]\n";
    }
  }
  for (int i = 0; i < N_UNIT; ++i) {
    std::string filename = save_dir + "/";
    std::string suffix = "_" + std::to_string(i) + ".dat";
    Tn[i].save((filename + "T" + suffix).c_str());
    eTt[i].save((filename + "Et" + suffix).c_str());
    eTr[i].save((filename + "Er" + suffix).c_str());
    eTb[i].save((filename + "Eb" + suffix).c_str());
    eTl[i].save((filename + "El" + suffix).c_str());
    C1[i].save((filename + "C1" + suffix).c_str());
    C2[i].save((filename + "C2" + suffix).c_str());
    C3[i].save((filename + "C3" + suffix).c_str());
    C4[i].save((filename + "C4" + suffix).c_str());
  }
  if (mpirank == 0) {
    for (int i = 0; i < N_UNIT; ++i) {
      std::ofstream ofs(save_dir + "/lambda_" + std::to_string(i) + ".dat");
      // max_digits10 round-trips a double exactly; the default 6 significant
      // digits silently truncated the Schmidt weights on every checkpoint.
      ofs << std::scientific
          << std::setprecision(std::numeric_limits<double>::max_digits10);
      for (int j = 0; j < nleg; ++j) {
        for (int k = 0; k < lattice.virtual_dims[i][j]; ++k) {
          ofs << lambda_tensor[i][j][k] << std::endl;
        }
      }
    }
  }
  if (finfo.enabled) {
    // The virtual-bond parity ledger is mutable state (the simple update
    // rewrites it through svd_trunc), so it has to travel with the tensors:
    // reloading with a stale ledger changes the measured energy without any
    // error message.
    //
    // The CTM environment tensors saved above carry no information in fermion
    // mode: save_tensors() runs before measure() (main.cpp, run_groundstate)
    // and the fermionic CTM is rebuilt from scratch inside measure().
    save_fermion_parity(save_dir);
  }
  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "Tensors saved in " << save_dir << std::endl;
  }
}

template <class ptensor>
void iTPS<ptensor>::save_fermion_parity(std::string const &save_dir) const {
  if (mpirank != 0) {
    return;
  }
  std::string filename = save_dir + "/fermion.dat";
  std::ofstream ofs(filename.c_str());
  constexpr int fermion_format_version = 1;
  ofs << fermion_format_version << " # Fermion_Format_Version\n";
  ofs << N_UNIT << " # N_UNIT\n";
  ofs << lattice.LX << " " << lattice.LY << " # L_sub\n";
  ofs << lattice.skew << " # skew\n";
  auto write_parity = [&ofs](tenes::fermion::parity_vector const &p) {
    for (std::size_t i = 0; i < p.size(); ++i) {
      ofs << (p[i] ? 1 : 0) << " ";
    }
  };
  for (int i = 0; i < N_UNIT; ++i) {
    write_parity(finfo.phys[i]);
    ofs << "# parity of the physical leg of Tn[" << i << "]\n";
    for (int leg = 0; leg < nleg; ++leg) {
      write_parity(finfo.virt[i][leg]);
      ofs << "# parity of the virtual leg " << leg << " of Tn[" << i << "]\n";
    }
  }
}

template <class ptensor>
void iTPS<ptensor>::load_tensors() {
  std::string const &load_dir = peps_parameters.tensor_load_dir;

  if (!util::isdir(load_dir)) {
    std::string msg = load_dir + " does not exists.";
    throw tenes::load_error(msg);
  }

  int tensor_format_version = 0;
  if (mpirank == 0) {
    std::string filename = load_dir + "/params.dat";
    std::string line;
    if (util::path_exists(filename)) {
      std::ifstream ifs(filename.c_str());
      std::getline(ifs, line);
      tensor_format_version = std::stoi(util::drop_comment(line));
    }
  }
  bcast(tensor_format_version, 0, comm);

  load_fermion_ledger(load_dir);

  if (tensor_format_version == 0) {
    load_tensors_v0();
  } else if (tensor_format_version == 1) {
    load_tensors_v1();
  } else {
    std::stringstream ss;
    ss << "ERROR: Unknown saved tensor format version: "
       << tensor_format_version;
    throw tenes::load_error(ss.str());
  }

  validate_loaded_fermion_tensors();
}

template <class ptensor>
void iTPS<ptensor>::load_fermion_ledger(std::string const &load_dir) {
  const std::string filename = load_dir + "/fermion.dat";
  int exists = 0;
  if (mpirank == 0) {
    exists = util::path_exists(filename) ? 1 : 0;
  }
  bcast(exists, 0, comm);

  if (!finfo.enabled) {
    if (exists != 0) {
      throw tenes::load_error(
          "ERROR: " + filename +
          " exists, i.e. the saved tensors come from a fermionic run, but "
          "parameter.general.fermion is false.\n"
          "HINT: set fermion = true, or load tensors saved by a "
          "non-fermionic run.");
    }
    return;
  }
  if (exists == 0) {
    throw tenes::load_error(
        "ERROR: cannot find " + filename +
        ".\n"
        "The saved tensors do not carry the fermionic parity ledger of the "
        "virtual bonds, which fermion mode needs in order to interpret them.\n"
        "HINT: they were saved by a non-fermionic run, or by a version of "
        "TeNeS that could not save fermionic tensors.");
  }

  std::string content;
  if (mpirank == 0) {
    std::ifstream ifs(filename.c_str());
    std::stringstream ss;
    ss << ifs.rdbuf();
    content = ss.str();
  }
  bcast(content, 0, comm);

  std::istringstream iss(content);
  std::string line;
  auto next_line = [&iss, &line, &filename]() -> std::string {
    if (!std::getline(iss, line)) {
      throw tenes::load_error("ERROR: " + filename + " ends unexpectedly");
    }
    return util::drop_comment(line);
  };
  std::string last_int_line;
  auto next_ints = [&next_line, &filename,
                    &last_int_line]() -> std::vector<int> {
    last_int_line = next_line();
    const auto words = util::split(last_int_line);
    std::vector<int> ret;
    ret.reserve(words.size());
    for (const auto &w : words) {
      try {
        ret.push_back(std::stoi(w));
      } catch (const std::exception &) {
        throw tenes::load_error("ERROR: cannot parse " + filename +
                                " line as integers: \"" + last_int_line + "\"");
      }
    }
    return ret;
  };
  auto next_scalar = [&next_ints, &filename, &last_int_line]() -> int {
    const auto values = next_ints();
    try {
      return values.at(0);
    } catch (const std::exception &) {
      throw tenes::load_error("ERROR: expected an integer in " + filename +
                              " line: \"" + last_int_line + "\"");
    }
  };

  const int version = next_scalar();
  if (version != 1) {
    std::stringstream ss;
    ss << "ERROR: " << filename << " has fermion format version " << version
       << " but this version of TeNeS supports only version 1";
    throw tenes::load_error(ss.str());
  }
  const int loaded_N_UNIT = next_scalar();
  if (loaded_N_UNIT != N_UNIT) {
    std::stringstream ss;
    ss << "ERROR: N_UNIT is " << N_UNIT << " but " << filename << " has "
       << loaded_N_UNIT;
    throw tenes::load_error(ss.str());
  }
  const auto lsub = next_ints();
  const int loaded_skew = next_scalar();
  if (lsub.size() != 2 || lsub[0] != lattice.LX || lsub[1] != lattice.LY ||
      loaded_skew != lattice.skew) {
    std::stringstream ss;
    ss << "ERROR: the unit cell of the saved tensors (L_sub = [";
    for (std::size_t i = 0; i < lsub.size(); ++i) {
      ss << (i == 0 ? "" : ", ") << lsub[i];
    }
    ss << "], skew = " << loaded_skew << ") differs from the input (L_sub = ["
       << lattice.LX << ", " << lattice.LY << "], skew = " << lattice.skew
       << ").\n"
       << "HINT: the parity ledger is indexed by (site, leg), so it only means "
          "the same thing on the same lattice.";
    throw tenes::load_error(ss.str());
  }

  auto to_parity = [](std::vector<int> const &v) {
    tenes::fermion::parity_vector p(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
      p[i] = (v[i] != 0);
    }
    return p;
  };

  std::vector<std::array<tenes::fermion::parity_vector, 4>> virt(N_UNIT);
  for (int i = 0; i < N_UNIT; ++i) {
    const auto phys = to_parity(next_ints());
    if (phys != finfo.phys[i]) {
      std::stringstream ss;
      ss << "ERROR: the physical parity of the tensor " << i << " in "
         << filename << " differs from tensor.unitcell.parity in the input";
      throw tenes::load_error(ss.str());
    }
    for (int leg = 0; leg < nleg; ++leg) {
      auto p = to_parity(next_ints());
      if (static_cast<int>(p.size()) != lattice.virtual_dims[i][leg]) {
        std::stringstream ss;
        ss << "ERROR: the virtual dimension of the leg " << leg
           << " of the tensor " << i << " is " << lattice.virtual_dims[i][leg]
           << " but the saved tensors have " << p.size() << ".\n"
           << "HINT: fermion mode cannot change virtual_dim on restart, "
              "because resizing a leg does not preserve its even/odd blocks. "
              "Keep virtual_dim as it was, or start without tensor_load.";
        throw tenes::load_error(ss.str());
      }
      virt[i][leg] = p;
    }
  }
  for (int i = 0; i < N_UNIT; ++i) {
    finfo.virt[i] = virt[i];
  }
  try {
    tenes::fermion::validate_neighbor_consistency(finfo, lattice);
  } catch (const std::exception &e) {
    throw tenes::load_error(
        std::string(e.what()) + "\nHINT: the parity ledger in " + filename +
        " disagrees between the two ends of a bond, so it does not describe "
        "the lattice being loaded into.");
  }
}

template <class ptensor>
void iTPS<ptensor>::validate_loaded_fermion_tensors() const {
  if (!finfo.enabled) {
    return;
  }
  for (int i = 0; i < N_UNIT; ++i) {
    const auto ft = tenes::fermion::wrap_Tn(Tn[i], finfo, i);
    // parity_violation scans the process-local slice; max_abs may already
    // reduce internally, but the extra allreduce is harmless.
    std::vector<double> reduced{tenes::fermion::parity_violation(ft),
                                tenes::fermion::max_abs(ft)};
    tenes::allreduce_max(reduced, comm);
    const double violation = reduced[0];
    const double scale = std::max(reduced[1], 1.0e-300);
    if (violation > 1.0e-12 * scale) {
      std::stringstream ss;
      ss << "ERROR: the loaded tensor " << i
         << " breaks fermion parity under the loaded parity ledger "
         << "(max violating amplitude " << violation << ", max amplitude "
         << reduced[1] << ").\n"
         << "HINT: the tensors and " << peps_parameters.tensor_load_dir
         << "/fermion.dat do not belong together.";
      throw tenes::load_error(ss.str());
    }
  }
}

template <class ptensor>
void load_tensor(ptensor &A, std::string const &name,
                 std::string const &directory, int iunit) {
  std::string filename =
      directory + "/" + name + "_" + std::to_string(iunit) + ".dat";
  if (!util::path_exists(filename)) {
    throw tenes::load_error("ERROR: cannot find a tensor file: " + filename);
  }
  ptensor temp;
  temp.load(filename.c_str());
  if (A.rank() != temp.rank()) {
    std::stringstream ss;
    ss << "ERROR: rank mismatch in load_tensor: ";
    ss << name << "[" << iunit << "] has " << A.rank() << " legs, but";
    ss << "loaded one has " << temp.rank() << " legs." << std::endl;
    ss << "HINT: check the calculation mode. The number of legs differs "
          "between ground state calculation and finite temperature "
          "calculation.";
    throw tenes::load_error(ss.str());
  }
  A = resize_tensor(temp, A.shape());
}

template <class ptensor>
void iTPS<ptensor>::load_tensors_v1() {
  std::string const &load_dir = peps_parameters.tensor_load_dir;

  int loaded_CHI = 1;
  std::vector<std::vector<int>> loaded_shape(N_UNIT,
                                             std::vector<int>(nleg + 1));
  if (mpirank == 0) {
    std::string filename = load_dir + "/params.dat";
    std::string line;
    std::ifstream ifs(filename.c_str());

    std::getline(ifs, line);
    const int format_version = std::stoi(util::drop_comment(line));
    if (format_version != 1) {
      std::stringstream ss;
      ss << "ERROR: " << filename << " has format version " << format_version
         << " but load_tensors_v1 supports only version 1";
      throw tenes::load_error(ss.str());
    }

    std::getline(ifs, line);
    const int loaded_N_UNIT = std::stoi(util::drop_comment(line));
    if (N_UNIT != loaded_N_UNIT) {
      std::stringstream ss;
      ss << "ERROR: N_UNIT is " << N_UNIT << " but loaded N_UNIT has "
         << loaded_N_UNIT << std::endl;
      throw tenes::load_error(ss.str());
    }

    std::getline(ifs, line);
    loaded_CHI = std::stoi(util::drop_comment(line));
    if (CHI != static_cast<std::size_t>(loaded_CHI)) {
      if (peps_parameters.print_level >= PrintLevel::info) {
        std::cout << "WARNING: parameters.ctm.dimension is " << CHI
                  << " but loaded tensors have CHI = " << loaded_CHI
                  << std::endl;
      }
    }

    for (int i = 0; i < N_UNIT; ++i) {
      std::getline(ifs, line);
      const auto shape = util::split(util::drop_comment(line));
      for (int j = 0; j < nleg; ++j) {
        loaded_shape[i][j] = std::stoi(shape[j]);
        const int vd_param = lattice.virtual_dims[i][j];
        if (vd_param != loaded_shape[i][j]) {
          if (peps_parameters.print_level >= PrintLevel::info) {
            std::cout << "WARNING: virtual dimension of the leg " << j
                      << " of the tensor " << i << " is " << vd_param
                      << " but loaded tensor has " << loaded_shape[i][j]
                      << std::endl;
          }
        }
      }
      loaded_shape[i][nleg] = std::stoi(shape[nleg]);
      const int pdim = lattice.physical_dims[i];
      if (pdim != loaded_shape[i][nleg]) {
        std::stringstream ss;
        ss << "ERROR: dimension of the physical bond of the tensor " << i
           << " is " << pdim << " but loaded tensor has "
           << loaded_shape[i][nleg] << std::endl;
        throw tenes::load_error(ss.str());
      }
    }
  }
  for (int i = 0; i < N_UNIT; ++i) {
    bcast(loaded_shape[i], 0, comm);
  }

  // #define LOAD_TENSOR_(A, name)                      \
  //   do {                                             \
  //     ptensor temp;                                  \
  //     temp.load((filename + name + suffix).c_str()); \
  //     A = resize_tensor(temp, A.shape());            \
  //   } while (false)

  for (int i = 0; i < N_UNIT; ++i) {
    std::string filename = load_dir + "/";
    std::string suffix = "_" + std::to_string(i) + ".dat";

    load_tensor(Tn[i], "T", load_dir, i);
    load_tensor(eTl[i], "El", load_dir, i);
    load_tensor(eTt[i], "Et", load_dir, i);
    load_tensor(eTr[i], "Er", load_dir, i);
    load_tensor(eTb[i], "Eb", load_dir, i);
    load_tensor(C1[i], "C1", load_dir, i);
    load_tensor(C2[i], "C2", load_dir, i);
    load_tensor(C3[i], "C3", load_dir, i);
    load_tensor(C4[i], "C4", load_dir, i);
    // LOAD_TENSOR_(Tn[i], "T");
    // LOAD_TENSOR_(eTl[i], "El");
    // LOAD_TENSOR_(eTt[i], "Et");
    // LOAD_TENSOR_(eTr[i], "Er");
    // LOAD_TENSOR_(eTb[i], "Eb");
    // LOAD_TENSOR_(C1[i], "C1");
    // LOAD_TENSOR_(C2[i], "C2");
    // LOAD_TENSOR_(C3[i], "C3");
    // LOAD_TENSOR_(C4[i], "C4");
  }
  // #undef LOAD_TENSOR_

  std::vector<double> ls;
  if (mpirank == 0) {
    for (int i = 0; i < N_UNIT; ++i) {
      std::string lambda_filename =
          load_dir + "/lambda_" + std::to_string(i) + ".dat";
      std::ifstream ifs(lambda_filename.c_str());
      for (int j = 0; j < nleg; ++j) {
        for (int k = 0; k < loaded_shape[i][j]; ++k) {
          double temp = 0.0;
          if (!(ifs >> temp)) {
            throw tenes::load_error(
                "ERROR: failed to read lambda values from " + lambda_filename);
          }
          ls.push_back(temp);
        }
      }
    }
  }
  bcast(ls, 0, comm);
  int index = 0;
  for (int i = 0; i < N_UNIT; ++i) {
    const auto vdim = lattice.virtual_dims[i];
    for (int j = 0; j < nleg; ++j) {
      lambda_tensor[i][j].clear();
      for (int k = 0; k < loaded_shape[i][j]; ++k) {
        lambda_tensor[i][j].push_back(ls[index]);
        ++index;
      }
      lambda_tensor[i][j].resize(vdim[j]);
    }
  }
}

template <class ptensor>
void iTPS<ptensor>::load_tensors_v0() {
  using mptensor::Shape;
  std::string const &load_dir = peps_parameters.tensor_load_dir;

  // load from the checkpoint
  if (!util::isdir(load_dir)) {
    std::string msg = load_dir + " does not exists.";
    throw tenes::load_error(msg);
  }
  for (int i = 0; i < N_UNIT; ++i) {
    std::string filename = load_dir + "/";
    std::string suffix = "_" + std::to_string(i) + ".dat";
    auto load = [&filename, &suffix](ptensor &A, const char *name) {
      std::string path = filename + name + suffix;
      if (!util::path_exists(path)) {
        throw tenes::load_error("ERROR: cannot find a tensor file: " + path);
      }
      A.load(path.c_str());
    };
    load(Tn[i], "T");
    load(eTt[i], "Et");
    load(eTr[i], "Er");
    load(eTb[i], "Eb");
    load(eTl[i], "El");
    load(C1[i], "C1");
    load(C2[i], "C2");
    load(C3[i], "C3");
    load(C4[i], "C4");
  }
  std::vector<double> ls;
  if (mpirank == 0) {
    for (int i = 0; i < N_UNIT; ++i) {
      const auto vdim = lattice.virtual_dims[i];
      std::string lambda_filename =
          load_dir + "/lambda_" + std::to_string(i) + ".dat";
      std::ifstream ifs(lambda_filename.c_str());
      for (int j = 0; j < nleg; ++j) {
        for (int k = 0; k < vdim[j]; ++k) {
          double temp = 0.0;
          if (!(ifs >> temp)) {
            throw tenes::load_error(
                "ERROR: failed to read lambda values from " + lambda_filename);
          }
          ls.push_back(temp);
        }
      }
    }
  }
  bcast(ls, 0, comm);
  int index = 0;
  for (int i = 0; i < N_UNIT; ++i) {
    const auto vdim = lattice.virtual_dims[i];
    for (int j = 0; j < nleg; ++j) {
      for (int k = 0; k < vdim[j]; ++k) {
        lambda_tensor[i][j][k] = ls[index];
        ++index;
      }
    }
  }

  // overwrite dimensions
  const Shape Cshape = C1[0].shape();
  if (CHI != Cshape[0]) {
    if (peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "WARNING: parameters.ctm.dimension is " << CHI
                << " but loaded tensors have CHI = " << Cshape[0] << std::endl;
    }
  }
  for (int i = 0; i < N_UNIT; ++i) {
    const Shape Tshape = Tn[i].shape();
    const int pdim = lattice.physical_dims[i];
    if (static_cast<std::size_t>(pdim) != Tshape[4]) {
      std::stringstream ss;
      ss << "ERROR: dimension of the physical bond of the tensor " << i
         << " is " << pdim << " but loaded tensor has " << Tshape[4]
         << std::endl;
      throw tenes::input_error(ss.str());
    }

    for (int l = 0; l < nleg; ++l) {
      const int vd_param = lattice.virtual_dims[i][l];
      const int vd_loaded = Tshape[l];
      if (vd_param != vd_loaded) {
        if (peps_parameters.print_level >= PrintLevel::info) {
          std::cout << "WARNING: virtual dimension of the leg " << l
                    << " of the tensor " << i << " is " << vd_param
                    << " but loaded tensor has " << vd_loaded << std::endl;
        }
      }
    }
  }
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace tenes::itps
