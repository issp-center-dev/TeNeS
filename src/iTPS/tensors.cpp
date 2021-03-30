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
#include <sys/stat.h>
#include <algorithm>
#include <complex>
#include <random>
#include <string>
#include <array>

#include "iTPS.hpp"

#include "../tensor.hpp"

#include "../operator.hpp"
#include "../PEPS_Parameters.hpp"
#include "../printlevel.hpp"
#include "../timer.hpp"
#include "../util/file.hpp"
#include "../util/string.hpp"

namespace tenes {

template <class ptensor>
void iTPS<ptensor>::initialize_tensors() {
  using mptensor::Shape;

  Tn.clear();
  eTt.clear();
  eTr.clear();
  eTb.clear();
  eTl.clear();
  C1.clear();
  C2.clear();
  C3.clear();
  C4.clear();
  lambda_tensor.clear();

  for (int i = 0; i < N_UNIT; ++i) {
    const auto pdim = lattice.physical_dims[i];
    const auto vdim = lattice.virtual_dims[i];

    Tn.push_back(ptensor(Shape(vdim[0], vdim[1], vdim[2], vdim[3], pdim)));
    eTt.push_back(ptensor(Shape(CHI, CHI, vdim[1], vdim[1])));
    eTr.push_back(ptensor(Shape(CHI, CHI, vdim[2], vdim[2])));
    eTb.push_back(ptensor(Shape(CHI, CHI, vdim[3], vdim[3])));
    eTl.push_back(ptensor(Shape(CHI, CHI, vdim[0], vdim[0])));
    C1.push_back(ptensor(Shape(CHI, CHI)));
    C2.push_back(ptensor(Shape(CHI, CHI)));
    C3.push_back(ptensor(Shape(CHI, CHI)));
    C4.push_back(ptensor(Shape(CHI, CHI)));

    std::vector<std::vector<double>> lambda(nleg);
    for (int j = 0; j < nleg; ++j) {
      lambda[j] = std::vector<double>(vdim[j], 1.0);
    }
    lambda_tensor.push_back(lambda);

    ptensor id(mptensor::Shape(pdim, pdim));
    for (int j = 0; j < pdim; ++j) {
      for (int k = 0; k < pdim; ++k) {
        id.set_value(mptensor::Index(j, k), (j == k ? 1.0 : 0.0));
      }
    }
    op_identity.push_back(id);
  }

  std::mt19937 gen(peps_parameters.seed);
  // use another rng for backward compatibility
  std::mt19937 gen_im(peps_parameters.seed * 11 + 137);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  int nr;

  if (peps_parameters.tensor_load_dir.empty()) {
    mptensor::Index index;
    for (int i = 0; i < lattice.N_UNIT; ++i) {
      const auto pdim = lattice.physical_dims[i];
      const auto vdim = lattice.virtual_dims[i];

      const size_t ndim = vdim[0] * vdim[1] * vdim[2] * vdim[3] * pdim;
      std::vector<double> ran_re(ndim);
      std::vector<double> ran_im(ndim);

      for (int j = 0; j < ndim; j++) {
        ran_re[j] = dist(gen);
        ran_im[j] = dist(gen_im);
      }
      auto &dir = lattice.initial_dirs[i];
      std::vector<double> dir_im(pdim);
      if (std::all_of(dir.begin(), dir.end(),
                      [=](double x) { return x == 0.0; })) {
        // random
        dir.resize(pdim);
        for (int j = 0; j < pdim; ++j) {
          dir[j] = dist(gen);
          dir_im[j] = dist(gen_im);
        }
      }

      for (int n = 0; n < Tn[i].local_size(); ++n) {
        index = Tn[i].global_index(n);
        if (index[0] == 0 && index[1] == 0 && index[2] == 0 && index[3] == 0) {
          auto v = std::complex<double>(dir[index[4]], dir_im[index[4]]);
          Tn[i].set_value(index, to_tensor_type(v));
        } else {
          nr = index[0] + index[1] * vdim[0] + index[2] * vdim[0] * vdim[1] +
               index[3] * vdim[0] * vdim[1] * vdim[2] +
               index[4] * vdim[0] * vdim[1] * vdim[2] * vdim[3];
          auto v =
              lattice.noises[i] * std::complex<double>(ran_re[nr], ran_im[nr]);
          Tn[i].set_value(index, to_tensor_type(v));
        }
      }
    }
  } else {
    load_tensors();
    if (peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "Tensors loaded from " << peps_parameters.tensor_load_dir
                << std::endl;
    }
  }  // end of else part of if(load_dir.empty())
}

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
      for (int j = 0; j < nleg; ++j) {
        for (int k = 0; k < lattice.virtual_dims[i][j]; ++k) {
          ofs << lambda_tensor[i][j][k] << std::endl;
        }
      }
    }
  }
  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "Tensors saved in " << save_dir << std::endl;
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
  if (tensor_format_version == 0) {
    load_tensors_v0();
  } else if (tensor_format_version == 1) {
    load_tensors_v1();
  } else {
    std::stringstream ss;
    ss << "ERROR: Unknown checkpoint format version: " << tensor_format_version;
    throw tenes::load_error(ss.str());
  }
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
    if (CHI != loaded_CHI) {
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

#define LOAD_TENSOR_(A, name)                      \
  do {                                             \
    ptensor temp;                                  \
    temp.load((filename + name + suffix).c_str()); \
    A = resize_tensor(temp, A.shape());            \
  } while (false)

  for (int i = 0; i < N_UNIT; ++i) {
    std::string filename = load_dir + "/";
    std::string suffix = "_" + std::to_string(i) + ".dat";

    LOAD_TENSOR_(Tn[i], "T");
    LOAD_TENSOR_(eTl[i], "El");
    LOAD_TENSOR_(eTt[i], "Et");
    LOAD_TENSOR_(eTr[i], "Er");
    LOAD_TENSOR_(eTb[i], "Eb");
    LOAD_TENSOR_(C1[i], "C1");
    LOAD_TENSOR_(C2[i], "C2");
    LOAD_TENSOR_(C3[i], "C3");
    LOAD_TENSOR_(C4[i], "C4");
  }
#undef LOAD_TENSOR_

  std::vector<double> ls;
  if (mpirank == 0) {
    for (int i = 0; i < N_UNIT; ++i) {
      std::ifstream ifs(load_dir + "/lambda_" + std::to_string(i) + ".dat");
      for (int j = 0; j < nleg; ++j) {
        for (int k = 0; k < loaded_shape[i][j]; ++k) {
          double temp = 0.0;
          ifs >> temp;
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
    Tn[i].load((filename + "T" + suffix).c_str());
    eTt[i].load((filename + "Et" + suffix).c_str());
    eTr[i].load((filename + "Er" + suffix).c_str());
    eTb[i].load((filename + "Eb" + suffix).c_str());
    eTl[i].load((filename + "El" + suffix).c_str());
    C1[i].load((filename + "C1" + suffix).c_str());
    C2[i].load((filename + "C2" + suffix).c_str());
    C3[i].load((filename + "C3" + suffix).c_str());
    C4[i].load((filename + "C4" + suffix).c_str());
  }
  std::vector<double> ls;
  if (mpirank == 0) {
    for (int i = 0; i < N_UNIT; ++i) {
      const auto vdim = lattice.virtual_dims[i];
      std::ifstream ifs(load_dir + "/lambda_" + std::to_string(i) + ".dat");
      for (int j = 0; j < nleg; ++j) {
        for (int k = 0; k < vdim[j]; ++k) {
          double temp = 0.0;
          ifs >> temp;
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
    if (pdim != Tshape[4]) {
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

}  // end of namespace tenes
