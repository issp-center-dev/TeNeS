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
#include <complex>
#include <random>
#include <string>
#include <array>
#include <cstdlib>

#include "iTPS.hpp"

#include "../tensor.hpp"

#include "../printlevel.hpp"
#include "../util/file.hpp"
#include "../util/string.hpp"

using std::size_t;

namespace tenes {
namespace itps {

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

    Tn.push_back(
        ptensor(comm, Shape(vdim[0], vdim[1], vdim[2], vdim[3], pdim)));
    eTt.push_back(ptensor(comm, Shape(CHI, CHI, vdim[1], vdim[1])));
    eTr.push_back(ptensor(comm, Shape(CHI, CHI, vdim[2], vdim[2])));
    eTb.push_back(ptensor(comm, Shape(CHI, CHI, vdim[3], vdim[3])));
    eTl.push_back(ptensor(comm, Shape(CHI, CHI, vdim[0], vdim[0])));
    C1.push_back(ptensor(comm, Shape(CHI, CHI)));
    C2.push_back(ptensor(comm, Shape(CHI, CHI)));
    C3.push_back(ptensor(comm, Shape(CHI, CHI)));
    C4.push_back(ptensor(comm, Shape(CHI, CHI)));

    std::vector<std::vector<double>> lambda(nleg);
    for (int j = 0; j < nleg; ++j) {
      lambda[j] = std::vector<double>(vdim[j], 1.0);
    }
    lambda_tensor.push_back(lambda);

    ptensor id(comm, mptensor::Shape(pdim, pdim));
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
void iTPS<ptensor>::initialize_tensors_density() {
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

    Tn.push_back(
        ptensor(comm, Shape(vdim[0], vdim[1], vdim[2], vdim[3], pdim, pdim)));
    eTt.push_back(ptensor(comm, Shape(CHI, CHI, vdim[1])));
    eTr.push_back(ptensor(comm, Shape(CHI, CHI, vdim[2])));
    eTb.push_back(ptensor(comm, Shape(CHI, CHI, vdim[3])));
    eTl.push_back(ptensor(comm, Shape(CHI, CHI, vdim[0])));
    C1.push_back(ptensor(comm, Shape(CHI, CHI)));
    C2.push_back(ptensor(comm, Shape(CHI, CHI)));
    C3.push_back(ptensor(comm, Shape(CHI, CHI)));
    C4.push_back(ptensor(comm, Shape(CHI, CHI)));

    std::vector<std::vector<double>> lambda(nleg);
    for (int j = 0; j < nleg; ++j) {
      lambda[j] = std::vector<double>(vdim[j], 1.0);
    }
    lambda_tensor.push_back(lambda);

    ptensor id(comm, mptensor::Shape(pdim, pdim));
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

  if (peps_parameters.tensor_load_dir.empty()) {
    mptensor::Index index;
    for (int i = 0; i < lattice.N_UNIT; ++i) {
      const auto pdim = lattice.physical_dims[i];
      const auto vdim = lattice.virtual_dims[i];

      for (int n = 0; n < Tn[i].local_size(); ++n) {
        index = Tn[i].global_index(n);
        if (index[0] == 0 && index[1] == 0 && index[2] == 0 && index[3] == 0 &&
            index[4] == index[5]) {
          Tn[i].set_value(index, to_tensor_type(1.0));
        } else {
          Tn[i].set_value(index, to_tensor_type(0.0));
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

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
