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
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <functional>

#include "iTPS.hpp"

#include "../tensor.hpp"

#include "../Lattice.hpp"
#include "../operator.hpp"
#include "../PEPS_Parameters.hpp"
#include "../printlevel.hpp"
#include "../timer.hpp"
#include "../correlation_length.hpp"

namespace tenes {

template <class ptensor>
std::vector<typename iTPS<ptensor>::transfer_matrix_eigenvalues_type>
iTPS<ptensor>::measure_transfer_matrix_eigenvalues() {
  // res[id][0]: direction
  // res[id][1]: coord
  // res[id][2]: value
  // res[id][3]: eigvals
  std::vector<transfer_matrix_eigenvalues_type> res;

  std::array<std::function<void(ptensor &, ptensor const &, int)>, 2> matvec;

  std::mt19937 gen(peps_parameters.seed * 137 + 31415);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  std::shared_ptr<CorrelationLengthCalculator<ptensor>> clength;
  if (peps_parameters.MeanField_Env) {
    clength = std::make_shared<CorrelationLengthCalculator_mf<ptensor>>(
        lattice, Tn, lambda_tensor);
  } else {
    clength = std::make_shared<CorrelationLengthCalculator_ctm<ptensor>>(
        lattice, Tn, C1, C2, C3, C4, eTl, eTt, eTr, eTb);
  }
  for (int dir = 0; dir < 2; ++dir) {
    int W = dir == 0 ? LY : LX;
    for (int fixed = 0; fixed < W; ++fixed) {
      auto eigvals = clength->eigenvalues(dir, fixed, clength_param, gen);
      res.push_back(std::make_tuple(dir, fixed, eigvals));
    }
  }

  return res;
}

template <class ptensor>
void iTPS<ptensor>::save_correlation_length(
    std::vector<typename iTPS<ptensor>::transfer_matrix_eigenvalues_type> const
        &lambdas) {
  if (mpirank != 0) {
    return;
  }
  std::string filename = outdir + "/correlation_length.dat";
  if (peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "    Save correlation length to " << filename << std::endl;
  }
  std::ofstream ofs(filename.c_str());
  ofs << std::scientific
      << std::setprecision(std::numeric_limits<double>::max_digits10);
  ofs << "#  $1: direction 0: +x, 1: +y\n";
  ofs << "#  $2: y (dir=0) or x (dir=1) coorinates\n";
  ofs << "#  $3: correlation length xi = 1/e_1 \n";
  ofs << "# $4-: eigenvalues e_i = -log|t_i/t_0|\n";
  ofs << "#      where i > 0 and t_i is i-th largest eigenvalue of T\n";
  ofs << std::endl;

  for (const auto &lambda : lambdas) {
    int dir, x;
    std::vector<std::complex<double>> eigvals;
    std::tie(dir, x, eigvals) = lambda;

    int L = 1;
    if (dir == 0) {
      L = LX;
    } else {
      L = 0;
      int x = 0;
      do {
        L += LY;
        x = (x + LX + lattice.skew) % LX;
      } while (x != 0);
    }

    const double e0 = std::abs(eigvals[0]);
    const double e1 = std::abs(eigvals[1]) / e0;
    const double correlation_length = -L / std::log(e1);
    ofs << dir << " " << x << " " << correlation_length;

    for (size_t i = 1; i < eigvals.size(); ++i) {
      const double e = std::abs(eigvals[i]) / e0;
      ofs << " " << -std::log(e) / L;
    }
    ofs << std::endl;
  }
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // end of namespace tenes
