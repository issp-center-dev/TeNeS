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

#include "SquareLattice.hpp"

#include <fstream>      // for operator<<, basic_ostream, char_traits, endl
#include <tuple>        // for tuple, make_tuple, tie

#include "exception.hpp"

namespace tenes {

SquareLattice::SquareLattice(int X, int Y, int skew)
    : LX(X),
      LY(Y),
      N_UNIT(LX * LY),
      skew(skew),
      physical_dims(N_UNIT, -1),
      virtual_dims(N_UNIT, std::array<int, 4>{-1, -1, -1, -1}),
      initial_dirs(N_UNIT, std::vector<double>(1)),
      noises(N_UNIT, 0.0) {
  if (X <= 0) {
    throw tenes::input_error("Lattice.X should be positive");
  }
  if (Y <= 0) {
    throw tenes::input_error("Lattice.Y should be positive");
  }
  calc_neighbors();
}

void SquareLattice::calc_neighbors() {
  N_UNIT = LX * LY;
  Tensor_list.assign(LX, std::vector<int>(LY));
  initial_dirs.assign(N_UNIT, std::vector<double>(1));
  noises.assign(N_UNIT, 0.0);
  NN_Tensor.resize(N_UNIT);
  for (int ix = 0; ix < LX; ++ix) {
    for (int iy = 0; iy < LY; ++iy) {
      const int i = ix + iy * LX;
      Tensor_list[ix][iy] = i;
    }
  }
  for (int iy = 0; iy < LY; ++iy) {
    for (int ix = 0; ix < LX; ++ix) {
      const int i = ix + iy * LX;
      NN_Tensor[i][0] = Tensor_list[(ix - 1 + LX) % LX][iy];
      NN_Tensor[i][1] = Tensor_list[ix][(iy + 1) % LY];
      NN_Tensor[i][2] = Tensor_list[(ix + 1) % LX][iy];
      NN_Tensor[i][3] = Tensor_list[ix][(iy - 1 + LY) % LY];
    }
  }
  if (skew != 0) {
    for (int ix = 0; ix < LX; ++ix) {
      NN_Tensor[ix][3] = Tensor_list[(ix + LX + skew) % LX][LY - 1];
      NN_Tensor[N_UNIT + ix - LX][1] = Tensor_list[(ix + LX - skew) % LX][0];
    }
  }
  logical_check();
}

int SquareLattice::other(int index, int dx, int dy) const {
  while (dx != 0) {
    if (dx > 0) {
      index = right(index);
      --dx;
    } else {
      index = left(index);
      ++dx;
    }
  }

  while (dy != 0) {
    if (dy > 0) {
      index = top(index);
      --dy;
    } else {
      index = bottom(index);
      ++dy;
    }
  }
  return index;
}

void SquareLattice::save(const char* filename, bool append) {
  std::ofstream ofs;
  if (append) {
    ofs.open(filename, std::ios::out | std::ios::app);
    ofs << std::endl;
  } else {
    ofs.open(filename, std::ios::out);
  }

  ofs << "Lsub = [ " << LX << " , " << LY << " ]" << std::endl;
  ofs << "skew = " << skew << std::endl;
}

void SquareLattice::Bcast(MPI_Comm comm, int root) {
  int irank;
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);
  std::vector<int> params_int(3);
  std::vector<std::vector<double>> init_dirs(N_UNIT);
  std::vector<double> ns;

  int ldof = 0;

  if (irank == root) {
    params_int[0] = LX;
    params_int[1] = LY;
    params_int[2] = N_UNIT;

    init_dirs.assign(initial_dirs.begin(), initial_dirs.end());

    MPI_Bcast(&params_int.front(), 3, MPI_INT, 0, comm);

    for (int i = 0; i < N_UNIT; ++i) {
      ldof = init_dirs[i].size();
      MPI_Bcast(&ldof, 1, MPI_INT, 0, comm);
      MPI_Bcast(&init_dirs[i].front(), ldof, MPI_DOUBLE, 0, comm);
    }

    ns.assign(noises.begin(), noises.end());
    MPI_Bcast(&ns.front(), N_UNIT, MPI_DOUBLE, 0, comm);

    std::vector<int> phys(physical_dims.begin(), physical_dims.end());
    MPI_Bcast(&phys.front(), N_UNIT, MPI_INT, 0, comm);

    for (int i = 0; i < N_UNIT; ++i) {
      std::array<int, 4> vs = virtual_dims[i];
      MPI_Bcast(vs.data(), 4, MPI_INT, 0, comm);
    }

  } else {
    MPI_Bcast(&params_int.front(), 3, MPI_INT, 0, comm);
    LX = params_int[0];
    LY = params_int[1];
    N_UNIT = params_int[2];

    for (int i = 0; i < N_UNIT; ++i) {
      MPI_Bcast(&ldof, 1, MPI_INT, 0, comm);
      init_dirs[i].resize(ldof);
      MPI_Bcast(&init_dirs[i].front(), ldof, MPI_DOUBLE, 0, comm);
    }

    ns.resize(N_UNIT);
    MPI_Bcast(&ns.front(), N_UNIT, MPI_DOUBLE, 0, comm);

    std::vector<int> phys(N_UNIT);
    MPI_Bcast(&phys.front(), N_UNIT, MPI_INT, 0, comm);
    physical_dims.assign(phys.begin(), phys.end());
    for (int i = 0; i < N_UNIT; ++i) {
      std::array<int, 4> vs;
      MPI_Bcast(vs.data(), 4, MPI_INT, 0, comm);
      virtual_dims[i] = vs;
    }
  }
  calc_neighbors();
  initial_dirs.assign(init_dirs.begin(), init_dirs.end());
  noises.assign(ns.begin(), ns.end());
}

void SquareLattice::check_dims() const {
  std::vector<std::tuple<int, int, int, int>> fails;
  for (int i = 0; i < N_UNIT; ++i) {
    if (virtual_dims[i][0] != virtual_dims[left(i)][2]) {
      fails.push_back(std::make_tuple(i, 0, left(i), 2));
    }
    if (virtual_dims[i][1] != virtual_dims[top(i)][3]) {
      fails.push_back(std::make_tuple(i, 1, top(i), 3));
    }
  }
  if (!fails.empty()) {
    int i, j;
    int i_leg, j_leg;
    std::stringstream ss;
    ss << "virtual dimension mismatch :\n";
    for (auto const& tpl : fails) {
      std::tie(i, i_leg, j, j_leg) = tpl;
      ss << "  " << i_leg << " bond of " << i << " site and " << j_leg
         << " bond of " << j << " site\n";
    }
    throw tenes::input_error(ss.str());
  }
}

void SquareLattice::logical_check() const {
#ifndef NDEBUG
  for (int i = 0; i < N_UNIT; ++i) {
    // index <-> coords
    {
      int x = this->x(i);
      int y = this->y(i);
      assert(i == this->index(x, y));
    }

    // neighbors
    {
      int left = this->left(i);
      assert(i == this->right(left));
    }

    {
      int right = this->right(i);
      assert(i == this->left(right));
    }

    {
      int top = this->top(i);
      assert(i == this->bottom(top));
    }

    {
      int bottom = this->bottom(i);
      assert(i == this->top(bottom));
    }

    // loops
    { assert(top(left(i)) == left(top(i))); }
  }
#endif
}

}  // end of namespace tenes
