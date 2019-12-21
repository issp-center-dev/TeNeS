#include <cassert>
#include <fstream>

#include "Lattice.hpp"

Lattice::Lattice(int X, int Y) : LX(X), LY(Y), N_UNIT(LX * LY), initial_dirs(N_UNIT, -1) {
  assert(X > 0);
  assert(Y > 0);

  reset();
}
void Lattice::reset() {
  N_UNIT = LX * LY;
  Tensor_list.assign(LX, std::vector<int>(LY));
  initial_dirs.assign(N_UNIT, -1);
  NN_Tensor.resize(N_UNIT);
  for (int ix = 0; ix < LX; ++ix) {
    for (int iy = 0; iy < LY; ++iy) {
      const int i = ix + iy * LX;
      Tensor_list[ix][iy] = i;
    }
  }
  for (int ix = 0; ix < LX; ++ix) {
    for (int iy = 0; iy < LY; ++iy) {
      const int i = ix + iy * LX;
      NN_Tensor[i][0] = Tensor_list[(ix - 1 + LX) % LX][iy];
      NN_Tensor[i][1] = Tensor_list[ix][(iy + 1) % LY];
      NN_Tensor[i][2] = Tensor_list[(ix + 1) % LX][iy];
      NN_Tensor[i][3] = Tensor_list[ix][(iy - 1 + LY) % LY];
    }
  }
}
void Lattice::reset(int X, int Y) {
  assert(X > 0);
  assert(Y > 0);
  LX = X;
  LY = Y;
  reset();
}

void Lattice::save(const char *filename, bool append) {
  std::ofstream ofs;
  if (append) {
    ofs.open(filename, std::ios::out | std::ios::app);
  } else {
    ofs.open(filename, std::ios::out);
  }

  ofs << "Lsub = [ " << LX << " , " << LY << " ]" << std::endl;
}

void Lattice::Bcast(MPI_Comm comm, int root) {
  int irank;
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);
  std::vector<int> params_int(3);
  std::vector<int> init_dirs;

  if (irank == root) {
    params_int[0] = LX;
    params_int[1] = LY;
    params_int[2] = N_UNIT;

    init_dirs.assign(initial_dirs.begin(), initial_dirs.end());

    MPI_Bcast(&params_int.front(), 3, MPI_INT, 0, comm);
    MPI_Bcast(&init_dirs.front(), N_UNIT, MPI_INT, 0, comm);
  } else {
    MPI_Bcast(&params_int.front(), 3, MPI_INT, 0, comm);
    initial_dirs.resize(N_UNIT);
    MPI_Bcast(&init_dirs.front(), N_UNIT, MPI_INT, 0, comm);

    LX = params_int[0];
    LY = params_int[1];
    N_UNIT = params_int[2];
  }
  reset();
  initial_dirs.assign(init_dirs.begin(), init_dirs.end());
}
