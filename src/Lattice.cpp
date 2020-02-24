#include "Lattice.hpp"

#include <cassert>
#include <fstream>

namespace tenes {

Lattice::Lattice(int X, int Y, int skew)
    : LX(X),
      LY(Y),
      N_UNIT(LX * LY),
      skew(skew),
      physical_dims(N_UNIT, -1),
      virtual_dims(N_UNIT, std::array<int, 4>{-1,-1,-1,-1}),
      initial_dirs(N_UNIT, std::vector<double>(1)),
      noises(N_UNIT, 0.0) {
  assert(X > 0);
  assert(Y > 0);

  calc_neighbors();
}

void Lattice::calc_neighbors() {
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
  if(skew!=0){
    for (int ix = 0; ix < LX; ++ix){
      NN_Tensor[ix][3] = Tensor_list[(ix+LX+skew)%LX][LY-1];
      NN_Tensor[N_UNIT+ix-LX][1] = Tensor_list[(ix+LX-skew)%LX][0];
    }
  }
  logical_check();
}

int Lattice::other(int index, int dx, int dy) const {
  while(dx != 0){
    if(dx>0){
      index = right(index);
      --dx;
    }else{
      index = left(index);
      ++dx;
    }
  }

  while(dy != 0){
    if(dy>0){
      index = top(index);
      --dy;
    }else{
      index = bottom(index);
      ++dy;
    }
  }
  return index;
}

void Lattice::reset(int X, int Y) {
  assert(X > 0);
  assert(Y > 0);
  LX = X;
  LY = Y;
  calc_neighbors();
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

    for(int i=0; i < N_UNIT; ++i){
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
    for(int i=0; i < N_UNIT; ++i){
      std::array<int, 4> vs;
      MPI_Bcast(vs.data(), 4, MPI_INT, 0, comm);
      virtual_dims[i] = vs;
    }
  }
  calc_neighbors();
  initial_dirs.assign(init_dirs.begin(), init_dirs.end());
  noises.assign(ns.begin(), ns.end());
}

void Lattice::check_dims() const{
  for(int i=0; i<N_UNIT; ++i){
    assert(virtual_dims[i][0] == virtual_dims[left(i)][2]);
    assert(virtual_dims[i][1] == virtual_dims[top(i)][3]);
    assert(virtual_dims[i][2] == virtual_dims[right(i)][0]);
    assert(virtual_dims[i][3] == virtual_dims[bottom(i)][1]);
  }
}

void Lattice::logical_check() const{
#ifndef NDEBUG
  for(int i=0; i<N_UNIT; ++i){
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
    {
      assert(top(left(i)) == left(top(i)));
    }
  }
#endif
}


}  // end of namespace tenes
