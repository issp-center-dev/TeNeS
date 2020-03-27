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

#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <array>
#include <vector>

#include "mpi.hpp"

namespace tenes {

// Lattice setting

/*
 * axis:
 *
 *  y
 *  ^
 *  |
 *  .->x
 *
 *  order:
 *
 *  2 3
 *  0 1
 *
 * edge index:
 *
 *   1
 *  0.2
 *   3
 *
 */
class Lattice {
 public:
  int LX;
  int LY;
  int N_UNIT;

  int skew;

  std::vector<std::vector<int>> Tensor_list;
  std::vector<std::array<int, 4>> NN_Tensor;

  std::vector<int> physical_dims;
  std::vector<std::array<int, 4>> virtual_dims;

  std::vector<std::vector<double>> initial_dirs;
  std::vector<double> noises;

  Lattice(int X, int Y, int skew=0);

  int x(int index) const { return index % LX; }
  int y(int index) const { return index / LX; }
  int index(int x, int y) const {
    int X = (x < 0) ? (x % LX + LX) : x % LX;  // c++11 requires neg%pos is neg
    int Y = (y < 0) ? (y % LY + LY) : y % LY;
    return X + Y * LX;
  }

  int neighbor(int index, int direction) const { return NN_Tensor[index][direction]; }
  int left(int index) const { return neighbor(index, 0); }
  int right(int index) const { return neighbor(index, 2); }
  int top(int index) const { return neighbor(index, 1); }
  int bottom(int index) const { return neighbor(index, 3); }

  int other(int index, int dx, int dy) const;

  void reset(int X, int Y);
  void calc_neighbors();

  void save(const char *filename, bool append = false);

  void save_append(const char *filename) { save(filename, true); }

  void Bcast(MPI_Comm comm, int root = 0);

  void check_dims() const;

private:
  void logical_check() const;
};

/*

class Lattice_skew : public Lattice {
public:
  int LX_ori;
  int LY_ori;
  int LX_diff;
  int LY_diff;

  std::vector<std::vector<int> > Tensor_position;

  Lattice_skew() {
    LX_ori = 2;
    LY_ori = 2;
    LX = LX_ori;
    LY = LY_ori;
    LX_diff = LX_ori;
    LY_diff = LY_ori;
    N_UNIT = LX_ori * LY_ori;

    Tensor_position =
        std::vector<std::vector<int> >(N_UNIT, std::vector<int>(2));
    Tensor_list = std::vector<std::vector<int> >(LX, std::vector<int>(LY));
    NN_Tensor = std::vector<std::vector<int> >(N_UNIT, std::vector<int>(4));
  }

  void read_parameters(const char *filename) {
    std::ifstream input_file;
    input_file.open(filename, std::ios::in);
    std::string reading_line_buffer;

    while (!input_file.eof()) {
      std::getline(input_file, reading_line_buffer);
      std::stringstream buf(reading_line_buffer);
      std::vector<std::string> result;
      while (buf >> reading_line_buffer) {
        result.push_back(reading_line_buffer);
      }

      if (result.size() > 1) {
        if (result[0].compare("LX_ori") == 0) {
          std::istringstream is(result[1]);
          is >> LX_ori;
        } else if (result[0].compare("LY_ori") == 0) {
          std::istringstream is(result[1]);
          is >> LY_ori;
        } else if (result[0].compare("LX") == 0) {
          std::istringstream is(result[1]);
          is >> LX;
        } else if (result[0].compare("LY") == 0) {
          std::istringstream is(result[1]);
          is >> LY;
        } else if (result[0].compare("LX_diff") == 0) {
          std::istringstream is(result[1]);
          is >> LX_diff;
        } else if (result[0].compare("LY_diff") == 0) {
          std::istringstream is(result[1]);
          is >> LY_diff;
        } else if (result[0].compare("N_UNIT") == 0) {
          std::istringstream is(result[1]);
          is >> N_UNIT;
        }
        // std::cout<< "## input data: "<<result[0]<<" =
        // "<<result[1]<<std::endl;
      }
    }
  }
  void output_parameters(const char *filename, bool append) {
    std::ofstream ofs;
    if (append) {
      ofs.open(filename, std::ios::out | std::ios::app);
    } else {
      ofs.open(filename, std::ios::out);
    }

    ofs << "LX_ori " << LX_ori << std::endl;
    ofs << "LY_ori " << LY_ori << std::endl;
    ofs << "LX " << LX << std::endl;
    ofs << "LY " << LY << std::endl;
    ofs << "LX_diff " << LX_diff << std::endl;
    ofs << "LY_diff " << LY_diff << std::endl;
    ofs << "N_UNIT " << N_UNIT << std::endl;
  }

  void output_parameters(const char *filename) {
    output_parameters(filename, false);
  }

  void output_parameters_append(const char *filename) {
    output_parameters(filename, true);
  }

  void Bcast_parameters(MPI_Comm comm) {
    int irank;
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
    std::vector<int> params_int(7);

    if (irank == 0) {
      params_int[0] = LX_ori;
      params_int[1] = LY_ori;
      params_int[2] = LX;
      params_int[3] = LY;
      params_int[4] = LX_diff;
      params_int[5] = LY_diff;
      params_int[6] = N_UNIT;

      MPI_Bcast(&params_int.front(), 7, MPI_INT, 0, comm);
    } else {
      MPI_Bcast(&params_int.front(), 7, MPI_INT, 0, comm);

      LX_ori = params_int[0];
      LY_ori = params_int[1];
      LX = params_int[2];
      LY = params_int[3];
      LX_diff = params_int[4];
      LY_diff = params_int[5];
      N_UNIT = params_int[6];
    }
  }
};

*/

}  // end of namespace tenes

#endif  // LATTICE_HPP
