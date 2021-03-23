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

#ifndef TENES_SRC_LATTICE_HPP_
#define TENES_SRC_LATTICE_HPP_

#include <array>
#include <vector>

#include "mpi.hpp"

namespace tenes {

//! Lattice (unit cell)
/*!
 * Geometries:
 *
 * axis:
 *
 * @verbatim
   y
   ^
   |
   .->x
   @endverbatim
 *
 * order of sites:
 *
 * @verbatim
    2 3
    0 1
   @endverbatim
 *
 * order of bonds or direction:
 *
 * @verbatim
    1
   0.2
    3
  @endverbatim
 *
 */
class Lattice {
 public:
  int LX;
  int LY;
  int N_UNIT;

  /*!
   * @brief Skew boundary condition
   *
   * <tt>T(x, y) = T(x + skew, y + LY)</tt>
   */
  int skew;

  std::vector<int> physical_dims;
  std::vector<std::array<int, 4>> virtual_dims;

  std::vector<std::vector<double>> initial_dirs;
  std::vector<double> noises;

  Lattice(int X, int Y, int skew = 0);

  //! x coordinate
  int x(int index) const { return index % LX; }
  //! y coordinate
  int y(int index) const { return index / LX; }
  //! indexing w/o boundary calculation
  int index_fast(int x, int y) const { return Tensor_list[x][y]; }
  //! indexing w/ boundary calculation
  int index(int x, int y) const {
    int y_offset = 0;
    if (y >= LY) {
      y_offset = y / LY;
    } else if (y < 0) {
      y_offset = (y + 1) / LY - 1;
    }
    y -= LY * y_offset;
    x -= skew * y_offset;
    x %= LX;
    if (x < 0) {
      x += LX;
    }

    return Tensor_list[x][y];
  }

  //! neighbor site
  int neighbor(int index, int direction) const {
    return NN_Tensor[index][direction];
  }

  //! left neighbor site
  int left(int index) const { return neighbor(index, 0); }
  //! right neighbor site
  int right(int index) const { return neighbor(index, 2); }
  //! top neighbor site
  int top(int index) const { return neighbor(index, 1); }
  //! bottom neighbor site
  int bottom(int index) const { return neighbor(index, 3); }

  /*! other site index
   *
   *  @param[in] index site index
   *  @param[in] dx displacement along x
   *  @param[in] dy displacement along y
   */
  int other(int index, int dx, int dy) const;

  void save(const char *filename, bool append = false);
  void save_append(const char *filename) { save(filename, true); }

  /*! @brief broadcast lattice
   *
   *  @param[in] comm communicator
   *  @param[in] root root rank
   */
  void Bcast(MPI_Comm comm, int root = 0);

  /*! @brief check consistency of bond dimensions
   *
   *  @throw tenes::input_error
   *
   */
  void check_dims() const;

 private:
  std::vector<std::vector<int>> Tensor_list;
  std::vector<std::array<int, 4>> NN_Tensor;
  void logical_check() const;
  void calc_neighbors();
};

}  // end of namespace tenes

#endif  // TENES_SRC_LATTICE_HPP_
