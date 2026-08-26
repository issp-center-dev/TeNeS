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

#ifndef TENES_SRC_SQUARELATTICE_HPP_
#define TENES_SRC_SQUARELATTICE_HPP_

#include <array>
#include <vector>

#include "mpi.hpp"

namespace tenes {

/*! @brief Square Lattice (unit cell)
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
class SquareLattice {
 public:
  int LX;             //!< unit-cell width
  int LX_noskew;      //!< unit-cell width before skew expansion
  int LY;             //!< unit-cell height
  int LY_noskew;      //!< unit-cell height before skew expansion
  int N_UNIT;         //!< number of sites in the unit cell
  int N_UNIT_noskew;  //!< number of sites before skew expansion

  /*!
   * @brief Skew boundary condition
   *
   * <tt>T(x, y) = T(x + skew, y + LY)</tt>
   */
  int skew;

  //! Physical (local Hilbert space) dimension of each site.
  std::vector<int> physical_dims;
  //! Virtual bond dimensions of each site, indexed by [site][direction].
  std::vector<std::array<int, 4>> virtual_dims;

  //! Initial state amplitude of each site in the physical basis
  //! (tensor.unitcell.initial_state); all zero means random.
  std::vector<std::vector<double>> initial_dirs;
  //! Amplitude of the random noise added to each site's initial tensor
  //! (tensor.unitcell.noise).
  std::vector<double> noises;

  //! Build an X x Y unit cell with the given skew boundary condition.
  SquareLattice(int X, int Y, int skew = 0);

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

  /*! parity of site
   *
   *  @param[in] index site index
   *
   *  @return parity (1 or -1)
   *
   *  @note
   *  parity of the site 0 is +1
   *
   */
  int parity(int index) const { return parities[index]; }

  //! Append the lattice description to a file (rank 0 only).
  void save(const char *filename, bool append = false);
  //! save() with append = true.
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
  std::vector<std::vector<int>> Tensor_list;  //!< index of site [x][y]
  std::vector<std::array<int, 4>>
      NN_Tensor;              //!< index of nearest neighbor site [site][bond]
  std::vector<int> parities;  //!< sublattice parity (+1/-1) of each site
  //! Validate the lattice sizes; throws tenes::input_error on nonsense.
  void logical_check() const;
  //! Fill Tensor_list, NN_Tensor, and parities from the sizes.
  void calc_neighbors();
};

}  // end of namespace tenes

#endif  // TENES_SRC_SQUARELATTICE_HPP_
