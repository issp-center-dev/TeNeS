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

//! @file
//! The folding-geometry and full-update test cases, compiled as one
//! translation unit.  They are written as one: full_update_env.cpp uses the
//! fixtures of fold_geometry.cpp (fg_*), full_update_bond.cpp uses both
//! (fub_*), decomposition_diagnostics.cpp and full_update_realctm.cpp build on
//! those in turn, and ctm_phase.cpp reaches into all of them.  Separating them
//! would mean declaring every fixture in a header of its own, which buys
//! nothing: this unit is already the longest pole of the executable, and the
//! point of the split is to fill a CI runner's other cores, not to make each
//! unit minimal.

#include "../test_fermion_common.hpp"

#include "fold_geometry.cpp"
#include "full_update_env.cpp"
#include "full_update_bond.cpp"
#include "decomposition_diagnostics.cpp"
#include "full_update_realctm.cpp"
#include "ctm_phase.cpp"
