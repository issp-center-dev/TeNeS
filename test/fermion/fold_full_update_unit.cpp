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
//! would mean declaring every fixture in a header of its own, and the two that
//! would be worth separating (fold_geometry at 5.4 s and ctm_phase at 4.2 s)
//! would still leave FreeFermionFull as the longest test in a parallel ctest,
//! so the work would buy nothing.

// Its own executable, so ctest can run it next to the rest of the fermion
// tests rather than after them: these cases are 12 s of the 18 s the single
// test_fermion_layer binary used to take, and a ctest test is one process.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../test_fermion_common.hpp"

#include "fold_geometry.cpp"
#include "full_update_env.cpp"
#include "full_update_bond.cpp"
#include "decomposition_diagnostics.cpp"
#include "full_update_realctm.cpp"
#include "ctm_phase.cpp"
