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
//! The r = 2 convention tests and the mean-field measurement tests, compiled
//! as one translation unit: mf_measure.cpp measures through the r2_* and r4_*
//! fixtures that r2_convention.cpp sets up.

#include "../test_fermion_common.hpp"

#include "r2_convention.cpp"
#include "mf_measure.cpp"
