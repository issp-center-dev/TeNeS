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
//! The fermion guard tests and the fast-full-update tests, compiled as one
//! translation unit: fast_full_update.cpp captures the solver's streams with
//! the fgd_stream_capture of fermion_guards.cpp.

#include "../test_fermion_common.hpp"

#include "fermion_guards.cpp"
#include "fast_full_update.cpp"
