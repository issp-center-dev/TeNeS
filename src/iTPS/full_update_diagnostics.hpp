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

/*! @file
 *  @brief User-facing text for a full update that could not continue.
 *
 *  A failed decomposition ends the run, and the run has written nothing
 *  yet: save_tensors() only follows a completed optimize(). The message is
 *  therefore all the user gets, so it names what failed, what the graded
 *  decomposition saw, which of the solver's own warnings to look back at,
 *  and which input.toml keys to change.
 */

#ifndef TENES_SRC_ITPS_FULL_UPDATE_DIAGNOSTICS_HPP_
#define TENES_SRC_ITPS_FULL_UPDATE_DIAGNOSTICS_HPP_

#include <sstream>
#include <string>

namespace tenes {
namespace itps {

//! What a full-update failure message should send the reader after first.
enum class full_update_failure_lead {
  state,         //!< The run's own numerics: the default.
  library,       //!< An info no conforming LAPACK can return.
  process_grid,  //!< ScaLAPACK reported that the MPI ranks disagreed.
};

/*!
 * @brief Assemble the error text for a failed fermionic full update.
 *
 * @param[in] what Which step gave up, e.g. "balancing SVD".
 * @param[in] detail One line of evidence, typically
 *        tenes::fermion::decomposition_diagnostics::describe().
 * @param[in] lead What to send the reader after first.
 * @return The message to hand tenes::runtime_error.
 */
inline std::string fermion_full_update_failure_message(
    const std::string &what, const std::string &detail,
    full_update_failure_lead lead = full_update_failure_lead::state) {
  std::ostringstream ss;
  ss << "fermion full update: " << what << " failed.\n"
     << "  " << detail << "\n";
  if (lead == full_update_failure_lead::process_grid) {
    ss << "That is ScaLAPACK's MIN(M,N)+1: pdgesvd found that the MPI ranks\n"
       << "did not agree on this block's singular values. It is a property "
          "of\n"
       << "spreading the decomposition over a process grid, not of the "
          "state;\n"
       << "nearly degenerate singular values are where ranks stop agreeing,\n"
       << "and the graded decomposition's parity blocks often have them.\n"
       << "Blocks of at most TENES_FERMION_LOCAL_DECOMP_MAX elements (4096 "
          "by\n"
       << "default) are already factorized on every rank instead of over the\n"
       << "grid, so this block was larger than that. What to do:\n"
       << "  - raise TENES_FERMION_LOCAL_DECOMP_MAX past this block's "
          "element\n"
       << "    count, so it is factorized on one rank too\n"
       << "  - rerun with fewer MPI processes; one rank cannot disagree "
          "with\n"
       << "    itself\n"
       << "  - check that every rank is the same CPU model with the same "
          "math\n"
       << "    library dispatch (heterogeneous nodes produce exactly this)\n"
       << "The settings below are unlikely to help.\n";
  } else if (lead == full_update_failure_lead::library) {
    ss << "The info above is outside what the LAPACK documentation allows "
          "for a\n"
       << "block that size, so no state can have produced it: check the "
          "LAPACK\n"
       << "and BLAS this binary is linked against (and try another "
          "implementation\n"
       << "or version) before anything else. The settings below are\n"
       << "unlikely to help; they are listed only in case the library\n"
       << "turns out to be innocent.\n";
  } else {
    ss << "This is a condition of the run, not an internal error: the "
          "two-site\n"
       << "state became numerically intractable. Look just above this line "
          "for\n"
       << "\"CTM did not converge\" or \"kept an empty parity sector\" - "
          "those\n"
       << "are usually the real cause.\n";
  }
  ss << "What to change, in order:\n"
     << "  1. parameter.ctm.iteration_max, parameter.ctm.convergence_epsilon\n"
     << "     - let the corner transfer matrices converge first\n"
     << "  2. parameter.ctm.dimension - chi >= D*D is recommended\n"
     << "  3. parameter.simple_update.num_step - more simple update first, so\n"
     << "     the full update starts from a better state\n"
     << "  4. parameter.full_update.tau - a smaller imaginary-time step\n"
     << "This run wrote no tensors: they are saved only after the whole\n"
     << "optimization finishes. Set parameter.general.tensor_save in a run "
        "that\n"
     << "does finish, and a later run can resume from it through\n"
     << "parameter.general.tensor_load instead of repeating the simple "
        "update.";
  return ss.str();
}

/*!
 * @brief Name the bond and the sweep step a full-update error came from.
 *
 * @param[in] source Site the gate is anchored at.
 * @param[in] target Its neighbour across the bond.
 * @param[in] source_leg Leg of @p source pointing at @p target, in the
 *        SquareLattice order (0 left, 1 top, 2 right, 3 bottom).
 * @param[in] step_index Zero-based index of the full-update sweep.
 * @param[in] nsteps Total number of sweeps configured.
 * @return A clause to append to the message, reported one-based.
 */
inline std::string full_update_bond_context(int source, int target,
                                            int source_leg, int step_index,
                                            int nsteps) {
  static const char *const leg_name[4] = {"left", "top", "right", "bottom"};
  const char *const dir =
      (source_leg >= 0 && source_leg < 4) ? leg_name[source_leg] : "unknown";
  std::ostringstream ss;
  ss << "\nIt happened at the bond from site " << source << " to its " << dir
     << " neighbour site " << target << ", full-update step "
     << (step_index + 1) << "/" << nsteps << ".";
  return ss.str();
}

}  // namespace itps
}  // namespace tenes

#endif  // TENES_SRC_ITPS_FULL_UPDATE_DIAGNOSTICS_HPP_
