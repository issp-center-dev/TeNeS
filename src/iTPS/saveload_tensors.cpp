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

#define _USE_MATH_DEFINES
#include <algorithm>
#include <array>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <cstdlib>

#include "iTPS.hpp"

#include <mptensor/file_io/file_io.hpp>

#include "../tensor.hpp"

#include "../fermion/fermion_info.hpp"
#include "../fermion/fops.hpp"
#include "../printlevel.hpp"
#include "../util/datetime.hpp"
#include "../util/file.hpp"
#include "../util/string.hpp"

using std::size_t;

namespace tenes::itps {

namespace {

/*! @brief Integer-per-line reader over a metadata file held whole in a string.
 *
 * The point is that every rank parses the same bytes. Parsing on rank 0 and
 * broadcasting the outcome piece by piece is what let a short or malformed
 * file abort rank 0 inside a collective, leaving the others waiting in the
 * next broadcast, and it is what let a shape loop run past the end of a
 * checkpoint that holds fewer sites than the input asks for.
 */
class metadata_reader {
 public:
  metadata_reader(std::string const &content, std::string const &filename)
      : iss_(content), filename_(filename) {}

  //! Next line with its comment dropped; throws at end of file.
  std::string next_line() {
    if (!std::getline(iss_, line_)) {
      throw tenes::load_error("ERROR: " + filename_ + " ends unexpectedly");
    }
    return util::drop_comment(line_);
  }

  //! The integers on the next line; throws if any word is not an integer.
  std::vector<int> next_ints() {
    last_ = next_line();
    const auto words = util::split(last_);
    std::vector<int> ret;
    ret.reserve(words.size());
    for (const auto &w : words) {
      try {
        ret.push_back(std::stoi(w));
      } catch (const std::exception &) {
        throw tenes::load_error("ERROR: cannot parse " + filename_ +
                                " line as integers: \"" + last_ + "\"");
      }
    }
    return ret;
  }

  //! The first integer on the next line; throws if the line has none.
  int next_scalar() {
    const auto values = next_ints();
    if (values.empty()) {
      throw tenes::load_error("ERROR: expected an integer in " + filename_ +
                              " line: \"" + last_ + "\"");
    }
    return values[0];
  }

 private:
  std::istringstream iss_;
  std::string line_;
  std::string last_;
  std::string filename_;
};

/*! @brief Whether @p name is one of the files a checkpoint is made of.
 *
 * Used to decide what a save may delete from its destination: files of this
 * shape that the save did not write are leftovers of an earlier run, anything
 * else belongs to the user and is left alone.
 */
bool is_checkpoint_name(const std::string &name) {
  if (name == "params.dat" || name == "fermion.dat") {
    return true;
  }
  auto all_digits = [](const std::string &s) {
    return !s.empty() && s.find_first_not_of("0123456789") == std::string::npos;
  };
  std::string stem = name;
  // mptensor writes one <base>.<rank>.bin and <base>.<rank>.idx per process.
  for (const std::string &suffix : {std::string(".bin"), std::string(".idx")}) {
    if (stem.size() > suffix.size() &&
        stem.compare(stem.size() - suffix.size(), suffix.size(), suffix) == 0) {
      const std::size_t dot = stem.rfind('.', stem.size() - suffix.size() - 1);
      if (dot == std::string::npos) {
        return false;
      }
      if (!all_digits(
              stem.substr(dot + 1, stem.size() - suffix.size() - dot - 1))) {
        return false;
      }
      stem = stem.substr(0, dot);
      break;
    }
  }
  if (stem.size() < 5 || stem.compare(stem.size() - 4, 4, ".dat") != 0) {
    return false;
  }
  stem = stem.substr(0, stem.size() - 4);
  const std::size_t underscore = stem.rfind('_');
  if (underscore == std::string::npos ||
      !all_digits(stem.substr(underscore + 1))) {
    return false;
  }
  const std::string prefix = stem.substr(0, underscore);
  for (const char *known :
       {"T", "El", "Et", "Er", "Eb", "C1", "C2", "C3", "C4", "lambda"}) {
    if (prefix == known) {
      return true;
    }
  }
  return false;
}

//! Working directory a save builds its checkpoint in, inside the destination.
//! TeNeS-specific because a save removes it unconditionally: the earlier name,
//! ".tmp", took a directory the user kept in the destination.
constexpr const char *save_work_dir_name = ".tenes-save-tmp";
//! Marker a save leaves while it moves files into place; see save_tensors().
//! TeNeS-specific for the same reason as save_work_dir_name.
constexpr const char *save_marker_name = ".tenes-save-incomplete";
//! How much of a marker's body a refused load reads, broadcasts and shows.
constexpr std::size_t marker_body_limit = 64 * 1024;

/*! @brief Move every entry of @p work_dir into @p dest, then confirm it is
 *         empty.
 *
 * Used both to finish a move an interrupted save left behind and to put a new
 * checkpoint in place. Stops at the first rename that fails and does not put
 * back what already moved: a rollback that fails halfway leaves a worse mess
 * than the marker describes.
 *
 * Success also requires @p work_dir to be empty afterwards. A listing that
 * came back short -- which entries_with_prefix() reports as a shorter list,
 * not as an error -- would otherwise pass for "everything moved", and the
 * caller would go on to sweep the destination and remove the files it missed.
 *
 * @param[out] moved_names names moved into @p dest
 * @param[out] failed what could not be moved, for a warning
 * @return true if everything moved and @p work_dir is empty
 */
bool move_all_into(const std::string &work_dir, const std::string &dest,
                   std::set<std::string> &moved_names, std::string &failed) {
  for (const auto &path : util::entries_with_prefix(work_dir, "")) {
    const std::string name = util::basename(path);
    if (!util::rename(path, dest + "/" + name)) {
      failed = name;
      return false;
    }
    moved_names.insert(name);
  }
  if (!util::is_empty_directory(work_dir)) {
    failed = "the files a listing of " + work_dir + " did not return";
    return false;
  }
  return true;
}

//! The whole of the file at @p path into @p content; false if it cannot be
//! read.
bool read_whole_file(const std::string &path, std::string &content) {
  std::ifstream ifs(path, std::ios::binary);
  if (!ifs) {
    return false;
  }
  content.assign(std::istreambuf_iterator<char>(ifs),
                 std::istreambuf_iterator<char>());
  return !ifs.bad();
}

/*! @brief Whether the files mptensor's save() wrote for @p t on this rank are
 *         there to their last byte.
 *
 * mptensor writes them without looking at its streams, so a write that failed
 * -- a full disk, a quota, a file-size limit -- used to leave a checkpoint
 * that said it was saved and loaded another state. A write that keeps failing
 * until the file is closed loses only the end of the file. One whose error
 * clears before the close -- space freed on a shared file system, say -- can
 * instead leave part of its buffer written twice, because an ofstream writes
 * what it still holds once more when it is closed. The checks below tell both
 * apart from a whole file:
 *
 * - @c <base>.<rank>.bin holds exactly local_size() elements.
 * - @c <base>.<rank>.idx gives local_size, the local rows r and the local
 *   columns c on its first three lines, then a header line and the global
 *   index of each row, and the same for the columns. save_index() breaks the
 *   line after every B indices and after each list, so the file holds exactly
 *   7 + r / B + c / B newlines, and a cut of any length leaves fewer. That
 *   includes a cut of the final empty line alone, which on a rank holding no
 *   elements is everything that follows the last header. It also holds
 *   exactly 8 + r + c words -- the six header words and values, the two list
 *   headers and the indices -- which a buffer written twice changes even when
 *   the newlines happen to match: a sweep over every failure point of such
 *   files found no file that kept both counts.
 * - @c <base>, which rank 0 alone writes, is exactly 8 lines.
 *
 * The counts follow mptensor's file format. Should a new version change it,
 * every save reports a failure, which ctest shows at once, instead of a cut
 * file passing for a whole one.
 */
template <class tensor>
bool saved_completely(const tensor &t, const std::string &base) {
  const int rank = t.get_comm_rank();
  const auto count_newlines = [](const std::string &s) {
    return static_cast<std::size_t>(std::count(s.begin(), s.end(), '\n'));
  };

  std::error_code ec;
  const auto data_size = std::filesystem::file_size(
      mptensor::io_helper::binary_filename(base, rank), ec);
  if (ec || data_size != sizeof(typename tensor::value_type) * t.local_size()) {
    return false;
  }

  std::string index;
  if (!read_whole_file(mptensor::io_helper::index_filename(base, rank),
                       index)) {
    return false;
  }
  std::istringstream header(index);
  std::string key_size, key_rows, key_cols;
  std::size_t local_size = 0, rows = 0, cols = 0;
  if (!(header >> key_size >> local_size >> key_rows >> rows >> key_cols >>
        cols) ||
      key_size != "local_size=" || key_rows != "local_n_row=" ||
      key_cols != "local_n_col=" || local_size != t.local_size()) {
    return false;
  }
  using matrix = std::decay_t<decltype(t.get_matrix())>;
  const std::size_t per_line =
      matrix::matrix_type_tag == MATRIX_TYPE_TAG_SCALAPACK ? 16 : 10;
  std::istringstream all_words(index);
  std::size_t words = 0;
  for (std::string word; all_words >> word;) {
    ++words;
  }
  if (count_newlines(index) != 7 + rows / per_line + cols / per_line ||
      words != 8 + rows + cols) {
    return false;
  }

  if (rank == 0) {
    std::string base_file;
    if (!read_whole_file(base, base_file) || count_newlines(base_file) != 8) {
      return false;
    }
  }
  return true;
}

}  // namespace

template <class ptensor>
bool iTPS<ptensor>::save_tensors() const {
  std::string const &save_dir = peps_parameters.tensor_save_dir;
  if (save_dir.empty()) {
    return true;
  }

  /* A checkpoint is built inside the destination, in a working directory
   * named by save_work_dir_name, and moved into place one file at a time when
   * it is complete.
   *
   * Crash safety: a run killed while writing -- a job hitting its wall clock,
   * say -- leaves the destination holding the previous checkpoint, because
   * nothing has been moved yet. The window in which the destination can hold a
   * mixture is the move, which is metadata only.
   *
   * The destination directory itself is never renamed or removed. An earlier
   * design swapped it wholesale, which is how "tensor_save = ." came to empty
   * the working directory and how a symlinked destination came to be replaced
   * by a real one. Here a destination that is a symbolic link is followed like
   * any other path, whatever else the user keeps in the destination survives,
   * and only the destination needs to be writable -- not its parent.
   */
  const std::string work_dir = save_dir + "/" + save_work_dir_name;
  const std::string marker = save_dir + "/" + save_marker_name;

  // 0: go ahead; 1: an interrupted move could not be finished;
  // 2: the working directory could not be created.
  int prepared = 0;
  std::string repair_failed;
  if (mpirank == 0) {
    bool repaired = true;
    if (util::path_exists(marker)) {
      /* A save before this one died partway through moving its files. Do
       * what the marker tells a person to do -- move the rest in, then delete
       * the marker -- before anything else. Clearing the working directory as
       * a mere leftover would destroy the files that repair needs and leave
       * the marker pointing at nothing, and a failure later in this save
       * would make that permanent.
       */
      std::set<std::string> moved_in;
      repaired = !util::path_exists(work_dir) ||
                 move_all_into(work_dir, save_dir, moved_in, repair_failed);
      if (repaired && !util::remove_all(marker)) {
        // Going on would write a checkpoint that the marker refuses, and
        // remove the working directory the marker points at.
        repaired = false;
        repair_failed = "the deletion of " + marker;
      }
    }
    if (!repaired) {
      prepared = 1;
    } else {
      // With no marker, a working directory was left by a save that died
      // before it moved anything: the destination is whole, so it can go.
      util::remove_all(work_dir);
      prepared = util::mkdir(work_dir) ? 0 : 2;
    }
  }
  // Doubles as the barrier that keeps the other ranks from writing into a
  // directory rank 0 has not made yet.
  bcast(prepared, 0, comm);
  if (prepared != 0) {
    // Not an exception: save_tensors() runs before measure(), so throwing here
    // would throw away a finished optimization over a directory permission.
    if (mpirank == 0) {
      if (prepared == 1) {
        std::cerr << "WARNING: " << marker
                  << " shows that an earlier save was interrupted while moving "
                     "its files, and finishing that move failed at "
                  << repair_failed << ", so no checkpoint was saved.\n"
                  << "  To finish the move by hand, move everything left in "
                  << work_dir << " into " << save_dir << ", then delete "
                  << marker << "." << std::endl;
      } else {
        std::cerr << "WARNING: cannot create the working directory " << work_dir
                  << " (" << save_dir
                  << " has to be writable), so no checkpoint was saved."
                  << std::endl;
      }
    }
    return false;
  }

  double wrote = 1.0;
  // Rank 0 alone. Every rank used to open params.dat and write the same
  // bytes into it at the same time.
  if (mpirank == 0) {
    // metadata
    std::string filename = work_dir + "/params.dat";
    std::ofstream ofs(filename.c_str());

    // Version 2 records whether the run was fermionic. Before it, the loader
    // inferred that from the presence of fermion.dat, which made a leftover
    // ledger change how a checkpoint was read.
    constexpr int tensor_format_version = 2;
    ofs << tensor_format_version << " # Format_Version\n";
    ofs << (finfo.enabled ? 1 : 0) << " # Fermion\n";
    ofs << N_UNIT << " # N_UNIT\n";
    ofs << CHI << " # CHI\n";
    for (int i = 0; i < N_UNIT; ++i) {
      for (int j = 0; j < nleg; ++j) {
        ofs << lattice.virtual_dims[i][j] << " ";
      }
      ofs << lattice.physical_dims[i] << " # Shape of Tn[" << i << "]\n";
    }
    ofs.flush();
    if (!ofs) {
      wrote = 0.0;
    }
  }
  // The CTM environment a fermionic run holds cannot be written out as it
  // stands: update_CTM() builds it through Calc_CTM_Environment_density, so
  // the edge tensors are folded to the single-layer form (CHI, CHI, D*D),
  // while initialize_tensors() allocates the double-layer (CHI, CHI, D, D)
  // that load_tensor() then insists on. A run that took a full-update step
  // therefore used to write a checkpoint it could not read back.
  //
  // Saving the environment in the shape a fresh run starts with costs nothing,
  // because the loaded environment is never used: Calc_CTM_Environment and its
  // density twin default to initialize = true and update_CTM() calls them that
  // way in both branches, so the environment is rebuilt from the site tensors
  // whichever mode we are in. Zero rather than whatever the allocator left
  // behind, so that the checkpoint of a given state is the same file every
  // time.
  const auto save_tensor = [&wrote](const ptensor &t, const std::string &path) {
    t.save(path);
    if (!saved_completely(t, path)) {
      wrote = 0.0;
    }
  };
  const auto save_placeholder = [this, &save_tensor](
                                    const mptensor::Shape &shape,
                                    const std::string &path) {
    ptensor t(comm, shape);
    for (size_t n = 0; n < t.local_size(); ++n) {
      t.set_value(t.global_index(n), typename ptensor::value_type(0.0));
    }
    save_tensor(t, path);
  };
  for (int i = 0; i < N_UNIT; ++i) {
    std::string filename = work_dir + "/";
    std::string suffix = "_" + std::to_string(i) + ".dat";
    save_tensor(Tn[i], filename + "T" + suffix);
    if (finfo.enabled) {
      const auto vdim = lattice.virtual_dims[i];
      save_placeholder(mptensor::Shape(CHI, CHI, vdim[1], vdim[1]),
                       filename + "Et" + suffix);
      save_placeholder(mptensor::Shape(CHI, CHI, vdim[2], vdim[2]),
                       filename + "Er" + suffix);
      save_placeholder(mptensor::Shape(CHI, CHI, vdim[3], vdim[3]),
                       filename + "Eb" + suffix);
      save_placeholder(mptensor::Shape(CHI, CHI, vdim[0], vdim[0]),
                       filename + "El" + suffix);
      for (const char *name : {"C1", "C2", "C3", "C4"}) {
        save_placeholder(mptensor::Shape(CHI, CHI), filename + name + suffix);
      }
    } else {
      save_tensor(eTt[i], filename + "Et" + suffix);
      save_tensor(eTr[i], filename + "Er" + suffix);
      save_tensor(eTb[i], filename + "Eb" + suffix);
      save_tensor(eTl[i], filename + "El" + suffix);
      save_tensor(C1[i], filename + "C1" + suffix);
      save_tensor(C2[i], filename + "C2" + suffix);
      save_tensor(C3[i], filename + "C3" + suffix);
      save_tensor(C4[i], filename + "C4" + suffix);
    }
  }
  if (mpirank == 0) {
    for (int i = 0; i < N_UNIT; ++i) {
      std::ofstream ofs(work_dir + "/lambda_" + std::to_string(i) + ".dat");
      // max_digits10 round-trips a double exactly; the default 6 significant
      // digits silently truncated the Schmidt weights on every checkpoint.
      ofs << std::scientific
          << std::setprecision(std::numeric_limits<double>::max_digits10);
      for (int j = 0; j < nleg; ++j) {
        for (int k = 0; k < lattice.virtual_dims[i][j]; ++k) {
          ofs << lambda_tensor[i][j][k] << std::endl;
        }
      }
      ofs.flush();
      if (!ofs) {
        wrote = 0.0;
      }
    }
  }
  if (finfo.enabled) {
    // The virtual-bond parity ledger is mutable state (the simple update
    // rewrites it through svd_trunc), so it has to travel with the tensors:
    // reloading with a stale ledger changes the measured energy without any
    // error message.
    //
    // The CTM environment saved above is a placeholder, see the comment there.
    if (!save_fermion_parity(work_dir)) {
      wrote = 0.0;
    }
  }

  // Every file above is covered: the ones TeNeS writes by their streams, the
  // tensor files by saved_completely().
  std::vector<double> wrote_everywhere{wrote};
  allreduce_min(wrote_everywhere, comm);
  if (wrote_everywhere[0] == 0.0) {
    if (mpirank == 0) {
      util::remove_all(work_dir);
      std::cerr << "WARNING: failed to write the checkpoint files for "
                << save_dir << ", so no new checkpoint was saved." << std::endl;
    }
    return false;
  }

  /* Move the finished checkpoint into place, rank 0 alone: rename touches
   * metadata, so there is nothing for the other ranks to do and no race to
   * arrange. The marker goes in first and comes out after the sweep, so while
   * it exists the working directory does exist too -- which is what makes its
   * one repair instruction correct at every point the move can stop.
   */
  int moved = 0;
  std::string failed_name;
  if (mpirank == 0) {
    bool marker_written = false;
    {
      std::ofstream ofs(marker.c_str());
      ofs << "# TeNeS: while this file is here, the checkpoint in this "
             "directory\n"
             "# cannot be loaded. A save was interrupted before it finished "
             "moving\n"
             "# its files into place, so the directory may hold a mixture of "
             "two runs.\n"
             "#\n"
             "# To repair: move the contents of the working directory below "
             "into this\n"
             "# directory, then delete this file.\n"
             "#\n";
      // Absolute: the marker sits inside the destination, so a relative
      // spelling would read as relative to the destination itself.
      ofs << "work_dir = " << util::absolute_path(work_dir) << "\n";
      ofs << "started  = " << util::datetime() << "\n";
      ofs.flush();
      marker_written = static_cast<bool>(ofs);
    }

    std::set<std::string> written_names;
    if (!marker_written) {
      // Moving without the marker would leave nothing to say the destination
      // is half moved if this run died now. Nothing has moved yet, so stop
      // while the destination is still whole.
      util::remove_all(marker);
      util::remove_all(work_dir);
      moved = 2;
    } else if (move_all_into(work_dir, save_dir, written_names, failed_name)) {
      // Files of this run's own shape that it did not write are leftovers of
      // an earlier, differently shaped run. Everything else is the user's.
      for (const auto &path : util::entries_with_prefix(save_dir, "")) {
        const std::string name = util::basename(path);
        if (name == save_work_dir_name || name == save_marker_name ||
            written_names.count(name) != 0) {
          continue;
        }
        if (is_checkpoint_name(name)) {
          util::remove_all(path);
        }
      }
      if (util::remove_all(marker)) {
        util::remove_all(work_dir);
        moved = 1;
      } else {
        // Everything is in place, but the marker refuses the load. The
        // working directory, empty now, stays with it, so that the marker
        // never points at nothing.
        moved = 3;
      }
    }
  }
  bcast(moved, 0, comm);
  bcast(failed_name, 0, comm);

  if (moved != 1) {
    if (mpirank == 0) {
      if (moved == 2) {
        std::cerr << "WARNING: could not write " << marker
                  << ", so nothing was moved into " << save_dir
                  << " and no new checkpoint was saved." << std::endl;
      } else if (moved == 3) {
        std::cerr << "WARNING: moved the checkpoint into " << save_dir
                  << ", but could not delete " << marker
                  << ". The checkpoint is complete, but it cannot be loaded "
                     "until that file is deleted."
                  << std::endl;
      } else {
        std::cerr << "WARNING: could not move " << failed_name << " into "
                  << save_dir
                  << ", so that directory may now hold a mixture of this run "
                     "and the previous one, and its checkpoint cannot be "
                     "loaded while "
                  << marker << " is there.\n"
                  << "  To finish the move by hand, move everything left in "
                  << work_dir << " into " << save_dir << ", then delete "
                  << marker << ". A later save into " << save_dir
                  << " does the same before it writes anything." << std::endl;
      }
    }
    return false;
  }

  if (mpirank == 0 && peps_parameters.print_level >= PrintLevel::info) {
    std::cout << "Tensors saved in " << save_dir << std::endl;
  }
  return true;
}

template <class ptensor>
bool iTPS<ptensor>::save_fermion_parity(std::string const &save_dir) const {
  if (mpirank != 0) {
    return true;
  }
  std::string filename = save_dir + "/fermion.dat";
  std::ofstream ofs(filename.c_str());
  constexpr int fermion_format_version = 1;
  ofs << fermion_format_version << " # Fermion_Format_Version\n";
  ofs << N_UNIT << " # N_UNIT\n";
  ofs << lattice.LX << " " << lattice.LY << " # L_sub\n";
  ofs << lattice.skew << " # skew\n";
  auto write_parity = [&ofs](tenes::fermion::parity_vector const &p) {
    for (std::size_t i = 0; i < p.size(); ++i) {
      ofs << (p[i] ? 1 : 0) << " ";
    }
  };
  for (int i = 0; i < N_UNIT; ++i) {
    write_parity(finfo.phys[i]);
    ofs << "# parity of the physical leg of Tn[" << i << "]\n";
    for (int leg = 0; leg < nleg; ++leg) {
      write_parity(finfo.virt[i][leg]);
      ofs << "# parity of the virtual leg " << leg << " of Tn[" << i << "]\n";
    }
  }
  ofs.flush();
  return static_cast<bool>(ofs);
}

template <class ptensor>
void iTPS<ptensor>::load_tensors() {
  std::string const &load_dir = peps_parameters.tensor_load_dir;

  if (!util::isdir(load_dir)) {
    // "or cannot be read": the existence tests report false rather than
    // throwing, so that they cannot abort one rank inside a collective,
    // and a directory without search permission is then indistinguishable
    // from a missing one.
    std::string msg =
        load_dir + " does not exist, or is not a readable directory.";
    throw tenes::load_error(msg);
  }

  /* A save that stopped partway through moving its files leaves this behind.
   * The directory may hold a mixture of two runs, which would load without
   * complaint and give wrong physics, so refuse until a person has dealt with
   * it. The marker carries its own repair instructions; pass them on rather
   * than rebuilding them here.
   */
  const std::string marker = load_dir + "/" + save_marker_name;
  int interrupted = 0;
  std::string marker_body;
  if (mpirank == 0) {
    if (util::path_exists(marker)) {
      interrupted = 1;
      // The head only: the body is for a person to read, and it goes to every
      // rank through bcast(std::string&), whose length is an int.
      std::ifstream ifs(marker.c_str());
      marker_body.resize(marker_body_limit);
      ifs.read(&marker_body[0],
               static_cast<std::streamsize>(marker_body_limit));
      marker_body.resize(static_cast<std::size_t>(ifs.gcount()));
      if (ifs.good() && ifs.peek() != std::char_traits<char>::eof()) {
        marker_body += "\n[... the rest of this file is not shown ...]\n";
      }
    }
  }
  bcast(interrupted, 0, comm);
  bcast(marker_body, 0, comm);
  if (interrupted != 0) {
    std::stringstream ss;
    ss << "ERROR: " << marker
       << " exists: a save into this directory was interrupted before it "
          "finished moving its files into place, so the checkpoint here may "
          "be a mixture of two runs.\n";
    if (!marker_body.empty()) {
      ss << marker_body;
    } else {
      ss << "HINT: move the contents of " << load_dir << "/"
         << save_work_dir_name << " into " << load_dir << ", then delete "
         << marker << ".";
    }
    throw tenes::load_error(ss.str());
  }

  // No params.dat at all is how the pre-versioning format announces itself;
  // one that is there but whose first line is not a number is a damaged
  // checkpoint, and saying so beats std::stoi's bare "stoi". The status
  // travels with the version so that every rank decides the same way: a
  // throw on rank 0 alone would leave the others in the broadcast below.
  const std::string params_file = load_dir + "/params.dat";
  int tensor_format_version = 0;
  int version_readable = 1;
  if (mpirank == 0) {
    if (util::path_exists(params_file)) {
      std::ifstream ifs(params_file.c_str());
      std::string line;
      version_readable = 0;
      if (std::getline(ifs, line)) {
        try {
          tensor_format_version = std::stoi(util::drop_comment(line));
          version_readable = 1;
        } catch (const std::exception &) {
          version_readable = 0;
        }
      }
    }
  }
  bcast(tensor_format_version, 0, comm);
  bcast(version_readable, 0, comm);
  if (version_readable == 0) {
    throw tenes::load_error(
        "ERROR: cannot read the format version from " + params_file +
        ".\n"
        "HINT: the file is empty or damaged; its first line has to be the "
        "saved tensor format version.");
  }

  if (tensor_format_version == 0) {
    std::vector<std::vector<int>> current_shape(N_UNIT,
                                                std::vector<int>(nleg + 1));
    for (int i = 0; i < N_UNIT; ++i) {
      for (int leg = 0; leg < nleg; ++leg) {
        current_shape[i][leg] = lattice.virtual_dims[i][leg];
      }
      current_shape[i][nleg] = lattice.physical_dims[i];
    }
    load_fermion_ledger(load_dir, current_shape, false, std::nullopt);
    load_tensors_v0();
  } else if (tensor_format_version == 1 || tensor_format_version == 2) {
    load_tensors_versioned(tensor_format_version);
  } else {
    std::stringstream ss;
    ss << "ERROR: " << params_file << " has an unknown saved tensor format "
       << "version: " << tensor_format_version;
    throw tenes::load_error(ss.str());
  }

  validate_loaded_fermion_tensors();
  ctm_valid_ = false;
}

template <class ptensor>
void iTPS<ptensor>::load_fermion_ledger(
    std::string const &load_dir,
    std::vector<std::vector<int>> const &saved_shape, bool validate_saved_shape,
    std::optional<bool> recorded_kind) {
  const std::string filename = load_dir + "/fermion.dat";

  if (recorded_kind.has_value()) {
    // The checkpoint says what it is, so the ledger file's presence decides
    // nothing: a stale one beside a non-fermionic checkpoint is just litter.
    if (*recorded_kind != finfo.enabled) {
      const std::string params_file = load_dir + "/params.dat";
      throw tenes::load_error("ERROR: " + params_file +
                              " records that the saved tensors come "
                              "from a " +
                              (*recorded_kind ? "fermionic" : "non-fermionic") +
                              " run, but parameter.general.fermion is " +
                              (finfo.enabled ? "true" : "false") +
                              ".\n"
                              "HINT: set fermion = " +
                              (*recorded_kind ? "true" : "false") +
                              ", or load a checkpoint saved by a " +
                              (finfo.enabled ? "fermionic" : "non-fermionic") +
                              " run.");
    }
    if (!finfo.enabled) {
      return;
    }
  }

  int exists = 0;
  if (mpirank == 0) {
    exists = util::path_exists(filename) ? 1 : 0;
  }
  bcast(exists, 0, comm);

  if (!finfo.enabled) {
    if (exists != 0) {
      throw tenes::load_error(
          "ERROR: " + filename +
          " exists, i.e. the saved tensors come from a fermionic run, but "
          "parameter.general.fermion is false.\n"
          "HINT: set fermion = true, or load tensors saved by a "
          "non-fermionic run.");
    }
    return;
  }
  if (exists == 0) {
    // "or cannot be read": path_exists() answers false rather than throwing
    // when the filesystem cannot say, so an unreadable directory arrives here
    // looking exactly like a missing file. Naming a cause outright -- the
    // earlier wording asserted the tensors came from a non-fermionic run --
    // then misdiagnoses a permission problem.
    throw tenes::load_error(
        "ERROR: cannot read " + filename +
        ".\n"
        "Fermion mode needs the fermionic parity ledger of the virtual bonds "
        "in order to interpret the saved tensors.\n"
        "HINT: the file is missing, because the tensors were saved by a "
        "non-fermionic run or by a version of TeNeS that could not save "
        "fermionic tensors; or it is there but cannot be read, because its "
        "directory cannot be searched.");
  }

  std::string content;
  if (mpirank == 0) {
    std::ifstream ifs(filename.c_str());
    std::stringstream ss;
    ss << ifs.rdbuf();
    content = ss.str();
  }
  bcast(content, 0, comm);

  metadata_reader reader(content, filename);
  auto next_ints = [&reader]() { return reader.next_ints(); };
  auto next_scalar = [&reader]() { return reader.next_scalar(); };

  const int version = next_scalar();
  if (version != 1) {
    std::stringstream ss;
    ss << "ERROR: " << filename << " has fermion format version " << version
       << " but this version of TeNeS supports only version 1";
    throw tenes::load_error(ss.str());
  }
  const int loaded_N_UNIT = next_scalar();
  if (loaded_N_UNIT != N_UNIT) {
    std::stringstream ss;
    ss << "ERROR: N_UNIT is " << N_UNIT << " but " << filename << " has "
       << loaded_N_UNIT;
    throw tenes::load_error(ss.str());
  }
  const auto lsub = next_ints();
  const int loaded_skew = next_scalar();
  if (lsub.size() != 2 || lsub[0] != lattice.LX || lsub[1] != lattice.LY ||
      loaded_skew != lattice.skew) {
    std::stringstream ss;
    ss << "ERROR: the unit cell of the saved tensors (L_sub = [";
    for (std::size_t i = 0; i < lsub.size(); ++i) {
      ss << (i == 0 ? "" : ", ") << lsub[i];
    }
    ss << "], skew = " << loaded_skew << ") differs from the input (L_sub = ["
       << lattice.LX << ", " << lattice.LY << "], skew = " << lattice.skew
       << ").\n"
       << "HINT: the parity ledger is indexed by (site, leg), so it only means "
          "the same thing on the same lattice.";
    throw tenes::load_error(ss.str());
  }

  auto to_parity = [](std::vector<int> const &v) {
    tenes::fermion::parity_vector p(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
      p[i] = (v[i] != 0);
    }
    return p;
  };

  std::vector<std::array<tenes::fermion::parity_vector, 4>> virt(N_UNIT);
  for (int i = 0; i < N_UNIT; ++i) {
    const auto phys = to_parity(next_ints());
    if (phys != finfo.phys[i]) {
      std::stringstream ss;
      ss << "ERROR: the physical parity of the tensor " << i << " in "
         << filename << " differs from tensor.unitcell.parity in the input";
      throw tenes::load_error(ss.str());
    }
    for (int leg = 0; leg < nleg; ++leg) {
      auto p = to_parity(next_ints());
      const int saved_dim = saved_shape[i][leg];
      if (!validate_saved_shape &&
          lattice.virtual_dims[i][leg] != static_cast<int>(p.size())) {
        std::stringstream ss;
        ss << "ERROR: the virtual dimension of the leg " << leg
           << " of the tensor " << i << " is " << lattice.virtual_dims[i][leg]
           << " but the saved tensors have " << p.size() << ".\n"
           << "HINT: legacy fermion tensor checkpoints cannot change "
              "virtual_dim on restart. Keep virtual_dim as it was, or start "
              "without tensor_load.";
        throw tenes::load_error(ss.str());
      }
      if (validate_saved_shape && static_cast<int>(p.size()) != saved_dim) {
        std::stringstream ss;
        ss << "ERROR: the virtual parity ledger of the leg " << leg
           << " of the tensor " << i << " in " << filename << " has "
           << p.size() << " entries but the saved tensor has dimension "
           << saved_dim << ".\n"
           << "HINT: the checkpoint is inconsistent; fermion.dat and "
              "params.dat do not describe the same saved tensors.";
        throw tenes::load_error(ss.str());
      }
      const int loaded_dim =
          validate_saved_shape ? saved_dim : static_cast<int>(p.size());
      if (lattice.virtual_dims[i][leg] < loaded_dim) {
        std::stringstream ss;
        ss << "ERROR: the virtual dimension of the leg " << leg
           << " of the tensor " << i << " is " << lattice.virtual_dims[i][leg]
           << " but the saved tensors have " << loaded_dim << ".\n"
           << "HINT: fermion mode can expand virtual_dim on restart, but "
              "cannot shrink it. Keep virtual_dim at least as large as the "
              "saved tensors, or start without tensor_load.";
        throw tenes::load_error(ss.str());
      }
      p = tenes::fermion::extend_parity(
          p, static_cast<std::size_t>(lattice.virtual_dims[i][leg]));
      virt[i][leg] = p;
    }
  }
  for (int i = 0; i < N_UNIT; ++i) {
    finfo.virt[i] = virt[i];
  }
  try {
    tenes::fermion::validate_neighbor_consistency(finfo, lattice);
  } catch (const std::exception &e) {
    throw tenes::load_error(
        std::string(e.what()) + "\nHINT: the parity ledger in " + filename +
        " disagrees between the two ends of a bond, so it does not describe "
        "the lattice being loaded into.");
  }
}

template <class ptensor>
void iTPS<ptensor>::validate_loaded_fermion_tensors() const {
  if (!finfo.enabled) {
    return;
  }
  for (int i = 0; i < N_UNIT; ++i) {
    const auto ft = tenes::fermion::wrap_Tn(Tn[i], finfo, i);
    // parity_violation scans the process-local slice; max_abs may already
    // reduce internally, but the extra allreduce is harmless.
    std::vector<double> reduced{tenes::fermion::parity_violation(ft),
                                tenes::fermion::max_abs(ft)};
    tenes::allreduce_max(reduced, comm);
    const double violation = reduced[0];
    const double scale = std::max(reduced[1], 1.0e-300);
    if (violation > 1.0e-12 * scale) {
      std::stringstream ss;
      ss << "ERROR: the loaded tensor " << i
         << " breaks fermion parity under the loaded parity ledger "
         << "(max violating amplitude " << violation << ", max amplitude "
         << reduced[1] << ").\n"
         << "HINT: the tensors and " << peps_parameters.tensor_load_dir
         << "/fermion.dat do not belong together.";
      throw tenes::load_error(ss.str());
    }
  }
}

template <class ptensor>
void load_tensor(ptensor &A, std::string const &name,
                 std::string const &directory, int iunit) {
  std::string filename =
      directory + "/" + name + "_" + std::to_string(iunit) + ".dat";
  if (!util::path_exists(filename)) {
    throw tenes::load_error("ERROR: cannot read a tensor file: " + filename +
                            " (it is missing, or its directory cannot be "
                            "searched)");
  }
  // On A's communicator: the default constructor would put temp on
  // MPI_COMM_WORLD, and load() does not change it, so a library caller
  // running the solver on a subset of the processes would get every loaded
  // tensor distributed over all of them.
  ptensor temp(A.get_comm());
  temp.load(filename.c_str());
  if (A.rank() != temp.rank()) {
    std::stringstream ss;
    ss << "ERROR: rank mismatch in load_tensor: ";
    ss << name << "[" << iunit << "] has " << A.rank() << " legs, but ";
    ss << "loaded one has " << temp.rank() << " legs." << std::endl;
    ss << "HINT: check the calculation mode. The number of legs differs "
          "between ground state calculation and finite temperature "
          "calculation.";
    throw tenes::load_error(ss.str());
  }
  A = resize_tensor(temp, A.shape());
}

template <class ptensor>
void iTPS<ptensor>::load_tensors_versioned(int expected_version) {
  std::string const &load_dir = peps_parameters.tensor_load_dir;
  const std::string params_file = load_dir + "/params.dat";

  std::string content;
  if (mpirank == 0) {
    std::ifstream ifs(params_file.c_str());
    std::stringstream ss;
    ss << ifs.rdbuf();
    content = ss.str();
  }
  bcast(content, 0, comm);
  metadata_reader reader(content, params_file);

  const int format_version = reader.next_scalar();
  if (format_version != expected_version) {
    std::stringstream ss;
    ss << "ERROR: " << params_file << " has format version " << format_version
       << " but was dispatched as version " << expected_version;
    throw tenes::load_error(ss.str());
  }

  // Version 2 records the kind of run. Before it, the loader had to infer that
  // from the presence of fermion.dat.
  std::optional<bool> recorded_kind;
  if (format_version >= 2) {
    const int kind = reader.next_scalar();
    if (kind != 0 && kind != 1) {
      std::stringstream ss;
      ss << "ERROR: " << params_file << " records the run kind as " << kind
         << "; it has to be 0 (not fermionic) or 1 (fermionic)";
      throw tenes::load_error(ss.str());
    }
    recorded_kind = (kind != 0);
  }

  // Before the shapes, not after: there is one shape line per saved site, so
  // reading N_UNIT of them out of a checkpoint that holds fewer is a read
  // past the end of the file.
  const int loaded_N_UNIT = reader.next_scalar();
  if (N_UNIT != loaded_N_UNIT) {
    std::stringstream ss;
    ss << "ERROR: N_UNIT is " << N_UNIT << " but " << params_file
       << " was saved with N_UNIT = " << loaded_N_UNIT;
    throw tenes::load_error(ss.str());
  }

  const int loaded_CHI = reader.next_scalar();
  if (CHI != static_cast<std::size_t>(loaded_CHI)) {
    if (mpirank == 0 && peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "WARNING: parameters.ctm.dimension is " << CHI
                << " but loaded tensors have CHI = " << loaded_CHI << std::endl;
    }
  }

  std::vector<std::vector<int>> loaded_shape(N_UNIT,
                                             std::vector<int>(nleg + 1));
  for (int i = 0; i < N_UNIT; ++i) {
    const auto shape = reader.next_ints();
    if (shape.size() < static_cast<std::size_t>(nleg) + 1) {
      std::stringstream ss;
      ss << "ERROR: " << params_file << " gives " << shape.size()
         << " dimensions for the tensor " << i << ", but " << (nleg + 1)
         << " are needed";
      throw tenes::load_error(ss.str());
    }
    for (int j = 0; j < nleg; ++j) {
      loaded_shape[i][j] = shape[j];
      const int vd_param = lattice.virtual_dims[i][j];
      if (vd_param != loaded_shape[i][j]) {
        if (mpirank == 0 && peps_parameters.print_level >= PrintLevel::info) {
          std::cout << "WARNING: virtual dimension of the leg " << j
                    << " of the tensor " << i << " is " << vd_param
                    << " but loaded tensor has " << loaded_shape[i][j]
                    << std::endl;
        }
      }
    }
    loaded_shape[i][nleg] = shape[nleg];
    const int pdim = lattice.physical_dims[i];
    if (pdim != loaded_shape[i][nleg]) {
      std::stringstream ss;
      ss << "ERROR: dimension of the physical bond of the tensor " << i
         << " is " << pdim << " but " << params_file << " has "
         << loaded_shape[i][nleg];
      throw tenes::load_error(ss.str());
    }
  }

  load_fermion_ledger(load_dir, loaded_shape, true, recorded_kind);

  // #define LOAD_TENSOR_(A, name)                      \
  //   do {                                             \
  //     ptensor temp;                                  \
  //     temp.load((filename + name + suffix).c_str()); \
  //     A = resize_tensor(temp, A.shape());            \
  //   } while (false)

  for (int i = 0; i < N_UNIT; ++i) {
    std::string filename = load_dir + "/";
    std::string suffix = "_" + std::to_string(i) + ".dat";

    load_tensor(Tn[i], "T", load_dir, i);
    load_tensor(eTl[i], "El", load_dir, i);
    load_tensor(eTt[i], "Et", load_dir, i);
    load_tensor(eTr[i], "Er", load_dir, i);
    load_tensor(eTb[i], "Eb", load_dir, i);
    load_tensor(C1[i], "C1", load_dir, i);
    load_tensor(C2[i], "C2", load_dir, i);
    load_tensor(C3[i], "C3", load_dir, i);
    load_tensor(C4[i], "C4", load_dir, i);
    // LOAD_TENSOR_(Tn[i], "T");
    // LOAD_TENSOR_(eTl[i], "El");
    // LOAD_TENSOR_(eTt[i], "Et");
    // LOAD_TENSOR_(eTr[i], "Er");
    // LOAD_TENSOR_(eTb[i], "Eb");
    // LOAD_TENSOR_(C1[i], "C1");
    // LOAD_TENSOR_(C2[i], "C2");
    // LOAD_TENSOR_(C3[i], "C3");
    // LOAD_TENSOR_(C4[i], "C4");
  }
  // #undef LOAD_TENSOR_

  std::vector<double> ls;
  std::string lambda_error;
  // Rank 0 reads and the others take the result, so a read failure has to
  // travel to them as data. Throwing here on rank 0 alone would leave every
  // other rank in the broadcast below.
  if (mpirank == 0) {
    lambda_error = [&]() -> std::string {
      for (int i = 0; i < N_UNIT; ++i) {
        std::string lambda_filename =
            load_dir + "/lambda_" + std::to_string(i) + ".dat";
        std::ifstream ifs(lambda_filename.c_str());
        for (int j = 0; j < nleg; ++j) {
          for (int k = 0; k < loaded_shape[i][j]; ++k) {
            double temp = 0.0;
            if (!(ifs >> temp)) {
              return "ERROR: failed to read lambda values from " +
                     lambda_filename;
            }
            ls.push_back(temp);
          }
        }
      }
      return std::string();
    }();
  }
  bcast(lambda_error, 0, comm);
  if (!lambda_error.empty()) {
    throw tenes::load_error(lambda_error);
  }
  bcast(ls, 0, comm);
  int index = 0;
  for (int i = 0; i < N_UNIT; ++i) {
    const auto vdim = lattice.virtual_dims[i];
    for (int j = 0; j < nleg; ++j) {
      lambda_tensor[i][j].clear();
      for (int k = 0; k < loaded_shape[i][j]; ++k) {
        lambda_tensor[i][j].push_back(ls[index]);
        ++index;
      }
      lambda_tensor[i][j].resize(vdim[j]);
    }
  }
}

template <class ptensor>
void iTPS<ptensor>::load_tensors_v0() {
  using mptensor::Shape;
  std::string const &load_dir = peps_parameters.tensor_load_dir;

  // load from the checkpoint
  if (!util::isdir(load_dir)) {
    // "or cannot be read": the existence tests report false rather than
    // throwing, so that they cannot abort one rank inside a collective,
    // and a directory without search permission is then indistinguishable
    // from a missing one.
    std::string msg =
        load_dir + " does not exist, or is not a readable directory.";
    throw tenes::load_error(msg);
  }
  for (int i = 0; i < N_UNIT; ++i) {
    std::string filename = load_dir + "/";
    std::string suffix = "_" + std::to_string(i) + ".dat";
    auto load = [&filename, &suffix](ptensor &A, const char *name) {
      std::string path = filename + name + suffix;
      if (!util::path_exists(path)) {
        throw tenes::load_error(
            "ERROR: cannot read a tensor file: " + path +
            " (it is missing, or its directory cannot be searched)");
      }
      A.load(path.c_str());
    };
    load(Tn[i], "T");
    load(eTt[i], "Et");
    load(eTr[i], "Er");
    load(eTb[i], "Eb");
    load(eTl[i], "El");
    load(C1[i], "C1");
    load(C2[i], "C2");
    load(C3[i], "C3");
    load(C4[i], "C4");
  }
  std::vector<double> ls;
  std::string lambda_error;
  // Rank 0 reads and the others take the result, so a read failure has to
  // travel to them as data. Throwing here on rank 0 alone would leave every
  // other rank in the broadcast below.
  if (mpirank == 0) {
    lambda_error = [&]() -> std::string {
      for (int i = 0; i < N_UNIT; ++i) {
        const auto vdim = lattice.virtual_dims[i];
        std::string lambda_filename =
            load_dir + "/lambda_" + std::to_string(i) + ".dat";
        std::ifstream ifs(lambda_filename.c_str());
        for (int j = 0; j < nleg; ++j) {
          for (int k = 0; k < vdim[j]; ++k) {
            double temp = 0.0;
            if (!(ifs >> temp)) {
              return "ERROR: failed to read lambda values from " +
                     lambda_filename;
            }
            ls.push_back(temp);
          }
        }
      }
      return std::string();
    }();
  }
  bcast(lambda_error, 0, comm);
  if (!lambda_error.empty()) {
    throw tenes::load_error(lambda_error);
  }
  bcast(ls, 0, comm);
  int index = 0;
  for (int i = 0; i < N_UNIT; ++i) {
    const auto vdim = lattice.virtual_dims[i];
    for (int j = 0; j < nleg; ++j) {
      for (int k = 0; k < vdim[j]; ++k) {
        lambda_tensor[i][j][k] = ls[index];
        ++index;
      }
    }
  }

  // overwrite dimensions
  const Shape Cshape = C1[0].shape();
  if (CHI != Cshape[0]) {
    if (peps_parameters.print_level >= PrintLevel::info) {
      std::cout << "WARNING: parameters.ctm.dimension is " << CHI
                << " but loaded tensors have CHI = " << Cshape[0] << std::endl;
    }
  }
  for (int i = 0; i < N_UNIT; ++i) {
    const Shape Tshape = Tn[i].shape();
    const int pdim = lattice.physical_dims[i];
    if (static_cast<std::size_t>(pdim) != Tshape[4]) {
      std::stringstream ss;
      ss << "ERROR: dimension of the physical bond of the tensor " << i
         << " is " << pdim << " but loaded tensor has " << Tshape[4]
         << std::endl;
      throw tenes::input_error(ss.str());
    }

    for (int l = 0; l < nleg; ++l) {
      const int vd_param = lattice.virtual_dims[i][l];
      const int vd_loaded = Tshape[l];
      if (vd_param != vd_loaded) {
        if (peps_parameters.print_level >= PrintLevel::info) {
          std::cout << "WARNING: virtual dimension of the leg " << l
                    << " of the tensor " << i << " is " << vd_param
                    << " but loaded tensor has " << vd_loaded << std::endl;
        }
      }
    }
  }
}

// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace tenes::itps
