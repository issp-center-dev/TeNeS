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
 *  @brief Exception types thrown by TeNeS, so callers can tell input
 *         mistakes from internal bugs.
 */

#ifndef TENES_SRC_EXCEPTION_HPP_
#define TENES_SRC_EXCEPTION_HPP_

#include <stdexcept>

namespace tenes {

//! Internal inconsistency: a bug in TeNeS, not in the input.
class logic_error : public std::logic_error {
 public:
  //! Construct with the given message.
  explicit logic_error(const std::string &what_arg)
      : std::logic_error(what_arg) {}
  //! @copybrief logic_error(const std::string &)
  explicit logic_error(const char *what_arg) : std::logic_error(what_arg) {}
};
//! A requested combination of features is not implemented.
class unimplemented_error : public tenes::logic_error {
 public:
  //! Construct with the given message.
  explicit unimplemented_error(const std::string &what_arg)
      : tenes::logic_error(what_arg) {}
  //! @copybrief unimplemented_error(const std::string &)
  explicit unimplemented_error(const char *what_arg)
      : tenes::logic_error(what_arg) {}
};

//! Failure at run time whose cause is outside TeNeS itself.
class runtime_error : public std::runtime_error {
 public:
  //! Construct with the given message.
  explicit runtime_error(const std::string &what_arg)
      : std::runtime_error(what_arg) {}
  //! @copybrief runtime_error(const std::string &)
  explicit runtime_error(const char *what_arg) : std::runtime_error(what_arg) {}
};
//! The input file is invalid; the message explains what to fix.
class input_error : public runtime_error {
 public:
  //! Construct with the given message.
  explicit input_error(const std::string &what_arg) : runtime_error(what_arg) {}
  //! @copybrief input_error(const std::string &)
  explicit input_error(const char *what_arg) : runtime_error(what_arg) {}
};
//! Saved tensors could not be loaded (missing, corrupt, or mismatched).
class load_error : public runtime_error {
 public:
  //! Construct with the given message.
  explicit load_error(const std::string &what_arg) : runtime_error(what_arg) {}
  //! @copybrief load_error(const std::string &)
  explicit load_error(const char *what_arg) : runtime_error(what_arg) {}
};

}  // end of namespace tenes

#endif  // TENES_SRC_EXCEPTION_HPP_
