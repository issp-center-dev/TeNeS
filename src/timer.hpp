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

#ifndef TENES_SRC_TIMER_HPP_
#define TENES_SRC_TIMER_HPP_

#include <chrono>
#include <cstddef>
#include <map>
#include <string>
#include <utility>

#include "mpi.hpp"

namespace tenes {

template <class CT = std::chrono::high_resolution_clock>
class Timer {
 public:
  using clock_type = CT;

  Timer() : start(clock_type::now()) {}
  void reset() { start = clock_type::now(); }
  double elapsed() const {
    const auto elapsed_time =
        std::chrono::duration_cast<std::chrono::nanoseconds>(clock_type::now() -
                                                             start);
    return elapsed_time.count() * 1.0e-9;
  }

 private:
  typename clock_type::time_point start;
};

//! @brief Process-wide accumulator of named timers
//!
//! Names are hierarchical, slash-separated (e.g. "contract/itps_ctm/2x2").
class TimerRegistry {
 public:
  struct Entry {
    std::size_t count = 0;
    double sum = 0.0;  //!< accumulated time in seconds
  };

  void add(std::string const &name, double seconds) {
    Entry &e = entries_[name];
    e.count += 1;
    e.sum += seconds;
  }

  std::map<std::string, Entry> const &entries() const { return entries_; }

  //! Registry used by instrumentation points across the process.
  static TimerRegistry &instance() {
    static TimerRegistry registry;
    return registry;
  }

 private:
  std::map<std::string, Entry> entries_;
};

//! RAII helper: measures from construction to destruction.
class ScopedTimer {
 public:
  explicit ScopedTimer(std::string name,
                       TimerRegistry &registry = TimerRegistry::instance())
      : name_(std::move(name)), registry_(registry) {}
  ~ScopedTimer() { registry_.add(name_, timer_.elapsed()); }
  ScopedTimer(ScopedTimer const &) = delete;
  ScopedTimer &operator=(ScopedTimer const &) = delete;

 private:
  std::string name_;
  TimerRegistry &registry_;
  Timer<> timer_;
};

//! Cross-rank view of one timer (times in seconds).
struct TimerAggregate {
  std::size_t count = 0;
  double sum = 0.0;       //!< value on the calling rank
  double max_rank = 0.0;  //!< maximum over MPI ranks
  double min_rank = 0.0;  //!< minimum over MPI ranks
};

/*! @brief Aggregate timer sums across MPI ranks (collective call)
 *
 * Assumes every rank recorded the same set of names
 * (tensor operations are collective).
 */
std::map<std::string, TimerAggregate> aggregate_timers(
    TimerRegistry const &registry, MPI_Comm comm);

//! Render aggregated timers and metadata as a JSON document.
std::string timers_to_json(std::map<std::string, TimerAggregate> const &timers,
                           std::string const &tenes_version, int mpi_size,
                           int omp_threads);

}  // end of namespace tenes

#endif  // TENES_SRC_TIMER_HPP_
