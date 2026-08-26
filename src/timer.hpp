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
 *  @brief Wall-clock timing: the ::tenes::Timer stopwatch, the named
 *         ::tenes::TimerRegistry the benchmark instrumentation reports
 *         into, and the MPI aggregation that renders timers.json.
 */

#ifndef TENES_SRC_TIMER_HPP_
#define TENES_SRC_TIMER_HPP_

#include <chrono>
#include <cstddef>
#include <map>
#include <string>
#include <utility>

#include "mpi.hpp"

namespace tenes {

//! Stopwatch: measures wall-clock seconds since construction or reset().
template <class CT = std::chrono::high_resolution_clock>
class Timer {
 public:
  using clock_type = CT;  //!< clock the stopwatch reads

  //! Start measuring now.
  Timer() : start(clock_type::now()) {}
  //! Restart measuring from now.
  void reset() { start = clock_type::now(); }
  //! Seconds since construction or the last reset().
  double elapsed() const {
    const auto elapsed_time =
        std::chrono::duration_cast<std::chrono::nanoseconds>(clock_type::now() -
                                                             start);
    return elapsed_time.count() * 1.0e-9;
  }

 private:
  typename clock_type::time_point start;  //!< when measuring began
};

//! @brief Process-wide accumulator of named timers
//!
//! Names are hierarchical, slash-separated (e.g. "contract/itps_ctm/2x2").
class TimerRegistry {
 public:
  //! Accumulated state of one named timer.
  struct Entry {
    std::size_t count = 0;  //!< number of recorded measurements
    double sum = 0.0;       //!< accumulated time in seconds
  };

  //! Record one measurement under the given name.
  void add(std::string const &name, double seconds) {
    Entry &e = entries_[name];
    e.count += 1;
    e.sum += seconds;
  }

  //! All recorded timers, keyed by name.
  std::map<std::string, Entry> const &entries() const { return entries_; }

  //! Registry used by instrumentation points across the process.
  static TimerRegistry &instance() {
    static TimerRegistry registry;
    return registry;
  }

 private:
  std::map<std::string, Entry> entries_;  //!< recorded timers by name
};

//! RAII helper: measures from construction to destruction.
class ScopedTimer {
 public:
  //! Start timing; the destructor records the elapsed time under name.
  explicit ScopedTimer(std::string name,
                       TimerRegistry &registry = TimerRegistry::instance())
      : name_(std::move(name)), registry_(registry) {}
  ~ScopedTimer() { registry_.add(name_, timer_.elapsed()); }
  ScopedTimer(ScopedTimer const &) = delete;
  ScopedTimer &operator=(ScopedTimer const &) = delete;

 private:
  std::string name_;         //!< timer name to record under
  TimerRegistry &registry_;  //!< registry receiving the measurement
  Timer<> timer_;            //!< the running stopwatch
};

//! Cross-rank view of one timer (times in seconds).
struct TimerAggregate {
  std::size_t count = 0;  //!< number of measurements on the calling rank
  double sum = 0.0;       //!< value on the calling rank
  double max_rank = 0.0;  //!< maximum over MPI ranks
  double min_rank = 0.0;  //!< minimum over MPI ranks
};

/*! @brief Aggregate timer sums across MPI ranks (collective call)
 *
 * Assumes every rank recorded the same set of names
 * (tensor operations are collective).
 * If ranks disagree on the recorded names, falls back to local values
 * (max_rank == min_rank == sum) with a warning on rank 0.
 */
std::map<std::string, TimerAggregate> aggregate_timers(
    TimerRegistry const &registry, MPI_Comm comm);

//! Render aggregated timers and metadata as a JSON document.
std::string timers_to_json(std::map<std::string, TimerAggregate> const &timers,
                           std::string const &tenes_version, int mpi_size,
                           int omp_threads);

}  // end of namespace tenes

#endif  // TENES_SRC_TIMER_HPP_
