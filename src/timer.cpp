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

#include "timer.hpp"

#include <iomanip>
#include <sstream>
#include <vector>

#include "mpi.hpp"

namespace tenes {

std::map<std::string, TimerAggregate> aggregate_timers(
    TimerRegistry const &registry, MPI_Comm comm) {
  std::vector<double> sums;
  sums.reserve(registry.entries().size());
  for (auto const &kv : registry.entries()) {
    sums.push_back(kv.second.sum);
  }
  std::vector<double> maxs = sums;
  std::vector<double> mins = sums;
  allreduce_max(maxs, comm);
  allreduce_min(mins, comm);

  std::map<std::string, TimerAggregate> result;
  std::size_t i = 0;
  for (auto const &kv : registry.entries()) {
    result[kv.first] =
        TimerAggregate{kv.second.count, sums[i], maxs[i], mins[i]};
    ++i;
  }
  return result;
}

namespace {
std::string escape_json(std::string const &s) {
  std::string ret;
  for (char c : s) {
    if (c == '"' || c == '\\') {
      ret += '\\';
    }
    ret += c;
  }
  return ret;
}
}  // namespace

std::string timers_to_json(std::map<std::string, TimerAggregate> const &timers,
                           std::string const &tenes_version, int mpi_size,
                           int omp_threads) {
  std::ostringstream oss;
  oss << std::setprecision(17);
  oss << "{\n";
  oss << "  \"meta\": {\n";
  oss << "    \"tenes_version\": \"" << escape_json(tenes_version) << "\",\n";
  oss << "    \"mpi_size\": " << mpi_size << ",\n";
  oss << "    \"omp_threads\": " << omp_threads << "\n";
  oss << "  },\n";
  oss << "  \"timers\": {\n";
  bool first = true;
  for (auto const &kv : timers) {
    if (!first) {
      oss << ",\n";
    }
    first = false;
    oss << "    \"" << escape_json(kv.first) << "\": {"
        << "\"count\": " << kv.second.count << ", \"sum\": " << kv.second.sum
        << ", \"max_rank\": " << kv.second.max_rank
        << ", \"min_rank\": " << kv.second.min_rank << "}";
  }
  oss << "\n  }\n}\n";
  return oss.str();
}

}  // namespace tenes
