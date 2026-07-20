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

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <string>

#include "../src/mpi.hpp"
#include "../src/timer.hpp"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  doctest::Context context(argc, argv);
  const int res = context.run();
  MPI_Finalize();
  return res;
}

TEST_CASE("TimerRegistry::add accumulates count and sum") {
  tenes::TimerRegistry registry;
  registry.add("contract/itps_ctm/2x2", 1.5);
  registry.add("contract/itps_ctm/2x2", 0.5);
  registry.add("phase/simple_update", 2.0);

  auto const &entries = registry.entries();
  REQUIRE(entries.size() == 2);
  CHECK(entries.at("contract/itps_ctm/2x2").count == 2);
  CHECK(entries.at("contract/itps_ctm/2x2").sum == doctest::Approx(2.0));
  CHECK(entries.at("phase/simple_update").count == 1);
  CHECK(entries.at("phase/simple_update").sum == doctest::Approx(2.0));
}

TEST_CASE("ScopedTimer records on destruction") {
  tenes::TimerRegistry registry;
  {
    tenes::ScopedTimer timer("scope/test", registry);
  }
  auto const &entries = registry.entries();
  REQUIRE(entries.count("scope/test") == 1);
  CHECK(entries.at("scope/test").count == 1);
  CHECK(entries.at("scope/test").sum >= 0.0);
}

TEST_CASE("TimerRegistry::instance returns a singleton") {
  auto &a = tenes::TimerRegistry::instance();
  auto &b = tenes::TimerRegistry::instance();
  CHECK(&a == &b);
}

TEST_CASE("aggregate_timers with a single process") {
  tenes::TimerRegistry registry;
  registry.add("a", 1.0);
  registry.add("a", 2.0);
  registry.add("b", 4.0);

  auto agg = tenes::aggregate_timers(registry, MPI_COMM_WORLD);
  REQUIRE(agg.size() == 2);
  CHECK(agg.at("a").count == 2);
  CHECK(agg.at("a").sum == doctest::Approx(3.0));
  CHECK(agg.at("a").max_rank == doctest::Approx(3.0));
  CHECK(agg.at("a").min_rank == doctest::Approx(3.0));
  CHECK(agg.at("b").sum == doctest::Approx(4.0));
}

TEST_CASE("timers_to_json renders the expected fields") {
  std::map<std::string, tenes::TimerAggregate> timers;
  timers["total"] = tenes::TimerAggregate{1, 1.5, 1.5, 1.5};
  timers["contract/itps_ctm/2x2"] = tenes::TimerAggregate{48, 4.5, 5.0, 4.0};
  const std::string json = tenes::timers_to_json(timers, "2.2-dev", 1, 8);
  CHECK(json.find("\"tenes_version\": \"2.2-dev\"") != std::string::npos);
  CHECK(json.find("\"mpi_size\": 1") != std::string::npos);
  CHECK(json.find("\"omp_threads\": 8") != std::string::npos);
  CHECK(json.find("\"total\": {\"count\": 1, \"sum\": 1.5") !=
        std::string::npos);
  CHECK(json.find("\"contract/itps_ctm/2x2\": {\"count\": 48") !=
        std::string::npos);
}
