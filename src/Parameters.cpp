#include <mpi.h>
#include <toml11/toml.hpp>
#include "util/find_or.hpp"
#include "Parameters.hpp"

Parameters::Parameters(const char *filename): Parameters(toml::parse(filename)){}
Parameters::Parameters(toml::Table data){
  set(data);
}

void Parameters::set(const char *filename){ set(toml::parse(filename)); }
void Parameters::set(toml::Table data){
  toml::Table param = toml::find<toml::Table>(data, "parameter");

  tau_simple = toml::find<double>(param, "tau_simple");
  num_simple_step = toml::find<int>(param, "num_simple_step");

  tau_full = toml::find<double>(param, "tau_full");
  num_full_step = util::find_or(param, "num_full_step", 0);
}

