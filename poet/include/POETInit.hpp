#pragma once

#include <string>
#include <vector>

struct POETInitCell {
  std::vector<std::string> solutions;
  std::vector<std::string> solution_primaries;
  std::vector<std::string> exchanger;
  std::vector<std::string> kinetics;
  std::vector<std::string> equilibrium;
  std::vector<std::string> surface_comps;
  std::vector<std::string> surface_charges;
};

struct POETConfig {
  std::string database;
  std::string input_script;
  POETInitCell cell;
};