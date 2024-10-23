#include "SolutionWrapper.hpp"
#include "NameDouble.h"
#include <set>
#include <vector>

SolutionWrapper::SolutionWrapper(
    cxxSolution *soln, const std::vector<std::string> &_solution_order)
    : solution(soln), solution_order(_solution_order.begin() + NUM_ESSENTIALS,
                                     _solution_order.end()) {
  this->num_elements = _solution_order.size();

  auto &totals = solution->Get_totals();
}

void SolutionWrapper::get(std::span<LDBLE> &data) const {
  data[0] = solution->Get_total_h();
  data[1] = solution->Get_total_o();
  data[2] = solution->Get_cb();
  data[3] = solution->Get_total("H(0)");
  data[4] = solution->Get_total("O(0)");

  std::size_t i = NUM_ESSENTIALS;
  for (const auto &tot_name : solution_order) {
    auto it = solution->Get_totals().find(tot_name);
    if (it == solution->Get_totals().end()) {
      data[i++] = 0.0;
      continue;
    }
    data[i++] = it->second > 1e-25 ? it->second : 0.;
  }
}

void SolutionWrapper::set(const std::span<LDBLE> &data) {
  std::size_t i = NUM_ESSENTIALS;
  cxxNameDouble new_totals;

  const double &total_h = data[0];
  const double &total_o = data[1];
  const double &cb = data[2];

  new_totals["H(0)"] = data[3];
  new_totals["O(0)"] = data[4];

  for (const auto &tot_name : solution_order) {
    const double value = data[i++];

    if (value < 1E-25) {
      continue;
    }
    new_totals[tot_name] = value;
  }

  this->solution->Update(total_h, total_o, cb, new_totals);
}

std::vector<std::string>
SolutionWrapper::names(cxxSolution *solution,
                       std::vector<std::string> &solution_order) {
  std::vector<std::string> names;

  names.insert(names.end(), ESSENTIALS.begin(), ESSENTIALS.end());

  std::set<std::string> names_set;
  for (const auto &name : solution->Get_totals()) {
    names_set.insert(name.first);
  }

  for (const auto &to_erase : ESSENTIALS) {
    // don't care if the element was not found
    names_set.erase(to_erase);
  }

  names.insert(names.end(), names_set.begin(), names_set.end());
  return names;
}