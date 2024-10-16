#include "SolutionWrapper.hpp"
#include "NameDouble.h"
#include <vector>

SolutionWrapper::SolutionWrapper(
    cxxSolution *soln, const std::vector<std::string> &solution_order_)
    : solution(soln), solution_order(solution_order_) {
  this->num_elements = solution_order.size() + NUM_ESSENTIALS;

  auto &totals = solution->Get_totals();
}

void SolutionWrapper::get(std::span<LDBLE> &data) const {
  data[0] = solution->Get_total_h();
  data[1] = solution->Get_total_o();
  data[2] = solution->Get_cb();

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
  for (const auto &tot_name : solution_order) {
    const double value = data[i++];

    if (value < 1E-25) {
      continue;
    }
    new_totals[tot_name] = value;
  }

  this->solution->Update(data[0], data[1], data[2], new_totals);
}

std::vector<std::string>
SolutionWrapper::names(cxxSolution *solution,
                       std::vector<std::string> &solution_order) {
  std::vector<std::string> names;

  names.push_back("H");
  names.push_back("O");
  names.push_back("Charge");

  for (const auto &[tot_name, value] : solution->Get_totals()) {
    if (tot_name == "H(0)" || tot_name == "O(0)") {
      continue;
    }

    solution_order.push_back(tot_name);
  }

  names.insert(names.end(), solution_order.begin(), solution_order.end());
  return names;
}