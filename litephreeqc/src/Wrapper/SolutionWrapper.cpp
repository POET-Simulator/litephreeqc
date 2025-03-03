/*
 * This project is subject to the original PHREEQC license. `litephreeqc` is a
 * version of the PHREEQC code that has been modified to be used as a library.
 *
 * It adds a C++ interface on top of the original PHREEQC code, with small
 * changes to the original code base.
 *
 * Authors of Modifications:
 * - Max Luebke (mluebke@uni-potsdam.de) - University of Potsdam
 * - Marco De Lucia (delucia@gfz.de) - GFZ Helmholz Centre for Geosciences
 *
 */

#include "SolutionWrapper.hpp"
#include "NameDouble.h"
#include <set>
#include <vector>

SolutionWrapper::SolutionWrapper(
    cxxSolution *soln, const std::vector<std::string> &_solution_order,
    bool with_redox)
    : solution(soln), solution_order(_solution_order.begin() + NUM_ESSENTIALS,
                                     _solution_order.end()),
      _with_redox(with_redox) {
  this->num_elements = _solution_order.size();

  auto &totals = solution->Get_totals();
}

void SolutionWrapper::get(std::span<LDBLE> &data) const {
  data[0] = solution->Get_total_h();
  data[1] = solution->Get_total_o();
  data[2] = solution->Get_cb();
  data[3] = solution->Get_tc();
  data[4] = solution->Get_patm();
  data[5] = solution->Get_soln_vol();
  data[6] = solution->Get_ph();
  data[7] = solution->Get_pe();

  const cxxNameDouble &totals =
      (_with_redox ? solution->Get_totals()
                   : solution->Get_totals().Simplify_redox());

  std::size_t i = NUM_ESSENTIALS;
  for (const auto &tot_name : solution_order) {
    auto it = totals.find(tot_name);
    if (it == totals.end()) {
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
  const double &tc = data[3];
  const double &patm = data[4];

  for (const auto &tot_name : solution_order) {
    const double value = data[i++];

    if (value < 1E-25) {
      continue;
    }
    new_totals[tot_name] = value;
  }

  this->solution->Update(total_h, total_o, cb, tc, patm,
                         _with_redox ? new_totals
                                     : new_totals.Simplify_redox());
}

std::vector<std::string>
SolutionWrapper::names(cxxSolution *solution, bool include_h0_o0,
                       std::vector<std::string> &solution_order,
                       bool with_redox) {
  std::vector<std::string> names;

  names.insert(names.end(), ESSENTIALS.begin(), ESSENTIALS.end());

  if (include_h0_o0) {
    names.push_back("H(0)");
    names.push_back("O(0)");
  }

  std::set<std::string> names_set;

  const cxxNameDouble &totals =
      (with_redox ? solution->Get_totals()
                  : solution->Get_totals().Simplify_redox());

  for (const auto &[name, _] : totals) {
    // Skip redox states of H and O
    if (name == "H(0)" || name == "O(0)") {
      continue;
    }
    names_set.insert(name);
  }

  names.insert(names.end(), names_set.begin(), names_set.end());
  return names;
}
