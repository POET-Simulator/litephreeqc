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

std::size_t SolutionWrapper::countActiveOptionals() const {
  std::size_t count = 0;
  for (const auto &flag : OPTIONAL_ESSENTIAL_FLAGS) {
    if (PhreeqcMatrix::hasFlag(_optional_essentials, flag)) {
      count++;
    }
  }
  return count;
}

SolutionWrapper::SolutionWrapper(
    cxxSolution *soln, const std::vector<std::string> &_solution_order,
    bool with_redox, PhreeqcMatrix::OptionalEssentials optional_essentials)
    : solution(soln),
      solution_order(_solution_order.begin() +
                         static_cast<std::ptrdiff_t>(NUM_MANDATORY_ESSENTIALS) +
                         static_cast<std::ptrdiff_t>([optional_essentials]() {
                           std::size_t count = 0;
                           for (const auto &flag : OPTIONAL_ESSENTIAL_FLAGS) {
                             if (PhreeqcMatrix::hasFlag(optional_essentials, flag)) {
                               count++;
                             }
                           }
                           return count;
                         }()),
                     _solution_order.end()),
      _with_redox(with_redox), _optional_essentials(optional_essentials) {
  this->num_elements = _solution_order.size();
}

void SolutionWrapper::get(std::span<LDBLE> &data) const {
  using OE = PhreeqcMatrix::OptionalEssentials;
  std::size_t idx = 0;

  // Mandatory essentials (always indices 0-4)
  data[idx++] = solution->Get_total_h();
  data[idx++] = solution->Get_total_o();
  data[idx++] = solution->Get_cb();
  data[idx++] = solution->Get_tc();
  data[idx++] = solution->Get_patm();

  // Optional essentials (conditional)
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::SolVol))
    data[idx++] = solution->Get_soln_vol();
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::pH))
    data[idx++] = solution->Get_ph();
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::pe))
    data[idx++] = solution->Get_pe();
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::MassH2O))
    data[idx++] = solution->Get_mass_water();
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::Viscosity))
    data[idx++] = solution->Get_viscosity();
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::Density))
    data[idx++] = solution->Get_density();

  const cxxNameDouble &totals =
      (_with_redox ? solution->Get_totals()
                   : solution->Get_totals().Simplify_redox());

  // Totals start at dynamic idx
  for (const auto &tot_name : solution_order) {
    auto it = totals.find(tot_name);
    if (it == totals.end()) {
      data[idx++] = 0.0;
      continue;
    }
    data[idx++] = it->second > 1e-25 ? it->second : 0.;
  }
}

void SolutionWrapper::set(const std::span<LDBLE> &data) {
  using OE = PhreeqcMatrix::OptionalEssentials;
  std::size_t idx = 0;
  cxxNameDouble new_totals;

  // Mandatory essentials (always first)
  const double total_h = data[idx++];
  const double total_o = data[idx++];
  const double cb = data[idx++];
  const double tc = data[idx++];
  const double patm = data[idx++];

  // Optional essentials (conditional) - use defaults if not present
  double ph = 7.0;
  double pe = 4.0;
  double massh2o = total_o/55.5;
  double soln_vol = massh2o;
  double viscosity = 1.0;
  double density = 1.0;

  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::SolVol))
    soln_vol = data[idx++];
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::pH))
    ph = data[idx++];
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::pe))
    pe = data[idx++];
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::MassH2O))
    massh2o = data[idx++];
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::Viscosity))
    viscosity = data[idx++];
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::Density))
    density = data[idx++];

  for (const auto &tot_name : solution_order) {
    const double value = data[idx++];

    if (value < 1E-25) {
      continue;
    }
    new_totals[tot_name] = value;
  }

  this->solution->Update(total_h, total_o, cb, tc, patm,
                         _with_redox ? new_totals
                                     : new_totals.Simplify_redox());

  // Set optional properties directly on the solution object after Update(),
  // since Update() does not handle them. Without this, the solution retains
  // stale values from the previous cell processed by this Engine.
  // Only overwrite with positive values; if transport zeroed them out,
  // keep the fallback from Update() (e.g. mass_water = o_tot/55.55).
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::SolVol) && soln_vol > 0)
    this->solution->Set_soln_vol(soln_vol);
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::pH))
    this->solution->Set_ph(ph);
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::pe))
    this->solution->Set_pe(pe);
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::MassH2O) && massh2o > 0)
    this->solution->Set_mass_water(massh2o);
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::Viscosity) && viscosity > 0)
    this->solution->Set_viscosity(viscosity);
  if (PhreeqcMatrix::hasFlag(_optional_essentials, OE::Density) && density > 0)
    this->solution->Set_density(density);
}

std::vector<std::string>
SolutionWrapper::names(cxxSolution *solution, bool include_h0_o0,
                       std::vector<std::string> &solution_order, bool with_redox,
                       PhreeqcMatrix::OptionalEssentials optional_essentials) {
  std::vector<std::string> names;

  // Add mandatory essentials
  names.insert(names.end(), MANDATORY_ESSENTIALS.begin(),
               MANDATORY_ESSENTIALS.end());

  // Add only the enabled optional essentials
  for (std::size_t i = 0; i < OPTIONAL_ESSENTIAL_NAMES.size(); ++i) {
    if (PhreeqcMatrix::hasFlag(optional_essentials, OPTIONAL_ESSENTIAL_FLAGS[i])) {
      names.push_back(OPTIONAL_ESSENTIAL_NAMES[i]);
    }
  }

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
