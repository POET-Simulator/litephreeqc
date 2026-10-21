/*
 * Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
 *                       Marco De Lucia (GFZ German Research Centre for Geosciences)
 *
 * SPDX-License-Identifier: EUPL-1.2
 *
 * This file is part of litephreeqc, a C++ interface library on top of PHREEQC.
 */
#include "Phreeqc.h"
#include <set>

const std::set<std::string> to_ignore = {
    "H", "O", "Charge", "tc", "patm", "SolVol", "pH", "pe", "MassH2O", "Viscosity", "Density", "H(0)", "O(0)"};

std::vector<std::string> Phreeqc::find_all_valence_states(
    const std::vector<std::string> &solution_names) {
  std::vector<std::string> solution_with_valences;
  solution_with_valences.reserve(solution_names.size());

  // to avoid duplicates store already evaluated master species
  std::set<std::string> master_species_found;

  auto is_ignored = [](const std::string &in) {
    return (to_ignore.find(in) != to_ignore.end());
  };

  for (std::size_t i = 0; i < solution_names.size(); i++) {

    if (is_ignored(solution_names[i])) {
      solution_with_valences.push_back(solution_names[i]);
      continue;
    }

    const auto master_primary =
        master_bsearch_primary(solution_names[i].c_str());

    assert(master_primary != NULL);

    const std::string master_species(master_primary->elt->name);

    // in case we already found valences for this master species we skip it
    if (master_species_found.find(master_species) !=
        master_species_found.end()) {
      continue;
    }

    master_species_found.insert(master_species);

    std::size_t start_idx = 0;
    if (master_primary->number >= 0 &&
        master_primary->number < static_cast<int>(master.size()) &&
        master[master_primary->number] == master_primary) {
      start_idx = static_cast<std::size_t>(master_primary->number);
    } else {
      while (start_idx < master.size() && master[start_idx] != master_primary) {
        start_idx++;
      }
    }

    bool has_valences = false;
    std::size_t last_valence = start_idx;

    while (last_valence + 1 < master.size()) {
      const auto next_master = master[last_valence + 1];
      if (next_master->elt && next_master->elt->primary == master_primary) {
        has_valences = true;
        last_valence++;
      } else {
        break;
      }
    }

    // in case we found valences we add them to the solution
    if (has_valences) {
      // skip the master species and only add the valences
      for (std::size_t k = start_idx + 1; k <= last_valence; k++) {
        solution_with_valences.push_back(master[k]->elt->name);
      }
      continue;
    }

    // otherwise we just add the master species without any valences
    solution_with_valences.push_back(master_species);
  }

  return solution_with_valences;
}
