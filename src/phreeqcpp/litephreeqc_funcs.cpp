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

#include "Phreeqc.h"
#include <set>

const std::set<std::string> to_ignore = {
    "H", "O", "Charge", "tc", "patm", "SolVol", "pH", "pe", "H(0)", "O(0)"};

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

    bool has_valences = false;
    std::size_t inner_loop_j = 0;
    std::size_t last_valence = inner_loop_j;

    // we HAVE to assume master species are already sorted!
    // loop over all known master species
    for (inner_loop_j = 0; inner_loop_j < master.size(); inner_loop_j++) {
      std::string curr_master_species(master[inner_loop_j]->elt->name);

      if (curr_master_species == master_species) {
        last_valence = inner_loop_j;

        while (last_valence < master.size() - 1) {
          std::string next_master_species(master[last_valence + 1]->elt->name);

          // check if the next name starts with the current master species
          if (next_master_species.compare(0, master_species.size(),
                                          master_species) != 0) {
            break;
          }

          // check if the next character is an opening parenthesis
          if (next_master_species[master_species.size()] != '(') {
            break;
          }

          // if this is all true we have a valence state
          has_valences = true;
          last_valence++;
        }
        break;
      }
    }

    // in case we found valences we add them to the solution
    if (has_valences) {
      // skip the master species and only add the valences
      for (std::size_t k = inner_loop_j + 1; k <= last_valence; k++) {
        solution_with_valences.push_back(master[k]->elt->name);
      }
      continue;
    }

    // otherwise we just add the master species without any valences
    solution_with_valences.push_back(master_species);
  }

  return solution_with_valences;
}
