#include "Phreeqc.h"

std::vector<std::string> Phreeqc::find_all_valence_states(
    const std::vector<std::string> &&solution_names, const std::size_t offset) {
  std::vector<std::string> solution_with_valences(
      solution_names.begin(), solution_names.begin() + offset);

  // to avoid duplicates store already evaluated master species
  std::set<std::string> master_species_found;

  for (std::size_t i = offset; i < solution_names.size(); i++) {
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