#include "EquilibriumWrapper.hpp"
#include "ExchangeWrapper.hpp"
#include "KineticWrapper.hpp"
#include "PhreeqcInit.hpp"
#include "SolutionWrapper.hpp"
#include "SurfaceWrapper.hpp"
#include "global_structures.h"

#include <cassert>
#include <cstddef>
#include <set>
#include <span>
#include <string>
#include <vector>

std::vector<std::string> PhreeqcInit::dump_solution_names(int cell_number) {
  constexpr std::size_t SKIP_H_O_CB = 3;

  std::vector<std::string> placeholder;

  return find_all_valence_states(
      SolutionWrapper::names(this->Get_solution(cell_number), placeholder),
      SKIP_H_O_CB);
}

void PhreeqcInit::dump_reactant_names(
    int cell_number, const std::vector<std::string> union_sol_names,
    essential_names &names) {
  // Exchange
  if (this->Get_exchange(cell_number) != NULL) {
    auto *exchange = this->Get_exchange(cell_number);

    std::vector<std::string> exch_formulas;
    names[POET_EXCH] = ExchangeWrapper::names(exchange, exch_formulas);

    this->exchanger[cell_number] = exch_formulas;
  }

  // Kinetics
  if (this->Get_kinetic(cell_number) != NULL) {
    auto *kinetic = this->Get_kinetic(cell_number);

    std::vector<std::string> kin_comps;

    names[POET_KIN] = KineticWrapper::names(kinetic, kin_comps);

    this->kinetics[cell_number] = kin_comps;
  }

  // PPassemblage
  if (this->Get_equilibrium(cell_number) != NULL) {
    auto *equilibrium = this->Get_equilibrium(cell_number);

    std::vector<std::string> equ_names;

    names[POET_EQUIL] = EquilibriumWrapper::names(equilibrium, equ_names);

    this->equilibrium[cell_number] = equ_names;
  }

  // Surface
  if (this->Get_surface(cell_number) != NULL) {
    auto *surface = this->Get_surface(cell_number);

    if (this->surface_primaries.empty()) {
      // this is fixed! Always add H and O
      this->surface_primaries.insert("H");
      this->surface_primaries.insert("O");

      for (std::size_t i = 3; i < union_sol_names.size(); i++) {
        const auto master_primary = this->PhreeqcPtr->master_bsearch_primary(
            union_sol_names[i].c_str());
        if (master_primary != NULL) {
          this->surface_primaries.insert(master_primary->elt->name);
        }
      }
    }

    std::vector<std::string> comp_formulas;
    std::vector<std::string> charge_names;

    names[POET_SURF] = SurfaceWrapper::names(surface, this->surface_primaries,
                                             comp_formulas, charge_names);

    this->surface_comps[cell_number] = comp_formulas;
    this->surface_charge[cell_number] = charge_names;
  }
}

std::vector<std::string> PhreeqcInit::find_all_valence_states(
    const std::vector<std::string> &solution_names, const std::size_t offset) {
  std::vector<std::string> solution_with_valences(
      solution_names.begin(), solution_names.begin() + offset);

  // to avoid duplicates store already evaluated master species
  std::set<std::string> master_species_found;

  for (std::size_t i = offset; i < solution_names.size(); i++) {
    const auto master_primary =
        this->PhreeqcPtr->master_bsearch_primary(solution_names[i].c_str());

    assert(master_primary != NULL);

    const std::string master_species(master_primary->elt->name);

    // in case we already found valences for this master species we skip it
    if (master_species_found.contains(master_species)) {
      continue;
    }

    master_species_found.insert(master_species);

    bool has_valences = false;
    std::size_t inner_loop_j = 0;
    std::size_t last_valence = inner_loop_j;

    // we HAVE to assume master species are already sorted!
    // loop over all known master species
    for (inner_loop_j = 0; inner_loop_j < this->PhreeqcPtr->master.size();
         inner_loop_j++) {
      std::string curr_master_species(
          this->PhreeqcPtr->master[inner_loop_j]->elt->name);

      if (curr_master_species == master_species) {
        last_valence = inner_loop_j;

        while (last_valence < this->PhreeqcPtr->master.size() - 1) {
          std::string next_master_species(
              this->PhreeqcPtr->master[last_valence + 1]->elt->name);

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
        solution_with_valences.push_back(
            this->PhreeqcPtr->master[k]->elt->name);
      }
      continue;
    }

    // otherwise we just add the master species without any valences
    solution_with_valences.push_back(master_species);
  }

  return solution_with_valences;
}

std::vector<double>
PhreeqcInit::get_essential_values_init(std::size_t cell_number,
                                       const std::vector<std::string> &order) {
  std::vector<double> essentials;

  // Solutions
  if (this->Get_solution(cell_number) != NULL) {
    SolutionWrapper solution(this->Get_solution(cell_number), order);
    std::vector<double> sol_values(solution.size());
    std::span<double> sol_span{sol_values};

    solution.get(sol_span);

    essentials.insert(essentials.end(), sol_values.begin(), sol_values.end());
  }

  // Exchange
  if (this->Get_exchange(cell_number) != NULL) {
    ExchangeWrapper exchange(this->Get_exchange(cell_number),
                             this->exchanger[cell_number]);

    std::vector<double> exch_values(exchange.size());
    std::span<double> exch_span{exch_values};

    exchange.get(exch_span);

    essentials.insert(essentials.end(), exch_values.begin(), exch_values.end());
  }

  // Kinetics
  if (this->Get_kinetic(cell_number) != NULL) {
    KineticWrapper kinetic(this->Get_kinetic(cell_number),
                           this->kinetics[cell_number]);

    std::vector<double> kin_values(kinetic.size());
    std::span<double> kin_span{kin_values};

    kinetic.get(kin_span);

    essentials.insert(essentials.end(), kin_values.begin(), kin_values.end());
  }

  // PPassemblage
  if (this->Get_equilibrium(cell_number) != NULL) {
    EquilibriumWrapper equilibrium(this->Get_equilibrium(cell_number),
                                   this->equilibrium[cell_number]);

    std::vector<double> equ_values(equilibrium.size());
    std::span<double> equ_span{equ_values};

    equilibrium.get(equ_span);

    essentials.insert(essentials.end(), equ_values.begin(), equ_values.end());
  }

  // Surface
  if (this->Get_surface(cell_number) != NULL) {

    SurfaceWrapper surface(
        this->Get_surface(cell_number), this->surface_primaries,
        this->surface_comps[cell_number], this->surface_charge[cell_number]);

    std::vector<double> surf_values(surface.size());
    std::span<double> surf_span{surf_values};

    surface.get(surf_span);

    essentials.insert(essentials.end(), surf_values.begin(), surf_values.end());
  }

  return essentials;
}

std::map<int, std::string> PhreeqcInit::raw_dumps() {
  std::map<int, std::string> dumps;

  this->SetDumpStringOn(true);

  for (const auto &[sol_id, unused] : this->raw_initials) {
    std::string call_string =
        "DUMP\n -cells " + std::to_string(sol_id) + "\nEND\n";
    this->RunString(call_string.c_str());
    dumps[sol_id] = this->GetDumpString();
  }

  this->SetDumpStringOn(false);

  return dumps;
}