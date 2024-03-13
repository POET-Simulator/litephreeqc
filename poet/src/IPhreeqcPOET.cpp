#include "IPhreeqcPOET.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <string>
#include <vector>

static inline std::vector<std::string>
unionStringVectors(const std::vector<std::string> &a,
                   const std::vector<std::string> &b) {
  std::vector<std::string> result;

  std::set_union(a.begin(), a.end(), b.begin(), b.end(),
                 std::back_inserter(result));

  return result;
}

static std::vector<double>
createConcVector(const std::vector<std::string> &conc_names,
                 const std::vector<double>::const_iterator &conc_it,
                 const std::vector<std::string> &new_names) {
  std::vector<double> conc_vec;

  for (const auto &i : new_names) {
    auto it = std::find(conc_names.begin(), conc_names.end(), i);

    auto conc_index = std::distance(conc_names.begin(), it);

    if (it != conc_names.end()) {
      conc_vec.push_back(*(conc_it + conc_index));
    } else {
      conc_vec.push_back(std::numeric_limits<double>::quiet_NaN());
    }
  }

  return conc_vec;
}

void IPhreeqcPOET::valuesFromModule(const std::string &module_name,
                                    int cell_number, essential_names &names,
                                    std::vector<double> &values) {
  std::size_t dest_module_i = 0;
  std::vector<double> to_insert;
  if (module_name == "exchange") { // 1
    this->Get_exchange(cell_number)->dump_essential_names(names[POET_EXCH]);
    this->Get_exchange(cell_number)->get_essential_values(to_insert);
    dest_module_i = 1;
  } else if (module_name == "kinetics") { // 2
    this->Get_kinetic(cell_number)->dump_essential_names(names[POET_KIN]);
    this->Get_kinetic(cell_number)->get_essential_values(to_insert);
    dest_module_i = 2;
  } else if (module_name == "surface") { // 4
    this->Get_surface(cell_number)->dump_essential_names(names[POET_SURF]);
    this->Get_surface(cell_number)->get_essential_values(to_insert);
    dest_module_i = 4;
  }

  std::size_t offset = 0;
  for (std::size_t i = 0; i < dest_module_i; i++) {
    offset += names[i].size();
  }

  values.insert(values.begin() + offset, to_insert.begin(), to_insert.end());
}

void IPhreeqcPOET::resolveSolutionUseKW(
    const std::vector<IPhreeqc::SolutionMapping> &unresolved,
    std::map<int, std::pair<essential_names, std::vector<double>>>
        &mapped_values) {

  for (const auto &input : unresolved) {
    if (mapped_values.find(input.module_n) != mapped_values.end()) {
      continue;
    }

    essential_names new_conc_names;
    new_conc_names[0] = mapped_values[input.sol_n].first[0];

    const auto &curr_sol_vec = mapped_values[input.sol_n].second;
    std::vector<double> new_conc_values(
        curr_sol_vec.begin(), curr_sol_vec.begin() + new_conc_names[0].size());

    valuesFromModule(input.module_name, input.module_n, new_conc_names,
                     new_conc_values);

    mapped_values[input.module_n] =
        std::make_pair(new_conc_names, new_conc_values);
  }
}

void IPhreeqcPOET::parseInitValues() {

  std::map<int, std::pair<essential_names, std::vector<double>>> init_values;

  for (const auto &[id, val] : this->PhreeqcPtr->Get_Rxn_solution_map()) {
    // A key less than zero indicates an internal solution
    if (id < 0) {
      continue;
    }

    essential_names curr_conc_names = this->dump_essential_names(id);

    init_values[id] = std::make_pair(
        curr_conc_names, this->get_essential_values(id, curr_conc_names[0]));
  }

  const auto unresolved_modules = this->getSolutionMapping();
  resolveSolutionUseKW(unresolved_modules, init_values);

  std::vector<int> ids_to_erase;

  // loop over found initial values and erase those that are not in the ids
  for (const auto &[id, values] : init_values) {
    // find key in vector of ids
    // auto it = std::find(ids.begin(), ids.end(), id);

    // if (it == ids.end()) {
    //   ids_to_erase.push_back(id);
    //   continue;
    // }

    // create a union of all known concentration names
    for (std::size_t i = 0; i < 5; i++) {
      this->initial_names[i] =
          unionStringVectors(this->initial_names[i], values.first[i]);

      if (i == 1) {
        for (auto &specie : this->initial_names[i]) {
          specie = subExchangeName(specie);
        }
      }
    }
  }

  for (const auto &key : ids_to_erase) {
    init_values.erase(key);
  }

  // create a vector of the initial concentrations with NaNs for missing
  // concentrations
  for (const auto &[key, val] : init_values) {
    std::vector<double> combined_conc_vec;

    for (std::size_t i = 0, offset = 0; i < 5; i++) {
      std::vector<double> union_vec = createConcVector(
          val.first[i], val.second.begin() + offset, this->initial_names[i]);

      combined_conc_vec.insert(combined_conc_vec.end(), union_vec.begin(),
                               union_vec.end());

      offset += val.first[i].size();
    }

    this->initial_concentrations.push_back(combined_conc_vec);
    this->solution_ids.push_back(key);
  }

  // const auto unresolved_modules = this->getSolutionMapping();

  // std::map<int, std::pair<essential_names, std::vector<double>>>
  // mapped_values;

  // for (const auto &val : unresolved_modules) {
  //   const auto &module_n = val.module_n;

  //   if (init_values.find(module_n) != init_values.end()) {
  //     break;
  //   }

  //   if (std::find(this->solution_ids.begin(), this->solution_ids.end(),
  //                 module_n) == this->solution_ids.end()) {
  //     const int solution_index =
  //         std::find(solution_ids.begin(), solution_ids.end(), val.sol_n) -
  //         solution_ids.begin();

  //     essential_names solution_only;
  //     solution_only[0] = init_values[val.sol_n].first[0];

  //     mapped_values[module_n] = std::make_pair(
  //         solution_only, init_values[val.sol_n].second);
  //   }
  // }
}

IPhreeqcPOET::essential_names
IPhreeqcPOET::dump_essential_names(std::size_t cell_number) {

  essential_names eNames;

  // Solutions
  if (this->Get_solution(cell_number) != NULL) {
    std::vector<std::string> &eSolNames = eNames[POET_SOL];
    this->Get_solution(cell_number)->dump_essential_names(eSolNames);
  }

  // Exchange
  if (this->Get_exchange(cell_number) != NULL) {
    std::vector<std::string> &eExchNames = eNames[POET_EXCH];
    this->Get_exchange(cell_number)->dump_essential_names(eExchNames);
  }

  // Kinetics
  if (this->Get_kinetic(cell_number) != NULL) {
    std::vector<std::string> &eKinNames = eNames[POET_KIN];
    this->Get_kinetic(cell_number)->dump_essential_names(eKinNames);
  }

  // PPassemblage
  if (this->Get_equilibrium(cell_number) != NULL) {
    std::vector<std::string> &eEquNames = eNames[POET_EQUIL];
    this->Get_equilibrium(cell_number)->dump_essential_names(eEquNames);
  }

  // Surface
  if (this->Get_surface(cell_number) != NULL) {
    std::vector<std::string> &eSurfNames = eNames[POET_SURF];
    this->Get_surface(cell_number)->dump_essential_names(eSurfNames);
  }

  return eNames;
}

std::vector<double>
IPhreeqcPOET::get_essential_values(std::size_t cell_number,
                                   const std::vector<std::string> &order) {
  std::vector<double> essentials;

  // Solutions
  if (this->Get_solution(cell_number) != NULL) {
    std::vector<double> sol_values;

    this->Get_solution(cell_number)->get_essential_values(sol_values, order);
    essentials.insert(essentials.end(), sol_values.begin(), sol_values.end());
  }

  // Exchange
  if (this->Get_exchange(cell_number) != NULL) {
    std::vector<double> exch_values;

    this->Get_exchange(cell_number)->get_essential_values(exch_values);
    essentials.insert(essentials.end(), exch_values.begin(), exch_values.end());
  }

  // Kinetics
  if (this->Get_kinetic(cell_number) != NULL) {
    std::vector<double> kin_values;

    this->Get_kinetic(cell_number)->get_essential_values(kin_values);
    essentials.insert(essentials.end(), kin_values.begin(), kin_values.end());
  }

  // PPassemblage
  if (this->Get_equilibrium(cell_number) != NULL) {
    std::vector<double> equ_values;

    this->Get_equilibrium(cell_number)->get_essential_values(equ_values);
    essentials.insert(essentials.end(), equ_values.begin(), equ_values.end());
  }

  // Surface
  if (this->Get_surface(cell_number) != NULL) {
    std::vector<double> surf_values;

    this->Get_surface(cell_number)->get_essential_values(surf_values);
    essentials.insert(essentials.end(), surf_values.begin(), surf_values.end());
  }

  return essentials;
}

void IPhreeqcPOET::set_essential_values(std::size_t cell_number,
                                        const std::vector<std::string> &order,
                                        std::vector<double> &values) {

  auto dump_it = values.begin();

  // Solutions
  if (this->Get_solution(cell_number) != NULL) {
    this->Get_solution(cell_number)->set_essential_values(dump_it, order);
  }

  // Exchange
  if (this->Get_exchange(cell_number) != NULL) {
    this->Get_exchange(cell_number)->set_essential_values(dump_it);
  }

  // Kinetics
  if (this->Get_kinetic(cell_number) != NULL) {
    this->Get_kinetic(cell_number)->set_essential_values(dump_it);
  }

  // PPassemblage
  if (this->Get_equilibrium(cell_number) != NULL) {
    this->Get_equilibrium(cell_number)->set_essential_values(dump_it);
  }

  // Surface
  if (this->Get_surface(cell_number) != NULL) {
    this->Get_surface(cell_number)->set_essential_values(dump_it);
  }
}