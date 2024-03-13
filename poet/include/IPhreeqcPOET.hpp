#pragma once

#include <IPhreeqc.hpp>

#include <Exchange.h>
#include <PPassemblage.h>
#include <Phreeqc.h>
#include <Solution.h>
#include <Surface.h>
#include <array>
#include <cstddef>
#include <cxxKinetics.h>
#include <map>
#include <string>
#include <vector>

class IPhreeqcPOET : public IPhreeqc {
public:
  std::vector<double>
  get_essential_values(std::size_t cell_number,
                       const std::vector<std::string> &order);

  void set_essential_values(std::size_t cell_number,
                            const std::vector<std::string> &order,
                            std::vector<double> &values);

  IPhreeqcPOET(const std::string &database, const std::string &input_script)
      : IPhreeqc() {
    this->LoadDatabase(database.c_str());
    this->RunFile(input_script.c_str());

    this->parseInitValues();
  }

  std::vector<std::string> getInitNames() const {
    std::vector<std::string> names;

    for (auto &sol : this->initial_names) {
      names.insert(names.end(), sol.begin(), sol.end());
    }

    return names;
  }

  std::vector<std::vector<double>> getInitValues() const {
    return this->initial_concentrations;
  }

  std::vector<int> getSolutionIds() const { return this->solution_ids; }

  std::map<int, std::string> raw_dumps() {
    std::map<int, std::string> dumps;

    this->SetDumpStringOn(true);

    for (const auto &sol_id : this->solution_ids) {
      std::string call_string = "DUMP\n -cells " + std::to_string(sol_id);
      this->RunString(call_string.c_str());
      dumps[sol_id] = this->GetDumpString();
    }

    return dumps;
  }

private:
  using essential_names = std::array<std::vector<std::string>, 5>;

  essential_names initial_names;
  std::vector<std::vector<double>> initial_concentrations;
  std::vector<int> solution_ids;

  void parseInitValues();

  essential_names dump_essential_names(std::size_t cell_number);

  cxxSolution *Get_solution(std::size_t n) {
    return Utilities::Rxn_find(this->PhreeqcPtr->Get_Rxn_solution_map(), n);
  }

  cxxExchange *Get_exchange(std::size_t n) {
    return Utilities::Rxn_find(this->PhreeqcPtr->Get_Rxn_exchange_map(), n);
  }

  cxxKinetics *Get_kinetic(std::size_t n) {
    return Utilities::Rxn_find(this->PhreeqcPtr->Get_Rxn_kinetics_map(), n);
  }

  cxxPPassemblage *Get_equilibrium(std::size_t n) {
    return Utilities::Rxn_find(this->PhreeqcPtr->Get_Rxn_pp_assemblage_map(),
                               n);
  }

  cxxSurface *Get_surface(std::size_t n) {
    return Utilities::Rxn_find(this->PhreeqcPtr->Get_Rxn_surface_map(), n);
  }

  void resolveSolutionUseKW(
      const std::vector<IPhreeqc::SolutionMapping> &unresolved,
      std::map<int, std::pair<essential_names, std::vector<double>>>
          &mapped_values);

  void valuesFromModule(const std::string &module_name, int cell_number,
                        essential_names &names, std::vector<double> &values);

  std::string subExchangeName(std::string name) {
    for (const auto &species : this->PhreeqcPtr->Get_species_list()) {
      const std::string &species_name = species.s->name;
      // check if `name` is a prefix of `species_name`
      if (species_name.compare(0, name.size(), name) == 0) {
        return species_name;
      }
    }
    return name;
  };
};