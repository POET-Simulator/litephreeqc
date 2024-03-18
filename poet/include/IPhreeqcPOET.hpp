#pragma once

#include <IPhreeqc.hpp>

#include <Exchange.h>
#include <PPassemblage.h>
#include <Phreeqc.h>
#include <Solution.h>
#include <Surface.h>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cxxKinetics.h>
#include <map>
#include <string>
#include <sys/types.h>
#include <utility>
#include <vector>

enum { POET_SOL = 0, POET_EXCH, POET_KIN, POET_EQUIL, POET_SURF };

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

    this->parseInit();
  }

  IPhreeqcPOET(const std::string &database, const std::string &input_script,
               const std::vector<std::string> &solutionInitVector,
               std::uint32_t n_cells)
      : IPhreeqc(), n_cells(n_cells), solutionInitVector(solutionInitVector) {
    this->LoadDatabase(database.c_str());
    this->RunFile(input_script.c_str());

    if (n_cells > 1) {
      std::string copy_string =
          "COPY cell 1 2-" + std::to_string(n_cells) + "\n";
    }
  }

  struct PhreeqcMat {
    std::vector<std::string> names;
    std::vector<int> ids;
    std::vector<std::vector<double>> values;
  };

  PhreeqcMat getPhreeqcMat();

  PhreeqcMat getPhreeqcMat(const std::vector<int> &ids);

  std::map<int, std::string> raw_dumps() {
    std::map<int, std::string> dumps;

    this->SetDumpStringOn(true);

    for (const auto &[sol_id, unused] : this->raw_initials) {
      std::string call_string = "DUMP\n -cells " + std::to_string(sol_id);
      this->RunString(call_string.c_str());
      dumps[sol_id] = this->GetDumpString();
    }

    this->SetDumpStringOn(false);

    return dumps;
  }

  using essential_names = std::array<std::vector<std::string>, 5>;
  using ModulesArray = std::array<std::uint32_t, 5>;

  ModulesArray getModuleSizes() const {
    ModulesArray module_sizes;
    for (std::uint8_t i = 0; i < 5; i++) {
      module_sizes[i] = this->initial_names[i].size();
    }

    return module_sizes;
  }

private:
  // required only for simulation
  essential_names initial_names;
  std::uint32_t n_cells;
  std::vector<std::string> solutionInitVector;

  void parseInitValues();
  void parseInit();

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

  using RawMap = std::map<int, std::pair<essential_names, std::vector<double>>>;

  essential_names union_raws(const RawMap &raws);

  std::vector<std::vector<double>>
  conc_from_essentials(const RawMap &raws, const essential_names &names);

  // void valuesFromModule(const std::string &module_name, int cell_number,
  //                       essential_names &names, std::vector<double> &values);

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

  RawMap raw_initials;
};