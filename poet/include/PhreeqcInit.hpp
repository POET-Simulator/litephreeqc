#pragma once

#include <IPhreeqc.hpp>

#include <Phreeqc.h>
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "Exchange.h"
#include "PPassemblage.h"
#include "Solution.h"
#include "Surface.h"
#include "cxxKinetics.h"

enum { POET_SOL = 0, POET_EXCH, POET_KIN, POET_EQUIL, POET_SURF };

class PhreeqcInit : public IPhreeqc {
public:
  PhreeqcInit(const std::string &database, const std::string &input_script);

  struct PhreeqcMat {
    std::vector<std::string> names;
    std::vector<int> ids;
    std::vector<std::vector<double>> values;
  };

  PhreeqcMat getPhreeqcMat() const { return pqc_mat; };

  std::map<int, std::string> raw_dumps();

  using essential_names = std::array<std::vector<std::string>, 5>;
  using ModulesArray = std::array<std::uint32_t, 5>;

  ModulesArray getModuleSizes(const std::vector<int> &cell_ids);

  std::vector<std::string> getSolutionPrimaries() {
    return std::vector<std::string>(this->surface_primaries.begin(),
                                    this->surface_primaries.end());
  }

  std::vector<std::string> getSolutionNames(int cell_id) {
    return this->raw_initials[cell_id].first[POET_SOL];
  }

  std::vector<std::string> getExchanger(int id) {
    if (this->exchanger.contains(id)) {
      return this->exchanger[id];
    }

    return {};
  }

  std::vector<std::string> getKineticsNames(int id) {
    if (this->kinetics.contains(id)) {
      return this->kinetics[id];
    }

    return {};
  }

  std::vector<std::string> getEquilibriumNames(int id) {
    if (this->equilibrium.contains(id)) {
      return this->equilibrium[id];
    }

    return {};
  }

  std::vector<std::string> getSurfaceCompNames(int id) {
    if (this->surface_comps.contains(id)) {
      return this->surface_comps[id];
    }

    return {};
  }

  std::vector<std::string> getSurfaceChargeNames(int id) {
    if (this->surface_charge.contains(id)) {
      return this->surface_charge[id];
    }

    return {};
  }

private:
  // required only for simulation
  essential_names initial_names;

  void parseInitValues();
  void parseInit();

  std::vector<std::string> dump_solution_names(int cell_number);
  void dump_reactant_names(int cell_number,
                           const std::vector<std::string> union_sol_names,
                           essential_names &names);

  std::vector<std::string>
  find_all_valence_states(const std::vector<std::string> &solution_names,
                          const std::size_t offset);

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

  PhreeqcMat pqc_mat;

  PhreeqcMat buildPhreeqcMat();

  std::map<int, std::vector<std::string>> exchanger;
  std::map<int, std::vector<std::string>> kinetics;
  std::map<int, std::vector<std::string>> equilibrium;
  std::map<int, std::vector<std::string>> surface_comps;
  std::map<int, std::vector<std::string>> surface_charge;

  std::set<std::string> surface_primaries;

  using RawMap = std::map<int, std::pair<essential_names, std::vector<double>>>;

  essential_names union_raws(const RawMap &raws);

  std::vector<std::vector<double>>
  conc_from_essentials(const RawMap &raws, const essential_names &names);

  // void valuesFromModule(const std::string &module_name, int cell_number,
  //                       essential_names &names, std::vector<double> &values);

  std::string subExchangeName(std::string name);

  std::vector<double>
  get_essential_values_init(std::size_t cell_number,
                            const std::vector<std::string> &order);

  RawMap raw_initials;
};