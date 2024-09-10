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

/**
 * @brief Class for initializing Phreeqc and create input for POET
 *
 * This class is used to initialize Phreeqc by reading a database and an input
 * Phreeqc script. It then extracts the 'plain' initial values for each cell and
 * module, which are then used for further processing in POET.
 */
class PhreeqcInit : public IPhreeqc {
public:
  /**
   * @brief Construct a new Phreeqc Init object
   *
   * @param database String containing the content of the database file
   * @param input_script String containing the content of the input script
   */
  PhreeqcInit(const std::string &database, const std::string &input_script);

  /**
   * Hold column names, indexes and values for each cell defined in the input
   * script
   */
  struct PhreeqcMat {
    std::vector<std::string> names;
    std::vector<int> ids;
    std::vector<std::vector<double>> values;
  };

  /**
   * @brief Get the PhreeqcMatrix
   *
   * @return PhreeqcMat Value matrix of the parsed Phreeqc input script
   */
  PhreeqcMat getPhreeqcMat() const { return pqc_mat; };

  /**
   * @brief Get the *_RAW input scripts for each cell
   *
   * @return std::map<int, std::string> Contains the raw input scripts for each
   * cell
   */
  std::map<int, std::string> raw_dumps();

  using essential_names = std::array<std::vector<std::string>, 5>;
  using ModulesArray = std::array<std::uint32_t, 5>;

  /**
   * @brief Get the actual reactant sizes for given cells
   *
   * @param cell_ids Vector of cell ids to return reactant sizes for
   * @return ModulesArray Combined sizes of all cells for each module
   */
  ModulesArray getModuleSizes(const std::vector<int> &cell_ids);

  /**
   * @brief Get solution primaries found by PhreeqcInit
   *
   * @return std::vector<std::string> Vector containing only the solution
   * primary names of the solution defined by the input script
   */
  std::vector<std::string> getSolutionPrimaries() {
    return std::vector<std::string>(this->surface_primaries.begin(),
                                    this->surface_primaries.end());
  }

  /**
   * @brief Get the solution names for a given cell
   *
   * @param cell_id ID of the cell to get the solution names for
   * @return std::vector<std::string> Whole vector of solution names for the
   * cell
   */
  std::vector<std::string> getSolutionNames(int cell_id) {
    return this->raw_initials[cell_id].first[POET_SOL];
  }

  /**
   * @brief Get the exchange names for a given cell
   *
   * @param cell_id ID of the cell to get the exchange names for
   * @return std::vector<std::string> Whole vector of exchange names for the
   * cell. Empty if no exchange is defined.
   */
  std::vector<std::string> getExchanger(int id) {
    if (contains_id(this->exchanger, id)) {
      return this->exchanger[id];
    }

    return {};
  }

  /**
   * @brief Get the kinetics names for a given cell
   *
   * @param cell_id ID of the cell to get the kinetics names for
   * @return std::vector<std::string> Whole vector of kinetics names for the
   * cell. Empty if no kinetics are defined.
   */
  std::vector<std::string> getKineticsNames(int id) {
    if (contains_id(this->kinetics, id)) {
      return this->kinetics[id];
    }

    return {};
  }

  /**
   * @brief Get the equilibrium names for a given cell
   *
   * @param cell_id ID of the cell to get the equilibrium names for
   * @return std::vector<std::string> Whole vector of equilibrium names for the
   * cell. Empty if no equilibrium is defined.
   */
  std::vector<std::string> getEquilibriumNames(int id) {
    if (contains_id(this->equilibrium, id)) {
      return this->equilibrium[id];
    }

    return {};
  }

  /**
   * @brief Get the surface component names for a given cell
   *
   * @param cell_id ID of the cell to get the surface component names for
   * @return std::vector<std::string> Whole vector of surface component names
   * for the cell. Empty if no surface is defined.
   */
  std::vector<std::string> getSurfaceCompNames(int id) {
    if (contains_id(this->surface_comps, id)) {
      return this->surface_comps[id];
    }

    return {};
  }

  /**
   * @brief Get the surface charge names for a given cell
   *
   * @param cell_id ID of the cell to get the surface charge names for
   * @return std::vector<std::string> Whole vector of surface charge names for
   * the cell. Empty if no surface is defined.
   */
  std::vector<std::string> getSurfaceChargeNames(int id) {
    if (contains_id(this->surface_charge, id)) {
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

  bool contains_id(const std::map<int, std::vector<std::string>> &map,
                   int id) const {
#if __cplusplus >= 202002L
    return map.contains(id);
#else
    return map.find(id) != map.end();
#endif
  }

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