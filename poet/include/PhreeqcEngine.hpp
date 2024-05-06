
#pragma once

#include "EquilibriumWrapper.hpp"
#include "ExchangeWrapper.hpp"
#include "KineticWrapper.hpp"
#include "POETInit.hpp"
#include "SolutionWrapper.hpp"
#include "SurfaceWrapper.hpp"

#include <IPhreeqc.hpp>

#include <Exchange.h>
#include <PPassemblage.h>
#include <Phreeqc.h>
#include <Solution.h>
#include <Surface.h>
#include <cstddef>
#include <cxxKinetics.h>
#include <memory>
#include <span>
#include <string>
#include <vector>

/**
 * @brief Class for running Phreeqc wrappped in POET
 *
 * Direct interface to Phreeqc, without utilizing *_MODIFY keywords/scripts to
 * set new values. Use with already initialized Phreeqc config.
 *
 */
class PhreeqcEngine : public IPhreeqc {
public:
  /**
   * @brief Construct a new Phreeqc Engine object
   *
   * Construct a new Phreeqc Engine object by previously initialized POETConfig.
   *
   * @param config Holds the configuration for the Phreeqc engine.
   */
  PhreeqcEngine(const POETConfig &config) : IPhreeqc() {
    this->LoadDatabaseString(config.database.c_str());
    this->RunString(config.input_script.c_str());

    this->init_wrappers(config.cell);
  }

  /**
   * @brief Siimulate a cell for a given time step
   *
   * @param cell_values Vector containing the input values for the cell
   * (**including the ID**). Output values are written back in place to this
   * vector!
   * @param time_step Time step to simulate in seconds
   */
  void runCell(std::vector<double> &cell_values, double time_step);

private:
  void init_wrappers(const POETInitCell &cell);
  void run(double time_step);

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

  void get_essential_values(std::span<double> &data);

  void set_essential_values(const std::span<double> &data);

  std::unique_ptr<SolutionWrapper> solutionWrapperPtr;
  std::unique_ptr<ExchangeWrapper> exchangeWrapperPtr;
  std::unique_ptr<KineticWrapper> kineticsWrapperPtr;
  std::unique_ptr<EquilibriumWrapper> equilibriumWrapperPtr;
  std::unique_ptr<SurfaceWrapper> surfaceWrapperPtr;

  bool has_exchange = false;
  bool has_kinetics = false;
  bool has_equilibrium = false;
  bool has_surface = false;
};