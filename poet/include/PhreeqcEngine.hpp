
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
#include <iomanip>
#include <iostream>
#include <memory>
#include <ostream>
#include <span>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <vector>

class PhreeqcEngine : public IPhreeqc {
public:
  PhreeqcEngine(const POETConfig &config) : IPhreeqc() {
    this->LoadDatabaseString(config.database.c_str());
    this->RunString(config.input_script.c_str());

    this->init_wrappers(config.cell);
  }

  void runCell(std::vector<double> &cell_values, double time_step);

  void run(int start_cell, int end_cell) {
    const std::string runs_string = "RUN_CELLS\n -cells " +
                                    std::to_string(start_cell) + "-" +
                                    std::to_string(end_cell) + "\nEND\n";
    this->RunString(runs_string.c_str());
  }

  void run(int start_cell, int end_cell, double time_step) {
    std::stringstream time_ss;
    time_ss << std::fixed << std::setprecision(20) << time_step;

    const std::string runs_string =
        "RUN_CELLS\n -cells " + std::to_string(start_cell) + "-" +
        std::to_string(end_cell) + "\n -time_step " + time_ss.str() + "\nEND\n";
    this->RunString(runs_string.c_str());
  }

private:
  void init_wrappers(const POETInitCell &cell);

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