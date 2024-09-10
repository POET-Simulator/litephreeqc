#include "PhreeqcEngine.hpp"
#include <cstddef>
#include <iomanip>
#include <span>
#include <sstream>

#include <IPhreeqc.hpp>
#include <Phreeqc.h>

#include "Wrapper/EquilibriumWrapper.hpp"
#include "Wrapper/ExchangeWrapper.hpp"
#include "Wrapper/KineticWrapper.hpp"
#include "Wrapper/SolutionWrapper.hpp"
#include "Wrapper/SurfaceWrapper.hpp"
class PhreeqcEngine::Impl : public IPhreeqc {
public:
  Impl(const POETConfig &config);
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

PhreeqcEngine::PhreeqcEngine(const POETConfig &config)
    : impl(std::make_unique<Impl>(config)) {}

PhreeqcEngine::Impl::Impl(const POETConfig &config) {
  this->LoadDatabaseString(config.database.c_str());
  this->RunString(config.input_script.c_str());

  this->init_wrappers(config.cell);
}

PhreeqcEngine::~PhreeqcEngine() = default;

void PhreeqcEngine::runCell(std::vector<double> &cell_values,
                            double time_step) {
  // skip ID
  std::span<double> cell_data{cell_values.begin() + 1, cell_values.end()};

  this->impl->set_essential_values(cell_data);
  this->impl->run(time_step);
  this->impl->get_essential_values(cell_data);
}

void PhreeqcEngine::Impl::run(double time_step) {
  std::stringstream time_ss;
  time_ss << std::fixed << std::setprecision(20) << time_step;

  const std::string runs_string =
      "RUN_CELLS\n -cells 1\n -time_step " + time_ss.str() + "\nEND\n";
  this->RunString(runs_string.c_str());
}

void PhreeqcEngine::Impl::init_wrappers(const POETInitCell &cell) {

  // Solutions
  this->solutionWrapperPtr =
      std::make_unique<SolutionWrapper>(this->Get_solution(1), cell.solutions);

  if (this->Get_exchange(1) != nullptr) {
    this->exchangeWrapperPtr = std::make_unique<ExchangeWrapper>(
        this->Get_exchange(1), cell.exchanger);
    this->has_exchange = true;
  }

  if (this->Get_kinetic(1) != nullptr) {
    this->kineticsWrapperPtr =
        std::make_unique<KineticWrapper>(this->Get_kinetic(1), cell.kinetics);

    this->has_kinetics = true;
  }

  if (this->Get_equilibrium(1) != nullptr) {
    this->equilibriumWrapperPtr = std::make_unique<EquilibriumWrapper>(
        this->Get_equilibrium(1), cell.equilibrium);

    this->has_equilibrium = true;
  }

  if (this->Get_surface(1) != nullptr) {
    std::set<std::string> primaries(cell.solution_primaries.begin(),
                                    cell.solution_primaries.end());
    this->surfaceWrapperPtr = std::make_unique<SurfaceWrapper>(
        this->Get_surface(1), primaries, cell.surface_comps,
        cell.surface_charges);

    this->has_surface = true;
  }
}

void PhreeqcEngine::Impl::get_essential_values(std::span<double> &data) {

  this->solutionWrapperPtr->get(data);

  std::size_t offset = this->solutionWrapperPtr->size();

  if (this->has_exchange) {
    std::span<double> exch_span{
        data.subspan(offset, this->exchangeWrapperPtr->size())};
    this->exchangeWrapperPtr->get(exch_span);

    offset += this->exchangeWrapperPtr->size();
  }

  if (this->has_kinetics) {
    std::span<double> kin_span{
        data.subspan(offset, this->kineticsWrapperPtr->size())};
    this->kineticsWrapperPtr->get(kin_span);

    offset += this->kineticsWrapperPtr->size();
  }

  if (this->has_equilibrium) {
    std::span<double> equ_span{
        data.subspan(offset, this->equilibriumWrapperPtr->size())};
    this->equilibriumWrapperPtr->get(equ_span);

    offset += this->equilibriumWrapperPtr->size();
  }

  if (this->has_surface) {
    std::span<double> surf_span{
        data.subspan(offset, this->surfaceWrapperPtr->size())};
    this->surfaceWrapperPtr->get(surf_span);
  }
}

void PhreeqcEngine::Impl::set_essential_values(const std::span<double> &data) {

  this->solutionWrapperPtr->set(data);
  this->PhreeqcPtr->initial_solutions_poet(1);

  std::size_t offset = this->solutionWrapperPtr->size();

  if (this->has_exchange) {
    std::span<double> exch_span{
        data.subspan(offset, this->exchangeWrapperPtr->size())};
    this->exchangeWrapperPtr->set(exch_span);

    offset += this->exchangeWrapperPtr->size();
  }

  if (this->has_kinetics) {
    std::span<double> kin_span{
        data.subspan(offset, this->kineticsWrapperPtr->size())};
    this->kineticsWrapperPtr->set(kin_span);

    offset += this->kineticsWrapperPtr->size();
  }

  if (this->has_equilibrium) {
    std::span<double> equ_span{
        data.subspan(offset, this->equilibriumWrapperPtr->size())};
    this->equilibriumWrapperPtr->set(equ_span);

    offset += this->equilibriumWrapperPtr->size();
  }

  if (this->has_surface) {
    std::span<double> surf_span{
        data.subspan(offset, this->surfaceWrapperPtr->size())};
    this->surfaceWrapperPtr->set(surf_span);
  }
}