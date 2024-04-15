#include "PhreeqcEngine.hpp"
#include <cstddef>
#include <span>

void PhreeqcEngine::init_wrappers(const POETInitCell &cell) {

  // TODO: Implement the rest of the wrappers. Currently supported: EXCHANGE

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

void PhreeqcEngine::get_essential_values(std::span<double> &data) {

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

  // // Exchange
  // if (this->Get_exchange(cell_number) != NULL) {
  //   std::vector<double> exch_values(this->exchangeWrapperPtr->size());

  //   std::span<double> exch_span{exch_values};

  //   this->exchangeWrapperPtr->get(exch_span);

  //   // this->Get_exchange(cell_number)->get_essential_values(exch_values);
  //   essentials.insert(essentials.end(), exch_values.begin(),
  //   exch_values.end());
  // }

  // // Kinetics
  // if (this->Get_kinetic(cell_number) != NULL) {
  //   std::vector<double> kin_values;

  //   this->Get_kinetic(cell_number)->get_essential_values(kin_values);
  //   essentials.insert(essentials.end(), kin_values.begin(),
  //   kin_values.end());
  // }

  // // PPassemblage
  // if (this->Get_equilibrium(cell_number) != NULL) {
  //   std::vector<double> equ_values;

  //   this->Get_equilibrium(cell_number)->get_essential_values(equ_values);
  //   essentials.insert(essentials.end(), equ_values.begin(),
  //   equ_values.end());
  // }

  // // Surface
  // if (this->Get_surface(cell_number) != NULL) {
  //   std::vector<double> surf_values(this->surfaceWrapperPtr->size());
  //   std::span<double> surf_span{surf_values};

  //   this->surfaceWrapperPtr->get(surf_span);

  //   // this->Get_surface(cell_number)->get_essential_values(surf_values);
  //   essentials.insert(essentials.end(), surf_values.begin(),
  //   surf_values.end());
  // }

  // return essentials;
}

void PhreeqcEngine::set_essential_values(const std::span<double> &data) {

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

  // // Solutions
  // if (this->Get_solution(cell_number) != NULL) {
  //   this->Get_solution(cell_number)->set_essential_values(dump_it, order);
  //   this->PhreeqcPtr->initial_solutions_poet(cell_number);
  // }

  // // Exchange
  // if (this->Get_exchange(cell_number) != NULL) {
  //   auto &exch_pointer = this->exchangeWrapperPtr;
  //   std::span<double> exch_span{dump_it, exch_pointer->size()};
  //   exch_pointer->set(exch_span);

  //   dump_it += exch_pointer->size();

  //   // this->Get_exchange(cell_number)->set_essential_values(dump_it);
  //   // this->PhreeqcPtr->Rxn_new_exchange.insert(cell_number);
  //   // this->Get_exchange(cell_number)->Set_new_def(true);
  //   // this->Get_exchange(cell_number)->Set_solution_equilibria(true);
  //   // this->Get_exchange(cell_number)->Set_n_solution(cell_number);
  //   // this->PhreeqcPtr->initial_exchangers(FALSE);
  // }

  // // Kinetics
  // if (this->Get_kinetic(cell_number) != NULL) {
  //   this->Get_kinetic(cell_number)->set_essential_values(dump_it);
  // }

  // // PPassemblage
  // if (this->Get_equilibrium(cell_number) != NULL) {
  //   this->Get_equilibrium(cell_number)->set_essential_values(dump_it);
  // }

  // // Surface
  // if (this->Get_surface(cell_number) != NULL) {
  //   auto *tmp = this->Get_surface(cell_number);
  //   auto &surf_pointer = this->surfaceWrapperPtr;
  //   std::span<double> surf_span{dump_it, surf_pointer->size()};
  //   surf_pointer->set(surf_span);

  //   dump_it += surf_pointer->size();

  //   // this->PhreeqcPtr->Rxn_new_surface.insert(cell_number);

  //   // this->Get_surface(cell_number)->set_essential_values(dump_it);
  //   // this->PhreeqcPtr->Rxn_new_surface.insert(cell_number);
  //   // this->PhreeqcPtr->initial_surfaces(TRUE);
  // }
}

void PhreeqcEngine::runCell(std::vector<double> &cell_values,
                            double time_step) {
  // skip ID
  std::span<double> cell_data{cell_values.begin() + 1, cell_values.end()};

  this->set_essential_values(cell_data);
  this->run(1, 1, time_step);
  this->get_essential_values(cell_data);
}