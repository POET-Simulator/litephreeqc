/*
 * This project is subject to the original PHREEQC license. `litephreeqc` is a
 * version of the PHREEQC code that has been modified to be used as a library.
 *
 * It adds a C++ interface on top of the original PHREEQC code, with small
 * changes to the original code base.
 *
 * Authors of Modifications:
 * - Max Luebke (mluebke@uni-potsdam.de) - University of Potsdam
 * - Marco De Lucia (delucia@gfz.de) - GFZ Helmholz Centre for Geosciences
 *
 */

#include "PhreeqcEngine.hpp"
#include <cstddef>
#include <iomanip>
#include <regex>
#include <span>
#include <sstream>

#include <IPhreeqc.hpp>
#include <Phreeqc.h>
#include <string>
#include <vector>

#include "PhreeqcKnobs.hpp"
#include "PhreeqcMatrix.hpp"
#include "Wrapper/EquilibriumWrapper.hpp"
#include "Wrapper/ExchangeWrapper.hpp"
#include "Wrapper/KineticWrapper.hpp"
#include "Wrapper/SolutionWrapper.hpp"
#include "Wrapper/SurfaceWrapper.hpp"
class PhreeqcEngine::Impl : public IPhreeqc {
public:
  Impl(const PhreeqcMatrix &pqc_mat, const int cell_id);
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
  struct InitCell {
    std::vector<std::string> solutions;
    bool with_redox;
    std::vector<std::string> exchanger;
    std::vector<std::string> kinetics;
    std::vector<std::string> equilibrium;
    std::vector<std::string> surface_comps;
    std::vector<std::string> surface_charges;
    std::vector<std::string> solution_primaries;
  };
  void init_wrappers(const InitCell &cell);
};

PhreeqcEngine::PhreeqcEngine(const PhreeqcMatrix &pqc_mat, const int cell_id)
    : impl(std::make_unique<Impl>(pqc_mat, cell_id)) {}

PhreeqcEngine::~PhreeqcEngine() = default;

static inline std::string
replaceRawKeywordID(const std::string &raw_dump_string) {
  std::regex re(R"((RAW\s+)(\d+))");
  return std::regex_replace(raw_dump_string, re, "RAW 1");
}

PhreeqcEngine::Impl::Impl(const PhreeqcMatrix &pqc_mat, const int cell_id) {

  if (!pqc_mat.checkIfExists(cell_id)) {
    throw std::invalid_argument("Cell ID does not exist in PhreeqcMatrix");
  }

  this->LoadDatabaseString(pqc_mat.getDatabase().c_str());

  pqc_mat.getKnobs().writeKnobs(this->PhreeqcPtr);

  const std::string pqc_string =
      replaceRawKeywordID(pqc_mat.getDumpStringsPQI(cell_id));

  this->RunString(pqc_string.c_str());

  InitCell cell = {pqc_mat.getSolutionNames(),
                   pqc_mat.withRedox(),
                   pqc_mat.getExchanger(cell_id),
                   pqc_mat.getKineticsNames(cell_id),
                   pqc_mat.getEquilibriumNames(cell_id),
                   pqc_mat.getSurfaceCompNames(cell_id),
                   pqc_mat.getSurfaceChargeNames(cell_id),
                   pqc_mat.getSolutionPrimaries()};

  this->init_wrappers(cell);
}

void PhreeqcEngine::runCell(std::vector<double> &cell_values,
                            double time_step) {

  if (time_step < 0) {
    throw std::invalid_argument("Time step must be positive");
  }

  // ID is already skipped by PhreeqcRunner, so no need to start ahead
  std::span<double> cell_data{cell_values.begin(), cell_values.end()};

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

  if (this->GetErrorStringLineCount() > 0) {
    std::cerr << ":: Error in Phreeqc script: " << this->GetErrorString()
              << "\n";
    throw std::runtime_error("Phreeqc script error");
  }
}

void PhreeqcEngine::Impl::init_wrappers(const InitCell &cell) {

  // Solutions
  this->solutionWrapperPtr = std::make_unique<SolutionWrapper>(
      this->Get_solution(1), cell.solutions, cell.with_redox);

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
  // this->PhreeqcPtr->initial_solutions_poet(1);

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
