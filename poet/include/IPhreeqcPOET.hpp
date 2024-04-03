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
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <utility>
#include <vector>

enum { POET_SOL = 0, POET_EXCH, POET_KIN, POET_EQUIL, POET_SURF };

class IPhreeqcPOET : public IPhreeqc {
public:
  IPhreeqcPOET(const std::string &database, const std::string &input_script)
      : IPhreeqc() {
    this->LoadDatabaseString(database.c_str());
    this->RunString(input_script.c_str());

    this->parseInit();
  }

  IPhreeqcPOET(const std::string &database, const std::string &input_script,
               const std::vector<std::string> &solutionInitVector,
               std::uint32_t n_cells)
      : IPhreeqc(), n_cells(n_cells), solutionInitVector(solutionInitVector) {
    this->LoadDatabaseString(database.c_str());
    this->RunString(input_script.c_str());

    if (n_cells > 1) {
      std::string copy_string =
          "COPY cell 1 2-" + std::to_string(n_cells) + "\n";
      this->RunString(copy_string.c_str());
    }
  }

  void queueCell(std::vector<double> &values) {
    if (this->queued_cell_pointer > this->n_cells) {
      return;
    }

    setCell(this->queued_cell_pointer++, values);
  }

  void dequeueCells(std::vector<std::vector<double>> &values) {
    if (this->queued_cell_pointer == 1) {
      return;
    }

    for (std::size_t i = 1; i < this->queued_cell_pointer; i++) {
      std::vector<double> cell_values;
      getCell(i, cell_values);
      values.push_back(cell_values);
    }

    this->queued_cell_pointer = 1;
  }

  void runQueuedCells(double time_step) {
    if (this->queued_cell_pointer == 1) {
      return;
    }
    run(1, this->queued_cell_pointer - 1, time_step);
  }

  void setCell(int cell_number, std::vector<double> &values) {
    this->set_essential_values(cell_number, this->solutionInitVector, values);
  }

  void getCell(int cell_number, std::vector<double> &values) {
    values = this->get_essential_values(cell_number, this->solutionInitVector);
  }

  void run(int start_cell, int end_cell, double time_step) {
    std::stringstream time_ss;
    time_ss << std::fixed << std::setprecision(20) << time_step;

    const std::string runs_string =
        "RUN_CELLS\n -cells " + std::to_string(start_cell) + "-" +
        std::to_string(end_cell) + "\n -time_step " + time_ss.str();
    this->RunString(runs_string.c_str());
  }

  struct PhreeqcMat {
    std::vector<std::string> names;
    std::vector<int> ids;
    std::vector<std::vector<double>> values;
  };

  PhreeqcMat getPhreeqcMat();

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

  ModulesArray getModuleSizes(const std::vector<int> &cell_ids);

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

  std::vector<double>
  get_essential_values(std::size_t cell_number,
                       const std::vector<std::string> &order);

  void set_essential_values(std::size_t cell_number,
                            const std::vector<std::string> &order,
                            std::vector<double> &values);

  RawMap raw_initials;

  std::size_t queued_cell_pointer = 1;
};