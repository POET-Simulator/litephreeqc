#include <IPhreeqcPOET.hpp>

IPhreeqcPOET::essential_names
IPhreeqcPOET::dump_essential_names(std::size_t cell_number) {

  essential_names eNames;

  // Solutions
  if (this->Get_solution(cell_number) != NULL) {
    std::vector<std::string> &eSolNames = eNames[POET_SOL];
    this->Get_solution(cell_number)->dump_essential_names(eSolNames);
  }

  // Exchange
  if (this->Get_exchange(cell_number) != NULL) {
    std::vector<std::string> &eExchNames = eNames[POET_EXCH];
    this->Get_exchange(cell_number)->dump_essential_names(eExchNames);
  }

  // Kinetics
  if (this->Get_kinetic(cell_number) != NULL) {
    std::vector<std::string> &eKinNames = eNames[POET_KIN];
    this->Get_kinetic(cell_number)->dump_essential_names(eKinNames);
  }

  // PPassemblage
  if (this->Get_equilibrium(cell_number) != NULL) {
    std::vector<std::string> &eEquNames = eNames[POET_EQUIL];
    this->Get_equilibrium(cell_number)->dump_essential_names(eEquNames);
  }

  // Surface
  if (this->Get_surface(cell_number) != NULL) {
    std::vector<std::string> &eSurfNames = eNames[POET_SURF];
    this->Get_surface(cell_number)->dump_essential_names(eSurfNames);
  }

  return eNames;
}

std::vector<double>
IPhreeqcPOET::get_essential_values(std::size_t cell_number,
                                   const std::vector<std::string> &order) {
  std::vector<double> essentials;

  // Solutions
  if (this->Get_solution(cell_number) != NULL) {
    std::vector<double> sol_values;

    this->Get_solution(cell_number)->get_essential_values(sol_values, order);
    essentials.insert(essentials.end(), sol_values.begin(), sol_values.end());
  }

  // Exchange
  if (this->Get_exchange(cell_number) != NULL) {
    std::vector<double> exch_values;

    this->Get_exchange(cell_number)->get_essential_values(exch_values);
    essentials.insert(essentials.end(), exch_values.begin(), exch_values.end());
  }

  // Kinetics
  if (this->Get_kinetic(cell_number) != NULL) {
    std::vector<double> kin_values;

    this->Get_kinetic(cell_number)->get_essential_values(kin_values);
    essentials.insert(essentials.end(), kin_values.begin(), kin_values.end());
  }

  // PPassemblage
  if (this->Get_equilibrium(cell_number) != NULL) {
    std::vector<double> equ_values;

    this->Get_equilibrium(cell_number)->get_essential_values(equ_values);
    essentials.insert(essentials.end(), equ_values.begin(), equ_values.end());
  }

  // Surface
  if (this->Get_surface(cell_number) != NULL) {
    std::vector<double> surf_values;

    this->Get_surface(cell_number)->get_essential_values(surf_values);
    essentials.insert(essentials.end(), surf_values.begin(), surf_values.end());
  }

  return essentials;
}

void IPhreeqcPOET::set_essential_values(std::size_t cell_number,
                                        const std::vector<std::string> &order,
                                        std::vector<double> &values) {

  auto dump_it = values.begin();

  // Solutions
  if (this->Get_solution(cell_number) != NULL) {
    this->Get_solution(cell_number)->set_essential_values(dump_it, order);
  }

  // Exchange
  if (this->Get_exchange(cell_number) != NULL) {
    this->Get_exchange(cell_number)->set_essential_values(dump_it);
  }

  // Kinetics
  if (this->Get_kinetic(cell_number) != NULL) {
    this->Get_kinetic(cell_number)->set_essential_values(dump_it);
  }

  // PPassemblage
  if (this->Get_equilibrium(cell_number) != NULL) {
    this->Get_equilibrium(cell_number)->set_essential_values(dump_it);
  }

  // Surface
  if (this->Get_surface(cell_number) != NULL) {
    this->Get_surface(cell_number)->set_essential_values(dump_it);
  }
}