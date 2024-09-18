#pragma once

#include "PhreeqcMatrix.hpp"
#include <memory>
#include <vector>

/**
 * @brief Class for running Phreeqc wrappped in POET
 *
 * Direct interface to Phreeqc, without utilizing Phreeqc's internal parser
 * to set/get values. To start a simulation, the interpreter is used internally.
 * No need for the user to directly interact with the interpreter.
 *
 */
class PhreeqcEngine {
public:
  /**
   * @brief Construct a new Phreeqc Engine object
   *
   * Construct a new Phreeqc Engine object by previously initialized
   * PhreeqcMatrix.
   *
   * @param pqc_mat PhreeqcMatrix initialized with the *full* Phreeqc script and
   * database.
   * @param cell_id ID of the cell (user id from Phreeqc script) to simulate
   */
  PhreeqcEngine(const PhreeqcMatrix &pqc_mat, const int cell_id);

  /**
   * @brief Destroy the Phreeqc Engine object
   *
   */
  ~PhreeqcEngine();

  /**
   * @brief Siimulate a cell for a given time step
   *
   * @param cell_values Vector containing the input values for the cell
   * (*including the ID*). Output values are written back in place to this
   * vector.
   * @param time_step Time step to simulate in seconds
   */
  void runCell(std::vector<double> &cell_values, double time_step);

private:
  class Impl;
  std::unique_ptr<Impl> impl;
};