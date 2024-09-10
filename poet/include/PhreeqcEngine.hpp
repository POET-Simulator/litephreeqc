#pragma once

#include "POETInit.hpp"

#include <memory>
#include <vector>

/**
 * @brief Class for running Phreeqc wrappped in POET
 *
 * Direct interface to Phreeqc, without utilizing *_MODIFY keywords/scripts
 to
 * set new values. Use with already initialized Phreeqc config.
 *
 */
class PhreeqcEngine {
public:
  /**
   * @brief Construct a new Phreeqc Engine object
   *
   * Construct a new Phreeqc Engine object by previously initialized
   POETConfig.
   *
   * @param config Holds the configuration for the Phreeqc engine.
   */
  PhreeqcEngine(const POETConfig &config);

  /**
   * @brief Destroy the Phreeqc Engine object
   *
   */
  ~PhreeqcEngine();

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
  class Impl;
  std::unique_ptr<Impl> impl;
};