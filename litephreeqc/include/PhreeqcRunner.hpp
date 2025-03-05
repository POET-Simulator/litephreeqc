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

#pragma once

#include "PhreeqcEngine.hpp"
#include "PhreeqcMatrix.hpp"
#include <cstddef>
#include <memory>
#include <unordered_map>
#include <vector>

/**
 * @class PhreeqcRunner
 * @brief Manages the execution of Phreeqc simulations.
 *
 * The PhreeqcRunner class is responsible for managing and running Phreeqc
 * simulations. It initializes a set of PhreeqcEngine instances based on the
 * provided PhreeqcMatrix and provides functionality to run simulations and
 * retrieve the number of engines.
 *
 * @note Copy and move operations are deleted to prevent unintended behavior.
 */
class PhreeqcRunner {
public:
  PhreeqcRunner(const PhreeqcRunner &) = delete;
  PhreeqcRunner(PhreeqcRunner &&) = delete;
  PhreeqcRunner &operator=(const PhreeqcRunner &) = delete;
  PhreeqcRunner &operator=(PhreeqcRunner &&) = delete;

  /**
   * @brief Constructs a PhreeqcRunner object with the given PhreeqcMatrix.
   *
   * Constructs a PhreeqcRunner object with the given PhreeqcMatrix. For each
   *  cell found in the PhreeqcMatrix, a PhreeqcEngine is created and stored in
   *  an internal map.
   *
   * @param matrix A reference to a PhreeqcMatrix object used to initialize the
   * PhreeqcRunner.
   */
  PhreeqcRunner(const PhreeqcMatrix &matrix);
  ~PhreeqcRunner() = default;

  /**
   * @brief Runs the simulation with the given input and output data for a
   * specified time step.
   *
   * @param simulationInOut A reference to a 2D vector containing the input data
   * for the simulation. The vector will be modified to contain the output data
   * after the simulation.
   * @param time_step The time step for the simulation.
   */
  void run(std::vector<std::vector<double>> &simulationInOut,
           const double time_step);

  /**
   * @brief Runs the simulation with the given input and output data.
   *
   * @param simulationInOut A reference to a 2D vector containing the simulation
   * input and output data.
   * @param time_step The time step for the simulation.
   * @param to_ignore A vector of indices specifying which elements to ignore
   * during the simulation.
   */
  void run(std::vector<std::vector<double>> &simulationInOut,
           const double time_step, const std::vector<std::size_t> &to_ignore);

  /**
   * @brief Returns the number of engines currently stored.
   *
   * This function provides the count of engines that are currently
   * stored in the _engineStorage container.
   *
   * @return std::size_t The number of engines.
   */
  std::size_t numEngines() const { return _engineStorage.size(); }

private:
  std::unordered_map<int, std::unique_ptr<PhreeqcEngine>> _engineStorage;
  std::vector<double> _buffer;
};