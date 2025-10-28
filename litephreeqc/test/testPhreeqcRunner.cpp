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

#include "PhreeqcRunner.hpp"
#include "utils.hpp"

#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <testInput.hpp>
#include <vector>

const std::string test_database = readFile(barite_test::database);
const std::string test_script = readFile(barite_test::script);

constexpr std::size_t num_cells = 10;

POET_TEST(PhreeqcRunnerConstructor) {
  PhreeqcMatrix pqc_mat(test_database, test_script);
  EXPECT_NO_THROW(PhreeqcRunner tmp(pqc_mat));
}

POET_TEST(PhreeqcRunnerSimulation) {
  PhreeqcMatrix pqc_mat(test_database, test_script);
  const auto subsetted_pqc_mat = pqc_mat.subset({2, 3});
  PhreeqcRunner runner(subsetted_pqc_mat);

  const auto stl_mat = subsetted_pqc_mat.get();
  const auto matrix_values = stl_mat.values;
  const auto num_columns = stl_mat.names.size();

  std::vector<std::vector<double>> simulationInOut;

  for (std::size_t index = 0; index < num_cells; ++index) {
    if (index < static_cast<std::size_t>(num_cells / 2)) {
      simulationInOut.push_back(std::vector<double>(
          matrix_values.begin(), matrix_values.begin() + num_columns));
    } else {
      simulationInOut.push_back(std::vector<double>(
          matrix_values.begin() + num_columns, matrix_values.end()));
    }
  }
  EXPECT_NO_THROW(runner.run(simulationInOut, 100));

  constexpr std::size_t half_cells = num_cells / 2;
  constexpr int expected_value_first_half = 2;
  constexpr int expected_value_second_half = 3;

  for (std::size_t cell_index = 0; cell_index < simulationInOut.size();
       ++cell_index) {
    const bool is_first_half = cell_index < half_cells;
    if (is_first_half) {
      EXPECT_EQ(simulationInOut[cell_index][0], expected_value_first_half);
      EXPECT_TRUE(std::isnan(simulationInOut[cell_index][14]));
    } else {
      EXPECT_EQ(simulationInOut[cell_index][0], expected_value_second_half);
      EXPECT_FALSE(std::isnan(simulationInOut[cell_index][14]));
    }

    EXPECT_NEAR(simulationInOut[cell_index][1], 111, 1);
    EXPECT_NEAR(simulationInOut[cell_index][2], 55, 1);
  }
}

POET_TEST(PhreeqcRunnerUnknownID) {
  PhreeqcMatrix pqc_mat(test_database, test_script);
  PhreeqcRunner runner(pqc_mat);

  std::vector<std::vector<double>> simulationInOut;

  for (std::size_t index = 0; index < num_cells; ++index) {
    simulationInOut.push_back(std::vector<double>(pqc_mat.get().names.size()));
  }

  simulationInOut[0][0] = 1000;

  EXPECT_THROW(runner.run(simulationInOut, 100), std::out_of_range);
}

POET_TEST(PhreeqcRunnerSimulationWithIgnoredCells) {
  PhreeqcMatrix pqc_mat(test_database, test_script);
  const auto subsetted_pqc_mat = pqc_mat.subset({2, 3});
  PhreeqcRunner runner(subsetted_pqc_mat);

  const auto stl_mat = subsetted_pqc_mat.get();
  const auto matrix_values = stl_mat.values;
  const auto num_columns = stl_mat.names.size();

  std::vector<double> second_line(matrix_values.begin() + num_columns,
                                  matrix_values.end());
  std::vector<std::vector<double>> simulationInOut;

  simulationInOut.push_back(second_line);

  EXPECT_NO_THROW(runner.run(simulationInOut, 10000, {0}));

  for (std::size_t i = 0; i < num_columns; ++i) {
    EXPECT_DOUBLE_EQ(simulationInOut[0][i], second_line[i]);
  }
}
