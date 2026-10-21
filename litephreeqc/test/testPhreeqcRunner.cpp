/*
 * Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
 *                       Marco De Lucia (GFZ German Research Centre for Geosciences)
 *
 * SPDX-License-Identifier: EUPL-1.2
 *
 * This file is part of litephreeqc, a C++ interface library on top of PHREEQC.
 */
#include "PhreeqcRunner.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <string>
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
      EXPECT_TRUE(std::isnan(simulationInOut[cell_index][17]));  // Barite_kin (index shifted by 2 for Viscosity, Density)
    } else {
      EXPECT_EQ(simulationInOut[cell_index][0], expected_value_second_half);
      EXPECT_FALSE(std::isnan(simulationInOut[cell_index][17]));  // Barite_kin
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

// The text column state_SO is NaN in every cell, and cell 1 lacks the PURE
// columns of cell 2: the runner must take the values of a cell by column, not
// by skipping NaNs (which shifted the SELECTED_OUTPUT values after state_SO)
POET_TEST(PhreeqcRunnerSelectedOutputTextColumn) {
  const std::string script = R"(SELECTED_OUTPUT
-high_precision true
-reset true
SOLUTION 1
pH 7.0
SOLUTION 2
units mol/kgw
pH 7
Na 0.1
Cl 0.1 charge
PURE 2
Calcite 0.0 1
RUN_CELLS
-cells 1 2
END)";
  PhreeqcMatrix pqc_mat(readFile(base_test::phreeqc_database), script);
  const auto exported = pqc_mat.get();
  const auto &names = exported.names;
  const std::size_t num_columns = names.size();
  auto col = [&names](const std::string &name) {
    return static_cast<std::size_t>(std::distance(
        names.begin(), std::find(names.begin(), names.end(), name)));
  };
  const std::size_t state = col("state_SO");
  const std::size_t ph = col("pH");
  const std::size_t ph_so = col("pH_SO");
  const std::size_t calcite = col("Calcite_eq");
  ASSERT_LT(state, num_columns);
  ASSERT_LT(ph, num_columns);
  ASSERT_LT(ph_so, num_columns);
  ASSERT_LT(calcite, num_columns);

  const auto cols_1 = pqc_mat.getCellColumns(1);
  const auto cols_2 = pqc_mat.getCellColumns(2);
  auto has = [](const std::vector<std::size_t> &cols, std::size_t j) {
    return std::find(cols.begin(), cols.end(), j) != cols.end();
  };
  EXPECT_TRUE(has(cols_1, state));
  EXPECT_TRUE(has(cols_2, state));
  EXPECT_FALSE(has(cols_1, calcite));
  EXPECT_TRUE(has(cols_2, calcite));

  std::vector<std::vector<double>> simulationInOut = {
      {exported.values.begin(), exported.values.begin() + num_columns},
      {exported.values.begin() + num_columns, exported.values.end()}};
  ASSERT_TRUE(std::isnan(simulationInOut[1][state]));

  PhreeqcRunner runner(pqc_mat);
  ASSERT_NO_THROW(runner.run(simulationInOut, 0));

  EXPECT_TRUE(std::isnan(simulationInOut[0][calcite]));
  EXPECT_FALSE(std::isnan(simulationInOut[1][calcite]));
  for (const auto &row : simulationInOut) {
    EXPECT_TRUE(std::isnan(row[state]));
    EXPECT_NEAR(row[ph_so], row[ph], 1e-8);
  }

  // a row of the wrong length is an error, not a misaligned cell
  simulationInOut[0].pop_back();
  EXPECT_THROW(runner.run(simulationInOut, 0), std::invalid_argument);
}
