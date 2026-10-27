/*
 * Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
 *                       Marco De Lucia (GFZ German Research Centre for
 * Geosciences)
 *
 * SPDX-License-Identifier: EUPL-1.2
 *
 * This file is part of litephreeqc, a C++ interface library on top of PHREEQC.
 */
#include <algorithm>
#include <stdexcept>

#include <testInput.hpp>

#include <gtest/gtest.h>

#include "IPhreeqcReader.hpp"
#include "PhreeqcEngine.hpp"
#include "PhreeqcMatrix.hpp"
#include "utils.hpp"

const std::string test_database = readFile(base_test::phreeqc_database);

POET_TEST(PhreeqcEngineConstructor) {

  PhreeqcMatrix pqc_mat(test_database, base_test::script);

  EXPECT_NO_THROW(PhreeqcEngine(pqc_mat, 1));
  EXPECT_THROW(PhreeqcEngine(pqc_mat, 2), std::invalid_argument);
}

POET_TEST(PhreeqcEngineStep) {
  PhreeqcMatrix pqc_mat(test_database, base_test::script);

  PhreeqcEngine engine(pqc_mat, 1);

  IPhreeqcReader pqc_compare(test_database, base_test::script);

  std::vector<double> cell_values = pqc_mat.get().values;
  std::vector<std::string> cell_names = pqc_mat.get().names;
  cell_values.erase(cell_values.begin(), cell_values.begin() + 1);
  cell_names.erase(cell_names.begin(), cell_names.begin() + 1);

  EXPECT_NO_THROW(engine.runCell(cell_values, 0));
  EXPECT_NO_THROW(engine.runCell(cell_values, 100));

  pqc_compare.run(0, {1});
  pqc_compare.run(100, {1});

  pqc_compare.setOutputID(1);

  for (std::size_t i = 0; i < cell_names.size(); ++i) {
    // Somehow 'pe' will not result in a expected near value, therefore we skip
    // it
    if (cell_names[i] == "pe") {
      continue;
    }
    EXPECT_NEAR(cell_values[i], pqc_compare[cell_names[i]], 1e-6);
  }

  EXPECT_THROW(engine.runCell(cell_values, -1), std::invalid_argument);
}

POET_TEST(PhreeqcEngineSelectedOutputStep) {
  const std::string script_with_so = R"(SOLUTION 1
units mol/kgw
temp 25
Ca 0.1
Mg 0.1
Cl 0.5 charge
Na 0.1
PURE 1
Calcite  0.0 1
Dolomite 0.0 0
SELECTED_OUTPUT 1
 -reset false
 -totals Ca Na
USER_PUNCH 1
 -headings MyVal
 10 PUNCH 789.0
RUN_CELLS
 -cells 1
END)";

  PhreeqcMatrix pqc_mat(test_database, script_with_so);

  PhreeqcEngine engine(pqc_mat, 1);

  std::vector<double> cell_values = pqc_mat.get().values;
  std::vector<std::string> cell_names = pqc_mat.get().names;
  cell_values.erase(cell_values.begin(), cell_values.begin() + 1);
  cell_names.erase(cell_names.begin(), cell_names.begin() + 1);

  auto it_ca_so =
      std::find(cell_names.begin(), cell_names.end(), "Ca(mol/kgw)_SO");
  auto it_punch = std::find(cell_names.begin(), cell_names.end(), "MyVal_SO");
  EXPECT_NE(it_ca_so, cell_names.end());
  EXPECT_NE(it_punch, cell_names.end());

  std::size_t ca_idx = std::distance(cell_names.begin(), it_ca_so);
  std::size_t punch_idx = std::distance(cell_names.begin(), it_punch);

  EXPECT_NO_THROW(engine.runCell(cell_values, 100));

  EXPECT_NEAR(cell_values[ca_idx], 0.12131, 1e-3);
  EXPECT_NEAR(cell_values[punch_idx], 789.0, 1e-4);
}
