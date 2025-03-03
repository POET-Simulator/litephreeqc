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

#include <stdexcept>

#include <testInput.hpp>

#include <gtest/gtest.h>

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

  std::vector<double> cell_values = pqc_mat.get().values;
  cell_values.erase(cell_values.begin(), cell_values.begin() + 1);

  EXPECT_NO_THROW(engine.runCell(cell_values, 0));
  EXPECT_NO_THROW(engine.runCell(cell_values, 100));

  for (std::size_t i = 0; i < cell_values.size(); ++i) {
    // skip Charge, H(0) and O(0)
    if (i >= 2 && i <= 4) {
      continue;
    }
    EXPECT_NEAR(cell_values[i], base_test::expected_values[i + 1],
                base_test::expected_errors[i + 1]);
  }

  EXPECT_THROW(engine.runCell(cell_values, -1), std::invalid_argument);
}