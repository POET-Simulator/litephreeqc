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

  EXPECT_NO_THROW(engine.runCell(cell_values, 0));
  EXPECT_NO_THROW(engine.runCell(cell_values, 100));

  for (std::size_t i = 0; i < cell_values.size(); ++i) {
    // ignore H(0) and O(0)
    if (i == 4 || i == 5) {
      continue;
    }
    EXPECT_NEAR(cell_values[i], base_test::expected_values[i],
                base_test::expected_errors[i]);
  }
}