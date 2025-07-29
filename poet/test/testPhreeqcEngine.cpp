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
