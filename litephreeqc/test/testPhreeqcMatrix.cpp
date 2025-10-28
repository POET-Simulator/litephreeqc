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

#include <cmath>
#include <gtest/gtest.h>
#include <string>
#include <vector>

#include <testInput.hpp>

#include "IPhreeqcReader.hpp"
#include "PhreeqcMatrix.hpp"
#include "utils.hpp"

const std::string base_db = readFile(base_test::phreeqc_database);

POET_TEST(PhreeqcInit) {
  EXPECT_NO_THROW(PhreeqcMatrix(base_db, base_test::script));
}

POET_TEST(PhreeqcMatrixOneSolution) {
  PhreeqcMatrix pqc_mat(base_db, base_test::script);
  const auto ids = pqc_mat.getIds();
  EXPECT_EQ(ids.size(), 1);
  EXPECT_EQ(ids[0], 1);

  EXPECT_TRUE(pqc_mat.checkIfExists(1));
  EXPECT_FALSE(pqc_mat.checkIfExists(2));

  PhreeqcMatrix::STLExport exported_init = pqc_mat.get();
  // ID + H,O,Charge + 6 Solutions + 4 Equil incl. params
  EXPECT_EQ(exported_init.names.size(), 19);

  IPhreeqcReader pqc_compare(base_db, base_test::script);
  pqc_compare.setOutputID(1);

  EXPECT_EQ(exported_init.names, base_test::expected_names);
  EXPECT_EQ(exported_init.values[0], 1);
  for (std::size_t i = 1; i < exported_init.values.size(); ++i) {
    EXPECT_NEAR(exported_init.values[i], pqc_compare[exported_init.names[i]],
                1e-7);
    // EXPECT_NEAR(exported_init.values[i], base_test::expected_values[i],
    //             base_test::expected_errors[i]);
  }

  auto dumps = pqc_mat.getDumpStringsPQI();
  EXPECT_EQ(dumps.size(), 1);
  EXPECT_GT(dumps[1].size(), 1);

  const auto kinetics = pqc_mat.getKineticsNames(1);
  EXPECT_EQ(kinetics.size(), 0);

  const auto equilibrium = pqc_mat.getEquilibriumNames(1);
  EXPECT_EQ(equilibrium.size(), 2);

  const auto expected_equilibrium =
      std::vector<std::string>({"Calcite", "Dolomite"});
  EXPECT_EQ(equilibrium, expected_equilibrium);
}

POET_TEST(PhreeqcMatrixBracketOperator) {
  PhreeqcMatrix pqc_mat(base_db, base_test::script);
  IPhreeqcReader pqc_compare(base_db, base_test::script);

  pqc_compare.setOutputID(1);

  EXPECT_NO_THROW(pqc_mat(1, "H"));
  EXPECT_NEAR(pqc_mat(1, "H"), pqc_compare["H"], 1e-7);
  EXPECT_ANY_THROW(pqc_mat(1, "J"));
  EXPECT_ANY_THROW(pqc_mat(2, "H"));
}

const std::string barite_db = readFile(barite_test::database);
const std::string barite_script = readFile(barite_test::script);

POET_TEST(PhreeqcMatrixMultiSolution) {
  PhreeqcMatrix pqc_mat(barite_db, barite_script);

  IPhreeqcReader pqc_compare(barite_db, barite_script);

  const auto ids = pqc_mat.getIds();
  EXPECT_EQ(ids.size(), 4);
  EXPECT_EQ(ids[0], 1);
  EXPECT_EQ(ids[1], 2);
  EXPECT_EQ(ids[2], 3);
  EXPECT_EQ(ids[3], 4);

  PhreeqcMatrix::STLExport exported = pqc_mat.get();

  EXPECT_EQ(exported.names, barite_test::expected_names);

  pqc_compare.setOutputID(1);

  for (std::size_t i = 1; i < exported.names.size(); i++) {
    if (i > 13 && i < 18) {
      EXPECT_TRUE(std::isnan(exported.values[i]));
      continue;
    }
    EXPECT_NEAR(exported.values[i], pqc_compare[exported.names[i]], 1e-7);
    // EXPECT_NEAR(exported.values[i], barite_test::expected_values_line_one[i],
    //             barite_test::expected_errors[i]);
  }

  auto dumps = pqc_mat.getDumpStringsPQI();
  EXPECT_EQ(dumps.size(), 4);

  const auto kinetics_sol_1 = pqc_mat.getKineticsNames(1);
  EXPECT_EQ(kinetics_sol_1.size(), 0);

  const auto equilibrium_sol_1 = pqc_mat.getEquilibriumNames(1);
  const auto expected_equilibrium_sol_1 =
      std::vector<std::string>({"Celestite"});
  EXPECT_EQ(equilibrium_sol_1, expected_equilibrium_sol_1);

  const auto equilibrium_sol_2 = pqc_mat.getEquilibriumNames(2);
  EXPECT_EQ(equilibrium_sol_2.size(), 0);

  const auto kinetics_sol_2 = pqc_mat.getKineticsNames(2);
  const auto expected_kinetics_sol_2 = std::vector<std::string>({"Celestite"});
  EXPECT_EQ(kinetics_sol_2, expected_kinetics_sol_2);

  const auto kinetics_sol_3 = pqc_mat.getKineticsNames(3);
  const auto expected_kinetics_sol_3 =
      std::vector<std::string>({"Barite", "Celestite"});
  EXPECT_EQ(kinetics_sol_3, expected_kinetics_sol_3);

  EXPECT_EQ(pqc_mat.getKineticsNames(4).size(), 0);
  EXPECT_EQ(pqc_mat.getEquilibriumNames(4).size(), 0);
}

POET_TEST(PhreeqcMatrixCtor) {
  PhreeqcMatrix pqc_mat(barite_db, barite_script);
  PhreeqcMatrix pqc_mat_copy(pqc_mat);
  PhreeqcMatrix pqc_mat_move(std::move(pqc_mat_copy));

  IPhreeqcReader pqc_compare(barite_db, barite_script);

  const auto ids = pqc_mat_move.getIds();
  EXPECT_EQ(ids.size(), 4);
  EXPECT_EQ(ids[0], 1);
  EXPECT_EQ(ids[1], 2);
  EXPECT_EQ(ids[2], 3);
  EXPECT_EQ(ids[3], 4);

  PhreeqcMatrix::STLExport exported = pqc_mat_move.get();

  EXPECT_EQ(exported.names, barite_test::expected_names);

  pqc_compare.setOutputID(1);

  for (std::size_t i = 1; i < exported.names.size(); i++) {
    if (i > 13 && i < 18) {
      EXPECT_TRUE(std::isnan(exported.values[i]));
      continue;
    }
    EXPECT_NEAR(exported.values[i], pqc_compare[exported.names[i]], 1e-7);
    // EXPECT_NEAR(exported.values[i], barite_test::expected_values_line_one[i],
    //             barite_test::expected_errors[i]);
  }
}

POET_TEST(PhreeqcMatrixOperator) {
  PhreeqcMatrix pqc_mat(barite_db, barite_script);
  PhreeqcMatrix pqc_mat_copy = pqc_mat;

  IPhreeqcReader pqc_compare(barite_db, barite_script);

  const auto ids = pqc_mat_copy.getIds();

  EXPECT_EQ(ids.size(), 4);
  EXPECT_EQ(ids[0], 1);
  EXPECT_EQ(ids[1], 2);
  EXPECT_EQ(ids[2], 3);
  EXPECT_EQ(ids[3], 4);

  PhreeqcMatrix::STLExport exported = pqc_mat_copy.get();

  EXPECT_EQ(exported.names, barite_test::expected_names);

  pqc_compare.setOutputID(1);

  for (std::size_t i = 1; i < exported.names.size(); i++) {
    if (i > 13 && i < 18) {
      EXPECT_TRUE(std::isnan(exported.values[i]));
      continue;
    }
    EXPECT_NEAR(exported.values[i], pqc_compare[exported.names[i]], 1e-7);
    // EXPECT_NEAR(exported.values[i], barite_test::expected_values_line_one[i],
    //             barite_test::expected_errors[i]);
  }
}

POET_TEST(PhreeqcMatrixRvalueManipulation) {
  PhreeqcMatrix pqc_mat(barite_db, barite_script);

  PhreeqcMatrix pqc_erased = pqc_mat.erase({1});

  const std::vector<int> expected_ids_erased = {2, 3, 4};

  PhreeqcMatrix::STLExport exported_erased = pqc_erased.get();

  EXPECT_EQ(pqc_erased.getIds(), expected_ids_erased);
  EXPECT_EQ(exported_erased.names, barite_test::expected_names_erased);
  EXPECT_EQ(exported_erased.values.size(),
            exported_erased.names.size() * expected_ids_erased.size());

  PhreeqcMatrix pqc_mat_subset = pqc_mat.subset({1});

  const std::vector<int> expected_ids_subset = {1};

  PhreeqcMatrix::STLExport exported_subset = pqc_mat_subset.get();

  EXPECT_EQ(pqc_mat_subset.getIds(), expected_ids_subset);
  EXPECT_EQ(exported_subset.names, barite_test::expected_names_subset);

  pqc_mat = pqc_mat.subset({1});

  exported_subset = pqc_mat.get();

  EXPECT_EQ(pqc_mat_subset.getIds(), expected_ids_subset);
  EXPECT_EQ(exported_subset.names, barite_test::expected_names_subset);
}

POET_TEST(PhreeqcMatrixColumnMajorExport) {
  PhreeqcMatrix pqc_mat(barite_db, barite_script);

  pqc_mat = pqc_mat.subset({2, 3});

  PhreeqcMatrix::STLExport exported =
      pqc_mat.get(PhreeqcMatrix::VectorExportType::COLUMN_MAJOR);

  const auto ids = pqc_mat.getIds();
  EXPECT_EQ(ids.size(), 2);
  EXPECT_EQ(ids[0], 2);
  EXPECT_EQ(ids[1], 3);

  EXPECT_EQ(exported.names, barite_test::expected_names_erased);
  EXPECT_EQ(exported.values.size(), exported.names.size() * ids.size());

  EXPECT_EQ(exported.values[0], 2);
  EXPECT_EQ(exported.values[1], 3);
}

POET_TEST(PhreeqcMatrixWithoutRedoxAndH0O0) {
  PhreeqcMatrix pqc_mat(barite_db, barite_script, false, false);

  const std::vector<std::string> expected_names_without_redox = {
      "H",  "O",  "Charge", "tc", "patm", "SolVol",
      "pH", "pe", "Ba",     "Cl", "S",    "Sr",
  };

  EXPECT_EQ(expected_names_without_redox, pqc_mat.getSolutionNames());
}
