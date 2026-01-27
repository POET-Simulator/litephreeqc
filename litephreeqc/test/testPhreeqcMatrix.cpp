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
  // ID + H,O,Charge,tc,patm,SolVol,pH,pe,MassH2O,Viscosity,Density + 6 Solutions + 4 Equil incl. params
  EXPECT_EQ(exported_init.names.size(), 22);

  IPhreeqcReader pqc_compare(base_db, base_test::script);
  pqc_compare.setOutputID(1);

  EXPECT_EQ(exported_init.names, base_test::expected_names);
  EXPECT_EQ(exported_init.values[0], 1);
  for (std::size_t i = 1; i < exported_init.values.size(); ++i) {
    // Skip Viscosity and Density as IPhreeqcReader doesn't provide them
    if (exported_init.names[i] == "Viscosity" ||
        exported_init.names[i] == "Density") {
      continue;
    }
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
    if (i > 16 && i < 21) {
      EXPECT_TRUE(std::isnan(exported.values[i]));
      continue;
    }
    // Skip Viscosity and Density as IPhreeqcReader doesn't provide them
    if (exported.names[i] == "Viscosity" || exported.names[i] == "Density") {
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
    if (i > 16 && i < 21) {
      EXPECT_TRUE(std::isnan(exported.values[i]));
      continue;
    }
    // Skip Viscosity and Density as IPhreeqcReader doesn't provide them
    if (exported.names[i] == "Viscosity" || exported.names[i] == "Density") {
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
    if (i > 16 && i < 21) {
      EXPECT_TRUE(std::isnan(exported.values[i]));
      continue;
    }
    // Skip Viscosity and Density as IPhreeqcReader doesn't provide them
    if (exported.names[i] == "Viscosity" || exported.names[i] == "Density") {
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
      "H",  "O",       "Charge", "tc", "patm", "SolVol", "pH",
      "pe", "MassH2O", "Viscosity", "Density", "Ba",     "Cl", "S",    "Sr",
  };

  EXPECT_EQ(expected_names_without_redox, pqc_mat.getSolutionNames());
}

POET_TEST(PhreeqcMatrixOptionalEssentials) {
  using OE = PhreeqcMatrix::OptionalEssentials;

  // Test with only pH, pe, and Viscosity selected
  PhreeqcMatrix pqc_mat(barite_db, barite_script, false, true,
                        OE::pH | OE::pe | OE::Viscosity);

  // Verify the optionalEssentials accessor returns the correct value
  EXPECT_EQ(pqc_mat.optionalEssentials(), OE::pH | OE::pe | OE::Viscosity);

  // Get solution names and verify they contain only the selected optional
  // essentials
  const auto solution_names = pqc_mat.getSolutionNames();

  // Should contain mandatory essentials: H, O, Charge, tc, patm
  EXPECT_NE(std::find(solution_names.begin(), solution_names.end(), "H"),
            solution_names.end());
  EXPECT_NE(std::find(solution_names.begin(), solution_names.end(), "O"),
            solution_names.end());
  EXPECT_NE(std::find(solution_names.begin(), solution_names.end(), "Charge"),
            solution_names.end());
  EXPECT_NE(std::find(solution_names.begin(), solution_names.end(), "tc"),
            solution_names.end());
  EXPECT_NE(std::find(solution_names.begin(), solution_names.end(), "patm"),
            solution_names.end());

  // Should contain selected optional essentials: pH, pe, Viscosity
  EXPECT_NE(std::find(solution_names.begin(), solution_names.end(), "pH"),
            solution_names.end());
  EXPECT_NE(std::find(solution_names.begin(), solution_names.end(), "pe"),
            solution_names.end());
  EXPECT_NE(std::find(solution_names.begin(), solution_names.end(), "Viscosity"),
            solution_names.end());

  // Should NOT contain excluded optional essentials: SolVol, MassH2O, Density
  EXPECT_EQ(std::find(solution_names.begin(), solution_names.end(), "SolVol"),
            solution_names.end());
  EXPECT_EQ(std::find(solution_names.begin(), solution_names.end(), "MassH2O"),
            solution_names.end());
  EXPECT_EQ(std::find(solution_names.begin(), solution_names.end(), "Density"),
            solution_names.end());

  // Verify getMatrixOutOnly returns only enabled optionals (plus tc, patm)
  const auto out_only = pqc_mat.getMatrixOutOnly();
  EXPECT_NE(std::find(out_only.begin(), out_only.end(), "tc"), out_only.end());
  EXPECT_NE(std::find(out_only.begin(), out_only.end(), "patm"), out_only.end());
  EXPECT_NE(std::find(out_only.begin(), out_only.end(), "pH"), out_only.end());
  EXPECT_NE(std::find(out_only.begin(), out_only.end(), "pe"), out_only.end());
  EXPECT_NE(std::find(out_only.begin(), out_only.end(), "Viscosity"),
            out_only.end());
  EXPECT_EQ(std::find(out_only.begin(), out_only.end(), "SolVol"),
            out_only.end());
  EXPECT_EQ(std::find(out_only.begin(), out_only.end(), "MassH2O"),
            out_only.end());
  EXPECT_EQ(std::find(out_only.begin(), out_only.end(), "Density"),
            out_only.end());

  // Verify getMatrixTransported excludes enabled optionals
  const auto transported = pqc_mat.getMatrixTransported();
  EXPECT_EQ(std::find(transported.begin(), transported.end(), "tc"),
            transported.end());
  EXPECT_EQ(std::find(transported.begin(), transported.end(), "patm"),
            transported.end());
  EXPECT_EQ(std::find(transported.begin(), transported.end(), "pH"),
            transported.end());
  EXPECT_EQ(std::find(transported.begin(), transported.end(), "pe"),
            transported.end());
  EXPECT_EQ(std::find(transported.begin(), transported.end(), "Viscosity"),
            transported.end());
  // H, O, Charge should be in transported
  EXPECT_NE(std::find(transported.begin(), transported.end(), "H"),
            transported.end());
  EXPECT_NE(std::find(transported.begin(), transported.end(), "O"),
            transported.end());
  EXPECT_NE(std::find(transported.begin(), transported.end(), "Charge"),
            transported.end());

  // Test with no optional essentials
  PhreeqcMatrix pqc_mat_none(barite_db, barite_script, false, true, OE::None);
  const auto solution_names_none = pqc_mat_none.getSolutionNames();

  // Should only have mandatory essentials + chemical species (no optional
  // essentials)
  EXPECT_EQ(
      std::find(solution_names_none.begin(), solution_names_none.end(), "pH"),
      solution_names_none.end());
  EXPECT_EQ(
      std::find(solution_names_none.begin(), solution_names_none.end(), "pe"),
      solution_names_none.end());
  EXPECT_EQ(std::find(solution_names_none.begin(), solution_names_none.end(),
                      "SolVol"),
            solution_names_none.end());
  EXPECT_EQ(std::find(solution_names_none.begin(), solution_names_none.end(),
                      "MassH2O"),
            solution_names_none.end());
  EXPECT_EQ(std::find(solution_names_none.begin(), solution_names_none.end(),
                      "Viscosity"),
            solution_names_none.end());
  EXPECT_EQ(std::find(solution_names_none.begin(), solution_names_none.end(),
                      "Density"),
            solution_names_none.end());

  // Test with all optional essentials (default behavior)
  PhreeqcMatrix pqc_mat_all(barite_db, barite_script, false, true, OE::All);
  const auto solution_names_all = pqc_mat_all.getSolutionNames();

  // Should have all optional essentials
  EXPECT_NE(
      std::find(solution_names_all.begin(), solution_names_all.end(), "pH"),
      solution_names_all.end());
  EXPECT_NE(
      std::find(solution_names_all.begin(), solution_names_all.end(), "pe"),
      solution_names_all.end());
  EXPECT_NE(
      std::find(solution_names_all.begin(), solution_names_all.end(), "SolVol"),
      solution_names_all.end());
  EXPECT_NE(
      std::find(solution_names_all.begin(), solution_names_all.end(), "MassH2O"),
      solution_names_all.end());
  EXPECT_NE(std::find(solution_names_all.begin(), solution_names_all.end(),
                      "Viscosity"),
            solution_names_all.end());
  EXPECT_NE(std::find(solution_names_all.begin(), solution_names_all.end(),
                      "Density"),
            solution_names_all.end());
}
